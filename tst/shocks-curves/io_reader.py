"""Run: python3 tst/shocks-curves/io_reader.py (standard library only)."""
import csv
import io
from pathlib import Path
import sys
import tempfile

sys.path.insert(0, str(Path(__file__).resolve().parents[2]/'vis/python'))
from curve_shock_history import read_history, read_suspended, read_trace
from shock_history import HistoryError, parse_event


NODE = 9007199254740993
SEQUENCE = 9007199254740995


def observation(time, generation, inside, valid=True):
    return list(map(str, [time, generation, int(valid), int(inside),
                         .25, .5, .75, 0, 0, 0,
                         2 if inside else 1, 2 if inside else 1, 1,
                         2, 2, 1, 0, 0, 1, 1, 1, 2, 2, 1]))


def row(sequence=SEQUENCE):
    core = list(map(str, [NODE, sequence, 0, 24]))
    for values in (observation(0, 0, False), observation(1, 1, True),
                   observation(2, 2, True), observation(3, 3, False),
                   observation(1, 1, True)):
        core += values
    assert len(core) == 124
    return core+list(map(str, [0, 3, .125, .125, .125, 4, .0625, .0625, .0625]))


def csv_line(values):
    out = io.StringIO()
    csv.writer(out, lineterminator='\n').writerow(values)
    return out.getvalue()


def write_chunk(root, generation, rows=(), suspended=False):
    prefix = 'suspended' if suspended else 'cycle'
    kind = 'SUSPENDED' if suspended else 'EVENTS'
    path = root/f'{prefix}{generation}.rank0.csv'
    path.write_text(f'# AKCURVE-SHOCK-{kind} 1 test {generation} {len(rows)}\n'
                    '# exact-config\n'+''.join(csv_line(r) for r in rows)+f'# END {len(rows)}\n')
    return path


def rejected(call):
    try:
        call()
    except HistoryError:
        return
    raise AssertionError('malformed history accepted')


def main():
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        (root/'manifest.rank0').write_text('AKCURVE-SHOCK-EVENTS 1\ntest\nexact-config\n')
        for generation in range(4):
            write_chunk(root, generation, [row()] if generation == 3 else [])
        events = read_history(root)
        assert len(events) == 1
        event = events[0]
        assert event['node_id'] == event['particle_id'] == NODE
        assert event['curve_id'] == NODE >> 32 and event['sequence'] == SEQUENCE
        assert event['identity'] == ('curve_node', 'test', NODE, SEQUENCE)
        assert event['identity'] != parse_event(row()[:124], 'test', 3)['identity']
        assert event['birth'] == 0 and event['last_level'] == 4
        assert event['first_spacing'] == [.125]*3
        core = parse_event(row()[:124], 'test', 3)
        assert all(event[key] == value for key, value in core.items() if key != 'identity')
        assert read_history(root, through_cycle=2) == []
        assert read_suspended(root) == []
        censored = row(SEQUENCE+1)
        censored[2:4] = ['1', '93']  # left/right censored, merged, cadence/no plateau
        censored[4:28] = observation(-1, 0, False, valid=False)
        censored[76:100] = observation(-1, 0, False, valid=False)
        censored[124] = '1'
        write_chunk(root, 3, [censored], suspended=True)
        views = read_suspended(root)
        assert len(views) == 1 and views[0]['sequence'] == SEQUENCE+1
        assert not views[0]['pre']['valid'] and not views[0]['post']['valid']
        assert len(read_history(root)) == 1  # suspended views never become events
        trace = list(map(str, [NODE, SEQUENCE+1, 1, 21]))+observation(2, 2, True)
        trace += ['1', '4', '.0625', '.0625', '.0625']
        path = root/'trace2.rank0.csv'
        path.write_text(csv_line(trace))
        samples = read_trace(path)
        assert samples[NODE]['birth'] == 1 and samples[NODE]['phase'] == 1
        assert samples[NODE]['spacing'] == [.0625]*3
        birthpath = root/'birthtrace2.rank0.csv'
        birthpath.write_text(csv_line(trace))
        assert read_trace(birthpath) == samples
        path.write_text(csv_line(trace)*2)
        rejected(lambda: read_trace(path))
        path.write_text(csv_line(trace[:-1]))
        rejected(lambda: read_trace(path))
        for index, value in [(0, '0'), (1, str(1 << 64)), (3, '2048'),
                             (124, 'nan'), (124, '2'), (125, '-1'),
                             (126, '0'), (130, 'inf')]:
            broken = row();broken[index] = value
            write_chunk(root, 3, [broken])
            rejected(lambda: read_history(root))
        write_chunk(root, 3, [row()[:-1]])
        rejected(lambda: read_history(root))
        write_chunk(root, 3, [row(), row()])
        rejected(lambda: read_history(root))
        write_chunk(root, 3, [row()])
        (root/'cycle1.rank0.csv').unlink()
        rejected(lambda: read_history(root))
        write_chunk(root, 1)
        duplicate = root/'manifest.rank00'
        duplicate.write_bytes((root/'manifest.rank0').read_bytes())
        rejected(lambda: read_history(root))
        duplicate.unlink()
        write_chunk(root, 3, [row()], suspended=True)
        rejected(lambda: read_suspended(root))
        assert len(read_history(root)) == 1
    with tempfile.TemporaryDirectory() as directory:
        root=Path(directory)
        # User parameter text must never select the ID encoding.
        legacy='curve-node-history-v1 HASH parameter 15: curve-id=8/24'
        (root/'manifest.rank0').write_text('AKCURVE-SHOCK-EVENTS 1\ntest\n'+legacy+'\n')
        for generation in range(4):
            path=write_chunk(root,generation,[row()] if generation==3 else [])
            path.write_text(path.read_text().replace('exact-config',legacy))
        assert read_history(root)[0]['curve_id']==NODE>>32
        current='curve-node-history-v2 HASH'
        (root/'manifest.rank0').write_text('AKCURVE-SHOCK-EVENTS 1\ntest\n'+current+'\n')
        for path in root.glob('cycle*.csv'):path.write_text(path.read_text().replace(legacy,current))
        rejected(lambda:read_history(root))  # old wide ID cannot alias into current namespace
        compact=row();compact[0]=str((1<<32)-1)
        path=write_chunk(root,3,[compact]);path.write_text(path.read_text().replace('exact-config',current))
        assert read_history(root)[0]['curve_id']==255
        # Variable-width IDs are never guessed from their numeric value.
        (root/'manifest.rank0').write_text('AKCURVE-SHOCK-EVENTS 1\ntest\ncurve-node-history-v3 4 HASH\n')
        trace=list(map(str,[(15<<28)|1,0,1,21]))+observation(2,2,True)
        trace+=['1','4','.0625','.0625','.0625']
        path=root/'trace2.rank0.csv';path.write_text(csv_line(trace))
        assert read_trace(path)[(15<<28)|1]['curve_id']==15
        rejected(lambda:read_trace(path,curve_id_bits=8))
        (root/'manifest.rank0').unlink()
        rejected(lambda:read_trace(path))
        assert read_trace(path,curve_id_bits=4)[(15<<28)|1]['curve_id']==15
        # Empty histories must validate layout metadata too.
        empty=root/'empty';empty.mkdir()
        (empty/'manifest.rank0').write_text('AKCURVE-SHOCK-EVENTS 1\ntest\ncurve-node-history-v3 31 HASH\n')
        path=write_chunk(empty,0,[])
        path.write_text(path.read_text().replace('exact-config','curve-node-history-v3 31 HASH'))
        rejected(lambda:read_history(empty))
        rejected(lambda:read_suspended(empty))
    print('Curve IO current uint32/historical uint64, metadata, censoring and adversarial checks PASS')


if __name__ == '__main__':
    main()
