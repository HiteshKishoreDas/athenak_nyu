"""One analysis entry point for native particle, curve and encounter outputs.

read_run() returns independent sampled streams; it never resamples them onto a
shared clock. Particle and curve-node namespaces remain separate. Legacy files
and opt-in native binary tables use the same existing semantic readers.
"""
from pathlib import Path
import re

import curve_shock_history
import material_curves
import shock_history
import tracer_binary
from tracer_thermo_history import read_curves, read_particles


def read_run(directory, thermo_stem=None, shock_directory='shock_events', include_particles=True):
    """Read independent streams; include_particles=False avoids unrelated particle IO."""
    root = Path(directory)
    if not root.is_dir():
        raise ValueError('run directory does not exist')
    result = dict(particles={},curve_nodes={},particle_times={},curve_times={},
                  geometry=[],particle_events=[],curve_events=[],
                  particle_suspended=[],curve_suspended=[],available=[])
    paths = sorted((root/'prtcl_thermo_history').glob('*.thp')) if include_particles else []
    if thermo_stem is not None:
        paths = [p for p in paths if p.stem == thermo_stem]
    if len(paths)>1:
        raise ValueError('multiple particle thermo streams; select thermo_stem')
    if paths:
        result['particles'],result['particle_times'] = read_particles(paths[0],allow_identical_replay=True)
        result['available'].append('particles')
    curve_thermo = root/'curve_thermo_history'
    if curve_thermo.is_dir():
        result['curve_nodes'],result['curve_times'] = read_curves(curve_thermo,stem=thermo_stem)
        result['available'].append('curve_nodes')
    for carrier, module, where in (('particle',shock_history,root/shock_directory),
                                  ('curve',curve_shock_history,root/shock_directory/'curves')):
        if carrier=='particle' and not include_particles:
            continue
        manifests = list(where.glob('manifest.rank*'))
        chunks = [p for pattern in ('cycle*.csv','cycle*.tbin','trace*.csv','trace*.tbin',
                                  'birthtrace*.csv','birthtrace*.tbin','suspended*.csv','suspended*.tbin')
                  for p in where.glob(pattern)]
        if chunks and not manifests:
            raise ValueError('orphan '+carrier+' shock chunks: missing manifests')
        if manifests:
            result[carrier+'_events'] = module.read_history(where)
            result[carrier+'_suspended'] = module.read_suspended(where)
            result['available'].append(carrier+'_events')
    # Native geometry is per rank; legacy geometry is one collective snapshot.
    groups = {}
    for folder in (root,root/'curves'):
        pattern='*.curves.*.rank*.tbin' if folder==root else '*.[0-9]*.rank*.tbin'
        for path in sorted(folder.glob(pattern)):
            if folder!=root and not re.fullmatch(r'.+\.[0-9]{8,}\.rank[0-9]+\.tbin',path.name):
                continue
            stem, sep, suffix = path.name.rpartition('.rank')
            if not sep or not suffix[:-5].isdigit():
                raise ValueError('malformed native geometry filename')
            state = tracer_binary.read_geometry(path)
            if int(suffix[:-5]) != state['rank']:
                raise ValueError('geometry filename rank mismatch')
            key = str(folder/stem)
            ranks = groups.setdefault(key,{})
            if state['rank'] in ranks:
                raise ValueError('duplicate geometry rank')
            ranks[state['rank']] = state
        pattern = '*.curves.*.txt' if folder == root else '*.txt'
        for path in sorted(folder.glob(pattern)):
            result['geometry'].append(dict(path=str(path),**material_curves.read(path)))
    for key,ranks in sorted(groups.items()):
        state = ranks[min(ranks)]
        if set(ranks) != set(range(state['nranks'])):
            raise ValueError('missing geometry rank')
        for other in ranks.values():
            if any(other[k] != state[k] for k in ('cycle','time','domain','dimensions','nranks')):
                raise ValueError('geometry rank metadata mismatch')
        merged = dict(state,path=key+'.rank*.tbin')
        merged['curves'] = sorted([c for s in ranks.values() for c in s['curves']],key=lambda c:c['id'])
        merged['curve_count'] = len(merged['curves'])
        merged['node_count'] = sum(s['node_count'] for s in ranks.values())
        result['geometry'].append(merged)
    result['geometry'].sort(key=lambda s:(s['cycle'],s['path']))
    if result['geometry']:
        result['available'].append('geometry')
    return result
