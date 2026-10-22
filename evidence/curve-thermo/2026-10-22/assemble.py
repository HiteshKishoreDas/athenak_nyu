"""Seal this exact local thermo campaign, source, binaries, failed fixtures and plots."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import tarfile

ROOT=Path(__file__).resolve().parents[3]


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024),b''):h.update(chunk)
    return h.hexdigest()


def write(path,value):path.write_text(json.dumps(value,indent=2)+'\n')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args(); out=args.output.resolve()
    assert (out/'index.html').is_file() and not (out/'manifest.json').exists()
    runs=ROOT/'build-lagrangian/evidence'
    reports={name:json.loads((runs/name/'report.json').read_text()) for name in
             ('curve-thermo-reviewed-code','thermo-reviewed-LT','thermo-reviewed-CSH')}
    for name,count in zip(reports,(18,14,20)):
        assert len(reports[name]['checks'])==count
        assert all(c['status']=='pass' for c in reports[name]['checks'])
    parity=json.loads((runs/'thermo-reviewed-parity.json').read_text())
    assert parity['status']=='pass' and len(parity['cases'])==53
    commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    cmake=(ROOT/'CMakeLists.txt').read_text()
    section=cmake[cmake.index('set(shock_source_files'):cmake.index('set(shock_source_hashes')]
    sources=re.findall(r'src/[A-Za-z0-9_./]+',section)
    assert len(set(sources))==len(sources)
    hashes={p:sha(ROOT/p) for p in sources}
    bindings={}
    for directory,problem in [('build-lagrangian','built_in_pgens'),
                              ('build-lagrangian-prescribed','material_curve_test')]:
        build=ROOT/directory
        computed=hashlib.sha256(('problem:'+problem+';'+''.join(p+':'+hashes[p]+';' for p in sources)).encode()).hexdigest()
        header=(build/'shock_version.hpp').read_text()
        assert re.search(r'SHOCK_IMPLEMENTATION_HASH "([a-f0-9]{64})"',header)[1]==computed
        binary=(build/'src/athena').read_bytes(); assert computed.encode() in binary
        bindings[directory]=dict(problem=problem,compiled_source_fingerprint=computed,
                                  executable_sha256=sha(build/'src/athena'))
        destination=out/'runtime'/directory;destination.mkdir(parents=True)
        for path in ['src/athena','CMakeCache.txt','config.hpp','shock_version.hpp',
                     'src/CMakeFiles/athena.dir/flags.make']:
            shutil.copy2(build/path,destination/Path(path).name)
    for r in reports.values():
        for row in r['runs']:
            assert row['executable_sha256'] in {b['executable_sha256'] for b in bindings.values()}
    write(out/'source-binding.json',dict(source_commit=commit,source_order=sources,
        sources_sha256=hashes,builds=bindings,kokkos_commit=subprocess.check_output(
            ['git','rev-parse','HEAD'],cwd=ROOT/'kokkos',text=True).strip()))
    source_paths=[]
    for prefix in ['src','kokkos','tst','inputs','vis/python','evidence/shocks/reference']:
        source_paths.extend(p for p in (ROOT/prefix).rglob('*') if p.is_file()
            and '.git' not in p.parts and '__pycache__' not in p.parts)
    source_paths.extend(ROOT/p for p in ['CMakeLists.txt','config.hpp.in','LICENSE',
        'docs/curve-thermo-history.md','evidence/curve-thermo/2026-10-22/reviews.md',
        'evidence/curve-thermo/2026-10-22/ledger.md','evidence/curve-thermo/2026-10-22/assemble.py'])
    source_paths=sorted(set(source_paths))
    with tarfile.open(out/'source.tar.gz','w:gz') as archive:
        for p in source_paths:archive.add(p,arcname=str(p.relative_to(ROOT)),recursive=False)
    write(out/'source-files.json',{str(p.relative_to(ROOT)):sha(p) for p in source_paths})
    raw_hashes={}
    directories=['curve-thermo-reviewed-code','thermo-reviewed-LT','thermo-reviewed-CSH',
                 'curve-thermo-dt-repro','curve-thermo-overdue-repro']
    with tarfile.open(out/'raw-native.tar.gz','w:gz') as archive:
        for directory in directories:
            for p in sorted((runs/directory).rglob('*')):
                if p.is_file():
                    relative=str(p.relative_to(runs));raw_hashes[relative]=sha(p)
                    archive.add(p,arcname=relative,recursive=False)
    # Verify archived bytes, not just the source inputs to the packer.
    with tarfile.open(out/'raw-native.tar.gz') as archive:
        for member in archive:
            assert hashlib.sha256(archive.extractfile(member).read()).hexdigest()==raw_hashes[member.name]
    write(out/'raw-files.json',raw_hashes)
    for name in reports:write(out/(name+'-report.json'),reports[name])
    write(out/'detector-parity.json',parity)
    for suffix in ['builtin-final','prescribed-final','built-final2','prescribed-final2']:
        for stage in ['configure','build']:
            p=Path('/private/tmp')/('curve-thermo-'+suffix+'-'+stage+'.log')
            if p.exists():shutil.copy2(p,out/p.name)
    shutil.copy2(ROOT/'docs/curve-thermo-history.md',out/'README.md')
    write(out/'verification.json',dict(native_checks=52,detector_parity_cases=53,
        source_fingerprints_recomputed=True,archived_raw_files_verified=len(raw_hashes),
        source_files=len(source_paths),scientific_approval='pending',hardware='local Serial/double'))
    write(out/'manifest.json',dict(source_commit=commit,files={str(p.relative_to(out)):sha(p)
        for p in sorted(out.rglob('*')) if p.is_file() and p.name!='manifest.json'}))
    print(json.dumps(dict(source_commit=commit,source_files=len(source_paths),
        archived_raw_files=len(raw_hashes),manifest_sha256=sha(out/'manifest.json'))))


if __name__=='__main__':main()
