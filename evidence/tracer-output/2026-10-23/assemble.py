"""Seal exact source, runtimes, final outputs and retained failures beside plots."""
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
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):h.update(block)
    return h.hexdigest()


def write(path,value):path.write_text(json.dumps(value,indent=2)+'\n')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();out=args.output.resolve()
    if not (out/'index.html').exists() or (out/'manifest.json').exists():
        raise ValueError('use an existing, unsealed visualization directory')
    commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    reports={}
    runpaths={'output':ROOT/'build-lagrangian/evidence/native-output-qualified-final-20261023',
        'thermo':ROOT/'build-lagrangian/evidence/native-output-frozen-CT',
        'curves':ROOT/'build-lagrangian/evidence/native-output-frozen-CSH',
        'lagrangian':ROOT/'build-lagrangian/evidence/native-output-frozen-LT',
        'setup':ROOT/'build-lagrangian-prescribed/evidence/native-output-frozen-setup'}
    for name,count in [('output',21),('thermo',18),('curves',20),('lagrangian',14),('setup',9)]:
        reports[name]=json.loads((runpaths[name]/'report.json').read_text())
        assert len(reports[name]['checks'])==count
        assert all(c['status']=='pass' for c in reports[name]['checks'])
        write(out/(name+'-report.json'),reports[name])
    cmake=(ROOT/'CMakeLists.txt').read_text()
    section=cmake[cmake.index('set(shock_source_files'):cmake.index('set(shock_source_hashes')]
    sources=re.findall(r'src/[A-Za-z0-9_./]+',section)
    assert len(set(sources))==len(sources)
    hashes={p:sha(ROOT/p) for p in sources};bindings={}
    for directory,problem in [('build-lagrangian','built_in_pgens'),
                              ('build-lagrangian-prescribed','material_curve_test')]:
        build=ROOT/directory
        fingerprint=hashlib.sha256(('problem:'+problem+';'+''.join(
            p+':'+hashes[p]+';' for p in sources)).encode()).hexdigest()
        assert re.search(r'SHOCK_IMPLEMENTATION_HASH "([a-f0-9]{64})"',
                         (build/'shock_version.hpp').read_text())[1]==fingerprint
        assert fingerprint.encode() in (build/'src/athena').read_bytes()
        bindings[directory]=dict(problem=problem,compiled_source_fingerprint=fingerprint,
                                  executable_sha256=sha(build/'src/athena'))
        destination=out/'runtime'/directory;destination.mkdir(parents=True)
        for name in ['src/athena','CMakeCache.txt','config.hpp','shock_version.hpp',
                     'src/CMakeFiles/athena.dir/flags.make']:
            shutil.copy2(build/name,destination/Path(name).name)
    for report in reports.values():
        for run in report['runs']:
            assert run['executable_sha256'] in {b['executable_sha256'] for b in bindings.values()}
    write(out/'source-binding.json',dict(source_commit=commit,source_order=sources,
        sources_sha256=hashes,builds=bindings,kokkos_commit=subprocess.check_output(
            ['git','rev-parse','HEAD'],cwd=ROOT/'kokkos',text=True).strip()))
    source_paths=[]
    for prefix in ['src','kokkos','tst','inputs','vis/python','evidence/shocks/reference']:
        source_paths.extend(p for p in (ROOT/prefix).rglob('*') if p.is_file()
            and '.git' not in p.parts and '__pycache__' not in p.parts)
    source_paths.extend(ROOT/p for p in ['CMakeLists.txt','config.hpp.in','LICENSE',
        'docs/tracer-output-unification.md','evidence/tracer-output/2026-10-23/ledger.md',
        'evidence/tracer-output/2026-10-23/reviews.md','evidence/tracer-output/2026-10-23/assemble.py'])
    source_paths=sorted(set(source_paths))
    with tarfile.open(out/'source.tar.gz','w:gz') as archive:
        for p in source_paths:archive.add(p,arcname=str(p.relative_to(ROOT)),recursive=False)
    write(out/'source-files.json',{str(p.relative_to(ROOT)):sha(p) for p in source_paths})
    failures={'initial-native-smoke':ROOT/'build-lagrangian/evidence/native-output-smoke-20261023',
              'startup-wrong-python':ROOT/'build-lagrangian-prescribed/evidence/native-output-final-setup',
              'birthtrace-harness-assumption':ROOT/'build-lagrangian/evidence/native-output-qualified-20261023'}
    raw_hashes={}
    with tarfile.open(out/'raw-native.tar.gz','w:gz') as archive:
        for label,directory in {**runpaths,**failures}.items():
            for path in sorted(directory.rglob('*')):
                if path.is_file():
                    relative=label+'/'+str(path.relative_to(directory));raw_hashes[relative]=sha(path)
                    archive.add(path,arcname=relative,recursive=False)
    for filename,expected in [('raw-native.tar.gz',raw_hashes),
                              ('source.tar.gz',{str(p.relative_to(ROOT)):sha(p) for p in source_paths})]:
        with tarfile.open(out/filename) as archive:
            names=[]
            for member in archive:
                names.append(member.name)
                assert hashlib.sha256(archive.extractfile(member).read()).hexdigest()==expected[member.name]
            assert set(names)==set(expected) and len(names)==len(expected)
    write(out/'raw-files.json',raw_hashes)
    for path in Path('/private/tmp').glob('native-output-frozen-*.log'):
        shutil.copy2(path,out/path.name)
    for name in ['native-output-qualified-final-check.log','native-output-final-setup.log']:
        path=Path('/private/tmp')/name
        if path.exists():shutil.copy2(path,out/path.name)
    shutil.copy2(ROOT/'docs/tracer-output-unification.md',out/'README.md')
    shutil.copy2(ROOT/'evidence/tracer-output/2026-10-23/reviews.md',out/'reviews.md')
    write(out/'verification.json',dict(native_check_groups=82,
        source_fingerprints_recomputed=True,archived_raw_files_verified=len(raw_hashes),
        archived_source_files_verified=len(source_paths),scientific_approval='pending',
        hardware='local Serial/double; MPI/GPU/cluster not qualified'))
    page=out/'index.html'
    page.write_text(page.read_text().replace('<h1>AthenaK: native particle, curve and shock outputs</h1>',
        '<h1>AthenaK: native particle, curve and shock outputs</h1>'
        '<p>Final integrated evidence: 82 native check groups passed (21 new output checks, '
        '18 thermo, 20 curve, 14 continuous-tracer and 9 startup checks). '
        '<a href="README.md">Setup and file-format documentation</a> · '
        '<a href="source-binding.json">Exact code/build binding</a> · '
        '<a href="verification.json">Bundle verification</a>.</p>'))
    write(out/'manifest.json',dict(source_commit=commit,files={str(p.relative_to(out)):sha(p)
        for p in sorted(out.rglob('*')) if p.is_file() and p.name!='manifest.json'}))
    print(json.dumps(dict(source_commit=commit,source_files=len(source_paths),
        archived_raw_files=len(raw_hashes),manifest_sha256=sha(out/'manifest.json'))))


if __name__=='__main__':main()
