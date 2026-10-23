"""Native startup-summary checks. Uses existing prescribed hydro fixtures and readers."""
import argparse
import html
import json
from pathlib import Path
import re
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[2]
sys.path[:0]=[str(ROOT/'tst/shocks-mc'),str(ROOT/'vis/python')]
from common import edit_input, provenance, sha, thermo
import material_curves


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--executable',type=Path,required=True)
    parser.add_argument('--baseline',type=Path,required=True)
    parser.add_argument('--workdir',type=Path,required=True)
    args=parser.parse_args();work=args.workdir.resolve();work.mkdir(parents=True,exist_ok=False)
    report=dict(runs=[],checks=[],source_sha256=sha(ROOT/'src/driver/driver.cpp'),
                runner_sha256=sha(__file__),hardware='local Serial/double')
    seeds={'particles/particle_type':'lagrangian_tracer','particles/pusher':'lagrangian_tracer',
        'particles/random_seed':11,'tracer_seed_initial/id':1,'tracer_seed_initial/count_per_event':128,
        'tracer_seed_initial/start_time':0,'tracer_seed_initial/seed':903,
        'tracer_seed_initial/region':'box','tracer_seed_initial/weight':'volume',
        'tracer_seed_initial/x1min':-.4,'tracer_seed_initial/x1max':.4,
        'tracer_seed_initial/x2min':-.4,'tracer_seed_initial/x2max':.4}
    def save():(work/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    def require(test,message):
        if not test:raise ValueError(message)
    def check(name,action):
        try:report['checks'].append(dict(id=name,status='pass',**(action() or {})))
        except Exception as exc:report['checks'].append(dict(id=name,status='fail',error=repr(exc)))
        save()
    def run(name,changes=None,fixture='smooth',baseline=False,particles=True,restart=None):
        directory=work/name;directory.mkdir()
        options={'shock_tracking/enabled':'false','time/nlim':0,'output1/dt':-1,
                 'output3/file_type':'prtcl_thermo_history','output3/source':'both' if particles else 'curves',
                 'output3/variables':'rho,p,v1,v2,v3','output3/dcycle':1}
        if particles:options.update(seeds)
        options.update(changes or {})
        input=directory/'input.athinput'
        input.write_text(edit_input((ROOT/('inputs/tests/shocks-curves/'+fixture+'.athinput')).read_text(),options))
        binary=(args.baseline if baseline else args.executable).resolve()
        command=[str(binary),'-r',str(restart)] if restart else [str(binary),'-i',str(input)]
        with (directory/'run.log').open('w') as log:
            status=subprocess.run(command,cwd=directory,stdout=log,stderr=subprocess.STDOUT,timeout=45).returncode
        report['runs'].append(dict(name=name,command=command,exit=status,
            input_used=not bool(restart),input_sha256=sha(input),
            restart_sha256=sha(restart) if restart else None,log_sha256=sha(directory/'run.log'),
            **provenance(binary)));save()
        require(status==0,name+': native failure; inspect run.log')
        log=(directory/'run.log').read_text()
        if baseline:return directory,''
        if name=='no_tracers':
            require('Tracer setup summary:' not in log,'summary printed for absent populations')
            return directory,''
        require(log.count('Tracer setup summary:')==1 and log.count('End tracer setup summary')==1,
                'summary missing or repeated')
        text=log.split('Tracer setup summary:',1)[1].split('End tracer setup summary',1)[0]
        if (directory/'prtcl_thermo_history').exists():
            rows,times=thermo(directory)
            count=sum(g==min(times) for g,p in rows)
            match=re.search(r'Particles: type=.*(?:actual|restored)=(\d+)',text)
            require(match and int(match[1])==count,'printed particle count differs from native records')
        snapshots=list(directory.glob('*.curves.initial.txt')) or list(directory.glob('*.curves.final.txt'))
        if snapshots:
            state=material_curves.read(snapshots[0])
            actual=sum(len(c['nodes']) for c in state['curves'])
            match=re.search(r'Curves: .*nodes(?:_actual)?=(\d+)',text)
            require(match and int(match[1])==actual,'printed node count differs from native geometry')
        (directory/'summary.txt').write_text('Tracer setup summary:'+text+'End tracer setup summary\n')
        return directory,text

    cases={}
    def counts(text):
        particle=re.search(r'requested_at_start=(\d+), actual=(\d+)',text)
        nodes=re.search(r'nodes_requested=(\d+), nodes_actual=(\d+)',text)
        return tuple(map(int,particle.groups()+nodes.groups()))
    def fresh():
        directory,text=run('fresh',{'time/nlim':4,'output1/dt':1})
        require('requested_at_start=128, actual=128' in text,'particle startup counts wrong')
        require('nodes_requested=32, nodes_actual=32' in text,'fixed node counts wrong')
        require('region=box, weight=volume' in text and 'cell centers' in text and
                'draws with replacement' in text and 'adaptation=false' in text,'placement/adaptation missing')
        require('source=both, every 1 cycles' in text and 'particle_shocks=false' in text,'history summary wrong')
        old,_=run('baseline',{'time/nlim':4,'output1/dt':1},baseline=True)
        require(next((directory/'prtcl_thermo_history').glob('*.thp')).read_bytes()==
                next((old/'prtcl_thermo_history').glob('*.thp')).read_bytes(),'summary changed particle history bytes')
        for filename in ['material.curves.initial.txt','material.curves.final.txt']:
            require((directory/filename).read_bytes()==(old/filename).read_bytes(),'summary changed curve state')
        cases['fresh']=counts(text)
        return dict(particle_history_and_curve_state_byte_equal=True)
    check('SUM01_fresh_and_passivity',fresh)

    def adaptive():
        directory,text=run('adaptive',fixture='amr')
        state=material_curves.read(directory/'material.curves.initial.txt')
        nodes=[n for c in state['curves'] for n in c['nodes']]
        requested=sum(n['anchor'] for n in nodes)
        require(len(nodes)>requested,'fixture did not insert nodes at startup')
        require(f'nodes_requested={requested}, nodes_actual={len(nodes)}' in text,'adapted counts wrong')
        require('adaptation=true' in text and 'max_spacing=auto(mesh)' in text and
                'chord_error=auto(spacing)' in text and 'unresolved_curves=0' in text,'adaptive settings missing')
        cases['adaptive']=counts(text)
        return dict(requested_nodes=requested,actual_nodes=len(nodes))
    check('SUM02_initial_adaptation',adaptive)

    def empty_region():
        _,text=run('empty_region',{'tracer_seed_initial/target':'rho','tracer_seed_initial/target_min':2})
        require('requested_at_start=128, actual=0' in text and 'filter=rho min=2' in text,
                'empty eligible selection hidden/misreported')
        cases['empty region']=counts(text)
    check('SUM03_no_eligible_cells',empty_region)

    def delayed():
        _,text=run('delayed',{'tracer_seed_initial/start_time':.2})
        require('requested_at_start=0, actual=0' in text and 'requested/event=128' in text and
                'start=0.2' in text,'future injection counted as missing startup particles')
        cases['future birth']=counts(text)
    check('SUM04_future_birth',delayed)

    def circle():
        _,text=run('circle',{'curve1/type':'circle','curve1/count':128,
            'curve1/x1':-.23,'curve1/x2':-.23,'curve1/radius':.14})
        require('type=circle, closed=true, placement=uniform angles' in text,'loop placement missing')
        require('nodes_requested=128, nodes_actual=128' in text,'circle counts wrong')
        cases['circle']=counts(text)
    check('SUM05_loop_placement',circle)

    def polyline():
        _,text=run('polyline',{'curve1/type':'polyline','curve1/coordinates':'lifted',
            'curve1/points':'-0.4 -0.3 0; 0 0 0; 0.4 0.3 0'})
        require('placement=explicit lifted vertices' in text and 'nodes_requested=3, nodes_actual=3' in text,
                'explicit polyline vertices not counted correctly')
    check('SUM06_polyline_placement',polyline)

    def restart():
        checkpoint=sorted((work/'fresh/rst').glob('*.rst'))[-1]
        directory,text=run('restart',restart=checkpoint)
        require('restart, time=' in text and 'restored=128' in text and 'Curves: restored=1, nodes=32' in text,
                'restored counts wrong')
        require('requested_at_start=' not in text and 'seed requests not reapplied' in text,
                'restart described as fresh seeding')
        require((directory/'material.curves.final.txt').read_bytes()==
                (work/'fresh/material.curves.final.txt').read_bytes(),'restart summary changed geometry')
    check('SUM07_restart_restored_counts',restart)

    def histories():
        _,text=run('histories',{'shock_tracking/enabled':'true','problem/flow':'shock_tube',
                               'output2/dcycle':0,'output2/dt':-1},fixture='planar')
        require('particle_shocks=true, curve_shocks=true, particle_shock_trace=true, curve_shock_trace=true' in text,
                'enabled shock histories not reported')
        _,off=run('disabled_histories',{'output3/dcycle':0})
        require('thermo histories: disabled' in off and 'prtcl_thermo_history' not in off,
                'disabled output configuration reported active')
        _,terminal=run('terminal_history',{'output3/dcycle':-1})
        require('source=both, initial/final only' in terminal,'terminal-only history cadence misreported')
    check('SUM08_enabled_histories',histories)

    def empty_and_separate():
        _,text=run('empty',{'material_curves/ncurves':0,'tracer_seed_initial/count_per_event':0})
        require('requested_at_start=0, actual=0' in text and 'nodes_requested=0, nodes_actual=0' in text,
                'empty populations not reported correctly')
        _,text=run('curve_only',particles=False)
        require('Particles: disabled' in text and 'source=curves' in text,'curve-only summary wrong')
        _,text=run('particles_only',{'material_curves/enabled':'false','output3/source':'particles'})
        require('Curves: disabled' in text and 'source=particles' in text,'particles-only summary wrong')
        run('no_tracers',{'material_curves/enabled':'false','output3/dcycle':0},particles=False)
    check('SUM09_empty_separate_and_absent',empty_and_separate)
    report['counts']=cases;save()
    print(json.dumps(report['checks'],indent=2))
    if not all(c['status']=='pass' for c in report['checks']):return 1

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,3,figsize=(15,5))
    states=material_curves.read(work/'circle/material.curves.initial.txt')
    cloud,_=thermo(work/'circle'); cloud=[s for (g,p),s in cloud.items() if g==0]
    axes[0].scatter([s['x'][0] for s in cloud],[s['x'][1] for s in cloud],s=12,alpha=.6,label='128 regular particles')
    nodes=states['curves'][0]['nodes'];closed=nodes+[nodes[0]]
    axes[0].plot([n['x'][0] for n in closed],[n['x'][1] for n in closed],'.-',ms=3,color='tab:orange',label='128 loop nodes')
    axes[0].set(xlabel=r'$x$ (code units)',ylabel=r'$y$ (code units)',title='Native initial placement',aspect='equal')
    axes[0].legend(fontsize=8)
    labels=list(cases);x=list(range(len(labels)))
    for ax,a,b,title in [(axes[1],0,1,'Regular particles at startup'),(axes[2],2,3,'Curve nodes after startup adaptation')]:
        ax.bar([i-.18 for i in x],[cases[k][a] for k in labels],width=.36,label='Requested')
        ax.bar([i+.18 for i in x],[cases[k][b] for k in labels],width=.36,label='Actual')
        ax.set(xticks=x,xticklabels=labels,ylabel='Count',title=title)
        ax.tick_params(axis='x',rotation=35);ax.legend(fontsize=8)
    fig.suptitle('Startup summary checked against actual native populations')
    fig.text(.5,.01,'Requested particles refer to births due at startup. No eligible cells gives 128 requested / 0 actual; future births give 0 / 0.\n'
             'Adaptation increases the initial curve from 8 to 15 nodes. Printing leaves particle-history bytes and curve state unchanged.',ha='center',fontsize=9)
    fig.tight_layout(rect=(0,.11,1,.94));fig.savefig(work/'startup-summary.png',dpi=170);plt.close(fig)
    summaries=''.join('<details><summary>'+html.escape(p.parent.name)+'</summary><pre>'+html.escape(p.read_text())+'</pre></details>'
                      for p in sorted(work.glob('*/summary.txt')))
    (work/'index.html').write_text('<!doctype html><meta charset="utf-8"><title>Tracer startup summary</title>'
        '<style>body{font:17px system-ui;max-width:1300px;margin:30px auto}img{width:100%}pre{white-space:pre-wrap;background:#f4f6f8;padding:15px}details{margin:15px}</style>'
        '<h1>Tracer startup summary</h1><p>Actual local Serial/double native runs; logging-only change. No scientific or MPI/GPU qualification added.</p>'
        '<img src="startup-summary.png">'+summaries)
    return 0


if __name__=='__main__':raise SystemExit(main())
