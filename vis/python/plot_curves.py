#!/usr/bin/env python3
"""Plot a curve's geometry, thermo profile, node history and shock evidence.

Uses existing supplied readers, with particle loading disabled. Samples
retain their own clocks. No resampling, parent-history inheritance or inferred
shock crossing times. Requires the existing NumPy/Matplotlib environment.
"""
import argparse
import json
from pathlib import Path

from tracer_history import read_run
from material_curves import node_curve_id


def select_curve(data, curve_id=None, node_id=None, field='pressure'):
    """Select exact integer identities and independent recorded streams."""
    if node_id is not None:
        try:
            inferred = node_curve_id(node_id)
        except ValueError as exc:
            raise ValueError('node ID must be an exact valid current uint32 or historical uint64 ID') from exc
        if curve_id is not None and inferred != curve_id:
            raise ValueError('node ID belongs to a different curve')
        curve_id = inferred
    if curve_id is not None and not 0 < curve_id < 1 << 32:
        raise ValueError('curve ID must be in [1, 4294967295]')
    available = {c['id'] for s in data['geometry'] for c in s['curves']}
    available.update(o['curve_id'] for o in data['curve_nodes'].values())
    available.update(e['curve_id'] for key in ('curve_events','curve_suspended') for e in data[key])
    if not available:
        raise ValueError('no curve geometry, thermo or encounter records found')
    if curve_id is None:
        if len(available)!=1:
            raise ValueError('choose --curve-id from '+str(sorted(available)))
        curve_id = next(iter(available))
    if curve_id not in available:
        raise ValueError('curve ID not found; available curves: '+str(sorted(available)))
    frames = {}
    for state in data['geometry']:
        for curve in state['curves']:
            if curve['id']!=curve_id:
                continue
            clock = state['cycle'],state['time']
            if clock in frames and json.dumps(frames[clock]['curve'],sort_keys=True)!=json.dumps(curve,sort_keys=True):
                raise ValueError('conflicting same-clock geometry snapshots for selected curve')
            frames[clock] = dict(cycle=state['cycle'],time=state['time'],curve=curve,
                                 domain=state['domain'],dimensions=state['dimensions'])
    frames = [frame for _,frame in sorted(frames.items())]
    samples = [o for o in data['curve_nodes'].values() if o['curve_id']==curve_id]
    fields = sorted({name for o in samples for name in o['fields']})
    if samples and field not in fields:
        raise ValueError('field not saved for this curve; available fields: '+', '.join(fields))
    events = {key:[e for e in data[key] if e['curve_id']==curve_id]
              for key in ('curve_events','curve_suspended')}
    ids = {n['id'] for f in frames for n in f['curve']['nodes']}
    ids.update(o['node_id'] for o in samples)
    ids.update(e['node_id'] for records in events.values() for e in records)
    if node_id is None:
        node_id = min(ids)
    if node_id not in ids:
        raise ValueError('node ID not found in selected curve')
    history = sorted((o for o in samples if o['node_id']==node_id),key=lambda o:o['cycle'])
    latest=max((o['cycle'] for o in samples),default=None)
    profile = sorted((o for o in samples if o['cycle']==latest),key=lambda o:o['a'])
    return dict(curve_id=curve_id,node_id=node_id,field=field,fields=fields,
                frames=frames,history=history,profile=profile,events=events,
                tracking_available='curve_events' in data['available'])


def curve_points(frame):
    """Close a lifted curve with its stored winding, never nearest-image wrapping."""
    points = [list(n['x']) for n in frame['curve']['nodes']]
    if frame['curve']['closed']:
        points.append([points[0][d]+frame['curve']['winding'][d]*
                       (frame['domain'][2*d+1]-frame['domain'][2*d]) for d in range(3)])
    return points


def plot_curve(selected, output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    output = Path(output)
    if output.exists():
        raise ValueError('output already exists; choose a new --output path')
    frames,history,profile = (selected[k] for k in ('frames','history','profile'))
    curve_id,node_id,field = (selected[k] for k in ('curve_id','node_id','field'))
    three_d = any(f['dimensions'][2]>1 for f in frames)
    fig = plt.figure(figsize=(13,9),layout='constrained')
    ax = fig.add_subplot(2,2,1,projection='3d' if three_d else None)
    if frames:
        for index in sorted({0,len(frames)//2,len(frames)-1}):
            frame=frames[index]; points=curve_points(frame)
            axes=list(zip(*points))
            ax.plot(*axes[:3 if three_d else 2],'.-',ms=3,
                    label=f"cycle {frame['cycle']}, t={frame['time']:.5g}")
        ax.legend(fontsize=8)
    else:
        ax.text(.1,.5,'No geometry snapshots saved',transform=ax.transAxes)
    ax.set(xlabel=r'lifted $x_1$ [code units]',ylabel=r'lifted $x_2$ [code units]',
           title='Stored geometry; edges follow connectivity')
    if three_d:
        ax.set_zlabel(r'lifted $x_3$ [code units]')
    else:
        ax.set_aspect('equal',adjustable='datalim')
    ax=fig.add_subplot(2,2,2)
    if profile:
        ax.plot([o['a'] for o in profile],[o['fields'][field] for o in profile],'.-')
        ax.set_title(f"Curve profile at cycle {profile[0]['cycle']}, t={profile[0]['time']:.5g}")
    else:
        ax.text(.1,.5,'No curve thermo history saved',transform=ax.transAxes)
        ax.set_title('Curve thermo profile unavailable')
    ax.set(xlabel=r'immutable material label $a$ [code units]',ylabel=field+' [code units]')
    ax=fig.add_subplot(2,2,3)
    if history:
        ax.plot([o['time'] for o in history],[o['fields'][field] for o in history],'.-',ms=3)
        birth=history[0]['birth']
        ax.axvline(birth,color='gray',ls=':',label=f'actual birth t={birth:.5g}')
        ax.legend(fontsize=8)
    else:
        ax.text(.1,.5,'No thermo samples for this node',transform=ax.transAxes)
    ax.set(xlabel='simulation time [code units]',ylabel=field+' [code units]',
           title='Selected node: saved samples; dotted line marks birth')
    ax=fig.add_subplot(2,2,4)
    labels=[]
    for key,terminal in [('curve_events',False),('curve_suspended',True)]:
        for e in selected['events'][key]:
            if e['node_id']!=node_id:
                continue
            row=len(labels)
            label=('terminal view' if terminal else 'complete' if e['status']==0 else 'partial')
            labels.append(f"{label}, seq {e['sequence']}, flags {e['flags']}")
            color='tab:purple' if terminal else 'tab:blue' if e['status']==0 else 'tab:orange'
            first,last=e['first']['time'],e['last']['time']
            ax.plot([first,last],[row,row],'.-',lw=4,color=color)
            for outside,inside in [(e['pre'],e['first']),(e['post'],e['last'])]:
                if outside['valid']:
                    ax.plot([outside['time'],inside['time']],[row,row],':',color=color)
            if not e['pre']['valid']:
                ax.annotate('entry unknown',(first,row),xytext=(0,8),textcoords='offset points',fontsize=7)
    ax.set_yticks(range(len(labels)),labels,fontsize=7)
    if not labels:
        ax.text(.1,.5,'No encounters for this node' if selected['tracking_available'] else
                'Shock history not present',transform=ax.transAxes)
    ax.set(xlabel='simulation time [code units]',title='Inside evidence and dotted timing brackets')
    fig.suptitle(f'Curve {curve_id} · node {node_id} · {field}',fontsize=14)
    fig.supxlabel('Independent output clocks · markers are saved samples · no inferred crossing time',fontsize=9)
    output.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(output,dpi=160)
    plt.close(fig)
    return dict(curve_id=curve_id,node_id=node_id,field=field,geometry_clocks=len(frames),
                node_samples=len(history),first_sample_time=history[0]['time'] if history else None,
                node_birth=history[0]['birth'] if history else None,
                emitted_node_events=sum(e['node_id']==node_id for e in selected['events']['curve_events']),
                terminal_node_views=sum(e['node_id']==node_id for e in selected['events']['curve_suspended']),
                geometry_3d=three_d,output=str(output.resolve()))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory',type=Path,help='native run directory')
    parser.add_argument('--curve-id',type=int,help='exact integer curve ID; required for multiple curves')
    parser.add_argument('--node-id',type=int,help='exact uint64 node ID; also identifies its curve')
    parser.add_argument('--field',default='pressure',help='canonical saved thermo field name')
    parser.add_argument('--thermo-stem',help='basename.output_id if several thermo streams exist')
    parser.add_argument('--shock-directory',default='shock_events')
    parser.add_argument('-o','--output',type=Path,required=True,help='new figure filename')
    args=parser.parse_args()
    try:
        # ponytail: loads curve histories into RAM; stream by cycle if histories outgrow memory.
        data=read_run(args.directory,args.thermo_stem,args.shock_directory,include_particles=False)
        selected=select_curve(data,args.curve_id,args.node_id,args.field)
        print(json.dumps(plot_curve(selected,args.output),indent=2))
    except (ValueError,OSError) as exc:
        parser.exit(1,str(exc)+'\n')


if __name__=='__main__':main()
