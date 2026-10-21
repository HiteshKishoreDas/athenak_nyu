"""Synchronized material-curve + shock showcase from native primitive snapshots.

Uses existing native readers and the exact supplied-state C++ detector. Saved
node/particle flags remain authoritative and must match every reconstructed cell.
No new physical setup, interpolated movie frames, or synthetic particle paths.
"""
import argparse
from collections import Counter
import csv
import json
from pathlib import Path
import subprocess
import sys
import textwrap

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D
import numpy as np
from PIL import Image

ROOT=Path(__file__).resolve().parents[2]
sys.path[:0]=[str(ROOT/'vis/python'),str(ROOT/'tst/shocks-mc'),str(ROOT/'tst/shocks-curves')]
from common import sha
from acceptance import snapshots, trace_rows
from bin_convert import read_binary
from shock_history import parse_observation, read_history, read_suspended
import curve_shock_history as curve_history

BLUE,ORANGE,CYAN,RED,PURPLE='#2166ac','#e47b13','#059ba7','#c03259','#8172b3'
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})


def parameters(path):
    result={};block=None
    for line in Path(path).read_text().splitlines():
        line=line.split('#',1)[0].strip()
        if line.startswith('<') and line.endswith('>'):block=line[1:-1]
        elif '=' in line and block:
            key,value=line.split('=',1);result[block+'/'+key.strip()]=value.strip()
    return result


def particle_frame(path,generation):
    result={}
    for row in csv.reader(Path(path).read_text().splitlines()):
        assert len(row)==28,'native trace schema changed'
        p=int(row[0]);o=parse_observation(row[4:])
        assert p not in result and o['generation']==generation
        result[p]=o
    return result


def detector_frame(data,config,executable,directory):
    assert data['n_mbs']==1 and data['Nx1']==data['Nx2']==data['Nx3']==32
    assert data['var_names']==['dens','velx','vely','velz','eint']
    fields=[np.asarray(data['mb_data'][k])[0] for k in data['var_names']]
    assert all(f.dtype==np.float64 and np.isfinite(f).all() for f in fields)
    assert np.all(fields[0]>0) and np.all(fields[-1]>0)
    nx,ny,nz=data['Nx1'],data['Nx2'],data['Nx3']
    prefix=directory/f"cycle{data['cycle']}";fixture=prefix.with_suffix('.txt')
    values=[float(config.get(key,default)) for key,default in [
        ('hydro/gamma',1.4),('shock_tracking/min_mach',1.1),
        ('shock_tracking/mach_max',1000),('shock_tracking/grad_floor',1e-30),
        ('shock_tracking/div_limit',0)]]
    width=int(config['shock_tracking/width'])
    values += [width,float(config.get('shock_tracking/density_floor',0)),
               float(config.get('shock_tracking/temperature_floor',0)),1/nx,1/ny,1/nz]
    with fixture.open('w') as f:
        f.write(f'{nx} {ny} {nz}\n')
        f.write(' '.join(str(v) if i==5 else format(v,'.17g') for i,v in enumerate(values))+'\n')
        np.savetxt(f,np.stack(fields,axis=-1).transpose(2,1,0,3).reshape(-1,5),fmt='%.17g')
    path=prefix.with_suffix('.csv');command=[str(executable),str(fixture)]
    with path.open('w') as f:
        result=subprocess.run(command,stdout=f,stderr=subprocess.PIPE,text=True,timeout=30)
    assert result.returncode==0,result.stderr
    numeric=np.loadtxt(path,delimiter=',',comments='#')
    assert numeric.shape==(nx*ny*nz,30)
    expected=np.indices((nx,ny,nz)).reshape(3,-1).T
    assert np.array_equal(numeric[:,:3],expected),'detector coordinate order changed'
    valid=numeric[:,3].reshape(nx,ny,nz).astype(bool)
    inside=numeric[:,5].reshape(nx,ny,nz).astype(bool)
    np.savez_compressed(prefix.with_suffix('.npz'),valid=valid,inside=inside)
    return fields[0],valid,inside,dict(command=command,fixture_sha256=sha(fixture),
                                     signals_sha256=sha(path),flagged_cells=int(inside.sum()))


def wrapped_line(points):
    wrapped=np.mod(points,1);result=[]
    for i,x in enumerate(wrapped):
        if i and np.any(np.abs(x-wrapped[i-1])>.5):result.append([np.nan]*3)
        result.append(x)
    return np.asarray(result)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run',type=Path,required=True)
    parser.add_argument('--detector',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();run=args.run.resolve();out=args.output.resolve()
    out.mkdir(parents=True,exist_ok=False);frames_dir=out/'frames';frames_dir.mkdir()
    diagnostics=out/'detector';diagnostics.mkdir()
    execution=json.loads((run/'execution.json').read_text());assert execution['status']=='pass'
    assert sha(run/'input.athinput')==execution['input_sha256']
    assert sha(run/'run.log')==execution['log_sha256']
    config=parameters(run/'input.athinput')
    assert config['particles/particle_type']==config['particles/pusher']=='lagrangian_tracer'
    assert config['shock_tracking/trace']=='true'
    for d in range(1,4):
        assert config[f'mesh/ix{d}_bc']==config[f'mesh/ox{d}_bc']=='periodic'
    states=snapshots(run);observed=trace_rows(run/'shock_events/curves',curve_history)
    ids=sorted({p for g,p in observed});generations=sorted(states)
    assert len(ids)==16 and all(set(p for q,p in observed if q==g)==set(ids) for g in generations)
    clocks={g:states[g]['time'] for g in generations}
    events=curve_history.read_history(run/'shock_events/curves')
    particle_events=read_history(run/'shock_events');terminal=read_suspended(run/'shock_events')
    assert len(events)==28 and len(particle_events)==6904
    assert all(e['status']==0 for e in events+particle_events)
    matrix=np.array([[(3 if observed[g,p]['phase']==0 else 1) if observed[g,p]['inside']
                          else 0 if observed[g,p]['valid'] else 2
                      for g in generations] for p in ids])
    initial_length=states[0]['curves'][0]['length'];frames=[];bindings=[]
    counts=[Counter(e['emission_generation'] for e in records) for records in [events,particle_events]]
    cumulative=[np.cumsum([c[g] for g in generations]) for c in counts]
    for field_file in sorted((run/'bin').glob('*.bin')):
        data=read_binary(field_file);g=data['cycle'];state=states[g];time=clocks[g]
        assert abs(data['time']-time)<=5e-7,'rounded binary clock inconsistent with native cycle'
        rho,valid,inside,binding=detector_frame(data,config,args.detector.resolve(),diagnostics)
        node_rows=[observed[g,p] for p in ids]
        particle_path=run/f'shock_events/trace{g}.rank0.csv'
        particles=particle_frame(particle_path,g);assert len(particles)==4096
        matches=0
        for o in node_rows+list(particles.values()):
            assert o['time']==time
            index=tuple(np.floor(np.mod(o['x'],1)*32).astype(int))
            assert bool(valid[index])==o['valid'],'cell/record validity mismatch'
            assert bool(inside[index])==o['inside'],'cell/record membership mismatch'
            matches+=1
        points=np.array([n['x'] for n in state['curves'][0]['nodes']])
        assert [n['id'] for n in state['curves'][0]['nodes']]==ids
        assert all(np.max(np.abs(np.mod(n['x'],1)-o['x']))<=2e-12
                   for n,o in zip(state['curves'][0]['nodes'],node_rows))
        fraction=inside.mean(axis=2).T
        frames.append(dict(cycle=g,time=time,column=rho.mean(axis=0),fraction=fraction,
            points=points,node_rows=node_rows,particles=particles,
            inside_nodes=sum(o['inside'] for o in node_rows),
            length_ratio=state['curves'][0]['length']/initial_length,
            node_events=sum(e['emission_generation']<=g for e in events),
            particle_events=sum(e['emission_generation']<=g for e in particle_events)))
        binding.update(cycle=g,time=time,native_field=str(field_file),field_sha256=sha(field_file),
            particle_trace_sha256=sha(particle_path),membership_matches=matches)
        bindings.append(binding)
        print(f'cycle {g}: {matches} exact mask/flag matches; {frames[-1]["inside_nodes"]}/16 nodes inside',flush=True)
    vmin=min(np.log10(f['column']).min() for f in frames)
    vmax=max(np.log10(f['column']).max() for f in frames)
    centers=(np.arange(32)+.5)/32
    legend=[Line2D([],[],color=CYAN,lw=2,label='Material curve'),
            Line2D([],[],marker='o',ls='',mfc=ORANGE,mec='white',label='Inside; episode active'),
            Line2D([],[],marker='o',ls='',mfc=PURPLE,mec='white',label='Inside; episode suppressed'),
            Line2D([],[],marker='o',ls='',mfc=BLUE,mec='white',label='Node outside shock')]

    def spatial(ax,frame,zoom=False):
        column=np.log10(frame['column']);fraction=frame['fraction'];points=frame['points']
        if zoom:
            image=ax.imshow(np.tile(column,(3,3)),origin='lower',extent=(-1,2,-1,2),
                            cmap='Greys',vmin=vmin,vmax=vmax,interpolation='nearest')
            grid=np.concatenate([centers-1,centers,centers+1]);field=np.tile(fraction,(3,3))
            line=points;center=points[:,:2].mean(axis=0)
            span=max(.14,float(np.ptp(points[:,:2],axis=0).max())*1.45)
            ax.set_xlim(center[0]-span/2,center[0]+span/2)
            ax.set_ylim(center[1]-span/2,center[1]+span/2)
        else:
            image=ax.imshow(column,origin='lower',extent=(0,1,0,1),cmap='Greys',
                            vmin=vmin,vmax=vmax,interpolation='nearest')
            grid=centers;field=fraction;line=wrapped_line(points)
            sample=[o for p,o in sorted(frame['particles'].items()) if p%16==0]
            xy=np.array([o['x'] for o in sample]);inside_=[o['inside'] for o in sample]
            ax.scatter(xy[:,0],xy[:,1],c=[ORANGE if flag else BLUE for flag in inside_],
                       s=7,alpha=.3,edgecolors='none')
            ax.set_xlim(0,1);ax.set_ylim(0,1)
        levels=[l for l in (.1,.3,.5) if field.min()<l<field.max()]
        if levels:ax.contour(grid,grid,field,levels=levels,colors=ORANGE,linewidths=.7,alpha=.6)
        ax.plot(line[:,0],line[:,1],color='white',lw=4,zorder=4)
        ax.plot(line[:,0],line[:,1],color=CYAN,lw=2,zorder=5)
        nodes=points if zoom else np.mod(points,1)
        colors=[(PURPLE if o['phase']==0 else ORANGE) if o['inside']
                else (BLUE if o['valid'] else RED) for o in frame['node_rows']]
        ax.scatter(nodes[:,0],nodes[:,1],c=colors,s=62 if zoom else 25,
                   edgecolors='white',linewidths=.8,zorder=6)
        if zoom:
            for p,(x,y,z) in enumerate(points):
                if p in (0,3,7,11,15):ax.annotate(str(p+1),(x,y),xytext=(4,5),
                    textcoords='offset points',fontsize=8,color='#152b40',zorder=7)
        ax.set_aspect('equal');ax.set(xlabel=r'$x$ [code units]',ylabel=r'$y$ [code units]')
        return image

    caption=('Same native forced 32³ run: 4096 continuous particles + 16 material-curve nodes. '
             'Background = density integrated along z; orange contours = fractions 0.1/0.3/0.5 of z-cells flagged as shocks. '
             'Projected features may lie at different depths. Orange/purple nodes have actual 3D inside flags; purple means the existing episode filter suppressed activation. '
             'All fields, geometry and flags share the displayed native cycle. Nine stored frames; no temporal interpolation. '
             'Short engineering demonstration; turbulence convergence and scientific approval remain pending.')
    image_paths=[]
    for frame in frames:
        fig,ax=plt.subplots(2,2,figsize=(13.5,9),gridspec_kw={'height_ratios':[1.65,1]})
        fig.suptitle('Material curves meet shocks — both features in one native run',
                     fontsize=17,fontweight='bold',x=.055,ha='left')
        im=spatial(ax[0,0],frame)
        fig.colorbar(im,ax=ax[0,0],label=r'$\log_{10}\Sigma$ [code units]',fraction=.045,pad=.02)
        ax[0,0].set_title(f"Full box, x–y projection | t={frame['time']:.6f}, cycle {frame['cycle']}")
        spatial(ax[0,1],frame,True)
        ax[0,1].set_title(f"Curve zoom | {frame['inside_nodes']}/16 raw inside; "
                         f"{frame['node_events']} complete passages")
        ax[0,1].legend(handles=legend,fontsize=8,loc='upper left')
        ax[0,1].text(.02,.02,f"3D length / initial length = {frame['length_ratio']:.2f}",
            transform=ax[0,1].transAxes,bbox=dict(facecolor='white',alpha=.85,edgecolor='none'),fontsize=9)
        ax[1,0].pcolormesh([clocks[g] for g in generations],range(1,17),matrix,shading='nearest',
            cmap=ListedColormap(['#e1eaf2',ORANGE,RED,PURPLE]),vmin=0,vmax=3)
        ax[1,0].axvline(frame['time'],color='#1b273b',lw=1.2)
        ax[1,0].set(xlabel=r'$t$ [code units]',ylabel='Curve node order',
                    title='Node history | orange = active inside; purple = suppressed inside')
        for total,n,label,color in zip(cumulative,[16,4096],['Curve nodes (16)','Particles (4096)'],[CYAN,BLUE]):
            ax[1,1].step([clocks[g] for g in generations],total/n,where='post',color=color,label=label)
            index=generations.index(frame['cycle'])
            ax[1,1].scatter(frame['time'],total[index]/n,color=color,s=30,zorder=5)
        ax[1,1].axvline(frame['time'],color='#1b273b',lw=1)
        ax[1,1].set(xlabel=r'$t$ [code units]',ylabel='Complete passages per initial carrier',
                    title=f"Committed records: {frame['node_events']} curve-node / {frame['particle_events']} particle")
        ax[1,1].legend(fontsize=9);ax[1,1].grid(alpha=.15)
        fig.text(.045,.015,textwrap.fill(caption,175),fontsize=9,va='bottom')
        fig.tight_layout(rect=(.015,.12,.995,.94),w_pad=3,h_pad=1.5)
        path=frames_dir/f"frame_{len(image_paths):03}.png";fig.savefig(path,dpi=130);plt.close(fig)
        image_paths.append(path)
    selected=max(range(len(frames)),key=lambda i:(0<frames[i]['inside_nodes']<16 and frames[i]['node_events']>0,
                                                       frames[i]['inside_nodes']))
    (out/'showcase.png').write_bytes(image_paths[selected].read_bytes())
    images=[Image.open(p).convert('RGB') for p in image_paths]
    images[0].save(out/'showcase.gif',save_all=True,append_images=images[1:],
                   duration=[1100]+[800]*(len(images)-2)+[1800],loop=0)
    for im in images:im.close()
    subprocess.run(['ffmpeg','-y','-framerate','1.5','-i',str(frames_dir/'frame_%03d.png'),
                    '-c:v','libx264','-pix_fmt','yuv420p','-vf','pad=ceil(iw/2)*2:ceil(ih/2)*2',
                    str(out/'showcase.mp4')],stdout=subprocess.DEVNULL,stderr=(out/'ffmpeg.log').open('w'),check=True)
    fig,axes=plt.subplots(2,2,figsize=(11,9))
    for ax,frame in zip(axes.flat,[frames[i] for i in [1,3,4,8]]):
        spatial(ax,frame,True);ax.set_title(f"t={frame['time']:.4f} | {frame['inside_nodes']}/16 inside | "
            f"{frame['node_events']} passages | 3D L/L₀={frame['length_ratio']:.2f}",fontsize=10)
    fig.legend(handles=legend,loc='upper center',ncol=3,bbox_to_anchor=(.5,.95),fontsize=10)
    fig.suptitle('The same material curve deforms and retains shock encounters',fontsize=15,fontweight='bold')
    fig.text(.04,.015,textwrap.fill(caption,155),fontsize=9)
    fig.tight_layout(rect=(.02,.15,.995,.91));fig.savefig(out/'curve_evolution.png',dpi=150);plt.close(fig)
    fig,(ax,bars)=plt.subplots(1,2,figsize=(13.5,6),gridspec_kw={'width_ratios':[3,1]})
    strength=[e['representative']['mach_pressure'] for e in events]
    assert all(np.isfinite(strength))
    norm=plt.Normalize(min(strength),max(strength));color_map=plt.get_cmap('viridis')
    node_index={p:i+1 for i,p in enumerate(ids)}
    for e in events:
        y=node_index[e['node_id']]
        for a,b in [('pre','first'),('last','post')]:
            assert e[a]['valid'] and e[b]['valid']
            ax.plot([e[a]['time'],e[b]['time']],[y,y],color='#97a8b8',lw=2)
        ax.plot([e['first']['time'],e['last']['time']],[y,y],
                color=color_map(norm(e['representative']['mach_pressure'])),lw=7,solid_capstyle='butt')
        if e['first']['time']==e['last']['time']:
            ax.scatter(e['first']['time'],y,s=20,marker='s',
                       color=color_map(norm(e['representative']['mach_pressure'])),zorder=5)
        ax.text((e['first']['time']+e['last']['time'])/2,y+.12,str(e['sequence']+1),
                ha='center',fontsize=8,color='#173047')
    ax.set(xlabel=r'$t$ [code units]',ylabel='Curve node order',yticks=range(1,17),
           title='Per-node episodes: color = strength evidence; numeral = episode number',ylim=(.5,16.7))
    fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=color_map),ax=ax,
                 label='Representative pressure-jump Mach proxy',pad=.02,fraction=.04)
    counts=Counter(e['node_id'] for e in events)
    bars.barh(range(1,17),[counts[p] for p in ids],color=CYAN)
    bars.set(xlabel='Complete passages',ylabel='Curve node order',yticks=range(1,17),
             title='28 records retained after exit',ylim=(.5,16.7))
    fig.suptitle('A material curve carries separate shock histories on its nodes',fontsize=16,fontweight='bold')
    fig.text(.04,.02,textwrap.fill('Thick segments span first-to-last accepted inside observations; squares mark single-sample episodes. '
        'Raw inside flags can be suppressed by the inherited upstream-order filter. Thin grey extensions are '
        'entry/exit sampling brackets, not exact crossings. Colors are saved strength estimates, not calibrated plateau '
        'truth. Counts are per node; they do not identify a unique shock surface along the whole curve.',165),fontsize=9)
    fig.tight_layout(rect=(.01,.16,.995,.92));fig.savefig(out/'node_passages.png',dpi=150);plt.close(fig)

    metadata=dict(status='pass',native_run=str(run),simulation_source_commit='1c325e5a59d64f743ffdaa4f8678a4805b12608f',
        input_sha256=sha(run/'input.athinput'),detector_executable_sha256=sha(args.detector),
        script_sha256=sha(__file__),frames=bindings,matched_memberships=sum(b['membership_matches'] for b in bindings),
        node_events=len(events),particle_events=len(particle_events),particle_terminal_censored=len(terminal),
        final_length_ratio=frames[-1]['length_ratio'],poster_cycle=frames[selected]['cycle'],
        raw_inside_node_samples=sum(o['inside'] for o in observed.values()),
        suppressed_inside_node_samples=sum(o['inside'] and o['phase']==0 for o in observed.values()),
        direction_ambiguous_node_events=sum(bool(e['flags']&256) for e in events),
        scientific_approval='pending',physics_changes=False,
        limitations=['3D projection mixes depths','9 stored fluid times; no interpolated movie frames',
                     'Short uniform-mesh CPU engineering run; no converged-turbulence/scientific approval'])
    (out/'report.json').write_text(json.dumps(metadata,indent=2)+'\n')
    labels=[f"t={f['time']:.6f}, cycle {f['cycle']}" for f in frames]
    (out/'index.html').write_text('''<!doctype html><meta charset="utf-8"><title>Material curves + shock histories</title>
<style>body{font:17px system-ui;margin:25px auto;max-width:1500px;padding:0 20px;color:#163047}img,video{width:100%}input{width:75%}button{font-size:17px;padding:8px}p{line-height:1.5}.note{background:#fff0dd;padding:16px}</style>
<h1>Material curves + shock histories, together</h1><p>Native forced 32³ box: 4096 continuous velocity-following particles and 16 curve nodes. The cyan curve moves with the flow. Orange/purple nodes are inside the numerical shock region; purple marks flags suppressed by the existing episode filter.</p>
<p class="note">Saved native cycles, synchronized by cycle ID. Background and contours project a 3D volume along z, so apparent overlaps can have different depths. Individual node colors use actual 3D samples. Scientific approval remains pending.</p>
<button id="play">Pause</button> <input id="slide" type="range" min="0" max="8" value="0"><span id="clock"></span>
<img id="frame" src="frames/frame_000.png"><p>Top: fluid and detected shocks with the material curve, then a close view. Bottom: each node’s saved flags and completed histories per initial carrier. Curves and particles have separate histories, sampled from the same fluid/detector.</p>
<h2>Curve deformation and shock passages</h2><img src="curve_evolution.png"><p>Four actual saved states. L/L₀ uses full 3D curve length. Nodes retain their IDs and event sequences; counted records are completed episodes, not orange samples.</p>
<h2>Histories retained on each curve node</h2><img src="node_passages.png"><p>Episode numbers follow each node’s saved sequence. Thin extensions show timing brackets. These are node events, not a count of unique shock surfaces. Of 775 raw inside node samples, 548 are suppressed by the inherited upstream-order filter; one completed event has a direction-ambiguity flag. These classifications are not scientifically approved.</p><p><a href="showcase.gif">Animated GIF</a> · <a href="showcase.mp4">MP4</a> · <a href="report.json">Source and mask/flag checks</a></p>
<script>const labels='''+json.dumps(labels)+''';let playing=true;const s=document.getElementById('slide');function show(i){s.value=i;document.getElementById('frame').src='frames/frame_'+String(i).padStart(3,'0')+'.png';document.getElementById('clock').textContent=labels[i];}s.oninput=()=>{playing=false;document.getElementById('play').textContent='Play';show(+s.value)};document.getElementById('play').onclick=()=>{playing=!playing;document.getElementById('play').textContent=playing?'Pause':'Play'};setInterval(()=>{if(playing)show((+s.value+1)%9)},1000);show(0);</script>''')
    print(json.dumps({k:v for k,v in metadata.items() if k!='frames'},indent=2))


if __name__=='__main__':main()
