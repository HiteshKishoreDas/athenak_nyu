"""Compact selected raw MC histories and analytic front overlays, no time estimator."""
import argparse
import csv
import json
import math
from pathlib import Path
import shutil

from common import ROOT, sha
from integration import traces
from shock_history import read_history, read_suspended


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--runs',type=Path,required=True)
    parser.add_argument('--report',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--plot',action='store_true')
    args=parser.parse_args();output=args.output.resolve();output.mkdir(parents=True,exist_ok=True)
    report=json.loads(args.report.read_text())
    selected=[]
    for label in ('mach1.2','mach1.5','mach5','boost2'):
        source=args.runs.resolve()/label/'shock_events';events=read_history(source)
        suspended=read_suspended(source);observed=traces(source)
        reentries={e['particle_id'] for e in events if e['sequence']>0}
        ids=sorted(reentries|{e['particle_id'] for e in suspended}|{0,1,2})
        speed=2.5 if label=='boost2' else .5
        rows=[]
        for (g,p),o in sorted(observed.items()):
            if p not in ids:continue
            front=.3+speed*o['time']
            rows.append([g,o['time'],p,o['phase'],int(o['valid']),int(o['inside']),
                *o['x'],*o['image'],o['rho'],o['pressure'],o['mach_pressure'],
                o['mach_temperature'],*o['normal'],front,o['x'][0]-front,o['sequence']])
        path=output/(label+'_selected.csv')
        with path.open('w') as file:
            writer=csv.writer(file);writer.writerow(['generation','time','particle_id','phase','valid','inside',
                'x','y','z','image_x','image_y','image_z','rho','pressure','mach_pressure','mach_temperature',
                'normal_x','normal_y','normal_z','analytic_front_x','x_minus_analytic_front','sequence'])
            writer.writerows(rows)
        transfer_path=output/(label+'_selected_transfers.csv')
        transfer_records=0
        with transfer_path.open('w') as file:
            writer=csv.writer(file);writer.writerow(['generation','rank','t_start','t_end','random_seed',
                'tag','source_gid','ip','jp','kp','start_mass','p_xminus','p_xplus','p_yminus','p_yplus',
                'p_zminus','p_zplus','source_x','source_y','source_z','draw','direction'])
            for source_transfer in sorted(source.glob('transfer*.rank*.csv')):
                lines=source_transfer.read_text().splitlines();head=lines[0].split()
                rank=int(source_transfer.name.split('.rank')[1].split('.')[0])
                for values in csv.reader(lines[1:-1]):
                    if int(float(values[0])) in ids:
                        writer.writerow([head[4],rank,*head[6:9],*values]);transfer_records+=1
        records=[e for e in events if e['particle_id'] in ids]
        item={'label':label,'selected_ids':ids,'original_raw_directory':str(source),
              'compiled_source_fingerprint':report.get('compiled_source_fingerprint'),
              'input_sha256':sha(source.parent/'input.athinput'),
              'selected_csv_sha256':sha(path),'rows':len(rows),
              'transfer_csv_sha256':sha(transfer_path),'transfer_rows':transfer_records,
              'events':records,'suspended':suspended,
              'scope':'cell-center MC samples and known Eulerian planar front; no continuous material trajectory or exact crossing time'}
        (output/(label+'_records.json')).write_text(json.dumps(item,indent=2,allow_nan=True)+'\n')
        target=output/'nonempty_chunks'/label;target.mkdir(parents=True,exist_ok=True)
        digests={}
        for chunk in source.glob('cycle*.rank0.csv'):
            digests[chunk.name]=sha(chunk)
            if int(chunk.read_text().splitlines()[0].split()[-1]):shutil.copy2(chunk,target/chunk.name)
        for chunk in source.glob('suspended*.rank0.csv'):shutil.copy2(chunk,target/chunk.name)
        shutil.copy2(source/'manifest.rank0',target/'manifest.rank0')
        (target/'all_committed_chunk_sha256.json').write_text(json.dumps(digests,indent=2)+'\n')
        selected.append(item)
        if args.plot:
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
            fig,(ax,dist)=plt.subplots(2,1,figsize=(8,7),sharex=True)
            for pid in ids:
                history=[r for r in rows if r[2]==pid]
                line=ax.step([r[1] for r in history],[r[6] for r in history],where='post',
                             linewidth=.9,label='ID '+str(pid))[0]
                hits=[r for r in history if r[5]]
                ax.scatter([r[1] for r in hits],[r[6] for r in hits],s=7,color=line.get_color())
                dist.step([r[1] for r in history],[r[20] for r in history],where='post',linewidth=.9)
            end=max(r[1] for r in rows)
            ax.plot([0,end],[.3,.3+speed*end],'k--',label='analytic Eulerian front')
            dist.axhline(0,color='black',ls='--',lw=1)
            ax.set_ylabel('sampled cell-center x');dist.set_ylabel('x − analytic front');dist.set_xlabel('code time')
            ax.set_title(label+': MC jumps; dots show unreduced shock-region membership')
            ax.legend(fontsize=8,ncol=4);fig.tight_layout()
            fig.savefig(output/(label+'_front.png'),dpi=180);fig.savefig(output/(label+'_front.pdf'))
            plt.close(fig)
    (output/'manifest.json').write_text(json.dumps({'source_report':str(args.report.resolve()),
        'source_report_sha256':sha(args.report),'selected_cases':selected},indent=2,allow_nan=True)+'\n')


if __name__=='__main__':main()
