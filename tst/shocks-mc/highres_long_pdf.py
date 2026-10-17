"""Adjacent captions for the larger native forced-box run; reuse vector PDF tools."""
import argparse,json,sys
from pathlib import Path
from tempfile import TemporaryDirectory
from pypdf import PdfReader,PdfWriter,Transformation
from reportlab.pdfgen import canvas
p=argparse.ArgumentParser();p.add_argument('packet',type=Path);p.add_argument('--root',type=Path,required=True);a=p.parse_args();root=a.packet.resolve();sys.path.insert(0,str(a.root/'tst/shocks-mc'))
from explained_pdf import paragraph,table,TITLE,SUB,BODY,SMALL,GRAY
r=json.loads((root/'report.json').read_text());assert r['status']=='completed_scientific_review_pending';writer=PdfWriter()
def text(c,x,y,value,style=BODY,width=530):
 item=paragraph(value,style);w,h=item.wrap(width,700);item.drawOn(c,x,y-h);return y-h-9
def frame(c,num,title):
 c.setFillColor(GRAY);c.setFont('Helvetica',8);c.drawString(36,766,'ATHENAK | 128-CUBED LONG FORCED BOX | 14 OCTOBER 2026');text(c,36,741,title,TITLE)
 c.setFillColor(GRAY);c.setFont('Helvetica',8);c.drawString(36,24,'Local CPU numerical evidence | Scientific approval pending');c.drawRightString(576,24,str(num))
with TemporaryDirectory(prefix='highres-pdf-') as tmp:
 temp=Path(tmp);cover=temp/'cover.pdf';c=canvas.Canvas(str(cover),pagesize=(612,792));frame(c,1,'Higher resolution, longer particle shock histories');y=654
 evo=r['field_evolution'];last=evo[-1];oc=r['outcomes']
 rows=[['Measurement','Actual recorded result'],['Grid / particles','128^3 / 8,192'],['End time / generations',f'{r["end_time"]:.6g} / {r["through_generation"]}'],['Complete / partial / open',f'{oc.get("complete",0):,} / {oc.get("partial",0):,} / {oc.get("open",0):,}'],['Complete episodes / particle',f'{r["complete_per_particle"]:.4f}'],['Final volume rms Mach',f'{last["mach_volume"]:.4f}'],['Final supersonic volume fraction',f'{100*last["supersonic_volume_fraction"]:.2f}%'],['Validated per-step observations',f'{r["counters"]["sampled"]:,}'],['Native run wall time',f'{r["execution"]["wall_seconds"]/60:.1f} minutes']]
 item=table(rows,[235,295]);w,h=item.wrap(530,700);item.drawOn(c,36,y-h);y-=h+18
 y=text(c,36,y,'Experiment',SUB);y=text(c,36,y,'Periodic unit cube, ideal Newtonian gas gamma=1.4, native RK2/PLM/HLLC and CFL0.30. Initial velocity rms Mach5, continuous native edot50 forcing, the same seeds and8,192 mass-flux MC tracers as the earlier64-cubed run. Detector width8 and8ghost cells preserve the physical endpoint half-span1/16box. This run evolves to t=1, compared with the earlier t=0.3. No production detector, particle or forcing machinery changed.')
 y=text(c,36,y,'What passed',SUB);y=text(c,36,y,'Exact frozen ShockIt mask/strength parity for new width8 fixtures, including an independent planar Mach3 jump oracle. Full per-step traces retain all particle IDs, raw shock flags, phases and precise observation clocks. Canonical event sequences are unique and contiguous; native counters agree. Positive finite checkpoint fluid states conserve periodic mass to1e-12. Same-binary zero-step checkpoint recovery matches every particle thermodynamic sample; plotted native detector slices match live particle mask/strength evidence to the unchanged1e-12 tolerance.')
 y=text(c,36,y,'Limits and reproducibility',SUB);y=text(c,36,y,'The inherited RK2 forcing thermal oracle remains a documented failure. Adiabatic heating can reduce the flow Mach despite continuous forcing. These are sampled numerical encounters, not certified physical shock counts or exact crossing times. Stationarity, ensembles, convergence, RUSTY/Frontier performance and human scientific approval remain open. The README records exact commands, source/executable/input hashes, output-cadence recovery, full trace validation and raw-data paths.',SMALL)
 assert y>42,y;c.save();writer.add_page(PdfReader(cover).pages[0])
 for num,f in enumerate(r['figures'],2):
  pagefile=temp/(f['name']+'.pdf');c=canvas.Canvas(str(pagefile),pagesize=(612,792));frame(c,num,f'{num-1:02d}. '+f['title'])
  plot=PdfReader(root/(f['name']+'.pdf')).pages[0];scale=min(540/float(plot.mediabox.width),390/float(plot.mediabox.height));height=float(plot.mediabox.height)*scale;bottom=654-height
  y=text(c,36,bottom-18,f['caption']);y=text(c,36,y,'Data: '+f['name']+'.pdf / .svg / .png; episodes.csv, evolution.csv, selected-traces.csv and report.json.',SMALL);assert y>42,(f['name'],y)
  c.save();page=PdfReader(pagefile).pages[0];page.merge_transformed_page(plot,Transformation().scale(scale).translate(36,bottom));writer.add_page(page)
 writer.add_metadata({'/Title':'AthenaK:128-cubed longer forced turbulence and particle shock histories','/Subject':'Exact numerical evidence; scientific approval pending'})
 with (root/'highres-long-report.pdf').open('wb') as f:writer.write(f)
print(json.dumps({'pdf':str(root/'highres-long-report.pdf'),'pages':len(writer.pages),'vector_plots_preserved':True}))
