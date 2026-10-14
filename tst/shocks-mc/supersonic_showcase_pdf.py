"""Compose vector showcase figures with adjacent explanations; bundled PDF Python."""
import argparse
import json
from pathlib import Path
from tempfile import TemporaryDirectory

from pypdf import PdfReader,PdfWriter,Transformation
from reportlab.pdfgen import canvas
from explained_pdf import paragraph,table,TITLE,SUB,BODY,SMALL,BLUE,GRAY


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('packet',type=Path);a=p.parse_args();root=a.packet.resolve();r=json.loads((root/'report.json').read_text());writer=PdfWriter()
    def text(c,x,y,value,style=BODY,width=530):
        item=paragraph(value,style);w,h=item.wrap(width,700);item.drawOn(c,x,y-h);return y-h-9
    def frame(c,num,title):
        c.setFillColor(GRAY);c.setFont('Helvetica',8);c.drawString(36,766,'ATHENAK | FORCED SUPERSONIC SHOWCASE | 13 OCTOBER 2026');text(c,36,741,title,TITLE)
        c.setFillColor(GRAY);c.setFont('Helvetica',8);c.drawString(36,24,'Evolving adiabatic flow | Physical strength calibration and scientific review pending');c.drawRightString(576,24,str(num))
    with TemporaryDirectory(prefix='showcase-pdf-') as tmp:
        temp=Path(tmp);cover=temp/'cover.pdf';c=canvas.Canvas(str(cover),pagesize=(612,792));frame(c,1,'Shock encounters in a forced supersonic box')
        y=text(c,36,656,'Native fluid-following Monte Carlo particles retain repeated numerical shock-region encounters, timing brackets, strength evidence and censored outcomes. This report showcases actual local simulation data; it is not a certified physical shock catalogue.')
        rows=[['Case','Complete','Open','Complete / particle','Final rms Mach']]
        for name,s in r['cases'].items():rows.append([name,str(s['completed']),str(s['open']),f'{s["emitted_per_particle"]:.3f}',f'{s["field_evolution"][-1]["mach_volume"]:.3f}'])
        item=table(rows,[125,70,70,130,135]);w,h=item.wrap(530,700);item.drawOn(c,36,y-h);y-=h+24
        y=text(c,36,y,'Experiment',SUB);y=text(c,36,y,'Periodic unit cube,64³ cells,8,192 mass-seeded native MC particles; uniform initial density 1 and pressure 1/1.4; initial Fourier velocity rms Mach 5. Newtonian ideal gas gamma=1.4, RK2/PLM/HLLC/CFL.30; continuous native edot=50 forcing, modes 1..3, npeak=2, solenoidal fraction 0.5, OU correlation=0.1, update=0.01. End time 0.30 equals three initial driving-scale turnover times (Ldrive=0.5, initial vrms=5); the later flow evolves and this is not a stationarity assertion.')
        y=text(c,36,y,'What is checked',SUB);y=text(c,36,y,'Five simulations succeeded. Tracking-on/off entire particle thermo bytes and every saved primitive/conserved/force array agree exactly. All particle IDs survive; event IDs are unique and sequences contiguous; tracker counters agree with strict event readers. Periodic mass error is at most1e-12. Three saved volumes are evaluated through the unchanged native detector fixture with genuine periodic neighbors. Binary times are rounded; plotting and sampling align by generation and full-double native particle clocks.')
        y=text(c,36,y,'Paper-use boundary',SUB);y=text(c,36,y,'Use the separate vector PDF/SVG figures and CSV data as a reproducible methods/implementation demonstration. Native RK2 forcing has a separately measured gas-heating bias. There is no cooling, stationary ensemble, accepted strength error budget, physical passage truth catalogue, independent-seed uncertainty estimate or cluster qualification. Captions state these limits next to each figure. All flags and open episodes remain in the data.')
        assert y>38;c.save();writer.add_page(PdfReader(cover).pages[0])
        for i,f in enumerate(r['figures'],2):
            pagefile=temp/(f['name']+'.pdf');c=canvas.Canvas(str(pagefile),pagesize=(612,792));frame(c,i,f'{i-1:02d}. '+f['title'])
            plot=PdfReader(root/(f['name']+'.pdf')).pages[0];scale=min(540/float(plot.mediabox.width),390/float(plot.mediabox.height));height=float(plot.mediabox.height)*scale
            # Caption sits directly under the actual vector figure, not in an appendix.
            bottom=654-height;y=text(c,36,bottom-18,f['caption'],BODY)
            y=text(c,36,y,'Data: '+f['name']+'.pdf / .svg / .png; episodes.csv, evolution.csv and report.json.',SMALL)
            assert y>42,(f['name'],y);c.save();page=PdfReader(pagefile).pages[0];page.merge_transformed_page(plot,Transformation().scale(scale).translate(36,bottom));writer.add_page(page)
        methods=temp/'methods.pdf';c=canvas.Canvas(str(methods),pagesize=(612,792));frame(c,12,'Methods, definitions and reproduction');y=654
        for title,body in [
            ('Event and time definitions','Membership uses the unreduced native ShockIt-compatible shock region, min_mach=1.1 and width 4 at 64³. Temperature gradients orient nearest-axis endpoints; first valid inside evidence is the representative strength. Latest valid outside pre-entry and first inside observations bracket entry; last inside and first valid post observations bracket exit. Two consecutive valid outside observations close/rearm; brief reentries merge and are flagged. Native MC motion follows density fluxes through stochastic parent-cell hops, not integrated smooth trajectories.'),
            ('Statistics and censoring','Counts include every particle and every unique canonical identity. Completed duration/span bounds are[last-first,post-pre]. Consecutive uncensored completed gaps are[max(0,next.pre-previous.post),max(0,next.first-previous.last)]. No midpoint/exact physical time is fabricated. Terminal-open views are separately retained; they do not contribute finite complete durations. All no_plateau/cadence_limited and other flags are available in the CSV/report. Rates are first-inside recording-time rates, not physical crossing rates.'),
            ('Comparison interpretation','HalfCFL changes fluid evolution, RNG draws and observation cadence.32³ width 4 doubles physical endpoint half-span;32³ width 2 matches64³ width 4 half-span. Both 32³ runs have identical native fluid/particle data, isolating their width effect. A single same-seed comparison is sensitivity evidence, not physical convergence or uncertainty. Source bias may affect turbulent temperature/strengths and is not equated with conservative shock heating.'),
            ('Exact code and environment','Feature coread4225edea78d255b67df97fe8f593c0c268fa1e; compiled fingerprint'+r['compiled_source_fingerprint']+'. Serial CPU/double macOS local run; no production source or diagnostic thresholds were changed to make the plots. Execution/report JSON bind exact executable, commands, source/configuration and logs. Raw simulation path is recorded in README. Parent forcing qualification retains its failing independent thermodynamic oracle.'),
            ('Reproduce','Run the showcase script with --executable build-release/src/athena --detector build-release/shock_detector_test --workdir /tmp/showcase-new --output /tmp/showcase-plots-new. Both directories must be new. The existing NumPy/Matplotlib environment generates10 figure families, each PDF/SVG/PNG. Run supersonic_showcase_pdf.py on that packet using bundled ReportLab/pypdf Python to generate this report. The README gives exact paths and commands. Production-scale/RUSTY/Frontier and human scientific approval remain separate gates.')]:
            y=text(c,36,y,title,SUB);y=text(c,36,y,body,SMALL)
        assert y>40;c.save();writer.add_page(PdfReader(methods).pages[0])
        writer.add_metadata({'/Title':'AthenaK: forced supersonic turbulence and particle shock encounters','/Author':'AthenaK showcase evidence','/Subject':'Numerical implementation showcase; scientific review pending'})
        with (root/'supersonic-showcase.pdf').open('wb') as f:writer.write(f)
    print(json.dumps({'pdf':str(root/'supersonic-showcase.pdf'),'pages':len(writer.pages),'figure_pages':10,'vector_plots_preserved':True}))


if __name__=='__main__':main()
