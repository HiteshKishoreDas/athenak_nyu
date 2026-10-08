"""Explain an existing laptop plot packet in a PDF, retaining vector plot pages.

Use the bundled document Python (ReportLab/pypdf), separately from the simulator.
"""
import argparse
from html import escape
import json
from pathlib import Path
from tempfile import TemporaryDirectory

from pypdf import PdfReader, PdfWriter, Transformation
from reportlab.lib import colors
from reportlab.lib.enums import TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.pdfgen import canvas
from reportlab.platypus import (SimpleDocTemplate, Paragraph, Spacer, Table,
                               TableStyle, PageBreak, Preformatted)

BLUE = colors.HexColor('#18486b')
GRAY = colors.HexColor('#526371')
BODY = ParagraphStyle('body', fontName='Helvetica', fontSize=10, leading=14,
                      spaceAfter=9, textColor=colors.HexColor('#1d2935'))
SMALL = ParagraphStyle('small', parent=BODY, fontSize=8.5, leading=11)
TITLE = ParagraphStyle('title', parent=BODY, fontName='Helvetica-Bold',
                       fontSize=23, leading=29, spaceAfter=16, textColor=BLUE)
HEADING = ParagraphStyle('heading', parent=TITLE, fontSize=16, leading=21,
                         spaceAfter=14)
SUB = ParagraphStyle('sub', parent=HEADING, fontSize=11, leading=15, spaceAfter=7)
MONO = ParagraphStyle('mono', fontName='Courier', fontSize=7.4, leading=10,
                      spaceAfter=10, alignment=TA_LEFT)


def paragraph(text, style=BODY):
    return Paragraph(escape(str(text)), style)


def table(rows, widths, header=True):
    value = Table([[paragraph(v, SMALL) for v in row] for row in rows],
                  colWidths=widths, repeatRows=int(header), hAlign='LEFT')
    commands = [('VALIGN', (0, 0), (-1, -1), 'TOP'),
                ('LEFTPADDING', (0, 0), (-1, -1), 8),
                ('RIGHTPADDING', (0, 0), (-1, -1), 8),
                ('TOPPADDING', (0, 0), (-1, -1), 7),
                ('BOTTOMPADDING', (0, 0), (-1, -1), 7),
                ('LINEBELOW', (0, 0), (-1, -1), .4, colors.HexColor('#dbe1e6'))]
    if header:
        commands += [('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#e8f0f6'))]
    value.setStyle(TableStyle(commands))
    return value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--packet', type=Path, required=True)
    parser.add_argument('--plots', type=Path, help='original 28-page vector plot PDF')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    packet = args.packet.resolve()
    report = json.loads((packet / 'report.json').read_text())
    integration = json.loads((packet / 'integration-report.json').read_text())
    physical = json.loads((packet / 'qualification-report.json').read_text())
    plots = PdfReader(args.plots or packet / 'plots-only.pdf')
    figures = report['figures']
    assert len(plots.pages) == len(figures) == 28
    assert integration['status'] == 'passed_supported_engineering_checks'
    assert len(integration['runs']) == 20 and all(r['passed'] for r in integration['runs'])
    assert len(physical['runs']) == 21 and all(r['exit'] == 0 for r in physical['runs'])
    assert report['compiled_source_fingerprint'] == physical['compiled_source_fingerprint']
    measured = {m['label']: m for m in physical['measurements']}
    width = A4[0] - 96
    story = []

    def text(s):
        story.append(paragraph(s))

    def section(s):
        if story:
            story.append(PageBreak())
        story.append(paragraph(s, HEADING))

    def sub(s):
        story.append(paragraph(s, SUB))

    def add_table(rows, fractions, header=True):
        story.append(table(rows, [width * f for f in fractions], header))
        story.append(Spacer(1, 10))

    story.append(paragraph('AthenaK shock encounters\nImplementation and test report', TITLE))
    text('Recorded local evidence | 8 October 2026 | Branch: codex/shock-multiphase')
    sub('What was implemented')
    text('On-the-fly shock detection and persistent per-particle shock-encounter histories. '
         'The live C++/Kokkos detector samples the completed fluid state; the tracker groups '
         'particle observations into episodes; migration, output and restart preserve those episodes.')
    add_table([['Component', 'Purpose'],
        ['Live mesh detector', 'Port the frozen ShockIt numerical contract and retain the unreduced shock region.'],
        ['Encounter tracker', 'Record entry/exit brackets, positions, Mach estimates, normals and validity/censoring flags.'],
        ['Persistence', 'Move complete typed state with each particle; use bounded event output and source-bound restart sidecars.']], [.28, .72])
    sub('Current result')
    text('Local engineering checks pass. The simple Mach 2 pilot records 128 completed episodes '
         'for 128 particles in each of eight seeds. Strong and boosted fixtures expose additional '
         'episodes and strength bias. Human scientific approval and R1 release acceptance remain pending.')
    sub('How to read this report')
    text('Pages 2-8 explain the test setup, methods, measurements, persistence, performance and '
         'reproduction. The following 28 figure pages retain the original vector plots and add '
         'plain-language reading notes and case-specific results.')
    text('The existing lagrangian_mc tracer makes stochastic parent-cell hops. Positions are cell '
         'centers. Lines in trajectory plots are sample sequences, and timing brackets describe '
         'neighboring shock-region observations. Neither is an exact continuous material crossing.')
    text('Evidence source: production commit 3de18d97; selected base multiphase_turb '
         '526a187f938e38634e89a88192ebe9731c4d36bb. This report adds explanation to recorded results; '
         'no simulations or numerical thresholds were changed to produce it.')

    section('2. The primary physical experiment')
    text('A translating planar shock is evolved in a thin two-dimensional domain. Independent '
         'Rankine-Hugoniot jump relations define its initial states and known speed. This is '
         'genuine 2D fluid output, but the planar setup does not qualify curved or intersecting shocks.')
    add_table([['Setting', 'Pilot configuration'],
        ['Domain and grid', 'x: 0 to 1; y: 0 to 0.1; 128 x 4 cells; no AMR.'],
        ['Physics and numerics', 'Newtonian ideal-gas hydro, gamma=1.4; RK2; PLM; HLLC; CFL=0.3.'],
        ['Shock', 'Mach 2; initial front x=0.3; speed=0.5; x_front(t)=0.3+0.5t.'],
        ['Particles and duration', '128 tracers seeded in x=0.45 to 0.75; t_limit=0.25.'],
        ['Boundaries', 'Outflow along x; periodic along y.'],
        ['Diagnostic settings', 'width=4 cells; min_mach=1.1; capacity=4096; trace=true. These are pilot settings, not approved scientific tolerances.']], [.30, .70])
    sub('Independent conservation reference')
    text('For upstream Mach M and gamma g, the density ratio is '
         '(g+1) M^2 / ((g-1) M^2+2), and the pressure ratio is '
         '(2g M^2-(g-1))/(g+1). The fixture checks mass, momentum and energy '
         'conservation in the shock frame before running the simulation.')
    text('At M=2, g=1.4: density ratio=2.6666667 and pressure ratio=4.5. '
         'Upstream density and pressure are both 1. The downstream values are 2.6666667 '
         'and 4.5. Velocities are -1.8664319 upstream and -0.3874120 downstream.')
    text('The pressure and density plots therefore have a reference independent of the detector. '
         'P/rho is a temperature proxy in code units. The Eulerian density-midpoint location is '
         'a front diagnostic; it is not used to invent particle crossing times.')
    sub('Why particles encounter this front')
    text('The front moves right while the upstream fluid and sampled tracers move predominantly '
         'left. Their histories meet the resolved shock region. Staying inside that region for '
         'many timesteps should create one numerical episode, not one event per timestep.')

    section('3. Engineering tests: what each checks')
    text('The recent integration suite has 20 subprocess runs. Some must succeed; others must '
         'abort with a specific diagnostic. All produced the expected outcome. Test IDs refer '
         'to the implementation plan and recorded checks.')
    add_table([['Test / IDs', 'Method and recorded result'],
        ['Non-interference\nCFG01, PHY08, SYS02', 'Run identical uniform and planar problems with tracking on/off. Particle observations are exact matches; 17-digit conserved slices are byte-identical.'],
        ['Sampling\nINT06, INT08', 'Compare final-step shock samples with native particle thermo history. Generation, time, position and density agree exactly.'],
        ['Negative controls\nCFG06, PHY09', 'Uniform flow, contact at constant pressure/velocity, and no particles. All emit zero events.'],
        ['Local migration\nINT01, IO01', 'Use multiple local mesh blocks. All initial IDs survive; 113 particles change parent block. This is not an inter-rank runtime test.'],
        ['Birth/censoring\nEVT02, EVT05', 'Seed initially or later inside the shock region. Missing valid pre-entry evidence is marked left-censored.'],
        ['All-phase restart\nRST01, RST07', 'Checkpoint outside, active and pending-exit phases. Resume and compare exact particle rows and canonical events.'],
        ['Damage/replay\nRST03/04, IO08', 'Remove/corrupt sidecars and damage retained history. Expected explicit rejection; immutable replay checks prevent conflicting event timelines.'],
        ['Capacity / configuration\nIO03/04, CFG', 'Force capacity=1 and an overflow; request unsupported RK3. Both explicitly reject instead of silently dropping state or accepting unqualified integration.']], [.30, .70])
    text('The native encounter unit also passes state-transition checks. Separate earlier evidence '
         'contains 12 explicit-rejection runs, 4 replay runs, 18 actual-method boundary cases '
         '(with some simplified harness context), and 13 actual production domain-exit records. '
         'These categories must not be summed as independent physical shock passages.')

    section('4. Detector parity and independent transport')
    sub('ShockIt parity: fidelity of the numerical port')
    text('53 fixtures compare the native detector against the frozen ShockIt local_pass reference. '
         'Shock-region masks must match exactly. Floating quantities use unchanged absolute '
         'and relative tolerances of 1e-12. All 53 pass.')
    text('The fixtures exercise constant fields, ideal shock jumps, gradients/normals, estimator '
         'thresholds and invalid inputs. Matching the reference validates the port; it does '
         'not independently certify that ShockIt or the encounter association is scientifically correct.')
    sub('Uniform transport: an independent physical probability model')
    text('For resolved positive uniform velocity u, the expected cell-hop probability per step '
         'is p=u*dt/dx. Summing the Bernoulli hops gives independent reference mean sum(p) '
         'and variance sum(p*(1-p)). This uses physical flux probabilities, not the detector or RNG implementation.')
    o = report['independent_transport']
    add_table([['Quantity', 'Reference', 'Measured'],
        ['Particles / steps', '', f"{o['particles']} tracers / {o['steps']} steps"],
        ['Mean cell hops', f"{o['expected_hops']:.6f}", f"{o['observed_hops_mean']:.6f}"],
        ['Hop variance', f"{o['expected_variance']:.6f}", f"{o['observed_hops_variance']:.6f}"],
        ['Mean / variance deviation', 'Existing +/-6 standard-error diagnostic', f"{o['mean_standard_errors']:.2f} / {o['variance_standard_errors']:.2f} standard errors"]], [.35, .30, .35])
    text('Both diagnostics pass. TRC11 also checks 119,808 actual transfer records: source parent, '
         'step start/end times, probabilities, random draw, selected direction and realized hop. '
         'Maximum probability discrepancy is 6.11e-16 against the unchanged 1e-12 tolerance.')
    text('This fixed-seed ensemble check is not proof of RNG independence at research scale. '
         'The inherited additive RNG key has alias/correlation questions that remain open. '
         'There is no deterministic per-particle material-path oracle.')

    section('5. Physical recovery and strength limitations')
    text('All 21 held-out fixtures executed. They vary eight Mach 2 seeds, Mach=1.2/1.5/3/5, '
         'gamma=5/3, stationary/boosted/reflected fronts, nx=64/256, CFL=0.15/0.45 and one '
         'thin-3D case. Run success is distinct from physical acceptance.')
    groups = [('Mach 2: each of eight seeds', 'seed11'), ('Mach 1.2', 'mach1.2'),
              ('Mach 1.5', 'mach1.5'), ('Mach 3', 'mach3'), ('Mach 5', 'mach5'),
              ('Boosted Mach 2', 'boost2')]
    rows = [['Fixture', 'Complete', 'Unique IDs', 'Suspended']]
    for title, key in groups:
        m = measured[key]
        rows.append([title, m['completed_events'], m['emitted_particles'], m['terminal_suspended']])
    rows.append(['Other gamma/speed/direction/resolution/CFL/thin-3D cases', '128 each', '128 each', '0'])
    add_table(rows, [.49, .17, .17, .17])
    text('Mach 5 gives six additional episodes for IDs 12, 43, 47, 57, 74 and 115. '
         'Boosted Mach 2 gives five additional episodes for IDs 3, 52, 67, 80 and 102, '
         'plus one suspended encounter. These are preserved limitations against the simple '
         'single-front expectation; stochastic returns and weak region lobes need scientific interpretation.')
    text('The boosted case adds 2 to the fluid and front velocities, retaining Mach 2. '
         'Its differing counts and bias expose frame dependence that needs qualification.')
    rows = [['Case', 'Pressure-Mach error range']]
    for title, key in [('Mach 2, seed 11', 'seed11'), ('Mach 1.2', 'mach1.2'),
                       ('Mach 5, all episodes', 'mach5'), ('Boosted Mach 2', 'boost2')]:
        low, high = measured[key]['pressure_mach_error_range']
        rows.append([title, f'{low:+.4f} to {high:+.4f}'])
    add_table(rows, [.49, .51])
    text('Error means estimated Mach minus prescribed Mach. Strength uses the first valid inside '
         'sample, with fixed endpoints that are not certified plateaus. Weak secondary episodes '
         'remain in the Mach 5 range. No outliers or tolerances were adjusted for acceptance.')

    section('6. Encounter timing, output and restart')
    text('An episode begins from valid unreduced-region membership, with an upstream-state '
         'entry guard where prior evidence exists. A first outside observation starts pending '
         'exit; two consecutive valid outside samples close/rearm. Earlier re-entry merges '
         'and flags flicker. Invalid samples cannot close an episode.')
    text('For valid pre/post evidence, entry bracket=[pre.time, first.time] and exit '
         'bracket=[last.time, post.time]. These bound adjacent sampled membership changes, '
         'not an exact material shock crossing. Cadence limitations, uncertain endpoints and '
         'direction ambiguity are explicit flags. Birth inside is left-censored; an unfinished '
         'run leaves a separately readable suspended right-censored view.')
    add_table([['Restart phase', 'Generation', 'Exact resumed rows', 'Exact events'],
        *[[r['phase'], integration['checkpoint_generations'][r['phase']],
           f"{r['rows_exact']:,}", r['events_exact']] for r in report['restart_results']]],
        [.27, .19, .30, .24])
    text('The strict reader detects no duplicate canonical event identities in these runs. '
         'Event keys are (lineage, unsigned 64-bit particle ID, unsigned 64-bit sequence). '
         'Complete encounter state migrates/compacts with the particle. Actual outflow records '
         'are emitted before particle deletion.')
    text('Per-generation event chunks have bounded capacity and immutable commits. Restart '
         'pairs the native checkpoint with a full encounter sidecar and requires the retained '
         'history prefix. Checks bind implementation, configuration, time, count, native layout '
         'and partition; missing, changed or incompatible data reject explicitly.')
    text('This is same-source, same-layout restart evidence. It is not an all-rank atomic '
         'transaction or power-loss/storage-fault qualification. MPI builds/singleton checks '
         'pass; the local two-rank launcher failed before AthenaK because of runtime '
         'compatibility. True distributed migration/restart and target filesystems remain untested.')

    section('7. Measured laptop performance')
    text('Uniform flow, 1,024 tracers, t_limit=0.1; two meshes and three modes. Three repetitions '
         'per mode/mesh, with execution order rotated. All 18 benchmark runs completed. '
         'Times below are whole-process medians, including startup and ordinary outputs.')
    values = {(r['nx'], r['mode']): r for r in report['benchmark_summary']}
    rows = [['Mesh', 'Tracking off', 'Tracking on', 'Tracking + trace']]
    for nx in (64, 256):
        cells = [f'{nx} x 4']
        for mode in ('off', 'tracking', 'trace'):
            v = values[nx, mode]
            cells.append(f"{v['median_whole_process_seconds']:.5f} s" +
                         (f" ({v['ratio_to_tracking_off']:.2f}x)" if mode != 'off' else ''))
        rows.append(cells)
    add_table(rows, [.17, .23, .30, .30])
    text('Tracking costs about 3.1-3.6x in these tiny jobs. Full observation/transfer traces '
         'cost about 25-27x. Durable empty event chunks are still committed every generation; '
         'trace mode writes every sampled observation and transfer. These costs must be '
         'visible when choosing operational output settings.')
    text('Output-cost figures report actual shock diagnostic bytes and file counts, excluding '
         'fluid/thermo output. Resident memory was not measured. Short local jobs do not '
         'establish sustained production overhead, scaling, MPI/GPU performance or costs on RUSTY/Frontier.')
    sub('Reproduction evidence and report scope')
    text('The initial fresh engineering and physical suites took 13.83 and 28.05 seconds. '
         'Their raw evidence occupied about 1.2 GB. The plotted workflow was independently '
         'repeated without reuse in a fresh directory: all 21 physical measurement sets and '
         '27 PNGs reproduced exactly, including output cost. The runtime figure retained fresh measurements.')
    text('The first plot attempt hit a Matplotlib layout error on an all-zero restart axis. '
         'Only the display scale/label coordinates were fixed; exact-equality assertions '
         'were retained. A packaging count was corrected from 26 to 27 matching figures '
         'because the output-cost plot also reproduced. These were reporting issues, not numerical fixes.')
    sub('Remaining gates')
    text('Human review and acceptable recovery/duplicate/timing/strength thresholds; RNG and '
         'frame dependence; sustained encounters and turbulence; true inter-rank operation; '
         'threaded/GPU/single precision and cluster budgets. AMR, cooling/source terms, '
         'MHD and relativistic physics require separate qualification. Material loops/strings are excluded.')

    section('8. Run the tests and locate the evidence')
    text('Run from the development worktree using the existing Release executable. Use new '
         'work and output directories each time. The existing ShockIt Python environment '
         'contains NumPy, Matplotlib and h5py; the PDF authoring step uses the separate '
         'bundled document runtime, with ReportLab and pypdf.')
    commands = '''cd /Users/hitesh/.codex/worktrees/shock-multiphase/athenak_nyu
MPLCONFIGDIR=/tmp/athenak-mpl \\
/Users/hitesh/hitesh/git/shockit/.venv/bin/python \\
  tst/shocks-mc/laptop_visuals.py \\
  --executable build-release/src/athena \\
  --workdir /tmp/athenak-report-run-new \\
  --output /tmp/athenak-report-plots-new
open /tmp/athenak-report-plots-new/index.html'''
    story.append(Preformatted(commands, MONO))
    text('The existing plot runner generates the raw tests, gallery and plot-only PDF. To '
         'produce this explained report from a new packet, run the following separate step:')
    story.append(Preformatted('''/Users/hitesh/.cache/codex-runtimes/\\
codex-primary-runtime/dependencies/python/bin/python3 \\
  tst/shocks-mc/explained_pdf.py \\
  --packet /tmp/athenak-report-plots-new \\
  --plots /tmp/athenak-report-plots-new/all-plots.pdf \\
  --output /tmp/athenak-report-explained.pdf''', MONO))
    text('That PDF step consumes existing results; it does not rerun the simulation. For a '
         'fresh build: cmake -S . -B build-release -DCMAKE_BUILD_TYPE=Release '
         '-DAthena_SHOCK_TESTS=ON; then cmake --build build-release -j8.')
    sub('Source and evidence identity')
    text('Core implementation fingerprint:')
    story.append(Preformatted(report['compiled_source_fingerprint'], MONO))
    text('Executable SHA256:')
    story.append(Preformatted(report['executable_sha256'], MONO))
    for line in ['Pilot input: inputs/tests/shock_mc.athinput',
                 'Test runners: tst/shocks-mc/integration.py, qualification.py, laptop_visuals.py',
                 'Recorded packet: evidence/shocks-mc/laptop-visuals/',
                 'Parity/persistence review: evidence/shocks-mc/review3/',
                 'Schema/support: docs/shocks-mc/schema.md, README.md, human-review.md',
                 'Recorded versions: Python 3.14.0, NumPy 2.5.0, Matplotlib 3.11.0, h5py 3.16.0']:
        story.append(paragraph(line, SMALL))
    text('Scientific sign-off: pending. CPU Serial/double/RK2/uniform Cartesian ideal-gas '
         'hydro is the measured local scope. Earlier evidence retains its own source identity.')

    total = 8 + len(figures)
    def page_frame(c, doc):
        c.setFont('Helvetica', 8)
        c.setFillColor(GRAY)
        c.drawString(48, A4[1] - 29, 'ATHENAK  |  SHOCK ENCOUNTER TEST REPORT')
        c.drawString(48, 27, 'Scientific review pending  |  Recorded local CPU evidence')
        c.drawRightString(A4[0] - 48, 27, f'{doc.page} / {total}')
    output = args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(prefix='athenak-explained-pdf-') as directory:
        front = Path(directory) / 'explanation.pdf'
        SimpleDocTemplate(str(front), pagesize=A4, leftMargin=48, rightMargin=48,
                          topMargin=55, bottomMargin=48, title='AthenaK shock-encounter test report',
                          author='AthenaK implementation evidence').build(
                              story, onFirstPage=page_frame, onLaterPages=page_frame)
        intro = PdfReader(front)
        assert len(intro.pages) == 8, f'Explanation overflow: {len(intro.pages)} pages'
        writer = PdfWriter()
        writer.append(intro)
        for i, (figure, source) in enumerate(zip(figures, plots.pages)):
            pw, ph = float(source.mediabox.width), float(source.mediabox.height)
            page_width, page_height = max(pw + 72, A4[0]), ph + 200
            page = Path(directory) / f'figure-{i}.pdf'
            c = canvas.Canvas(str(page), pagesize=(page_width, page_height))
            c.setFillColor(BLUE)
            c.setFont('Helvetica-Bold', 14)
            c.drawString(36, page_height - 40, f'Figure {i + 1:02d}  |  ' + figure['name'].replace('_', ' '))
            c.setFillColor(GRAY)
            c.setFont('Helvetica', 9)
            c.drawString(36, page_height - 61, 'Read the plotted quantities with the notes below; code time and spatial code units apply.')
            note = figure['caption']
            name = figure['name']
            if name.endswith('_parents') or name.endswith('_strength_time'):
                key = name.removesuffix('_parents').removesuffix('_strength_time')
                m = measured[key]
                note += (f" Case: input Mach={m['expected_input_mach']}; {m['completed_events']} complete "
                         f"episodes / {m['emitted_particles']} unique particles; {m['terminal_suspended']} suspended. "
                         'Each trace is a selected ID, not an ensemble-averaged material path.')
            if name == 'episode_counts':
                note += ' All eight Mach 2 seeds recover 128/128; Mach 5 gives 134/128; boosted Mach 2 gives 132/127 plus one suspended.'
            if name == 'restart_agreement':
                note += ' All plotted differences are exactly zero. The 1e-14 display unit is not a tolerance; comparisons use exact equality.'
            if name.startswith('2d_'):
                note += ' Top panel is the initial state (t=0); later panels use the actual recorded generation/time. All four y rows are native grid values.'
            if name == 'timing_bracket_widths':
                note += ' Coincident curves can overlap: halving CFL and doubling x resolution give similar observation cadences in these fixtures.'
            if name == 'uniform_transport':
                note += ' The reference is a sum of per-step Bernoulli probabilities. The measured mean/variance deviations are -0.89/-2.40 standard errors.'
            p = paragraph(note, SMALL)
            _, height = p.wrap(page_width - 72, 90)
            assert height <= 84, (name, height)
            p.drawOn(c, 36, 111 - height)
            c.setFillColor(GRAY)
            c.setFont('Helvetica', 8)
            c.drawString(36, 24, 'Scientific review pending  |  Original vector plot retained')
            c.drawRightString(page_width - 36, 24, f'{9 + i} / {total}')
            c.showPage()
            c.save()
            explained = PdfReader(page).pages[0]
            explained.merge_transformed_page(source, Transformation().translate(tx=36, ty=120))
            writer.add_page(explained)
        writer.add_metadata({'/Title': 'AthenaK shock encounters - implementation and test report',
                             '/Author': 'AthenaK implementation evidence',
                             '/Subject': 'Engineering methods, measured outcomes, physical limitations and 28 explained plots'})
        with output.open('wb') as stream:
            writer.write(stream)
    check = PdfReader(output)
    assert len(check.pages) == total
    contents = '\n'.join(page.extract_text() or '' for page in check.pages)
    for text_check in ('Rankine-Hugoniot', 'Scientific', '3de18d97', '134', 'same-layout', '119,808'):
        assert text_check in contents, text_check
    print(json.dumps({'output': str(output), 'pages': total, 'explanation_pages': 8,
                      'explained_vector_figures': len(figures), 'scientific_signoff': 'pending'}))


if __name__ == '__main__':
    main()
