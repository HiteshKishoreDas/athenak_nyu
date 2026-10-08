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


def figure_description(name, report, measured, packet):
    """Test setup, plot reading, measured result and caveat for this exact figure."""
    pilot = ('Mach 2 planar shock, gamma=1.4, RK2/PLM/HLLC, 128 x 4 cells, '
             '128 MC tracers, CFL=0.3, t_limit=0.25. The known front is x=0.3+0.5t.')
    endpoints = ('Fixed sampling endpoints are not certified undisturbed plateaus. '
                 'Strength and physical passage acceptance remain pending.')
    if name.startswith('2d_'):
        field = {'2d_dens': 'density', '2d_pressure': 'pressure',
                 '2d_temperature': 'temperature proxy P/rho'}[name]
        plateau = {'2d_dens': 'density 1 upstream and 2.6667 downstream',
                   '2d_pressure': 'pressure 1 upstream and 4.5 downstream',
                   '2d_temperature': 'P/rho 1 upstream and 1.6875 downstream'}[name]
        return (pilot + ' This page shows native ' + field + ' snapshots.',
                'Read top to bottom from t=0 to t=0.25; the middle panels show their actual '
                'recorded times. Color is ' + field + '. Dashed red marks the analytic front; '
                'white dots are sampled particles inside the unreduced shock region.',
                'The independent jump conditions prescribe ' + plateau + '. The front moves '
                'right to x=0.425 at t=0.25. Synchronization checks separately match particle '
                'time, cell position and density exactly to native thermo observations.',
                'The y direction is stretched to display four real grid rows. White dots are '
                'particle samples, not a full mesh mask. This planar fixture does not qualify '
                'curved, interacting or turbulent 2D shocks.')
    if name == 'rh_profiles':
        return (pilot + ' Initial states satisfy mass, momentum and energy conservation in the shock frame.',
                'The four panels compare numerical density, pressure, x velocity and P/rho '
                'with dashed independent Rankine-Hugoniot plateau values at t=0.25.',
                'The reference front is x=0.425. Density ratio is 8/3 and pressure ratio is 4.5; '
                'velocities are -1.86643 upstream and -0.38741 downstream. The resolved '
                'profiles show finite shock thickness and small downstream oscillations.',
                'This checks the fluid against a physical reference independent of ShockIt. '
                'It does not establish a plateau-certified particle Mach estimator or an exact '
                'particle crossing time.')
    if name == 'shock_space_time':
        return (pilot + ' Native grid snapshots are written throughout the run.',
                'Horizontal coordinate is x; vertical coordinate is simulation time. Color '
                'shows density. The red dashed line is the prescribed Eulerian front; white '
                'dots mark synchronized inside-region particle observations.',
                'The density jump follows the translating front, and the sampled membership '
                'appears around the resolved shock. At the end the analytic front is x=0.425. '
                'The matching particle timestamp/position/density tests passed.',
                'Space-time interpolation of sparse grid outputs is a visualization. It is '
                'not a continuous material trajectory or a fitted crossing-time oracle.')
    if name == 'front_position':
        import csv
        rows = list(csv.DictReader((packet / 'data/front-position.csv').open()))
        errors = [float(r['error_cells']) for r in rows]
        return (pilot + ' The fluid front is measured by linearly interpolating its density midpoint.',
                'Left: midpoint location versus the known front. Right: midpoint minus '
                'analytic position, divided by the x cell width; zero means coincident locations.',
                f'The recorded midpoint offsets range from {min(errors):.3f} to '
                f'{max(errors):.3f} cells across {len(rows)} snapshots. Initial cell placement '
                'contributes an offset; the midpoint is required to have a unique crossing.',
                'This is an Eulerian front diagnostic, not a particle timing estimate. No '
                'physical acceptance tolerance is inferred from the displayed offset.')
    if name.endswith('_parents') or name.endswith('_strength_time'):
        key = name.removesuffix('_parents').removesuffix('_strength_time')
        m = measured[key]
        speed = 2.5 if key == 'boost2' else 0 if key == 'stationary' else -.5 if key == 'reflected' else .5
        origin = .7 if key == 'reflected' else .3
        nx = 64 if key == 'nx64' else 256 if key == 'nx256' else 128
        duration = .45 if key == 'mach1.2' else .3 if key == 'mach1.5' else .25
        setup = (f"Input Mach {m['expected_input_mach']:g}, {nx} x 4 cells, 128 MC tracers, "
                 f'gamma=1.4, CFL=0.3, t_limit={duration:g}. Analytic front: '
                 f'x={origin:g}+({speed:g})t. The plotted IDs are selected samples.')
        low, high = m['pressure_mach_error_range']
        result = (f"{m['completed_events']} completed numerical episodes for "
                  f"{m['emitted_particles']} unique particles; {m['terminal_suspended']} suspended. "
                  f'First-inside pressure-Mach error is {low:+.4f} to {high:+.4f} '
                  '(estimate minus input Mach).')
        if key in ('mach5', 'boost2'):
            result += f" {len(m['particles_with_multiple_episodes'])} particles have a second episode; scientific interpretation remains pending."
        if name.endswith('_parents'):
            read = ('Top: sampled cell-center x versus time, with inside-region dots and the '
                    'known dashed front. Bottom: parent x minus analytic front x. Stair steps '
                    'show discrete cell hops; selected extra-episode/censored IDs are included.')
            limit = ('MC sample sequences are not continuous material paths. A numerical '
                     're-entry need not represent a distinct physical shock. No exact crossing '
                     'time is assigned, and a suspended encounter is retained at cutoff.')
        else:
            read = ('Top: outside/active/pending-exit tracker state. Lower panels: pressure '
                    'and temperature Mach for every valid inside observation. Dashed lines '
                    'mark input Mach; gaps mean no valid inside strength sample.')
            limit = ('The lower panels include follow-up samples, whereas each emitted event '
                     'uses its first valid inside sample as representative. ' + endpoints)
        return setup, read, result, limit
    if name == 'episode_counts':
        return ('21 held-out runs: eight Mach 2 seeds, four other Mach numbers, gamma, front '
                'speed/direction, resolution, CFL and one thin-3D variation; 128 initial particles each.',
                'Blue bars count completed numerical episodes; orange adds partial episodes. '
                'Black dots count unique IDs with events; red crosses count suspended views. '
                'The dashed reference is 128 initially seeded particles.',
                'Every Mach 2 seed gives 128/128. Mach 5 gives 134 episodes/128 IDs. Boosted '
                'Mach 2 gives 132/127 plus one suspended. Mach 1.5 gives 127 plus one suspended. '
                'Other held-out cases give 128/128.',
                'Successful execution is not scientific acceptance. Extra episodes violate '
                'the simple single-front counting expectation; unfinished encounters are '
                'censored evidence, not discarded failures.')
    if name == 'mach_errors':
        return ('First-inside representative strengths from all emitted events in the 21 '
                'held-out cases, compared with independently prescribed input Mach.',
                'Top boxes use pressure Mach; bottom boxes use temperature Mach. The y axis '
                'is estimated minus input Mach. Zero is unbiased. Boxes, whiskers and outliers '
                'retain every emitted event, including weak secondary detections.',
                'Pressure-Mach errors: seed-11 Mach 2, -0.0163 to +0.0118; Mach 1.2, -0.1031 '
                'to -0.0157; Mach 5, -3.9657 to +0.0594; boosted Mach 2, -0.1091 to -0.0771.',
                'A box plot diagnoses bias/spread, not an approved error budget. No outliers '
                'were removed, and the full vertical range intentionally exposes strong-case '
                'failures. ' + endpoints)
    if name == 'timing_bracket_widths':
        widths = {}
        for key in ('nx64', 'seed11', 'nx256', 'cfl0.15', 'cfl0.45'):
            values = [e['entry_region_bracket'][1]-e['entry_region_bracket'][0]
                      for e in measured[key]['event_evidence']]
            widths[key] = max(values)
        return ('Compare x resolutions 64/128/256 and CFL 0.15/0.3/0.45 for the same Mach 2 '
                'planar test, measuring each event\'s recorded entry and exit bracket width.',
                'Left: entry widths; right: exit widths. Each cumulative curve shows the '
                'fraction of events below a width. A left shift means tighter observation '
                'brackets. Coincident curves can overlap.',
                f"Largest recorded entry widths: nx64 {widths['nx64']:.6g}, baseline "
                f"{widths['seed11']:.6g}, nx256 {widths['nx256']:.6g} code-time units. "
                'All five cases record 128 episodes for 128 IDs.',
                'Widths describe adjacent mask observations, not errors against an exact '
                'material crossing. Halving CFL and doubling resolution can yield similar '
                'cadence; that does not prove converged physical passage histories.')
    if name == 'encounter_brackets':
        return ('Selected particle IDs in baseline Mach 2, Mach 5 and boosted Mach 2, including '
                'IDs with rearmed episodes and terminal suspended evidence.',
                'Blue is the first episode; red is a later episode for the same ID. Pale '
                'thick ends are entry/exit brackets. Solid segments join first and last inside '
                'observations; a black cross denotes a suspended terminal view.',
                'Baseline has one event per ID. Mach 5 has six IDs with a second episode; '
                'the boosted case has five and one suspended encounter. The tracker preserves '
                'these distinct identities and timing evidence rather than hiding them.',
                'A solid segment is not exact physical residence: flicker can be merged, and '
                'MC returns are stochastic. Brackets limit sampled region transitions, not '
                'exact material shock crossings.')
    if name == 'uniform_transport':
        o = report['independent_transport']
        return ('Uniform periodic flow, velocity u=1, 64 x 4 cells, 1,024 tracers, t_limit=0.25. '
                'The independent reference uses per-step Bernoulli probabilities p=u*dt/dx.',
                'Left: measured net-hop histogram versus the probability distribution of the '
                'Bernoulli sum. Right: mean and variance differences measured in their '
                'predeclared standard errors, with the existing +/-6 diagnostic lines.',
                f"Expected mean {o['expected_hops']:.4f}, measured {o['observed_hops_mean']:.4f}; "
                f"expected variance {o['expected_variance']:.4f}, measured {o['observed_hops_variance']:.4f}. "
                'Mean/variance deviations are -0.89/-2.40 standard errors; both pass.',
                'This is a fixed-seed ensemble transport check, not proof of RNG independence '
                'at research scale. The inherited additive RNG key has alias/correlation '
                'questions; no deterministic per-particle path oracle is used.')
    if name == 'negative_controls':
        return ('Uniform flow, a contact with density ratio 2:1 at constant pressure and '
                'velocity, and a zero-particle run. A contact should not be mistaken for a shock.',
                'Left: sampled density versus parent-cell x at the first and last generations '
                'for uniform/contact runs. Right: emitted event counts for all three controls.',
                'All three emit zero shock episodes. Tracking-on/off uniform and planar runs '
                'also preserve exact native particle observations and byte-identical '
                '17-digit conserved slices in the separate non-interference checks.',
                'These controls check specific false detections and empty-population behavior. '
                'They do not bound the false-positive rate in turbulence or complex shock geometry.')
    if name == 'restart_agreement':
        return ('Checkpoint the baseline run in outside, active and pending-exit phases '
                '(generations 1, 61 and 68), then resume with identical source/configuration/partition.',
                'Each panel plots maximum absolute position and density differences from '
                'uninterrupted particle observations. The display unit is 1e-14 code units; '
                'comparisons use exact equality, not a 1e-14 tolerance.',
                'All differences are exactly zero. Outside/active/pending resumes match '
                '41,600 / 33,920 / 33,024 rows, respectively, and each preserves all 128 '
                'canonical events with no duplicate identities.',
                'Requires the matching implementation, native ABI, rank count and mesh '
                'partition plus retained event prefix. No cross-rank, all-rank atomic '
                'transaction or power-loss/storage-fault qualification is claimed.')
    if name == 'runtime_cost':
        return ('Uniform flow, 1,024 particles, t_limit=0.1, grids 64 x 4 and 256 x 4. '
                'Three trials per off/tracking/trace mode; order is rotated.',
                'Left: median whole-process wall seconds, with dots for actual repeats. '
                'Right: runtime divided by the matching tracking-off median. Blue is off, '
                'orange tracking with trace off, and green full diagnostic traces.',
                '64 x 4: off 0.01163 s, tracking 0.04171 s (3.59x), trace 0.28880 s (24.84x). '
                '256 x 4: off 0.04174 s, tracking 0.12862 s (3.08x), trace 1.10772 s (26.54x).',
                'Startup, ordinary output and durable empty event chunks are included. '
                'Full traces are expensive small-job qualification diagnostics. These '
                'timings do not establish sustained RUSTY/Frontier or MPI/GPU performance.')
    if name == 'output_cost':
        import statistics
        summary = {}
        for nx in (64, 256):
            rows = [r for r in report['benchmark_runs'] if r['mesh'][0] == nx and r['mode'] == 'trace']
            summary[nx] = (statistics.median(r['shock_output_bytes'] for r in rows)/1e6,
                           int(statistics.median(r['shock_output_files'] for r in rows)))
        return ('The same two-mesh, 1,024-particle benchmark, inspecting shock diagnostic '
                'directories after each off/tracking/trace run.',
                'Left: bytes written in decimal MB; right: file count. Bars use the same '
                'three diagnostic modes as the runtime plot. Counts include immutable '
                'per-generation chunks and optional observation/transfer CSVs.',
                f"Trace mode writes {summary[64][0]:.2f} MB / {summary[64][1]} files on 64 x 4, "
                f"and {summary[256][0]:.2f} MB / {summary[256][1]} files on 256 x 4. "
                'Tracking-off writes no shock diagnostic directory.',
                'Fluid and legacy thermo outputs are excluded from these totals. Empty '
                'event generations still require durable files. Resident memory and '
                'production-filesystem behavior were not measured.')
    raise ValueError('Missing local explanation: ' + name)


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
    text('Pages 2-29 pair each original vector plot with its own test setup, reading guide, '
         'measured result and limitation directly above it. General methods, persistence '
         'details and reproducible run instructions follow in the appendix on pages 30-36.')
    text('The existing lagrangian_mc tracer makes stochastic parent-cell hops. Positions are cell '
         'centers. Lines in trajectory plots are sample sequences, and timing brackets describe '
         'neighboring shock-region observations. Neither is an exact continuous material crossing.')
    text('Evidence source: production commit 3de18d97; selected base multiphase_turb '
         '526a187f938e38634e89a88192ebe9731c4d36bb. This report adds explanation to recorded results; '
         'no simulations or numerical thresholds were changed to produce it.')

    section('Appendix A. The primary physical experiment')
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

    section('Appendix B. Engineering tests: what each checks')
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

    section('Appendix C. Detector parity and independent transport')
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

    section('Appendix D. Physical recovery and strength limitations')
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

    section('Appendix E. Encounter timing, output and restart')
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

    section('Appendix F. Measured laptop performance')
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

    section('Appendix G. Run the tests and locate the evidence')
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
        number = 1 if doc.page == 1 else doc.page + len(figures)
        c.drawRightString(A4[0] - 48, 27, f'{number} / {total}')
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
        writer.add_page(intro.pages[0])
        for i, (figure, source) in enumerate(zip(figures, plots.pages)):
            pw, ph = float(source.mediabox.width), float(source.mediabox.height)
            page_width = max(pw + 72, A4[0])
            column_width = (page_width - 96) / 2
            setup, read, result, caveat = figure_description(
                figure['name'], report, measured, packet)
            columns = [[('Test setup', setup), ('How to read the plot', read)],
                       [('Measured result', result), ('What this establishes', caveat)]]
            blocks = []
            for column in columns:
                parts = []
                for label, content in column:
                    for value, style in ((label, SUB), (content, BODY)):
                        p = paragraph(value, style)
                        _, height = p.wrap(column_width, 1000)
                        parts.append((p, height, 6 if style is SUB else 10))
                blocks.append(parts)
            explanation_height = max(sum(h + gap for p, h, gap in parts)
                                     for parts in blocks)
            page_height = ph + 133 + explanation_height
            page = Path(directory) / f'figure-{i}.pdf'
            c = canvas.Canvas(str(page), pagesize=(page_width, page_height))
            c.setFillColor(BLUE)
            c.setFont('Helvetica-Bold', 14)
            c.drawString(36, page_height - 34, f'Test {i + 1:02d}  |  ' +
                         figure['name'].replace('_', ' '))
            for column, parts in enumerate(blocks):
                y = page_height - 65
                x = 36 + column * (column_width + 24)
                for p, height, gap in parts:
                    y -= height
                    p.drawOn(c, x, y)
                    y -= gap
                assert y >= ph + 48, 'Explanation overlaps its plot'
            c.setFillColor(GRAY)
            c.setFont('Helvetica', 8)
            c.drawString(36, 24, 'Scientific review pending  |  Original vector plot retained')
            c.drawRightString(page_width - 36, 24, f'{2 + i} / {total}')
            c.showPage()
            c.save()
            explained = PdfReader(page).pages[0]
            explained.merge_transformed_page(
                source, Transformation().translate(tx=36, ty=48))
            writer.add_page(explained)
        for page in intro.pages[1:]:
            writer.add_page(page)
        writer.add_metadata({'/Title': 'AthenaK shock encounters - implementation and test report',
                             '/Author': 'AthenaK implementation evidence',
                             '/Subject': '28 plots with adjacent test setup, reading guide, measured results and limitations'})
        with output.open('wb') as stream:
            writer.write(stream)
    check = PdfReader(output)
    assert len(check.pages) == total
    contents = '\n'.join(page.extract_text() or '' for page in check.pages)
    for text_check in ('Rankine-Hugoniot', 'Scientific', '3de18d97', '134', 'same-layout', '119,808'):
        assert text_check in contents, text_check
    print(json.dumps({'output': str(output), 'pages': total, 'explanation_pages': 8,
                      'explained_vector_figures': len(figures), 'plot_local_explanations': 28,
                      'overview_pages': 1, 'appendix_pages': 7, 'scientific_signoff': 'pending'}))


if __name__ == '__main__':
    main()
