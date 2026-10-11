#!/usr/bin/env python3
"""Bounded forcing pilot; no stationarity/reference-closure claim from elapsed time.

Requires a serial double -DPROBLEM=sgs_verify binary with stage-consistent forcing.
Uses native-double checkpoint inventories, fixed driver inputs and explicit on/off
lanes. Every attempt is retained. Numerical thresholds are inherited unchanged.
"""
import argparse
import json
from pathlib import Path
import subprocess
import time

import verify


class ForcingVerification(verify.Verification):
    def run_forcing(self, name, enabled, sol_fraction, seed=42, tlim=1., rst=None,
                    oracle=False, nscalars=2, ndim=3, overrides=None, timeout=60,
                    retain_snapshots=True):
        config = verify.parameters((verify.SOURCE/'inputs/hydro/sgs_verify.athinput').read_text())
        config.update({'job/basename':name, 'problem/mode':'uniform',
            'problem/velocity':'0', 'problem/k0':'1e-6', 'problem/pressure0':'1',
            'hydro/gamma':'1.6666666666666667', 'hydro/nscalars':str(nscalars),
            'sgs/enabled':str(enabled).lower(), 'time/tlim':str(tlim),
            'time/cfl_number':'0.1', 'time/nlim':'5000', 'time/ndiag':'100',
            'output1/dt':'0.1', 'output1/data_format':'%24.17e',
            'output2/dt':'-1', 'output3/dt':'-1', 'output4/dt':'0.25',
            'turb_driving/type':'hydro', 'turb_driving/stage_consistent':'true',
            'turb_driving/dedt':'0' if oracle else '0.1',
            'turb_driving/tcorr':'0.1', 'turb_driving/dt_update':'0.01',
            'turb_driving/rseed':str(seed), 'turb_driving/nlow':'1',
            'turb_driving/nhigh':'3', 'turb_driving/sol_fraction':str(sol_fraction),
            'turb_driving/normalization':'edot'})
        if oracle:
            config['problem/forcing_oracle'] = 'true'
        for d in (1,2,3):
            config[f'mesh/nx{d}'] = '16' if d <= ndim else '1'
            config[f'meshblock/nx{d}'] = '8' if d <= ndim else '1'
        for key, value in (overrides or {}).items():
            if value is None:
                config.pop(key, None)
            else:
                config[key] = value
        wd = self.output/name
        wd.mkdir()
        inp = wd/'input.athinput'
        inp.write_text(verify.input_text(config))
        cmd = [str(self.exe), '-i', str(inp)]
        if rst:
            cmd = [str(self.exe), '-r', str(rst), f'job/basename={name}',
                   f'time/tlim={tlim}']
            if overrides and 'time/cfl_number' in overrides:
                cmd.append(f"time/cfl_number={config['time/cfl_number']}")
        row = dict(name=name, enabled=enabled, sol_fraction=sol_fraction,
                   seed=seed, command=cmd, input_sha256=verify.sha(inp), status='running',
                   timeout_seconds=timeout)
        self.results.append(row)
        self.save()
        started = time.monotonic()
        try:
            with (wd/'run.log').open('w') as log:
                proc = subprocess.run(cmd, cwd=wd, stdout=log, stderr=subprocess.STDOUT,
                                      timeout=timeout)
            row['returncode'] = proc.returncode
            if proc.returncode:
                raise RuntimeError(f'execution returned {proc.returncode}')
            snapshots = []
            first = last = None
            for path in sorted((wd/'rst').glob('*.rst')):
                snap = verify.restart(path)
                if (not snap['finite'] or snap['min_density'] <= 0 or
                        snap['min_thermal'] <= 0 or
                        (snap['min_sgs'] is not None and snap['min_sgs'] < 0)):
                    raise RuntimeError('inadmissible stored reservoir')
                if first is None:
                    first = snap
                last = snap
                if retain_snapshots:
                    snapshots.append(snap)
            if first is None:
                raise RuntimeError('no native checkpoints')
            if not retain_snapshots:
                snapshots = [first, last]
            hist = verify.history(wd/f'{name}.hydro.hst')
            row.update(actual_time=last['time'], cycles=last['cycle'])
            if verify.scaled_diff(last['time'],tlim) > 1e-12:
                raise RuntimeError('did not complete tlim within resource cap')
            row['checks'] = [dict(id='forcing.admissibility', passed=True)]
            initial, final = hist[0], hist[-1]
            work = final['force-work']-initial['force-work']
            # Histories supply stage-integrated applied forcing; native states supply E+K.
            expected_delta = work
            E0 = first['totals'][4]+(first['totals'][-1] if enabled else 0)
            E1 = last['totals'][4]+(last['totals'][-1] if enabled else 0)
            # On restart, history and checkpoint cadences can start at different times.
            # Compare final absolute inventories to checkpoint-header initial counters.
            if rst:
                original = verify.restart(rst)
                E0 = original['totals'][4]+(original['totals'][-1] if enabled else 0)
                work0 = float(original['config']['turb_driving/integrated_work'])
                expected_delta = final['force-work']-work0
                first = original
            self.check(row,'P24.actual_work',verify.scaled_diff(E1-E0,expected_delta),
                       verify.LIMITS['represented_energy_scaled_absolute'],
                       'represented inventory change equals actually applied forcing work')
            for v in (1,2,3):
                impulse0 = float(first['config'].get(f'turb_driving/integrated_momentum{v}','0'))
                impulse = final[f'force-mom{v}']-impulse0
                self.check(row,f'P24.momentum{v}',verify.scaled_diff(
                    last['totals'][v]-first['totals'][v],impulse),
                    verify.LIMITS['mass_momentum_scaled_absolute'],
                    'independent periodic momentum change versus applied impulse')
            self.check(row,'forcing.mass',verify.scaled_diff(first['totals'][0],last['totals'][0]),
                       verify.LIMITS['mass_momentum_scaled_absolute'],'periodic mass')
            for v in range(5,5+int(first['config'].get('hydro/nscalars','0'))):
                self.check(row,f'forcing.scalar{v-5}',verify.scaled_diff(
                    first['totals'][v],last['totals'][v]),
                    verify.LIMITS['scalar_mass_scaled_absolute'],'distinct passive scalar inventory')
            self.check(row,'forcing.history.thermal',verify.scaled_diff(
                last['thermal'],final.get('thermal-E',final['tot-E']-
                                         sum(final[f'{v}-KE'] for v in (1,2,3)))),
                verify.LIMITS['history_vs_restart_scaled_absolute'],'thermal inventory')
            if enabled:
                source0 = float(first['config'].get('sgs/integrated_production','0'))-float(
                    first['config'].get('sgs/integrated_dissipation','0'))
                source1 = final['sgs-production-integral']-final['sgs-dissipation-integral']
                self.check(row,'forcing.SGS_transfer',verify.scaled_diff(
                    last['totals'][-1]-first['totals'][-1],source1-source0),
                    verify.LIMITS['paired_transfer_scaled_absolute'],
                    'SGS production/dissipation ledger remains separate from resolved forcing')
                for label,value in [('sgs-E',last['totals'][-1]),
                                    ('represented-E',E1),('thermal-E',last['thermal'])]:
                    self.check(row,'forcing.history.'+label,verify.scaled_diff(value,final[label]),
                               verify.LIMITS['history_vs_restart_scaled_absolute'],
                               'online history versus independent native-double inventory')
            self.check(row,'forcing.history.resolved',verify.scaled_diff(last['totals'][4],final['tot-E']),
                       verify.LIMITS['history_vs_restart_scaled_absolute'],'resolved energy inventory')
            self.check(row,'forcing.history.kinetic',verify.scaled_diff(last['kinetic'],
                       sum(final[f'{v}-KE'] for v in (1,2,3))),
                       verify.LIMITS['history_vs_restart_scaled_absolute'],'kinetic inventory')
            row['final_inventories'] = dict(resolved=last['totals'][4],thermal=last['thermal'],
                kinetic=last['kinetic'],sgs=last['totals'][-1] if enabled else 0.,
                represented=E1,applied_work=final['force-work'])
            # A trend measures evidence sufficiency; it is never a stationarity certificate.
            tail = [h for h in hist if h['time'] >= 0.5*tlim]
            if len(tail) >= 3:
                ts = [h['time'] for h in tail]
                kes = [sum(h[f'{v}-KE'] for v in (1,2,3)) for h in tail]
                tm,km = sum(ts)/len(ts),sum(kes)/len(kes)
                slope = sum((t-tm)*(k-km) for t,k in zip(ts,kes))/sum((t-tm)**2 for t in ts)
                row['late_kinetic_trend'] = dict(samples=len(tail),slope=slope,
                    fractional_drift=slope*(ts[-1]-ts[0])/km if km else None,
                    window=[ts[0],ts[-1]],
                    qualification='pilot trend only; no independence/stationarity claim')
            row['status'] = 'pass' if all(c['passed'] for c in row['checks']) else 'fail'
            self.states[name] = (first,last,snapshots,hist)
        except subprocess.TimeoutExpired:
            row.update(status='fail',error=f'{timeout}-second cap exceeded; retained as incomplete')
        except Exception as exc:
            row.update(status='fail',error=str(exc))
        row['seconds'] = time.monotonic()-started
        self.save()
        print(name,row['status'],row.get('error',''),flush=True)
        return row

    def forcing_restart(self, name, original_name, enabled, sol_fraction, tlim):
        if original_name not in self.states:
            self.results.append(dict(name=name,status='blocked',error='pilot unavailable'))
            self.save()
            return
        paths = sorted((self.output/original_name/'rst').glob('*.rst'))
        chosen = next((p for p in paths if .2*tlim <= verify.restart(p)['time'] <= .7*tlim),None)
        if chosen is None:
            self.results.append(dict(name=name,status='blocked',error='no normal intermediate checkpoint'))
            self.save()
            return
        self.run_forcing(name,enabled,sol_fraction,tlim=tlim,rst=chosen)
        if name not in self.states:
            return
        a,b = self.states[original_name][1],self.states[name][1]
        checkpoint = verify.restart(chosen)
        row = dict(name='forcing.restart:'+name,status='running',
            checkpoint=str(chosen), checkpoint_time=checkpoint['time'],
            checkpoint_cycle=checkpoint['cycle'],
            checkpoint_ou_updates=checkpoint['forcing_metadata']['n_updates'],
            final_times=[a['time'],b['time']], final_cycles=[a['cycle'],b['cycle']],
            final_ou_updates=[a['forcing_metadata']['n_updates'],
                              b['forcing_metadata']['n_updates']])
        self.results.append(row)
        self.check(row,'forcing.restart.fields',verify.state_diff(a,b),
                   verify.LIMITS['checkpoint_restart_scaled_linf'],'native double active-cell state')
        diff = max(verify.scaled_diff(x,y) for key in a['forcing']
                   for x,y in zip(a['forcing'][key],b['forcing'][key]))
        self.check(row,'forcing.restart.modal_OU',diff,verify.LIMITS['checkpoint_restart_scaled_linf'],
                   'Fourier OU amplitudes, independent of normalized applied acceleration')
        self.check(row,'forcing.restart.metadata',
                   int(a['forcing_metadata']['fields'] != b['forcing_metadata']['fields']),
                   0, 'identical native modal driver metadata')
        self.check(row,'forcing.restart.OU_updates',
                   abs(a['forcing_metadata']['n_updates']-b['forcing_metadata']['n_updates']),
                   0, 'exact OU update count in native metadata')
        self.check(row,'forcing.restart.RNG',int(a['rng_state'] != b['rng_state']),0,
                   'native RNG state bytes')
        for label in ('force-work','force-mom1','force-mom2','force-mom3'):
            self.check(row,'forcing.restart.'+label,verify.scaled_diff(
                self.states[original_name][3][-1][label],self.states[name][3][-1][label]),
                verify.LIMITS['checkpoint_restart_scaled_linf'],'precise integrated forcing metadata')
        row['status'] = 'pass' if all(c['passed'] for c in row['checks']) else 'fail'
        self.save()

    def reject_contracts(self, original_name, suffix=''):
        """Reject legacy SGS forcing and every immutable restart override sampled."""
        input_path = self.output/original_name/'input.athinput'
        checkpoint = sorted((self.output/original_name/'rst').glob('*.rst'))[1]
        cases = [('legacy_forcing', ['-i',str(input_path),
                  'turb_driving/stage_consistent=false'],
                  'SGS requires stage_consistent = true')]
        for key,value in [('dedt','0.2'),('kpeak','3'),('rseed','43'),
                          ('stage_consistent','false'),('dt_update','0.02'),
                          ('integrated_work','0'),('new_parameter','1')]:
            diagnostic = ("Parameter 'new_parameter' in block 'turb_driving' on command line not found"
                          if key == 'new_parameter' else
                          'SGS restart mode, coefficients, forcing and scalar layout must match')
            cases.append(('restart_override_'+key,
                ['-r',str(checkpoint),f'turb_driving/{key}={value}'],
                diagnostic))
        for name,args,expected in cases:
            name += suffix
            wd = self.output/name
            wd.mkdir()
            command = [str(self.exe)]+args
            row = dict(name=name,command=command,status='running')
            self.results.append(row)
            self.save()
            started = time.monotonic()
            try:
                with (wd/'run.log').open('w') as log:
                    proc = subprocess.run(command,cwd=wd,stdout=log,
                        stderr=subprocess.STDOUT,timeout=20)
                output = (wd/'run.log').read_text()
                row.update(returncode=proc.returncode,
                    status='pass' if proc.returncode != 0 and expected in output else 'fail',
                    expected_diagnostic=expected)
            except subprocess.TimeoutExpired:
                row.update(status='fail',error='rejection exceeded 20-second cap')
            row['seconds'] = time.monotonic()-started
            self.save()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--exe',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--tlim',type=float,default=1.)
    parser.add_argument('--oracle-only',action='store_true')
    args = parser.parse_args()
    if not 0 < args.tlim <= 2:
        parser.error('initial pilot tlim must be in (0,2]')
    suite = ForcingVerification(args.exe,args.output)
    if args.oracle_only:
        suite.run_forcing('forcing_oracle',True,1.,tlim=.001,oracle=True)
    else:
        for fraction,lane in ((1.,'sol'),(0.,'comp')):
            for enabled in (True,False):
                name = f'{lane}_sgs{int(enabled)}_seed42'
                suite.run_forcing(name,enabled,fraction,tlim=args.tlim)
        suite.run_forcing('forcing_12scalars',True,1.,tlim=.05,nscalars=12)
        suite.forcing_restart('sol_sgs1_restart','sol_sgs1_seed42',True,1.,args.tlim)
        suite.forcing_restart('sol_sgs0_restart','sol_sgs0_seed42',False,1.,args.tlim)
        suite.forcing_restart('comp_sgs1_restart','comp_sgs1_seed42',True,0.,args.tlim)
        suite.reject_contracts('sol_sgs1_seed42')
    summary = {status:sum(row['status']==status for row in suite.results)
               for status in ('pass','fail','blocked')}
    print(json.dumps(summary,indent=2))
    return int(summary['fail'] > 0)


if __name__ == '__main__':
    raise SystemExit(main())
