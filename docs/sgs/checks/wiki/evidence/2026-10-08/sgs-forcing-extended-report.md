## Bounded extension to t=2

The same frozen binary and all production source hashes match the clean t=1 packet. This fresh extension used the unchanged four 16³, seed42 pilot lanes, three normal intermediate restarts and eight rejection checks. **18 passed, 0 failed, 0 blocked** under 5000-cycle/60-second-per-case caps; aggregate execution time was 61.6 seconds. No source or harness was edited for this extension. Exact commands/input/source hashes and measured checks are in `/private/tmp/sgs-forcing-extended`.

| Lane | t=2 actual work | Final kinetic | Final thermal | Final SGS | Late KE fit drift, t≈1–2 | Fit fractional rate: earlier → extended | Late KE range / mean |
|---|---:|---:|---:|---:|---:|---:|---:|
| sol_sgs1_seed42 | 0.1997975026 | 0.04207958783 | 1.657422124 | 0.0002967907684 | 0.17450 | 0.44176 → 0.17453 | 0.17820 |
| sol_sgs0_seed42 | 0.1994670879 | 0.04217417149 | 1.657292916 | 0 | 0.17551 | 0.44574 → 0.17581 | 0.17887 |
| comp_sgs1_seed42 | 0.1974086899 | 0.01768824726 | 1.679611282 | 0.0001101601999 | 0.12878 | 0.52100 → 0.12883 | 0.91088 |
| comp_sgs0_seed42 | 0.192847412 | 0.01734053015 | 1.675506882 | 0 | 0.13980 | 0.44766 → 0.14014 | 0.91958 |

Fractional fit rate has units of inverse code time and divides the fitted fractional drift by its actual window length. The earlier t=1 window is roughly 0.5–1 (six samples); this t=2 window is roughly 1–2 (eleven samples). Rates decrease, but the fitted t=2 drift is still 12.9–17.6%. Compressive kinetic energy fluctuates strongly: late-window range is 91–92% of the window mean, so the positive fitted slope does not imply monotonic growth. No stationarity acceptance threshold was invented after observing these values.

All applied-work/represented-energy checks pass, with worst residual 6.384e-16; momentum/actual-impulse worst residual 1.055e-17. Both passive scalars are exactly conserved in every physical lane. Three continued/restarted endpoint states, raw OU memory, RNG bytes and integrated work/impulse counters again match exactly. The nominal requested injection remains 0.1, while actually applied work differs between SGS-on/off lanes; a physical equal-work comparison is still unqualified.

Only seed42 was executed. A three-seed campaign, 32³ expansion, long-time independent statistics and a validated turbulence reference remain not run: the M7 reference and stationarity/statistics gates are unmet. The additional run is a bounded local pilot, and scientific M7/M8 acceptance is unchanged.

Extended manifest commit: `14620bb33d580246e32d80f4bd584c6c021bbd1e`; executable SHA-256 `f8b29a46a2e7d6014a1e73227a322dcac6dce056208b25688a5a8b955830c3c5`.
