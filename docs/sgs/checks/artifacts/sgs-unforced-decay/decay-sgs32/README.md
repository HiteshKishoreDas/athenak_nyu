# 32³ unforced SGS decay comparison — 2026-10-13

Completed two new native runs, 32³ SGS on and off, through 10 initial crossing times T. Revalidated/reused the matching 256³ and 128³ SGS-off references; no new reference simulations. Same Gaussian low-mode seed-42 initial field, Mach 9, periodic unit cube, uniform initial density, gamma 1.001, RK2/LLF/PPM4, CFL 0.2, 8 OpenMP threads, no forcing, original SGS coefficients. Near-isothermal ideal EOS and local undeveloped initial conditions differ from the imported developed turbulence in the published decay test. No solver changes or coefficient tuning.

The existing harness run() supports this single-resolution execution; its command-line entrypoint expects two coarse resolutions, so run32.py calls run() with one requested resolution. Fresh outputs are required: run32.py refuses to overwrite the existing native folder. To rerun, copy run32.py to a fresh campaign folder and invoke it with OMP_NUM_THREADS=8 OMP_PROC_BIND=false, using Python with numpy and matplotlib. Then copy analyze32.py alongside it and run that. Feature checkout and reference paths are recorded in the scripts. Source branch codex/multiphase-schmidt-federrath-hd, commit ab54c37deba469b47897966fb511abeccd924947. The executable and all 324 native source file hashes match the retained manifest; source worktree remains clean.

## Results and limits

Over 1–10T, the ratio of integrated mean SGS energy to filtered 256³ unresolved energy is 0.758196. At 10T it is 0.355263; against filtered 128³ it is 0.576079. Model energy starts close to the reference and declines more rapidly relative to it. The late deficit persists against both numerical references. Maximum fixed-filter reference sensitivity over 1–10T is 0.391717, failing the frozen 10% screen. This is not a converged-reference demonstration of absolute model error or full paper reproduction.

Integrated total kinetic energy relative to full 256³ reference: SGS on 0.875674; off 0.860772. Integrated resolved-only kinetic energy relative to filtered reference: on 0.905411; off 1.078661. These metrics give mixed performance, not uniform SGS benefit.

## How to read the plots

energy-and-performance.png: top left directly compares K from removed reference velocities with evolved K; top right shows their ratio and reference resolution sensitivity. Bottom panels compare total and resolved-only kinetic energies, separately, with the appropriate reference. Grey 0–1T startup is excluded from integrated scores; initial K matching was imposed.

budgets-and-reference.png: represented-energy drift, thermal-pressure-work-corrected heating, and the fixed-filter convergence screen. Numerical/RK residual heating is 32.99% of corrected heating; that fraction is not a measured fraction of missing SGS energy.

direct-1T.png, direct-4T.png, direct-10T.png: left shows RMS velocity removed by filtering the fine state to 32³, sqrt(2 K_removed/rho_filtered). Middle shows equivalent speed from the corresponding SGS-on cells, sqrt(2 K_model/rho_model), with identical map colour limits and same central cell plane. Right directly compares the two energy densities in every corresponding coarse cell; red diagonal indicates equality. K_removed = 1/2 average(rho_fine |u_fine-u_Favre|²); its independent kinetic-energy subtraction identity passes at all three times. All six direct input hashes match the archived inventories. Actual snapshot times are retained in summary.json: the largest mismatch is approximately 0.003978T. No spatial interpolation. Late differences include independently evolved flow divergence; 10T reference is mildly subsonic. Scalar SGS energy provides an RMS speed, not fluctuation directions.

report.pdf contains all five figures, 16 panels. 310/310 numerical checks passed, including initial roundtrip, binary/history agreement, represented-energy conservation and SGS transfer accounting. These checks validate numerical bookkeeping, not closure accuracy. Raw native results and full provenance are under native/; summary.json and verification.json retain metrics and post-run checks.
