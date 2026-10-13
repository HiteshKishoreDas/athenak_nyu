# SGS equation fidelity audit evidence

The audit finds no equation, normalization or RK2 source-weight defect in the supported baseline. The native solver and coefficients are unchanged. FIDELITY_AUDIT.md maps equations to source and identifies numerical differences; checks-summary.json and verification.json record the passing checks and source provenance. The newly installed repository test is tst/sgs/verify_specific_diffusion.py on feature commit 470406af. It verifies transport of K/rho with varying density using one-step paired controls.

The authors arXiv manuscript supplies the audited equations. Journal and NASA ADS final-PDF endpoints returned HTTP 403, so the final typeset version was not independently compared. This evidence does not establish published forced-turbulence benchmark reproduction.

The verification subdirectories contain native inputs, checkpoints, histories, logs and detailed results. The preinstallation artifact diffusion run is retained separately; aggregate counts include only verify_specific_diffusion_repo. forcing_oracle is the analytic constant-acceleration control. closure-oracle.log records 22575 passing compiled assertions. fidelity-checks.png and fidelity-checks.pdf show four deterministic numerical checks; they do not depict turbulent closure accuracy.
