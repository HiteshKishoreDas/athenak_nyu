# SGS unresolved-energy spatial check, 2026-10-13

Postprocessing only: existing 10T snapshots of 64³/128³ SGS-on and 256³ SGS-off. No new simulations or solver edits. Reference density and unresolved kinetic energy use the existing conservative finite-volume/Favre restriction. K is energy per volume; sigma is sqrt(2K/rho).

energy-distributions.png compares volume PDFs of log10 K and density-conditioned unresolved velocities (own-run density; volume-cell median and 25–75 percentiles). energy-maps.png compares the same central plane with shared log10 K color limits. metrics.json reports whole-volume statistics, same-cell Pearson correlations, and recovery of reference cells in its highest-K decile by the model highest-K decile. Independent random decile sets would recover 10 percent. All raw input hashes and final mean ratios match the archived campaign.

The 256³ numerical reference is not converged at these filter scales. Late pointwise comparisons combine model differences with independently evolved flow differences; they are not an isolated closure-accuracy test. A similar PDF does not establish correct spatial placement, and initial matching was imposed. Figures describe 10T, when the reference has become mildly subsonic; they do not establish agreement throughout supersonic decay.

Reproduction: run analyze.py using Python with numpy and matplotlib and the existing feature checkout and raw snapshots at the paths in the script. The script uses the saved reference and model fields; it does not launch AthenaK.
