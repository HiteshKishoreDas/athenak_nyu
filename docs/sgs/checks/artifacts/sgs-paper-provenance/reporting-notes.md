# Paper-test scope and provenance

The table reference is the solenoidal64-cubed row in Section4/Table4 of
https://arxiv.org/html/1010.4492 . Forcing setup reference:
https://arxiv.org/pdf/0905.1060 pp3-4.

The32-cubed pilot is below the paper resolution floor. The64-cubed campaign
uses the same realization to4T and a native restart to10T. Ending at4T
truncates the last step; no bitwise identity to an unrun fresh10T is claimed.
All applied counters remain absolute; same-time duplicates must agree before
removing the repeated boundary. Source/executable/input fingerprints are in
individual manifests and the campaign protocol. Solver or forcing changes
are not made between the two64-cubed segments.

Native accel_rms fixes instantaneous physical acceleration magnitude, unlike
published stochastic modal amplitude. Native spectral envelope equivalence
is unproven; Table4 does not explicitly give V/c0 for solenoidal forcing.
V/c0=5.5 is declared, not inferred from table results. These departures prevent
calling this a literal Enzo reproduction or treating table differences as a
unique implementation error. No closure coefficient was tuned.

Original32/64 pilot runner uses uncentered velocity second moment. The final
combined diagnostics reread raw states and center the covariance; the earlier
pilot statistic plots are descriptive. Native Mach/K/epsilon are unchanged.
The original pilot runner is preserved with matching manifest hash.

The first32 pilot command exited1 after numerical PASS because the plot module
was still being written. Its native result remainsPASS and plots were then
regenerated. The first signed-forcing rejection test lacked explicit max_k input
keys; its failedattempt is retained and corrected tests passed21/21.

Original64 full numericalrunPASS; initialreport failed on advanced4.1T-vs4T comparison. Process enumeration was unavailable and no signal was sent. Corrected raw analysisPASS uses exact zero-duration load probe, and full-window statistics retain all absolute counters.
