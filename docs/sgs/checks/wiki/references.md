---
title: AthenaK subgrid turbulence — literature and code landscape
class: append-only
created: 2026-10-06
---

# Literature and code landscape

## 2026-10-06 — exploratory digest

All entries were retrieved this session on 2026-10-06 from primary research,
publisher or upstream code pages. Reading depth is explicit. This is a targeted
landscape, not an exhaustive prior-art sweep or full bibliography audit.

| Reference | Connection | Reading depth |
|---|---|---|
| W. Schmidt & C. Federrath, [A Fluid-Dynamical Subgrid Scale Model for Highly Compressible Astrophysical Turbulence](https://arxiv.org/abs/1010.4492), preprint 2010; A&A 528, A106 (2011), DOI 10.1051/0004-6361/201015630 | Compressible SGS energy model; mixed closure; filtered-box and evolved-model tests | arXiv metadata/abstract and [full-text passages](https://arxiv.org/html/1010.4492) on filtering, energy, effective pressure and mixed closure. Publication context also cross-checked in the 2025 review. |
| W. Schmidt, A. S. Almgren, H. Braun, J. F. Engels, J. C. Niemeyer, J. Schulz, R. R. Mekuria, A. J. Aspden & J. B. Bell, [Cosmological Fluid Mechanics with Adaptively Refined Large Eddy Simulations](https://arxiv.org/abs/1309.3996), preprint 2013; MNRAS 440, 3051–3077 (2014) | Conservative resolved/SGS compensation across AMR levels and forced-turbulence tests | Metadata/abstract; [Section 2 and Appendix A.2.2 passages](https://arxiv.org/html/1309.3996) on restriction, refinement, bulk-flow contamination and positivity. Publisher search record supplied journal metadata; full publisher retrieval failed. |
| D. G. Vlaykov, P. Grete, W. Schmidt & D. R. G. Schleicher (2016), [A nonlinear structural subgrid-scale closure for compressible MHD Part I: derivation and energy dissipation properties](https://arxiv.org/abs/1606.02064), Physics of Plasmas 23, 062316 | Compressible MHD stresses, electromotive force and energy transfers | arXiv metadata/abstract; derivation not independently audited. |
| P. Grete, D. G. Vlaykov, W. Schmidt & D. R. G. Schleicher (2016), [A nonlinear structural subgrid-scale closure for compressible MHD Part II: a priori comparison on turbulence simulation data](https://arxiv.org/abs/1606.01573), Physics of Plasmas 23, 062317 | Filtered-reference validation across flow conditions | arXiv metadata/abstract; no reproduction of their tests or correlation measurements. |
| P. Grete, D. G. Vlaykov, W. Schmidt & D. R. G. Schleicher (2017), [Comparative statistics of selected subgrid-scale models in large eddy simulations of decaying, supersonic MHD turbulence](https://arxiv.org/abs/1703.00858), Physical Review E 95, 033206 | A posteriori comparison: numerical damping and explicit filter affect whether closure changes statistics | arXiv abstract via search; [publisher metadata/abstract](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.95.033206) via search. Full text not inspected. |
| W. Schmidt-Brückner (2025), [Large eddy simulations in astrophysics](https://doi.org/10.1007/s41115-025-00023-1), Living Reviews in Computational Astrophysics 11, 2 | Updated map of explicit/implicit LES, HD and MHD models, mixing and numerical-diffusion limitations | Publisher metadata/abstract and selected passages on filtering, adaptive energy exchange and applications. This is the major 2025 revision, not the original 2015 version. |

### AthenaK infrastructure inspected

- [Upstream README](https://github.com/IAS-Astrophysics/athenak): Kokkos,
  block-based AMR and HD/MHD infrastructure.
- [Source directory](https://github.com/IAS-Astrophysics/athenak/tree/main/src)
  and [source terms](https://github.com/IAS-Astrophysics/athenak/tree/main/src/srcterms):
  module listings; not a repository-wide source audit.
- [Turbulence driver](https://github.com/IAS-Astrophysics/athenak/blob/main/src/srcterms/turb_driver.cpp):
  selected mode initialization, task integration and forcing-update passages.
- [Fluid problem generators](https://github.com/IAS-Astrophysics/athenak/tree/main/src/pgen/fluids)
  and [hydro inputs](https://github.com/IAS-Astrophysics/athenak/tree/main/inputs/hydro):
  listings show `turb.cpp` and `turb.athinput`. The attempted input-file retrieval
  failed; its contents were not inspected.
- [Hydrodynamics documentation](https://ias-astrophysics.github.io/athenak-docs/hydro.html)
  and [AMR documentation](https://ias-astrophysics.github.io/athenak-docs/amr.html):
  pages retrieved, without an end-to-end verification against a pinned code revision.

GitHub pages are mutable upstream snapshots and were not pinned to a commit.
Selected source passages confirm starting infrastructure; they do not establish
correct forcing normalization through regridding or absence of an SGS implementation.
No local AthenaK checkout was inspected or changed.

### Search coverage

- `Schmidt Federrath subgrid scale model highly compressible astrophysical turbulence 2011`
- `Schmidt adaptive mesh refinement large eddy simulations turbulence energy conservation 2014`
- `AthenaK github turbulence forcing hydro diffusion viscosity`
- `"AthenaK" "subgrid" turbulence`
- `"AthenaK" "SGS"`
- `"AthenaK" "subgrid turbulence"`
- `Schmidt Niemeyer Hillebrandt localized subgrid scale model turbulence 2006 I 0601499`
- `"A nonlinear structural subgrid-scale closure for compressible MHD" arxiv Grete Vlaykov`
- `"Comparative statistics of selected subgrid-scale models" arxiv`

General web search supplied leads; only primary pages support the subject's
scientific claims. NASA ADS, SciSpace and Semantic Scholar were not independently
searched. No full citation-graph, branch, fork or unpublished-work search was done.
No exact AthenaK implementation was established by these queries, which does not
prove one is absent. Cooling/multiphase and observational extensions remain
proposed connections; their exact prior-art coverage was not established.

## 2026-10-08 — local AthenaK inspection and implementation priority

The operator chose to implement the Schmidt–Federrath method before extensions.
The [session record](sessions/2026-10-08-implementation-first.md) identifies the
published 2011 HD model as the baseline and separates that decision from the
assistant's proposed uniform-grid validation sequence.

Read-only local inspection identified `/Users/hitesh/hitesh/git/athenak_nyu`,
branch `gotham_development`, HEAD `07c529202befbeb6ef5e694d95286929618ff7c1`.
Selected hydro, EOS, diffusion/source, turbulence, output and build/test files were
read. Existing scalar transport and boundary infrastructure suggest an integration
path; scalar positivity and resolved/SGS energy coupling need explicit handling.
The first passive scalar also has an existing metallicity role in cooling code.

This supersedes the earlier statement that no local checkout had been inspected.
It does not establish absence of SGS support across all AthenaK forks or branches.
The operator's target checkout remains unconfirmed. No code changes, build,
simulation, performance measurement or new literature verification were performed.

## 2026-10-08 — agentic development planning: primary text and code recheck

Reopened [Schmidt & Federrath (2011), §§2–4](https://arxiv.org/html/1010.4492).
Checked stress normalization/sign, the trace-pressure relation, production,
dissipation and SGS-specific-energy transport, paired resolved/SGS equations,
coefficient values and the energy-evolving pseudo-isothermal test setup. These
support the equation contract in [development-plan.md](development-plan.md);
discrete algorithms, tolerance targets and test designs are proposals, not
algorithms independently reproduced from the paper.

Reopened [Schmidt et al. (2014), §2 and Appendix A.2](https://arxiv.org/html/1309.3996)
on conservative inter-level compensation and positivity limitations. This supports
the later AMR work boundary; no transfer rule was selected or implemented.

Rechecked the candidate local commit `07c529202befbeb6ef5e694d95286929618ff7c1`.
Inspected hydro task/fallback logic, restart counts, output configuration, global
timestep selection in `src/mesh/mesh.cpp`, turbulence-driving settings and the
existing `tst/scripts/hydro/hydro_linwave.py` harness example. The code map remains
specific to that candidate revision; target selection and all execution coverage
remain pending. No solver build, numerical test or performance run was performed.
