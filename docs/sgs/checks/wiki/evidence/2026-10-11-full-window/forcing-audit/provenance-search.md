# Bounded search for the original Schmidt–Federrath Enzo setup

Checked 2026-10-11. Scope: public primary sources, author repositories,
paper metadata, and Enzo archive metadata. No simulations, native code changes,
new dependencies, coefficient tuning, or messages to authors.

## Result

This search did **not** recover an identifiable source revision, input file,
forcing realization, or numerical reference dataset for the Enzo 1.5
hydrodynamic LES reported by Schmidt & Federrath (2011), especially Table 4.
The sources found establish related conventions, but do not establish exact
forcing equivalence. This is a bounded negative search result, not proof that
the original material does not exist.

The strongest newly found provenance evidence is the official Enzo release
history: Enzo 2.5, released in March 2016, lists stochastic forcing as a new
feature. The public forcing file history begins in October 2014. Consequently,
the modern public implementation cannot simply be assigned to the 2011 runs.
[Official release notes](https://enzo-project.org/ReleaseNotes.html),
[October 2014 introduction commit](https://github.com/enzo-project/enzo-dev/commit/d368b31235441402748ec7bf9ac819b6a4af7040).

## What remains identifiable

| Quantity | Evidence and status |
| --- | --- |
| Original solver version | Enzo 1.5, stated in SF §4. |
| Integral length and time | Half the box size; force correlation time is T=L/V. Confirmed in SF §4. |
| Pure solenoidal mixture | ζ=1. Confirmed. |
| Pure solenoidal Table 4 V/c0 | Not explicitly specified by the inspected SF text. Table 5 is a different mixed forcing experiment. |
| Original absolute OU forcing amplitude | Not recovered. Schmidt et al. (2009) defines a related normalization, but no original Table 4 input or source revision connects it unambiguously to these runs. |
| Original radial forcing envelope/bandwidth | Not recovered from a run input or the original implementation. Later code and related simulations are insufficient provenance. |
| Gaussian comparison filter length | SF §4.2 gives ΔG=16Δ=L/32 for its 1024³ ILES versus 64³ LES comparison. |
| Exact Gaussian kernel/standard deviation | SF text inspected here does not write the kernel or convert ΔG to σ. No original filtering script was recovered. Do not equate ΔG with σ without an additional convention. |
| Original reference fields or time series | Not recovered. Published figures and table values are available; modern Enzo test reference files are a different experiment. |

These paper statements are from the [author-submitted SF manuscript, §§4 and
4.2](https://arxiv.org/html/1010.4492). Its [arXiv record](https://arxiv.org/abs/1010.4492)
lists only v1 (October 2010). The publisher HTML endpoint returned 403 in this
check, so this search does not independently verify differences between the
submitted manuscript and the final journal text.

The three meanings of “Gaussian width” must remain separate: the spatial
comparison filter, a radial forcing bandwidth, and the variance of Gaussian
OU innovations. Recovering one does not recover the other two.

## Exact discovery queries and their results

Search queries are reproduced verbatim. Results were followed only when they
led to primary material relevant to this provenance question.

1. `"Schmidt" "Federrath" "2011" "Enzo" "code" "forcing"`
2. `"Schmidt" "Enzo 1.5" turbulence`
3. `"1010.4492" "data"`
4. `"StochasticForcing" "Schmidt" "2005" Enzo`

This batch found current Enzo parameter/test documentation and related
turbulence papers. No original SF run archive appeared. Current
[test documentation](https://enzo.readthedocs.io/en/latest/testproblems/index.html)
attributes the stochastic forcing test to Philipp Grete (2014) and describes
the later SGS model. It is not an identifiable SF2011 Table 4 input.

5. `"Schmidt" "Federrath" "2011" turbulence "supplementary"`
6. `"Enzo 1.5" "Schmidt" source code`
7. `"Schmidt Federrath" "Gaussian" "filter" 2011`
8. `"Wolfram Schmidt" "Enzo" code download`

This batch found Schmidt's official Hamburg group page and related publications.
The group page links a public Enzo fork explicitly associated with modern
disk galaxies and an MHD SGS model; it provides no link identified as the
2011 hydrodynamic turbulence archive.
[Author group page](https://www.sternwarte.uni-hamburg.de/en/hs/group---schmidt.html).

9. `Enzo "1.5" download source archive lca ucsd`
10. `"Schmidt" "Federrath" "A106" "data" turbulence`
11. `site:ita.uni-heidelberg.de/~chfeder "Schmidt" "2011" "SGS"`
12. `"FEARLESS" "code" "download" Enzo Schmidt`

This batch found the [ANU institutional paper record](https://researchportalplus.anu.edu.au/en/publications/a-fluid-dynamical-subgrid-scale-model-for-highly-compressible-ast/)
and an author conference talk. The institutional record links the article DOI;
no original source/input/data supplement was identified in its returned record.
Unrelated FEARLESS products and publications were discarded.

13. `enzo 1.5 source download`
    — domains: `enzo-project.org`, `github.com`, `lca.ucsd.edu`
14. `Schmidt Federrath 2011 A106 SGS data`
    — domains: `mso.anu.edu.au`, `ita.uni-heidelberg.de`, `zenodo.org`
15. `Schmidt 2006 localised subgrid scale Gaussian filter exp delta`
    — domain: `arxiv.org`

This final batch found the official Enzo project/release pages, modern
repositories, and the cited 2006 theory paper. It produced no identifiable
original SF archive. Search-engine results are discovery evidence only, not
evidence for numerical parameter values.

## Primary pages and repository checks

- [SF arXiv abstract](https://arxiv.org/abs/1010.4492): PDF, experimental HTML,
  and TeX-source links; no explicit ancillary-data entry was present in the
  returned page. Generic arXiv code/data service toggles are not an author
  dataset or original run input. The TeX archive was not extracted.
- [Publisher HTML](https://www.aanda.org/articles/aa/full_html/2011/04/aa15630-10/aa15630-10.html):
  access failed with 403. No supplementary archive can be excluded behind
  that endpoint.
- The LCA site named by the paper, `http://lca.ucsd.edu/`, was inaccessible
  through the browser tool. Historical web snapshots and private author
  storage were not searched.
- [Federrath's official homepage](https://www.mso.anu.edu.au/~chfeder/) links
  [TurbGen](https://www.mso.anu.edu.au/~chfeder/codes/turbgen/turbgen.html),
  GitHub and Bitbucket. The code page describes a portable OU generator and
  cites the 2010 forcing paper and a 2022 software record. It supplies no
  identified Enzo 1.5 Table 4 archive. GitHub API reports the
  `chfeder/turbulence_generator` repository creation as 2022-04-25; this is
  the hosting date, not proof of the source algorithm's age.
- The author's Bitbucket repository listing returned no readable entries.
  The historical Enzo stochastic-forcing pull request linked by the release
  notes also failed to load. Neither failure is evidence of absence.

The following GitHub API endpoints were read directly, without cloning,
building, or saving foreign source trees:

| Exact endpoint | Result |
| --- | --- |
| `https://api.github.com/repos/enzo-project/enzo-dev/tags?per_page=100` | 32 tags; release tags run from Enzo 2.0 to 2.6.1. No tag identified as Enzo 1.5. |
| `https://api.github.com/repos/enzo-project/enzo-dev/tags?per_page=100&page=2` | Empty list, completing this tag inventory. |
| `https://api.github.com/users/wolfram-schmidt/repos?per_page=100` | Seven public repositories: amrex, athenapk, enzo-dev, enzo-e, isolated-galaxy, masterfiles, Nyx. No separately named original SF2011 repository. |
| `https://api.github.com/repos/wolfram-schmidt/enzo-dev/branches?per_page=100` | cloud-wind, diskgalaxy, diskgalaxy-sgs, main, stable. No branch explicitly identified as the 2011 Enzo 1.5 run source. |
| `https://api.github.com/repos/wolfram-schmidt/enzo-dev/commits?path=src/enzo/StochasticForcing_Init.C&per_page=100` | Six entries; oldest is d368b31235441402748ec7bf9ac819b6a4af7040, 2014-10-15, “Actually committed the StochasticForcing files...”. This is the history of that path on the default branch, not a search of every renamed file or detached historical commit. |
| `https://api.github.com/users/chfeder/repos?per_page=100` | Eight public repositories, including turbulence_generator and cfflash; none explicitly named as the original SF2011 Enzo source/input/data archive. |
| `https://api.github.com/repos/enzo-project/enzo-dev/git/trees/df9c286b47a49806f494222c579242550601d642?recursive=1` | Complete Enzo 2.0 release tree (`truncated=false`). Case-insensitive path search for `stochastic`, `drivenflow`, `schmidt`, `federrath`, `gaussian`, `fearless` found no paths. This is a filename search, not a whole-tree content audit. |
| `https://api.github.com/repos/wolfram-schmidt/enzo-dev/git/trees/839511f7a2d1f601cc5994eff8bca526289c043a?recursive=1` | Complete diskgalaxy-sgs branch tree. Contains StochasticForcing implementation and `run/MHD/3D/StochasticForcing` input plus Mach/time/vorticity reference files. No provenance connecting these modern files to Table 4 was found. |

The inspected modern main tree is pinned at
`18812cdde8851e07592733bda26f7575554ba52c`. Its forcing files describe
adaptation to Enzo 2.4 in 2014; their older author headers do not by themselves
establish identical code, parameters, or normalization in 2011.
[Pinned modern initialization source](https://github.com/enzo-project/enzo-dev/blob/18812cdde8851e07592733bda26f7575554ba52c/src/enzo/StochasticForcing_Init.C).

## Consequence for the present validation

Keep the ongoing forcing protocol frozen and finish the prescribed 2T–10T
averages. Compare with the paper as a **related setup with unestablished exact
forcing/filter equivalence**. Do not set acceleration RMS to 6, change c0,
or change the radial envelope merely to move Table 4 residuals toward zero.

The existing forcing-only oracle measures a real statistical difference
between instantaneous RMS normalization and modal OU normalization. It does
not prove which implementation generated Table 4, or identify an original
absolute target. The related Schmidt2009 normalization remains a conditional
reference, not an established SF2011 input.

Closing the original provenance gap requires an identifiable original Enzo
1.5 patch/revision and input, original filtering script/kernel definition,
or explicit clarification with run provenance from the authors. No request
to the authors was sent. Stop this search here; further generic searches or
tuning against published aggregate values would not establish equivalence.
