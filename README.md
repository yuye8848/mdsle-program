# MD–SLE Shared Program — Ménière's Disease × Systemic Lupus Erythematosus

**Live resource: https://mdsle.erbihou.org/**

This repository hosts the companion website and data archive of a study that asked whether Ménière's disease (MD) and systemic lupus erythematosus (SLE) — two diseases not usually studied together — share a transcriptional program in blood. They do: a 123-gene program (65-gene core) that is stable under resampling, localises to the myeloid / dendritic-cell axis in SLE single-cell data, validates with a fixed scoring formula across independent cohorts, and sits at the top of a graded specificity map against neighbouring autoimmune diseases.

梅尼埃病与系统性红斑狼疮共享的血液转录程序：123 个共享基因（核心 65 个）、多层证据矩阵、三法计算扰动裁决与 14 队列特异性图谱，附在线查询与本地打分脚本。网站为中英双语，默认英文。

## Contents

- **The program:** 123 direction-consistent shared genes; core of 65 (36 up-regulated, 29 down-regulated), median resampling inclusion frequency 0.855
- **Evidence matrix:** 123 genes × 29 columns of evidence — independent-cohort differential expression, WGCNA co-expression, single-cell localisation, in-silico perturbation
- **Validation & specificity:** fixed-formula scores over 14 cohorts (MD, SLE, systemic sclerosis, psoriasis, rheumatoid arthritis, primary Sjögren's, Crohn's disease, ulcerative colitis, multiple sclerosis, and acute unilateral vestibulopathy as a vestibular control); MD/SLE cohorts score highest (mean S = 1.556 vs 0.592 for neighbours and the control)
- **Perturbation:** three in-silico methods (scTenifoldKnk, CellOracle, Geneformer) adjudicated under a two-of-three rule nominate **CEBPB** as the consensus regulator of the core up-regulated program (tested in the SLE myeloid compartment)

## Website features

- **Query** — look up any of the 123 genes and see its full evidence card, or any of the 14 perturbation targets and its three-method verdicts
- **Compute** — paste a gene list and check its overlap and direction concordance with the program, computed locally in the browser
- **Download** — the evidence matrix, perturbation-target table, specificity map, and `score_program.R`, the script that reproduces the fixed-formula validation scores on your own expression data

## Data availability

The evidence matrix and perturbation-target table are archived on Zenodo (CC BY 4.0): **DOI [10.5281/zenodo.23681592](https://doi.org/10.5281/zenodo.23681592)** (record v3; matrix schema v2; concept DOI 10.5281/zenodo.23587696). The same files are available from the website's download page.

## Citation

If you use this resource in your work, please cite:

> Lei P. Multi-layer evidence matrix for the 123 shared genes of Ménière's disease and systemic lupus erythematosus. Zenodo, 2026. DOI 10.5281/zenodo.23681592. Website: https://mdsle.erbihou.org/

Companion manuscript: Lei P, Dong W, Tan J, Yang X. *A Shared Blood Immune Transcriptional Program in Ménière's Disease and Systemic Lupus Erythematosus: Computational Perturbation and Cross-Cohort Specificity Validation.* Submitted to PLOS ONE (2026).

## Related resources

- [PolypAtlas](https://polyp.erbihou.org/) — curated multi-omics resource for chronic rhinosinusitis with nasal polyps (CRSwNP)
- [RhinitisAtlas](https://rhinitis.erbihou.org/) — unified public-transcriptome resource for allergic rhinitis

## Repository layout

- This repository is the static website itself, served via GitHub Pages at the address above.

## Contact

Peiliang Lei (of@live.com), Department of Otolaryngology, The Third Hospital of Shijiazhuang, Shijiazhuang, Hebei, China.
