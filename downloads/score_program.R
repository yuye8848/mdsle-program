#!/usr/bin/env Rscript
## score_program.R - locked-formula scoring of the MD-SLE shared core program
##
## Companion script to the multi-layer evidence matrix
## (Zenodo DOI: 10.5281/zenodo.23680382, version 2).
##
## The formula was locked before validation and is used unchanged across all
## cohorts of the study: ssGSEA (GSVA v2, ssgseaParam) of the 36 core-up genes
## and the 29 core-down genes, scored separately per sample. The cohort-level
## specificity statistic is
##     S = Cliff's delta(core_up) - Cliff's delta(core_down)   (range -2..2)
## with an optional label-permutation p value.
##
## Reference S values from the 14-cohort specificity map (same code path):
##   MD discovery (GSE202657) 2.000 | MD myeloid (GSE244979) 1.917
##   SLE discovery (GSE50772) 1.880 | SLE reference (GSE61635) 1.747
##   SLE (GSE72509) 1.162           | SLE validation (GSE138458) 0.632
##   Neighbour diseases (8 cohorts): -0.105 (MS) .. 0.980 (systemic sclerosis)
##
## Usage:
##   source("score_program.R")
##   sc <- score_program(expr)        # expr: genes x samples matrix,
##                                    # rownames = HGNC symbols, log-scale
##   S_statistic(sc, groups)          # groups: "Case" / "Control" per sample
## Requires: GSVA (>= 2.0). Input should be normalised log-scale expression
## (arrays as deposited; counts as log2(CPM + 1); RPKM as log2(x + 1)).

CORE_UP <- c(
  "B3GNT5", "KLHL2", "NFIL3", "NSUN7", "OSM", "SAT1",
  "STEAP4", "TKT", "TLR2", "TRIOBP", "ACSL1", "BASP1",
  "BCL2A1", "CHMP1B", "CHMP2A", "CLEC4D", "CLIC1", "DNTTIP1",
  "GLIPR2", "IDI1", "JDP2", "KRT23", "MND1", "MNDA",
  "MSRB1", "NAMPT", "NCAPG", "PTEN", "RALB", "RGS2",
  "RILPL2", "RRAGD", "RRP12", "SERPINB1", "TNFAIP3", "VMP1")

CORE_DOWN <- c(
  "AGFG2", "PHC1", "PLXDC1", "ATXN2L", "BCL11A", "BCL9L",
  "CBLL1", "CCDC142", "CD28", "GATAD1", "GPRASP1", "HOOK1",
  "ID3", "LCOR", "NCR3", "NOG", "PPP1R3E", "S1PR1",
  "SGPP1", "SLC6A16", "SS18", "ST6GALNAC1", "TMC8", "TMEM41A",
  "TMEM63A", "TTN", "ZC3H14", "ZFP90", "ZNF793")

score_program <- function(expr) {
  stopifnot(requireNamespace("GSVA", quietly = TRUE))
  p <- GSVA::ssgseaParam(as.matrix(expr),
                         list(core_up = CORE_UP, core_down = CORE_DOWN))
  GSVA::gsva(p)
}

cliffs_delta <- function(x, y) {
  r <- rank(c(x, y)); n1 <- length(x); n2 <- length(y)
  (2 * sum(r[seq_len(n1)]) - n1 * (n1 + n2 + 1)) / (n1 * n2)
}

S_statistic <- function(scores, groups, nperm = 0, seed = 11) {
  up <- scores["core_up", ]; dn <- scores["core_down", ]
  ca <- groups == "Case"; ct <- groups == "Control"
  out <- list(S = cliffs_delta(up[ca], up[ct]) - cliffs_delta(dn[ca], dn[ct]),
              delta_up = cliffs_delta(up[ca], up[ct]),
              delta_down = cliffs_delta(dn[ca], dn[ct]))
  if (nperm > 0) {
    set.seed(seed); cnt <- 0
    for (i in seq_len(nperm)) {
      p <- sample(groups)
      cnt <- cnt + (cliffs_delta(up[p == "Case"], up[p == "Control"]) -
                    cliffs_delta(dn[p == "Case"], dn[p == "Control"]) >= out$S)
    }
    out$perm_p <- (cnt + 1) / (nperm + 1)
  }
  out
}
