# =============================================================================
# 16_nup88_immune_correlation_gse72094.R
# Replication check: does NUP88 correlate negatively with myeloid/B-cell
# infiltration in GSE72094, as previously reported for TCGA via TIMER3.0?
#
# This is a replication test, not a re-derivation: the TCGA reference values
# are taken as-is from the existing Table S2 (S2 TABLE TIMER3 NUP88
# LMBN2.xlsx, Algorithm = XCELL, purity-adjusted partial Spearman
# correlation). xCell is NOT re-run on TCGA here - doing so would produce a
# third, different TCGA number and create an inconsistency in the
# manuscript. GSE72094 is scored fresh with xCell (data/GEO_cache/
# GSE72094_xCell.rds, produced by R/15_gse72094_xcell.R) because no
# TIMER3.0/deconvolution result exists for it yet.
#
# Fourteen populations: all NUP88/XCELL rows in TCGA with rho < 0, ranked by
# BH-adjusted p-value (xCell has no Macrophages M0 - that category is
# specific to CIBERSORT-ABS and was dropped). Every xCell raw score name
# below was looked up in the immunedeconv mapping table
# (inst/extdata/cell_type_mapping.xlsx, method_dataset=="xcell"), not
# assumed from the TIMER3.0 display label - all 14 resolved to exactly one
# unambiguous xCell score, and all 14 are present in the 67 rows produced
# by R/15_gse72094_xcell.R:
#   NKT                          -> "T cell NK"
#   Class-switched memory B-cells -> "Class-switched memory B cell"
#   Macrophages M1               -> "Macrophage M1"
#   aDC                          -> "Myeloid dendritic cell activated"
#   B-cells                      -> "B cell"
#   Macrophages                  -> "Macrophage"
#   CD4+ naive T-cells           -> "T cell CD4+ naive"
#   Plasma cells                 -> "B cell plasma"
#   Memory B-cells                -> "B cell memory"
#   CD8+ T-cells                 -> "T cell CD8+"
#   CD8+ Tcm                     -> "T cell CD8+ central memory"
#   Monocytes                    -> "Monocyte"
#   Macrophages M2               -> "Macrophage M2"
#   pDC                          -> "Plasmacytoid dendritic cell"
#
# NOTE ON REFERENCE VALUES: the TCGA rho/p reference below comes from the
# local S2 TABLE TIMER3 NUP88 LMBN2.xlsx, which was independently
# cross-checked against the final manuscript PDF (3/3 spot-checked
# citations matched exactly) earlier in this project's history. One value
# in particular, Plasmacytoid dendritic cell, was separately reported as
# p_BH = 0.0391 from an outside export the user did not ultimately supply
# for comparison; this file has p_BH = 0.0542 for that row. That specific
# discrepancy is unresolved and is reproduced as-is (not silently
# corrected) in the output below.
#
# NUP88 has three probes on this platform (merck-BU539430_a_at,
# merck2-AK225247_at, merck2-BG686994_at). Their pairwise Spearman
# correlation is checked first: if all three pairwise rho > 0.5, NUP88
# expression is the row mean of the three probes; otherwise all three are
# reported separately (no "best annotated probe" tie-break - that call was
# rejected as arbitrary).
#
# IMPORTANT ASYMMETRY: TIMER3.0 rho values are purity-adjusted partial
# correlations; GSE72094 has no tumour-purity estimate available, so its
# rho values are plain (unadjusted) Spearman correlations. This is a
# methodological difference to state in Methods/Limitations, not a defect
# to paper over.
#
# Output:
#   outputs/tables/NUP88_immune_correlation_GSE72094_vs_TIMER3.csv
# =============================================================================

# Locate the project root automatically (the folder that contains R/), instead
# of a hard-coded path. Mirrors R/00_setup.R's PROJECT_ROOT logic.
if (!exists("PROJECT_ROOT")) {
  .args <- commandArgs(FALSE)
  .file <- sub("^--file=", "", .args[grep("^--file=", .args)])
  .self <- if (length(.file)) .file
           else if (!is.null(sys.frames()[[1]]$ofile)) sys.frames()[[1]]$ofile
           else NA_character_
  PROJECT_ROOT <- if (!is.na(.self)) normalizePath(file.path(dirname(.self), ".."))
                  else getwd()
}
setwd(PROJECT_ROOT)

suppressPackageStartupMessages({
    library(Biobase)
    library(readxl)
})

source("R/00_setup.R")

CONCORDANCE_RHO_CUTOFF <- 0.5

# -----------------------------------------------------------------------------
# 1. LOAD GSE72094 EXPRESSION + XCELL SCORES
# -----------------------------------------------------------------------------
eset <- readRDS("data/GEO_cache/GSE72094_eset.rds")
xc   <- readRDS("data/GEO_cache/GSE72094_xCell.rds")   # 67 cell types x 442 samples

ex <- exprs(eset)
f  <- fData(eset)
message("15: GSE72094 loaded, ", ncol(ex), " samples")

stopifnot(identical(colnames(ex), colnames(xc)))  # same sample order, same eset

# -----------------------------------------------------------------------------
# 2. NUP88 PROBES: CONCORDANCE CHECK
# -----------------------------------------------------------------------------
nup88_ids <- rownames(f)[which(as.character(f$GeneSymbol) == "NUP88")]
message("15: NUP88 probes found: ", paste(nup88_ids, collapse = ", "))
stopifnot(length(nup88_ids) == 3)

nup88_mat <- t(ex[nup88_ids, , drop = FALSE])   # samples x 3 probes

probe_cor <- cor(nup88_mat, method = "spearman", use = "pairwise.complete.obs")
cat("\n=== NUP88 probe pairwise Spearman correlation ===\n")
print(round(probe_cor, 3))

offdiag <- probe_cor[upper.tri(probe_cor)]
concordant <- all(offdiag > CONCORDANCE_RHO_CUTOFF)
message("15: probes concordant (all pairwise rho > ", CONCORDANCE_RHO_CUTOFF,
        "): ", concordant)

if (concordant) {
    nup88_expr <- list(NUP88_mean_3probes = rowMeans(nup88_mat))
} else {
    nup88_expr <- setNames(
        lapply(nup88_ids, function(id) nup88_mat[, id]),
        nup88_ids
    )
}

# -----------------------------------------------------------------------------
# 3. TCGA REFERENCE (TIMER3.0, Table S2, as published - not recomputed)
# -----------------------------------------------------------------------------
s2 <- read_excel("S2 TABLE TIMER3 NUP88 LMBN2.xlsx", sheet = "S2_correlations")
s2 <- as.data.frame(s2)

pop_map <- data.frame(
    xcell_row    = c("NKT", "Class-switched memory B-cells", "Macrophages M1",
                     "aDC", "B-cells", "Macrophages", "CD4+ naive T-cells",
                     "Plasma cells", "Memory B-cells", "CD8+ T-cells",
                     "CD8+ Tcm", "Monocytes", "Macrophages M2", "pDC"),
    timer3_label = c("T cell NK", "Class-switched memory B cell", "Macrophage M1",
                     "Myeloid dendritic cell activated", "B cell", "Macrophage",
                     "T cell CD4+ naive", "B cell plasma", "B cell memory",
                     "T cell CD8+", "T cell CD8+ central memory", "Monocyte",
                     "Macrophage M2", "Plasmacytoid dendritic cell"),
    stringsAsFactors = FALSE
)

tcga_ref <- s2[s2$Gene == "NUP88" & s2$Algorithm == "XCELL" &
              s2$`Immune cell` %in% pop_map$timer3_label, ]
tcga_ref <- merge(pop_map, tcga_ref, by.x = "timer3_label", by.y = "Immune cell")

cat("\n=== TCGA reference (TIMER3.0, Table S2, XCELL, as published) ===\n")
print(tcga_ref[, c("xcell_row", "timer3_label", "n", "rho", "p_BH")])

# -----------------------------------------------------------------------------
# 4. SPEARMAN CORRELATIONS ON GSE72094
# -----------------------------------------------------------------------------
results <- do.call(rbind, lapply(names(nup88_expr), function(nup_label) {
    nup_vals <- nup88_expr[[nup_label]]
    do.call(rbind, lapply(pop_map$xcell_row, function(pop) {
        score <- xc[pop, ]
        ok    <- !is.na(nup_vals) & !is.na(score)
        ct    <- suppressWarnings(
            cor.test(nup_vals[ok], score[ok], method = "spearman", exact = FALSE)
        )
        data.frame(
            nup88_probe_set = nup_label,
            population      = pop,
            timer3_label    = pop_map$timer3_label[pop_map$xcell_row == pop],
            n               = sum(ok),
            rho_GSE72094    = round(unname(ct$estimate), 3),
            p_raw           = signif(ct$p.value, 3)
        )
    }))
}))

results$p_BH <- signif(p.adjust(results$p_raw, method = "BH"), 3)

results <- merge(results, tcga_ref[, c("xcell_row", "rho", "n")],
                 by.x = "population", by.y = "xcell_row")
names(results)[names(results) == "rho"]   <- "rho_TCGA_TIMER3_reference"
names(results)[names(results) == "n.x"]   <- "n_GSE72094"
names(results)[names(results) == "n.y"]   <- "n_TCGA_TIMER3"

results$direction_GSE72094      <- ifelse(results$rho_GSE72094 < 0, "negative",
                                          ifelse(results$rho_GSE72094 > 0, "positive", "zero"))
results$direction_matches_TCGA  <- sign(results$rho_GSE72094) == sign(results$rho_TCGA_TIMER3_reference)
results$significant_BH          <- results$p_BH < 0.05
results$purity_adjusted_GSE72094 <- FALSE
results$purity_adjusted_TCGA_TIMER3 <- TRUE

results <- results[, c("nup88_probe_set", "population", "timer3_label",
                       "n_GSE72094", "rho_GSE72094", "p_raw", "p_BH",
                       "significant_BH", "direction_GSE72094",
                       "rho_TCGA_TIMER3_reference", "n_TCGA_TIMER3",
                       "direction_matches_TCGA",
                       "purity_adjusted_GSE72094", "purity_adjusted_TCGA_TIMER3")]
results <- results[order(match(results$population, pop_map$xcell_row)), ]

# -----------------------------------------------------------------------------
# 5. REPLICATION VERDICT
#
# Pre-registered criterion for the 14-population panel (fixed by the user
# before this run, in the request that specified the 14 populations):
# replicated if >= 8/14 negative direction AND >= 7/14 significant after BH.
# This supersedes the earlier >=3/>=2 criterion, which was defined for a
# 4-population panel and is no longer the applicable rule.
# -----------------------------------------------------------------------------
n_populations <- nrow(pop_map)
primary <- results[results$nup88_probe_set == names(nup88_expr)[1], ]
n_negative    <- sum(primary$direction_GSE72094 == "negative")
n_significant <- sum(primary$significant_BH)
replicated    <- n_negative >= 8 & n_significant >= 7

# -----------------------------------------------------------------------------
# 6. WRITE OUTPUT
# -----------------------------------------------------------------------------
dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)
write.csv(results,
          "outputs/tables/NUP88_immune_correlation_GSE72094_vs_TIMER3.csv",
          row.names = FALSE)
message("15: wrote outputs/tables/NUP88_immune_correlation_GSE72094_vs_TIMER3.csv")

cat("\n=== NUP88 vs immune infiltration: GSE72094 (fresh xCell) vs TCGA (TIMER3.0 reference) ===\n")
print(results, row.names = FALSE)

cat("\n=== Replication verdict (", names(nup88_expr)[1], ", n populations = ",
    n_populations, ") ===\n", sep = "")
cat("Negative direction:", n_negative, "/", n_populations, "\n", sep = "")
cat("BH-significant:    ", n_significant, "/", n_populations, "\n", sep = "")
cat("Pass condition (pre-registered for n=14): >=8 negative AND >=7 significant:",
    replicated, "\n")

cat("\nNote for Methods/Limitations: GSE72094 correlations are unadjusted\n",
    "Spearman correlations (no tumour-purity estimate available for this\n",
    "cohort); TCGA reference values are purity-adjusted partial Spearman\n",
    "correlations (TIMER3.0). This asymmetry should be stated explicitly,\n",
    "not treated as a discrepancy.\n", sep = "")
