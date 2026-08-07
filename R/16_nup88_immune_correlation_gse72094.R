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
# Four populations only (xCell has no Macrophages M0 - that category is
# specific to CIBERSORT-ABS and was dropped):
#   aDC             -> TIMER3.0 label "Myeloid dendritic cell activated"
#                      (confirmed via immunedeconv source, not by intuition:
#                      inst/extdata/cell_type_mapping.xlsx maps xCell's aDC,
#                      not cDC or generic DC, to this label)
#   Macrophages     -> "Macrophage"
#   Macrophages M1  -> "Macrophage M1"
#   B-cells         -> "B cell"
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

setwd("~/Desktop/luad-pipeline")

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
    xcell_row       = c("aDC", "Macrophages", "Macrophages M1", "B-cells"),
    timer3_label    = c("Myeloid dendritic cell activated", "Macrophage",
                        "Macrophage M1", "B cell"),
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
# 5. REPLICATION VERDICT (fixed criterion: >=3/4 negative direction,
#    >=2/4 significant after BH, evaluated on the primary probe-set only)
# -----------------------------------------------------------------------------
primary <- results[results$nup88_probe_set == names(nup88_expr)[1], ]
n_negative    <- sum(primary$direction_GSE72094 == "negative")
n_significant <- sum(primary$significant_BH)
replicated    <- n_negative >= 3 & n_significant >= 2

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

cat("\n=== Replication verdict (", names(nup88_expr)[1], ", n populations = 4) ===\n", sep = "")
cat("Negative direction:", n_negative, "/4\n")
cat("BH-significant:    ", n_significant, "/4\n")
cat("Replicated (>=3/4 negative AND >=2/4 significant):", replicated, "\n")

cat("\nNote for Methods/Limitations: GSE72094 correlations are unadjusted\n",
    "Spearman correlations (no tumour-purity estimate available for this\n",
    "cohort); TCGA reference values are purity-adjusted partial Spearman\n",
    "correlations (TIMER3.0). This asymmetry should be stated explicitly,\n",
    "not treated as a discrepancy.\n", sep = "")
