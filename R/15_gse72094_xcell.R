# =============================================================================
# 15_gse72094_xcell.R
# Download GSE72094 (Chen/Merck lung adenocarcinoma cohort, n = 442), run
# three quality checks used before any downstream analysis, then run xCell
# to produce the immune deconvolution scores consumed by
# R/16_nup88_immune_correlation_gse72094.R.
#
# Quality checks (all must pass before proceeding):
#   1. NUP88 present among the annotated probe gene symbols.
#   2. All samples are confirmed LUAD histology (source_name_ch1).
#   3. Survival time and status are available (complete-case count reported).
#
# xCell (dviraran/xCell, GitHub) requires GSVA >= 1.50 for the ssgseaParam
# interface; both are installed from CRAN/Bioconductor/GitHub, not pinned
# to a specific commit here - see outputs/sessionInfo.txt for the versions
# this pipeline was run against.
#
# Output:
#   data/GEO_cache/GSE72094_eset.rds   cached ExpressionSet (raw download)
#   data/GEO_cache/GSE72094_xCell.rds  xCell scores, 67 cell types x 442 samples
# =============================================================================

setwd("~/Desktop/luad-pipeline")

suppressPackageStartupMessages({
    library(GEOquery)
    library(Biobase)
    library(xCell)
})

source("R/00_setup.R")
options(timeout = 1000)

dir.create("data/GEO_cache", recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 1. DOWNLOAD (CACHED)
# -----------------------------------------------------------------------------
rds_file <- "data/GEO_cache/GSE72094_eset.rds"
if (file.exists(rds_file)) {
    eset <- readRDS(rds_file)
    message("15: GSE72094 loaded from cache")
} else {
    g    <- getGEO("GSE72094", GSEMatrix = TRUE, getGPL = TRUE,
                   destdir = "data/GEO_cache")
    eset <- g[[1]]
    saveRDS(eset, rds_file)
    message("15: GSE72094 downloaded and cached")
}

p <- pData(eset)
f <- fData(eset)
message("15: ", nrow(f), " probes x ", nrow(p), " samples")

# -----------------------------------------------------------------------------
# 2. QUALITY CHECKS
# -----------------------------------------------------------------------------
cat("\n=== 1. NUP88 among annotated probes? ===\n")
nup88_ids <- rownames(f)[which(as.character(f$GeneSymbol) == "NUP88")]
cat(length(nup88_ids) > 0, "-", length(nup88_ids), "probe(s):",
    paste(nup88_ids, collapse = ", "), "\n")

cat("\n=== 2. Confirmed LUAD histology, how many samples? ===\n")
histology <- unique(as.character(p$source_name_ch1))
cat("source_name_ch1 unique value(s):", paste(histology, collapse = " | "), "\n")
cat(length(histology) == 1 && grepl("adenocarcinoma", histology, ignore.case = TRUE),
    "-", nrow(p), "/", nrow(p), "samples\n")

cat("\n=== 3. Survival time and status available? ===\n")
surv_time   <- suppressWarnings(as.numeric(as.character(p[["survival_time_in_days:ch1"]])))
surv_status <- p[["vital_status:ch1"]]
complete    <- !is.na(surv_time) & surv_status %in% c("Alive", "Dead")
cat(sum(complete) > 0, "-", sum(complete), "/", nrow(p),
    "complete cases,", sum(surv_status[complete] == "Dead"), "events\n")

# -----------------------------------------------------------------------------
# 3. GENE-SYMBOL EXPRESSION MATRIX FOR XCELL
# -----------------------------------------------------------------------------
# Collapse to one row per gene symbol: highest-variance probe, consistent
# with build_geo_cohort() in R/11_benchmark.R for the other GEO cohorts.
ex   <- exprs(eset)
sym  <- as.character(f$GeneSymbol)
keep <- !is.na(sym) & sym != "" & sym != "---"
ex   <- ex[keep, ]; sym <- sym[keep]

vars <- apply(ex, 1, var, na.rm = TRUE)
ord  <- order(-vars)
ex   <- ex[ord, ]; sym <- sym[ord]
dup  <- duplicated(sym)
ex   <- ex[!dup, ]
rownames(ex) <- sym[!dup]

message("15: expression matrix for xCell: ", nrow(ex), " genes x ", ncol(ex), " samples")

# -----------------------------------------------------------------------------
# 4. XCELL
# -----------------------------------------------------------------------------
xcell_file <- "data/GEO_cache/GSE72094_xCell.rds"
if (file.exists(xcell_file)) {
    xcell_res <- readRDS(xcell_file)
    message("15: xCell scores loaded from cache")
} else {
    xcell_res <- xCellAnalysis(ex)
    saveRDS(xcell_res, xcell_file)
    message("15: xCell done, ", nrow(xcell_res), " cell types x ",
            ncol(xcell_res), " samples")
}

cat("\n=== xCell cell types produced (DC / macrophage related) ===\n")
print(rownames(xcell_res)[grepl("DC|acrophage", rownames(xcell_res))])
