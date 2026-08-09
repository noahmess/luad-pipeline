# =============================================================================
# 11_benchmark.R
# Comparative benchmark of Score A against four published signatures across
# four cohorts.
#
# Output:
#   outputs/tables/Benchmark_main_3cohorts.csv
#     Score A, LASSO-Cox, Li2023_TME6, Zhou2022, CCP_like
#     on TCGA held-out (n = 118), GSE31210 (n = 204), GSE50081 (n = 181).
#   outputs/tables/Benchmark_suppl_TCGAfull.csv
#     Same signatures on TCGA full (n = 394).
#   BENCHMARK_ENV.RData
#     Saved environment with cohorts and helper functions for downstream use
#     (12_sensitivity_LOO.R).
#
# Metrics per cohort:
#   HR per 1 SD (continuous Cox), Cox p-value
#   Harrell C-index with 95% CI
#   time-dependent AUC at 12, 36, 60 months (timeROC, Blanche 2013)
#
# Cohorts:
#   TCGA full     : n = 394, harmonised survival cohort
#   TCGA train    : 70% partition stratified on events, seed = 42
#   TCGA held-out : 30% partition, used as primary internal validation
#   GSE31210      : n = 204, Japanese cohort, 30 events
#   GSE50081      : n = 181, Canadian cohort, 75 events
#
# Note on TCGA expression. The TCGA expression matrix is built from the
# DESeq2-normalised counts (norm_tumor) produced by 04_cohort_harmonization,
# log2-transformed. Using raw counts here would change the variance structure
# and shift Score A statistics away from the manuscript values.
# =============================================================================

setwd("~/Desktop/luad-pipeline")

suppressPackageStartupMessages({
    library(survival)
    library(glmnet)
    library(timeROC)
    library(GEOquery)
    library(Biobase)
    library(SummarizedExperiment)
})

source("R/00_setup.R")
source("R/10_signatures.R")

set.seed(42)
options(timeout = 1000)

# -----------------------------------------------------------------------------
# 1. HELPER FUNCTIONS
# -----------------------------------------------------------------------------

.standardise <- function(x) {
    s <- apply(x, 2, function(col) {
        m <- mean(col, na.rm = TRUE); sdv <- sd(col, na.rm = TRUE)
        if (is.na(sdv) || sdv == 0) rep(0, length(col)) else (col - m) / sdv
    })
    matrix(s, nrow = nrow(x), dimnames = dimnames(x))
}

compute_score <- function(sig, expr, train_expr = NULL, train_surv = NULL) {
    if (isTRUE(sig$refit)) {
        return(.fit_lasso_score(sig, expr, train_expr, train_surv))
    }
    # Apply alias map: rename platform columns from alias to canonical symbol
    expr <- .apply_alias_map(expr, sig$genes)
    present <- intersect(sig$genes, colnames(expr))
    if (length(present) == 0) return(rep(NA_real_, nrow(expr)))
    X <- as.matrix(expr[, present, drop = FALSE])
    w <- sig$coefs[present]
    if (isTRUE(sig$zscore)) X <- .standardise(X)
    as.numeric(X %*% w) / ifelse(isTRUE(sig$zscore), length(present), 1)
}

.apply_alias_map <- function(expr, target_genes) {
    expr_cols <- colnames(expr)
    for (canonical in names(SIGNATURE_ALIAS_MAP)) {
        if (canonical %in% target_genes && !(canonical %in% expr_cols)) {
            aliases <- SIGNATURE_ALIAS_MAP[[canonical]]
            found_alias <- intersect(aliases, expr_cols)
            if (length(found_alias) > 0) {
                idx <- match(found_alias[1], colnames(expr))
                colnames(expr)[idx] <- canonical
            }
        }
    }
    expr
}

.fit_lasso_score <- function(sig, expr, train_expr, train_surv) {
    if (is.null(train_expr) || is.null(train_surv))
        stop("LASSO refit requires train_expr and train_surv")
    candidate_genes <- intersect(PANEL_9, colnames(train_expr))
    Xtr <- as.matrix(train_expr[, candidate_genes, drop = FALSE])
    Xtr <- .standardise(Xtr)
    ytr <- Surv(train_surv[, "time"], train_surv[, "status"])
    cv  <- cv.glmnet(Xtr, ytr, family = "cox", alpha = 1, nfolds = 10)
    coefs <- as.numeric(coef(cv, s = "lambda.min"))
    names(coefs) <- candidate_genes
    selected <- coefs[coefs != 0]
    if (length(selected) == 0) return(rep(NA_real_, nrow(expr)))
    # Apply to test expression
    test_genes <- intersect(names(selected), colnames(expr))
    if (length(test_genes) == 0) return(rep(NA_real_, nrow(expr)))
    Xte <- as.matrix(expr[, test_genes, drop = FALSE])
    Xte <- .standardise(Xte)
    as.numeric(Xte %*% selected[test_genes])
}

cindex_ci <- function(score, surv) {
    ok <- !is.na(score)
    if (sum(ok) < 5 || length(unique(score[ok])) < 2)
        return(c(C = NA, lower = NA, upper = NA, se = NA))
    df <- data.frame(time = surv[ok, "time"], status = surv[ok, "status"],
                     s = score[ok])
    fit <- tryCatch(coxph(Surv(time, status) ~ s, data = df),
                    error = function(e) NULL)
    if (is.null(fit)) return(c(C = NA, lower = NA, upper = NA, se = NA))
    cc <- summary(fit)$concordance
    C  <- unname(cc[1]); se <- unname(cc[2])
    c(C = C, lower = C - 1.96 * se, upper = C + 1.96 * se, se = se)
}

hr_pvalue <- function(score, surv) {
    ok <- !is.na(score)
    if (sum(ok) < 5 || length(unique(score[ok])) < 2)
        return(c(HR = NA, HR_low = NA, HR_high = NA, p = NA))
    df <- data.frame(time = surv[ok, "time"], status = surv[ok, "status"],
                     s = scale(score[ok])[, 1])
    fit <- tryCatch(coxph(Surv(time, status) ~ s, data = df),
                    error = function(e) NULL)
    if (is.null(fit)) return(c(HR = NA, HR_low = NA, HR_high = NA, p = NA))
    s <- summary(fit)
    c(HR = s$conf.int[1], HR_low = s$conf.int[3],
      HR_high = s$conf.int[4], p = s$coefficients[5])
}

auc_times <- function(score, surv, times = c(12, 36, 60)) {
    ok <- !is.na(score)
    if (sum(ok) < 10)
        return(setNames(rep(NA, length(times)), paste0("AUC_", times, "m")))
    res <- tryCatch(
        timeROC(T = surv[ok, "time"], delta = surv[ok, "status"],
                marker = score[ok], cause = 1, times = times, iid = FALSE),
        error = function(e) NULL)
    if (is.null(res))
        return(setNames(rep(NA, length(times)), paste0("AUC_", times, "m")))
    setNames(as.numeric(res$AUC), paste0("AUC_", times, "m"))
}

# -----------------------------------------------------------------------------
# 2. BUILD TCGA COHORT FROM THE HARMONISED RData
# -----------------------------------------------------------------------------

message("11: building TCGA cohort")
load("data/TCGA_LUAD_cohort.RData")
# Provides norm_tumor (genes x 540, DESeq2-normalised counts), expr9, expr4_v,
# time_v, status_v, surv_obj, etc.

load("data/TCGA_LUAD_data.RData")
# Provides data and gene_symbols (used to name the rows of norm_tumor).

# Build a (samples x genes) expression matrix on the 394 survival patients,
# using the DESeq2-normalised counts already produced by 04_cohort_harmonization.
# This keeps Score A and all benchmark statistics aligned with the manuscript.
ids394   <- rownames(expr4_v)
gene_symbols <- as.character(rowData(data)$gene_name)
expr_sym <- norm_tumor
rownames(expr_sym) <- gene_symbols
# Keep one row per symbol, the one with the highest mean expression
expr_sym <- expr_sym[order(-rowMeans(expr_sym)), ]
expr_sym <- expr_sym[!duplicated(rownames(expr_sym)), ]
expr_tcga_log <- t(log2(expr_sym[, ids394] + 1))

cohorts <- list()
cohorts$TCGA <- list(expr = expr_tcga_log, surv = surv_obj)

# Partition 70/30 stratified on events
events     <- surv_obj[, "status"]
idx_events <- which(events == 1)
idx_censor <- which(events == 0)
train_idx  <- sort(c(sample(idx_events, round(0.7 * length(idx_events))),
                     sample(idx_censor, round(0.7 * length(idx_censor)))))
test_idx   <- setdiff(seq_along(events), train_idx)

cohorts$TCGA_train   <- list(expr = expr_tcga_log[train_idx, ],
                             surv = surv_obj[train_idx, ])
cohorts$TCGA_heldout <- list(expr = expr_tcga_log[test_idx, ],
                             surv = surv_obj[test_idx, ])

message("11: TCGA train ", length(train_idx), " patients, ",
        sum(events[train_idx]), " events")
message("11: TCGA held-out ", length(test_idx), " patients, ",
        sum(events[test_idx]), " events")

# -----------------------------------------------------------------------------
# 3. BUILD GEO COHORTS WITH CACHING
# -----------------------------------------------------------------------------

dir.create("data/GEO_cache", recursive = TRUE, showWarnings = FALSE)

safe_get_geo <- function(gse_id) {
    rds_file <- file.path("data/GEO_cache", paste0(gse_id, "_eset.rds"))
    if (file.exists(rds_file)) return(readRDS(rds_file))
    g    <- getGEO(gse_id, GSEMatrix = TRUE, getGPL = TRUE,
                   destdir = "data/GEO_cache")
    eset <- g[[1]]
    saveRDS(eset, rds_file)
    eset
}

build_geo_cohort <- function(gse_id, time_col, status_col, time_unit,
                             exclusion_col = NULL) {
    eset  <- safe_get_geo(gse_id)
    expr  <- exprs(eset); pheno <- pData(eset); feat <- fData(eset)
    sym   <- sapply(strsplit(as.character(feat[["Gene Symbol"]]),
                             " /// ", fixed = TRUE), `[`, 1)
    keep  <- !is.na(sym) & sym != "" & sym != "---"
    expr  <- expr[keep, ]; sym <- sym[keep]
    # Keep one probe per gene: the one with the highest variance
    vars  <- apply(expr, 1, var, na.rm = TRUE)
    ord   <- order(-vars); expr <- expr[ord, ]; sym <- sym[ord]
    dup   <- duplicated(sym); expr <- expr[!dup, ]
    rownames(expr) <- sym[!dup]
    # Survival
    time_v <- as.numeric(as.character(pheno[[time_col]]))
    if (time_unit == "days")  time_v <- time_v / 30.44
    if (time_unit == "years") time_v <- time_v * 12
    status_v <- as.integer(pheno[[status_col]] == "dead")
    keep_pat <- !is.na(time_v) & !is.na(status_v) & time_v > 0
    if (!is.null(exclusion_col)) {
        excl <- pheno[[exclusion_col]]
        keep_pat <- keep_pat & (is.na(excl) | excl == "none")
    }
    list(expr = t(expr[, keep_pat, drop = FALSE]),
         surv = Surv(time_v[keep_pat], status_v[keep_pat]))
}

message("11: building GSE31210")
cohorts$GSE31210 <- build_geo_cohort(
    "GSE31210",
    time_col      = "days before death/censor:ch1",
    status_col    = "death:ch1",
    time_unit     = "days",
    exclusion_col = "exclude for prognosis analysis due to incomplete resection or adjuvant therapy:ch1"
)
message("11: GSE31210 ", nrow(cohorts$GSE31210$expr), " patients, ",
        sum(cohorts$GSE31210$surv[, 2]), " events")

message("11: building GSE50081")
cohorts$GSE50081 <- build_geo_cohort(
    "GSE50081",
    time_col   = "survival time:ch1",
    status_col = "status:ch1",
    time_unit  = "years"
)
message("11: GSE50081 ", nrow(cohorts$GSE50081$expr), " patients, ",
        sum(cohorts$GSE50081$surv[, 2]), " events")

# -----------------------------------------------------------------------------
# 4. RUN THE BENCHMARK
# -----------------------------------------------------------------------------

benchmark_row <- function(sig_name, sig, cohort_name) {
    ch    <- cohorts[[cohort_name]]
    train <- if (isTRUE(sig$refit)) cohorts$TCGA_train else NULL
    score <- compute_score(sig, ch$expr,
                           train_expr = if (!is.null(train)) train$expr else NULL,
                           train_surv = if (!is.null(train)) train$surv else NULL)
    hr  <- hr_pvalue(score, ch$surv)
    ci  <- cindex_ci(score, ch$surv)
    auc <- auc_times(score, ch$surv)
    # Gene coverage
    expr_a <- .apply_alias_map(ch$expr, sig$genes)
    n_found <- length(intersect(sig$genes, colnames(expr_a)))
    n_total <- length(sig$genes)
    data.frame(
        cohort    = cohort_name,
        signature = sig_name,
        n_genes   = paste0(n_found, "/", n_total),
        HR        = round(unname(hr["HR"]), 3),
        HR_low    = round(unname(hr["HR_low"]), 3),
        HR_high   = round(unname(hr["HR_high"]), 3),
        p         = signif(unname(hr["p"]), 3),
        C         = round(unname(ci["C"]), 3),
        C_low     = round(unname(ci["lower"]), 3),
        C_high    = round(unname(ci["upper"]), 3),
        AUC_12m   = round(unname(auc["AUC_12m"]), 3),
        AUC_36m   = round(unname(auc["AUC_36m"]), 3),
        AUC_60m   = round(unname(auc["AUC_60m"]), 3),
        notes     = NA_character_
    )
}

sig_names <- names(signatures)

# Main benchmark on three cohorts
main_cohorts <- c("TCGA_heldout", "GSE31210", "GSE50081")
benchmark_main <- do.call(rbind, lapply(main_cohorts, function(cn) {
    do.call(rbind, lapply(sig_names, function(sn) {
        message("11: ", sn, " on ", cn)
        benchmark_row(sn, signatures[[sn]], cn)
    }))
}))

# Supplementary benchmark on TCGA full
benchmark_suppl <- do.call(rbind, lapply(sig_names, function(sn) {
    message("11: ", sn, " on TCGA full")
    benchmark_row(sn, signatures[[sn]], "TCGA")
}))

# -----------------------------------------------------------------------------
# 5. WRITE OUTPUTS
# -----------------------------------------------------------------------------

dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)

write.csv(benchmark_main,
          "outputs/tables/Benchmark_main_3cohorts.csv",
          row.names = FALSE)
write.csv(benchmark_suppl,
          "outputs/tables/Benchmark_suppl_TCGAfull.csv",
          row.names = FALSE)

# Save the environment for 12_sensitivity_LOO.R to reuse without recomputing
save(cohorts, compute_score, cindex_ci, hr_pvalue, auc_times,
     .standardise, .apply_alias_map, .fit_lasso_score,
     SIGNATURE_ALIAS_MAP,
     file = "BENCHMARK_ENV.RData")

message("11: wrote outputs/tables/Benchmark_main_3cohorts.csv")
message("11: wrote outputs/tables/Benchmark_suppl_TCGAfull.csv")
message("11: saved BENCHMARK_ENV.RData")

cat("\n=== Main benchmark, three cohorts ===\n")
print(benchmark_main)

cat("\n=== Supplementary benchmark, TCGA full ===\n")
print(benchmark_suppl)
