# =============================================================================
# 13_random_signature_benchmark.R
# Permutation test: is Score A better than chance, and better than a
# proliferation-matched null?
#
# Score A (LMNB2 + TMPO + NDC1 + LBR, equal-weight Z-score, see
# R/10_signatures.R) is compared against N_PERM random 4-gene panels, drawn
# from two different null pools, on four cohorts:
#
#   TCGA          n = 394  DISCOVERY cohort - Score A's genes (incl. the
#                          exclusion of LMNB1) were selected on this exact
#                          cohort in R/05_survival_cox.R. Estimates here are
#                          optimistic / circular, not a validation.
#   TCGA_heldout  n = 118  30% partition of TCGA, internal validation.
#   GSE31210      n = 204  external validation, Japanese cohort.
#   GSE50081      n = 181  external validation, Canadian cohort.
#
# Two null pools per cohort:
#   protein_coding          any protein-coding gene expressed in that cohort
#                            (excluding Score A's own genes).
#   proliferation_matched   protein-coding genes whose |Spearman rho| with
#                            MKI67 exceeds 0.40 in that cohort (excluding
#                            MKI67 itself and Score A's own genes). This is
#                            the harder null: it asks whether Score A beats
#                            a generic proliferation signal, not just any
#                            random gene set.
#
# Each random panel is scored with the exact same method as Score A (mean
# Z-score, no refitting), so the only thing that differs is gene
# composition. Empirical p-value = proportion of random panels that match
# or beat Score A, out of N_PERM + 1 (observed counted as one draw).
#
# Output:
#   outputs/tables/Random_signature_benchmark_detail.csv
#     One row per random panel (all cohorts x both null pools) + one row
#     per cohort for the observed Score A.
#   outputs/tables/Random_signature_benchmark_summary.csv
#     Observed vs random mean/SD, empirical p (C-index and Cox p), per
#     cohort x null pool x metric. gene_selection_cohort flags TCGA (the
#     discovery cohort) so optimistic and independent estimates aren't
#     read off the same line.
#
# Prerequisite:
#   BENCHMARK_ENV.RData (produced by R/11_benchmark.R) for the four
#   cohorts (expr, surv) and the cindex_ci / hr_pvalue / compute_score
#   helpers. data/TCGA_LUAD_data.RData for GENCODE gene_type annotation.
#   If missing, this script aborts with a clear message.
# =============================================================================

setwd("~/Desktop/luad-pipeline")

suppressPackageStartupMessages({
    library(survival)
})

source("R/00_setup.R")
source("R/10_signatures.R")

N_PERM       <- 1000
MKI67_RHO_CUTOFF <- 0.40
set.seed(SEED)

# -----------------------------------------------------------------------------
# 1. LOAD BENCHMARK ENVIRONMENT
# -----------------------------------------------------------------------------
if (!file.exists("BENCHMARK_ENV.RData")) {
    stop("BENCHMARK_ENV.RData not found. Run R/11_benchmark.R first.")
}
load("BENCHMARK_ENV.RData")
message("13: loaded BENCHMARK_ENV.RData")

cohort_names <- c("TCGA", "TCGA_heldout", "GSE31210", "GSE50081")
# The cohort on which Score A's genes were selected (R/05_survival_cox.R,
# fixed n = 394 cohort). Everything else is independent validation.
discovery_cohort <- "TCGA"

n_sig_genes <- length(GENES_SCORE_A)

# -----------------------------------------------------------------------------
# 2. PROTEIN-CODING ANNOTATION (universal HGNC symbol list, GENCODE/GDC)
# -----------------------------------------------------------------------------
# Restrict every random pool to protein-coding genes only. Score A is made
# of four protein-coding genes; comparing it against a pool that also
# contains pseudogenes, lncRNAs and other non-coding biotypes (noisier,
# often near-zero counts) would bias random panels downward and make the
# permutation test look more favourable than a fair like-for-like
# comparison.
message("13: loading gene biotype annotation")
.ge <- new.env()
load("data/TCGA_LUAD_data.RData", envir = .ge)
row_info <- SummarizedExperiment::rowData(get("data", envir = .ge))
protein_coding_symbols <- unique(
    row_info$gene_name[row_info$gene_type == "protein_coding"]
)
rm(.ge)

# -----------------------------------------------------------------------------
# 3. HELPERS
# -----------------------------------------------------------------------------
make_sig <- function(genes) {
    list(genes  = genes,
         coefs  = setNames(rep(1, length(genes)), genes),
         zscore = TRUE,
         refit  = FALSE)
}

score_panel <- function(genes, expr, surv) {
    sig   <- make_sig(genes)
    score <- compute_score(sig, expr)
    hr    <- hr_pvalue(score, surv)
    ci    <- cindex_ci(score, surv)
    data.frame(
        genes   = paste(genes, collapse = "|"),
        n_genes = length(genes),
        HR      = round(unname(hr["HR"]), 3),
        HR_low  = round(unname(hr["HR_low"]), 3),
        HR_high = round(unname(hr["HR_high"]), 3),
        p       = signif(unname(hr["p"]), 3),
        C       = round(unname(ci["C"]), 3),
        C_low   = round(unname(ci["lower"]), 3),
        C_high  = round(unname(ci["upper"]), 3)
    )
}

# Build both null pools for one cohort's expression matrix.
build_pools <- function(expr) {
    gene_var <- apply(expr, 2, var, na.rm = TRUE)
    expressed <- colnames(expr)[!is.na(gene_var) & gene_var > 0]
    coding    <- intersect(expressed, protein_coding_symbols)

    pool_protein_coding <- setdiff(coding, GENES_SCORE_A)

    pool_prolif <- character(0)
    if ("MKI67" %in% colnames(expr)) {
        rho <- suppressWarnings(
            cor(expr[, coding, drop = FALSE], expr[, "MKI67"],
                method = "spearman", use = "pairwise.complete.obs")[, 1]
        )
        prolif_genes <- names(rho)[!is.na(rho) & abs(rho) > MKI67_RHO_CUTOFF]
        pool_prolif  <- setdiff(prolif_genes, c(GENES_SCORE_A, "MKI67"))
    }

    list(protein_coding = pool_protein_coding,
        proliferation_matched = pool_prolif)
}

# Run N_PERM random panels from one pool and return detail + summary rows.
run_null <- function(pool, expr, surv, cohort, null_type, obs_C, obs_p) {
    if (length(pool) < n_sig_genes) {
        message("13:   ", cohort, " / ", null_type,
                ": pool too small (", length(pool), " genes), skipped")
        return(list(detail = NULL, summary = NULL))
    }
    panels <- replicate(N_PERM, sample(pool, n_sig_genes), simplify = FALSE)
    res <- do.call(rbind, lapply(seq_along(panels), function(i) {
        if (i %% 250 == 0)
            message("13:   ", cohort, " / ", null_type,
                    ": permutation ", i, "/", N_PERM)
        row <- score_panel(panels[[i]], expr, surv)
        row$cohort    <- cohort
        row$null_type <- null_type
        row$type      <- "random"
        row$perm      <- i
        row
    }))

    perm_C <- res$C
    perm_p <- res$p
    p_emp_C <- (1 + sum(perm_C >= obs_C, na.rm = TRUE)) /
               (1 + sum(!is.na(perm_C)))
    p_emp_p <- (1 + sum(perm_p <= obs_p, na.rm = TRUE)) /
               (1 + sum(!is.na(perm_p)))

    summary_rows <- data.frame(
        cohort      = cohort,
        null_type   = null_type,
        metric      = c("C_index", "Cox_p"),
        observed    = c(obs_C, obs_p),
        random_mean = c(mean(perm_C, na.rm = TRUE), mean(perm_p, na.rm = TRUE)),
        random_sd   = c(sd(perm_C, na.rm = TRUE),   sd(perm_p, na.rm = TRUE)),
        random_min  = c(min(perm_C, na.rm = TRUE),  min(perm_p, na.rm = TRUE)),
        random_max  = c(max(perm_C, na.rm = TRUE),  max(perm_p, na.rm = TRUE)),
        n_perm      = c(sum(!is.na(perm_C)), sum(!is.na(perm_p))),
        empirical_p = c(p_emp_C, p_emp_p)
    )
    list(detail = res, summary = summary_rows)
}

# -----------------------------------------------------------------------------
# 4. RUN PER COHORT
# -----------------------------------------------------------------------------
detail_rows  <- list()
summary_rows <- list()

for (cohort in cohort_names) {
    ch   <- cohorts[[cohort]]
    expr <- ch$expr
    surv <- ch$surv
    message("13: === ", cohort, " (n = ", nrow(expr), ") ===")

    observed <- score_panel(GENES_SCORE_A, expr, surv)
    observed$cohort    <- cohort
    observed$null_type <- NA_character_
    observed$type      <- "observed"
    observed$perm      <- 0L
    detail_rows[[paste0(cohort, "_observed")]] <- observed

    pools <- build_pools(expr)

    for (null_type in names(pools)) {
        out <- run_null(pools[[null_type]], expr, surv, cohort, null_type,
                        obs_C = observed$C, obs_p = observed$p)
        if (!is.null(out$detail)) {
            detail_rows[[paste0(cohort, "_", null_type)]] <- out$detail
            summary_rows[[paste0(cohort, "_", null_type)]] <- out$summary
        }
    }
}

benchmark_detail  <- do.call(rbind, detail_rows)
benchmark_summary <- do.call(rbind, summary_rows)

benchmark_detail$gene_selection_cohort  <- benchmark_detail$cohort == discovery_cohort
benchmark_summary$gene_selection_cohort <- benchmark_summary$cohort == discovery_cohort

benchmark_detail <- benchmark_detail[, c("cohort", "gene_selection_cohort",
                                         "null_type", "type", "perm",
                                         "genes", "n_genes",
                                         "HR", "HR_low", "HR_high", "p",
                                         "C", "C_low", "C_high")]
benchmark_summary <- benchmark_summary[, c("cohort", "gene_selection_cohort",
                                           "null_type", "metric", "observed",
                                           "random_mean", "random_sd",
                                           "random_min", "random_max",
                                           "n_perm", "empirical_p")]

# -----------------------------------------------------------------------------
# 5. WRITE OUTPUTS
# -----------------------------------------------------------------------------
dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)

old_file <- "outputs/tables/Random_signature_benchmark_TCGA.csv"
if (file.exists(old_file)) {
    file.remove(old_file)
    message("13: removed superseded ", old_file)
}

write.csv(benchmark_detail,
          "outputs/tables/Random_signature_benchmark_detail.csv",
          row.names = FALSE)
write.csv(benchmark_summary,
          "outputs/tables/Random_signature_benchmark_summary.csv",
          row.names = FALSE)

message("13: wrote outputs/tables/Random_signature_benchmark_detail.csv")
message("13: wrote outputs/tables/Random_signature_benchmark_summary.csv")

cat("\n=== Empirical significance, all cohorts x null pools ===\n")
print(benchmark_summary, row.names = FALSE)

cat("\nNote: gene_selection_cohort = TRUE (", discovery_cohort, ") is the ",
    "discovery cohort where Score A's genes were chosen - its p-values are ",
    "optimistic, not an independent validation.\n", sep = "")
