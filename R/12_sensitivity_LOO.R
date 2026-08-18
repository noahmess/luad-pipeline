# =============================================================================
# 12_sensitivity_LOO.R
# Leave-one-gene-out sensitivity analysis of Score A and nested likelihood-
# ratio test for the marginal contribution of NDC1.
#
# Output:
#   outputs/tables/Supplementary_Table_S3_LOO.csv
#     Five Score A variants (full + four leave-one-out) on four cohorts.
#   Console summary of the nested LRT on TCGA full.
#
# Prerequisite:
#   BENCHMARK_ENV.RData (produced by 11_benchmark.R) provides cohorts and
#   helper functions. If missing, this script aborts with a clear message.
#
# Cohorts tested:
#   TCGA full     : n = 394, harmonised survival cohort
#   TCGA held-out : n = 118, 30% partition
#   GSE31210      : n = 204
#   GSE50081      : n = 181
#
# Metrics per variant per cohort:
#   HR per 1 SD, Cox p-value, C-index with 95% CI,
#   time-dependent AUC at 12, 36, 60 months, delta_C vs full.
#
# Nested LRT on TCGA full:
#   M0 : Surv ~ LMNB2 + TMPO + LBR
#   M1 : Surv ~ LMNB2 + TMPO + LBR + NDC1
#   Reports chi-square, p, AIC, C-index, NDC1 coefficient with 95% CI.
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
    library(survival)
    library(timeROC)
})

# -----------------------------------------------------------------------------
# 1. LOAD BENCHMARK ENVIRONMENT
# -----------------------------------------------------------------------------
if (!file.exists("BENCHMARK_ENV.RData")) {
    stop("BENCHMARK_ENV.RData not found. Run R/11_benchmark.R first.")
}
load("BENCHMARK_ENV.RData")
message("12: loaded BENCHMARK_ENV.RData")

# -----------------------------------------------------------------------------
# 2. DEFINE LEAVE-ONE-GENE-OUT SIGNATURES
# -----------------------------------------------------------------------------
GENES_FULL <- c("LMNB2", "TMPO", "NDC1", "LBR")

make_sig <- function(genes) {
    list(genes  = genes,
         coefs  = setNames(rep(1, length(genes)), genes),
         zscore = TRUE,
         refit  = FALSE)
}

signatures_loo <- list(
    Score_A_full     = make_sig(GENES_FULL),
    Score_A_no_LMNB2 = make_sig(setdiff(GENES_FULL, "LMNB2")),
    Score_A_no_TMPO  = make_sig(setdiff(GENES_FULL, "TMPO")),
    Score_A_no_NDC1  = make_sig(setdiff(GENES_FULL, "NDC1")),
    Score_A_no_LBR   = make_sig(setdiff(GENES_FULL, "LBR"))
)

# -----------------------------------------------------------------------------
# 3. BUILD THE SENSITIVITY TABLE
# -----------------------------------------------------------------------------
sens_row <- function(sig, sig_name, cohort_name) {
    ch    <- cohorts[[cohort_name]]
    score <- compute_score(sig, ch$expr)
    hr    <- hr_pvalue(score, ch$surv)
    ci    <- cindex_ci(score, ch$surv)
    auc   <- auc_times(score, ch$surv)
    data.frame(
        cohort       = cohort_name,
        signature    = sig_name,
        gene_removed = if (sig_name == "Score_A_full") "none"
                       else sub("Score_A_no_", "", sig_name),
        n_genes      = length(sig$genes),
        HR           = round(unname(hr["HR"]), 3),
        HR_low       = round(unname(hr["HR_low"]), 3),
        HR_high      = round(unname(hr["HR_high"]), 3),
        p            = signif(unname(hr["p"]), 3),
        C            = round(unname(ci["C"]), 3),
        C_low        = round(unname(ci["lower"]), 3),
        C_high       = round(unname(ci["upper"]), 3),
        AUC_12m      = round(unname(auc["AUC_12m"]), 3),
        AUC_36m      = round(unname(auc["AUC_36m"]), 3),
        AUC_60m      = round(unname(auc["AUC_60m"]), 3)
    )
}

cohorts_to_test <- c("TCGA", "TCGA_heldout", "GSE31210", "GSE50081")

sens_table <- do.call(rbind, lapply(cohorts_to_test, function(cn) {
    do.call(rbind, lapply(names(signatures_loo), function(sn) {
        message("12: ", sn, " on ", cn)
        sens_row(signatures_loo[[sn]], sn, cn)
    }))
}))

# Add delta_C versus the full Score A in the same cohort
sens_table$delta_C <- NA_real_
for (cn in cohorts_to_test) {
    ref_C <- sens_table$C[sens_table$cohort == cn &
                          sens_table$signature == "Score_A_full"]
    rows  <- sens_table$cohort == cn
    sens_table$delta_C[rows] <- round(sens_table$C[rows] - ref_C, 3)
}

# -----------------------------------------------------------------------------
# 4. NESTED LRT ON TCGA FULL
# -----------------------------------------------------------------------------
tcga_expr <- cohorts$TCGA$expr
tcga_surv <- cohorts$TCGA$surv
zof <- function(g) scale(tcga_expr[, g])[, 1]

df_tcga <- data.frame(
    time   = tcga_surv[, "time"],
    status = tcga_surv[, "status"],
    LMNB2  = zof("LMNB2"),
    TMPO   = zof("TMPO"),
    LBR    = zof("LBR"),
    NDC1   = zof("NDC1")
)

M0 <- coxph(Surv(time, status) ~ LMNB2 + TMPO + LBR,        data = df_tcga)
M1 <- coxph(Surv(time, status) ~ LMNB2 + TMPO + LBR + NDC1, data = df_tcga)

lrt_anova <- anova(M0, M1)
lrt_chisq <- round(lrt_anova[["Chisq"]][2], 3)
lrt_p     <- round(lrt_anova[["Pr(>|Chi|)"]][2], 4)
aic_M0    <- round(AIC(M0), 2)
aic_M1    <- round(AIC(M1), 2)
c_M0      <- round(summary(M0)$concordance[1], 3)
c_M1      <- round(summary(M1)$concordance[1], 3)
ndc1_hr   <- round(summary(M1)$conf.int["NDC1", c(1, 3, 4)], 3)
ndc1_p    <- round(summary(M1)$coefficients["NDC1", 5], 3)

# -----------------------------------------------------------------------------
# 5. WRITE OUTPUT
# -----------------------------------------------------------------------------
dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)
write.csv(sens_table,
          "outputs/tables/Supplementary_Table_S3_LOO.csv",
          row.names = FALSE)

message("12: wrote outputs/tables/Supplementary_Table_S3_LOO.csv")

cat("\n=== Leave-one-gene-out sensitivity (Supplementary Table S3) ===\n")
print(sens_table[, c("cohort", "signature", "gene_removed", "n_genes",
                     "HR", "p", "C", "C_low", "C_high",
                     "AUC_12m", "AUC_36m", "AUC_60m", "delta_C")],
      row.names = FALSE)

cat("\n=== Nested LRT on TCGA full ===\n")
cat("M0 : LMNB2 + TMPO + LBR\n")
cat("M1 : LMNB2 + TMPO + LBR + NDC1\n\n")
cat("LRT chi-square :", lrt_chisq, "\n")
cat("LRT p-value    :", lrt_p, "\n")
cat("AIC M0         :", aic_M0, "\n")
cat("AIC M1         :", aic_M1, "\n")
cat("Delta AIC      :", round(aic_M0 - aic_M1, 2),
    "(positive favours M1)\n")
cat("C-index M0     :", c_M0, "\n")
cat("C-index M1     :", c_M1, "\n")
cat("Delta C        :", round(c_M1 - c_M0, 3), "\n")
cat("NDC1 in M1     : HR =", ndc1_hr[1], ", 95% CI",
    ndc1_hr[2], "-", ndc1_hr[3], ", p =", ndc1_p, "\n")
