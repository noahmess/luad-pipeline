# =============================================================================
# incremental_analyses.R
# Three analyses to strengthen the Score A manuscript:
#   1. Incremental value: clinical model vs clinical + Score A (LR test + C-index)
#   2. Stage I stratification: does Score A work in stage I alone?
#   3. Harrell C-index: Score A vs clinical vs combined
#
# Requires pipeline objects in memory:
#   age_v, sex_v, stage_v, time_v, status_v, expr4_v (all length 394)
#
# Run:  source("/Users/nouha/Desktop/luad-pipeline/incremental_analyses.R")
# =============================================================================

suppressPackageStartupMessages({
  library(survival)
  library(survminer)
})

FIG_DIR <- file.path(PROJECT_ROOT, "outputs", "figures")
TAB_DIR <- file.path(PROJECT_ROOT, "outputs", "tables")

n <- length(time_v)
stopifnot(n == 394)

# ---- Compute Score A ----
scoreA <- rowMeans(scale(expr4_v))

# ---- Build clinical variables ----
# Stage as numeric ordinal (I=1, II=2, III=3, IV=4)
stage_num <- as.integer(factor(stage_v, levels = c("I","II","III","IV")))
sex_bin   <- as.integer(sex_v == "male")  # 1 = male

df <- data.frame(
  time    = time_v,
  status  = status_v,
  age     = age_v,
  sex     = sex_bin,
  stage   = stage_num,
  scoreA  = scoreA,
  stage_f = factor(stage_v, levels = c("I","II","III","IV"))
)

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════╗\n")
cat("║            INCREMENTAL VALUE & C-INDEX ANALYSES                ║\n")
cat("║            Harmonized cohort n = 394 (129 events)              ║\n")
cat("╚══════════════════════════════════════════════════════════════════╝\n")

# =====================================================================
# 1. INCREMENTAL VALUE: Clinical vs Clinical + Score A
# =====================================================================
cat("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("  1. INCREMENTAL VALUE OF SCORE A OVER CLINICAL VARIABLES\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

# Model 1: Clinical only (age + sex + stage)
cox_clin <- coxph(Surv(time, status) ~ age + sex + stage, data = df)
# Model 2: Clinical + Score A
cox_full <- coxph(Surv(time, status) ~ age + sex + stage + scoreA, data = df)

cat("Model 1 (Clinical only: age + sex + stage):\n")
print(summary(cox_clin)$coefficients)
cat("\nModel 2 (Clinical + Score A):\n")
print(summary(cox_full)$coefficients)

# Likelihood Ratio Test
lr_test <- anova(cox_clin, cox_full)
lr_p    <- lr_test[["Pr(>|Chi|)"]][2]
lr_chi  <- lr_test[["Chisq"]][2]
cat(sprintf("\nLikelihood Ratio Test (adding Score A): chi2 = %.2f, df = 1, p = %.4g\n",
            lr_chi, lr_p))

# C-index comparison
c_clin <- summary(cox_clin)$concordance
c_full <- summary(cox_full)$concordance
cat(sprintf("\nC-index (Clinical only):      %.3f (SE %.3f)\n", c_clin[1], c_clin[2]))
cat(sprintf("C-index (Clinical + ScoreA):  %.3f (SE %.3f)\n", c_full[1], c_full[2]))
cat(sprintf("Delta C-index:                +%.3f\n", c_full[1] - c_clin[1]))

# Score A alone
cox_scoreA <- coxph(Surv(time, status) ~ scoreA, data = df)
c_scoreA   <- summary(cox_scoreA)$concordance
cat(sprintf("C-index (Score A alone):      %.3f (SE %.3f)\n", c_scoreA[1], c_scoreA[2]))

# AIC comparison
cat(sprintf("\nAIC  Clinical only:    %.1f\n", AIC(cox_clin)))
cat(sprintf("AIC  Clinical + ScoreA: %.1f  (lower = better)\n", AIC(cox_full)))
cat(sprintf("Delta AIC:              %.1f\n", AIC(cox_clin) - AIC(cox_full)))

# Score A HR in the full model
hr_sa <- exp(coef(cox_full)["scoreA"])
ci_sa <- exp(confint(cox_full)["scoreA",])
p_sa  <- summary(cox_full)$coefficients["scoreA", "Pr(>|z|)"]
cat(sprintf("\nScore A in full model: HR = %.3f (%.3f - %.3f), p = %.4g\n",
            hr_sa, ci_sa[1], ci_sa[2], p_sa))

# =====================================================================
# 2. STAGE I STRATIFICATION
# =====================================================================
cat("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("  2. SCORE A STRATIFICATION IN STAGE I PATIENTS\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

stI <- df$stage_f == "I"
df_stI <- df[stI, ]
n_stI <- nrow(df_stI)
ev_stI <- sum(df_stI$status)
cat(sprintf("Stage I patients: n = %d, events = %d\n", n_stI, ev_stI))

# Cox continuous in stage I
cox_stI <- coxph(Surv(time, status) ~ scoreA, data = df_stI)
hr_stI  <- exp(coef(cox_stI))
ci_stI  <- exp(confint(cox_stI))
p_stI   <- summary(cox_stI)$coefficients[, "Pr(>|z|)"]
cat(sprintf("Cox continuous (Stage I): HR = %.3f (%.3f - %.3f), p = %.4g\n",
            hr_stI, ci_stI[1], ci_stI[2], p_stI))

# Cox adjusted for age + sex in stage I
cox_stI_adj <- coxph(Surv(time, status) ~ scoreA + age + sex, data = df_stI)
hr_stI_adj  <- exp(coef(cox_stI_adj)["scoreA"])
ci_stI_adj  <- exp(confint(cox_stI_adj)["scoreA",])
p_stI_adj   <- summary(cox_stI_adj)$coefficients["scoreA", "Pr(>|z|)"]
cat(sprintf("Cox adjusted age+sex (Stage I): HR = %.3f (%.3f - %.3f), p = %.4g\n",
            hr_stI_adj, ci_stI_adj[1], ci_stI_adj[2], p_stI_adj))

# KM median split in stage I
grp_stI <- factor(ifelse(df_stI$scoreA > median(df_stI$scoreA), "High", "Low"),
                  c("Low", "High"))
df_stI$grp <- grp_stI
sdiff_stI <- survdiff(Surv(time, status) ~ grp, data = df_stI)
lrp_stI   <- pchisq(sdiff_stI$chisq, 1, lower.tail = FALSE)
cat(sprintf("KM log-rank (Stage I): p = %.4g\n", lrp_stI))

# 5-year OS in stage I
km_stI <- survfit(Surv(time, status) ~ grp, data = df_stI)
s5_stI <- summary(km_stI, times = 60)
if (length(s5_stI$surv) == 2) {
  cat(sprintf("5-yr OS Stage I: Low = %.1f%%, High = %.1f%%\n",
              100*s5_stI$surv[1], 100*s5_stI$surv[2]))
}

# C-index in stage I
c_stI <- summary(cox_stI)$concordance
cat(sprintf("C-index (Score A, Stage I): %.3f (SE %.3f)\n", c_stI[1], c_stI[2]))

# KM plot for stage I
n_low <- sum(grp_stI == "Low"); n_high <- sum(grp_stI == "High")
pKM <- ggsurvplot(km_stI, data = df_stI, pval = TRUE, pval.size = 4.5,
  risk.table = TRUE, risk.table.height = 0.25,
  title = paste0("Score A in Stage I patients (n = ", n_stI, ")"),
  legend.labs = c(paste0("Low Risk (n=", n_low, ")"),
                  paste0("High Risk (n=", n_high, ")")),
  legend.title = "", xlab = "Time (months)",
  ylab = "Overall survival probability",
  palette = c("#2E86C1", "#E74C3C"),
  ggtheme = theme_classic(base_size = 12),
  surv.median.line = "hv", tables.theme = theme_cleantable())

png(file.path(FIG_DIR, "ScoreA_StageI_KM.png"), 7, 8, "in", res = 300)
print(pKM); dev.off()
tiff(file.path(FIG_DIR, "ScoreA_StageI_KM.tiff"), 7, 8, "in", res = 300)
print(pKM); dev.off()
cat("saved: ScoreA_StageI_KM.png + .tiff\n")

# =====================================================================
# 3. HARRELL C-INDEX SUMMARY
# =====================================================================
cat("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("  3. HARRELL C-INDEX COMPARISON\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

# Already computed above, consolidate
cox_stage  <- coxph(Surv(time, status) ~ stage, data = df)
c_stage    <- summary(cox_stage)$concordance

cat("Full cohort (n = 394):\n")
cat(sprintf("  Stage alone:            C = %.3f (SE %.3f)\n", c_stage[1], c_stage[2]))
cat(sprintf("  Score A alone:          C = %.3f (SE %.3f)\n", c_scoreA[1], c_scoreA[2]))
cat(sprintf("  Clinical (age+sex+st):  C = %.3f (SE %.3f)\n", c_clin[1], c_clin[2]))
cat(sprintf("  Clinical + Score A:     C = %.3f (SE %.3f)\n", c_full[1], c_full[2]))
cat(sprintf("\n  Gain from adding Score A to clinical model: +%.3f\n",
            c_full[1] - c_clin[1]))

cat(sprintf("\nStage I only (n = %d):\n", n_stI))
cat(sprintf("  Score A alone:          C = %.3f (SE %.3f)\n", c_stI[1], c_stI[2]))

# =====================================================================
# SUMMARY TABLE
# =====================================================================
cat("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("  PUBLICATION-READY SUMMARY\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

summary_df <- data.frame(
  Analysis = c(
    "Score A univariate (n=394)",
    "Score A adjusted age+sex+stage (n=394)",
    "Score A in Stage I univariate",
    "Score A in Stage I adjusted age+sex",
    "LR test: clinical vs clinical+ScoreA",
    "C-index clinical only",
    "C-index clinical + Score A",
    "C-index gain"
  ),
  Result = c(
    sprintf("HR=%.3f (%.3f-%.3f), p=%.4g", exp(coef(cox_scoreA)),
            exp(confint(cox_scoreA))[1], exp(confint(cox_scoreA))[2],
            summary(cox_scoreA)$coefficients[,"Pr(>|z|)"]),
    sprintf("HR=%.3f (%.3f-%.3f), p=%.4g", hr_sa, ci_sa[1], ci_sa[2], p_sa),
    sprintf("HR=%.3f (%.3f-%.3f), p=%.4g", hr_stI, ci_stI[1], ci_stI[2], p_stI),
    sprintf("HR=%.3f (%.3f-%.3f), p=%.4g", hr_stI_adj, ci_stI_adj[1], ci_stI_adj[2], p_stI_adj),
    sprintf("chi2=%.2f, p=%.4g", lr_chi, lr_p),
    sprintf("%.3f", c_clin[1]),
    sprintf("%.3f", c_full[1]),
    sprintf("+%.3f", c_full[1] - c_clin[1])
  ),
  stringsAsFactors = FALSE
)
print(summary_df, row.names = FALSE)
outpath <- file.path(TAB_DIR, "Incremental_analyses_summary.csv")
write.csv(summary_df, outpath, row.names = FALSE)
cat("\nSaved:", outpath, "\n")

cat("\n══════════════════════════════════════════════════════════════════\n")
cat("  DONE. Key results for the manuscript:\n")
cat("  - LR test p-value tells if Score A adds to clinical model\n")
cat("  - C-index gain quantifies the improvement\n")
cat("  - Stage I analysis shows clinical relevance for adjuvant decisions\n")
cat("══════════════════════════════════════════════════════════════════\n")
