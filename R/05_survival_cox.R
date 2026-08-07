# =============================================================================
# 05_survival_cox.R
# All survival modelling, on the fixed n = 394 cohort built in step 04.
#
#   1. Univariate Cox for each of the 9 panel genes (standardised expression).
#   2. Multivariate Cox for each gene, adjusted for age, sex and stage.
#   3. Score A = mean Z-score of LMNB2, TMPO, NDC1, LBR:
#        continuous HR, adjusted HR, median split + log-rank + 5-yr OS.
#   4. Per-gene median-split 5-yr overall survival (the four signature genes).
#   5. Schoenfeld test (proportional-hazards check) for the four genes.
#
# Exports to outputs/tables/:
#   Cox_univarie_n394.csv, Cox_multivarie_n394.csv,
#   Schoenfeld_results.csv, Comparison_Individual_vs_Signature.csv
# =============================================================================

if (!exists("RDATA")) source("R/00_setup.R")
.need("survival")
load(RDATA$cohort)  # expr4_v, expr9_h, time_v, status_v, age_v, sex_v, stage_v, surv_obj

stopifnot(length(time_v) == 394, length(status_v) == 394)
message("05: cohort = ", length(time_v), " patients, ",
        sum(status_v), " events")
set.seed(SEED)

# Helper: tidy a Cox fit row for a standardised covariate named `z`.
cox_row <- function(gene, fit) {
  s <- summary(fit)
  data.frame(Gene = gene, n = fit$n,
             HR      = round(s$coefficients["z", "exp(coef)"], 3),
             HR_low  = round(s$conf.int["z", "lower .95"], 3),
             HR_high = round(s$conf.int["z", "upper .95"], 3),
             p       = signif(s$coefficients["z", "Pr(>|z|)"], 3),
             stringsAsFactors = FALSE)
}

# ---- 1 & 2. Univariate and multivariate Cox, per panel gene -----------------
uni <- multi <- NULL
for (g in colnames(expr9_h)) {
  z <- scale(expr9_h[, g])[, 1]
  uni   <- rbind(uni,   cox_row(g, coxph(surv_obj ~ z)))
  dat   <- data.frame(z = z, age = age_v, sex = sex_v, stage = stage_v)
  multi <- rbind(multi, cox_row(g, coxph(surv_obj ~ z + age + sex + stage,
                                         data = dat)))
}
write.csv(uni,   file.path(PATHS$tables, "Cox_univarie_n394.csv"),   row.names = FALSE)
write.csv(multi, file.path(PATHS$tables, "Cox_multivarie_n394.csv"), row.names = FALSE)
message("05: univariate / multivariate Cox tables written.")

# ---- 3. Score A -------------------------------------------------------------
scoreA <- rowMeans(scale(expr4_v))                       # mean Z of the 4 genes
fitA   <- coxph(surv_obj ~ scoreA)
sA     <- summary(fitA)
datA   <- data.frame(scoreA = scoreA, age = age_v, sex = sex_v, stage = stage_v)
fitAm  <- coxph(surv_obj ~ scoreA + age + sex + stage, data = datA)
sAm    <- summary(fitAm)

grpA   <- factor(ifelse(scoreA > median(scoreA), "High", "Low"),
                 levels = c("Low", "High"))
fitKM  <- survfit(surv_obj ~ grpA)
lr     <- survdiff(surv_obj ~ grpA)
p_lr   <- 1 - pchisq(lr$chisq, df = 1)
os5A   <- summary(fitKM, times = 60)$surv               # Low, High

message(sprintf("05: Score A continuous HR = %.3f (%.3f-%.3f) p = %s",
                sA$conf.int[1, 1], sA$conf.int[1, 3], sA$conf.int[1, 4],
                signif(sA$coefficients[1, 5], 3)))
message(sprintf("05: Score A adjusted   HR = %.3f p = %s",
                sAm$conf.int["scoreA", 1],
                signif(sAm$coefficients["scoreA", 5], 3)))
message(sprintf("05: Score A 5-yr OS  Low = %.1f%%  High = %.1f%%  log-rank p = %s",
                os5A[1] * 100, os5A[2] * 100,
                formatC(p_lr, format = "e", digits = 2)))

# ---- 4. Per-gene median-split 5-yr OS (four signature genes) ----------------
comp <- data.frame(Marker = character(), HR = numeric(),
                   CI_low = numeric(), CI_high = numeric(), pvalue = numeric(),
                   Surv5y_High = numeric(), Surv5y_Low = numeric(),
                   stringsAsFactors = FALSE)
for (g in SIGNATURE_4) {
  z   <- scale(expr4_v[, g])[, 1]
  s   <- summary(coxph(surv_obj ~ z))
  grp <- factor(ifelse(expr4_v[, g] > median(expr4_v[, g]), "High", "Low"),
                levels = c("Low", "High"))
  os5 <- summary(survfit(surv_obj ~ grp), times = 60)$surv
  comp <- rbind(comp, data.frame(
    Marker = g,
    HR = round(s$conf.int[1, 1], 3),
    CI_low = round(s$conf.int[1, 3], 3),
    CI_high = round(s$conf.int[1, 4], 3),
    pvalue = signif(s$coefficients[1, 5], 3),
    Surv5y_High = round(os5[2] * 100, 1),
    Surv5y_Low  = round(os5[1] * 100, 1)))
}
# Append the composite Score A as the final row.
comp <- rbind(comp, data.frame(
  Marker = "Score A", HR = round(sA$conf.int[1, 1], 3),
  CI_low = round(sA$conf.int[1, 3], 3), CI_high = round(sA$conf.int[1, 4], 3),
  pvalue = signif(sA$coefficients[1, 5], 3),
  Surv5y_High = round(os5A[2] * 100, 1), Surv5y_Low = round(os5A[1] * 100, 1)))
write.csv(comp, file.path(PATHS$tables, "Comparison_Individual_vs_Signature.csv"),
          row.names = FALSE)

# ---- 5. Schoenfeld proportional-hazards test (four genes) -------------------
schoen <- NULL
for (g in SIGNATURE_4) {
  expr_g <- expr4_v[, g]
  zph <- cox.zph(coxph(surv_obj ~ expr_g + age_v + sex_v + stage_v))
  schoen <- rbind(schoen, data.frame(
    Gene = g,
    p_gene   = round(zph$table["expr_g", "p"], 4),
    p_global = round(zph$table["GLOBAL", "p"], 4)))
}
write.csv(schoen, file.path(PATHS$tables, "Schoenfeld_results.csv"),
          row.names = FALSE)
message("05: comparison + Schoenfeld tables written. Done.")
