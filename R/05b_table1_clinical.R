# =============================================================================
# table1_clinical.R
# Compute Table 1: Clinical characteristics of the harmonized cohort (n = 394).
#
# Requires the pipeline objects already in memory:
#   age_v, sex_v, stage_v, time_v, status_v  (all length 394)
#
# Run after the pipeline:
#   source("/Users/nouha/Desktop/luad-pipeline/table1_clinical.R")
# =============================================================================

suppressPackageStartupMessages(library(survival))

n <- length(time_v)
stopifnot(n == 394, length(age_v)==n, length(sex_v)==n,
          length(stage_v)==n, length(status_v)==n)

# ---- Age ----
age_med   <- median(age_v, na.rm=TRUE)
age_iqr   <- quantile(age_v, c(0.25, 0.75), na.rm=TRUE)
age_range <- range(age_v, na.rm=TRUE)

# ---- Sex ----
sex_tab <- table(sex_v)
sex_pct <- round(100*prop.table(sex_tab), 1)

# ---- Stage ----
stage_tab <- table(stage_v)
stage_pct <- round(100*prop.table(stage_tab), 1)

# ---- Events ----
n_death    <- sum(status_v == 1)
n_censored <- sum(status_v == 0)

# ---- Follow-up (reverse KM = median follow-up of censored) ----
fu_med   <- median(time_v, na.rm=TRUE)
fu_range <- range(time_v, na.rm=TRUE)
# Reverse Kaplan-Meier for median potential follow-up
fu_rev <- survfit(Surv(time_v, 1 - status_v) ~ 1)
fu_rev_med <- summary(fu_rev)$table["median"]

# ---- 5-year OS ----
km <- survfit(Surv(time_v, status_v) ~ 1)
s5 <- summary(km, times = 60)  # 60 months = 5 years
os5  <- round(100*s5$surv, 1)
ci5l <- round(100*s5$lower, 1)
ci5u <- round(100*s5$upper, 1)

# ---- Build and print ----
cat("\n")
cat("=======================================================================\n")
cat("  TABLE 1: Clinical characteristics of the harmonized cohort (n=394)\n")
cat("=======================================================================\n\n")

cat(sprintf("  Age, median (IQR)           %g (%g - %g)\n", age_med, age_iqr[1], age_iqr[2]))
cat(sprintf("  Age, range                  %g - %g\n", age_range[1], age_range[2]))
cat("\n")

cat("  Sex\n")
for (s in names(sex_tab))
  cat(sprintf("    %-26s %d (%s%%)\n", s, sex_tab[s], sex_pct[s]))
cat("\n")

cat("  Tumor stage\n")
for (s in names(stage_tab))
  cat(sprintf("    Stage %-21s %d (%s%%)\n", s, stage_tab[s], stage_pct[s]))
cat("\n")

cat(sprintf("  Follow-up, median (range)   %.1f months (%.1f - %.1f)\n", fu_med, fu_range[1], fu_range[2]))
cat(sprintf("  Reverse KM median follow-up %.1f months\n", fu_rev_med))
cat(sprintf("  Events (deaths)             %d (%.1f%%)\n", n_death, 100*n_death/n))
cat(sprintf("  Censored                    %d (%.1f%%)\n", n_censored, 100*n_censored/n))
cat(sprintf("  5-year overall survival     %.1f%% (95%% CI: %.1f - %.1f)\n", os5, ci5l, ci5u))

cat("\n=======================================================================\n")

# ---- Save as CSV for easy Word import ----
out <- data.frame(
  Variable = c("Age, median (IQR)",
               "Sex, Female", "Sex, Male",
               paste0("Stage ", names(stage_tab)),
               "Follow-up, median (months)",
               "Events (deaths)", "Censored",
               "5-year OS (95% CI)"),
  Value = c(sprintf("%g (%g - %g)", age_med, age_iqr[1], age_iqr[2]),
            sprintf("%d (%s%%)", sex_tab["female"], sex_pct["female"]),
            sprintf("%d (%s%%)", sex_tab["male"], sex_pct["male"]),
            sapply(names(stage_tab), function(s)
              sprintf("%d (%s%%)", stage_tab[s], stage_pct[s])),
            sprintf("%.1f (%.1f - %.1f)", fu_med, fu_range[1], fu_range[2]),
            sprintf("%d (%.1f%%)", n_death, 100*n_death/n),
            sprintf("%d (%.1f%%)", n_censored, 100*n_censored/n),
            sprintf("%.1f%% (%.1f - %.1f)", os5, ci5l, ci5u)),
  stringsAsFactors = FALSE
)
outpath <- file.path(PROJECT_ROOT, "outputs/tables/Table1_Clinical.csv")
write.csv(out, outpath, row.names=FALSE)
cat("\nSaved:", outpath, "\n")
print(out, row.names=FALSE)
