# =============================================================================
# 06_external_validation.R
# External validation of the four-gene Score A in two independent microarray
# cohorts cited in the Methods: GSE31210 (n = 226 used) and GSE50081 (n = 181).
#
# >>> RECONSTRUCTED MODULE <<<
# The original GEO download/processing code was NOT among the provided scripts;
# only the resulting objects (zscore_geo, risk_geo, surv_geo, ...) survived in
# the session. This file reconstructs the most likely pipeline from the
# manuscript Methods and the structure of those objects. The phenotype-column
# mapping and probe-to-gene selection MUST be checked against the original
# analysis before the numbers are quoted - GEO series label survival fields
# inconsistently. Treat the output as provisional until verified.
#
# Output: outputs/tables/External_validation_GEO.csv (Score A log-rank p / HR).
# =============================================================================

if (!exists("RDATA")) source("R/00_setup.R")
.need(c("GEOquery", "survival"))
set.seed(SEED)

# Best-probe-per-gene Score A on a GEO ExpressionSet: returns a numeric vector.
score_A_geo <- function(eset, gene_col = "Gene Symbol") {
  fd  <- Biobase::fData(eset)
  ex  <- Biobase::exprs(eset)
  symcol <- gene_col
  if (!symcol %in% colnames(fd))                       # fall back to any symbol-like column
    symcol <- grep("symbol", colnames(fd), ignore.case = TRUE, value = TRUE)[1]
  z <- sapply(SIGNATURE_4, function(g) {
    probes <- rownames(fd)[fd[[symcol]] == g]
    probes <- probes[probes %in% rownames(ex)]
    if (!length(probes)) return(rep(NA_real_, ncol(ex)))
    # pick the probe with the highest variance (most informative)
    best <- probes[which.max(apply(ex[probes, , drop = FALSE], 1, var))]
    scale(ex[best, ])[, 1]
  })
  rowMeans(z)
}

# Run a Score-A median-split survival test; returns one result row.
validate <- function(score, time, status, cohort) {
  ok  <- is.finite(score) & is.finite(time) & time > 0 & !is.na(status)
  grp <- factor(ifelse(score[ok] > median(score[ok]), "High", "Low"),
                levels = c("Low", "High"))
  so  <- survival::Surv(time[ok], status[ok])
  hr  <- summary(coxph(so ~ score[ok]))$conf.int
  lr  <- survdiff(so ~ grp)
  data.frame(Cohort = cohort, n = sum(ok),
             HR = round(hr[1, 1], 3),
             p_logrank = signif(1 - pchisq(lr$chisq, 1), 3))
}

res <- list()

# --- GSE31210 ---------------------------------------------------------------
# Phenotype fields (from the observed object): "days before death/censor:ch1"
# and "death:ch1" (dead/alive).
g1 <- tryCatch(getGEO("GSE31210", GSEMatrix = TRUE)[[1]], error = function(e) NULL)
if (!is.null(g1)) {
  ph     <- Biobase::pData(g1)
  status <- ifelse(ph[["death:ch1"]] == "dead", 1, 0)
  time   <- as.numeric(ph[["days before death/censor:ch1"]]) / 30.44
  res$GSE31210 <- validate(score_A_geo(g1), time, status, "GSE31210")
}

# --- GSE50081 ---------------------------------------------------------------
# Phenotype fields: "survival time:ch1" (years) and "status:ch1" (dead/alive).
g2 <- tryCatch(getGEO("GSE50081", GSEMatrix = TRUE)[[1]], error = function(e) NULL)
if (!is.null(g2)) {
  ph     <- Biobase::pData(g2)
  status <- ifelse(ph[["status:ch1"]] == "dead", 1, 0)
  time   <- as.numeric(ph[["survival time:ch1"]]) * 12        # years -> months
  res$GSE50081 <- validate(score_A_geo(g2), time, status, "GSE50081")
}

if (length(res)) {
  out <- do.call(rbind, res)
  write.csv(out, file.path(PATHS$tables, "External_validation_GEO.csv"),
            row.names = FALSE)
  print(out)
  message("06: external validation written (RECONSTRUCTED - verify).")
} else {
  message("06: GEO download failed or unavailable; skipped.")
}
