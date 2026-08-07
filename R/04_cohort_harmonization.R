# =============================================================================
# 04_cohort_harmonization.R
# The single source of truth for the two cohorts used in the paper.
#
#   EXPRESSION COHORT  (n = 540 tumours)
#       Size-factor-normalised, log2 expression of the tumour samples.
#       Used by every analysis that does NOT need survival covariates:
#       co-expression (Fig 4), immune correlation (Fig 8), PD-L1 (Fig 9),
#       Ki-67 (Fig S1). Kept at full size to maximise statistical power.
#
#   SURVIVAL COHORT    (n = 394 tumours, 129 events)
#       The subset of tumours with COMPLETE survival time + status + age +
#       sex + stage. Fixed once here so univariate and multivariate Cox run
#       on an identical sample (no n drift between models).
#
# Output: data/TCGA_LUAD_cohort.RData containing
#   norm_tumor                         : normalised counts, 540 tumour columns
#   expr9                              : 540 x 9   (PANEL_9), log2
#   expr_all                           : 540 x (PANEL_9 + 18 immune markers)
#   expr4_v, expr9_h                   : 394 x 4 and 394 x 9 (survival cohort)
#   time_v, status_v, age_v, sex_v, stage_v, surv_obj : 394-length covariates
# =============================================================================

if (!exists("RDATA")) source("R/00_setup.R")
.need(c("DESeq2", "SummarizedExperiment", "survival"))

if (file.exists(RDATA$cohort)) {
  message("04: cache found (", basename(RDATA$cohort), ") - skipping.")
} else {
  load(RDATA$tcga)  # data, gene_symbols, tumor_idx, normal_idx

  # ---------------------------------------------------------------------------
  # A. Normalised tumour expression (the n = 540 side)
  # ---------------------------------------------------------------------------
  counts <- SummarizedExperiment::assay(data, "unstranded")
  dds_t  <- DESeqDataSetFromMatrix(
    countData = counts[, tumor_idx],
    colData   = data.frame(row.names = colnames(counts)[tumor_idx],
                           cond = rep("T", length(tumor_idx))),
    design    = ~ 1
  )
  dds_t      <- estimateSizeFactors(dds_t)
  norm_tumor <- counts(dds_t, normalized = TRUE)        # genes x 540

  # Helper: log2 expression matrix (samples x genes) for a set of symbols.
  expr_of <- function(symbols) {
    idx <- match(symbols, gene_symbols)
    stopifnot(!anyNA(idx))
    m <- t(log2(norm_tumor[idx, , drop = FALSE] + 1))
    colnames(m) <- symbols
    m
  }

  expr9 <- expr_of(PANEL_9)                              # 540 x 9

  # Immune / proliferation markers (used by 07_immune and Figs 8, 9, S1).
  IMMUNE_MARKERS <- c("CD8A", "CD4", "FOXP3", "CD68", "CD163", "NOS2",
                      "CD19", "NCAM1", "CD274", "PDCD1", "CTLA4", "HAVCR2",
                      "LAG3", "IFNG", "GZMA", "PRF1", "MKI67", "PECAM1")
  found_markers <- intersect(IMMUNE_MARKERS, gene_symbols)
  expr_all <- expr_of(c(PANEL_9, found_markers))         # 540 x (9 + markers)

  # ---------------------------------------------------------------------------
  # B. Harmonised survival cohort (the n = 394 side)
  # ---------------------------------------------------------------------------
  clin <- as.data.frame(SummarizedExperiment::colData(data))[tumor_idx, ]

  # Overall survival in months: time to death, or last follow-up if alive.
  surv_time   <- as.numeric(clin$days_to_death)
  surv_status <- ifelse(clin$vital_status == "Dead", 1, 0)
  alive       <- which(is.na(surv_time) | surv_status == 0)
  surv_time[alive] <- as.numeric(clin$days_to_last_follow_up)[alive]
  surv_time_months <- surv_time / 30.44

  age   <- as.numeric(clin$age_at_index)
  sex   <- factor(clin$gender)
  stage <- collapse_stage(clin$ajcc_pathologic_stage)   # helper from 00_setup

  # Complete-case filter -> the fixed 394-patient survival cohort.
  valid <- which(surv_time_months > 0 & !is.na(surv_time_months) &
                 !is.na(surv_status) & !is.na(age) &
                 !is.na(sex) & !is.na(stage))

  expr4_v  <- expr9[valid, SIGNATURE_4, drop = FALSE]    # 394 x 4
  expr9_h  <- expr9[valid, , drop = FALSE]               # 394 x 9
  time_v   <- surv_time_months[valid]
  status_v <- surv_status[valid]
  age_v    <- age[valid]
  sex_v    <- sex[valid]
  stage_v  <- stage[valid]
  surv_obj <- survival::Surv(time_v, status_v)

  message("04: expression cohort n = ", nrow(expr9),
          " | survival cohort n = ", length(time_v),
          " (", sum(status_v), " events)")

  save(norm_tumor, expr9, expr_all, found_markers,
       expr4_v, expr9_h, time_v, status_v, age_v, sex_v, stage_v, surv_obj,
       file = RDATA$cohort)
  message("04: saved ", basename(RDATA$cohort))
}
