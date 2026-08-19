# =============================================================================
# 01_download_tcga.R
# Download TCGA-LUAD RNA-seq + clinical data with TCGAbiolinks and assemble the
# base objects used by every downstream script.
#
# Output: data/TCGA_LUAD_data.RData containing
#   data         : SummarizedExperiment (STAR - Counts, "unstranded" assay)
#   gene_symbols : character vector of HGNC symbols, one per row of `data`
#   tumor_idx    : column indices of the primary-tumour samples (~540)
#   normal_idx   : column indices of the adjacent-normal samples (~59)
#   sample_type  : raw TCGA sample-type label for every column
#
# This is the slow, network-bound step (~2.5 GB from the GDC, 601 files).
# It is cached: if TCGA_LUAD_data.RData already exists the script exits early.
# Delete the file to force a fresh download.
#
# Known issue: GDCdownload() can fail mid-extraction with
#   "tar: Error exit delayed from previous errors." /
#   "Erreur dans if (ret == 1) break : l'argument est de longueur nulle"
# on some configurations (observed with macOS bsdtar + TCGAbiolinks 2.40.0).
# See README.md, "Reproducing the main pipeline", for the diagnosis, what was
# ruled out, and a verified workaround for downstream scripts that only need
# gene_symbols/gene_type annotation rather than the raw count matrix.
# =============================================================================

if (!exists("RDATA")) source("R/00_setup.R")

if (file.exists(RDATA$tcga)) {
  message("01: cache found (", basename(RDATA$tcga), ") - skipping download.")
} else {
  .need(c("TCGAbiolinks", "SummarizedExperiment"))

  # --- Query the GDC for harmonised STAR gene-expression counts --------------
  query <- GDCquery(
    project       = "TCGA-LUAD",
    data.category = "Transcriptome Profiling",
    data.type     = "Gene Expression Quantification",
    workflow.type = "STAR - Counts"
  )
  GDCdownload(query, directory = file.path(PATHS$data, "GDCdata"))
  data <- GDCprepare(query, directory = file.path(PATHS$data, "GDCdata"))

  # --- Gene symbols ----------------------------------------------------------
  # rowData carries the HGNC symbol under `gene_name` for STAR-Counts objects.
  gene_symbols <- as.character(SummarizedExperiment::rowData(data)$gene_name)

  # --- Sample-type split: primary tumour vs adjacent normal ------------------
  sample_type <- as.character(SummarizedExperiment::colData(data)$sample_type)
  tumor_idx  <- which(sample_type %in% c("Primary Tumor", "Primary solid Tumor"))
  normal_idx <- which(sample_type %in% c("Solid Tissue Normal"))

  message("01: ", length(tumor_idx), " tumour + ",
          length(normal_idx), " normal samples downloaded.")

  # --- Patient-level pseudo-replication check (DE cohort) --------------------
  # TCGA barcodes: the first 12 characters (TCGA-XX-XXXX) identify the patient.
  # Pseudo-replication = one patient contributing more than one sample within a
  # condition group (tumour or normal) of the tumour-vs-normal DE comparison.
  # Report only - existing objects and results are left unchanged.
  .bc    <- colnames(data)
  .pat_t <- substr(.bc[tumor_idx],  1, 12)
  .pat_n <- substr(.bc[normal_idx], 1, 12)
  .dup_t <- table(.pat_t)[table(.pat_t) > 1]
  .dup_n <- table(.pat_n)[table(.pat_n) > 1]
  message(sprintf(
    "01: DE-cohort patient check | tumour: %d samples / %d patients (%d multi-sample) | normal: %d samples / %d patients (%d multi-sample) | %d patients paired T+N",
    length(tumor_idx), length(unique(.pat_t)), length(.dup_t),
    length(normal_idx), length(unique(.pat_n)), length(.dup_n),
    length(intersect(.pat_t, .pat_n))))
  if (length(.dup_t) == 0 && length(.dup_n) == 0) {
    message("01: no pseudo-replication - each patient contributes <= 1 sample ",
            "per condition group.")
  } else {
    warning(sprintf(paste0(
      "01: PSEUDO-REPLICATION in the DE cohort - %d tumour and %d normal ",
      "patients contribute multiple samples (%d extra tumour, %d extra normal ",
      "samples). DESeq2 in R/02 treats these as independent replicates; ",
      "consider collapsing to one sample per patient if this matters."),
      length(.dup_t), length(.dup_n),
      sum(.dup_t) - length(.dup_t), sum(.dup_n) - length(.dup_n)))
  }
  # Persist the report as new output tables (existing results untouched).
  dir.create(PATHS$tables, recursive = TRUE, showWarnings = FALSE)
  write.csv(
    data.frame(group = c("tumour", "normal"),
               n_samples  = c(length(tumor_idx), length(normal_idx)),
               n_patients = c(length(unique(.pat_t)), length(unique(.pat_n))),
               n_multisample_patients = c(length(.dup_t), length(.dup_n)),
               n_extra_samples = c(sum(.dup_t) - length(.dup_t),
                                   sum(.dup_n) - length(.dup_n))),
    file.path(PATHS$tables, "Patient_duplicate_check.csv"), row.names = FALSE)
  if (length(.dup_t) > 0)
    write.csv(
      data.frame(patient = names(.dup_t), n_tumour_samples = as.integer(.dup_t)),
      file.path(PATHS$tables, "Patient_duplicate_tumour_list.csv"),
      row.names = FALSE)

  save(data, gene_symbols, tumor_idx, normal_idx, sample_type,
       file = RDATA$tcga)
  message("01: saved ", basename(RDATA$tcga))
}
