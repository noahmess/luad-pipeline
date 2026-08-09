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

  save(data, gene_symbols, tumor_idx, normal_idx, sample_type,
       file = RDATA$tcga)
  message("01: saved ", basename(RDATA$tcga))
}
