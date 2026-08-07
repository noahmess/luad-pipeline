# =============================================================================
# 00_setup.R
# Global configuration for the LUAD nuclear-envelope pipeline.
# Sourced first by every other script (and by run_all.R).
#
# It does three things:
#   1. Fixes the random seed so jitter / label-repel placement is reproducible.
#   2. Declares the project paths and creates the output folders.
#   3. Declares the gene sets, colour palettes and small helper functions that
#      are shared across the analysis and figure scripts.
#
# Cohort logic (see README.md for the full rationale):
#   * n = 540  -> expression-only analyses (DEG, co-expression, immune, ROC).
#   * n = 394  -> survival analyses (Cox, Score A, Kaplan-Meier, Schoenfeld),
#                 i.e. patients with complete time/status/age/sex/stage.
# =============================================================================

# ---- 1. Reproducibility ------------------------------------------------------
# A single global seed. Stochastic steps (geom_jitter, ggrepel) call set.seed()
# again with this value immediately before they run, because some ggplot draw
# steps consume the RNG lazily.
SEED <- 42L
set.seed(SEED)

# ---- 2. Paths ----------------------------------------------------------------
# All paths are relative to the project root so the repository is portable.
# Run R from the project root, or set PROJECT_ROOT before sourcing.
if (!exists("PROJECT_ROOT")) PROJECT_ROOT <- getwd()

PATHS <- list(
  data    = file.path(PROJECT_ROOT, "data"),       # cached .RData intermediates
  figures = file.path(PROJECT_ROOT, "outputs", "figures"),
  tables  = file.path(PROJECT_ROOT, "outputs", "tables")
)
for (p in PATHS) dir.create(p, recursive = TRUE, showWarnings = FALSE)

# Named paths of the cached intermediates produced by the numbered scripts.
RDATA <- list(
  tcga       = file.path(PATHS$data, "TCGA_LUAD_data.RData"),        # 01
  deg        = file.path(PATHS$data, "TCGA_LUAD_DEG_results.RData"), # 02
  enrichment = file.path(PATHS$data, "TCGA_LUAD_enrichment.RData"),  # 03
  cohort     = file.path(PATHS$data, "TCGA_LUAD_cohort.RData"),      # 04
  immune     = file.path(PATHS$data, "TCGA_LUAD_immune.RData")       # 07
)

# ---- 3. Packages -------------------------------------------------------------
# Loaded quietly; install.packages()/BiocManager::install() is the user's job
# (documented in README.md). Sourcing fails early with a clear message if a
# package is missing, rather than deep inside an analysis.
.need <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing))
    stop("Missing packages: ", paste(missing, collapse = ", "),
         "\nInstall them (see README.md) before running the pipeline.",
         call. = FALSE)
  invisible(lapply(pkgs, function(p)
    suppressPackageStartupMessages(library(p, character.only = TRUE))))
}

# ---- 4. Gene sets ------------------------------------------------------------
# The 20 nuclear-envelope genes selected a priori, grouped by sub-compartment.
NE_GENES_20 <- c(
  # Nuclear lamina
  "LMNB1", "LMNB2", "LMNA", "TMPO", "LBR", "BANF1", "EMD",
  # LINC complex
  "SUN1", "SUN2", "SYNE1", "SYNE2",
  # Nuclear pore complex
  "NUP62", "NUP88", "NUP107", "NUP155", "NDC1", "POM121",
  # Inner nuclear membrane
  "TOR1AIP1", "LEMD3", "ANKLE2"
)

# The 9-gene up-regulated panel used for co-expression, immune and ROC analyses.
PANEL_9 <- c("LMNB2", "LMNB1", "TMPO", "NDC1", "NUP62",
             "BANF1", "NUP107", "NUP88", "LBR")

# The final four-gene prognostic signature (Score A).
# LMNB1 was excluded (multivariate Cox p = 0.0571, borderline).
SIGNATURE_4 <- c("LMNB2", "TMPO", "NDC1", "LBR")

# ---- 5. Palettes -------------------------------------------------------------
PAL <- list(
  surv      = c(Low = "#2C7FB8", High = "#D95F0E"),  # Kaplan-Meier groups
  deg       = c(Upregulated = "#E74C3C",
                Downregulated = "#2E86C1",
                NS = "grey60"),                       # volcano status
  nup88     = c(NUP88 = "#E24B4A", Other = "#AFA9EC"),# NUP88 vs rest
  stage     = c(I = "#4393C3", II = "#74C476",
                III = "#FD8D3C", IV = "#D73027"),
  diverging = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100)
)

# ---- 6. Helpers --------------------------------------------------------------
# Save a ggplot or recordedplot to both TIFF (300 dpi, LZW - Frontiers' format)
# and PNG (300 dpi, for review). Width/height are in inches.
save_fig <- function(plot, name, width, height) {
  tiff_path <- file.path(PATHS$figures, paste0(name, ".tiff"))
  png_path  <- file.path(PATHS$figures, paste0(name, ".png"))
  ggplot2::ggsave(tiff_path, plot, width = width, height = height,
                  dpi = 300, compression = "lzw")
  ggplot2::ggsave(png_path,  plot, width = width, height = height, dpi = 300)
  message("  saved: ", basename(tiff_path), " + ", basename(png_path))
  invisible(tiff_path)
}

# Collapse a TCGA AJCC stage string ("Stage IIIA") to a coarse factor I-IV.
collapse_stage <- function(x) {
  s <- gsub("Stage ", "", x)
  s <- gsub("[ABC]", "", s)
  s[!s %in% c("I", "II", "III", "IV")] <- NA
  factor(s, levels = c("I", "II", "III", "IV"))
}

message("00_setup.R loaded | seed = ", SEED, " | root = ", PROJECT_ROOT)
