# =============================================================================
# 10_signatures.R
# Definitions of all prognostic signatures compared in the benchmark.
#
# Five signatures are compared in the revision:
#   - Score A         : the four-gene coefficient-free score from this paper
#                       (LMNB2 + TMPO + NDC1 + LBR, mean of Z-scores)
#   - LASSO-Cox       : a regularised LASSO-Cox model refitted on TCGA train,
#                       glmnet alpha = 1, lambda.min, 10-fold CV
#   - Li2023_TME6     : Li et al. 2023, 6-gene tumour microenvironment signature
#                       (Sci Rep, DOI 10.1038/s41598-023-40980-2)
#   - Zhou2022        : Zhou et al. 2022, 3-gene signature
#                       (Transl Lung Cancer Res, DOI 10.21037/tlcr-22-444)
#   - CCP_like        : a 31-gene cell-cycle proliferation signature in the
#                       spirit of Cuzick 2011 and Wistuba 2013, equal-weight
#                       Z-score, used as a proliferation benchmark
#
# Each signature is a list with four fields:
#   genes   : character vector of HGNC gene symbols
#   coefs   : named numeric vector of weights aligned to genes
#   zscore  : logical, TRUE = compute mean Z-score, FALSE = use coefs directly
#   refit   : logical, TRUE = refit on TCGA training partition (LASSO only)
#
# Aliases. Some old gene symbols may appear in microarray annotations:
#   CDK1  <-> CDC2
#   PCLAF <-> KIAA0101
#   ORC6  <-> ORC6L
#   RSRC2 <-> C18orf24
# Resolution is handled in 11_benchmark.R when matching to platform symbols.
# =============================================================================

signatures <- list()

# -----------------------------------------------------------------------------
# Score A: the four-gene coefficient-free score from this paper
# -----------------------------------------------------------------------------
GENES_SCORE_A <- c("LMNB2", "TMPO", "NDC1", "LBR")

signatures$Score_A <- list(
    genes  = GENES_SCORE_A,
    coefs  = setNames(rep(1, length(GENES_SCORE_A)), GENES_SCORE_A),
    zscore = TRUE,
    refit  = FALSE
)

# -----------------------------------------------------------------------------
# LASSO-Cox: refitted on TCGA train partition (handled in 11_benchmark.R)
# Placeholder definition. Genes and coefs are filled at fit time.
# -----------------------------------------------------------------------------
signatures$LASSO_Cox <- list(
    genes  = character(0),
    coefs  = numeric(0),
    zscore = FALSE,
    refit  = TRUE
)

# -----------------------------------------------------------------------------
# Li 2023, 6-gene tumour microenvironment signature
# Coefficients from the published Cox model.
# -----------------------------------------------------------------------------
signatures$Li2023_TME6 <- list(
    genes  = c("PLK1", "LDHA", "FURIN", "FSCN1", "RAB27B", "MS4A1"),
    coefs  = c(
        PLK1   = 0.1835,
        LDHA   = 0.3800,
        FURIN  = 0.1700,
        FSCN1  = 0.1600,
        RAB27B = 0.1900,
        MS4A1  = -0.1400
    ),
    zscore = FALSE,
    refit  = FALSE
)

# -----------------------------------------------------------------------------
# Zhou 2022, 3-gene signature
# -----------------------------------------------------------------------------
signatures$Zhou2022 <- list(
    genes  = c("COL1A1", "GPX3", "PLEK2"),
    coefs  = c(
        COL1A1 = 0.1446053,
        GPX3   = -0.2426827,
        PLEK2  = 0.2697514
    ),
    zscore = FALSE,
    refit  = FALSE
)

# -----------------------------------------------------------------------------
# CCP-like, 31-gene cell-cycle proliferation panel
# Used as an equal-weight Z-score benchmark for proliferation-driven risk.
# -----------------------------------------------------------------------------
CCP_GENES <- c(
    "FOXM1",  "CDC20",  "CDKN3",  "CDK1",   "KIF11",  "PCLAF",  "NUSAP1",
    "CENPF",  "ASPM",   "BUB1B",  "RRM2",   "DLGAP5", "BIRC5",  "KIF20A",
    "PLK1",   "TOP2A",  "TK1",    "PBK",    "ASF1B",  "RSRC2",  "RAD54L",
    "PTTG1",  "CDCA3",  "MCM10",  "PRC1",   "DTL",    "CEP55",  "RAD51",
    "CENPM",  "CDCA8",  "ORC6"
)

signatures$CCP_like <- list(
    genes  = CCP_GENES,
    coefs  = setNames(rep(1, length(CCP_GENES)), CCP_GENES),
    zscore = TRUE,
    refit  = FALSE
)

# -----------------------------------------------------------------------------
# Alias map for old microarray annotations
# -----------------------------------------------------------------------------
SIGNATURE_ALIAS_MAP <- list(
    CDK1  = c("CDC2"),
    PCLAF = c("KIAA0101"),
    RSRC2 = c("C18orf24"),
    ORC6  = c("ORC6L")
)

message("10: signatures defined - ",
        paste(names(signatures), collapse = ", "))
