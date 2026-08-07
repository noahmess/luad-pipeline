# =============================================================================
# run_revision.R
# Master script for the revision analyses (benchmark + sensitivity).
#
# This script runs the three revision modules in order:
#   1. 10_signatures.R       : signature definitions
#   2. 11_benchmark.R        : five-signature benchmark on four cohorts
#   3. 12_sensitivity_LOO.R  : leave-one-gene-out + nested LRT
#
# The original pipeline (run_all.R) is not re-run here. The benchmark and
# sensitivity modules read the harmonised cohort from
# data/TCGA_LUAD_cohort.RData, which is produced by R/04_cohort_harmonization.R.
# If that file is missing, run run_all.R first.
#
# Outputs:
#   outputs/tables/Benchmark_main_3cohorts.csv
#   outputs/tables/Benchmark_suppl_TCGAfull.csv
#   outputs/tables/Supplementary_Table_S3_LOO.csv
#   BENCHMARK_ENV.RData
#
# Approximate runtime:
#   - first run, GEO downloads needed : 5 to 10 minutes
#   - subsequent runs, GEO cached     : under 2 minutes
# =============================================================================

setwd("~/Desktop/luad-pipeline")

t_start <- Sys.time()
message("==========================================")
message("LUAD revision pipeline started at ", format(t_start))
message("==========================================")

# Sanity check: the harmonised cohort must exist
if (!file.exists("data/TCGA_LUAD_cohort.RData")) {
    stop("data/TCGA_LUAD_cohort.RData is missing. Run run_all.R first.")
}

# Step 1: signatures
message("\n--- Step 1: loading signature definitions ---")
source("R/10_signatures.R")

# Step 2: benchmark
message("\n--- Step 2: benchmark of five signatures on four cohorts ---")
source("R/11_benchmark.R")

# Step 3: sensitivity LOO + nested LRT
message("\n--- Step 3: leave-one-gene-out sensitivity + nested LRT ---")
source("R/12_sensitivity_LOO.R")

t_end <- Sys.time()
elapsed <- round(as.numeric(difftime(t_end, t_start, units = "mins")), 2)

message("\n==========================================")
message("Revision pipeline finished at ", format(t_end))
message("Elapsed time: ", elapsed, " minutes")
message("==========================================")
message("Outputs:")
message("  outputs/tables/Benchmark_main_3cohorts.csv")
message("  outputs/tables/Benchmark_suppl_TCGAfull.csv")
message("  outputs/tables/Supplementary_Table_S3_LOO.csv")
message("  BENCHMARK_ENV.RData")
