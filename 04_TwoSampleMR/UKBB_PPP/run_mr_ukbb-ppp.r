# ==============================================================================
# Two-sample Mendelian randomisation: UKBB-PPP cis-pQTLs
#
# Runs two-sample MR between LD-clumped UKBB-PPP cis-pQTL instruments and
# a specified outcome GWAS.
#
# Instruments should previously have been prepared using:
#   P < 5e-8
#   LD r2 < 0.001
#   LD window = 10,000 kb
# ==============================================================================


library(TwoSampleMR)
library(data.table)


# ------------------------------------------------------------------------------
# User-defined paths
# ------------------------------------------------------------------------------

# Directory containing exposure-specific, LD-clumped UKBB-PPP instrument files
exposure_dir <- "/path/to/UKBB_PPP/cis_IVs_rsids_LDclmpd"

# Outcome GWAS
outcome_file <- "/path/to/outcome_gwas.txt"

# Name used to identify the outcome in MR results
outcome_name <- "outcome_name"

# Output directory
output_dir <- "/path/to/output"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------------------------
# Load outcome GWAS
# ------------------------------------------------------------------------------

outcome_raw <- fread(
  outcome_file,
  data.table = FALSE
)


# ==============================================================================
# Format outcome data
#
# Adjust column mappings below to match the harmonised GWAS files used in the
# analysis.
# ==============================================================================

outcome_data <- format_data(
  outcome_raw,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  eaf_col = "eaf",
  pval_col = "pval"
)

outcome_data$outcome <- outcome_name
outcome_data$id.outcome <- outcome_name


# ==============================================================================
# Identify exposure files
# ==============================================================================

exposure_files <- list.files(
  exposure_dir,
  pattern = "_LDclumped\\.txt$",
  full.names = TRUE
)

if (length(exposure_files) == 0) {
  stop("No LD-clumped UKBB-PPP exposure files found.")
}

message(
  "UKBB-PPP exposures identified: ",
  length(exposure_files)
)


# ==============================================================================
# Run two-sample MR
# ==============================================================================

mr_results <- vector(
  mode = "list",
  length = length(exposure_files)
)

for (i in seq_along(exposure_files)) {

  exposure_file <- exposure_files[i]

  exposure_name <- sub(
    "_LDclumped\\.txt$",
    "",
    basename(exposure_file)
  )

  message("")
  message(
    "Processing ",
    exposure_name,
    " (",
    i,
    "/",
    length(exposure_files),
    ")"
  )


  # ---------------------------------------------------------------------------
  # Load exposure instruments
  # ---------------------------------------------------------------------------

  exposure_raw <- fread(
    exposure_file,
    data.table = FALSE
  )


  # ---------------------------------------------------------------------------
  # Format exposure data
  #
  # These mappings correspond to the UKBB-PPP instrument data used in the
  # original workflow. Amend only if source column names differ.
  # ---------------------------------------------------------------------------

  exposure_data <- format_data(
    exposure_raw,
    type = "exposure",
    snp_col = "SNP",
    beta_col = "beta",
    se_col = "se",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    eaf_col = "eaf",
    pval_col = "pval"
  )

  exposure_data$exposure <- exposure_name
  exposure_data$id.exposure <- exposure_name


  # ---------------------------------------------------------------------------
  # Restrict outcome to exposure SNPs
  # ---------------------------------------------------------------------------

  outcome_subset <- outcome_data[
    outcome_data$SNP %in% exposure_data$SNP,
    ,
    drop = FALSE
  ]

  exposure_subset <- exposure_data[
    exposure_data$SNP %in% outcome_subset$SNP,
    ,
    drop = FALSE
  ]

  if (nrow(exposure_subset) == 0) {

    message("  No instruments present in outcome GWAS. Skipping.")
    next
  }


  # ---------------------------------------------------------------------------
  # Harmonise exposure and outcome
  # ---------------------------------------------------------------------------

  mr_data <- harmonise_data(
    exposure_dat = exposure_subset,
    outcome_dat = outcome_subset,
    action = 1
  )

  mr_data <- mr_data[
    mr_data$mr_keep,
    ,
    drop = FALSE
  ]

  if (nrow(mr_data) == 0) {

    message("  No instruments retained after harmonisation. Skipping.")
    next
  }


  # ---------------------------------------------------------------------------
  # Mendelian randomisation
  # ---------------------------------------------------------------------------

  mr_output <- mr(
    mr_data
  )

  if (is.null(mr_output) || nrow(mr_output) == 0) {
    next
  }


  # ---------------------------------------------------------------------------
  # Retain IVW result for downstream analysis
  # ---------------------------------------------------------------------------

  ivw_output <- mr_output[
    mr_output$method == "Inverse variance weighted",
    ,
    drop = FALSE
  ]

  if (nrow(ivw_output) == 0) {

    message(
      "  No IVW estimate available for ",
      exposure_name,
      "."
    )

    next
  }


  mr_results[[i]] <- data.frame(
    exposure = exposure_name,
    outcome = outcome_name,
    mr_ivw_beta = ivw_output$b[1],
    mr_ivw_se = ivw_output$se[1],
    mr_ivw_pval = ivw_output$pval[1],
    mr_ivw_nsnp = ivw_output$nsnp[1],
    stringsAsFactors = FALSE
  )

  message(
    "  IVW beta = ",
    signif(ivw_output$b[1], 4),
    "; P = ",
    signif(ivw_output$pval[1], 4),
    "; SNPs = ",
    ivw_output$nsnp[1]
  )
}


# ==============================================================================
# Combine results
# ==============================================================================

mr_results <- rbindlist(
  mr_results,
  use.names = TRUE,
  fill = TRUE
)


# ==============================================================================
# Save results
# ==============================================================================

timestamp <- format(
  Sys.time(),
  "%Y%m%d_%H%M%S"
)

output_file <- file.path(
  output_dir,
  paste0(
    "UKBB_PPP_",
    outcome_name,
    "_MR_",
    timestamp,
    ".txt"
  )
)

fwrite(
  mr_results,
  file = output_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# ==============================================================================
# Complete
# ==============================================================================

message("")
message("UKBB-PPP two-sample MR complete.")
message("Outcome: ", outcome_name)
message("Exposures analysed: ", nrow(mr_results))
message("Results: ", output_file)