# ==============================================================================
# Two-sample Mendelian randomisation: eQTLGen cis-eQTLs
#
# Runs two-sample MR between LD-clumped eQTLGen cis-eQTL instruments and a
# specified outcome GWAS.
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

# Directory containing gene-specific LD-clumped eQTLGen instruments
exposure_dir <- "/path/to/eqtlgen_clumped"

# Outcome GWAS
outcome_file <- "/path/to/outcome_gwas.txt"

# Outcome identifier
outcome_name <- "outcome_name"

# Output directory
output_dir <- "/path/to/output"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ==============================================================================
# 1. Load and format outcome GWAS
# ==============================================================================

outcome_raw <- fread(
  outcome_file,
  data.table = FALSE
)

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
# 2. Identify eQTLGen exposure files
# ==============================================================================

exposure_files <- list.files(
  exposure_dir,
  pattern = "_LDclumped\\.txt$",
  full.names = TRUE
)

if (length(exposure_files) == 0) {
  stop("No LD-clumped eQTLGen exposure files found.")
}

message(
  "eQTLGen exposures identified: ",
  length(exposure_files)
)


# ==============================================================================
# 3. Run two-sample MR
# ==============================================================================

mr_results <- vector(
  mode = "list",
  length = length(exposure_files)
)

for (i in seq_along(exposure_files)) {

  exposure_file <- exposure_files[i]

  gene <- sub(
    "_LDclumped\\.txt$",
    "",
    basename(exposure_file)
  )

  message("")
  message(
    "Processing ",
    gene,
    " (",
    i,
    "/",
    length(exposure_files),
    ")"
  )


  # ---------------------------------------------------------------------------
  # Load eQTLGen instruments
  # ---------------------------------------------------------------------------

  exposure_raw <- fread(
    exposure_file,
    data.table = FALSE
  )


  # ---------------------------------------------------------------------------
  # Format exposure
  #
  # IMPORTANT:
  # Replace the mappings below with the exact column names in the source
  # eQTLGen dataset if they differ.
  # ---------------------------------------------------------------------------

  exposure_data <- format_data(
    exposure_raw,
    type = "exposure",
    snp_col = "SNP",
    beta_col = "Zscore",
    effect_allele_col = "AlleleAssessed",
    other_allele_col = "OtherAllele",
    eaf_col = "MAF",
    pval_col = "Pvalue"
  )

  exposure_data$exposure <- gene
  exposure_data$id.exposure <- gene


  # ---------------------------------------------------------------------------
  # Restrict to instruments present in outcome
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
  # Retain IVW result
  # ---------------------------------------------------------------------------

  ivw_output <- mr_output[
    mr_output$method == "Inverse variance weighted",
    ,
    drop = FALSE
  ]

  if (nrow(ivw_output) == 0) {

    message(
      "  No IVW estimate available for ",
      gene,
      "."
    )

    next
  }


  mr_results[[i]] <- data.frame(
    exposure = gene,
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
# 4. Combine results
# ==============================================================================

mr_results <- rbindlist(
  mr_results,
  use.names = TRUE,
  fill = TRUE
)


# ==============================================================================
# 5. Save results
# ==============================================================================

timestamp <- format(
  Sys.time(),
  "%Y%m%d_%H%M%S"
)

output_file <- file.path(
  output_dir,
  paste0(
    "eQTLGen_",
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
message("eQTLGen two-sample MR complete.")
message("Outcome: ", outcome_name)
message("Exposures analysed: ", nrow(mr_results))
message("Results: ", output_file)