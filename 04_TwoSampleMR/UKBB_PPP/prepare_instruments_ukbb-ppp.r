# ==============================================================================
# UKBB-PPP cis-pQTL instrument preparation
#
# Prepares genome-wide significant, LD-independent cis-pQTL instruments from
# non-clumped UKBB-PPP cis-pQTL files.
#
# Instrument criteria:
#   P < 5e-8
#   LD clumping r2 < 0.001
#   LD clumping window = 10,000 kb
#
# LD clumping is performed locally using a European reference panel.
#
# ==============================================================================


library(data.table)
library(ieugwasr)
library(genetics.binaRies)


# ------------------------------------------------------------------------------
# User-defined paths
# ------------------------------------------------------------------------------

# Directory containing non-clumped UKBB-PPP cis-pQTL files
input_dir <- "/path/to/UKBB_PPP/cis_IVs_rsids"

# Directory for LD-clumped instruments
output_dir <- "/path/to/UKBB_PPP/cis_IVs_rsids_LDclmpd"

# PLINK-format European LD reference prefix (.bed/.bim/.fam)
eur_bfile <- "/path/to/EUR"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------------------------
# Instrument-selection parameters
# ------------------------------------------------------------------------------

P_THRESHOLD <- 5e-8
CLUMP_R2 <- 0.001
CLUMP_KB <- 10000


# ------------------------------------------------------------------------------
# Identify UKBB-PPP cis-pQTL files
# ------------------------------------------------------------------------------

pqtl_files <- list.files(
  input_dir,
  full.names = TRUE
)

if (length(pqtl_files) == 0) {
  stop("No UKBB-PPP cis-pQTL files found in input_dir.")
}

message(
  "UKBB-PPP files identified: ",
  length(pqtl_files)
)


# ==============================================================================
# Prepare and LD-clump each cis-pQTL exposure
# ==============================================================================

for (i in seq_along(pqtl_files)) {

  input_file <- pqtl_files[i]

  exposure_name <- tools::file_path_sans_ext(
    basename(input_file)
  )

  message("")
  message(
    "Processing ",
    exposure_name,
    " (",
    i,
    "/",
    length(pqtl_files),
    ")"
  )


  # ---------------------------------------------------------------------------
  # Load non-clumped cis-pQTL data
  # ---------------------------------------------------------------------------

  exposure_data <- fread(
    input_file,
    data.table = FALSE
  )

  if (!all(c("SNP", "log10p") %in% names(exposure_data))) {

    warning(
      "Required columns SNP/log10p absent from ",
      basename(input_file),
      ". Skipping."
    )

    next
  }


  # ---------------------------------------------------------------------------
  # Convert -log10(P) to P-value
  # ---------------------------------------------------------------------------

  exposure_data$pval <- 10^(
    -as.numeric(exposure_data$log10p)
  )


  # ---------------------------------------------------------------------------
  # Restrict to genome-wide significant cis-pQTLs
  # ---------------------------------------------------------------------------

  exposure_data <- exposure_data[
    !is.na(exposure_data$pval) &
      exposure_data$pval < P_THRESHOLD,
    ,
    drop = FALSE
  ]

  message(
    "  Genome-wide significant cis-pQTLs: ",
    nrow(exposure_data)
  )

  if (nrow(exposure_data) == 0) {

    message(
      "  No variants pass P < ",
      P_THRESHOLD,
      ". Skipping."
    )

    next
  }


  # ---------------------------------------------------------------------------
  # Construct ieugwasr clumping input
  # ---------------------------------------------------------------------------

  clump_input <- unique(
    data.frame(
      rsid = exposure_data$SNP,
      pval = exposure_data$pval,
      id = exposure_name,
      stringsAsFactors = FALSE
    )
  )


  # ---------------------------------------------------------------------------
  # LD clumping
  # ---------------------------------------------------------------------------

  clumped <- ld_clump(
    d = clump_input,
    clump_kb = CLUMP_KB,
    clump_r2 = CLUMP_R2,
    clump_p = P_THRESHOLD,
    plink_bin = genetics.binaRies::get_plink_binary(),
    bfile = eur_bfile
  )

  message(
    "  LD-independent instruments retained: ",
    nrow(clumped)
  )

  if (nrow(clumped) == 0) {
    next
  }


  # ---------------------------------------------------------------------------
  # Recover complete UKBB-PPP records for retained instruments
  # ---------------------------------------------------------------------------

  exposure_clumped <- exposure_data[
    exposure_data$SNP %in% clumped$rsid,
    ,
    drop = FALSE
  ]

  # Preserve the order returned by LD clumping
  exposure_clumped <- exposure_clumped[
    match(
      clumped$rsid,
      exposure_clumped$SNP
    ),
    ,
    drop = FALSE
  ]


  # ---------------------------------------------------------------------------
  # Save exposure-specific instrument file
  # ---------------------------------------------------------------------------

  output_file <- file.path(
    output_dir,
    paste0(
      exposure_name,
      "_LDclumped.txt"
    )
  )

  fwrite(
    exposure_clumped,
    file = output_file,
    sep = "\t",
    quote = FALSE,
    na = "NA"
  )
}


# ==============================================================================
# Complete
# ==============================================================================

message("")
message("UKBB-PPP instrument preparation complete.")
message("P threshold: ", P_THRESHOLD)
message("LD r2 threshold: ", CLUMP_R2)
message("LD window: ", CLUMP_KB, " kb")
message("Output directory: ", output_dir)