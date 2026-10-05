# ==============================================================================
# Build UKBB-PPP SNP-to-protein lookup
#
# Generates a SNP tagger from the full, non-LD-clumped UKBB-PPP cis-pQTL data.
# Each SNP is linked to the protein/exposure represented by the corresponding
# cis-pQTL file.
#
# This lookup is subsequently used to map SNPs identified by cross-trait
# SuSiE colocalisation to candidate protein exposures.
# ==============================================================================


library(data.table)


# ------------------------------------------------------------------------------
# User-defined paths
# ------------------------------------------------------------------------------

# Directory containing full, non-LD-clumped UKBB-PPP cis-pQTL files
input_dir <- "/path/to/UKBB_PPP/cis_IVs_rsids"

# Output file
output_file <- "/path/to/output/UKBB_PPP_snp_tagger.rds"


# ==============================================================================
# 1. Identify UKBB-PPP cis-pQTL files
# ==============================================================================

pqtl_files <- list.files(
  input_dir,
  full.names = TRUE
)

if (length(pqtl_files) == 0) {
  stop("No UKBB-PPP cis-pQTL files found in input_dir.")
}

message(
  "UKBB-PPP cis-pQTL files identified: ",
  length(pqtl_files)
)


# ==============================================================================
# 2. Extract SNP-to-exposure mappings
# ==============================================================================

snp_tagger_list <- vector(
  mode = "list",
  length = length(pqtl_files)
)

for (i in seq_along(pqtl_files)) {

  input_file <- pqtl_files[i]

  # Exposure/protein identifier is encoded in the source filename
  exposure <- tools::file_path_sans_ext(
    basename(input_file)
  )

  message(
    "Processing ",
    exposure,
    " (",
    i,
    "/",
    length(pqtl_files),
    ")"
  )


  # ---------------------------------------------------------------------------
  # Load non-clumped cis-pQTL data
  # ---------------------------------------------------------------------------

  pqtl_data <- fread(
    input_file,
    select = "SNP",
    data.table = FALSE
  )

  if (!"SNP" %in% names(pqtl_data)) {

    warning(
      "SNP column absent from ",
      basename(input_file),
      ". Skipping."
    )

    next
  }


  # ---------------------------------------------------------------------------
  # Associate each SNP with its molecular exposure
  # ---------------------------------------------------------------------------

  snp_tagger_list[[i]] <- data.frame(
    SNP = pqtl_data$SNP,
    exposure = exposure,
    stringsAsFactors = FALSE
  )
}


# ==============================================================================
# 3. Combine SNP-to-exposure mappings
# ==============================================================================

snp_tagger <- rbindlist(
  snp_tagger_list,
  use.names = TRUE,
  fill = TRUE
)

snp_tagger <- unique(
  snp_tagger,
  by = c("SNP", "exposure")
)

message("")
message(
  "Unique SNP-exposure mappings: ",
  nrow(snp_tagger)
)

message(
  "Unique SNPs: ",
  uniqueN(snp_tagger$SNP)
)

message(
  "Unique exposures: ",
  uniqueN(snp_tagger$exposure)
)


# ==============================================================================
# 4. Save full UKBB-PPP SNP tagger
# ==============================================================================

saveRDS(
  snp_tagger,
  file = output_file
)


# ==============================================================================
# Complete
# ==============================================================================

message("")
message("UKBB-PPP SNP tagger complete.")
message("Output: ", output_file)