# ==============================================================================
# Map cross-trait SuSiE SNPs to UKBB-PPP protein exposures
#
# Extracts SNPs identified as hit1 or hit2 by cross-trait SuSiE
# colocalisation and queries them against the full UKBB-PPP SNP-to-protein
# lookup.
#
# The resulting proteins are candidate molecular exposures for downstream
# Mendelian randomisation.
#
# ==============================================================================


library(data.table)


# ------------------------------------------------------------------------------
# User-defined paths
# ------------------------------------------------------------------------------

# Merged cross-trait SuSiE colocalisation results
coloc_file <- "/path/to/merged_coloc_susie_results.rds"

# Full UKBB-PPP SNP-to-protein lookup
snp_tagger_file <- "/path/to/UKBB_PPP_snp_tagger.rds"

# Output file
output_file <- "/path/to/output/UKBB_PPP_coloc_exposures.rds"


# ==============================================================================
# 1. Load cross-trait SuSiE results
# ==============================================================================

coloc_results <- readRDS(
  coloc_file
)

coloc_results <- as.data.table(
  coloc_results
)


# ==============================================================================
# 2. Extract SuSiE hit SNPs
# ==============================================================================

required_columns <- c(
  "hit1",
  "hit2"
)

missing_columns <- setdiff(
  required_columns,
  names(coloc_results)
)

if (length(missing_columns) > 0) {
  stop(
    "Missing SuSiE columns: ",
    paste(missing_columns, collapse = ", ")
  )
}


susie_snps <- unique(
  c(
    coloc_results$hit1,
    coloc_results$hit2
  )
)

susie_snps <- susie_snps[
  !is.na(susie_snps) &
    susie_snps != ""
]

message(
  "Unique SuSiE hit1/hit2 SNPs: ",
  length(susie_snps)
)


# ==============================================================================
# 3. Load full UKBB-PPP SNP tagger
# ==============================================================================

snp_tagger <- readRDS(
  snp_tagger_file
)

snp_tagger <- as.data.table(
  snp_tagger
)

required_columns <- c(
  "SNP",
  "exposure"
)

missing_columns <- setdiff(
  required_columns,
  names(snp_tagger)
)

if (length(missing_columns) > 0) {
  stop(
    "Missing SNP-tagger columns: ",
    paste(missing_columns, collapse = ", ")
  )
}


# ==============================================================================
# 4. Map SuSiE SNPs to UKBB-PPP protein exposures
# ==============================================================================

coloc_exposures <- snp_tagger[
  SNP %in% susie_snps
]

coloc_exposures <- unique(
  coloc_exposures,
  by = c("SNP", "exposure")
)


# ==============================================================================
# 5. Summarise mapping
# ==============================================================================

message(
  "SuSiE SNPs represented in UKBB-PPP: ",
  uniqueN(coloc_exposures$SNP)
)

message(
  "Candidate UKBB-PPP protein exposures: ",
  uniqueN(coloc_exposures$exposure)
)


# ==============================================================================
# 6. Save SNP-to-protein mappings
# ==============================================================================

saveRDS(
  coloc_exposures,
  file = output_file
)


# ==============================================================================
# Complete
# ==============================================================================

message("")
message("UKBB-PPP SuSiE SNP mapping complete.")
message("Output: ", output_file)