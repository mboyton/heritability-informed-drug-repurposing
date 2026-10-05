# ==============================================================================
# Build eQTLGen SNP-to-gene lookup
#
# Generates a SNP tagger from the full, non-LD-clumped eQTLGen cis-eQTL
# dataset. Each SNP is linked to its corresponding gene/exposure.
#
# This lookup is subsequently used to map SNPs identified by cross-trait
# SuSiE colocalisation to candidate gene-expression exposures.
#
# ==============================================================================


library(data.table)


# ------------------------------------------------------------------------------
# User-defined paths
# ------------------------------------------------------------------------------

# Full, non-LD-clumped eQTLGen cis-eQTL dataset
eqtlgen_file <- "/path/to/eqtlgen_maf_annotated.txt"

# Output file
output_file <- "/path/to/output/eQTLGen_snp_tagger.rds"


# ==============================================================================
# 1. Load full eQTLGen cis-eQTL dataset
# ==============================================================================

eqtlgen_cis <- fread(
  eqtlgen_file,
  select = c("SNP", "GeneSymbol")
)

required_columns <- c(
  "SNP",
  "GeneSymbol"
)

missing_columns <- setdiff(
  required_columns,
  names(eqtlgen_cis)
)

if (length(missing_columns) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}


# ==============================================================================
# 2. Construct SNP-to-gene lookup
# ==============================================================================

snp_tagger <- eqtlgen_cis[
  !is.na(SNP) &
    !is.na(GeneSymbol) &
    SNP != "" &
    GeneSymbol != "",
  .(
    SNP,
    exposure = GeneSymbol
  )
]


# Remove exact duplicate SNP-gene mappings
snp_tagger <- unique(
  snp_tagger,
  by = c("SNP", "exposure")
)


# ==============================================================================
# 3. Summarise lookup
# ==============================================================================

message(
  "Unique SNP-gene mappings: ",
  nrow(snp_tagger)
)

message(
  "Unique SNPs: ",
  uniqueN(snp_tagger$SNP)
)

message(
  "Unique genes: ",
  uniqueN(snp_tagger$exposure)
)


# ==============================================================================
# 4. Save full eQTLGen SNP tagger
# ==============================================================================

saveRDS(
  snp_tagger,
  file = output_file
)


# ==============================================================================
# Complete
# ==============================================================================

message("")
message("eQTLGen SNP tagger complete.")
message("Output: ", output_file)