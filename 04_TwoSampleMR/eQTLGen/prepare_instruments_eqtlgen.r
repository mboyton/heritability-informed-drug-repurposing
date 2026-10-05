# ==============================================================================
# eQTLGen cis-eQTL instrument preparation
#
# Prepares genome-wide significant, LD-independent cis-eQTL instruments from
# eQTLGen summary statistics for downstream two-sample Mendelian randomisation.
#
# Instrument criteria:
#   P < 5e-8
#   LD clumping r2 < 0.01
#   LD clumping window = 10,000 kb
#
# LD clumping is performed locally using a European reference panel.
# ==============================================================================


library(data.table)
library(ieugwasr)
library(genetics.binaRies)


# ------------------------------------------------------------------------------
# User-defined paths
# ------------------------------------------------------------------------------

# Full eQTLGen cis-eQTL dataset
eqtlgen_file <- "/path/to/eqtlgen_cis.txt"

# Output directory for gene-specific LD-clumped instruments
output_dir <- "/path/to/eqtlgen_clumped"

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
CLUMP_R2 <- 0.01
CLUMP_KB <- 10000


# ==============================================================================
# 1. Load eQTLGen cis-eQTL data
# ==============================================================================

eqtlgen <- fread(
  eqtlgen_file
)

# Expected minimum columns:
#   SNP
#   GeneSymbol
#   Pvalue
#
# Amend P-value column name below if required for the source eQTLGen file.

required_columns <- c(
  "SNP",
  "GeneSymbol",
  "Pvalue"
)

missing_columns <- setdiff(
  required_columns,
  names(eqtlgen)
)

if (length(missing_columns) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}


# ==============================================================================
# 2. Restrict to genome-wide significant cis-eQTLs
# ==============================================================================

eqtlgen <- eqtlgen[
  !is.na(Pvalue) &
    Pvalue < P_THRESHOLD
]

message(
  "Genome-wide significant cis-eQTL associations: ",
  nrow(eqtlgen)
)

if (nrow(eqtlgen) == 0) {
  stop("No eQTLGen associations pass the specified P-value threshold.")
}


# ==============================================================================
# 3. Identify genes with eligible cis-eQTL instruments
# ==============================================================================

genes <- sort(
  unique(eqtlgen$GeneSymbol)
)

genes <- genes[
  !is.na(genes) &
    genes != ""
]

message(
  "Genes with at least one eligible cis-eQTL: ",
  length(genes)
)


# ==============================================================================
# 4. LD-clump instruments separately for each gene
# ==============================================================================

for (i in seq_along(genes)) {

  gene <- genes[i]

  message("")
  message(
    "Processing ",
    gene,
    " (",
    i,
    "/",
    length(genes),
    ")"
  )


  # ---------------------------------------------------------------------------
  # Extract cis-eQTLs for this gene
  # ---------------------------------------------------------------------------

  gene_data <- eqtlgen[
    GeneSymbol == gene
  ]

  if (nrow(gene_data) == 0) {
    next
  }


  # ---------------------------------------------------------------------------
  # Construct ieugwasr clumping input
  # ---------------------------------------------------------------------------

  clump_input <- unique(
    data.frame(
      rsid = gene_data$SNP,
      pval = gene_data$Pvalue,
      id = gene,
      stringsAsFactors = FALSE
    )
  )

  # Remove incomplete records
  clump_input <- clump_input[
    !is.na(clump_input$rsid) &
      !is.na(clump_input$pval),
    ,
    drop = FALSE
  ]

  if (nrow(clump_input) == 0) {
    next
  }


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
    "  Eligible cis-eQTLs: ",
    nrow(clump_input)
  )

  message(
    "  LD-independent instruments retained: ",
    nrow(clumped)
  )

  if (nrow(clumped) == 0) {
    next
  }


  # ---------------------------------------------------------------------------
  # Recover complete eQTLGen records for retained instruments
  # ---------------------------------------------------------------------------

  gene_clumped <- gene_data[
    SNP %in% clumped$rsid
  ]

  gene_clumped <- gene_clumped[
    match(
      clumped$rsid,
      SNP
    )
  ]


  # ---------------------------------------------------------------------------
  # Save gene-specific instrument file
  # ---------------------------------------------------------------------------

  # Replace characters unsuitable for filenames
  gene_filename <- gsub(
    "[^A-Za-z0-9._-]",
    "_",
    gene
  )

  output_file <- file.path(
    output_dir,
    paste0(
      gene_filename,
      "_LDclumped.txt"
    )
  )

  fwrite(
    gene_clumped,
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
message("eQTLGen instrument preparation complete.")
message("P threshold: ", P_THRESHOLD)
message("LD r2 threshold: ", CLUMP_R2)
message("LD window: ", CLUMP_KB, " kb")
message("Output directory: ", output_dir)
