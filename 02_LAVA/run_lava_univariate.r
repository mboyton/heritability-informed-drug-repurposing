# ==============================================================================
# LAVA univariate local heritability analysis
#
# Runs univariate LAVA analysis across 2,495 approximately LD-independent
# genomic loci for a single phenotype.
#
# Requirements:
#   - LAVA
#   - data.table
#   - LAVA input-info file
#   - 1000 Genomes European reference data
#   - LAVA 2,495-locus definition file
# ==============================================================================

library(LAVA)
library(data.table)


# ------------------------------------------------------------------------------
# User-defined parameters
# ------------------------------------------------------------------------------

# Phenotype name as specified in the LAVA input-info file
phenotype <- "phenotype_name"

# LAVA input files
input_info_file <- "/path/to/input_info_file.txt"
sample_overlap_file <- NULL

# 1000 Genomes European reference prefix
ref_prefix <- "/path/to/g1000_eur"

# LAVA 2,495-locus definition file
loci_file <- "/path/to/loci_2495.txt"

# Output directory
output_dir <- "/path/to/output"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


# ==============================================================================
# 1. Process LAVA input
# ==============================================================================

input <- process.input(
  input.info.file = input_info_file,
  sample.overlap.file = sample_overlap_file,
  ref.prefix = ref_prefix,
  phenos = phenotype
)

loci <- read.loci(loci_file)

chromosomes <- sort(unique(loci$CHR))


# ==============================================================================
# 2. Run univariate analysis
# ==============================================================================

chromosome_results <- vector(
  mode = "list",
  length = length(chromosomes)
)

for (chr_index in seq_along(chromosomes)) {

  chr <- chromosomes[chr_index]

  message(
    "Analysing chromosome ",
    chr,
    " (",
    chr_index,
    "/",
    length(chromosomes),
    ")"
  )

  loci_chr <- loci[loci$CHR == chr, ]

  # Initialise chromosome-level results table
  univ_table <- data.frame(
    LOC = loci_chr$LOC,
    CHR = loci_chr$CHR,
    START = loci_chr$START,
    STOP = loci_chr$STOP,
    phen = NA_character_,
    h2.obs = NA_real_,
    p = NA_real_,
    n_snps = NA_integer_
  )

  for (i in seq_len(nrow(loci_chr))) {

    message(
      "  Processing locus ",
      i,
      "/",
      nrow(loci_chr),
      " (LOC ",
      loci_chr$LOC[i],
      ")"
    )

    locus <- process.locus(
      loci_chr[i, ],
      input
    )

    # process.locus() may return NULL where a locus cannot be analysed
    if (is.null(locus)) {
      next
    }

    univ_data <- run.univ(locus)

    if (nrow(univ_data) == 0) {
      next
    }

    # One phenotype is analysed per run
    if (nrow(univ_data) == 1) {

      univ_table$phen[i] <- univ_data$phen
      univ_table$h2.obs[i] <- univ_data$h2.obs
      univ_table$p[i] <- univ_data$p
      univ_table$n_snps[i] <- locus$n.snps

    }
  }

  chromosome_results[[chr_index]] <- univ_table

  # Save chromosome-level results
  chromosome_output <- file.path(
    output_dir,
    paste0(
      phenotype,
      "_LAVA_univ_chr_",
      chr,
      ".txt"
    )
  )

  write.table(
    univ_table,
    file = chromosome_output,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}


# ==============================================================================
# 3. Merge chromosome-level results
# ==============================================================================

univ_merged <- rbindlist(
  chromosome_results,
  use.names = TRUE,
  fill = TRUE
)

timestamp <- format(
  Sys.time(),
  "%Y%m%d_%H%M%S"
)

merged_output <- file.path(
  output_dir,
  paste0(
    phenotype,
    "_LAVA_univariate_",
    timestamp,
    ".txt"
  )
)

fwrite(
  univ_merged,
  file = merged_output,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# ==============================================================================
# Complete
# ==============================================================================

message("")
message("LAVA univariate analysis complete.")
message("Phenotype: ", phenotype)
message("Loci analysed: ", nrow(loci))
message("Results: ", merged_output)