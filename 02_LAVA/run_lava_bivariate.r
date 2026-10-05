# ==============================================================================
# LAVA bivariate local genetic correlation analysis
#
# Runs bivariate LAVA analysis between two phenotypes at loci showing
# univariate evidence (P < 0.05) for both phenotypes.
#
# Sample overlap is accounted for using cross-trait LDSC intercepts
# (gcov_int), with diagonal elements set to 1.
#
# Requirements:
#   - LAVA
#   - data.table
#   - LAVA univariate results
#   - matrix of LDSC cross-trait intercepts
#   - 1000 Genomes European reference data
#   - LAVA 2,495-locus definition file
# ==============================================================================

library(LAVA)
library(data.table)


# ------------------------------------------------------------------------------
# User-defined parameters
# ------------------------------------------------------------------------------

phenotype_1 <- "phenotype_1"
phenotype_2 <- "phenotype_2"

# Univariate LAVA significance threshold
univ_p_thresh <- 0.05

# Directory containing phenotype-specific univariate LAVA results
univariate_data_directory <- "/path/to/univariate/results"

# Matrix of LDSC cross-trait intercepts (gcov_int)
gcov_int_file <- "/path/to/gcov_int_matrix.txt"

# LAVA dependencies
input_info_file <- "/path/to/input_info_file.txt"
ref_prefix <- "/path/to/g1000_eur"
loci_file <- "/path/to/loci_2495.txt"

# Output directory
output_dir <- "/path/to/output"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


# ==============================================================================
# 1. Construct pair-specific sample-overlap matrix
# ==============================================================================

gcov_int_matrix <- as.matrix(
  read.table(
    gcov_int_file,
    header = TRUE,
    row.names = 1,
    check.names = FALSE
  )
)

phenotypes <- c(phenotype_1, phenotype_2)

if (!all(phenotypes %in% rownames(gcov_int_matrix)) ||
    !all(phenotypes %in% colnames(gcov_int_matrix))) {
  stop("One or both phenotypes are absent from the gcov_int matrix.")
}

gcov_int_subset <- gcov_int_matrix[
  phenotypes,
  phenotypes,
  drop = FALSE
]

# LAVA sample-overlap matrices require diagonal values of 1
diag(gcov_int_subset) <- 1

sample_overlap_file <- file.path(
  output_dir,
  paste0(
    phenotype_1,
    "_",
    phenotype_2,
    "_sample_overlap.txt"
  )
)

write.table(
  gcov_int_subset,
  file = sample_overlap_file,
  row.names = TRUE,
  col.names = NA,
  sep = "\t",
  quote = FALSE
)


# ==============================================================================
# 2. Load univariate LAVA results
# ==============================================================================

get_univariate_file <- function(phenotype) {

  phenotype_dir <- file.path(
    univariate_data_directory,
    phenotype
  )

  files <- list.files(
    phenotype_dir,
    pattern = "LAVA_univariate.*\\.txt$",
    full.names = TRUE
  )

  if (length(files) == 0) {
    stop("No merged LAVA univariate result found for: ", phenotype)
  }

  if (length(files) > 1) {
    files <- files[which.max(file.info(files)$mtime)]
    message("Multiple files found for ", phenotype, "; using most recent: ", files)
  }

  return(files)
}

pheno_1_univ <- fread(
  get_univariate_file(phenotype_1)
)

pheno_2_univ <- fread(
  get_univariate_file(phenotype_2)
)


# ==============================================================================
# 3. Process LAVA input
# ==============================================================================

input <- process.input(
  input.info.file = input_info_file,
  sample.overlap.file = sample_overlap_file,
  ref.prefix = ref_prefix,
  phenos = phenotypes
)

loci <- read.loci(loci_file)


# ==============================================================================
# 4. Identify loci with univariate evidence for both phenotypes
# ==============================================================================

pheno_1_sig_loci <- pheno_1_univ[
  !is.na(p) & p < univ_p_thresh,
  LOC
]

pheno_2_sig_loci <- pheno_2_univ[
  !is.na(p) & p < univ_p_thresh,
  LOC
]

common_loci <- intersect(
  pheno_1_sig_loci,
  pheno_2_sig_loci
)

message(
  "Loci with univariate P < ",
  univ_p_thresh,
  " for both phenotypes: ",
  length(common_loci),
  " of ",
  nrow(loci)
)


# ==============================================================================
# 5. Run bivariate LAVA
# ==============================================================================

bivar_table <- data.frame(
  LOC = loci$LOC,
  CHR = loci$CHR,
  START = loci$START,
  STOP = loci$STOP,
  phen1 = NA_character_,
  phen2 = NA_character_,
  rho = NA_real_,
  rho.lower = NA_real_,
  rho.upper = NA_real_,
  r2 = NA_real_,
  r2.lower = NA_real_,
  r2.upper = NA_real_,
  p = NA_real_
)

for (i in seq_along(common_loci)) {

  locus_id <- common_loci[i]

  message(
    "Processing locus ",
    i,
    "/",
    length(common_loci),
    " (LOC ",
    locus_id,
    ")"
  )

  # Match using the LAVA locus identifier rather than assuming
  # LOC corresponds to the dataframe row number.
  locus_row <- loci[loci$LOC == locus_id, ]

  if (nrow(locus_row) != 1) {
    warning("Could not uniquely identify LOC ", locus_id, ". Skipping.")
    next
  }

  locus <- process.locus(
    locus_row,
    input
  )

  if (is.null(locus)) {
    next
  }

  bivar_data <- run.bivar(
    locus,
    param.lim = 1000
  )

  if (is.null(bivar_data) || nrow(bivar_data) == 0) {
    next
  }

  result_row <- which(bivar_table$LOC == locus_id)

  bivar_table$phen1[result_row] <- bivar_data$phen1[1]
  bivar_table$phen2[result_row] <- bivar_data$phen2[1]
  bivar_table$rho[result_row] <- bivar_data$rho[1]
  bivar_table$rho.lower[result_row] <- bivar_data$rho.lower[1]
  bivar_table$rho.upper[result_row] <- bivar_data$rho.upper[1]
  bivar_table$r2[result_row] <- bivar_data$r2[1]
  bivar_table$r2.lower[result_row] <- bivar_data$r2.lower[1]
  bivar_table$r2.upper[result_row] <- bivar_data$r2.upper[1]
  bivar_table$p[result_row] <- bivar_data$p[1]
}


# ==============================================================================
# 6. Save results
# ==============================================================================

timestamp <- format(
  Sys.time(),
  "%Y%m%d_%H%M%S"
)

output_file <- file.path(
  output_dir,
  paste0(
    phenotype_1,
    "_",
    phenotype_2,
    "_LAVA_bivariate_",
    timestamp,
    ".txt"
  )
)

fwrite(
  bivar_table,
  file = output_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# ==============================================================================
# 7. Summary
# ==============================================================================

n_tested <- sum(!is.na(bivar_table$p))

n_nominal <- sum(
  bivar_table$p < 0.05,
  na.rm = TRUE
)

message("")
message("Bivariate LAVA analysis complete.")
message("Phenotype pair: ", phenotype_1, " / ", phenotype_2)
message("Loci passing univariate filter: ", length(common_loci))
message("Loci successfully tested: ", n_tested)
message("Loci with bivariate P < 0.05: ", n_nominal)
message("Results: ", output_file)