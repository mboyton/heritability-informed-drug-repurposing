# ==============================================================================
# SuSiE fine-mapping and cross-trait colocalisation
#
# For loci showing evidence of local genetic correlation in LAVA:
#   1. Extract GWAS summary statistics for both phenotypes
#   2. Harmonise alleles to the 1000 Genomes European reference
#   3. Generate an LD matrix using the European reference panel
#   4. Fine-map each phenotype using SuSiE
#   5. Perform cross-trait colocalisation using coloc.susie()
#
# Expected GWAS columns:
#   SNP, A1, A2, beta, pval
#
#
# Requirements:
#   - LAVA
#   - coloc
#   - ieugwasr
#   - genetics.binaRies
#   - stringr
#   - data.table
# ==============================================================================


library(LAVA)
library(coloc)
library(ieugwasr)
library(genetics.binaRies)
library(stringr)
library(data.table)


# ------------------------------------------------------------------------------
# User-defined parameters
# ------------------------------------------------------------------------------

phenotype_1 <- "phenotype_1"
phenotype_2 <- "phenotype_2"

analysis_id <- paste(
 phenotype_1,
 phenotype_2,
 sep = "_"
)

# LAVA bivariate threshold used to select loci for colocalisation
bivar_p_thresh <- 0.05

# GWAS summary statistics
phenotype_1_file <- "/path/to/phenotype_1_gwas.txt"
phenotype_2_file <- "/path/to/phenotype_2_gwas.txt"

# Bivariate LAVA results
bivar_results_file <- "/path/to/LAVA_bivariate_results.txt"

# LAVA dependencies
input_info_file <- "/path/to/input_info_file.txt"
sample_overlap_file <- "/path/to/sample_overlap_file.txt"
ref_prefix <- "/path/to/g1000_eur"
loci_file <- "/path/to/loci_2495.txt"

# 1000 Genomes European PLINK reference
eur_bfile <- "/path/to/EUR"

# Output directory
output_dir <- "/path/to/output"

dir.create(
 output_dir,
 recursive = TRUE,
 showWarnings = FALSE
)


# ==============================================================================
# Helper functions
# ==============================================================================


# ------------------------------------------------------------------------------
# Calculate SE from beta and a two-sided P-value
# ------------------------------------------------------------------------------

calculate_se <- function(beta, pval) {
 
 # QC FIX:
 # Explicitly return NA for P-values that cannot be used to reconstruct
 # a finite Wald Z statistic.
 valid <- (
  is.finite(beta) &
   is.finite(pval) &
   pval > 0 &
   pval < 1
 )
 
 se <- rep(
  NA_real_,
  length(beta)
 )
 
 z <- qnorm(
  1 - pval[valid] / 2
 )
 
 se[valid] <- abs(
  beta[valid] / z
 )
 
 se
}


# ------------------------------------------------------------------------------
# Identify strand-ambiguous palindromic SNPs
# ------------------------------------------------------------------------------


is_palindromic <- function(a1, a2) {
 
 allele_pair <- paste0(
  toupper(a1),
  toupper(a2)
 )
 
 allele_pair %in% c(
  "AT",
  "TA",
  "CG",
  "GC"
 )
}


# ------------------------------------------------------------------------------
# Harmonise GWAS alleles to reference alleles
# ------------------------------------------------------------------------------

# SNPs with directly matching alleles are retained.
# SNPs with reversed alleles are flipped and beta is multiplied by -1.
# Irreconcilable SNPs are removed.

harmonise_to_reference <- function(gwas, reference) {
 
 m <- match(
  gwas$SNP,
  reference$rsid
 )
 
 keep <- !is.na(m)
 
 gwas <- gwas[
  keep,
  ,
  drop = FALSE
 ]
 
 reference <- reference[
  m[keep],
  ,
  drop = FALSE
 ]
 
 
 # remove strand-ambiguous variants
 
 palindromic <- is_palindromic(
  reference$A1,
  reference$A2
 )
 
 n_palindromic <- sum(
  palindromic,
  na.rm = TRUE
 )
 
 if (any(palindromic, na.rm = TRUE)) {
  
  keep <- !palindromic
  
  gwas <- gwas[
   keep,
   ,
   drop = FALSE
  ]
  
  reference <- reference[
   keep,
   ,
   drop = FALSE
  ]
 }
 
 
 #  allele-matching logic
 
 is_match <- (
  gwas$A1 == reference$A1 &
   gwas$A2 == reference$A2
 )
 
 is_flip <- (
  gwas$A1 == reference$A2 &
   gwas$A2 == reference$A1
 )
 
 is_error <- !(
  is_match |
   is_flip
 )
 
 n_flip <- sum(
  is_flip,
  na.rm = TRUE
 )
 
 n_error <- sum(
  is_error,
  na.rm = TRUE
 )
 
 
 #  beta/allele flipping
 
 if (any(is_flip, na.rm = TRUE)) {
  
  old_a1 <- gwas$A1[
   is_flip
  ]
  
  gwas$A1[
   is_flip
  ] <- gwas$A2[
   is_flip
  ]
  
  gwas$A2[
   is_flip
  ] <- old_a1
  
  gwas$beta[
   is_flip
  ] <- -gwas$beta[
   is_flip
  ]
 }
 
 
 #  removal of incompatible alleles
 
 keep <- !is_error
 
 gwas <- gwas[
  keep,
  ,
  drop = FALSE
 ]
 
 reference <- reference[
  keep,
  ,
  drop = FALSE
 ]
 
 
 # Retain actual genomic position
 
 gwas$position <- reference$position
 
 
 list(
  data = gwas,
  n_flipped = n_flip,
  n_removed = n_error,
  n_palindromic_removed = n_palindromic
 )
}


# ==============================================================================
# 1. Load input data
# ==============================================================================

phenotype_1_raw <- fread(
 phenotype_1_file,
 data.table = FALSE
)

phenotype_2_raw <- fread(
 phenotype_2_file,
 data.table = FALSE
)


# Validate the minimum columns required by the original analysis.

required_gwas_columns <- c(
 "SNP",
 "A1",
 "A2",
 "beta",
 "pval"
)

missing_p1 <- setdiff(
 required_gwas_columns,
 names(phenotype_1_raw)
)

missing_p2 <- setdiff(
 required_gwas_columns,
 names(phenotype_2_raw)
)

if (length(missing_p1) > 0) {
 
 stop(
  "Phenotype 1 GWAS is missing required columns: ",
  paste(missing_p1, collapse = ", ")
 )
}

if (length(missing_p2) > 0) {
 
 stop(
  "Phenotype 2 GWAS is missing required columns: ",
  paste(missing_p2, collapse = ", ")
 )
}


# Standardise the fields used numerically later in the analysis.

phenotype_1_raw$beta <- as.numeric(
 phenotype_1_raw$beta
)

phenotype_1_raw$pval <- as.numeric(
 phenotype_1_raw$pval
)

phenotype_2_raw$beta <- as.numeric(
 phenotype_2_raw$beta
)

phenotype_2_raw$pval <- as.numeric(
 phenotype_2_raw$pval
)


# 1000 Genomes European reference BIM file

ref_snps_all <- fread(
 paste0(
  eur_bfile,
  ".bim"
 ),
 header = FALSE,
 data.table = FALSE
)

colnames(ref_snps_all) <- c(
 "chr",
 "rsid",
 "cm",
 "position",
 "A1",
 "A2"
)


# match() silently returns the first occurrence of a duplicated rsID.
# Exclude all duplicated rsIDs so that reference allele assignment is unique.

duplicate_ref_rsids <- unique(
 ref_snps_all$rsid[
  duplicated(ref_snps_all$rsid) |
   duplicated(
    ref_snps_all$rsid,
    fromLast = TRUE
   )
 ]
)

if (length(duplicate_ref_rsids) > 0) {
 
 warning(
  length(duplicate_ref_rsids),
  " duplicated rsIDs detected in the EUR reference and excluded."
 )
 
 ref_snps_all <- ref_snps_all[
  !(ref_snps_all$rsid %in% duplicate_ref_rsids),
  ,
  drop = FALSE
 ]
}


# ==============================================================================
# 2. Process LAVA inputs
# ==============================================================================

input <- process.input(
 input.info.file = input_info_file,
 sample.overlap.file = sample_overlap_file,
 ref.prefix = ref_prefix,
 phenos = c(
  phenotype_1,
  phenotype_2
 )
)

loci <- read.loci(
 loci_file
)

bivar_results <- fread(
 bivar_results_file,
 data.table = FALSE
)

# Retain loci showing nominal evidence of local genetic correlation

bivar_results <- bivar_results[
 !is.na(bivar_results$p) &
  bivar_results$p < bivar_p_thresh,
 ,
 drop = FALSE
]

loci_for_susie <- bivar_results$LOC

message(
 length(loci_for_susie),
 " LAVA loci with P < ",
 bivar_p_thresh,
 " selected for SuSiE colocalisation."
)


# ==============================================================================
# 3. Retrieve phenotype sample sizes
# ==============================================================================

input_info <- read.table(
 input_info_file,
 header = TRUE
)

phenotype_1_info <- input_info[
 input_info$phenotype == phenotype_1,
 ,
 drop = FALSE
]

phenotype_2_info <- input_info[
 input_info$phenotype == phenotype_2,
 ,
 drop = FALSE
]

if (
 nrow(phenotype_1_info) != 1 ||
 nrow(phenotype_2_info) != 1
) {
 
 stop(
  "Could not uniquely identify both phenotypes in input_info_file."
 )
}


phenotype_1_n <- (
 phenotype_1_info$cases +
  phenotype_1_info$controls
)

phenotype_2_n <- (
 phenotype_2_info$cases +
  phenotype_2_info$controls
)


# Explicit case fractions for the case-control coloc datasets.

phenotype_1_s <- (
 phenotype_1_info$cases /
  phenotype_1_n
)

phenotype_2_s <- (
 phenotype_2_info$cases /
  phenotype_2_n
)

if (
 !is.finite(phenotype_1_s) ||
 phenotype_1_s <= 0 ||
 phenotype_1_s >= 1
) {
 
 stop(
  "Invalid case fraction for phenotype 1."
 )
}

if (
 !is.finite(phenotype_2_s) ||
 phenotype_2_s <= 0 ||
 phenotype_2_s >= 1
) {
 
 stop(
  "Invalid case fraction for phenotype 2."
 )
}


# ==============================================================================
# 4. Initialise output
# ==============================================================================

timestamp <- format(
 Sys.time(),
 "%Y%m%d_%H%M%S"
)

output_log <- data.frame(
 
 locus_i = seq_along(
  loci_for_susie
 ),
 
 timestamp = NA_character_,
 
 locus = loci_for_susie,
 
 chr = NA_integer_,
 start = NA_integer_,
 stop = NA_integer_,
 
 locus_snps = NA_integer_,
 
 p1_snps = NA_integer_,
 p2_snps = NA_integer_,
 p1p2_snps = NA_integer_,
 
 ref_snps = NA_integer_,
 
 palindromic_snps_removed = NA_integer_,
 
 p1_ref_mismatch = NA_integer_,
 p1_ref_req_align = NA_integer_,
 
 p2_ref_mismatch = NA_integer_,
 p2_ref_req_align = NA_integer_,
 
 p1_beta_zero = NA_integer_,
 p2_beta_zero = NA_integer_,
 
 p1_betase_zero = NA_integer_,
 p2_betase_zero = NA_integer_,
 
 p1_snps_not_in_ldmatrix = NA_integer_,
 p2_snps_not_in_ldmatrix = NA_integer_,
 
 ldmatrix_nsnps = NA_integer_,
 
 p1_ld_alignment = NA_real_,
 p2_ld_alignment = NA_real_,
 
 status = "not_started",
 error_message = NA_character_,
 
 phenotype_1 = phenotype_1,
 phenotype_2 = phenotype_2,
 
 stringsAsFactors = FALSE
)

susie_results <- list()


# Output files

results_file <- file.path(
 output_dir,
 paste0(
  analysis_id,
  "_coloc_susie_",
  timestamp,
  ".rds"
 )
)

log_file <- file.path(
 output_dir,
 paste0(
  analysis_id,
  "_coloc_susie_log_",
  timestamp,
  ".txt"
 )
)


# ==============================================================================
# 5. Iterate over LAVA loci
# ==============================================================================

for (locus_i in seq_along(loci_for_susie)) {
 
 locus_id <- loci_for_susie[
  locus_i
 ]
 
 message("")
 message(
  "Processing locus ",
  locus_i,
  "/",
  length(loci_for_susie),
  " (LOC ",
  locus_id,
  ")"
 )
 
 
 # Catch unexpected locus-level failures so that one problematic locus does
 # not terminate analysis of all subsequent loci.
 
 tryCatch({
  
  output_log$status[
   locus_i
  ] <- "processing"
  
  
  # --------------------------------------------------------------------------
  # Identify LAVA locus
  # --------------------------------------------------------------------------
  
  locus_row <- loci[
   loci$LOC == locus_id,
   ,
   drop = FALSE
  ]
  
  if (nrow(locus_row) != 1) {
   
   stop(
    "Could not uniquely identify LOC ",
    locus_id,
    "."
   )
  }
  
  
  locus <- process.locus(
   locus_row,
   input
  )
  
  if (is.null(locus)) {
   
   stop(
    "LAVA could not process LOC ",
    locus_id,
    "."
   )
  }
  
  
  # Update locus-level log
  
  output_log$chr[
   locus_i
  ] <- locus_row$CHR
  
  output_log$start[
   locus_i
  ] <- locus_row$START
  
  output_log$stop[
   locus_i
  ] <- locus_row$STOP
  
  output_log$locus_snps[
   locus_i
  ] <- length(
   locus$snps
  )
  
  
  # --------------------------------------------------------------------------
  # Extract SNPs present in both GWAS datasets
  # --------------------------------------------------------------------------
  
  p1 <- phenotype_1_raw[
   phenotype_1_raw$SNP %in% locus$snps,
   ,
   drop = FALSE
  ]
  
  p2 <- phenotype_2_raw[
   phenotype_2_raw$SNP %in% locus$snps,
   ,
   drop = FALSE
  ]
  
  
  output_log$p1_snps[
   locus_i
  ] <- nrow(p1)
  
  output_log$p2_snps[
   locus_i
  ] <- nrow(p2)
  
  
  # Duplicate GWAS rsIDs make SNP matching/order ambiguous.
  # Remove all duplicated IDs within the locus.
  
  duplicate_p1 <- unique(
   p1$SNP[
    duplicated(p1$SNP) |
     duplicated(
      p1$SNP,
      fromLast = TRUE
     )
   ]
  )
  
  duplicate_p2 <- unique(
   p2$SNP[
    duplicated(p2$SNP) |
     duplicated(
      p2$SNP,
      fromLast = TRUE
     )
   ]
  )
  
  duplicate_gwas_snps <- unique(
   c(
    duplicate_p1,
    duplicate_p2
   )
  )
  
  if (length(duplicate_gwas_snps) > 0) {
   
   p1 <- p1[
    !(p1$SNP %in% duplicate_gwas_snps),
    ,
    drop = FALSE
   ]
   
   p2 <- p2[
    !(p2$SNP %in% duplicate_gwas_snps),
    ,
    drop = FALSE
   ]
  }
  
  
  common_snps <- intersect(
   p1$SNP,
   p2$SNP
  )
  
  output_log$p1p2_snps[
   locus_i
  ] <- length(
   common_snps
  )
  
  if (length(common_snps) == 0) {
   
   stop(
    "No overlapping GWAS SNPs at LOC ",
    locus_id,
    "."
   )
  }
  
  
  p1 <- p1[
   p1$SNP %in% common_snps,
   ,
   drop = FALSE
  ]
  
  p2 <- p2[
   p2$SNP %in% common_snps,
   ,
   drop = FALSE
  ]
  
  
  # --------------------------------------------------------------------------
  # Restrict to SNPs available in the 1000G European reference
  # --------------------------------------------------------------------------
  
  ref_snps <- ref_snps_all[
   ref_snps_all$rsid %in% common_snps,
   ,
   drop = FALSE
  ]
  
  output_log$ref_snps[
   locus_i
  ] <- nrow(
   ref_snps
  )
  
  if (nrow(ref_snps) == 0) {
   
   stop(
    "No SNPs available in EUR reference for LOC ",
    locus_id,
    "."
   )
  }
  
  
  p1 <- p1[
   p1$SNP %in% ref_snps$rsid,
   ,
   drop = FALSE
  ]
  
  p2 <- p2[
   p2$SNP %in% ref_snps$rsid,
   ,
   drop = FALSE
  ]
  
  
  # --------------------------------------------------------------------------
  # Harmonise both GWAS datasets to reference alleles
  # --------------------------------------------------------------------------
  
  p1_harmonised <- harmonise_to_reference(
   p1,
   ref_snps
  )
  
  p2_harmonised <- harmonise_to_reference(
   p2,
   ref_snps
  )
  
  
  output_log$p1_ref_mismatch[
   locus_i
  ] <- p1_harmonised$n_removed
  
  output_log$p1_ref_req_align[
   locus_i
  ] <- p1_harmonised$n_flipped
  
  output_log$p2_ref_mismatch[
   locus_i
  ] <- p2_harmonised$n_removed
  
  output_log$p2_ref_req_align[
   locus_i
  ] <- p2_harmonised$n_flipped
  
  
  # Palindromic SNPs are reference-defined and therefore identical between
  # the two harmonisation calls. Record the number once.
  
  output_log$palindromic_snps_removed[
   locus_i
  ] <- max(
   p1_harmonised$n_palindromic_removed,
   p2_harmonised$n_palindromic_removed
  )
  
  
  p1 <- p1_harmonised$data
  p2 <- p2_harmonised$data
  
  
  # --------------------------------------------------------------------------
  # Re-synchronise SNP sets after allele harmonisation
  # --------------------------------------------------------------------------
  
  common_snps <- intersect(
   p1$SNP,
   p2$SNP
  )
  
  if (length(common_snps) == 0) {
   
   stop(
    "No SNPs remain after allele harmonisation at LOC ",
    locus_id,
    "."
   )
  }
  
  
  p1 <- p1[
   match(
    common_snps,
    p1$SNP
   ),
   ,
   drop = FALSE
  ]
  
  p2 <- p2[
   match(
    common_snps,
    p2$SNP
   ),
   ,
   drop = FALSE
  ]
  
  
  # --------------------------------------------------------------------------
  # Remove SNPs with beta = 0
  # --------------------------------------------------------------------------
  
  p1_beta_zero <- p1$SNP[
   !is.na(p1$beta) &
    p1$beta == 0
  ]
  
  p2_beta_zero <- p2$SNP[
   !is.na(p2$beta) &
    p2$beta == 0
  ]
  
  
  output_log$p1_beta_zero[
   locus_i
  ] <- length(
   p1_beta_zero
  )
  
  output_log$p2_beta_zero[
   locus_i
  ] <- length(
   p2_beta_zero
  )
  
  
  zero_beta_snps <- unique(
   c(
    p1_beta_zero,
    p2_beta_zero
   )
  )
  
  if (length(zero_beta_snps) > 0) {
   
   p1 <- p1[
    !(p1$SNP %in% zero_beta_snps),
    ,
    drop = FALSE
   ]
   
   p2 <- p2[
    !(p2$SNP %in% zero_beta_snps),
    ,
    drop = FALSE
   ]
  }
  
  
  # --------------------------------------------------------------------------
  # Calculate standard errors
  # --------------------------------------------------------------------------
  
  # SE is reconstructed from beta and two-sided P-value.
  
  p1$beta_se <- calculate_se(
   p1$beta,
   p1$pval
  )
  
  p2$beta_se <- calculate_se(
   p2$beta,
   p2$pval
  )
  
  
  p1_bad <- p1$SNP[
   !is.finite(p1$beta_se) |
    p1$beta_se <= 0
  ]
  
  p2_bad <- p2$SNP[
   !is.finite(p2$beta_se) |
    p2$beta_se <= 0
  ]
  
  
  output_log$p1_betase_zero[
   locus_i
  ] <- length(
   p1_bad
  )
  
  output_log$p2_betase_zero[
   locus_i
  ] <- length(
   p2_bad
  )
  
  
  bad_snps <- unique(
   c(
    p1_bad,
    p2_bad
   )
  )
  
  if (length(bad_snps) > 0) {
   
   p1 <- p1[
    !(p1$SNP %in% bad_snps),
    ,
    drop = FALSE
   ]
   
   p2 <- p2[
    !(p2$SNP %in% bad_snps),
    ,
    drop = FALSE
   ]
  }
  
  
  # Re-synchronise SNP order
  
  common_snps <- intersect(
   p1$SNP,
   p2$SNP
  )
  
  p1 <- p1[
   match(
    common_snps,
    p1$SNP
   ),
   ,
   drop = FALSE
  ]
  
  p2 <- p2[
   match(
    common_snps,
    p2$SNP
   ),
   ,
   drop = FALSE
  ]
  
  
  if (length(common_snps) < 2) {
   
   stop(
    "Fewer than two SNPs remain at LOC ",
    locus_id,
    "."
   )
  }
  
  
  # --------------------------------------------------------------------------
  # Generate EUR LD matrix
  # --------------------------------------------------------------------------
  
  ld_matrix_raw <- ld_matrix(
   common_snps,
   plink_bin = genetics.binaRies::get_plink_binary(),
   bfile = eur_bfile
  )
  
  
  if (
   is.null(ld_matrix_raw) ||
   nrow(ld_matrix_raw) == 0
  ) {
   
   stop(
    "Could not construct LD matrix for LOC ",
    locus_id,
    "."
   )
  }
  
  
  # ieugwasr LD-matrix names contain SNP and allele information
  
  ld_names <- colnames(
   ld_matrix_raw
  )
  
  ld_snps <- str_extract(
   ld_names,
   "^[^_]+"
  )
  
  ld_a1 <- sub(
   "^.*_([^_]+)_([^_]+)$",
   "\\1",
   ld_names
  )
  
  ld_a2 <- sub(
   "^.*_([^_]+)_([^_]+)$",
   "\\2",
   ld_names
  )
  
  ld_info <- data.frame(
   SNP = ld_snps,
   A1 = ld_a1,
   A2 = ld_a2,
   stringsAsFactors = FALSE
  )
  
  
  # The LD matrix must contain one unique row/column per SNP.
  
  if (anyDuplicated(ld_info$SNP)) {
   
   stop(
    "Duplicate SNP identifiers returned in LD matrix at LOC ",
    locus_id,
    "."
   )
  }
  
  
  # --------------------------------------------------------------------------
  # Restrict GWAS datasets to SNPs returned in the LD matrix
  # --------------------------------------------------------------------------
  
  n_before_ld_p1 <- nrow(p1)
  n_before_ld_p2 <- nrow(p2)
  
  
  p1 <- p1[
   p1$SNP %in% ld_info$SNP,
   ,
   drop = FALSE
  ]
  
  p2 <- p2[
   p2$SNP %in% ld_info$SNP,
   ,
   drop = FALSE
  ]
  
  
  output_log$p1_snps_not_in_ldmatrix[
   locus_i
  ] <- (
   n_before_ld_p1 -
    nrow(p1)
  )
  
  output_log$p2_snps_not_in_ldmatrix[
   locus_i
  ] <- (
   n_before_ld_p2 -
    nrow(p2)
  )
  
  
  # --------------------------------------------------------------------------
  # Align GWAS SNP order to LD-matrix order
  # --------------------------------------------------------------------------
  
  ld_info <- ld_info[
   ld_info$SNP %in% p1$SNP &
    ld_info$SNP %in% p2$SNP,
   ,
   drop = FALSE
  ]
  
  
  p1 <- p1[
   match(
    ld_info$SNP,
    p1$SNP
   ),
   ,
   drop = FALSE
  ]
  
  p2 <- p2[
   match(
    ld_info$SNP,
    p2$SNP
   ),
   ,
   drop = FALSE
  ]
  
  
  ld_index <- match(
   ld_info$SNP,
   ld_snps
  )
  
  ld_matrix_final <- ld_matrix_raw[
   ld_index,
   ld_index,
   drop = FALSE
  ]
  
  
  # --------------------------------------------------------------------------
  # Harmonise LD matrix orientation to GWAS effect alleles
  # --------------------------------------------------------------------------
  
  is_ld_match <- (
   ld_info$A1 == p1$A1 &
    ld_info$A2 == p1$A2
  )
  
  is_ld_flip <- (
   ld_info$A1 == p1$A2 &
    ld_info$A2 == p1$A1
  )
  
  is_ld_error <- !(
   is_ld_match |
    is_ld_flip
  )
  
  
  if (any(is_ld_error)) {
   
   keep <- !is_ld_error
   
   p1 <- p1[
    keep,
    ,
    drop = FALSE
   ]
   
   p2 <- p2[
    keep,
    ,
    drop = FALSE
   ]
   
   ld_info <- ld_info[
    keep,
    ,
    drop = FALSE
   ]
   
   ld_matrix_final <- ld_matrix_final[
    keep,
    keep,
    drop = FALSE
   ]
   
   is_ld_flip <- is_ld_flip[
    keep
   ]
  }
  
  
  if (nrow(p1) < 2) {
   
   stop(
    "Fewer than two SNPs remain after LD-reference harmonisation at LOC ",
    locus_id,
    "."
   )
  }
  
  
  # Multiply rows/columns by -1 where LD-reference allele orientation
  # is reversed relative to the GWAS effect allele.
  
  multiplier <- ifelse(
   is_ld_flip,
   -1,
   1
  )
  
  ld_matrix_final <- (
   ld_matrix_final *
    outer(
     multiplier,
     multiplier
    )
  )
  
  
  rownames(
   ld_matrix_final
  ) <- p1$SNP
  
  colnames(
   ld_matrix_final
  ) <- p1$SNP
  
  
  output_log$ldmatrix_nsnps[
   locus_i
  ] <- nrow(
   ld_matrix_final
  )
  
  
  # --------------------------------------------------------------------------
  # Final alignment checks
  # --------------------------------------------------------------------------
  
  # In addition to SNP identity/order, explicitly require valid LD dimensions,
  # finite summary statistics, and finite LD values.
  
  stopifnot(
   identical(
    p1$SNP,
    p2$SNP
   ),
   identical(
    p1$SNP,
    rownames(ld_matrix_final)
   ),
   identical(
    p1$SNP,
    colnames(ld_matrix_final)
   ),
   all(
    p1$A1 == p2$A1
   ),
   all(
    p1$A2 == p2$A2
   ),
   nrow(ld_matrix_final) == nrow(p1),
   ncol(ld_matrix_final) == nrow(p1),
   all(is.finite(p1$beta)),
   all(is.finite(p2$beta)),
   all(is.finite(p1$beta_se)),
   all(is.finite(p2$beta_se)),
   all(is.finite(ld_matrix_final))
  )
  
  
  # --------------------------------------------------------------------------
  # Construct coloc/SuSiE datasets
  # --------------------------------------------------------------------------
  
  data_p1_susie <- list(
   
   beta = p1$beta,
   varbeta = p1$beta_se^2,
   snp = p1$SNP,
   type = "cc",
   LD = ld_matrix_final,
   N = phenotype_1_n,
   
   # Use actual reference genomic position rather than seq_len().
   position = p1$position,
   
   # Explicit case fraction for case-control data.
   s = phenotype_1_s
  )
  
  
  data_p2_susie <- list(
   
   beta = p2$beta,
   varbeta = p2$beta_se^2,
   snp = p2$SNP,
   type = "cc",
   LD = ld_matrix_final,
   N = phenotype_2_n,
   
   position = p2$position,
   
   s = phenotype_2_s
  )
  
  
  # --------------------------------------------------------------------------
  # Validate coloc/SuSiE datasets
  # --------------------------------------------------------------------------
  
  # Official coloc input validation prior to runsusie().
  
  coloc::check_dataset(
   data_p1_susie,
   req = "LD"
  )
  
  coloc::check_dataset(
   data_p2_susie,
   req = "LD"
  )
  
  
  # Assess consistency between the signs of GWAS effects and the signed
  # LD matrix. This is logged as a diagnostic rather than used as an
  # arbitrary exclusion threshold.
  
  alignment_p1 <- coloc::check_alignment(
   data_p1_susie,
   do_plot = FALSE
  )
  
  alignment_p2 <- coloc::check_alignment(
   data_p2_susie,
   do_plot = FALSE
  )
  
  
  output_log$p1_ld_alignment[
   locus_i
  ] <- alignment_p1
  
  output_log$p2_ld_alignment[
   locus_i
  ] <- alignment_p2
  
  
  # --------------------------------------------------------------------------
  # SuSiE fine-mapping and colocalisation
  # --------------------------------------------------------------------------
  
  
  susie_p1 <- coloc::runsusie(
   data_p1_susie
  )
  
  susie_p2 <- coloc::runsusie(
   data_p2_susie
  )
  
  coloc_result <- coloc::coloc.susie(
   susie_p1,
   susie_p2
  )
  
  
  susie_results[
   [as.character(locus_id)]
  ] <- coloc_result
  
  
  output_log$status[
   locus_i
  ] <- "success"
  
  
 }, error = function(e) {
  
  # Record the failure and continue to the next locus.
  
  output_log$status[
   locus_i
  ] <<- "failed"
  
  output_log$error_message[
   locus_i
  ] <<- conditionMessage(e)
  
  warning(
   "LOC ",
   locus_id,
   " failed: ",
   conditionMessage(e)
  )
 })
 
 
 # ---------------------------------------------------------------------------
 # Save progress after every locus
 # ---------------------------------------------------------------------------
 
 output_log$timestamp[
  locus_i
 ] <- format(
  Sys.time(),
  "%Y%m%d_%H%M%S"
 )
 
 
 saveRDS(
  susie_results,
  file = results_file
 )
 
 
 fwrite(
  output_log,
  file = log_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
 )
}


# ==============================================================================
# 6. Save reproducibility information
# ==============================================================================

# Record package/R versions used for the analysis.

session_info_file <- file.path(
 output_dir,
 paste0(
  analysis_id,
  "_sessionInfo_",
  timestamp,
  ".txt"
 )
)

writeLines(
 capture.output(
  sessionInfo()
 ),
 con = session_info_file
)


# ==============================================================================
# Complete
# ==============================================================================

message("")
message("SuSiE colocalisation workflow complete.")
message(
 "Phenotype pair: ",
 phenotype_1,
 " / ",
 phenotype_2
)
message(
 "LAVA loci considered: ",
 length(loci_for_susie)
)
message(
 "Successful loci: ",
 sum(
  output_log$status == "success",
  na.rm = TRUE
 )
)
message(
 "Failed loci: ",
 sum(
  output_log$status == "failed",
  na.rm = TRUE
 )
)
message(
 "Colocalisation results generated: ",
 length(susie_results)
)