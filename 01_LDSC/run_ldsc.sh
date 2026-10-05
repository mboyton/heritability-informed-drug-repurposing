#!/bin/bash
set -euo pipefail

# ==============================================================================
# LD Score Regression (LDSC)
#
# Performs:
#   1. GWAS summary-statistic munging
#   2. Cross-trait genetic correlation analysis
#   3. Collation of genetic correlation results
#
# Expected GWAS columns:
#   SNP, A1, A2, beta, pval, gwas_ntotal
#
# Requirements:
#   - LDSC
#   - LDSC conda environment
#   - HapMap3 SNP list (w_hm3.snplist)
#   - European LD-score reference files (eur_w_ld_chr/)
# ==============================================================================


# ------------------------------------------------------------------------------
# User-defined paths
# ------------------------------------------------------------------------------

# Directory containing input GWAS summary-statistic .txt files
GWAS_DIR="/path/to/gwas"

# Cloned LDSC repository
LDSC_DIR="/path/to/ldsc"

# HapMap3 SNP list
HAPMAP3_SNPLIST="${LDSC_DIR}/w_hm3.snplist"

# European LD-score reference panel
LD_REFERENCE="${LDSC_DIR}/eur_w_ld_chr/"


# ------------------------------------------------------------------------------
# Initialise LDSC environment
# ------------------------------------------------------------------------------

cd "$LDSC_DIR"

# Initialise conda for non-interactive shell use
source "$(conda info --base)/etc/profile.d/conda.sh"

if ! conda env list | awk '{print $1}' | grep -qx "ldsc"; then
    echo "LDSC conda environment not found."
    echo "Creating environment from environment.yml..."
    conda env create --file environment.yml
fi

conda activate ldsc


# ------------------------------------------------------------------------------
# Output directories
# ------------------------------------------------------------------------------

TIMESTAMP=$(date +"%Y%m%d_%H%M%S")

MUNGED_DIR="${GWAS_DIR}/ldsc_munged_${TIMESTAMP}"
RESULTS_DIR="${GWAS_DIR}/ldsc_results_${TIMESTAMP}"

mkdir -p "$MUNGED_DIR" "$RESULTS_DIR"


# ==============================================================================
# 1. Munge GWAS summary statistics
# ==============================================================================

echo "Munging GWAS summary statistics..."

for GWAS_FILE in "$GWAS_DIR"/*.txt; do

    BASENAME=$(basename "$GWAS_FILE" .txt)

    # Extract total GWAS sample size from the first data row
    GWAS_NTOTAL=$(awk '
        NR == 1 {
            for (i = 1; i <= NF; i++) {
                if ($i == "gwas_ntotal") col = i
            }
        }
        NR == 2 {
            if (col) print $col
            exit
        }
    ' "$GWAS_FILE")

    # Skip datasets without a valid integer sample size
    if ! [[ "$GWAS_NTOTAL" =~ ^[0-9]+$ ]]; then
        echo "WARNING: Invalid or missing gwas_ntotal for ${BASENAME}. Skipping."
        continue
    fi

    echo "Processing ${BASENAME} (N=${GWAS_NTOTAL})..."

    ./munge_sumstats.py \
        --sumstats "$GWAS_FILE" \
        --snp SNP \
        --a1 A1 \
        --a2 A2 \
        --signed-sumstats beta,0 \
        --p pval \
        --N "$GWAS_NTOTAL" \
        --chunksize 500000 \
        --merge-alleles "$HAPMAP3_SNPLIST" \
        --out "${MUNGED_DIR}/${BASENAME}_munged"

done

echo "GWAS munging complete."


# ==============================================================================
# 2. Cross-trait genetic correlation
# ==============================================================================

echo "Running cross-trait genetic correlations..."

MUNGED_FILES=("$MUNGED_DIR"/*_munged.sumstats.gz)

if [ ${#MUNGED_FILES[@]} -lt 2 ]; then
    echo "ERROR: At least two successfully munged GWAS datasets are required."
    exit 1
fi

# Construct comma-separated list of all munged datasets
ALL_FILES=$(IFS=,; echo "${MUNGED_FILES[*]}")

for FILE_1 in "${MUNGED_FILES[@]}"; do

    TRAIT=$(basename "$FILE_1" _munged.sumstats.gz)

    echo "Running LDSC genetic correlations for ${TRAIT}..."

    ./ldsc.py \
        --rg "${FILE_1},${ALL_FILES}" \
        --ref-ld-chr "$LD_REFERENCE" \
        --w-ld-chr "$LD_REFERENCE" \
        --out "${RESULTS_DIR}/${TRAIT}_vs_all"

done

echo "Genetic correlation analyses complete."


# ==============================================================================
# 3. Collate LDSC genetic correlation results
# ==============================================================================

SUMMARY_FILE="${RESULTS_DIR}/ldsc_genetic_correlations_${TIMESTAMP}.txt"

echo -e \
"p1\tp1_userfacinglabel\tp2\tp2_userfacinglabel\trg\tse\tz\tp\th2_obs\th2_obs_se\th2_int\th2_int_se\tgcov_int\tgcov_int_se" \
> "$SUMMARY_FILE"

for LOG_FILE in "$RESULTS_DIR"/*_vs_all.log; do

    awk '/^Summary of Genetic Correlation Results/,/^Analysis finished/' "$LOG_FILE" |
    awk 'NR > 2 && !/^Analysis finished/ && NF > 0 { print $0 }' |
    while read -r LINE; do

        P1=$(echo "$LINE" | awk '{print $1}')
        P2=$(echo "$LINE" | awk '{print $2}')

        # Generate concise phenotype labels from GWAS filenames
        P1_LABEL=$(echo "$P1" | sed -E 's|.*/([^/]+)_build.*|\1|')
        P2_LABEL=$(echo "$P2" | sed -E 's|.*/([^/]+)_build.*|\1|')

        RESULTS=$(echo "$LINE" |
            awk '{$1=$2=""; print $0}' |
            sed 's/^ *//')

        echo -e \
            "${P1}\t${P1_LABEL}\t${P2}\t${P2_LABEL}\t${RESULTS}" \
            >> "$SUMMARY_FILE"

    done

done

# Remove blank lines
sed -i '/^$/d' "$SUMMARY_FILE"

echo
echo "LDSC workflow complete."
echo "Munged GWAS: ${MUNGED_DIR}"
echo "LDSC results: ${RESULTS_DIR}"
echo "Combined results: ${SUMMARY_FILE}"