## ============================================================
## Prepare OmniPath and Open Targets annotation data
## ============================================================
##
## Outputs:
##   1. OmniPath signed/directed interaction network
##   2. Open Targets pharmacological annotation
##
## The curated Open Targets disease-to-phenotype mapping is
## supplied separately as:
##   annotated_opentargets_diseases_v0-1.csv
## ============================================================


library(otargen)
library(biomaRt)
library(data.table)
library(dplyr)
library(OmnipathR)


## ------------------------------------------------------------
## Input / output paths
## ------------------------------------------------------------

# Molecular exposures used to define the protein universe
exposure_file <- "/path/to/LDclmpd_cispQTL_UKBBPPP.txt"

# Manually curated mapping of Open Targets disease names
# to phenotypes included in the analysis
phenotype_mapping_file <-
  "/path/to/annotated_opentargets_diseases_v0-1.csv"

# Outputs
omnipath_output_file <-
  "/path/to/output/OmniPath_interactions.rds"

opentargets_output_file <-
  "/path/to/output/OpenTargets_annotation.rds"


## ============================================================
## 1. Import OmniPath interaction network
## ============================================================

omnipath_interactions <- import_omnipath_interactions()

message(
  "OmniPath interactions imported: ",
  nrow(omnipath_interactions)
)

# These fields are subsequently used for signed and directed
# network-mediated target prioritisation.
required_omnipath_columns <- c(
  "source_genesymbol",
  "target_genesymbol",
  "is_stimulation",
  "is_inhibition"
)

missing_columns <- setdiff(
  required_omnipath_columns,
  colnames(omnipath_interactions)
)

if (length(missing_columns) > 0) {
  stop(
    "Required OmniPath columns not found: ",
    paste(missing_columns, collapse = ", ")
  )
}

saveRDS(
  omnipath_interactions,
  file = omnipath_output_file
)


## ============================================================
## 2. Define protein universe for Open Targets annotation
## ============================================================

# Load molecular exposures
exposures <- fread(exposure_file)

if (!"id.exposure" %in% colnames(exposures)) {
  stop("Expected column 'id.exposure' not found in exposure file.")
}

protein_universe <- unique(
  c(
    omnipath_interactions$source_genesymbol,
    omnipath_interactions$target_genesymbol,
    exposures$id.exposure
  )
)

protein_universe <- protein_universe[
  !is.na(protein_universe) &
  protein_universe != ""
]

message(
  "Proteins in Open Targets query universe: ",
  length(protein_universe)
)


## ============================================================
## 3. Map HGNC symbols to Ensembl gene IDs
## ============================================================

mart <- useMart(
  "ensembl",
  dataset = "hsapiens_gene_ensembl"
)

mapping <- getBM(
  attributes = c(
    "hgnc_symbol",
    "ensembl_gene_id"
  ),
  filters = "hgnc_symbol",
  values = protein_universe,
  mart = mart
)

message(
  "HGNC symbols mapped to Ensembl IDs: ",
  length(unique(mapping$hgnc_symbol))
)


## ============================================================
## 4. Query Open Targets
## ============================================================

open_targets_list <- vector(
  "list",
  nrow(mapping)
)

for (i in seq_len(nrow(mapping))) {

  message(
    "Querying protein ",
    i,
    " of ",
    nrow(mapping),
    ": ",
    mapping$hgnc_symbol[i]
  )

  query_res <- knownDrugsQuery(
    mapping$ensembl_gene_id[i],
    cursor = NULL,
    freeTextQuery = NULL,
    size = 10
  )

  if (!is.null(query_res) && nrow(query_res) > 0) {

    colnames(query_res) <- paste(
      "OPENTARGETS",
      colnames(query_res),
      sep = "_"
    )

    query_res$query_protein <-
      mapping$hgnc_symbol[i]

    query_res$query_protein_ensid <-
      mapping$ensembl_gene_id[i]

    open_targets_list[[i]] <- query_res
  }
}

open_targets_list <- Filter(
  Negate(is.null),
  open_targets_list
)

if (length(open_targets_list) == 0) {
  stop("No Open Targets drug annotations were returned.")
}

open_targets_annotation <- bind_rows(
  open_targets_list
)


## ============================================================
## 5. Add curated phenotype mapping
## ============================================================

phenotype_mapping <- read.csv(
  phenotype_mapping_file,
  stringsAsFactors = FALSE
)

if (!all(
  c(
    "OPENTARGETS_disease.name",
    "relevant_phenotype"
  ) %in% colnames(phenotype_mapping)
)) {
  stop(
    "Required columns are missing from phenotype mapping file."
  )
}

# Retain disease terms manually identified as relevant to
# phenotypes included in the analysis.
phenotype_mapping <- phenotype_mapping[
  phenotype_mapping$relevant_phenotype == TRUE,
]

open_targets_annotation <- open_targets_annotation %>%
  left_join(
    phenotype_mapping,
    by = "OPENTARGETS_disease.name"
  )


## ============================================================
## 6. Extract Open Targets mechanism of action
## ============================================================

n <- nrow(open_targets_annotation)

MOA_actionType <- character(n)
MOA_targets <- character(n)

for (i in seq_len(n)) {

  entry <-
    open_targets_annotation$
      OPENTARGETS_drug.mechanismsOfAction.rows[[i]]

  if (
    !is.null(entry) &&
    is.data.frame(entry) &&
    nrow(entry) > 0
  ) {

    MOA_actionType[i] <-
      as.character(entry$actionType[1])

    MOA_targets[i] <-
      as.character(entry$targets[1])

  } else {

    MOA_actionType[i] <- NA
    MOA_targets[i] <- NA
  }
}

open_targets_annotation$MOA_actionType <-
  MOA_actionType

open_targets_annotation$MOA_targets <-
  MOA_targets


## ============================================================
## 7. Convert mechanism of action to MR-interpretable classes
## ============================================================

moa_map <- c(

  # Inhibitory pharmacology
  "INHIBITOR" = "inhibits",
  "ANTAGONIST" = "inhibits",
  "ANTISENSE INHIBITOR" = "inhibits",
  "RNAI INHIBITOR" = "inhibits",
  "BLOCKER" = "inhibits",
  "NEGATIVE MODULATOR" = "inhibits",
  "NEGATIVE ALLOSTERIC MODULATOR" = "inhibits",
  "INVERSE AGONIST" = "inhibits",

  # Activating pharmacology
  "AGONIST" = "activates",
  "PARTIAL AGONIST" = "activates",
  "ACTIVATOR" = "activates",
  "OPENER" = "activates",
  "RELEASING AGENT" = "activates",
  "POSITIVE MODULATOR" = "activates",
  "POSITIVE ALLOSTERIC MODULATOR" = "activates",
  "STABILISER" = "activates",

  # Mimetic pharmacology
  "EXOGENOUS PROTEIN" = "mimics",
  "EXOGENOUS GENE" = "mimics",
  "SUBSTRATE" = "mimics",

  # Ambiguous pharmacology
  "PROTEOLYTIC ENZYME" = "ambiguous",
  "HYDROLYTIC ENZYME" = "ambiguous",
  "DISRUPTING AGENT" = "ambiguous",
  "CROSS-LINKING AGENT" = "ambiguous",
  "BINDING AGENT" = "ambiguous",
  "MODULATOR" = "ambiguous",
  "OTHER" = "ambiguous",

  # Vaccine
  "VACCINE ANTIGEN" = "vaccine"
)

open_targets_annotation$effect_on_target <-
  unname(
    moa_map[
      open_targets_annotation$MOA_actionType
    ]
  )


## ============================================================
## 8. Summary
## ============================================================

message(
  "Open Targets annotation rows: ",
  nrow(open_targets_annotation)
)

message(
  "Unique queried proteins with annotations: ",
  length(
    unique(open_targets_annotation$query_protein)
  )
)

message(
  "MOA classifications:"
)

print(
  table(
    open_targets_annotation$effect_on_target,
    useNA = "ifany"
  )
)


## ============================================================
## 9. Save Open Targets annotation
## ============================================================

saveRDS(
  open_targets_annotation,
  file = opentargets_output_file
)

message(
  "Saved OmniPath interactions to: ",
  omnipath_output_file
)

message(
  "Saved Open Targets annotation to: ",
  opentargets_output_file
)