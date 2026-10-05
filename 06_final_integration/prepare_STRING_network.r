## ============================================================
## Prepare STRING protein interaction network
## ============================================================
##
## This script annotates STRING v12 protein-protein interactions
## with gene/protein names for use in the downstream target
## prioritisation pipeline.
##
## The STRING combined_score is retained without filtering.
## Confidence thresholds are applied downstream.
## ============================================================


library(data.table)


## ------------------------------------------------------------
## Input / output paths
## ------------------------------------------------------------

protein_info_file <- "/path/to/9606.protein.info.v12.0.txt"

protein_links_file <- "/path/to/9606.protein.links.v12.0.txt"

output_file <- "/path/to/output/STRING_v12_network_annotated.rds"


## ------------------------------------------------------------
## Load STRING v12 reference data
## ------------------------------------------------------------

protein_info <- fread(protein_info_file)

protein_links <- fread(protein_links_file)


## ------------------------------------------------------------
## Annotate STRING protein IDs with preferred names
## ------------------------------------------------------------

protein_links[
  ,
  protein1_name := protein_info[
    match(protein1, `#string_protein_id`),
    preferred_name
  ]
]

protein_links[
  ,
  protein2_name := protein_info[
    match(protein2, `#string_protein_id`),
    preferred_name
  ]
]


## ------------------------------------------------------------
## Basic checks
## ------------------------------------------------------------

message("STRING interactions: ", nrow(protein_links))

message(
  "Interactions with mapped protein names: ",
  sum(
    !is.na(protein_links$protein1_name) &
    !is.na(protein_links$protein2_name)
  )
)

message(
  "Combined score range: ",
  min(protein_links$combined_score, na.rm = TRUE),
  " - ",
  max(protein_links$combined_score, na.rm = TRUE)
)


## ------------------------------------------------------------
## Save annotated STRING network
## ------------------------------------------------------------

saveRDS(
  protein_links,
  file = output_file
)

message("Saved annotated STRING network to: ", output_file)