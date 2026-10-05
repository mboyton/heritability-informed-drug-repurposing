## ============================================================
## Cross-trait target prioritisation and drug repurposing
## ============================================================
##
##
## Integrates:
##   - LAVA local genetic correlation
##   - SuSiE colocalisation
##   - molecular QTL SNP-to-exposure mapping
##   - STRING protein interaction data
##   - two-sample MR
##   - Open Targets pharmacological annotation
##   - OmniPath signed/directed interactions
##
## Two analysis modes were used:
##
## High-confidence:
##   LAVA P <= 0.01
##   coloc PP.H4 >= 0.90
##   STRING combined score >= 900
##   MR concordance score >= 0.90
##   MR null-effect score >= 0.90
##   OmniPath network-mediated targets excluded
##
## Exploratory:
##   LAVA P <= 0.05
##   coloc PP.H4 >= 0.70
##   STRING combined score >= 800
##   MR concordance score >= 0.85
##   MR null-effect score >= 0.85
##   OmniPath network-mediated targets included
## ============================================================


library(data.table)
library(dplyr)


## ============================================================
## 1. Analysis parameters
## ============================================================

# Currently configured for high-confidence analysis.

LAVA_P_THRESHOLD <- 0.01          # 0.05 for exploratory
COLOC_PPH4_THRESHOLD <- 0.90      # 0.70 for exploratory
STRING_SCORE_THRESHOLD <- 900     # 800 for exploratory
MR_CONCORDANCE_THRESHOLD <- 0.90  # 0.85 for exploratory
MR_NULL_THRESHOLD <- 0.90         # 0.85 for exploratory

INCLUDE_STRING <- TRUE
INCLUDE_OMNIPATH <- FALSE         # TRUE for exploratory


## Molecular QTL source:
##   "pQTL" = UKBB-PPP
##   "eQTL" = eQTLGen

MOLECULAR_QTL <- "pQTL"


## ============================================================
## 2. Input files
## ============================================================

lava_file <-
  "/path/to/LAVA_bivariate_results.rds"

coloc_file <-
  "/path/to/merged_coloc_susie_results.rds"

string_file <-
  "/path/to/STRING_v12_network_annotated.rds"

omnipath_file <-
  "/path/to/OmniPath_interactions.rds"

opentargets_file <-
  "/path/to/OpenTargets_annotation.rds"


# Molecular-QTL-specific inputs

ukbb_snp_tagger_file <-
  "/path/to/UKBB_PPP_coloc_exposures.rds"

ukbb_mr_file <-
  "/path/to/UKBB_PPP_MR_results.rds"


eqtlgen_snp_tagger_file <-
  "/path/to/eQTLGen_coloc_exposures.rds"

eqtlgen_mr_file <-
  "/path/to/eQTLGen_MR_results.rds"


## Output

output_file <-
  "/path/to/output/target_prioritisation_results.rds"


## ============================================================
## 3. Load data
## ============================================================

lava_bivar <- as.data.table(
  readRDS(lava_file)
)

coloc_results <- as.data.table(
  readRDS(coloc_file)
)

string_network <- as.data.table(
  readRDS(string_file)
)

omnipath_db <- as.data.table(
  readRDS(omnipath_file)
)

opentargets_db <- as.data.table(
  readRDS(opentargets_file)
)


if (MOLECULAR_QTL == "pQTL") {

  snp_tagger <- as.data.table(
    readRDS(ukbb_snp_tagger_file)
  )

  mr_results <- as.data.table(
    readRDS(ukbb_mr_file)
  )

} else if (MOLECULAR_QTL == "eQTL") {

  snp_tagger <- as.data.table(
    readRDS(eqtlgen_snp_tagger_file)
  )

  mr_results <- as.data.table(
    readRDS(eqtlgen_mr_file)
  )

} else {

  stop(
    "MOLECULAR_QTL must be either 'pQTL' or 'eQTL'."
  )
}


## ============================================================
## 4. MR scoring functions
## ============================================================

compute_concordance_score <- function(
  beta1,
  se1,
  beta2,
  se2
) {

  p1_positive <- pnorm(
    0,
    mean = beta1,
    sd = se1,
    lower.tail = FALSE
  )

  p2_positive <- pnorm(
    0,
    mean = beta2,
    sd = se2,
    lower.tail = FALSE
  )

  p1_negative <- 1 - p1_positive
  p2_negative <- 1 - p2_positive

  # Probability that both effects have the same direction
  p1_positive * p2_positive +
    p1_negative * p2_negative
}


compute_null_effect_score <- function(
  beta1,
  se1,
  beta2,
  se2
) {

  # Two-sided p-values against the null
  pval1 <- 2 * (
    1 - pnorm(abs(beta1 / se1))
  )

  pval2 <- 2 * (
    1 - pnorm(abs(beta2 / se2))
  )

  nonnull1 <- 1 - pval1
  nonnull2 <- 1 - pval2

  # Median strength of evidence against the null
  median(
    c(nonnull1, nonnull2)
  )
}


## ============================================================
## 5. Identify phenotype combinations
## ============================================================

phenotype_pairs <- unique(
  coloc_results[
    ,
    .(
      phenotype_1 = phenotype_1,
      phenotype_2 = phenotype_2
    )
  ]
)


## ============================================================
## 6. Initialise outputs
## ============================================================

all_mr_results <- list()

all_direct_repurposing_hits <- list()

all_direct_multimodal_targets <- list()

all_network_repurposing_hits <- list()

all_network_multimodal_targets <- list()


## ============================================================
## 7. Analyse each phenotype pair
## ============================================================

for (pair_index in seq_len(nrow(phenotype_pairs))) {

  phenotype_1 <-
    phenotype_pairs$phenotype_1[pair_index]

  phenotype_2 <-
    phenotype_pairs$phenotype_2[pair_index]

  pair_id <- paste(
    phenotype_1,
    phenotype_2,
    sep = "__"
  )


  message(
    "Processing: ",
    phenotype_1,
    " vs ",
    phenotype_2
  )


  ## ----------------------------------------------------------
  ## 7.1 LAVA local genetic correlation
  ## ----------------------------------------------------------

  lava_pair <- lava_bivar[
    (phen1 == phenotype_1 &
       phen2 == phenotype_2) |
      (phen1 == phenotype_2 &
         phen2 == phenotype_1)
  ]


  if (nrow(lava_pair) == 0) {

    message(
      "No LAVA results for ",
      pair_id
    )

    next
  }


  significant_lava_loci <- lava_pair[
    p < LAVA_P_THRESHOLD,
    unique(LOC)
  ]


  if (length(significant_lava_loci) == 0) {
    next
  }


  ## ----------------------------------------------------------
  ## 7.2 Colocalisation
  ## ----------------------------------------------------------

  coloc_pair <- coloc_results[
    (
      phenotype_1 == ..phenotype_1 &
      phenotype_2 == ..phenotype_2
    ) |
      (
        phenotype_1 == ..phenotype_2 &
        phenotype_2 == ..phenotype_1
      )
  ]


  coloc_pair <- coloc_pair[
    PP.H4.abf >= COLOC_PPH4_THRESHOLD &
      LOC %in% significant_lava_loci
  ]


  if (nrow(coloc_pair) == 0) {
    next
  }


  coloc_snps <- unique(
    c(
      coloc_pair$hit1,
      coloc_pair$hit2
    )
  )

  coloc_snps <- coloc_snps[
    !is.na(coloc_snps) &
      coloc_snps != ""
  ]


  ## ----------------------------------------------------------
  ## 7.3 Map colocalisation SNPs to molecular exposures
  ## ----------------------------------------------------------

  tagged_exposures <- unique(
    snp_tagger[
      SNP %in% coloc_snps,
      exposure
    ]
  )


  tagged_exposures <- tagged_exposures[
    !is.na(tagged_exposures) &
      tagged_exposures != ""
  ]


  if (length(tagged_exposures) == 0) {
    next
  }


  ## ----------------------------------------------------------
  ## 7.4 STRING network expansion
  ## ----------------------------------------------------------

  candidate_exposures <- tagged_exposures


  if (INCLUDE_STRING) {

    string_subset <- string_network[
      combined_score >= STRING_SCORE_THRESHOLD &
        protein1_name %in% tagged_exposures
    ]


    string_interactors <- unique(
      string_subset$protein2_name
    )


    string_interactors <- string_interactors[
      !is.na(string_interactors) &
        string_interactors != ""
    ]


    candidate_exposures <- unique(
      c(
        candidate_exposures,
        string_interactors
      )
    )
  }


  ## ----------------------------------------------------------
  ## 7.5 Retrieve MR results for both phenotypes
  ## ----------------------------------------------------------

  mr_pheno_1 <- mr_results[
    outcome == phenotype_1 &
      exposure %in% candidate_exposures &
      !is.na(mr_ivw_beta)
  ]


  mr_pheno_2 <- mr_results[
    outcome == phenotype_2 &
      exposure %in% candidate_exposures &
      !is.na(mr_ivw_beta)
  ]


  if (
    nrow(mr_pheno_1) == 0 ||
    nrow(mr_pheno_2) == 0
  ) {
    next
  }


  mr_pair <- merge(
    mr_pheno_1,
    mr_pheno_2,
    by = "exposure",
    suffixes = c(
      "_pheno_1",
      "_pheno_2"
    )
  )


  if (nrow(mr_pair) == 0) {
    next
  }


  ## ----------------------------------------------------------
  ## 7.6 Score MR concordance
  ## ----------------------------------------------------------

  mr_pair <- mr_pair %>%
    rowwise() %>%
    mutate(

      concordance_score =
        compute_concordance_score(
          mr_ivw_beta_pheno_1,
          mr_ivw_se_pheno_1,
          mr_ivw_beta_pheno_2,
          mr_ivw_se_pheno_2
        ),

      null_effect_score =
        compute_null_effect_score(
          mr_ivw_beta_pheno_1,
          mr_ivw_se_pheno_1,
          mr_ivw_beta_pheno_2,
          mr_ivw_se_pheno_2
        )

    ) %>%
    ungroup()


  mr_pair$phenotype_1 <- phenotype_1
  mr_pair$phenotype_2 <- phenotype_2


  all_mr_results[[pair_id]] <-
    mr_pair


  ## ----------------------------------------------------------
  ## 7.7 Retain concordant MR effects
  ## ----------------------------------------------------------

  concordant_hits <- mr_pair %>%
    filter(
      concordance_score >=
        MR_CONCORDANCE_THRESHOLD,
      null_effect_score >=
        MR_NULL_THRESHOLD
    )


  if (nrow(concordant_hits) == 0) {
    next
  }


  ## ==========================================================
  ## 8. Direct pharmacological targets
  ## ==========================================================

  direct_repurposing_hits <- list()
  direct_multimodal_targets <- list()


  for (
    hit_index in seq_len(
      nrow(concordant_hits)
    )
  ) {

    hit <- concordant_hits[
      hit_index,
    ]


    mr_direction <- sign(
      hit$mr_ivw_beta_pheno_1
    )


    ## --------------------------------------------------------
    ## Protective exposure:
    ## increasing exposure predicts lower disease risk.
    ## Seek activating/mimetic pharmacology.
    ## --------------------------------------------------------

    if (mr_direction == -1) {

      candidate_drugs <- opentargets_db[
        query_protein == hit$exposure &
          effect_on_target %in%
            c(
              "activates",
              "mimics"
            ) &
          phenotype %in%
            c(
              phenotype_1,
              phenotype_2
            )
      ]

      desired_moa <- "agonist"
    }


    ## --------------------------------------------------------
    ## Risk-increasing exposure:
    ## increasing exposure predicts greater disease risk.
    ## Seek inhibitory pharmacology.
    ## --------------------------------------------------------

    if (mr_direction == 1) {

      candidate_drugs <- opentargets_db[
        query_protein == hit$exposure &
          effect_on_target %in%
            c(
              "inhibits",
              "vaccine"
            ) &
          phenotype %in%
            c(
              phenotype_1,
              phenotype_2
            )
      ]

      desired_moa <- "antagonist"
    }


    if (mr_direction == 0) {
      next
    }


    ## --------------------------------------------------------
    ## Identify cross-indication repurposing opportunities
    ## --------------------------------------------------------

    if (nrow(candidate_drugs) > 0) {

      repurposing_hits <-
        candidate_drugs %>%
        group_by(
          OPENTARGETS_drug.id
        ) %>%
        summarise(
          phenotypes_present =
            list(
              unique(phenotype)
            ),
          n_phenotypes =
            n_distinct(
              phenotype
            ),
          .groups = "drop"
        ) %>%
        filter(
          n_phenotypes == 1
        ) %>%
        left_join(
          candidate_drugs,
          by =
            "OPENTARGETS_drug.id"
        )


      if (nrow(repurposing_hits) > 0) {

        repurposing_hits$exposure <-
          hit$exposure

        repurposing_hits$phenotype_1 <-
          phenotype_1

        repurposing_hits$phenotype_2 <-
          phenotype_2

        direct_repurposing_hits[
          [length(
            direct_repurposing_hits
          ) + 1]
        ] <- repurposing_hits

      } else {

        direct_multimodal_targets[
          [length(
            direct_multimodal_targets
          ) + 1]
        ] <- data.frame(
          target = hit$exposure,
          moa = desired_moa,
          phenotype_1 =
            phenotype_1,
          phenotype_2 =
            phenotype_2
        )
      }

    } else {

      direct_multimodal_targets[
        [length(
          direct_multimodal_targets
        ) + 1]
      ] <- data.frame(
        target = hit$exposure,
        moa = desired_moa,
        phenotype_1 =
          phenotype_1,
        phenotype_2 =
          phenotype_2
      )
    }
  }


  if (
    length(
      direct_repurposing_hits
    ) > 0
  ) {

    all_direct_repurposing_hits[
      [pair_id]
    ] <- bind_rows(
      direct_repurposing_hits
    )
  }


  if (
    length(
      direct_multimodal_targets
    ) > 0
  ) {

    all_direct_multimodal_targets[
      [pair_id]
    ] <- bind_rows(
      direct_multimodal_targets
    )
  }


  ## ==========================================================
  ## 9. OmniPath network-mediated pharmacological targets
  ## ==========================================================

  if (!INCLUDE_OMNIPATH) {
    next
  }


  network_repurposing_hits <- list()

  network_multimodal_targets <- list()


  for (
    hit_index in seq_len(
      nrow(concordant_hits)
    )
  ) {

    hit <- concordant_hits[
      hit_index,
    ]

    exposure <- hit$exposure

    mr_direction <- sign(
      hit$mr_ivw_beta_pheno_1
    )


    ## --------------------------------------------------------
    ## MR direction < 0
    ## --------------------------------------------------------

    if (mr_direction == -1) {


      ## Exposure inhibits downstream target

      network_targets <- unique(
        omnipath_db[
          source_genesymbol ==
            exposure &
            is_inhibition == 1,
          target_genesymbol
        ]
      )


      if (length(network_targets) > 0) {

        drugs <- opentargets_db[
          query_protein %in%
            network_targets &
            effect_on_target %in%
              c(
                "inhibits",
                "vaccine"
              ) &
            phenotype %in%
              c(
                phenotype_1,
                phenotype_2
              )
        ]


        if (nrow(drugs) > 0) {

          drugs$exposure <- exposure

          drugs$mechanistic_link <-
            paste(
              exposure,
              "inhibits"
            )

          network_repurposing_hits[
            [length(
              network_repurposing_hits
            ) + 1]
          ] <- drugs

        } else {

          network_multimodal_targets[
            [length(
              network_multimodal_targets
            ) + 1]
          ] <- data.frame(
            target =
              network_targets,
            moa =
              "antagonist",
            phenotype_1 =
              phenotype_1,
            phenotype_2 =
              phenotype_2,
            mechanistic_link =
              paste(
                exposure,
                "inhibits"
              ),
            exposure =
              exposure
          )
        }
      }


      ## Upstream activator of exposure

      network_targets <- unique(
        omnipath_db[
          target_genesymbol ==
            exposure &
            is_stimulation == 1,
          source_genesymbol
        ]
      )


      if (length(network_targets) > 0) {

        drugs <- opentargets_db[
          query_protein %in%
            network_targets &
            effect_on_target %in%
              c(
                "activates",
                "mimics"
              ) &
            phenotype %in%
              c(
                phenotype_1,
                phenotype_2
              )
        ]


        if (nrow(drugs) > 0) {

          drugs$exposure <- exposure

          drugs$mechanistic_link <-
            paste(
              exposure,
              "is activated by"
            )

          network_repurposing_hits[
            [length(
              network_repurposing_hits
            ) + 1]
          ] <- drugs

        } else {

          network_multimodal_targets[
            [length(
              network_multimodal_targets
            ) + 1]
          ] <- data.frame(
            target =
              network_targets,
            moa =
              "agonist",
            phenotype_1 =
              phenotype_1,
            phenotype_2 =
              phenotype_2,
            mechanistic_link =
              paste(
                exposure,
                "is activated by"
              ),
            exposure =
              exposure
          )
        }
      }
    }


    ## --------------------------------------------------------
    ## MR direction > 0
    ## --------------------------------------------------------

    if (mr_direction == 1) {


      ## Exposure activates downstream target

      network_targets <- unique(
        omnipath_db[
          source_genesymbol ==
            exposure &
            is_stimulation == 1,
          target_genesymbol
        ]
      )


      if (length(network_targets) > 0) {

        drugs <- opentargets_db[
          query_protein %in%
            network_targets &
            effect_on_target %in%
              c(
                "inhibits",
                "vaccine"
              ) &
            phenotype %in%
              c(
                phenotype_1,
                phenotype_2
              )
        ]


        if (nrow(drugs) > 0) {

          drugs$exposure <- exposure

          drugs$mechanistic_link <-
            paste(
              exposure,
              "activates"
            )

          network_repurposing_hits[
            [length(
              network_repurposing_hits
            ) + 1]
          ] <- drugs

        } else {

          network_multimodal_targets[
            [length(
              network_multimodal_targets
            ) + 1]
          ] <- data.frame(
            target =
              network_targets,
            moa =
              "antagonist",
            phenotype_1 =
              phenotype_1,
            phenotype_2 =
              phenotype_2,
            mechanistic_link =
              paste(
                exposure,
                "activates"
              ),
            exposure =
              exposure
          )
        }
      }


      ## Upstream inhibitor of exposure

      network_targets <- unique(
        omnipath_db[
          target_genesymbol ==
            exposure &
            is_inhibition == 1,
          source_genesymbol
        ]
      )


      if (length(network_targets) > 0) {

        drugs <- opentargets_db[
          query_protein %in%
            network_targets &
            effect_on_target %in%
              c(
                "activates",
                "mimics"
              ) &
            phenotype %in%
              c(
                phenotype_1,
                phenotype_2
              )
        ]


        if (nrow(drugs) > 0) {

          drugs$exposure <- exposure

          drugs$mechanistic_link <-
            paste(
              exposure,
              "is inhibited by"
            )

          network_repurposing_hits[
            [length(
              network_repurposing_hits
            ) + 1]
          ] <- drugs

        } else {

          network_multimodal_targets[
            [length(
              network_multimodal_targets
            ) + 1]
          ] <- data.frame(
            target =
              network_targets,
            moa =
              "agonist",
            phenotype_1 =
              phenotype_1,
            phenotype_2 =
              phenotype_2,
            mechanistic_link =
              paste(
                exposure,
                "is inhibited by"
              ),
            exposure =
              exposure
          )
        }
      }
    }
  }


  if (
    length(
      network_repurposing_hits
    ) > 0
  ) {

    all_network_repurposing_hits[
      [pair_id]
    ] <- bind_rows(
      network_repurposing_hits
    )
  }


  if (
    length(
      network_multimodal_targets
    ) > 0
  ) {

    all_network_multimodal_targets[
      [pair_id]
    ] <- bind_rows(
      network_multimodal_targets
    )
  }
}


## ============================================================
## 10. Combine outputs
## ============================================================

mr_output <- bind_rows(
  all_mr_results,
  .id = "phenotype_pair"
)

direct_repurposing_output <- bind_rows(
  all_direct_repurposing_hits,
  .id = "phenotype_pair"
)

direct_multimodal_output <- bind_rows(
  all_direct_multimodal_targets,
  .id = "phenotype_pair"
)

network_repurposing_output <- bind_rows(
  all_network_repurposing_hits,
  .id = "phenotype_pair"
)

network_multimodal_output <- bind_rows(
  all_network_multimodal_targets,
  .id = "phenotype_pair"
)


## ============================================================
## 11. Save
## ============================================================

results <- list(

  parameters = list(
    molecular_qtl =
      MOLECULAR_QTL,
    lava_p_threshold =
      LAVA_P_THRESHOLD,
    coloc_pph4_threshold =
      COLOC_PPH4_THRESHOLD,
    string_score_threshold =
      STRING_SCORE_THRESHOLD,
    mr_concordance_threshold =
      MR_CONCORDANCE_THRESHOLD,
    mr_null_threshold =
      MR_NULL_THRESHOLD,
    include_string =
      INCLUDE_STRING,
    include_omnipath =
      INCLUDE_OMNIPATH
  ),

  mr_results =
    mr_output,

  direct_repurposing_hits =
    direct_repurposing_output,

  direct_multimodal_targets =
    direct_multimodal_output,

  network_repurposing_hits =
    network_repurposing_output,

  network_multimodal_targets =
    network_multimodal_output
)


saveRDS(
  results,
  file = output_file
)

message(
  "Target prioritisation complete."
)

message(
  "Results saved to: ",
  output_file
)