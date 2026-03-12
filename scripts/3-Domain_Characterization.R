#!/usr/bin/env Rscript

# =============================================================================
# DOMAIN-SPECIFIC SEVERITY SCORE CHARACTERIZATION
# PPMI Cohort - Unsupervised Transcriptomic Severity Score (UTSS)
# =============================================================================
#
# Description: This script processes the severity signals generated from
#              pathway-specific clustering to compute the Unsupervised
#              Transcriptomic Severity Score (UTSS) for motor and non-motor
#              domains. It aggregates signals across pathways, combines data
#              across visits, and integrates with clinical covariates.
#
#              The UTSS is calculated as the sum of severity signals across
#              pathways within each domain:
#              - UTSS_Motor: Sum of signals from motor-related covariates
#              - UTSS_NonMotor: Sum of signals from non-motor-related covariates
#
# Author: [Jaime Ñíguez Baeza]
# Date: 2026
# Corresponding author: [juanbot@um.es]
#
# Input files:
#   - Sev_score_detail_M*_pathways_Reactome_pathways_filtered_50_protein_coding.rds:
#     Severity signals for each visit (M0, M12, M24, M36)
#   - PPMI_label_forScript: Clinical labels with visit information
#   - PPMI_covariates_extra: Additional clinical covariates
#   - Covariates_new_all: Extended covariate set
#
# Output files:
#   - Covariates_with_UTSS_scores.rds: Combined dataset with UTSS scores
#     and clinical covariates
#
# Domain definitions (based on established clinical frameworks):
#   Motor domain: Hoehn & Yahr stage, UPDRS II/III, tremor, PIGD, Schwab & England
#   Non-motor domain: UPDRS I, cognitive impairment, hallucinations,
#                     depression, anxiety, apathy, RBD, ESS, MoCA
#
# Dependencies: tidyverse, dplyr, stringr
# =============================================================================

# =============================================================================
# 1. INITIAL SETUP AND DATA LOADING
# =============================================================================

# Load required libraries
library(tidyverse)
library(dplyr)
library(stringr)

# Set seed for reproducibility
set.seed(123)

# Define base directory (modify as needed)
base_dir <- "/home/jaimeniguez/Pipeline"

cat("========================================\n")
cat("DOMAIN-SPECIFIC SEVERITY SCORE CHARACTERIZATION\n")
cat("Unsupervised Transcriptomic Severity Score (UTSS)\n")
cat("========================================\n\n")

# -----------------------------------------------------------------------------
# 1.1 Load clinical labels
# -----------------------------------------------------------------------------
cat("1.1 Loading clinical labels...\n")
PPMI_label <- readRDS(file.path(base_dir, "PPMI_label_forScript"))
cat("  - Total samples in labels:", nrow(PPMI_label), "\n")
cat("  - Available visits:", paste(unique(PPMI_label$visit_name), collapse = ", "), "\n\n")

# -----------------------------------------------------------------------------
# 1.2 Load severity signals for each visit
# -----------------------------------------------------------------------------
cat("1.2 Loading severity signals for each visit...\n")

results_dir <- file.path(base_dir, "results_pathways")

# Load M0 signals
Sev_score_detail_M0 <- readRDS(
  file.path(results_dir, "Sev_score_detail_M0_pathways_Reactome_pathways_filtered_50_protein_coding.rds")
)
cat("  - M0 signals loaded:", nrow(Sev_score_detail_M0), "records\n")

# Load M12 signals
Sev_score_detail_M12 <- readRDS(
  file.path(results_dir, "Sev_score_detail_M12_pathways_Reactome_pathways_filtered_50_protein_coding.rds")
)
cat("  - M12 signals loaded:", nrow(Sev_score_detail_M12), "records\n")

# Load M24 signals
Sev_score_detail_M24 <- readRDS(
  file.path(results_dir, "Sev_score_detail_M24_pathways_Reactome_pathways_filtered_50_protein_coding.rds")
)
cat("  - M24 signals loaded:", nrow(Sev_score_detail_M24), "records\n")

# Load M36 signals
Sev_score_detail_M36 <- readRDS(
  file.path(results_dir, "Sev_score_detail_M36_pathways_Reactome_pathways_filtered_50_protein_coding.rds")
)
cat("  - M36 signals loaded:", nrow(Sev_score_detail_M36), "records\n\n")

# =============================================================================
# 2. GENERATION OF MOTOR AND NON-MOTOR SEVERITY SCORES
# =============================================================================
cat("========================================\n")
cat("2. GENERATION OF MOTOR AND NON-MOTOR SEVERITY SCORES\n")
cat("========================================\n\n")

#' Create severity score summary for a specific visit
#'
#' Aggregates severity signals across pathways for each patient and covariate,
#' then calculates total scores. This function transforms the detailed
#' pathway-level signals into patient-level scores per clinical covariate.
#'
#' @param Sev_score_detail_MSigDB Data frame with severity signals from clustering
#' @param PPMI_label Clinical labels with visit information
#' @param visit_label Visit identifier (e.g., "M0", "M12")
#' @return Data frame with severity scores per patient for the specified visit
create_score_summary <- function(Sev_score_detail_MSigDB, PPMI_label, visit_label) {
  
  # Convert rownames to "sample_id" column for joining
  PPMI_label <- PPMI_label %>%
    tibble::rownames_to_column(var = "sample_id")
  
  # Extract sample_ids for the specific visit
  sample_ids_visit <- PPMI_label %>%
    filter(visit_name == visit_label) %>%
    dplyr::select(sample_id, visit_name)
  
  # Sum binary scores (SevScore) by sample_id and Covariate
  # This aggregates across all pathways: each patient receives a count
  # of how many pathways identified them as being in the more severe cluster
  # for each clinical covariate
  score_orig <- Sev_score_detail_MSigDB %>%
    group_by(sample_id, Covariate) %>%
    summarise(score_sum = sum(SevScore, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(
      names_from = Covariate,
      values_from = score_sum,
      names_prefix = "Score_",
      values_fill = list(score_sum = 0)
    )
  
  # Add patients without any scores and fill with zeros
  # This ensures all patients from the visit are represented
  full_score_summary <- sample_ids_visit %>%
    left_join(score_orig, by = "sample_id") %>%
    mutate(across(where(is.numeric), ~replace_na(., 0)))
  
  return(full_score_summary)
}

# -----------------------------------------------------------------------------
# 2.1 Generate score summaries for each visit
# -----------------------------------------------------------------------------
cat("2.1 Generating score summaries for each visit...\n")

cat("  - Processing M0...\n")
Severity_score_summ_M0 <- create_score_summary(
  Sev_score_detail_M0, PPMI_label, visit_label = "M0"
) %>%
  mutate(
    Visit = 0,
    participant_id = as.numeric(sub("PP-(\\d+)-.*", "\\1", sample_id))
  )

cat("  - Processing M12...\n")
Severity_score_summ_M12 <- create_score_summary(
  Sev_score_detail_M12, PPMI_label, visit_label = "M12"
) %>%
  mutate(
    Visit = 12,
    participant_id = as.numeric(sub("PP-(\\d+)-.*", "\\1", sample_id))
  )

cat("  - Processing M24...\n")
Severity_score_summ_M24 <- create_score_summary(
  Sev_score_detail_M24, PPMI_label, visit_label = "M24"
) %>%
  mutate(
    Visit = 24,
    participant_id = as.numeric(sub("PP-(\\d+)-.*", "\\1", sample_id))
  )

cat("  - Processing M36...\n")
Severity_score_summ_M36 <- create_score_summary(
  Sev_score_detail_M36, PPMI_label, visit_label = "M36"
) %>%
  mutate(
    Visit = 36,
    participant_id = as.numeric(sub("PP-(\\d+)-.*", "\\1", sample_id))
  )

cat("✅ Score summaries generated for all visits\n\n")

# -----------------------------------------------------------------------------
# 2.2 Combine all visits
# -----------------------------------------------------------------------------
cat("2.2 Combining scores across all visits...\n")

list_of_scores <- list(
  Severity_score_summ_M0, 
  Severity_score_summ_M12, 
  Severity_score_summ_M24, 
  Severity_score_summ_M36
)

# Combine dataframes with full join to preserve all patients and visits
combined_severity_score <- purrr::reduce(list_of_scores, full_join, by = NULL)

# Replace NA with 0 for all scores (patients with no signals in a visit)
combined_severity_score[is.na(combined_severity_score)] <- 0

cat("  - Combined data dimensions:", dim(combined_severity_score)[1], "rows ×", 
    dim(combined_severity_score)[2], "columns\n")
cat("  - Unique patients:", length(unique(combined_severity_score$participant_id)), "\n")
cat("  - Visits:", paste(unique(combined_severity_score$Visit), collapse = ", "), "\n\n")

# =============================================================================
# 3. DOMAIN DEFINITION AND UTSS CALCULATION
# =============================================================================
cat("========================================\n")
cat("3. DOMAIN DEFINITION AND UTSS CALCULATION\n")
cat("========================================\n\n")

# -----------------------------------------------------------------------------
# 3.1 Define motor and non-motor domains
# -----------------------------------------------------------------------------
cat("3.1 Defining motor and non-motor domains based on clinical frameworks...\n")

# Motor domain covariates
# These variables reflect motor symptoms and functional impairment
Motor <- c(
  "code_upd2hy_hoehn_and_yahr_stage",      # Hoehn & Yahr stage
  "mds_updrs_part_ii_summary_score",       # UPDRS Part II (motor experiences of daily living)
  "mds_updrs_part_iii_summary_score",      # UPDRS Part III (motor examination)
  "Tremor_Right_UM",                        # Right upper limb tremor
  "Tremor_Left_UM",                         # Left upper limb tremor
  "Tremor_all_Um",                          # Overall tremor score
  "PIGD_UM",                                 # Postural instability and gait difficulty
  "mod_schwab_england_pct_adl_score"        # Schwab & England ADL scale
)

# Non-motor domain covariates
# These variables reflect cognitive, psychiatric, and autonomic symptoms
Non_Motor <- c(
  "code_upd2106_dopamine_dysregulation_syndrome_features",  # DDS features
  "mds_updrs_part_i_summary_score",                          # UPDRS Part I (non-motor experiences)
  "code_upd2102_hallucinations_and_psychosis",               # Hallucinations/psychosis
  "code_upd2101_cognitive_impairment",                       # Cognitive impairment
  "code_upd2103_depressed_mood",                             # Depression
  "code_upd2104_anxious_mood",                               # Anxiety
  "code_upd2105_apathy",                                      # Apathy
  "rbd_summary_score",                                        # RBD screening
  "ess_summary_score",                                        # Epworth Sleepiness Scale
  "moca_total_score"                                          # Montreal Cognitive Assessment
)

# Column names for the raw scores (prefixed with "Score_" from pivot_wider)
Motor_raw <- paste0("Score_", Motor)
Non_Motor_raw <- paste0("Score_", Non_Motor)

cat("  - Motor covariates:", length(Motor), "\n")
cat("    ", paste(Motor, collapse = ", "), "\n\n")
cat("  - Non-motor covariates:", length(Non_Motor), "\n")
cat("    ", paste(Non_Motor, collapse = ", "), "\n\n")

# -----------------------------------------------------------------------------
# 3.2 Calculate UTSS for motor and non-motor domains
# -----------------------------------------------------------------------------
cat("3.2 Calculating Unsupervised Transcriptomic Severity Scores (UTSS)...\n")

# Calculate UTSS for motor and non-motor domains by summing signals
# across all covariates within each domain
combined_severity_score_all_visits <- combined_severity_score %>%
  mutate(
    # Unsupervised Transcriptomic Severity Score (UTSS)
    UTSS_Motor = rowSums(dplyr::select(., any_of(Motor_raw)), na.rm = TRUE),
    UTSS_NonMotor = rowSums(dplyr::select(., any_of(Non_Motor_raw)), na.rm = TRUE)
  )

# Display summary statistics
cat("\nUTSS Summary Statistics:\n")
cat("  UTSS_Motor - Range:", 
    min(combined_severity_score_all_visits$UTSS_Motor, na.rm = TRUE), "to",
    max(combined_severity_score_all_visits$UTSS_Motor, na.rm = TRUE), "\n")
cat("  UTSS_NonMotor - Range:", 
    min(combined_severity_score_all_visits$UTSS_NonMotor, na.rm = TRUE), "to",
    max(combined_severity_score_all_visits$UTSS_NonMotor, na.rm = TRUE), "\n\n")

# -----------------------------------------------------------------------------
# 3.3 Generate cumulative scores to track progression
# -----------------------------------------------------------------------------
cat("3.3 Generating cumulative scores to track progression...\n")

# Calculate cumulative sums across visits for each patient
# This allows analysis of disease progression over time
combined_severity_score_all_visits_acc <- combined_severity_score_all_visits %>%
  arrange(participant_id, Visit) %>%
  group_by(participant_id) %>%
  mutate(across(
    c(UTSS_Motor, UTSS_NonMotor),
    cumsum,
    .names = "{.col}_Cumulative"
  )) %>%
  ungroup()

cat("✅ UTSS scores calculated for motor and non-motor domains\n")
cat("   - Cumulative scores generated for progression analysis\n\n")

# =============================================================================
# 4. INTEGRATION WITH ADDITIONAL CLINICAL COVARIATES
# =============================================================================
cat("========================================\n")
cat("4. INTEGRATION WITH ADDITIONAL CLINICAL COVARIATES\n")
cat("========================================\n\n")

# -----------------------------------------------------------------------------
# 4.1 Load and prepare additional covariates
# -----------------------------------------------------------------------------
cat("4.1 Loading additional clinical covariates...\n")

# Load main covariates dataset
Covariates <- as.data.frame(
  readRDS(file.path(base_dir, "metadata_reduced", "PPMI_covariates_extra"))
)
rownames(Covariates) <- Covariates$sample_id

cat("  - Main covariates loaded:", nrow(Covariates), "records\n")

# Load additional covariates (age at onset, prognosis, genetics, etc.)
cna <- readRDS(file.path(base_dir, "Pipeline", "Covariates_new_all"))
cat("  - Extended covariates loaded:", nrow(cna), "records\n")

# Merge with main covariates
Covariates <- Covariates %>%
  left_join(
    cna %>% dplyr::select(
      sample_id, 
      AGE_ONSET, 
      prognosis_status, 
      apoe, 
      mapt, 
      dementia_status, 
      MOCA_Status
    ), 
    by = "sample_id"
  )

# Propagate time-invariant variables within participants
# (variables that don't change over time are filled across visits)
Covariates <- Covariates %>%
  group_by(participant_id) %>%
  fill(AGE_ONSET, prognosis_status, apoe, mapt, .direction = "downup") %>%
  ungroup()

cat("  - Final covariates dimensions:", dim(Covariates)[1], "rows ×", 
    dim(Covariates)[2], "columns\n\n")

# -----------------------------------------------------------------------------
# 4.2 Filter covariates for patients with UTSS scores
# -----------------------------------------------------------------------------
cat("4.2 Filtering covariates for patients with UTSS scores...\n")

# Identify patients with UTSS scores
common_patients <- unique(combined_severity_score_all_visits$participant_id)
common_patients_id <- paste("PP-", common_patients, sep = "")

# Filter covariates to include only these patients and exclude M6 visit
Covariates_filtered <- Covariates %>% 
  filter(participant_id %in% common_patients_id & visit_name != "M6") %>%
  mutate(Visit = as.numeric(str_extract(visit_name, "\\d+")))

cat("  - Patients with UTSS scores:", length(common_patients), "\n")
cat("  - Covariates records after filtering:", nrow(Covariates_filtered), "\n")
cat("  - Visits in filtered covariates:", 
    paste(unique(Covariates_filtered$visit_name), collapse = ", "), "\n\n")

# -----------------------------------------------------------------------------
# 4.3 Merge UTSS scores with clinical covariates
# -----------------------------------------------------------------------------
cat("4.3 Merging UTSS scores with clinical covariates...\n")

# Define columns to merge
score_columns <- c(
  "UTSS_Motor", 
  "UTSS_NonMotor",
  "UTSS_Motor_Cumulative", 
  "UTSS_NonMotor_Cumulative"
)

# Merge datasets
Covariates_score <- Covariates_filtered %>%
  left_join(
    combined_severity_score_all_visits_acc %>% 
      dplyr::select("sample_id", all_of(score_columns)), 
    by = "sample_id"
  )

cat("  - Final merged dataset dimensions:", 
    dim(Covariates_score)[1], "rows ×", dim(Covariates_score)[2], "columns\n")

# Display summary of UTSS scores in the merged dataset
cat("\nUTSS scores in final dataset:\n")
print(summary(Covariates_score[, c("participant_id", "visit_name", score_columns)]))

# -----------------------------------------------------------------------------
# 4.4 Save final dataset
# -----------------------------------------------------------------------------
cat("\n4.4 Saving final dataset with UTSS scores...\n")

saveRDS(
  Covariates_score, 
  file.path(results_dir, "Covariates_with_UTSS_scores.rds")
)

cat("  - Dataset saved to:", file.path(results_dir, "Covariates_with_UTSS_scores.rds"), "\n")
cat("✅ Final dataset saved with UTSS motor and non-motor scores\n\n")

# =============================================================================
# 5. SUMMARY STATISTICS
# =============================================================================
cat("========================================\n")
cat("5. SUMMARY STATISTICS\n")
cat("========================================\n\n")

# Calculate overall statistics
total_patients <- length(unique(Covariates_score$participant_id))
total_visits <- nrow(Covariates_score)
visits_per_patient <- Covariates_score %>%
  group_by(participant_id) %>%
  summarise(n_visits = n())

cat("Final dataset summary:\n")
cat("  - Total unique patients:", total_patients, "\n")
cat("  - Total observations (patient-visits):", total_visits, "\n")
cat("  - Average visits per patient:", 
    round(mean(visits_per_patient$n_visits), 2), "\n")
cat("  - Visit distribution:\n")
print(table(Covariates_score$visit_name))

cat("\nUTSS score statistics by visit:\n")
Covariates_score %>%
  group_by(visit_name) %>%
  summarise(
    n_patients = n(),
    mean_UTSS_Motor = mean(UTSS_Motor, na.rm = TRUE),
    sd_UTSS_Motor = sd(UTSS_Motor, na.rm = TRUE),
    mean_UTSS_NonMotor = mean(UTSS_NonMotor, na.rm = TRUE),
    sd_UTSS_NonMotor = sd(UTSS_NonMotor, na.rm = TRUE)
  ) %>%
  print()

# =============================================================================
# 6. SESSION INFORMATION
# =============================================================================
cat("\n========================================\n")
cat("6. SESSION INFORMATION\n")
cat("========================================\n\n")

print(sessionInfo())

cat("\n========================================\n")
cat("Domain characterization completed successfully!\n")
cat("UTSS motor and non-motor scores have been generated and integrated with clinical covariates.\n")
cat("========================================\n")