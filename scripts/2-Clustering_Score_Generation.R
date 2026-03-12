#!/usr/bin/env Rscript

# =============================================================================
# PATHWAY-SPECIFIC TRANSCRIPTOMIC CLUSTERING AND SEVERITY SCORE GENERATION
# PPMI Cohort - Unsupervised Transcriptomic Severity Score (UTSS)
# =============================================================================
#
# Description: This script performs pathway-specific hierarchical clustering
#              of patients and generates severity scores for motor and 
#              non-motor domains based on pairwise cluster comparisons.
#
#              For each pathway (n=336), patients are clustered using Ward's
#              method with absolute Pearson correlation distance. Clusters
#              are compared pairwise across clinical covariates to identify
#              significant differences in disease severity. Patients in more
#              severe clusters receive +1 scores, which are summed across
#              pathways within motor and non-motor domains to generate UTSS.
#
# Author: [Jaime Ñíguez Baeza]
# Date: 2026
# Corresponding author: [juanbot@um.es]
#
# Input files:
#   - Subsets_pathways_Reactome: Reactome pathway gene sets
#   - pathways_filtered_similarity_50: Filtered pathway list (n=336)
#   - PPMI_label_forScript: Clinical covariates with labels
#   - PPMI_exp: Normalized expression matrix
#   - protein_coding_genes: List of protein-coding genes
#   - Covariates_new_all: Additional clinical covariates
#
# Output files:
#   - Sev_score_detail_M*_pathways_Reactome_pathways_filtered_50_*.rds:
#     Detailed severity signals per pathway, cluster, and patient
#   - Combined severity scores with motor/non-motor domains
#
# Clinical covariates used for severity assessment:
#   Risk_var: Covariates where higher values indicate greater severity
#   Inv_Var: Covariates where lower values indicate greater severity
#            (moca_total_score, mod_schwab_england_pct_adl_score)
#
# Dependencies: Various (see library calls)
# =============================================================================

# =============================================================================
# 1. INITIAL SETUP AND DATA LOADING
# =============================================================================

# Load required libraries
library(tidyverse)
library(dplyr)
library(ggplot2)
library(patchwork)
library(DT)

# Clustering and statistics
library(cluster)
library(fpc)
library(clValid)
library(clusterSim)
library(factoextra)
library(rstatix)
library(dendextend)
library(pheatmap)

# Machine learning (for compatibility, not actively used)
library(caret)
library(glmnet)
library(ranger)
library(randomForest)
library(MLmetrics)

# Parallel processing
library(foreach)
library(doParallel)

# Visualization
library(plotly)
library(circlize)
library(networkD3)
library(ggpubr)
library(RColorBrewer)
library(forcats)

# Gene annotation and pathways
library(biomaRt)
library(msigdbr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(gprofiler2)
library(stringr)

# Multivariate analysis
library(FactoMineR)
library(e1071)

# Set seed for reproducibility
set.seed(123)

# Define base directory (modify as needed)
base_dir <- "/home/jaimeniguez/Pipeline"

cat("========================================\n")
cat("Starting pathway clustering and severity score generation\n")
cat("========================================\n\n")

# -----------------------------------------------------------------------------
# 1.1 Load pathway gene sets
# -----------------------------------------------------------------------------
cat("1.1 Loading pathway gene sets...\n")

# Load Reactome pathways
Subsets_pathways <- readRDS(file.path(base_dir, "Subsets_pathways_Reactome"))

# Load filtered pathways (336 non-redundant pathways after clustering)
pathways_filtered_similarity <- readRDS(file.path(base_dir, "pathways_filtered_similarity_50"))
Subsets_pathways <- Subsets_pathways[pathways_filtered_similarity]

cat("  - Total pathways loaded:", length(Subsets_pathways), "\n")

# Display pathway size statistics
pathway_sizes <- sapply(Subsets_pathways, length)
cat("  - Pathway size range:", min(pathway_sizes), "to", max(pathway_sizes), "genes\n")
cat("  - Median pathway size:", median(pathway_sizes), "genes\n\n")

# -----------------------------------------------------------------------------
# 1.2 Load clinical and expression data
# -----------------------------------------------------------------------------
cat("1.2 Loading clinical and expression data...\n")

# Load clinical labels
PPMI_label <- readRDS(file.path(base_dir, "PPMI_label_forScript"))

# Load expression data
PPMI_exp <- readRDS("~/1aFase_ReducirGenes/PPMI_exp")

# Load protein-coding genes and filter expression matrix
protein_coding_genes <- readRDS(file.path(base_dir, "protein_coding_genes"))
PPMI_exp <- PPMI_exp[, protein_coding_genes]

cat("  - Clinical samples:", nrow(PPMI_label), "\n")
cat("  - Expression matrix dimensions:", dim(PPMI_exp)[1], "samples ×", 
    dim(PPMI_exp)[2], "protein-coding genes\n\n")

# =============================================================================
# 2. CLINICAL COVARIATES FOR SEVERITY ASSESSMENT
# =============================================================================
cat("========================================\n")
cat("2. CLINICAL COVARIATES FOR SEVERITY ASSESSMENT\n")
cat("========================================\n\n")

# Covariates where higher values indicate greater severity
Risk_var <- c(
  "mds_updrs_part_i_summary_score",
  "code_upd2101_cognitive_impairment",
  "code_upd2102_hallucinations_and_psychosis",
  "code_upd2103_depressed_mood",
  "code_upd2104_anxious_mood",
  "code_upd2105_apathy",
  "code_upd2106_dopamine_dysregulation_syndrome_features",
  "mds_updrs_part_ii_summary_score",
  "mds_updrs_part_iii_summary_score",
  "Tremor_Right_UM",
  "Tremor_Left_UM",
  "Tremor_all_Um",
  "PIGD_UM",
  "ess_summary_score",
  "moca_total_score",
  "Abeta",
  "p-Tau",
  "Tau",
  "rbd_summary_score",
  "code_upd2hy_hoehn_and_yahr_stage",
  "apoe",
  "mapt",
  "mod_schwab_england_pct_adl_score"
)

# Covariates where lower values indicate greater severity (inverse direction)
Inv_Var <- c("moca_total_score", "mod_schwab_england_pct_adl_score")

cat("2.1 Clinical covariates loaded:\n")
cat("  - Total risk covariates:", length(Risk_var), "\n")
cat("  - Inverse direction covariates:", length(Inv_Var), "\n\n")

# =============================================================================
# 3. PATHWAY-SPECIFIC CLUSTERING AND PAIRWISE COMPARISONS
# =============================================================================
cat("========================================\n")
cat("3. PATHWAY-SPECIFIC CLUSTERING AND PAIRWISE COMPARISONS\n")
cat("========================================\n\n")

#' Perform hierarchical clustering and pairwise cluster comparisons for pathways
#' 
#' For each pathway, patients are clustered using hierarchical clustering with
#' Ward's method and absolute Pearson correlation distance. Clusters are compared
#' pairwise across clinical covariates to identify significant differences in
#' disease severity. Patients in more severe clusters receive +1 scores.
#'
#' @param df_covariates Data frame with clinical covariates (samples × covariates)
#' @param df_transcript Expression matrix (samples × genes)
#' @param Subsets_pathways List of pathway gene sets
#' @param n_clust Number of clusters (k) for hierarchical clustering
#' @param use_random_pathways Logical; if TRUE, use random gene sets of same size
#' @param seed Random seed for reproducibility
#' @return Data frame with severity signals (pathway, cluster, patient, covariate)
dividing_data_cluster <- function(df_covariates, df_transcript, 
                                  Subsets_pathways, n_clust, 
                                  use_random_pathways = FALSE,
                                  seed = NULL) {
  
  # Initialize results dataframe
  Sev_score_detail <- data.frame()
  
  if (!is.null(seed)) set.seed(seed)
  
  pathways_to_process <- list()
  
  # Prepare pathways (original or random)
  if (use_random_pathways) {
    # Generate random gene sets of same size as original pathways
    all_valid_genes <- colnames(df_transcript)
    
    for (pathway_name in names(Subsets_pathways)) {
      original_genes <- Subsets_pathways[[pathway_name]][
        Subsets_pathways[[pathway_name]] %in% all_valid_genes
      ]
      n_genes_original <- length(original_genes)
      if (n_genes_original < 2) next
      random_genes <- sample(all_valid_genes, n_genes_original)
      random_pathway_name <- paste0("RANDOM_", pathway_name)
      pathways_to_process[[random_pathway_name]] <- random_genes
    }
    
  } else {
    # Use original pathways, keeping only genes present in expression matrix
    for (pathway_name in names(Subsets_pathways)) {
      valid_genes <- Subsets_pathways[[pathway_name]][
        Subsets_pathways[[pathway_name]] %in% colnames(df_transcript)
      ]
      if (length(valid_genes) >= 2) {
        pathways_to_process[[pathway_name]] <- valid_genes
      }
    }
  }
  
  if (length(pathways_to_process) == 0) {
    stop("No valid pathways to process")
  }
  
  cat("\nTotal pathways to process:", length(pathways_to_process), "\n")
  cat(rep("-", 50), "\n")
  
  pathway_counter <- 0
  
  # Process each pathway
  for (pathway_name in names(pathways_to_process)) {
    
    pathway_counter <- pathway_counter + 1
    
    # Extract expression data for this pathway
    valid_columns <- pathways_to_process[[pathway_name]]
    df_transcript_pathway <- df_transcript[, valid_columns, drop = FALSE]
    
    # Remove constant or NA columns
    is_constant_or_na <- apply(df_transcript_pathway, 2, 
                               function(x) all(is.na(x)) || sd(x, na.rm = TRUE) == 0)
    df_transcript_pathway <- df_transcript_pathway[, !is_constant_or_na, drop = FALSE]
    
    if (ncol(df_transcript_pathway) < 2) next
    
    # Calculate correlation-based distance matrix
    # Similarity: absolute Pearson correlation
    # Dissimilarity: 1 - |correlation|
    cor.pe <- cor(t(as.matrix(df_transcript_pathway)), method = "pearson")
    cor.pe[cor.pe >= 0.9999] <- 0.9999
    cor.pe[cor.pe <= -0.9999] <- -0.9999
    DissMatrix <- as.dist(1 - abs(cor.pe))
    
    tryCatch({
      
      # Hierarchical clustering with Ward's method
      hclust_object <- hclust(DissMatrix, method = "ward.D2")
      cluster_assignments <- cutree(hclust_object, k = n_clust)
      df_covariates$cluster <- cluster_assignments
      full_cluster <- unique(cluster_assignments)
      signals_found <- 0
      
      # Pairwise comparisons between all clusters
      for (i in 1:(length(full_cluster)-1)) {
        for (j in (i+1):length(full_cluster)) {
          
          cluster_i <- df_covariates[cluster_assignments == i, ]
          cluster_j <- df_covariates[cluster_assignments == j, ]
          variables <- setdiff(colnames(cluster_i), "cluster")
          
          for (variable in variables) {
            
            if (!variable %in% Risk_var) next
            
            # ================= NUMERIC VARIABLES =================
            if (is.numeric(df_covariates[[variable]])) {
              
              # Mann-Whitney U test
              wt <- wilcox.test(cluster_i[[variable]], cluster_j[[variable]],
                                exact = min(nrow(cluster_i), nrow(cluster_j)) < 20)
              
              # Bonferroni correction for multiple comparisons
              wt$p.value <- p.adjust(wt$p.value,
                                     method = "bonferroni",
                                     n = n_clust * length(variables))
              
              if (!is.na(wt$p.value) && wt$p.value <= 0.05) {
                
                m_i <- median(cluster_i[[variable]], na.rm = TRUE)
                m_j <- median(cluster_j[[variable]], na.rm = TRUE)
                
                # Determine which cluster has more severe profile
                if (variable %in% Inv_Var) {
                  # For inverse variables, lower values indicate more severity
                  worsening_cluster <- ifelse(m_i < m_j, i, j)
                } else {
                  # For regular variables, higher values indicate more severity
                  worsening_cluster <- ifelse(m_i > m_j, i, j)
                }
                
                target_cluster <- if (worsening_cluster == i) cluster_i else cluster_j
                target_vs <- if (worsening_cluster == i) j else i
                
                # Create severity signal record
                row <- data.frame(
                  Pathway = pathway_name,
                  Cluster = worsening_cluster,
                  k = n_clust,
                  n = nrow(target_cluster),
                  Covariate = variable,
                  sample_id = rownames(target_cluster),
                  p_value = wt$p.value,
                  VS = target_vs,
                  SevScore = 1,  # Binary score
                  Pathway_Type = ifelse(grepl("^RANDOM_", pathway_name),
                                        "Random", "Original")
                )
                
                Sev_score_detail <- rbind(Sev_score_detail, row)
                signals_found <- signals_found + 1
              }
              
              # ================= CATEGORICAL VARIABLES =================
            } else if (is.factor(df_covariates[[variable]])) {
              
              full_classes <- sort(unique(c(levels(df_covariates[[variable]]),
                                            as.character(cluster_i[[variable]]), 
                                            as.character(cluster_j[[variable]]))))
              is_inverse <- variable %in% Inv_Var
              
              for (c1 in full_classes) {
                for (c2 in full_classes) {
                  if (c1 >= c2) next
                  
                  # For categorical variables, determine which class represents more severity
                  upper_class <- ifelse(is_inverse,
                                        min(as.numeric(c1), as.numeric(c2)),
                                        max(as.numeric(c1), as.numeric(c2)))
                  lower_class <- ifelse(is_inverse,
                                        max(as.numeric(c1), as.numeric(c2)),
                                        min(as.numeric(c1), as.numeric(c2)))
                  
                  # Create contingency table
                  tab <- matrix(0, 2, 2,
                                dimnames = list(Class = c(upper_class, lower_class),
                                                Group = c(paste0("C", i), paste0("C", j))))
                  tab[1,1] <- sum(cluster_i[[variable]] == upper_class, na.rm = TRUE)
                  tab[2,1] <- sum(cluster_i[[variable]] == lower_class, na.rm = TRUE)
                  tab[1,2] <- sum(cluster_j[[variable]] == upper_class, na.rm = TRUE)
                  tab[2,2] <- sum(cluster_j[[variable]] == lower_class, na.rm = TRUE)
                  
                  if (all(rowSums(tab) > 0) && all(colSums(tab) > 0)) {
                    # Choose appropriate test based on expected frequencies
                    test_to_use <- if (any(chisq.test(tab)$expected < 5)) fisher.test else chisq.test
                    res <- test_to_use(tab)
                    
                    # Bonferroni correction
                    res$p.value <- p.adjust(res$p.value,
                                            method = "bonferroni",
                                            n = n_clust * length(variables))
                    
                    prop_upper_i <- tab[1,1] / sum(tab[,1])
                    prop_upper_j <- tab[1,2] / sum(tab[,2])
                    
                    if (!is.na(res$p.value) && res$p.value <= 0.05) {
                      
                      # Determine which cluster has higher proportion in more severe class
                      worsening_cluster <- if (prop_upper_i > prop_upper_j) i else j
                      target_cluster <- if (worsening_cluster == i) cluster_i else cluster_j
                      target_vs <- if (worsening_cluster == i) j else i
                      
                      # Create severity signal record
                      row <- data.frame(
                        Pathway = pathway_name,
                        Cluster = worsening_cluster,
                        k = n_clust,
                        n = nrow(target_cluster),
                        Covariate = variable,
                        sample_id = rownames(target_cluster),
                        p_value = res$p.value,
                        VS = target_vs,
                        SevScore = 1,  # Binary score
                        Pathway_Type = ifelse(grepl("^RANDOM_", pathway_name),
                                              "Random", "Original")
                      )
                      
                      Sev_score_detail <- rbind(Sev_score_detail, row)
                      signals_found <- signals_found + 1
                    }
                  }
                }
              }
            }
          }
        }
      }
      
    }, error = function(e) {
      warning("Error in pathway ", pathway_name, ": ", e$message)
    })
    
  } # End pathway loop
  
  cat("\nANALYSIS COMPLETED\n")
  cat("Total records:", nrow(Sev_score_detail), "\n")
  
  return(Sev_score_detail)
}

# =============================================================================
# 4. CLUSTERING ANALYSIS FOR EACH VISIT
# =============================================================================
cat("========================================\n")
cat("4. CLUSTERING ANALYSIS FOR EACH VISIT\n")
cat("========================================\n\n")

# The same clustering strategy is applied to all visits (M0, M12, M24, M36)
# For each visit, patients are clustered using k values from 2 to 20

# -----------------------------------------------------------------------------
# 4.1 Prepare data for M0 visit
# -----------------------------------------------------------------------------
cat("4.1 Preparing data for M0 visit...\n")

PPMI_label_M0 <- PPMI_label %>% filter(visit_name == "M0")
PPMI_label_M0$visit_name <- NULL

# Ensure consistent samples between clinical and expression data
common_ids <- intersect(rownames(PPMI_label_M0), rownames(PPMI_exp))

PPMI_label_M0 <- PPMI_label_M0[common_ids, ]
PPMI_exp_M0   <- PPMI_exp[common_ids, ]

cat("  - M0 samples:", nrow(PPMI_label_M0), "\n\n")

# Create results directory
results_dir <- file.path(base_dir, "results_pathways")
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

# -----------------------------------------------------------------------------
# 4.2 Perform clustering for M0 visit (k = 2 to 20)
# -----------------------------------------------------------------------------
cat("4.2 Performing clustering for M0 visit...\n")
cat("   Testing k values from 2 to 20\n\n")

Sev_score_detail_M0 <- data.frame()

for (n in 2:20) {
  cat("   Processing k =", n, "...\n")
  
  res_k <- dividing_data_cluster(
    df_covariates = PPMI_label_M0,
    df_transcript = PPMI_exp_M0,
    Subsets_pathways = Subsets_pathways,
    n_clust = n,
    use_random_pathways = FALSE,
    seed = NULL
  )
  
  Sev_score_detail_M0 <- rbind(Sev_score_detail_M0, res_k)
}

# Save results for M0
saveRDS(
  Sev_score_detail_M0,
  file.path(
    results_dir,
    "Sev_score_detail_M0_pathways_Reactome_pathways_filtered_50_protein_coding.rds"
  )
)
cat("\n✅ M0 results saved\n\n")

# -----------------------------------------------------------------------------
# 4.3 Perform clustering for M12, M24, and M36 visits
# -----------------------------------------------------------------------------
cat("4.3 Performing clustering for remaining visits...\n")
cat("   Following the same strategy (k = 2 to 20) for M12, M24, M36\n\n")

# M12 visit
cat("   Processing M12...\n")
PPMI_label_M12 <- PPMI_label %>% filter(visit_name == "M12")
PPMI_label_M12$visit_name <- NULL
common_ids <- intersect(rownames(PPMI_label_M12), rownames(PPMI_exp))
PPMI_label_M12 <- PPMI_label_M12[common_ids, ]
PPMI_exp_M12   <- PPMI_exp[common_ids, ]

Sev_score_detail_M12 <- data.frame()
for (n in 2:20) {
  res_k <- dividing_data_cluster(
    df_covariates = PPMI_label_M12,
    df_transcript = PPMI_exp_M12,
    Subsets_pathways = Subsets_pathways,
    n_clust = n,
    use_random_pathways = FALSE,
    seed = NULL
  )
  Sev_score_detail_M12 <- rbind(Sev_score_detail_M12, res_k)
}
saveRDS(Sev_score_detail_M12,
        file.path(results_dir, "Sev_score_detail_M12_pathways_Reactome_pathways_filtered_50_protein_coding.rds"))
cat("   ✅ M12 complete\n")

# M24 visit
cat("\n   Processing M24...\n")
PPMI_label_M24 <- PPMI_label %>% filter(visit_name == "M24")
PPMI_label_M24$visit_name <- NULL
common_ids <- intersect(rownames(PPMI_label_M24), rownames(PPMI_exp))
PPMI_label_M24 <- PPMI_label_M24[common_ids, ]
PPMI_exp_M24   <- PPMI_exp[common_ids, ]

Sev_score_detail_M24 <- data.frame()
for (n in 2:20) {
  res_k <- dividing_data_cluster(
    df_covariates = PPMI_label_M24,
    df_transcript = PPMI_exp_M24,
    Subsets_pathways = Subsets_pathways,
    n_clust = n,
    use_random_pathways = FALSE,
    seed = NULL
  )
  Sev_score_detail_M24 <- rbind(Sev_score_detail_M24, res_k)
}
saveRDS(Sev_score_detail_M24,
        file.path(results_dir, "Sev_score_detail_M24_pathways_Reactome_pathways_filtered_50_protein_coding.rds"))
cat("   ✅ M24 complete\n")

# M36 visit
cat("\n   Processing M36...\n")
PPMI_label_M36 <- PPMI_label %>% filter(visit_name == "M36")
PPMI_label_M36$visit_name <- NULL
common_ids <- intersect(rownames(PPMI_label_M36), rownames(PPMI_exp))
PPMI_label_M36 <- PPMI_label_M36[common_ids, ]
PPMI_exp_M36   <- PPMI_exp[common_ids, ]

Sev_score_detail_M36 <- data.frame()
for (n in 2:20) {
  res_k <- dividing_data_cluster(
    df_covariates = PPMI_label_M36,
    df_transcript = PPMI_exp_M36,
    Subsets_pathways = Subsets_pathways,
    n_clust = n,
    use_random_pathways = FALSE,
    seed = NULL
  )
  Sev_score_detail_M36 <- rbind(Sev_score_detail_M36, res_k)
}
saveRDS(Sev_score_detail_M36,
        file.path(results_dir, "Sev_score_detail_M36_pathways_Reactome_pathways_filtered_50_protein_coding.rds"))
cat("   ✅ M36 complete\n\n")

# =============================================================================
# 5. GENERATION OF MOTOR AND NON-MOTOR SEVERITY SCORES
# =============================================================================
cat("========================================\n")
cat("5. GENERATION OF MOTOR AND NON-MOTOR SEVERITY SCORES\n")
cat("========================================\n\n")

#' Create severity score summary for a specific visit
#'
#' Aggregates severity signals across pathways for each patient and covariate,
#' then calculates total scores.
#'
#' @param Sev_score_detail_MSigDB Data frame with severity signals
#' @param PPMI_label Clinical labels with visit information
#' @param visit_label Visit identifier (e.g., "M0", "M12")
#' @return Data frame with severity scores per patient
create_score_summary <- function(Sev_score_detail_MSigDB, PPMI_label, visit_label) {
  
  # Convert rownames to "sample_id" column
  PPMI_label <- PPMI_label %>%
    tibble::rownames_to_column(var = "sample_id")
  
  # Extract sample_ids for the specific visit
  sample_ids_visit <- PPMI_label %>%
    filter(visit_name == visit_label) %>%
    dplyr::select(sample_id, visit_name)
  
  # Sum binary scores (SevScore) by sample_id and Covariate
  score_orig <- Sev_score_detail_MSigDB %>%
    group_by(sample_id, Covariate) %>%
    summarise(score_sum = sum(SevScore, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(
      names_from = Covariate,
      values_from = score_sum,
      names_prefix = "Score_",
      values_fill = list(score_sum = 0)
    )
  
  # Add patients without scores and fill with zeros
  full_score_summary <- sample_ids_visit %>%
    left_join(score_orig, by = "sample_id") %>%
    mutate(across(where(is.numeric), ~replace_na(., 0)))
  
  return(full_score_summary)
}

# -----------------------------------------------------------------------------
# 5.1 Generate score summaries for each visit
# -----------------------------------------------------------------------------
cat("5.1 Generating score summaries for each visit...\n")

Severity_score_summ_M0 <- create_score_summary(Sev_score_detail_M0, PPMI_label, visit_label = "M0") %>%
  mutate(Visit = 0,
         participant_id = as.numeric(sub("PP-(\\d+)-.*", "\\1", sample_id)))

Severity_score_summ_M12 <- create_score_summary(Sev_score_detail_M12, PPMI_label, visit_label = "M12") %>%
  mutate(Visit = 12,
         participant_id = as.numeric(sub("PP-(\\d+)-.*", "\\1", sample_id)))

Severity_score_summ_M24 <- create_score_summary(Sev_score_detail_M24, PPMI_label, visit_label = "M24") %>%
  mutate(Visit = 24,
         participant_id = as.numeric(sub("PP-(\\d+)-.*", "\\1", sample_id)))

Severity_score_summ_M36 <- create_score_summary(Sev_score_detail_M36, PPMI_label, visit_label = "M36") %>%
  mutate(Visit = 36,
         participant_id = as.numeric(sub("PP-(\\d+)-.*", "\\1", sample_id)))

# -----------------------------------------------------------------------------
# 5.2 Combine all visits
# -----------------------------------------------------------------------------
cat("5.2 Combining scores across all visits...\n")

list_of_scores <- list(Severity_score_summ_M0, 
                       Severity_score_summ_M12, 
                       Severity_score_summ_M24, 
                       Severity_score_summ_M36)

# Combine dataframes with full join
combined_severity_score <- purrr::reduce(list_of_scores, full_join, by = NULL)

# Replace NA with 0 for all scores
combined_severity_score[is.na(combined_severity_score)] <- 0

cat("  - Combined data dimensions:", dim(combined_severity_score), "\n")

# -----------------------------------------------------------------------------
# 5.3 Define motor and non-motor domains
# -----------------------------------------------------------------------------
cat("5.3 Defining motor and non-motor domains...\n")

# Motor domain covariates (based on established clinical frameworks)
Motor <- c(
  "code_upd2hy_hoehn_and_yahr_stage", 
  "mds_updrs_part_ii_summary_score",
  "mds_updrs_part_iii_summary_score", 
  "Tremor_Right_UM", 
  "Tremor_Left_UM",
  "Tremor_all_Um", 
  "PIGD_UM",
  "mod_schwab_england_pct_adl_score"
)

# Non-motor domain covariates
Non_Motor <- c(
  "code_upd2106_dopamine_dysregulation_syndrome_features",
  "mds_updrs_part_i_summary_score", 
  "code_upd2102_hallucinations_and_psychosis",
  "code_upd2101_cognitive_impairment", 
  "code_upd2103_depressed_mood",
  "code_upd2104_anxious_mood", 
  "code_upd2105_apathy",
  "rbd_summary_score", 
  "ess_summary_score",
  "moca_total_score"
)

# Column names for scores
Motor_raw <- paste0("Score_", Motor)
Non_Motor_raw <- paste0("Score_", Non_Motor)

cat("  - Motor covariates:", length(Motor), "\n")
cat("  - Non-motor covariates:", length(Non_Motor), "\n")

# -----------------------------------------------------------------------------
# 5.4 Calculate motor and non-motor severity scores
# -----------------------------------------------------------------------------
cat("5.4 Calculating motor and non-motor severity scores...\n")

# Calculate UTSS for motor and non-motor domains
combined_severity_score_all_visits <- combined_severity_score %>%
  mutate(
    # Unsupervised Transcriptomic Severity Score (UTSS)
    UTSS_Motor = rowSums(dplyr::select(., any_of(Motor_raw)), na.rm = TRUE),
    UTSS_NonMotor = rowSums(dplyr::select(., any_of(Non_Motor_raw)), na.rm = TRUE)
  )

# Generate cumulative scores to track progression
combined_severity_score_all_visits_acc <- combined_severity_score_all_visits %>%
  arrange(participant_id, Visit) %>%
  group_by(participant_id) %>%
  mutate(across(
    c(UTSS_Motor, UTSS_NonMotor),
    cumsum,
    .names = "{.col}_Cumulative"
  )) %>%
  ungroup()

cat("✅ UTSS scores calculated for motor and non-motor domains\n\n")

# =============================================================================
# 6. INTEGRATION WITH ADDITIONAL CLINICAL COVARIATES
# =============================================================================
cat("========================================\n")
cat("6. INTEGRATION WITH ADDITIONAL CLINICAL COVARIATES\n")
cat("========================================\n\n")

# -----------------------------------------------------------------------------
# 6.1 Load and prepare additional covariates
# -----------------------------------------------------------------------------
cat("6.1 Loading additional clinical covariates...\n")

Covariates <- as.data.frame(readRDS(file.path(base_dir, "metadata_reduced", "PPMI_covariates_extra")))
rownames(Covariates) <- Covariates$sample_id

# Add additional covariates
cna <- readRDS(file.path(base_dir, "Pipeline", "Covariates_new_all"))

Covariates <- Covariates %>%
  left_join(cna %>% dplyr::select(sample_id, AGE_ONSET, prognosis_status, 
                                  apoe, mapt, dementia_status, MOCA_Status), 
            by = "sample_id")

# Propagate time-invariant variables within participants
Covariates <- Covariates %>%
  group_by(participant_id) %>%
  fill(AGE_ONSET, prognosis_status, apoe, mapt, .direction = "downup") %>%
  ungroup()

# -----------------------------------------------------------------------------
# 6.2 Filter covariates for patients with UTSS scores
# -----------------------------------------------------------------------------
cat("6.2 Filtering covariates for patients with UTSS scores...\n")

common_patients <- unique(combined_severity_score_all_visits$participant_id)
common_patients_id <- paste("PP-", common_patients, sep = "")

Covariates_filtered <- Covariates %>% 
  filter(participant_id %in% common_patients_id & visit_name != "M6") %>%
  mutate(Visit = as.numeric(str_extract(visit_name, "\\d+")))

cat("  - Patients with UTSS scores:", length(common_patients), "\n")
cat("  - Covariates records:", nrow(Covariates_filtered), "\n")

# -----------------------------------------------------------------------------
# 6.3 Merge UTSS scores with clinical covariates
# -----------------------------------------------------------------------------
cat("6.3 Merging UTSS scores with clinical covariates...\n")

score_columns <- c(
  "UTSS_Motor", "UTSS_NonMotor",
  "UTSS_Motor_Cumulative", "UTSS_NonMotor_Cumulative"
)

Covariates_score <- Covariates_filtered %>%
  left_join(combined_severity_score_all_visits_acc %>% 
              dplyr::select("sample_id", all_of(score_columns)), 
            by = "sample_id")

# Display summary
cat("\nFinal dataset summary:\n")
print(summary(Covariates_score[, c("participant_id", "visit_name", score_columns)]))

# Save final dataset
saveRDS(Covariates_score, file.path(results_dir, "Covariates_with_UTSS_scores.rds"))
cat("\n✅ Final dataset saved with UTSS motor and non-motor scores\n")

# =============================================================================
# 7. SESSION INFORMATION
# =============================================================================
cat("\n========================================\n")
cat("7. SESSION INFORMATION\n")
cat("========================================\n\n")

print(sessionInfo())

cat("\n========================================\n")
cat("Clustering and severity score generation completed successfully!\n")
cat("========================================\n")