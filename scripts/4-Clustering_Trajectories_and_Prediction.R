#!/usr/bin/env Rscript

# =============================================================================
# LONGITUDINAL DISEASE PROGRESSION TRAJECTORY ANALYSIS AND PREDICTION
# PPMI Cohort - Unsupervised Transcriptomic Severity Score (UTSS)
# =============================================================================
#
# Description: This script performs longitudinal analysis of disease progression
#              trajectories based on UTSS scores and builds predictive models
#              using baseline transcriptomic data.
#
#              Key analyses include:
#              1. Definition of progression trajectories using k-means clustering
#                 on baseline UTSS and longitudinal changes between visits
#              2. Prediction of progression trajectories from baseline gene
#                 expression using Random Forest models
#              3. Identification of progression-associated transcriptomic features
#                 through variable importance analysis and differential expression
#
# Author: [Jaime Ñíguez Baeza]
# Date: 2026
# Corresponding author: [juanbot@um.es]
#
# Input files:
#   - Covariates_with_UTSS_scores.rds: Combined dataset with UTSS scores
#   - PPMI_exp_M0: Baseline expression matrix
#   - Subsets_pathways_Reactome: Reactome pathway gene sets
#   - Sev_score_detail_M*: Severity signals for gene selection
#   - PPMI_raw_counts: Raw count data for differential expression
#
# Output files:
#   - progression_prediction_*.rds: Model results for each k value
#   - *_VarImp_elbow_results.csv: Variable importance rankings
#   - *_DEG_results.csv: Differential expression results
#   - PD_genes_*.csv: Parkinson's disease genes identified
#   - Variable_Importance_Results_elbow.RData: Complete importance results
#
# Dependencies: Various (see library calls)
# =============================================================================

# =============================================================================
# 1. INITIAL SETUP AND DATA LOADING
# =============================================================================

# Load required libraries
library(tidyverse)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(readxl)
library(readr)

# Machine learning
library(caret)
library(ranger)
library(randomForest)

# Differential expression
library(DESeq2)

# Gene annotation
library(biomaRt)

# Set seed for reproducibility
set.seed(202401201)

# Define base directory (modify as needed)
base_dir <- "/home/jaimeniguez/Pipeline"

cat("========================================\n")
cat("LONGITUDINAL DISEASE PROGRESSION TRAJECTORY ANALYSIS\n")
cat("========================================\n\n")

# -----------------------------------------------------------------------------
# 1.1 Load severity scores and expression data
# -----------------------------------------------------------------------------
cat("1.1 Loading severity scores and expression data...\n")

# Load combined severity scores with UTSS
severity_scores <- readRDS(file.path(base_dir, "results_pathways", "Covariates_with_UTSS_scores.rds"))
cat("  - Severity scores loaded:", nrow(severity_scores), "observations\n")

# Load baseline expression data
PPMI_exp_M0 <- readRDS(file.path(base_dir, "PPMI_exp_M0.rds"))  # Adjust path as needed
cat("  - Baseline expression loaded:", dim(PPMI_exp_M0)[1], "samples ×", 
    dim(PPMI_exp_M0)[2], "genes\n")

# Load pathways and signal genes
Subsets_pathways <- readRDS(file.path(base_dir, "Subsets_pathways_Reactome"))

# Get genes from pathways that generated signals
Sev_score_detail_M0 <- readRDS(file.path(results_dir, "Sev_score_detail_M0_pathways_Reactome_pathways_filtered_50_protein_coding.rds"))
Sev_score_detail_M12 <- readRDS(file.path(results_dir, "Sev_score_detail_M12_pathways_Reactome_pathways_filtered_50_protein_coding.rds"))
Sev_score_detail_M24 <- readRDS(file.path(results_dir, "Sev_score_detail_M24_pathways_Reactome_pathways_filtered_50_protein_coding.rds"))
Sev_score_detail_M36 <- readRDS(file.path(results_dir, "Sev_score_detail_M36_pathways_Reactome_pathways_filtered_50_protein_coding.rds"))

pathways_signal <- unique(c(
  Sev_score_detail_M0$Pathway, 
  Sev_score_detail_M12$Pathway, 
  Sev_score_detail_M24$Pathway, 
  Sev_score_detail_M36$Pathway
))

pathways_gene_signal <- Subsets_pathways[names(Subsets_pathways) %in% pathways_signal]
genes_signal <- unique(unlist(pathways_gene_signal))
genes_present <- genes_signal[genes_signal %in% colnames(PPMI_exp_M0)]

cat("  - Signal-generating pathways:", length(pathways_signal), "\n")
cat("  - Unique genes in these pathways:", length(genes_signal), "\n")
cat("  - Genes present in expression data:", length(genes_present), "\n\n")

# =============================================================================
# 2. DISEASE PROGRESSION TRAJECTORY ANALYSIS
# =============================================================================
cat("========================================\n")
cat("2. DISEASE PROGRESSION TRAJECTORY ANALYSIS\n")
cat("========================================\n\n")

#' Evaluate clustering of progression trajectories for a specific k value
#'
#' This function:
#' 1. Transforms UTSS values to percentile ranks at each visit
#' 2. Creates feature vectors combining baseline and longitudinal changes
#' 3. Performs k-means clustering to identify progression trajectories
#' 4. Trains Random Forest models to predict trajectories from baseline expression
#'
#' @param severity_scores Data frame with UTSS scores per patient and visit
#' @param expression_data Baseline expression matrix
#' @param score_type Domain to analyze ("Motor" or "Non_Motor")
#' @param k Number of clusters for trajectory definition
#' @param test_size Proportion of data for testing
#' @param n_folds Number of cross-validation folds
#' @param mtry_select Values of mtry to test in tuning
#' @param ntrees Number of trees in Random Forest
#' @param seed Random seed for reproducibility
#' @return List with clustering results, model, and performance metrics
evaluate_clusters_simple <- function(
    severity_scores,
    expression_data,
    score_type = c("Motor", "Non_Motor"),
    k = 3,
    test_size = 0.2,
    n_folds = 5,
    mtry_select = c(10, 30, 60),
    ntrees = 500,
    seed = 202401201
) {
  
  library(caret)
  library(ranger)
  library(dplyr)
  library(tidyr)
  
  score_type <- match.arg(score_type)
  
  # Select appropriate score column based on domain
  score_column <- switch(score_type,
                         "Motor" = "UTSS_Motor",
                         "Non_Motor" = "UTSS_NonMotor"
  )
  
  # ------------------------------------------------------------------
  # 2.1 TRANSFORM TO PERCENTILE RANKS AT EACH VISIT
  # ------------------------------------------------------------------
  cat("    Transforming to percentile ranks...\n")
  
  # Calculate percentile ranks for each visit separately
  severity_scores <- severity_scores %>%
    group_by(Visit) %>%
    mutate(
      !!paste0(score_column, "_percentile") := percent_rank(!!sym(score_column)) * 100
    ) %>%
    ungroup()
  
  # ------------------------------------------------------------------
  # 2.2 CREATE LONGITUDINAL FEATURE VECTORS
  # ------------------------------------------------------------------
  cat("    Creating longitudinal feature vectors...\n")
  
  percentile_column <- paste0(score_column, "_percentile")
  
  df_long <- severity_scores[, c("participant_id", "Visit", percentile_column)]
  df_long <- df_long[order(df_long$participant_id, df_long$Visit), ]
  
  # Pivot to wide format
  df_wide <- pivot_wider(
    df_long,
    names_from = Visit,
    values_from = all_of(percentile_column),
    names_prefix = "Visit_"
  )
  
  # Create feature vector: baseline + changes between visits
  # {S0, (S12-S0), (S24-S12), (S36-S24)}
  df_diff <- df_wide %>%
    mutate(
      S0 = Visit_0,
      diff_0_12 = Visit_12 - Visit_0,
      diff_12_24 = Visit_24 - Visit_12,
      diff_24_36 = Visit_36 - Visit_24
    ) %>%
    dplyr::select(participant_id, S0, diff_0_12, diff_12_24, diff_24_36)
  
  # Scale features for clustering
  data_scaled <- scale(df_diff[, -1])
  
  # ------------------------------------------------------------------
  # 2.3 K-MEANS CLUSTERING TO IDENTIFY TRAJECTORIES
  # ------------------------------------------------------------------
  cat("    Performing k-means clustering with k =", k, "...\n")
  
  set.seed(seed)
  km <- kmeans(data_scaled, centers = k, nstart = 25)
  
  cluster_df <- data.frame(
    participant_id = df_diff$participant_id,
    cluster = factor(make.names(km$cluster))
  )
  
  # ------------------------------------------------------------------
  # 2.4 MERGE WITH BASELINE EXPRESSION DATA
  # ------------------------------------------------------------------
  cat("    Merging with baseline expression data...\n")
  
  # Add participant_id to expression data
  PPMI_exp_M0_id <- expression_data %>%
    as.data.frame() %>%
    mutate(
      sample_id = rownames(.),
      participant_id = as.numeric(sub("PP-(\\d+)-.*", "\\1", sample_id))
    ) %>%
    dplyr::select(-sample_id)
  
  df_initial <- PPMI_exp_M0_id %>%
    inner_join(cluster_df, by = "participant_id") %>%
    dplyr::select(-participant_id)
  
  df_initial$cluster <- factor(df_initial$cluster)
  
  # ------------------------------------------------------------------
  # 2.5 TRAIN/TEST SPLIT (STRATIFIED BY TRAJECTORY)
  # ------------------------------------------------------------------
  cat("    Splitting data into training/test sets...\n")
  
  set.seed(seed)
  train_idx <- createDataPartition(df_initial$cluster, p = 1 - test_size, list = FALSE)
  
  train_data <- df_initial[train_idx, ]
  test_data <- df_initial[-train_idx, ]
  
  # ------------------------------------------------------------------
  # 2.6 REMOVE HIGHLY CORRELATED FEATURES (TRAINING DATA ONLY)
  # ------------------------------------------------------------------
  cat("    Removing highly correlated genes (r > 0.95)...\n")
  
  remove_correlated_genes <- function(data, threshold = 0.95) {
    gene_data <- data %>% dplyr::select(-cluster)
    
    if (ncol(gene_data) > 1) {
      cor_mat <- cor(gene_data, use = "complete.obs")
      high_cor <- which(abs(cor_mat) > threshold & upper.tri(cor_mat), arr.ind = TRUE)
      
      if (length(high_cor) > 0) {
        genes_to_remove <- unique(colnames(gene_data)[high_cor[, 2]])
        genes_keep <- setdiff(colnames(gene_data), genes_to_remove)
        result <- data %>% dplyr::select(all_of(c(genes_keep, "cluster")))
        return(result)
      }
    }
    return(data)
  }
  
  train_filtered <- remove_correlated_genes(train_data, 0.95)
  genes_keep <- setdiff(colnames(train_filtered), "cluster")
  test_filtered <- test_data %>% dplyr::select(all_of(c(genes_keep, "cluster")))
  
  cat("    Genes retained after correlation filtering:", length(genes_keep), "\n")
  
  # ------------------------------------------------------------------
  # 2.7 HYPERPARAMETER TUNING WITH 5-FOLD CROSS-VALIDATION
  # ------------------------------------------------------------------
  cat("    Tuning hyperparameters with 5-fold CV...\n")
  
  tune_grid <- expand.grid(
    mtry = mtry_select,
    splitrule = "gini",
    min.node.size = 10
  )
  
  train_control <- trainControl(
    method = "cv",
    number = n_folds,
    classProbs = TRUE,
    summaryFunction = multiClassSummary,
    savePredictions = "final",
    allowParallel = FALSE
  )
  
  set.seed(seed)
  cv_model <- train(
    x = train_filtered %>% dplyr::select(-cluster),
    y = train_filtered$cluster,
    method = "ranger",
    tuneGrid = tune_grid,
    trControl = train_control,
    num.trees = ntrees,
    importance = "permutation",
    metric = "Accuracy"
  )
  
  best_params <- cv_model$bestTune
  cat("    Best mtry:", best_params$mtry, "\n")
  
  # ------------------------------------------------------------------
  # 2.8 FINAL MODEL TRAINING
  # ------------------------------------------------------------------
  cat("    Training final Random Forest model...\n")
  
  final_model <- ranger(
    cluster ~ .,
    data = train_filtered,
    mtry = best_params$mtry,
    splitrule = best_params$splitrule,
    min.node.size = best_params$min.node.size,
    num.trees = ntrees,
    importance = "permutation",
    seed = seed
  )
  
  # ------------------------------------------------------------------
  # 2.9 TEST SET EVALUATION
  # ------------------------------------------------------------------
  cat("    Evaluating on test set...\n")
  
  test_pred <- predict(final_model, test_filtered)$predictions
  test_cm <- confusionMatrix(test_pred, test_filtered$cluster)
  
  # Extract metrics
  metrics <- list(
    OverallAccuracy = test_cm$overall["Accuracy"],
    BalancedAccuracy = if(is.matrix(test_cm$byClass)) {
      mean(test_cm$byClass[, "Balanced Accuracy"], na.rm = TRUE)
    } else {
      test_cm$byClass["Balanced Accuracy"]
    },
    Kappa = test_cm$overall["Kappa"],
    PValue = test_cm$overall["AccuracyPValue"],
    CV_Accuracy = cv_model$results %>%
      filter(mtry == best_params$mtry) %>%
      pull(Accuracy)
  )
  
  cat("    Test accuracy:", round(metrics$OverallAccuracy, 3), "\n")
  cat("    P-value vs No Information Rate:", round(metrics$PValue, 4), "\n")
  
  # ------------------------------------------------------------------
  # 2.10 RETURN RESULTS
  # ------------------------------------------------------------------
  list(
    ScoreType = score_type,
    K = k,
    BestParams = best_params,
    FinalModel = final_model,
    TestConfusionMatrix = test_cm,
    Metrics = metrics,
    VariableImportance = data.frame(
      Gene = names(final_model$variable.importance),
      Importance = final_model$variable.importance
    ) %>% arrange(desc(Importance)),
    TrainingSize = nrow(train_filtered),
    TestSize = nrow(test_filtered),
    GenesUsed = length(genes_keep),
    ClusteringResults = list(
      kmeans_obj = km,
      cluster_assignments = cluster_df,
      feature_matrix = df_diff
    )
  )
}

#' Run clustering and prediction for a range of k values
#'
#' @param severity_scores Data frame with UTSS scores
#' @param expression_data Baseline expression matrix
#' @param score_type Domain to analyze ("Motor" or "Non_Motor")
#' @param k_range Range of k values to test (typically 2-10)
#' @param n_folds Number of cross-validation folds
#' @param save_prefix Prefix for output files
#' @return List of results for each k value
run_all_clusters_simple <- function(
    severity_scores,
    expression_data,
    score_type = c("Motor", "Non_Motor"),
    k_range = 2:10,
    n_folds = 5,
    save_prefix = "progression_prediction"
) {
  
  score_type <- match.arg(score_type)
  results <- list()
  
  cat("\n", rep("-", 50), "\n", sep = "")
  cat("Running", score_type, "domain analysis\n")
  cat(rep("-", 50), "\n", sep = "")
  
  for (k in k_range) {
    cat("\n>>> Testing k =", k, "\n")
    
    result <- evaluate_clusters_simple(
      severity_scores = severity_scores,
      expression_data = expression_data,
      score_type = score_type,
      k = k,
      n_folds = n_folds
    )
    
    # Save individual results
    filename <- paste0(save_prefix, "_", score_type, "_k", k, ".rds")
    saveRDS(result, filename)
    
    # Print summary
    cat("\n  Summary for k =", k, ":\n")
    cat("  Accuracy:", round(result$Metrics$OverallAccuracy, 3), 
        "(Balanced:", round(result$Metrics$BalancedAccuracy, 3), ")\n")
    cat("  Kappa:", round(result$Metrics$Kappa, 3), 
        "P-value:", round(result$Metrics$PValue, 4), "\n")
    cat("  ---\n")
    
    results[[paste0("k", k)]] <- result
  }
  
  return(results)
}

# =============================================================================
# 3. EXECUTE TRAJECTORY ANALYSIS FOR BOTH DOMAINS
# =============================================================================
cat("\n========================================\n")
cat("3. EXECUTING TRAJECTORY ANALYSIS\n")
cat("========================================\n\n")

# Define prefix for saving results
prefix <- "progression_prediction"

# Run analysis for Non-Motor domain
cat("\n>>> NON-MOTOR DOMAIN TRAJECTORIES <<<\n")
non_motor_results <- run_all_clusters_simple(
  severity_scores = severity_scores,
  expression_data = PPMI_exp_M0[, genes_present],
  score_type = "Non_Motor",
  k_range = 2:10,
  n_folds = 5,
  save_prefix = prefix
)

# Run analysis for Motor domain
cat("\n>>> MOTOR DOMAIN TRAJECTORIES <<<\n")
motor_results <- run_all_clusters_simple(
  severity_scores = severity_scores,
  expression_data = PPMI_exp_M0[, genes_present],
  score_type = "Motor",
  k_range = 2:10,
  n_folds = 5,
  save_prefix = prefix
)

# Identify best k based on test accuracy
motor_performance <- data.frame(
  k = 2:10,
  Accuracy = sapply(motor_results, function(x) x$Metrics$OverallAccuracy),
  BalancedAcc = sapply(motor_results, function(x) x$Metrics$BalancedAccuracy),
  PValue = sapply(motor_results, function(x) x$Metrics$PValue)
)

non_motor_performance <- data.frame(
  k = 2:10,
  Accuracy = sapply(non_motor_results, function(x) x$Metrics$OverallAccuracy),
  BalancedAcc = sapply(non_motor_results, function(x) x$Metrics$BalancedAccuracy),
  PValue = sapply(non_motor_results, function(x) x$Metrics$PValue)
)

motor_best_k <- motor_performance[which.max(motor_performance$Accuracy), ]
non_motor_best_k <- non_motor_performance[which.max(non_motor_performance$Accuracy), ]

cat("\n=== OPTIMAL CLUSTER SOLUTIONS ===\n")
cat("Motor domain - Best k =", motor_best_k$k, 
    "(Accuracy:", round(motor_best_k$Accuracy, 3), 
    "P-value:", round(motor_best_k$PValue, 4), ")\n")
cat("Non-Motor domain - Best k =", non_motor_best_k$k, 
    "(Accuracy:", round(non_motor_best_k$Accuracy, 3), 
    "P-value:", round(non_motor_best_k$PValue, 4), ")\n\n")

# Store best models for downstream analysis
best_motor_model <- motor_results[[paste0("k", motor_best_k$k)]]
best_non_motor_model <- non_motor_results[[paste0("k", non_motor_best_k$k)]]

# =============================================================================
# 4. IDENTIFICATION OF PROGRESSION-ASSOCIATED TRANSCRIPTOMIC FEATURES
# =============================================================================
cat("========================================\n")
cat("4. IDENTIFICATION OF PROGRESSION-ASSOCIATED FEATURES\n")
cat("========================================\n\n")

# -----------------------------------------------------------------------------
# 4.1 Load Parkinson's disease risk gene lists
# -----------------------------------------------------------------------------
cat("4.1 Loading Parkinson's disease risk gene lists...\n")

# 1. Processed PD gene list
processed_PD_gene_list <- readxl::read_excel("./genes_risk_pd/processed_PD_gene_list.xls")

# 2. Nalls et al. 2019 (PMID: 31701892)
MikeNalls_PMID31701892_associations_export <- readxl::read_excel("./genes_risk_pd/MikeNalls_PMID31701892_associations_export.xls")

# 3. GWAS from Kim et al. 2023 (PMID: 38155330)
gwas_2023_Kim <- read_delim("~/Pipeline/gwas-association-downloaded_2026-02-05-pubmedId_38155330.tsv", delim = "\t")

# Extract gene symbols
processed_PD_gene <- processed_PD_gene_list$symbol
MikeNalls_genes <- MikeNalls_PMID31701892_associations_export$symbol
gwas_2023_Kim_genes <- gwas_2023_Kim$MAPPED_GENE

# Clean and process symbols
clean_gene_symbols <- function(gene_vector) {
  # Split multiple genes (separated by comma, semicolon, etc.)
  genes <- unlist(strsplit(as.character(gene_vector), "[,;|]"))
  # Remove spaces and convert to uppercase
  genes <- trimws(toupper(genes))
  # Remove empty entries
  genes <- genes[genes != "" & !is.na(genes)]
  # Remove duplicates
  unique(genes)
}

processed_PD_gene <- clean_gene_symbols(processed_PD_gene)
MikeNalls_genes <- clean_gene_symbols(MikeNalls_genes)
gwas_2023_Kim_genes <- clean_gene_symbols(gwas_2023_Kim_genes)

cat("  - Processed PD genes:", length(processed_PD_gene), "genes\n")
cat("  - Nalls et al. 2019:", length(MikeNalls_genes), "genes\n")
cat("  - Kim et al. 2023:", length(gwas_2023_Kim_genes), "genes\n\n")

# Load Ensembl mart for gene annotation
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# -----------------------------------------------------------------------------
# 4.2 Function to identify PD genes in importance lists
# -----------------------------------------------------------------------------
#' Identify Parkinson's disease genes in variable importance lists
#'
#' @param important_genes Vector of important genes (Ensembl IDs)
#' @param domain_name Domain being analyzed ("Motor" or "Non_Motor")
#' @param pd_gene_lists List of PD gene sets
#' @param min_intersect Minimum number of lists where gene must appear
#' @param return_all Return detailed information
#' @return Vector of highlighted genes or detailed list
get_pd_highlight_genes <- function(
    important_genes,
    domain_name,
    pd_gene_lists = list(
      processed_PD = processed_PD_gene,
      MikeNalls = MikeNalls_genes,
      Kim_2023 = gwas_2023_Kim_genes
    ),
    min_intersect = 1,
    return_all = FALSE
) {
  
  # Convert Ensembl IDs to gene symbols if needed
  if (any(grepl("^ENSG", important_genes[1:10]))) {
    cat("    Converting Ensembl IDs to gene symbols...\n")
    gene_map <- getBM(
      attributes = c("ensembl_gene_id", "hgnc_symbol"),
      filters = "ensembl_gene_id",
      values = important_genes,
      mart = mart
    )
    
    important_symbols <- sapply(important_genes, function(id) {
      symbol <- gene_map$hgnc_symbol[match(id, gene_map$ensembl_gene_id)]
      ifelse(!is.na(symbol) & symbol != "", symbol, id)
    })
  } else {
    important_symbols <- important_genes
  }
  
  # Convert to uppercase for comparison
  important_symbols <- toupper(important_symbols)
  
  # Calculate intersections
  intersections <- list()
  highlight_counts <- list()
  
  for (list_name in names(pd_gene_lists)) {
    pd_list <- toupper(pd_gene_lists[[list_name]])
    intersect_genes <- intersect(important_symbols, pd_list)
    intersections[[list_name]] <- intersect_genes
    highlight_counts[[list_name]] <- length(intersect_genes)
  }
  
  # Find genes appearing in multiple lists
  all_intersected <- unlist(intersections)
  gene_freq <- table(all_intersected)
  
  # Select genes based on frequency
  if (min_intersect > 1) {
    highlight_genes <- names(gene_freq[gene_freq >= min_intersect])
  } else {
    highlight_genes <- unique(all_intersected)
  }
  
  # Print summary
  cat("\n    PD GENES FOUND IN", domain_name, ":\n")
  for (list_name in names(pd_gene_lists)) {
    cat(sprintf("    %-20s: %4d genes\n", list_name, highlight_counts[[list_name]]))
  }
  
  cat("    Genes in ≥", min_intersect, "list(s):", length(highlight_genes), "\n")
  
  if (return_all) {
    return(list(
      highlight_genes = highlight_genes,
      intersections = intersections,
      counts = highlight_counts,
      gene_frequencies = gene_freq
    ))
  } else {
    return(highlight_genes)
  }
}

# -----------------------------------------------------------------------------
# 4.3 Function for variable importance analysis with elbow selection
# -----------------------------------------------------------------------------
#' Analyze variable importance from Random Forest models
#'
#' Uses curvature-based elbow method to select optimal number of genes
#'
#' @param results_object Model results from evaluate_clusters_simple
#' @param domain_name Domain being analyzed ("Motor" or "Non_Motor")
#' @param additional_vars Additional variables to show beyond elbow
#' @param highlight_genes Manually specified genes to highlight
#' @param pd_gene_lists List of PD gene sets
#' @param pd_min_freq Minimum frequency for PD gene highlighting
#' @param save_plot Save importance plot to file
#' @param output_dir Output directory for plots
#' @param file_prefix Prefix for output files
#' @return List with importance results
analyze_variable_importance <- function(
    results_object,
    domain_name = "Motor",
    additional_vars = 1000,
    highlight_genes = NULL,
    pd_gene_lists = list(
      processed_PD = processed_PD_gene,
      MikeNalls = MikeNalls_genes,
      Kim_2023 = gwas_2023_Kim_genes
    ),
    pd_min_freq = 1,
    save_plot = FALSE,
    output_dir = ".",
    file_prefix = "VarImp"
) {
  
  cat("\n", rep("-", 50), "\n", sep = "")
  cat("Analyzing variable importance for", domain_name, "domain\n")
  cat(rep("-", 50), "\n", sep = "")
  
  # ===============================
  # Extract variable importance
  # ===============================
  var_imp <- results_object$VariableImportance %>%
    as.data.frame() %>%
    filter(Importance > 0) %>%  # Only genes with positive importance
    arrange(desc(Importance)) %>%
    mutate(
      Importance_raw = Importance,
      Importance_plot = 100 * Importance / sum(Importance)
    )
  
  cat("\nGenes with positive importance:", nrow(var_imp), "/", 
      nrow(results_object$VariableImportance), 
      "(", round(nrow(var_imp)/nrow(results_object$VariableImportance)*100, 1), "%)\n")
  
  # ===============================
  # ELBOW DETECTION (MAXIMUM CURVATURE METHOD)
  # ===============================
  calculate_elbow_point_simple <- function(var_imp) {
    
    importances <- var_imp$Importance
    importances <- importances / sum(importances)
    cumulative_percentage <- cumsum(importances)
    
    n_vars <- length(importances)
    
    coords <- cbind(1:n_vars, cumulative_percentage)
    start <- coords[1, ]
    end <- coords[n_vars, ]
    
    # Calculate distance from each point to the line connecting start and end
    distances <- apply(coords, 1, function(p) {
      abs((end[2] - start[2]) * p[1] -
            (end[1] - start[1]) * p[2] +
            end[1] * start[2] -
            end[2] * start[1]) /
        sqrt((end[2] - start[2])^2 + (end[1] - start[1])^2)
    })
    
    # Maximum curvature point (traditional elbow)
    elbow_point <- which.max(distances)
    
    cat("\n=== SELECTION METHOD ===\n")
    cat("Method: Maximum curvature (traditional elbow)\n")
    cat("Selected cutoff point:", elbow_point, "genes\n")
    cat("Cumulative importance:", round(cumulative_percentage[elbow_point] * 100, 1), "%\n")
    
    list(
      elbow_index = elbow_point,
      cumulative_importance = cumulative_percentage[elbow_point]
    )
  }
  
  elbow <- calculate_elbow_point_simple(var_imp)
  
  important_vars <- var_imp[1:elbow$elbow_index, ]
  n_show <- min(elbow$elbow_index + additional_vars, nrow(var_imp))
  var_imp_filtered <- var_imp[1:n_show, ]
  
  # ===============================
  # IDENTIFY PD GENES AUTOMATICALLY
  # ===============================
  pd_genes_info <- get_pd_highlight_genes(
    important_genes = important_vars$Gene,
    domain_name = domain_name,
    pd_gene_lists = pd_gene_lists,
    min_intersect = pd_min_freq,
    return_all = TRUE
  )
  
  # Combine PD genes with manually specified genes
  all_highlight_genes <- unique(c(highlight_genes, pd_genes_info$highlight_genes))
  
  # ===============================
  # CONVERT ENSEMBL IDs TO SYMBOLS FOR PLOTTING
  # ===============================
  var_imp_filtered$Gene_Symbol <- var_imp_filtered$Gene
  
  # Map Ensembl IDs to symbols
  if (length(all_highlight_genes) > 0 || any(grepl("^ENSG", var_imp_filtered$Gene[1:10]))) {
    
    gene_map <- getBM(
      attributes = c("ensembl_gene_id", "hgnc_symbol"),
      filters = "ensembl_gene_id",
      values = var_imp_filtered$Gene,
      mart = mart
    )
    
    var_imp_filtered$Gene_Symbol <- gene_map$hgnc_symbol[
      match(var_imp_filtered$Gene, gene_map$ensembl_gene_id)
    ]
    
    var_imp_filtered$Gene_Symbol[
      is.na(var_imp_filtered$Gene_Symbol) | var_imp_filtered$Gene_Symbol == ""
    ] <- var_imp_filtered$Gene[is.na(var_imp_filtered$Gene_Symbol) | 
                                 var_imp_filtered$Gene_Symbol == ""]
    
    # Mark genes to highlight
    var_imp_filtered$Highlight <- 
      toupper(var_imp_filtered$Gene_Symbol) %in% toupper(all_highlight_genes)
    
    # Identify PD vs manual genes
    var_imp_filtered$Highlight_Type <- ifelse(
      toupper(var_imp_filtered$Gene_Symbol) %in% toupper(pd_genes_info$highlight_genes),
      "PD Gene",
      ifelse(toupper(var_imp_filtered$Gene_Symbol) %in% toupper(highlight_genes),
             "Manual", "None")
    )
    
  } else {
    var_imp_filtered$Highlight <- FALSE
    var_imp_filtered$Highlight_Type <- "None"
  }
  
  # ===============================
  # CREATE IMPORTANCE PLOT
  # ===============================
  y_cut <- nrow(var_imp_filtered) - elbow$elbow_index + 0.5
  
  pd_highlight_count <- sum(var_imp_filtered$Highlight_Type == "PD Gene")
  manual_highlight_count <- sum(var_imp_filtered$Highlight_Type == "Manual")
  
  p <- ggplot(
    var_imp_filtered,
    aes(
      x = Importance_plot,
      y = reorder(Gene_Symbol, Importance_raw),
      fill = Highlight_Type
    )
  ) +
    geom_col(width = 0.7, alpha = 0.85) +
    scale_fill_manual(
      values = c("None" = "steelblue", "PD Gene" = "red", "Manual" = "orange"),
      guide = "none"
    ) +
    geom_hline(
      yintercept = y_cut,
      color = "red",
      linewidth = 1.3,
      alpha = 0.9
    ) +
    geom_text(
      data = subset(var_imp_filtered, Highlight_Type != "None"),
      aes(label = Gene_Symbol, color = Highlight_Type),
      hjust = -0.1,
      size = 2,
      fontface = "bold"
    ) +
    scale_color_manual(
      values = c("PD Gene" = "red", "Manual" = "orange"),
      guide = "none"
    ) +
    labs(
      title = paste("Variable Importance –", domain_name, "Domain"),
      subtitle = paste0(
        "Elbow-selected variables: ", elbow$elbow_index,
        " (", round(elbow$cumulative_importance * 100, 1), "% cumulative importance)"
      ),
      x = "Relative Importance (%)",
      y = NULL
    ) +
    theme_minimal() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5)
    ) +
    coord_cartesian(clip = "off")
  
  # Save plot if requested
  if (save_plot) {
    file_name <- paste0(file_prefix, "_", domain_name, "_VarImp.png")
    ggsave(
      filename = file.path(output_dir, file_name),
      plot = p,
      width = 10,
      height = 8,
      dpi = 300
    )
    cat("\nPlot saved as:", file_name, "\n")
  }
  
  # ===============================
  # CONVERT IMPORTANT GENES TO SYMBOLS
  # ===============================
  convert_ensembl_to_symbols <- function(var_imp_data, mart = NULL) {
    if (is.null(mart)) {
      mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    }
    
    ensembl_ids <- var_imp_data$Gene
    gene_map <- getBM(
      attributes = c("ensembl_gene_id", "hgnc_symbol"),
      filters = "ensembl_gene_id",
      values = ensembl_ids,
      mart = mart
    )
    
    # Create new symbols
    new_genes <- sapply(ensembl_ids, function(id) {
      symbol <- gene_map$hgnc_symbol[match(id, gene_map$ensembl_gene_id)]
      ifelse(symbol != "" & !is.na(symbol), symbol, id)
    })
    
    var_imp_data$Gene_Symbol <- new_genes
    var_imp_data$Ensembl_ID <- ensembl_ids
    
    return(var_imp_data)
  }
  
  important_genes_symbols <- convert_ensembl_to_symbols(important_vars, mart)
  
  # ===============================
  # RETURN RESULTS
  # ===============================
  return(list(
    plot = p,
    important_vars = important_genes_symbols,
    elbow_point = elbow$elbow_index,
    cumulative_importance = elbow$cumulative_importance,
    selection_method = "Maximum curvature (traditional elbow)",
    pd_genes = pd_genes_info,
    highlighted_genes = list(
      pd_genes = pd_genes_info$highlight_genes,
      manual_genes = highlight_genes,
      all_highlighted = all_highlight_genes
    )
  ))
}

# -----------------------------------------------------------------------------
# 4.4 Function to summarize analysis results
# -----------------------------------------------------------------------------
summarize_analysis_results <- function(results, domain_name) {
  cat("\n", rep("=", 60), "\n", sep = "")
  cat("FINAL SUMMARY - ", domain_name, "\n", sep = "")
  cat(rep("=", 60), "\n", sep = "")
  
  cat("\n1. GENE SELECTION:\n")
  cat("   Method: Maximum curvature (traditional elbow)\n")
  cat("   Selected genes:", results$elbow_point, "\n")
  cat("   Cumulative importance:", round(results$cumulative_importance * 100, 1), "%\n")
  
  cat("\n2. PARKINSON'S DISEASE GENES IDENTIFIED:\n")
  for (list_name in names(results$pd_genes$counts)) {
    cat(sprintf("   %-15s: %4d genes\n", list_name, results$pd_genes$counts[[list_name]]))
  }
  
  cat("\n3. HIGHLIGHTED GENES:\n")
  cat("   PD genes highlighted:", length(results$highlighted_genes$pd_genes), "\n")
  if (length(results$highlighted_genes$manual_genes) > 0) {
    cat("   Manual genes highlighted:", length(results$highlighted_genes$manual_genes), "\n")
  }
  
  cat("\n4. TOP 10 IMPORTANT GENES:\n")
  top_genes <- head(results$important_vars[, c("Gene_Symbol", "Importance")], 10)
  for (i in 1:nrow(top_genes)) {
    is_pd <- top_genes$Gene_Symbol[i] %in% results$highlighted_genes$pd_genes
    pd_marker <- ifelse(is_pd, " [PD]", "")
    cat(sprintf("   %2d. %-15s: %.6f%s\n", 
                i, top_genes$Gene_Symbol[i], top_genes$Importance[i], pd_marker))
  }
  
  cat(rep("=", 60), "\n\n")
}

# =============================================================================
# 5. EXECUTE VARIABLE IMPORTANCE ANALYSIS
# =============================================================================
cat("\n========================================\n")
cat("5. EXECUTING VARIABLE IMPORTANCE ANALYSIS\n")
cat("========================================\n\n")

# ----------------------------------------------------------------------------
# MOTOR DOMAIN ANALYSIS
# ----------------------------------------------------------------------------
cat("\n>>> MOTOR DOMAIN VARIABLE IMPORTANCE <<<\n")

motor_importance <- analyze_variable_importance(
  results_object = best_motor_model,
  domain_name = "Motor",
  additional_vars = 100,
  highlight_genes = NULL,
  pd_min_freq = 1,
  save_plot = TRUE,
  output_dir = ".",
  file_prefix = "Motor_elbow"
)

print(motor_importance$plot)
summarize_analysis_results(motor_importance, "Motor")

# Save detailed results
write.csv(
  motor_importance$important_vars,
  "PPMI_motor_VarImp_elbow_results.csv", 
  row.names = FALSE
)
cat("Motor results saved to: PPMI_motor_VarImp_elbow_results.csv\n")

# Save PD genes specific to motor
if (length(motor_importance$highlighted_genes$pd_genes) > 0) {
  pd_motor_genes <- data.frame(
    Gene_Symbol = motor_importance$highlighted_genes$pd_genes,
    Domain = "Motor",
    Selection_Method = "Maximum curvature"
  )
  write.csv(pd_motor_genes, "PD_genes_Motor_elbow.csv", row.names = FALSE)
  cat("Motor PD genes saved to: PD_genes_Motor_elbow.csv\n")
}

# ----------------------------------------------------------------------------
# NON-MOTOR DOMAIN ANALYSIS
# ----------------------------------------------------------------------------
cat("\n>>> NON-MOTOR DOMAIN VARIABLE IMPORTANCE <<<\n")

non_motor_importance <- analyze_variable_importance(
  results_object = best_non_motor_model,
  domain_name = "Non-Motor",
  additional_vars = 100,
  highlight_genes = NULL,
  pd_min_freq = 1,
  save_plot = TRUE,
  output_dir = ".",
  file_prefix = "NonMotor_elbow"
)

print(non_motor_importance$plot)
summarize_analysis_results(non_motor_importance, "Non-Motor")

# Save detailed results
write.csv(
  non_motor_importance$important_vars,
  "PPMI_non_motor_VarImp_elbow_results.csv", 
  row.names = FALSE
)
cat("Non-Motor results saved to: PPMI_non_motor_VarImp_elbow_results.csv\n")

# Save PD genes specific to non-motor
if (length(non_motor_importance$highlighted_genes$pd_genes) > 0) {
  pd_non_motor_genes <- data.frame(
    Gene_Symbol = non_motor_importance$highlighted_genes$pd_genes,
    Domain = "Non-Motor",
    Selection_Method = "Maximum curvature"
  )
  write.csv(pd_non_motor_genes, "PD_genes_NonMotor_elbow.csv", row.names = FALSE)
  cat("Non-Motor PD genes saved to: PD_genes_NonMotor_elbow.csv\n")
}

# ----------------------------------------------------------------------------
# COMPARATIVE ANALYSIS
# ----------------------------------------------------------------------------
cat("\n", rep("=", 70), "\n", sep = "")
cat("FINAL COMPARATIVE ANALYSIS BETWEEN DOMAINS\n")
cat(rep("=", 70), "\n\n")

comparison_df <- data.frame(
  Domain = c("Motor", "Non-Motor"),
  Selection_Method = c("Maximum curvature", "Maximum curvature"),
  Selected_Genes = c(motor_importance$elbow_point, non_motor_importance$elbow_point),
  Cumulative_Importance = c(
    round(motor_importance$cumulative_importance * 100, 1),
    round(non_motor_importance$cumulative_importance * 100, 1)
  ),
  PD_Genes_Identified = c(
    length(motor_importance$highlighted_genes$pd_genes),
    length(non_motor_importance$highlighted_genes$pd_genes)
  ),
  Top_Gene = c(
    motor_importance$important_vars$Gene_Symbol[1],
    non_motor_importance$important_vars$Gene_Symbol[1]
  ),
  Top_Gene_Is_PD = c(
    motor_importance$important_vars$Gene_Symbol[1] %in% motor_importance$highlighted_genes$pd_genes,
    non_motor_importance$important_vars$Gene_Symbol[1] %in% non_motor_importance$highlighted_genes$pd_genes
  )
)

print(comparison_df)

# Save comparison
write.csv(comparison_df, "Comparison_Motor_vs_NonMotor_elbow.csv", row.names = FALSE)
cat("\nComparison saved to: Comparison_Motor_vs_NonMotor_elbow.csv\n")

# Create variables for downstream use
VarImpMotor <- motor_importance$important_vars$Gene_Symbol
VarImpNonMotor <- non_motor_importance$important_vars$Gene_Symbol

cat("\n=== VARIABLES CREATED FOR DOWNSTREAM ANALYSIS ===\n")
cat("VarImpMotor_elbow:", length(VarImpMotor), "genes\n")
cat("VarImpNonMotor_elbow:", length(VarImpNonMotor), "genes\n")

# Save all results
save(VarImpMotor, VarImpNonMotor, 
     motor_importance, non_motor_importance,
     file = "Variable_Importance_Results_elbow.RData")

cat("\nAll results saved to: Variable_Importance_Results_elbow.RData\n")

# =============================================================================
# 6. DIFFERENTIAL EXPRESSION ANALYSIS
# =============================================================================
cat("\n========================================\n")
cat("6. DIFFERENTIAL EXPRESSION ANALYSIS\n")
cat("========================================\n\n")

#' Perform differential expression analysis for progression trajectories
#'
#' Compares gene expression between trajectory groups using DESeq2
#'
#' @param data Data frame with expression counts and group labels
#' @param group_column Column name containing group labels
#' @param metadata_columns Columns to retain as metadata
#' @param visits Visits to include in analysis
#' @param base_column_name Base name for group column
#' @return List with differential expression results
analyze_differential_expression <- function(data, 
                                            group_column, 
                                            metadata_columns = c("sample_id"), 
                                            visits = c("M0"),
                                            base_column_name = "cluster") {
  
  cat("\nAnalyzing differential expression for", group_column, "\n")
  
  # 1. Filter data for selected visits
  filtered_data <- data[data$visit_name %in% visits, ]
  
  # Ensure group column is factor
  filtered_data[[group_column]] <- factor(filtered_data[[group_column]])
  
  # 2. Extract count matrix (assuming genes start from column 4)
  countData <- t(as.matrix(filtered_data[, 4:ncol(filtered_data)]))
  countData <- as.matrix(countData)
  
  # 3. Convert Ensembl IDs to gene symbols
  ensembl_ids <- rownames(countData)
  
  gene_map <- getBM(
    attributes = c("ensembl_gene_id", "hgnc_symbol"),
    filters = "ensembl_gene_id", 
    values = ensembl_ids,
    mart = mart
  )
  
  new_rownames <- sapply(ensembl_ids, function(id) {
    symbol <- gene_map$hgnc_symbol[match(id, gene_map$ensembl_gene_id)]
    ifelse(symbol != "" & !is.na(symbol), symbol, id)
  })
  
  rownames(countData) <- new_rownames
  countData <- as.matrix(countData)
  
  # 4. Prepare metadata for DESeq2
  group_column_clean <- gsub(paste0("^", base_column_name, "_"), "", group_column)
  metadata <- filtered_data[, c(metadata_columns, group_column)]
  colnames(metadata)[colnames(metadata) == group_column] <- "group"
  
  # 5. Create DESeqDataSet object
  dds <- DESeqDataSetFromMatrix(
    countData = countData,
    colData = metadata,
    design = ~ group
  )
  
  # 6. Determine number of groups and execute appropriate analysis
  num_groups <- length(levels(filtered_data[[group_column]]))
  
  if (num_groups == 2) {
    # Analysis for 2 groups
    cat("  Two-group comparison\n")
    dds <- DESeq(dds)
    
    # Get group levels
    group_levels <- levels(filtered_data[[group_column]])
    
    # Compare group 1 vs 2
    res <- results(dds, contrast = c("group", as.character(group_levels[1]), 
                                     as.character(group_levels[2])))
    
    # Convert to data.frame and add annotations
    res_df <- as.data.frame(res) %>%
      tibble::rownames_to_column("gene") %>%
      mutate(
        significant = ifelse(
          padj < 0.05 & abs(log2FoldChange) > 0.5, "DEG", "No DEG"
        ),
        overexpressed_in = case_when(
          log2FoldChange > 0.5 & padj < 0.05 ~ paste("Group", group_levels[1]),
          log2FoldChange < -0.5 & padj < 0.05 ~ paste("Group", group_levels[2]), 
          TRUE ~ "None"
        ),
        comparison = paste(group_levels[1], "vs", group_levels[2]),
        cluster_type = group_column_clean
      )
    
    return(list(results = res_df, analysis_type = "2_groups"))
    
  } else if (num_groups >= 3) {
    # Analysis for 3 or more groups
    cat("  Multi-group comparison (", num_groups, " groups)\n")
    
    # Helper function for DESeq and results
    get_DEGs <- function(dds, group_var, contrast, label, cluster_type){
      dds_tmp <- dds
      colData(dds_tmp)$grp <- factor(colData(dds_tmp)[[group_var]])
      design(dds_tmp) <- ~ grp
      dds_tmp <- DESeq(dds_tmp)
      
      res <- results(dds_tmp, contrast = c("grp", contrast[1], contrast[2]))
      res_df <- as.data.frame(res) %>%
        rownames_to_column("gene") %>%
        mutate(
          comparison = label,
          cluster_type = cluster_type,
          significant = ifelse(!is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 0.5, 
                               "DEG", "No DEG"),
          overexpressed_in = case_when(
            log2FoldChange > 0.5 & padj < 0.05 ~ contrast[1],
            log2FoldChange < -0.5 & padj < 0.05 ~ contrast[2],
            TRUE ~ "None"
          )
        )
      return(res_df)
    }
    
    # Prepare combined groups for comparisons
    group_levels <- levels(filtered_data[[group_column]])
    
    # Create all possible 1 vs rest combinations
    combined_results <- list()
    
    for (i in seq_along(group_levels)) {
      current_group <- group_levels[i]
      remaining_groups <- paste(group_levels[-i], collapse = "+")
      
      # Create grouped variable
      colData(dds)[[paste0("group_", current_group, "_vs_rest")]] <- 
        ifelse(colData(dds)$group == current_group, 
               as.character(current_group), 
               remaining_groups)
      
      # Execute analysis for this comparison
      res_temp <- get_DEGs(
        dds, 
        paste0("group_", current_group, "_vs_rest"), 
        c(as.character(current_group), remaining_groups),
        paste(current_group, "vs", remaining_groups),
        group_column_clean
      )
      
      combined_results[[i]] <- res_temp
    }
    
    # Combine all results
    res_all <- bind_rows(combined_results)
    
    return(list(results = res_all, analysis_type = "multi_groups"))
    
  } else {
    stop("At least 2 groups are required for analysis")
  }
}

# =============================================================================
# 7. PREPARE DATA FOR DIFFERENTIAL EXPRESSION
# =============================================================================
cat("\n========================================\n")
cat("7. PREPARING DATA FOR DIFFERENTIAL EXPRESSION\n")
cat("========================================\n\n")

# -----------------------------------------------------------------------------
# 7.1 Load raw counts and cluster assignments
# -----------------------------------------------------------------------------
cat("7.1 Loading raw counts and cluster assignments...\n")

# Load raw counts (not VST-transformed)
PPMI_raw_counts <- readRDS(file.path(base_dir, "1aFase_ReducirGenes", "PPMI_exp")) %>% 
  rownames_to_column(.data = ., "sample_id") %>%
  filter(sample_id %in% rownames(PPMI_exp_M0)) %>%
  dplyr::select(all_of(c("sample_id", genes_present)))

cat("  - Raw counts dimensions:", dim(PPMI_raw_counts), "\n")

# Get cluster assignments from best models
motor_cluster_assignments <- best_motor_model$ClusteringResults$cluster_assignments %>%
  rename(cluster_motor = cluster)

non_motor_cluster_assignments <- best_non_motor_model$ClusteringResults$cluster_assignments %>%
  rename(cluster_non_motor = cluster)

# Merge with severity scores to get visit information
Covariates_clusterMotor <- severity_scores %>%
  filter(Visit == 0) %>%
  left_join(motor_cluster_assignments, by = "participant_id") %>%
  mutate(cluster = cluster_motor) %>%
  dplyr::select(sample_id, cluster, visit_name)

Covariates_clusterNonMotor <- severity_scores %>%
  filter(Visit == 0) %>%
  left_join(non_motor_cluster_assignments, by = "participant_id") %>%
  mutate(cluster = cluster_non_motor) %>%
  dplyr::select(sample_id, cluster, visit_name)

cat("\n  - Motor clusters:", 
    paste(sort(unique(Covariates_clusterMotor$cluster)), collapse = ", "), "\n")
cat("  - Non-motor clusters:", 
    paste(sort(unique(Covariates_clusterNonMotor$cluster)), collapse = ", "), "\n")

# -----------------------------------------------------------------------------
# 7.2 Create expression datasets with cluster labels
# -----------------------------------------------------------------------------
cat("\n7.2 Creating expression datasets with cluster labels...\n")

genes_motor <- Covariates_clusterMotor %>%
  inner_join(PPMI_raw_counts, by = "sample_id")

genes_non_motor <- Covariates_clusterNonMotor %>%
  inner_join(PPMI_raw_counts, by = "sample_id")

cat("  - Motor dataset:", nrow(genes_motor), "samples,", 
    ncol(genes_motor) - 3, "genes\n")
cat("  - Non-motor dataset:", nrow(genes_non_motor), "samples,", 
    ncol(genes_non_motor) - 3, "genes\n")

# =============================================================================
# 8. EXECUTE DIFFERENTIAL EXPRESSION ANALYSIS
# =============================================================================
cat("\n========================================\n")
cat("8. EXECUTING DIFFERENTIAL EXPRESSION ANALYSIS\n")
cat("========================================\n\n")

# Motor domain analysis
cat("\n>>> MOTOR DOMAIN DIFFERENTIAL EXPRESSION <<<\n")
motor_DEG_results <- analyze_differential_expression(
  data = genes_motor,
  group_column = "cluster",
  visits = c("M0")
)

# Non-motor domain analysis
cat("\n>>> NON-MOTOR DOMAIN DIFFERENTIAL EXPRESSION <<<\n")
non_motor_DEG_results <- analyze_differential_expression(
  data = genes_non_motor,
  group_column = "cluster", 
  visits = c("M0")
)

# Extract significant DEGs
motor_df_DEG <- motor_DEG_results$results %>% 
  filter(significant == "DEG") %>%
  arrange(padj)

non_motor_df_DEG <- non_motor_DEG_results$results %>% 
  filter(significant == "DEG") %>%
  arrange(padj)

cat("\n=== DIFFERENTIAL EXPRESSION SUMMARY ===\n")
cat("Motor domain DEGs:", nrow(motor_df_DEG), "\n")
cat("Non-motor domain DEGs:", nrow(non_motor_df_DEG), "\n")

# Display top DEGs
cat("\nTop 10 Motor DEGs:\n")
print(head(motor_df_DEG[, c("gene", "log2FoldChange", "padj", "overexpressed_in")], 10))

cat("\nTop 10 Non-motor DEGs:\n")
print(head(non_motor_df_DEG[, c("gene", "log2FoldChange", "padj", "overexpressed_in")], 10))

# Save results
write.csv(
  motor_df_DEG %>% dplyr::select(-comparison, -cluster_type),
  "motor_DEG_results.csv", 
  row.names = FALSE
)

write.csv(
  non_motor_df_DEG %>% dplyr::select(-comparison, -cluster_type),
  "non_motor_DEG_results.csv", 
  row.names = FALSE
)

cat("\nDEG results saved to:\n")
cat("  - motor_DEG_results.csv\n")
cat("  - non_motor_DEG_results.csv\n")

# =============================================================================
# 9. SESSION INFORMATION
# =============================================================================
cat("\n========================================\n")
cat("9. SESSION INFORMATION\n")
cat("========================================\n\n")

print(sessionInfo())

cat("\n========================================\n")
cat("LONGITUDINAL TRAJECTORY ANALYSIS COMPLETED SUCCESSFULLY!\n")
cat("========================================\n")