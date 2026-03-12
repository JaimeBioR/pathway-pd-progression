#!/usr/bin/env Rscript

# =============================================================================
# TRANSCRIPTOMIC AND CLINICAL DATA PREPROCESSIONG FOR PARKINSON'S DISEASE ANALYSIS
# PPMI and PDBP Cohorts
# =============================================================================
#
# Description: This script performs comprehensive preprocessing of transcriptomic
#              and clinical data from the PPMI and PDBP cohorts, including:
#              - Metadata filtering and missing value imputation
#              - Gene filtering using elbow method
#              - Variance stabilizing transformation (DESeq2)
#              - Gene biotype stratification
#              - Pathway-based gene set selection and clustering
#
# Author: [Jaime Ñíguez Baeza]
# Date: 2026
# Corresponding author: [juanbot@um.es]
#
# Input files:
#   - PPMI_covariates_extra: Clinical metadata
#   - PPMI_exp: Raw expression matrix
#   - Subsets_pathways_Reactome: Reactome pathway gene sets
#   - mart: Ensembl mart object for gene annotation
#
# Output files:
#   - PPMI_exp_filtered_elbow: Filtered expression matrix
#   - PPMI_exp_vst: Variance-stabilized expression matrix
#   - Gene_biotypes_PPMI_and_pathways.csv: Gene annotation table
#   - clustering_full_results.rds: Complete pathway clustering results
#   - pathway_cluster_mapping.csv: Pathway-to-cluster assignments
#   - Various diagnostic plots (.png)
#
# Dependencies: tidyverse, DESeq2, biomaRt, msigdbr, ggplot2, patchwork
# =============================================================================

# =============================================================================
# 1. INITIAL SETUP AND DATA LOADING
# =============================================================================

# Load required libraries
library(tidyverse)
library(DESeq2)
library(biomaRt)
library(msigdbr)
library(ggplot2)
library(patchwork)  # For combining plots
library(DT)         # For interactive tables (optional)

# Set seed for reproducibility
set.seed(123)

# Define base directory (modify as needed)
base_dir <- getwd()  # Or specify your path

cat("========================================\n")
cat("Starting preprocessing pipeline\n")
cat("========================================\n\n")

# -----------------------------------------------------------------------------
# 1.1 Load clinical metadata
# -----------------------------------------------------------------------------
cat("Loading clinical metadata...\n")
Covariates <- as.data.frame(readRDS(file.path(base_dir, "metadata_reduced", "PPMI_covariates_extra")))
rownames(Covariates) <- Covariates$sample_id

# -----------------------------------------------------------------------------
# 1.2 Define analysis parameters
# -----------------------------------------------------------------------------
cat("Setting analysis parameters...\n")
Visit <- c("M0", "M12", "M24", "M36")           # Visits to include
cohort <- c("Sporadic")                          # Mutation status to include
Diagnosis <- c("Case")                           # Diagnosis status

# -----------------------------------------------------------------------------
# 1.3 Load pathway gene sets
# -----------------------------------------------------------------------------
cat("Loading pathway gene sets...\n")
Subsets_pathways <- readRDS(file.path(base_dir, "Pipeline", "Subsets_pathways_Reactome"))

# Calculate and display pathway size statistics
pathway_sizes <- sapply(Subsets_pathways, length)
max_pathway <- names(which.max(pathway_sizes))
max_size <- max(pathway_sizes)
min_pathway <- names(which.min(pathway_sizes))
min_size <- min(pathway_sizes)

cat("\nPathway size statistics:\n")
cat("  - Maximum size:", max_pathway, "(", max_size, "genes)\n")
cat("  - Minimum size:", min_pathway, "(", min_size, "genes)\n")
cat("  - Total pathways:", length(Subsets_pathways), "\n\n")

# =============================================================================
# 2. CLINICAL DATA PREPROCESSING
# =============================================================================
cat("========================================\n")
cat("2. CLINICAL DATA PREPROCESSING\n")
cat("========================================\n\n")

# -----------------------------------------------------------------------------
# 2.1 Filter metadata by visit, diagnosis, and mutation status
# -----------------------------------------------------------------------------
cat("2.1 Filtering metadata by visit, diagnosis, and mutation status...\n")

# Display available visits by mutation status
cat("Available visits by mutation status:\n")
Covariates %>%
  filter(!is.na(mutation)) %>%
  group_by(mutation) %>%
  summarise(visits = paste(unique(visit_name), collapse = ", "),
            n = n()) %>%
  arrange(mutation) %>%
  print()

# Propagate mutation status within participants
Covariates <- Covariates %>%
  group_by(participant_id) %>%
  mutate(mutation = ifelse(all(is.na(mutation)), NA, dplyr::first(na.omit(mutation)))) %>%
  ungroup()

# Identify participants with all required visits
valid_participants <- Covariates %>%
  group_by(participant_id) %>%
  filter(all(Visit %in% visit_name)) %>%
  pull(participant_id) %>%
  unique()

# Apply all filters
Covariates <- Covariates %>%
  filter(participant_id %in% valid_participants &
           visit_name %in% Visit &
           case_control_other_latest %in% Diagnosis &
           mutation %in% cohort)

cat("  - Participants retained:", length(unique(Covariates$participant_id)), "\n")
cat("  - Samples retained:", nrow(Covariates), "\n\n")

# -----------------------------------------------------------------------------
# 2.2 Missing value treatment
# -----------------------------------------------------------------------------
cat("2.2 Handling missing values...\n")

#' Calculate percentage of missing values per column
#' @param df Data frame
#' @return Vector with percentage of NAs per column
percent_na <- function(df) {
  colMeans(is.na(df)) * 100
}

# Remove columns with >10% missing values
PercentNA_Df <- as.data.frame(percent_na(Covariates), nm = "Percent.NAs")
PercentNA_names_clean <- rownames(subset(PercentNA_Df, Percent.NAs > 10))
Covariates_clean <- Covariates[, !colnames(Covariates) %in% PercentNA_names_clean]
Covariates <- Covariates_clean
cat("  - Removed", length(PercentNA_names_clean), "variables with >10% missing values\n")

# Impute remaining NAs
## Numeric columns: impute with median
numeric_columns <- Covariates[, sapply(Covariates, is.numeric)]
medians <- apply(numeric_columns, 2, function(x) median(x, na.rm = TRUE))

numeric_imputed <- as.data.frame(lapply(1:ncol(numeric_columns), function(i) {
  x <- numeric_columns[, i]
  x[is.na(x)] <- medians[i]
  return(x)
}), col.names = colnames(numeric_columns))

## Non-numeric columns: impute with mode
#' Calculate mode (most frequent value)
#' @param x Vector
#' @return Mode value
calculate_mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

non_numeric_columns <- Covariates[, !sapply(Covariates, is.numeric)]
non_numeric_imputed <- as.data.frame(lapply(1:ncol(non_numeric_columns), function(i) {
  x <- non_numeric_columns[, i]
  mode_val <- calculate_mode(x[!is.na(x)])
  x[is.na(x)] <- mode_val
  return(x)
}), col.names = colnames(non_numeric_columns))

# Combine imputed data
Covariates_imputed <- cbind(non_numeric_imputed, numeric_imputed)
Covariates_imputed <- Covariates_imputed[, colnames(Covariates)]
rownames(Covariates_imputed) <- Covariates_imputed$sample_id
Covariates <- Covariates_imputed

# Verify no missing values remain
na_check <- table(percent_na(Covariates) == 0)
cat("  - Variables with 0% missing values:", na_check["TRUE"], "\n")
cat("  - Variables with >0% missing values:", na_check["FALSE"], "\n\n")

# =============================================================================
# 3. TRANSCRIPTOMIC DATA PREPROCESSING
# =============================================================================
cat("========================================\n")
cat("3. TRANSCRIPTOMIC DATA PREPROCESSING\n")
cat("========================================\n\n")

# -----------------------------------------------------------------------------
# 3.1 Load and filter expression data
# -----------------------------------------------------------------------------
cat("3.1 Loading expression data...\n")
PPMI_exp <- readRDS(file.path(base_dir, "1aFase_ReducirGenes", "PPMI_exp"))
cat("  - Original dimensions:", dim(PPMI_exp)[1], "samples ×", dim(PPMI_exp)[2], "genes\n")

# Filter samples based on Covariates
PPMI_exp_case <- PPMI_exp %>%
  filter(row.names(.) %in% Covariates$sample_id)

# Merge with case information
PPMI_exp_case$sample_id <- rownames(PPMI_exp_case)
case <- Covariates[, c("case_control_other_latest", "sample_id")]
PPMI_exp_case <- merge(PPMI_exp_case, case, by = "sample_id")

# Set row names and remove non-expression columns
rownames(PPMI_exp_case) <- PPMI_exp_case$sample_id
PPMI_exp_case$sample_id <- NULL
PPMI_exp_case$case_control_other_latest <- NULL

PPMI_exp <- PPMI_exp_case
cat("  - Filtered dimensions:", dim(PPMI_exp)[1], "samples ×", dim(PPMI_exp)[2], "genes\n\n")

# -----------------------------------------------------------------------------
# 3.2 Gene filtering using elbow method
# -----------------------------------------------------------------------------
cat("3.2 Applying gene filtering with elbow method...\n")

#' Find optimal expression threshold using elbow method
#' 
#' @param expression_matrix Gene expression matrix (samples × genes)
#' @param min_samples_percent Minimum proportion of samples with expression > threshold
#' @return List with optimal threshold and columns to keep
find_optimal_expression_threshold <- function(expression_matrix, min_samples_percent = 0.8) {
  
  # Calculate number of columns retained for each threshold
  thresholds <- seq(0, 99, by = 1)
  num_columns <- sapply(thresholds, function(threshold) {
    high_expr_percentage <- colSums(expression_matrix > threshold, na.rm = TRUE) / nrow(expression_matrix)
    columns_to_keep <- names(high_expr_percentage[high_expr_percentage >= min_samples_percent & !is.na(high_expr_percentage)])
    length(columns_to_keep)
  })
  
  # Find elbow point using curvature method
  find_elbow_point <- function(x, y) {
    n <- length(x)
    coordinates <- cbind(x, y)
    start_point <- coordinates[1, ]
    end_point <- coordinates[n, ]
    
    distances <- apply(coordinates, 1, function(point) {
      numerator <- abs((end_point[2] - start_point[2]) * point[1] - 
                         (end_point[1] - start_point[1]) * point[2] + 
                         end_point[1] * start_point[2] - end_point[2] * start_point[1])
      denominator <- sqrt((end_point[2] - start_point[2])^2 + (end_point[1] - start_point[1])^2)
      numerator / denominator
    })
    
    which.max(distances)
  }
  
  elbow_index <- find_elbow_point(thresholds, num_columns)
  optimal_threshold <- thresholds[elbow_index]
  optimal_columns <- num_columns[elbow_index]
  
  # Get final columns to keep
  final_high_expr <- colSums(expression_matrix > optimal_threshold, na.rm = TRUE) / nrow(expression_matrix)
  final_columns <- names(final_high_expr[final_high_expr >= min_samples_percent & !is.na(final_high_expr)])
  
  # Plot results
  png("figures/gene_filtering_elbow.png", width = 800, height = 600, res = 120)
  plot(thresholds, num_columns, type = "l", lwd = 2, col = "steelblue",
       xlab = "Expression Threshold", ylab = "Number of Genes Retained",
       main = "Optimal Threshold Selection (Elbow Method)")
  abline(v = optimal_threshold, col = "red", lty = 2, lwd = 2)
  points(optimal_threshold, optimal_columns, col = "red", pch = 19, cex = 1.5)
  legend("topright", legend = paste("Optimal threshold:", optimal_threshold), 
         col = "red", lty = 2, lwd = 2)
  dev.off()
  
  return(list(
    optimal_threshold = optimal_threshold,
    num_columns_retained = optimal_columns,
    columns_to_keep = final_columns
  ))
}

# Apply filtering
# Uncomment to run filtering (commented out to use pre-computed results)
# optimal_results <- find_optimal_expression_threshold(PPMI_exp, min_samples_percent = 0.8)
# PPMI_exp_filtered <- PPMI_exp[, optimal_results$columns_to_keep]
# saveRDS(PPMI_exp_filtered, "PPMI_exp_filtered_elbow")

# Load pre-filtered data
PPMI_exp_filtered <- readRDS("PPMI_exp_filtered_elbow")
PPMI_exp <- PPMI_exp_filtered
cat("  - Genes retained after filtering:", ncol(PPMI_exp), "\n\n")

# -----------------------------------------------------------------------------
# 3.3 Variance stabilizing transformation
# -----------------------------------------------------------------------------
cat("3.3 Applying variance stabilizing transformation (VST)...\n")

# Uncomment to run VST
# PPMI_exp_vst <- vst(as.matrix(t(PPMI_exp)))
# saveRDS(PPMI_exp_vst, "PPMI_exp_vst")

# Load pre-computed VST data
PPMI_exp_vst <- readRDS("PPMI_exp_vst")
PPMI_exp <- as.data.frame(t(PPMI_exp_vst))
cat("  - VST complete. Final dimensions:", dim(PPMI_exp)[1], "samples ×", dim(PPMI_exp)[2], "genes\n\n")

# =============================================================================
# 4. GENE BIOTYPE STRATIFICATION
# =============================================================================
cat("========================================\n")
cat("4. GENE BIOTYPE STRATIFICATION\n")
cat("========================================\n\n")

cat("4.1 Annotating genes with biotype information...\n")

# Load Ensembl mart (pre-saved to avoid repeated downloads)
mart <- readRDS("mart")

# Get genes from expression matrix and pathways
genes_ppmi <- colnames(PPMI_exp)
genes_pathways <- unique(unlist(Subsets_pathways))
genes_intersect <- intersect(genes_ppmi, genes_pathways)

# Query biotype information from Ensembl
gene_info <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "gene_biotype"),
  filters = "ensembl_gene_id",
  values = unique(c(genes_ppmi, genes_pathways)),
  mart = mart
)
colnames(gene_info)[1] <- "Gene"

# Add pathway membership information
gene_info$Found_in_pathway <- gene_info$Gene %in% genes_pathways

# Count pathways per gene
pathway_counts <- table(unlist(Subsets_pathways))
pathway_counts <- as.data.frame(pathway_counts)
colnames(pathway_counts) <- c("Gene", "Pathway_count")
gene_info <- merge(gene_info, pathway_counts, by = "Gene", all.x = TRUE)
gene_info$Pathway_count[is.na(gene_info$Pathway_count)] <- 0
gene_info$In_PPMI_exp <- gene_info$Gene %in% genes_ppmi

# Display biotype distribution
cat("\nBiotype distribution in PPMI + pathways:\n")
biotype_table <- table(gene_info$gene_biotype[gene_info$In_PPMI_exp & gene_info$Found_in_pathway])
print(biotype_table)

# Save gene information
write.csv(gene_info, "Gene_biotypes_PPMI_and_pathways.csv", row.names = FALSE)
cat("\n  - Gene annotation saved to: Gene_biotypes_PPMI_and_pathways.csv\n")

# Create visualization of biotype distribution
biotype_df <- as.data.frame(biotype_table)
colnames(biotype_df) <- c("Gene_biotype", "Count")
biotype_df <- biotype_df %>%
  mutate(Legend_label = paste0(Gene_biotype, " (", Count, ")"))

# Bar plot
p_bar <- ggplot(biotype_df, aes(x = reorder(Gene_biotype, -Count), y = Count, fill = Legend_label)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(
    title = "Gene Biotype Distribution",
    x = "Gene biotype",
    y = "Number of Genes",
    fill = "Biotype (n)"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Pie chart
p_pie <- ggplot(biotype_df, aes(x = "", y = Count, fill = Legend_label)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  theme_void() +
  labs(title = "Biotype Gene Proportion",
       fill = "Biotype(n)")

# Save plots
ggsave("figures/gene_biotype_barplot.png", p_bar, width = 10, height = 6, dpi = 300)
ggsave("figures/gene_biotype_pie.png", p_pie, width = 8, height = 6, dpi = 300)

# Extract protein-coding genes for downstream analysis
protein_coding_genes <- gene_info[
  gene_info$gene_biotype == "protein_coding" & 
    gene_info$Found_in_pathway & 
    gene_info$In_PPMI_exp, 
]$Gene
saveRDS(protein_coding_genes, "protein_coding_genes")
cat("  - Protein-coding genes saved:", length(protein_coding_genes), "\n\n")

# =============================================================================
# 5. PATHWAY-BASED GENE SET SELECTION AND CLUSTERING
# =============================================================================
cat("========================================\n")
cat("5. PATHWAY-BASED GENE SET SELECTION\n")
cat("========================================\n\n")

cat("5.1 Loading Reactome pathways from MSigDB...\n")
# Uncomment to download pathways (commented out as likely pre-computed)
# geneset_filter <- msigdbr(species = "Homo sapiens") %>%
#   filter(gs_collection == "C2" & gs_subcollection %in% c("CP:REACTOME"))
# geneset_filter$gs_description_ont <- paste(geneset_filter$gs_subcollection, 
#                                            geneset_filter$gs_description, sep = ":")
# gene_list <- split(geneset_filter$ensembl_gene, geneset_filter$gs_description_ont)
# vector_lists <- lapply(gene_list, unique)
# saveRDS(vector_lists, "Subsets_pathways_Reactome")

cat("Example pathway:\n")
print(Subsets_pathways[123])

# =============================================================================
# 5.2 Pathway clustering to reduce redundancy
# =============================================================================
cat("\n5.2 Performing pathway clustering to select non-redundant pathways...\n")

#' Calculate distance matrix between pathways based on gene overlap
#' 
#' @param subsets_pathways List of pathway gene sets
#' @param genes_filter Optional vector of genes to filter pathways
#' @param min_genes Minimum number of genes required for a pathway
#' @return List with distance matrix and pathway information
calculate_pathway_distance_matrix <- function(subsets_pathways, genes_filter = NULL, min_genes = 5) {
  
  cat("  - Filtering pathways with at least", min_genes, "genes...\n")
  
  # Apply gene filter if provided
  if (!is.null(genes_filter)) {
    subsets_pathways_filtered <- lapply(subsets_pathways, function(pathway_genes) {
      intersect(pathway_genes, genes_filter)
    })
  } else {
    subsets_pathways_filtered <- subsets_pathways
  }
  
  # Filter by minimum size
  pathway_sizes <- sapply(subsets_pathways_filtered, length)
  valid_pathways <- which(pathway_sizes >= min_genes)
  
  cat("  - Pathways retained:", length(valid_pathways), "\n")
  
  if (length(valid_pathways) < 2) {
    stop("Need at least 2 pathways for clustering")
  }
  
  subsets_pathways_filtered <- subsets_pathways_filtered[valid_pathways]
  
  # Calculate similarity matrix (Szymkiewicz-Simpson coefficient)
  cat("  - Calculating similarity matrix...\n")
  
  calculate_ss_similarity <- function(genes1, genes2) {
    intersection <- length(intersect(genes1, genes2))
    min_size <- min(length(genes1), length(genes2))
    if (min_size == 0) return(0)
    return(intersection / min_size)
  }
  
  n_pathways <- length(subsets_pathways_filtered)
  pathway_names <- names(subsets_pathways_filtered)
  
  sim_matrix <- matrix(0, nrow = n_pathways, ncol = n_pathways)
  rownames(sim_matrix) <- pathway_names
  colnames(sim_matrix) <- pathway_names
  
  # Progress bar for large calculations
  pb <- txtProgressBar(min = 0, max = n_pathways, style = 3)
  for (i in 1:n_pathways) {
    for (j in i:n_pathways) {
      if (i == j) {
        sim_matrix[i, j] <- 1
      } else {
        sim_matrix[i, j] <- calculate_ss_similarity(
          subsets_pathways_filtered[[i]],
          subsets_pathways_filtered[[j]]
        )
        sim_matrix[j, i] <- sim_matrix[i, j]
      }
    }
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  # Convert to distance matrix
  dist_matrix <- as.dist(1 - sim_matrix)
  
  return(list(
    distance_matrix = dist_matrix,
    similarity_matrix = sim_matrix,
    pathways = pathway_names,
    pathway_genes = subsets_pathways_filtered,
    n_pathways = n_pathways,
    mean_similarity = mean(sim_matrix[upper.tri(sim_matrix)]),
    pathway_sizes = pathway_sizes[valid_pathways]
  ))
}

#' Perform hierarchical clustering of pathways
#' 
#' @param distance_data Output from calculate_pathway_distance_matrix
#' @param k Number of clusters
#' @param min_cluster_size Minimum cluster size to retain
#' @return List with clustering results
cluster_pathways <- function(distance_data, k = 20, min_cluster_size = 3) {
  
  dist_matrix <- as.dist(distance_data$distance_matrix)
  pathways <- names(distance_data$pathway_genes)
  
  if (length(pathways) < k) {
    k <- length(pathways) - 1
  }
  
  # Hierarchical clustering
  hc <- hclust(dist_matrix, method = "ward.D2")
  clusters <- cutree(hc, k = k)
  
  cluster_df <- data.frame(
    Pathway = names(clusters),
    Cluster = as.integer(clusters),
    Pathway_Size = sapply(names(clusters),
                          function(p) length(distance_data$pathway_genes[[p]]))
  )
  
  # Find representative pathway for each cluster
  sim_matrix <- 1 - as.matrix(dist_matrix)
  
  find_representative <- function(cluster_id) {
    cluster_pathways <- cluster_df$Pathway[cluster_df$Cluster == cluster_id]
    if (length(cluster_pathways) < min_cluster_size) return(NULL)
    
    avg_similarity <- sapply(cluster_pathways, function(p) {
      others <- setdiff(cluster_pathways, p)
      if (length(others) > 0) mean(sim_matrix[p, others], na.rm = TRUE) else 0
    })
    
    representative <- cluster_pathways[which.max(avg_similarity)]
    
    list(
      representative = representative,
      mean_similarity = max(avg_similarity),
      cluster_size = length(cluster_pathways),
      representative_size = cluster_df$Pathway_Size[cluster_df$Pathway == representative]
    )
  }
  
  # Compile representative information
  cluster_representatives <- list()
  for (i in 1:k) {
    rep_info <- find_representative(i)
    if (!is.null(rep_info)) {
      cluster_representatives[[i]] <- data.frame(
        Cluster_ID = i,
        Representative_Pathway = rep_info$representative,
        Cluster_Size = rep_info$cluster_size,
        Cluster_Mean_Similarity = rep_info$mean_similarity,
        Representative_Size = rep_info$representative_size
      )
    }
  }
  
  cluster_representatives_df <- do.call(rbind, cluster_representatives)
  cluster_representatives_df <- cluster_representatives_df[
    cluster_representatives_df$Cluster_Size >= min_cluster_size, ]
  
  # Re-number clusters
  if (nrow(cluster_representatives_df) > 0) {
    old_to_new <- setNames(1:nrow(cluster_representatives_df), 
                           cluster_representatives_df$Cluster_ID)
    
    filtered_clusters <- ifelse(
      clusters %in% cluster_representatives_df$Cluster_ID,
      old_to_new[as.character(clusters)],
      NA
    )
    names(filtered_clusters) <- names(clusters)
    
    cluster_df_filtered <- cluster_df[
      cluster_df$Cluster %in% cluster_representatives_df$Cluster_ID, ]
    cluster_df_filtered$Cluster <- old_to_new[as.character(cluster_df_filtered$Cluster)]
    
    cluster_representatives_df$Cluster_ID <- 1:nrow(cluster_representatives_df)
  }
  
  return(list(
    hclust_object = hc,
    clusters = filtered_clusters,
    cluster_info = cluster_df_filtered,
    cluster_representatives = cluster_representatives_df,
    similarity_matrix = sim_matrix,
    distance_matrix = dist_matrix,
    n_clusters_final = nrow(cluster_representatives_df)
  ))
}

# Execute pathway clustering
cat("\n5.3 Executing pathway clustering with automatic k selection...\n")

# Calculate distance matrix
distance_data <- calculate_pathway_distance_matrix(
  subsets_pathways = Subsets_pathways,
  genes_filter = NULL,
  min_genes = 5
)

# Evaluate different k values
k_range <- 10:500
results_k <- data.frame()

for (k in k_range) {
  if (k >= 3 && k <= length(distance_data$pathways)) {
    clustering_temp <- cluster_pathways(
      distance_data = distance_data,
      k = k,
      min_cluster_size = 1
    )
    
    if (!is.null(clustering_temp$cluster_representatives) && 
        nrow(clustering_temp$cluster_representatives) > 0) {
      
      mean_similarity <- mean(clustering_temp$cluster_representatives$Cluster_Mean_Similarity, 
                              na.rm = TRUE)
      
      results_k <- rbind(results_k, data.frame(
        k = k,
        n_clusters = nrow(clustering_temp$cluster_representatives),
        mean_similarity = mean_similarity,
        sd_similarity = sd(clustering_temp$cluster_representatives$Cluster_Mean_Similarity, 
                           na.rm = TRUE)
      ))
    }
  }
}

# Select optimal k (maximizing mean similarity)
optimal_k <- results_k$k[which.max(results_k$mean_similarity)]
optimal_ss <- max(results_k$mean_similarity)

cat("\nOptimal clustering parameters:\n")
cat("  - Optimal k:", optimal_k, "\n")
cat("  - Maximum average similarity:", round(optimal_ss, 3), "\n")

# Perform final clustering with optimal k
clustering_result <- cluster_pathways(
  distance_data = distance_data,
  k = optimal_k,
  min_cluster_size = 1
)

# Create pathway mapping
pathway_mapping <- clustering_result$cluster_info %>%
  left_join(clustering_result$cluster_representatives, 
            by = c("Cluster" = "Cluster_ID")) %>%
  dplyr::select(
    Original_Pathway = Pathway,
    Cluster_ID = Cluster,
    Representative_Pathway,
    Original_Pathway_Size = Pathway_Size,
    Cluster_Size,
    Cluster_Mean_Similarity
  )

# Save results
write.csv(pathway_mapping, "pathway_cluster_mapping.csv", row.names = FALSE)
write.csv(clustering_result$cluster_representatives, "cluster_representatives.csv", 
          row.names = FALSE)

saveRDS(list(
  clustering = clustering_result,
  k_evaluation = results_k,
  optimal_k = optimal_k
), "clustering_full_results.rds")

cat("\nPathway clustering complete:\n")
cat("  - Original pathways:", distance_data$n_pathways, "\n")
cat("  - Final clusters:", nrow(clustering_result$cluster_representatives), "\n")
cat("  - Reduction:", round(100 * (1 - nrow(clustering_result$cluster_representatives) / 
                                     distance_data$n_pathways), 1), "%\n")

# Create vector of representative pathways
create_cluster_vector <- function(pathway_mapping) {
  cluster_vector <- list()
  unique_clusters <- unique(pathway_mapping$Cluster_ID)
  
  for (cluster_id in unique_clusters) {
    representative <- pathway_mapping %>%
      filter(Cluster_ID == cluster_id) %>%
      pull(Representative_Pathway) %>%
      first()
    
    cluster_pathways <- pathway_mapping %>%
      filter(Cluster_ID == cluster_id) %>%
      pull(Original_Pathway)
    
    cluster_vector[[representative]] <- cluster_pathways
  }
  
  return(cluster_vector)
}

cluster_vector <- create_cluster_vector(pathway_mapping)
saveRDS(cluster_vector, "pathway_clusters.rds")

# =============================================================================
# 6. SESSION INFORMATION
# =============================================================================
cat("\n========================================\n")
cat("6. SESSION INFORMATION\n")
cat("========================================\n\n")

print(sessionInfo())

cat("\n========================================\n")
cat("Preprocessing pipeline completed successfully!\n")
cat("========================================\n")