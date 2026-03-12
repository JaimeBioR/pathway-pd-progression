#!/usr/bin/env Rscript

# =============================================================================
# MODEL VALIDATION ACROSS GENETIC COHORT
# PPMI Cohort - Unsupervised Transcriptomic Severity Score (UTSS)
# =============================================================================
#
# Description: This script validates the generalizability of UTSS by training
#              linear regression models in the sporadic PD cohort and applying
#              them to genetic PD cohorts (GBA, LRRK2, SNCA) and Healthy Controls.
#
#              Key analyses include:
#              1. Training linear models using clinical covariates as predictors
#                 and UTSS values as outcomes for motor and non-motor domains
#              2. Consistent scaling of predictors using training set parameters
#              3. Prediction of UTSS in genetic cohorts and healthy controls
#              4. Statistical comparison of predicted UTSS distributions between
#                 genetic cohorts and the sporadic PD reference cohort using
#                 Mann–Whitney U test
#
# Author: [Jaime Ñíguez Baeza]
# Date: 2026
# Corresponding author: [juanbot@um.es]
#
# Input files:
#   - Covariates_with_UTSS_scores.rds: Combined dataset with UTSS scores
#   - PPMI_covariates_extra: Clinical metadata for all cohorts
#
# Output files:
#   - PPMI_resultados_motor_score_clinical: Motor model results
#   - PPMI_resultados_nonmotor_score_clinical: Non-motor model results
#   - Figure_R3.png: Combined heatmap of coefficients with R² and p-values
#   - Figure_R5_n_Predicted.tiff: Predicted UTSS distributions across cohorts
#   - Model_summary_statistics_3sig.csv: Model statistics with 3 significant figures
#
# Dependencies: tidyverse, ggplot2, ggpubr, broom
# =============================================================================

# =============================================================================
# 1. INITIAL SETUP AND DATA LOADING
# =============================================================================

# Load required libraries
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(broom)
library(purrr)

# Set seed for reproducibility
set.seed(202401201)

# Define base directory (modify as needed)
base_dir <- "/home/jaimeniguez/Pipeline"

cat("========================================\n")
cat("MODEL VALIDATION ACROSS GENETIC AND EXTERNAL COHORTS\n")
cat("========================================\n\n")

# -----------------------------------------------------------------------------
# 1.1 Load data with UTSS scores
# -----------------------------------------------------------------------------
cat("1.1 Loading UTSS scores and clinical covariates...\n")

Covariates_score <- readRDS(file.path(base_dir, "results_pathways", "Covariates_with_UTSS_scores.rds"))
cat("  - Data loaded:", nrow(Covariates_score), "observations\n")

# -----------------------------------------------------------------------------
# 1.2 Define motor and non-motor covariates
# -----------------------------------------------------------------------------
cat("1.2 Defining clinical covariates for each domain...\n")

# Motor domain covariates
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

cat("  - Motor covariates:", length(Motor), "\n")
cat("  - Non-motor covariates:", length(Non_Motor), "\n\n")

# Dictionary for covariate labels in plots
covariate_labels <- c(
  # --- Motor ---
  "code_upd2hy_hoehn_and_yahr_stage" = "Hoehn & Yahr Stage",
  "mds_updrs_part_ii_summary_score" = "MDS-UPDRS-II",
  "mds_updrs_part_iii_summary_score" = "MDS-UPDRS-III",
  "Tremor_Right_UM" = "Right Tremor",
  "Tremor_Left_UM" = "Left Tremor",
  "Tremor_all_Um" = "Overall Tremor",
  "PIGD_UM" = "PIGD",
  "mod_schwab_england_pct_adl_score" = "Schwab & England ADL",
  
  # --- Non-Motor ---
  "code_upd2106_dopamine_dysregulation_syndrome_features" = "DDS",
  "mds_updrs_part_i_summary_score" = "MDS-UPDRS-I",
  "code_upd2102_hallucinations_and_psychosis" = "Hallucinations & Psychosis",
  "code_upd2101_cognitive_impairment" = "Cognitive Impairment",
  "code_upd2103_depressed_mood" = "Depressed Mood",
  "code_upd2104_anxious_mood" = "Anxious Mood",
  "code_upd2105_apathy" = "Apathy",
  "rbd_summary_score" = "RBD",
  "ess_summary_score" = "ESS",
  "moca_total_score" = "MoCA Score"
)

# =============================================================================
# 2. LINEAR REGRESSION MODELS WITH CONSISTENT SCALING
# =============================================================================
cat("========================================\n")
cat("2. LINEAR REGRESSION MODELS WITH CONSISTENT SCALING\n")
cat("========================================\n\n")

#' Multivariate linear regression analysis with consistent scaling
#'
#' Trains linear models using clinical covariates as predictors and UTSS as outcome.
#' All predictors are scaled using mean and SD from the training set, and these
#' scaling parameters are saved for consistent application to external cohorts.
#'
#' @param data Data frame with UTSS scores and clinical covariates
#' @param score_col Name of the UTSS column to predict
#' @param visits Vector of visit names to analyze
#' @param predictor_vars Vector of predictor variable names
#' @param alpha Significance threshold for p-values
#' @param cohort_name Name of the training cohort
#' @return List containing models, scaling parameters, and results
multivariate_lm_consistent <- function(data, 
                                       score_col = "UTSS_Motor", 
                                       visits = c("M0", "M12", "M24", "M36"),
                                       predictor_vars = NULL,
                                       alpha = 0.05,
                                       cohort_name = "training") {
  
  if (is.null(predictor_vars)) {
    predictor_vars <- names(data) %>%
      intersect(names(data)[sapply(data, is.numeric)]) %>%
      setdiff(c(score_col, "visit", "visit_name", "participant_id", "Order"))
  }
  
  # Lists to store results
  results <- list()
  models <- list()
  model_data <- list()
  training_predictions <- list()
  scale_params <- list()  # NEW: Save scaling parameters
  
  for (v in visits) {
    # Filter data for current visit
    df_visit <- data %>% 
      filter(visit_name == v) %>%
      dplyr::select(all_of(c(score_col, "participant_id", predictor_vars)))
    
    # ================================
    # MEDIAN IMPUTATION (TRAINING)
    # ================================
    for (var in predictor_vars) {
      if (is.numeric(df_visit[[var]])) {
        med <- median(df_visit[[var]], na.rm = TRUE)
        df_visit[[var]][is.na(df_visit[[var]])] <- med
      }
    }
    
    # Remove only if target is missing
    df_visit <- df_visit %>%
      filter(!is.na(.data[[score_col]]))
    
    # Check if sufficient data
    if (nrow(df_visit) < length(predictor_vars) + 10) {
      message(paste("Insufficient data for visit", v, "- skipping"))
      next
    }
    
    # ================================
    # SCALE PREDICTORS - SAVE PARAMETERS
    # ================================
    scale_info <- list()
    df_visit_scaled <- df_visit
    
    for (var in predictor_vars) {
      if (var %in% names(df_visit)) {
        # Calculate mean and standard deviation from training set
        mean_val <- mean(df_visit[[var]], na.rm = TRUE)
        sd_val <- sd(df_visit[[var]], na.rm = TRUE)
        
        # Save parameters
        scale_info[[var]] <- list(
          mean = mean_val,
          sd = sd_val,
          scaled = TRUE
        )
        
        # Apply scaling
        if (sd_val > 0) {
          df_visit_scaled[[var]] <- (df_visit[[var]] - mean_val) / sd_val
        } else {
          df_visit_scaled[[var]] <- 0  # If no variability
          scale_info[[var]]$scaled <- FALSE
        }
      }
    }
    
    # Save scaling parameters
    scale_params[[v]] <- scale_info
    
    # Create formula with all predictors
    formula_str <- paste(score_col, "~", paste(predictor_vars, collapse = " + "))
    formula_obj <- as.formula(formula_str)
    
    # Train multiple linear model with scaled data
    model <- lm(formula_obj, data = df_visit_scaled)
    
    # Predict on training data
    train_pred <- predict(model, newdata = df_visit_scaled, se.fit = TRUE)
    
    # Store training predictions
    training_predictions[[v]] <- data.frame(
      participant_id = df_visit$participant_id,
      visit = v,
      cohort = cohort_name,
      actual_score = df_visit[[score_col]],  # Keep original score
      predicted_score = train_pred$fit,
      se_fit = train_pred$se.fit,
      df = train_pred$df,
      residual_scale = train_pred$residual.scale,
      stringsAsFactors = FALSE
    )
    
    # Model summary
    summ <- summary(model)
    
    # Extract coefficients and statistics
    coef_df <- broom::tidy(model) %>%
      filter(term != "(Intercept)") %>%
      mutate(
        visit = v,
        significant = p.value < alpha,
        term = gsub("`", "", term)
      )
    
    # Store results
    results[[v]] <- list(
      visit = v,
      n_observations = nrow(df_visit),
      r_squared = summ$r.squared,
      adj_r_squared = summ$adj.r.squared,
      f_statistic = summ$fstatistic[1],
      f_pvalue = pf(summ$fstatistic[1], summ$fstatistic[2], summ$fstatistic[3], lower.tail = FALSE),
      coefficients = coef_df
    )
    
    models[[v]] <- model
    model_data[[v]] <- df_visit_scaled  # Save scaled data
  }
  
  # Combine all coefficients into one dataframe
  all_coefficients <- bind_rows(lapply(results, function(x) x$coefficients))
  all_training_predictions <- bind_rows(training_predictions)
  
  # Create summary by visit
  summary_df <- data.frame(
    visit = sapply(results, function(x) x$visit),
    n_obs = sapply(results, function(x) x$n_observations),
    r_squared = sapply(results, function(x) x$r_squared),
    adj_r_squared = sapply(results, function(x) x$adj_r_squared),
    f_statistic = sapply(results, function(x) x$f_statistic),
    f_pvalue = sapply(results, function(x) x$f_pvalue),
    stringsAsFactors = FALSE
  )
  
  return(list(
    models = models,                       # Full models by visit
    model_data = model_data,                # Scaled data used for training
    training_predictions = all_training_predictions, # Predictions on training cohort
    scale_parameters = scale_params,        # NEW: Scaling parameters by visit
    summary = summary_df,                   # Statistical summary
    coefficients = all_coefficients,        # Detailed coefficients
    detailed_results = results,              # Complete results by visit
    formula = formula_str,                   # Formula used
    variables = predictor_vars,              # Predictor variables
    target = score_col,                      # Target variable
    cohort_name = cohort_name                 # Training cohort name
  ))
}

# =============================================================================
# 3. VISUALIZATION FUNCTIONS
# =============================================================================
cat("========================================\n")
cat("3. VISUALIZATION FUNCTIONS\n")
cat("========================================\n\n")

#' Plot coefficients with points
plot_coefficients <- function(results) {
  
  coef_data <- results$coefficients
  
  ggplot(coef_data, aes(x = estimate, y = term, color = significant)) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
    geom_point(size = 2) +
    facet_grid(visit ~ ., scales = "free_y") +
    labs(title = "Linear model coefficients per visit",
         x = "Value",
         y = "Clinical covariates",
         color = "P-Value < 0.05") +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      axis.text.y = element_text(size = 5),
      axis.text.x = element_text(size = 9),
      strip.text = element_text(size = 10)
    )
}

#' Plot coefficients as heatmap
plot_coefficients_heatmap <- function(results) {
  
  coef_data <- results$coefficients %>%
    mutate(
      # Add asterisks based on significance level
      significance_stars = case_when(
        p.value < 0.001 ~ "***",
        p.value < 0.01 ~ "**",
        p.value < 0.05 ~ "*",
        p.value < 0.1 ~ ".",
        TRUE ~ ""
      ),
      # Create label with coefficient and asterisks
      label = sprintf("%.3f%s", estimate, significance_stars)
    )
  
  ggplot(coef_data, aes(x = visit, y = term, fill = estimate)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0, name = "Coefficients") +
    geom_text(aes(label = label), size = 5) +
    labs(title = "Coefficients per visit",
         x = NULL, y = NULL) +
    theme_minimal() +
    theme(
      aspect.ratio = 1.1,
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.caption = element_text(size = 9, hjust = 0, face = "italic")
    )
}

# =============================================================================
# 4. MOTOR AND NON-MOTOR MODEL TRAINING
# =============================================================================
cat("========================================\n")
cat("4. TRAINING MODELS FOR MOTOR AND NON-MOTOR DOMAINS\n")
cat("========================================\n\n")

#' Run complete motor analysis with consistent scaling
run_motor_analysis <- function(data, Motor_vars, Non_Motor_vars = NULL) {
  # Train model with main cohort
  model_results <- multivariate_lm_consistent(
    data, 
    score_col = "UTSS_Motor",
    predictor_vars = c(Motor_vars, Non_Motor_vars),
    cohort_name = "PPMI - Sporadic"
  )
  
  # Generate plots
  plot_motor <- plot_coefficients(model_results)
  heatmap_motor <- plot_coefficients_heatmap(model_results)
  
  # Display adjusted R²
  cat("Motor Score - Adjusted R² by visit:\n")
  print(model_results$summary$adj_r_squared)
  
  return(list(
    results = model_results,
    plot = plot_motor,
    heatmap = heatmap_motor
  ))
}

#' Run complete non-motor analysis with consistent scaling
run_non_motor_analysis <- function(data, Non_Motor_vars, Motor_vars = NULL) {
  # Train model with main cohort
  model_results <- multivariate_lm_consistent(
    data, 
    score_col = "UTSS_NonMotor",
    predictor_vars = c(Non_Motor_vars, Motor_vars),
    cohort_name = "PPMI - Sporadic"
  )
  
  # Generate plots
  plot_non_motor <- plot_coefficients(model_results)
  heatmap_non_motor <- plot_coefficients_heatmap(model_results)
  
  # Display adjusted R²
  cat("Non-Motor Score - Adjusted R² by visit:\n")
  print(model_results$summary$adj_r_squared)
  
  return(list(
    results = model_results,
    plot = plot_non_motor,
    heatmap = heatmap_non_motor
  ))
}

# =============================================================================
# 5. CREATE COMBINED FIGURE (FIGURE R3)
# =============================================================================
cat("========================================\n")
cat("5. CREATING COMBINED FIGURE (FIGURE R3)\n")
cat("========================================\n\n")

#' Create combined figure with R² and p-values (3 significant figures)
create_combined_figure <- function(motor_analysis, non_motor_analysis, save = TRUE) {
  
  # Function to format with 3 significant figures
  format_3sig <- function(x) {
    if (is.na(x)) return(NA)
    if (x == 0) return("0")
    
    # For very small numbers
    if (abs(x) < 0.001) {
      # Use scientific notation for very small numbers
      return(formatC(x, format = "e", digits = 2))
    }
    
    # Determine number of decimals based on magnitude
    magnitude <- floor(log10(abs(x)))
    
    if (magnitude >= 0) {
      # Numbers >= 1: 3 significant figures
      return(as.character(signif(x, 3)))
    } else {
      # Numbers < 1: calculate decimals for 3 significant figures
      decimals <- 2 - magnitude  # 3 significant figures
      return(as.character(round(x, decimals)))
    }
  }
  
  # Function to format p-value specifically
  format_pvalue <- function(p) {
    if (is.na(p)) return(NA)
    if (p < 0.001) {
      return("<0.001")
    } else {
      # For p-values, use 3 decimals or significant figures
      if (p >= 0.01) {
        return(sprintf("%.3f", p))
      } else {
        # For p between 0.001 and 0.01, show 3 decimals
        return(sprintf("%.3f", p))
      }
    }
  }
  
  # Function to format R² specifically
  format_r2 <- function(r2) {
    if (is.na(r2)) return(NA)
    if (r2 >= 0.01) {
      return(sprintf("%.3f", r2))
    } else {
      # Very small R²
      return(formatC(r2, format = "e", digits = 2))
    }
  }
  
  # Extract heatmaps
  heatmap_motor <- motor_analysis$heatmap
  heatmap_non_motor <- non_motor_analysis$heatmap
  
  # Rename Y axis according to dictionary
  heatmap_motor$data$term <- covariate_labels[heatmap_motor$data$term]
  heatmap_non_motor$data$term <- covariate_labels[heatmap_non_motor$data$term]
  
  # Remove significance stars (if any)
  heatmap_motor$data$significance_stars <- ""
  heatmap_non_motor$data$significance_stars <- ""
  
  # Create X axis labels with R² and model p-value (3 significant figures)
  motor_labels <- sapply(1:nrow(motor_analysis$results$summary), function(i) {
    visit <- motor_analysis$results$summary$visit[i]
    r2 <- motor_analysis$results$summary$adj_r_squared[i]
    pval <- motor_analysis$results$summary$f_pvalue[i]
    
    # Format with 3 significant figures
    r2_formatted <- format_r2(r2)
    p_formatted <- format_pvalue(pval)
    
    # Create label with line break
    paste0(visit, "\nR²=", r2_formatted, "\np=", p_formatted)
  })
  names(motor_labels) <- motor_analysis$results$summary$visit
  
  non_motor_labels <- sapply(1:nrow(non_motor_analysis$results$summary), function(i) {
    visit <- non_motor_analysis$results$summary$visit[i]
    r2 <- non_motor_analysis$results$summary$adj_r_squared[i]
    pval <- non_motor_analysis$results$summary$f_pvalue[i]
    
    # Format with 3 significant figures
    r2_formatted <- format_r2(r2)
    p_formatted <- format_pvalue(pval)
    
    # Create label with line break
    paste0(visit, "\nR²=", r2_formatted, "\np=", p_formatted)
  })
  names(non_motor_labels) <- non_motor_analysis$results$summary$visit
  
  # Apply labels to heatmap
  heatmap_motor <- heatmap_motor +
    scale_x_discrete(labels = motor_labels) +
    labs(title = "Motor") +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 8),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size = 12),
      plot.title = element_text(hjust = 0, face = "bold"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  heatmap_non_motor <- heatmap_non_motor +
    scale_x_discrete(labels = non_motor_labels) +
    labs(title = "Non-Motor") +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 8),
      axis.text.y = element_text(size = 12),
      axis.title.y = element_text(face = "bold"),
      plot.title = element_text(hjust = 0, face = "bold"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  # Save p-values separately
  Motor_significance <- heatmap_motor$data$p.value
  NonMotor_significance <- heatmap_non_motor$data$p.value
  
  # Combine heatmaps
  Figure_R3 <- ggarrange(
    heatmap_non_motor, heatmap_motor, 
    ncol = 2, nrow = 1,
    common.legend = FALSE,
    legend = "right",
    labels = c("A", "B"),
    font.label = list(size = 14, face = "bold"),
    label.x = 0,
    label.y = 1
  )
  
  # Add shared X axis title
  Figure_R3 <- annotate_figure(Figure_R3,
                               bottom = text_grob("Visit (with adjusted R² and model p-value, 3 significant figures)", 
                                                  face = "bold", size = 10))
  
  # Save if requested
  if (save) {
    ggsave("Figure_R3.png", Figure_R3, width = 18, height = 6, dpi = 300, bg = "white")
    
    # Save objects for reference
    save(heatmap_motor, heatmap_non_motor,
         Motor_significance, NonMotor_significance,
         file = "Figure_R3_objects.RData")
    
    # Save CSV with model results (with 3 significant figures)
    summary_combined <- data.frame(
      Model = rep(c("Motor", "Non-Motor"), each = nrow(motor_analysis$results$summary)),
      Visit = rep(motor_analysis$results$summary$visit, 2),
      N = c(motor_analysis$results$summary$n_obs, 
            non_motor_analysis$results$summary$n_obs),
      Adj_R2 = sapply(c(motor_analysis$results$summary$adj_r_squared,
                        non_motor_analysis$results$summary$adj_r_squared), 
                      function(x) format_r2(x)),
      F_statistic = sapply(c(motor_analysis$results$summary$f_statistic,
                             non_motor_analysis$results$summary$f_statistic),
                           function(x) format_3sig(x)),
      Model_P_value = sapply(c(motor_analysis$results$summary$f_pvalue,
                               non_motor_analysis$results$summary$f_pvalue),
                             function(x) format_pvalue(x)),
      stringsAsFactors = FALSE
    )
    
    write.csv(summary_combined, "Model_summary_statistics_3sig.csv", row.names = FALSE)
  }
  
  # Display summary in console with 3 significant figures
  cat("=== MODEL SUMMARY (3 significant figures) ===\n")
  
  cat("\nMotor Score Models:\n")
  motor_summary <- data.frame(
    Visit = motor_analysis$results$summary$visit,
    N = motor_analysis$results$summary$n_obs,
    Adj_R2 = sapply(motor_analysis$results$summary$adj_r_squared, format_r2),
    P_value = sapply(motor_analysis$results$summary$f_pvalue, format_pvalue),
    stringsAsFactors = FALSE
  )
  print(motor_summary)
  
  cat("\nNon-Motor Score Models:\n")
  nonmotor_summary <- data.frame(
    Visit = non_motor_analysis$results$summary$visit,
    N = non_motor_analysis$results$summary$n_obs,
    Adj_R2 = sapply(non_motor_analysis$results$summary$adj_r_squared, format_r2),
    P_value = sapply(non_motor_analysis$results$summary$f_pvalue, format_pvalue),
    stringsAsFactors = FALSE
  )
  print(nonmotor_summary)
  
  return(list(
    figure = Figure_R3,
    pvalues_motor = Motor_significance,
    pvalues_non_motor = NonMotor_significance,
    summary_motor = motor_summary,
    summary_non_motor = nonmotor_summary,
    labels_motor = motor_labels,
    labels_non_motor = non_motor_labels
  ))
}

# =============================================================================
# 6. PREDICTION ON GENETIC COHORTS WITH CONSISTENT SCALING
# =============================================================================
cat("========================================\n")
cat("6. PREDICTION ON GENETIC COHORTS WITH CONSISTENT SCALING\n")
cat("========================================\n\n")

#' Prepare genetic cohort data
prepare_genetic_cohort_data <- function(base_dir, target_cohorts = NULL, diagnoses = NULL) {
  
  # Load data
  Covariates_all <- as.data.frame(readRDS(file.path(base_dir, "metadata_reduced", "PPMI_covariates_extra")))
  
  # Default configuration
  if (is.null(target_cohorts)) {
    target_cohorts <- c("Healthy Control", "GBA - Aff", "LRRK2 - Aff", "LRRK2 - Unaff", 
                        "SNCA - Aff", "SNCA - Unaff", "GBA - Unaff")
  }
  
  if (is.null(diagnoses)) {
    diagnoses <- c("Case", "Control")
  }
  
  # Process mutation data
  Covariates_all <- Covariates_all %>%
    group_by(participant_id) %>%
    mutate(mutation = ifelse(all(is.na(mutation)), NA, dplyr::first(na.omit(mutation)))) %>%
    ungroup()
  
  # Filter data according to conditions
  Covariates_all <- Covariates_all %>%
    filter(case_control_other_latest %in% diagnoses &
             mutation %in% target_cohorts)
  
  # Visit statistics by mutation
  visit_stats <- Covariates_all %>%
    filter(!is.na(mutation)) %>%
    group_by(mutation) %>%
    summarise(visits = paste(unique(visit_name), collapse = ", "),
              n = n()) %>%
    arrange(mutation)
  
  print("Visit statistics by cohort:")
  print(visit_stats)
  
  return(Covariates_all)
}

#' Configure prediction cohorts
configure_prediction_cohorts <- function(data, reference_cohort = "PPMI - Sporadic", 
                                         target_cohorts = NULL) {
  
  if (is.null(target_cohorts)) {
    target_cohorts <- c("Healthy Control", "GBA - Aff", "LRRK2 - Aff", "SNCA - Aff")
  }
  
  # Normalize cohort names
  data_processed <- data %>%
    mutate(mutation = case_when(
      mutation == "Sporadic" ~ reference_cohort,
      TRUE ~ mutation
    ))
  
  # Filter data of interest
  data_filtered <- data_processed %>%
    filter(mutation %in% c(target_cohorts, reference_cohort))
  
  return(list(
    data = data_filtered,
    reference_cohort = reference_cohort,
    target_cohorts = target_cohorts
  ))
}

#' Verify consistent scaling
verify_consistent_scaling <- function(model_results) {
  cat("=== CONSISTENT SCALING VERIFICATION ===\n\n")
  
  for (v in names(model_results$scale_parameters)) {
    cat(paste("Visit:", v, "\n"))
    
    scale_info <- model_results$scale_parameters[[v]]
    
    for (var in names(scale_info)) {
      params <- scale_info[[var]]
      cat(paste("  ", var, ":\n"))
      cat(paste("    Mean (training):", round(params$mean, 4), "\n"))
      cat(paste("    SD (training):", round(params$sd, 4), "\n"))
      cat(paste("    Scaling applied:", params$scaled, "\n\n"))
    }
    
    # Verify scaled training data
    training_data <- model_results$model_data[[v]]
    if (!is.null(training_data)) {
      cat("  Scaled data statistics (training):\n")
      for (var in names(scale_info)) {
        if (var %in% names(training_data)) {
          mean_scaled <- mean(training_data[[var]], na.rm = TRUE)
          sd_scaled <- sd(training_data[[var]], na.rm = TRUE)
          cat(paste("    ", var, "- Mean:", round(mean_scaled, 4), 
                    ", SD:", round(sd_scaled, 4), "\n"))
        }
      }
    }
    cat("\n")
  }
  
  cat("✓ Verification complete.\n")
}

#' Predict on new cohort with consistent scaling
predict_on_cohort_scaled <- function(model_results, new_cohort, 
                                     cohort_name = "test", visits = NULL,
                                     max_na_percent = 5) {
  
  if (is.null(visits)) {
    visits <- names(model_results$models)
  }
  
  predictions <- list()
  
  # Function to calculate NA percentage
  percent_na <- function(df) {
    colMeans(is.na(df)) * 100
  }
  
  # Function to calculate mode
  calculate_mode <- function(x) {
    ux <- unique(na.omit(x))
    if (length(ux) == 0) return(NA)
    ux[which.max(tabulate(match(x, ux)))]
  }
  
  for (v in visits) {
    # Check if model exists for this visit
    if (!v %in% names(model_results$models)) {
      message(paste("No model for visit", v, "- skipping"))
      next
    }
    
    original_model <- model_results$models[[v]]
    scale_info <- model_results$scale_parameters[[v]]  # Get scaling parameters
    
    # Filter new cohort data for this visit
    df_visit_new <- new_cohort %>%
      filter(visit_name == v) %>%
      dplyr::select(all_of(c("participant_id", model_results$variables))) 
    
    # 1. Identify variables with many NAs in new data
    na_percent <- percent_na(df_visit_new)
    vars_with_many_nas <- names(na_percent[na_percent > max_na_percent])
    
    # Available variables for prediction (excluding problematic ones)
    available_vars <- setdiff(model_results$variables, vars_with_many_nas)
    
    if (length(vars_with_many_nas) > 0) {
      message(paste("Visit", v, ": Variables excluded due to >", max_na_percent, "% NAs:", 
                    paste(vars_with_many_nas, collapse = ", ")))
      message(paste("Visit", v, ": Available variables:", length(available_vars), "of", 
                    length(model_results$variables)))
    }
    
    # Check that we have at least some variables for prediction
    if (length(available_vars) == 0) {
      warning(paste("Visit", v, ": No variables available after excluding those with many NAs"))
      next
    }
    
    # 2. Retrain model only with available variables
    # Get original training data
    training_data <- model_results$model_data[[v]]
    
    # Filter only available variables
    training_data_adapted <- training_data %>%
      dplyr::select(all_of(c(model_results$target, available_vars))) %>%
      na.omit()
    
    # Retrain model
    formula_adapted <- as.formula(paste(model_results$target, "~", 
                                        paste(available_vars, collapse = " + ")))
    
    model_adapted <- lm(formula_adapted, data = training_data_adapted)
    
    # 3. Prepare data for prediction WITH CONSISTENT SCALING
    df_visit_for_prediction <- df_visit_new %>%
      dplyr::select(all_of(c("participant_id", available_vars)))
    
    # Apply same scaling as in training
    for (col_name in available_vars) {
      if (col_name %in% names(scale_info)) {
        scale_params <- scale_info[[col_name]]
        
        if (scale_params$scaled && scale_params$sd > 0) {
          # SCALE WITH THE SAME PARAMETERS FROM TRAINING
          df_visit_for_prediction[[col_name]] <- 
            (df_visit_for_prediction[[col_name]] - scale_params$mean) / scale_params$sd
        }
      }
      
      # Impute missing values after scaling (if any)
      if (any(is.na(df_visit_for_prediction[[col_name]]))) {
        if (is.numeric(df_visit_for_prediction[[col_name]])) {
          # Impute with 0 (scaled value corresponding to training mean)
          df_visit_for_prediction[[col_name]][is.na(df_visit_for_prediction[[col_name]])] <- 0
        } else {
          mode_val <- calculate_mode(df_visit_for_prediction[[col_name]])
          if (!is.na(mode_val)) {
            df_visit_for_prediction[[col_name]][is.na(df_visit_for_prediction[[col_name]])] <- mode_val
          }
        }
      }
    }
    
    # 4. Predict with adapted model
    tryCatch({
      predictions_obj <- predict(model_adapted, newdata = df_visit_for_prediction, se.fit = TRUE)
      
      predictions[[v]] <- data.frame(
        participant_id = df_visit_for_prediction$participant_id,
        visit = v,
        cohort = cohort_name,
        predicted_score = predictions_obj$fit,
        se_fit = predictions_obj$se.fit,
        df = predictions_obj$df,
        residual_scale = predictions_obj$residual.scale,
        n_variables_used = length(available_vars),
        variables_used = paste(available_vars, collapse = ", "),
        variables_excluded = if (length(vars_with_many_nas) > 0) {
          paste(vars_with_many_nas, collapse = ", ")
        } else {
          "None"
        },
        # NEW: Add scaling information
        scaling_applied = TRUE,
        scaling_parameters = paste("Training mean/SD for", 
                                   length(available_vars), "variables"),
        stringsAsFactors = FALSE
      )
      
      message(paste("✓ Visit", v, ": Successful predictions using consistent scaling for", 
                    length(available_vars), "variables for", 
                    nrow(df_visit_for_prediction), "participants"))
      
    }, error = function(e) {
      message(paste("✗ Error in visit", v, ":", e$message))
    })
  }
  
  # Combine all predictions
  if (length(predictions) > 0) {
    all_predictions <- bind_rows(predictions)
    return(all_predictions)
  } else {
    warning("Could not make predictions for any visit")
    return(NULL)
  }
}

#' Predict on all cohorts with consistent scaling
predict_all_cohorts_scaled <- function(model_results, data, visits = c("M0","M12","M24","M36"),
                                       max_na_percent = 10, target_cohorts = NULL,
                                       reference_cohort = "PPMI - Sporadic") {
  
  if (is.null(target_cohorts)) {
    target_cohorts <- c("Healthy Control", "GBA - Aff", "LRRK2 - Aff", "SNCA - Aff")
  }
  
  # Take only target cohorts present in data
  cohorts_present <- intersect(unique(data$mutation), c(target_cohorts, reference_cohort))
  
  predictions_by_cohort <- map(cohorts_present, function(c) {
    # Only predict for target cohorts (exclude reference if not in target_cohorts)
    if (!(c %in% target_cohorts) && c != reference_cohort) return(NULL)
    
    data_filtered <- data %>% filter(mutation == c, visit_name %in% visits)
    if (nrow(data_filtered) == 0) return(NULL)
    
    predict_on_cohort_scaled(
      model_results = model_results,
      new_cohort = data_filtered,
      cohort_name = c,
      visits = visits,
      max_na_percent = max_na_percent
    )
  })
  
  bind_rows(predictions_by_cohort)
}

#' Prepare combined predictions
prepare_combined_predictions <- function(model_results, external_predictions, 
                                         reference_cohort = "PPMI - Sporadic",
                                         plot_visits = c("M0", "M12", "M24", "M36")) {
  
  # Prepare training predictions
  train_preds <- model_results$training_predictions %>%
    dplyr::select(participant_id, visit, cohort, predicted_score, actual_score) %>%
    mutate(cohort = ifelse(cohort == "Sporadic", reference_cohort, cohort))
  
  # Combine all predictions
  all_preds <- bind_rows(
    train_preds,
    external_predictions %>% mutate(actual_score = NA)
  ) %>%
    filter(visit %in% plot_visits)
  
  return(all_preds)
}

# =============================================================================
# 7. VISUALIZATION FUNCTIONS FOR PREDICTIONS
# =============================================================================
cat("========================================\n")
cat("7. VISUALIZATION FUNCTIONS FOR PREDICTIONS\n")
cat("========================================\n\n")

#' General boxplot of predictions
plot_predictions_general <- function(all_preds, title = "Severity predictions") {
  
  ggplot(all_preds, aes(x = interaction(cohort, visit), y = predicted_score, fill = cohort)) +
    geom_boxplot(alpha = 0.6) +
    labs(
      title = title,
      x = "Cohort and Visit",
      y = "Predicted Score"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

#' Compare predictions vs actual values
plot_predictions_vs_actual <- function(all_preds, title = "Predictions vs Actual Values") {
  
  all_preds_long <- all_preds %>%
    pivot_longer(cols = c(predicted_score, actual_score),
                 names_to = "type",
                 values_to = "score")
  
  ggplot(all_preds_long, aes(x = cohort, y = score, fill = type)) +
    geom_boxplot(alpha = 0.6, position = position_dodge(width = 0.8)) +
    labs(
      title = title,
      x = "Cohort",
      y = "Score",
      fill = "Score Type"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

#' Plot by visit with statistical tests
plot_predictions_by_visit <- function(all_preds, reference_cohort, target_cohorts,
                                      visits = c("M0", "M12", "M24", "M36")) {
  
  unique_visits <- unique(all_preds$visit)
  
  for (v in unique_visits) {
    data_visit <- all_preds %>% filter(visit == v)
    
    # Cohorts in this visit
    cohorts_in_visit <- intersect(unique(data_visit$cohort), c(target_cohorts, reference_cohort))
    
    # Build comparison list
    comparisons_visit <- lapply(setdiff(cohorts_in_visit, reference_cohort), function(c) {
      c(c, reference_cohort)
    })
    
    if (length(comparisons_visit) == 0) {
      message(paste("No comparisons for visit", v, "- skipping test."))
      next
    }
    
    # Create individual plot
    p_individual <- ggplot(data_visit %>% filter(cohort %in% c(target_cohorts, reference_cohort)),
                           aes(x = cohort, y = predicted_score, fill = cohort)) +
      geom_boxplot(alpha = 0.7) +
      labs(
        title = paste("Motor Severity Prediction -", v),
        subtitle = paste("Comparisons: each cohort vs", reference_cohort),
        x = "Cohort",
        y = "Predicted Motor Score",
        caption = "Significance: *** p < 0.001; ** p < 0.01; * p < 0.05; ns = not significant"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.subtitle = element_text(size = 10, color = "gray40"),
        plot.caption = element_text(size = 9, color = "gray40")
      ) +
      stat_compare_means(
        comparisons = comparisons_visit,
        method = "wilcox.test",
        label = "p.signif"
      )
    
    # Save individual plot
    ggsave(
      filename = paste0("cohort_predict_comparison_Motor_", v, ".png"),
      plot = p_individual,
      width = 12,
      height = 8,
      dpi = 300
    )
    message(paste("Plot saved for visit:", v))
  }
}

# =============================================================================
# 8. FIGURE R5 - PREDICTION COMPARISON BY COHORT (WITH SAMPLE SIZES)
# =============================================================================
cat("========================================\n")
cat("8. CREATING FIGURE R5 - PREDICTION COMPARISON BY COHORT\n")
cat("========================================\n\n")

#' Create Figure R5 - Prediction comparison by cohort with sample sizes
create_figure_R5 <- function(all_preds_motor, all_preds_nonmotor, 
                             save = TRUE, width = 16, height = 16) {
  
  # Configuration
  cohort_order <- c("Healthy Control", "SNCA - Aff", "GBA - Aff", "LRRK2 - Aff", "PPMI - Sporadic")
  cohort_colors <- c(
    "Healthy Control" = "#1b9e77",
    "SNCA - Aff"      = "#d95f02",
    "GBA - Aff"       = "#7570b3",
    "LRRK2 - Aff"     = "#e7298a",
    "PPMI - Sporadic" = "#66a61e"
  )
  reference_cohort <- "PPMI - Sporadic"
  plot_visits <- c("M0", "M12", "M24", "M36")
  
  # Function to calculate sample sizes by visit and cohort
  calculate_n_by_cohort <- function(all_preds) {
    n_counts <- all_preds %>%
      filter(visit %in% plot_visits, cohort %in% cohort_order) %>%
      group_by(visit, cohort) %>%
      summarise(n = n(), .groups = "drop") %>%
      complete(visit = plot_visits, cohort = cohort_order, fill = list(n = 0)) %>%
      arrange(cohort, visit)
    
    return(n_counts)
  }
  
  # Calculate sample sizes for motor and non-motor
  n_counts_motor <- calculate_n_by_cohort(all_preds_motor)
  n_counts_nonmotor <- calculate_n_by_cohort(all_preds_nonmotor)
  
  # Function to create labels with n
  create_labels_with_n <- function(n_counts_visit) {
    labels <- character(length(cohort_order))
    
    for (i in seq_along(cohort_order)) {
      cohort_name <- cohort_order[i]
      n_val <- n_counts_visit$n[n_counts_visit$cohort == cohort_name]
      
      if (length(n_val) > 0 && n_val > 0) {
        # Format short names with n
        short_name <- case_when(
          cohort_name == "Healthy Control" ~ "HC",
          cohort_name == "SNCA - Aff" ~ "SNCA",
          cohort_name == "GBA - Aff" ~ "GBA",
          cohort_name == "LRRK2 - Aff" ~ "LRRK2",
          cohort_name == "PPMI - Sporadic" ~ "sPD",
          TRUE ~ cohort_name
        )
        
        labels[i] <- sprintf("%s\n(n=%d)", short_name, n_val)
      } else {
        labels[i] <- ""
      }
    }
    
    return(labels)
  }
  
  # Helper function for significance plots
  plot_significant_predictions <- function(visit_data, ref_cohort, colors, order, 
                                           n_counts_visit, ylab = "") {
    
    # Prepare labels with n for this visit
    labels_with_n <- create_labels_with_n(n_counts_visit)
    
    visit_data <- visit_data %>% 
      mutate(cohort = factor(cohort, levels = order))
    
    cohorts_to_compare <- setdiff(levels(visit_data$cohort), ref_cohort)
    
    # Calculate p-values and keep only significant
    signif_df <- lapply(cohorts_to_compare, function(c) {
      g1 <- visit_data$predicted_score[visit_data$cohort == c]
      g2 <- visit_data$predicted_score[visit_data$cohort == ref_cohort]
      if (length(g1) > 1 && length(g2) > 1) {
        p <- wilcox.test(g1, g2)$p.value
        if (!is.na(p) && p < 0.05) {
          y_pos <- max(g1, na.rm = TRUE) + 0.05 * (max(visit_data$predicted_score, na.rm = TRUE) -
                                                     min(visit_data$predicted_score, na.rm = TRUE))
          label <- ifelse(p < 0.001, "***",
                          ifelse(p < 0.01, "**", "*"))
          return(data.frame(cohort = c, ypos = y_pos, label = label))
        }
      }
      return(NULL)
    })
    
    signif_df <- do.call(rbind, signif_df)
    
    # Create base plot
    p <- ggplot(visit_data, aes(x = cohort, y = predicted_score, fill = cohort)) +
      geom_boxplot(alpha = 0.85, width = 0.65, color = "black",
                   outlier.shape = 16, outlier.size = 1.5) +
      scale_fill_manual(values = colors, drop = FALSE) +
      scale_x_discrete(labels = labels_with_n, drop = FALSE) +  # Use labels with n
      labs(x = NULL, y = ylab, title = unique(visit_data$visit)) +
      theme_minimal(base_size = 14) +
      theme(
        axis.text.x = element_text(size = 18, vjust = 0.5),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.6),
        plot.title = element_text(face = "bold", hjust = 0.5)
      )
    
    # Add significance labels
    if (!is.null(signif_df) && nrow(signif_df) > 0) {
      p <- p + geom_text(data = signif_df, aes(x = cohort, y = ypos, label = label),
                         inherit.aes = FALSE, size = 7.5)
    }
    
    return(p)
  }
  
  # Generate plots by visit
  motor_plots <- lapply(plot_visits, function(v) {
    data_visit <- all_preds_motor %>% filter(visit == v)
    n_counts_visit <- n_counts_motor %>% filter(visit == v)
    
    if (nrow(data_visit) == 0) return(NULL)
    
    plot_significant_predictions(
      data_visit, reference_cohort, cohort_colors, cohort_order,
      n_counts_visit, ylab = "Predicted Motor Score"
    )
  })
  
  nonmotor_plots <- lapply(plot_visits, function(v) {
    data_visit <- all_preds_nonmotor %>% filter(visit == v)
    n_counts_visit <- n_counts_nonmotor %>% filter(visit == v)
    
    if (nrow(data_visit) == 0) return(NULL)
    
    plot_significant_predictions(
      data_visit, reference_cohort, cohort_colors, cohort_order,
      n_counts_visit, ylab = "Predicted Non-Motor Score"
    )
  })
  
  # Filter null plots
  motor_plots <- motor_plots[!sapply(motor_plots, is.null)]
  nonmotor_plots <- nonmotor_plots[!sapply(nonmotor_plots, is.null)]
  
  # Combine plots
  col_nonmotor <- ggarrange(plotlist = nonmotor_plots, ncol = 1, nrow = length(nonmotor_plots))
  col_motor <- ggarrange(plotlist = motor_plots, ncol = 1, nrow = length(motor_plots))
  
  combined_plot <- ggarrange(col_nonmotor, col_motor, ncol = 2, nrow = 1, widths = c(1, 1))
  
  # Column titles
  title_nonmotor <- text_grob("Predicted Non-Motor Score", size = 20, face = "bold", hjust = 0.5)
  title_motor <- text_grob("Predicted Motor Score", size = 20, face = "bold", hjust = 0.5)
  
  titles_row <- ggarrange(title_nonmotor, title_motor, 
                          ncol = 2, widths = c(1, 1), 
                          labels = c("A", "B"),
                          font.label = list(size = 12),
                          label.x = 0.05, label.y = 0.75, hjust = 0, align = "hv")
  
  # Combine final figure
  final_plot <- ggarrange(titles_row, combined_plot, nrow = 2, heights = c(0.05, 0.8))
  
  # Save if requested
  if (save) {
    ggsave(filename = "Figure_R5_n_Predicted.tiff", plot = final_plot, 
           width = width, height = height, dpi = 300, bg = "white")
    
    # Save objects for reproducibility
    save(motor_plots, nonmotor_plots, combined_plot, final_plot,
         n_counts_motor, n_counts_nonmotor,
         calculate_n_by_cohort, create_labels_with_n,
         plot_significant_predictions, cohort_order, cohort_colors, 
         reference_cohort, plot_visits,
         file = "Figure_R5_con_n_objects.RData")
  }
  
  # Display sample size table
  cat("=== SAMPLE SIZES BY COHORT AND VISIT ===\n")
  
  cat("\nMotor Scores:\n")
  print(n_counts_motor %>% pivot_wider(names_from = visit, values_from = n))
  
  cat("\nNon-Motor Scores:\n")
  print(n_counts_nonmotor %>% pivot_wider(names_from = visit, values_from = n))
  
  return(final_plot)
}

# =============================================================================
# 9. COMPLETE PIPELINE WITH CONSISTENT SCALING
# =============================================================================
cat("========================================\n")
cat("9. COMPLETE PREDICTION PIPELINE WITH CONSISTENT SCALING\n")
cat("========================================\n\n")

#' Execute complete prediction pipeline with consistent scaling
execute_prediction_pipeline_scaled <- function(base_dir, Motor_vars, Non_Motor_vars, 
                                               motor_model_results = NULL,
                                               nonmotor_model_results = NULL) {
  
  # 1. Prepare genetic cohort data
  cat("=== PREPARING GENETIC COHORT DATA ===\n")
  cohort_data <- prepare_genetic_cohort_data(base_dir)
  
  # 2. Configure prediction cohorts
  config <- configure_prediction_cohorts(cohort_data)
  data_filtered <- config$data
  reference_cohort <- config$reference_cohort
  target_cohorts <- config$target_cohorts
  
  cat("Cohorts included in analysis:", paste(unique(data_filtered$mutation), collapse = ", "), "\n")
  
  # 3. Verify scaling (optional)
  if (!is.null(motor_model_results)) {
    cat("\n=== VERIFYING SCALING FOR MOTOR MODEL ===\n")
    verify_consistent_scaling(motor_model_results)
  }
  
  # 4. Execute predictions with consistent scaling
  if (!is.null(motor_model_results)) {
    cat("=== EXECUTING MOTOR PREDICTIONS (WITH CONSISTENT SCALING) ===\n")
    motor_predictions <- predict_all_cohorts_scaled(
      model_results = motor_model_results,
      data = data_filtered,
      target_cohorts = target_cohorts,
      reference_cohort = reference_cohort
    )
    
    all_preds_motor <- prepare_combined_predictions(
      motor_model_results, motor_predictions, reference_cohort
    )
  }
  
  # 5. Execute predictions for non-motor
  if (!is.null(nonmotor_model_results)) {
    cat("=== EXECUTING NON-MOTOR PREDICTIONS (WITH CONSISTENT SCALING) ===\n")
    nonmotor_predictions <- predict_all_cohorts_scaled(
      model_results = nonmotor_model_results,
      data = data_filtered,
      target_cohorts = target_cohorts,
      reference_cohort = reference_cohort
    )
    
    all_preds_nonmotor <- prepare_combined_predictions(
      nonmotor_model_results, nonmotor_predictions, reference_cohort
    )
  }
  
  # 6. Generate Figure R5
  if (!is.null(motor_model_results) && !is.null(nonmotor_model_results)) {
    cat("=== GENERATING FIGURE R5 WITH CONSISTENT PREDICTIONS ===\n")
    figure_R5 <- create_figure_R5(all_preds_motor, all_preds_nonmotor)
  }
  
  # Return results
  results <- list(
    cohort_data = data_filtered,
    reference_cohort = reference_cohort,
    target_cohorts = target_cohorts
  )
  
  if (!is.null(motor_model_results)) {
    results$motor_predictions <- motor_predictions
    results$all_preds_motor <- all_preds_motor
  }
  
  if (!is.null(nonmotor_model_results)) {
    results$nonmotor_predictions <- nonmotor_predictions
    results$all_preds_nonmotor <- all_preds_nonmotor
  }
  
  if (exists("figure_R5")) {
    results$figure_R5 <- figure_R5
  }
  
  cat("\n=== PIPELINE COMPLETED SUCCESSFULLY ===\n")
  cat("✓ Consistent scaling applied between training and prediction.\n")
  cat("✓ Valid and comparable predictions across cohorts.\n")
  
  return(results)
}

# =============================================================================
# 10. MAIN EXECUTION - COMPLETE PIPELINE WITH CONSISTENT SCALING
# =============================================================================
cat("========================================\n")
cat("10. MAIN EXECUTION - COMPLETE PIPELINE\n")
cat("========================================\n\n")

# Run main analyses first with consistent scaling
cat(">>> TRAINING MOTOR MODEL <<<\n")
motor_analysis <- run_motor_analysis(Covariates_score, Motor, Non_Motor_vars = NULL)

cat("\n>>> TRAINING NON-MOTOR MODEL <<<\n")
non_motor_analysis <- run_non_motor_analysis(Covariates_score, Non_Motor, Motor_vars = NULL)

# Save model results
saveRDS(motor_analysis, "PPMI_motor_score_clinical_results.rds")
saveRDS(non_motor_analysis, "PPMI_nonmotor_score_clinical_results.rds")
cat("\n✓ Model results saved\n")

# Generate Figure R3
cat("\n>>> GENERATING FIGURE R3 <<<\n")
figure_R3 <- create_combined_figure(
  motor_analysis,
  non_motor_analysis,
  save = TRUE
)

# Add title and save final version
F3 <- annotate_figure(figure_R3$figure,
                      left = text_grob("Clinical covariates", 
                                       face = "bold", 
                                       size = 12,
                                       rot = 90,
                                       vjust = 0.5,
                                       hjust = 0.5))

ggsave(filename = "F3.png", plot = F3, width = 18, height = 6, dpi = 300, bg = "white")
cat("✓ Figure R3 saved as F3.png\n")

# Execute prediction pipeline with consistent scaling
cat("\n>>> EXECUTING PREDICTION PIPELINE FOR GENETIC COHORTS <<<\n")
prediction_results <- execute_prediction_pipeline_scaled(
  base_dir = base_dir,
  Motor_vars = Motor,
  Non_Motor_vars = Non_Motor,
  motor_model_results = motor_analysis$results,
  nonmotor_model_results = non_motor_analysis$results
)

cat("\n=== PIPELINE COMPLETED SUCCESSFULLY ===\n")
cat("✓ Linear models trained on sporadic PD cohort\n")
cat("✓ Consistent scaling applied for external predictions\n")
cat("✓ UTSS predicted for genetic cohorts and healthy controls\n")
cat("✓ Statistical comparisons performed using Mann-Whitney U test\n")

# =============================================================================
# 11. SESSION INFORMATION
# =============================================================================
cat("\n========================================\n")
cat("11. SESSION INFORMATION\n")
cat("========================================\n\n")

print(sessionInfo())

cat("\n========================================\n")
cat("MODEL VALIDATION ACROSS GENETIC COHORTS COMPLETED SUCCESSFULLY!\n")
cat("========================================\n")