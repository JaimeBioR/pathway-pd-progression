# PPMI Transcriptomic Analysis Pipeline

This repository contains the complete analysis pipeline for the paper "Parkinson’s Disease motor and non-motor progression models emerge from pathway-level transcriptomics". The pipeline processes transcriptomic and clinical data from the Parkinson's Progression Markers Initiative (PPMI) to derive and validate severity scores for motor and non-motor domains.


## Overview

The pipeline consists of five main scripts that sequentially process data from raw expression values to validated severity scores:

    1-Preprocessing.R - Data preprocessing and quality control

    2-Clustering_Score_Generation.R - Pathway-specific clustering and severity signal generation

    3-Domain_Characterization.R - Calculation of UTSS for motor and non-motor domains

    4-Clustering_Trajectories_and_Prediction.R - Longitudinal trajectory analysis and prediction

    5-Replication_Genetic_Cohort.R - Validation across genetic cohorts and healthy controls
    
    
## Methodology

### Pathway-Specific Transcriptomic Clustering

For each of the 336 selected Reactome pathways, patients are clustered using hierarchical clustering with Ward's method and absolute Pearson correlation distance. Clusters are generated for k values ranging from 2 to 20.

### Severity Score Generation

Pathway-specific clusters are compared pairwise across clinical covariates to identify significant differences in disease severity. When significant differences are detected (Bonferroni-corrected p < 0.05), patients in the more severe cluster receive a +1 score. Scores are summed across pathways within motor and non-motor domains to generate the Unsupervised Transcriptomic Severity Score (UTSS).

### Longitudinal Trajectory Analysis

UTSS values are transformed to percentile ranks at each visit. Disease progression trajectories are defined using feature vectors combining baseline UTSS and longitudinal changes between visits. K-means clustering identifies distinct progression patterns.

### Prediction from Baseline Transcriptomics

Random Forest models are trained to predict progression trajectories from baseline gene expression, using UTSS pathway genes as predictors. Variable importance analysis identifies key genes associated with progression.

### Validation Across Cohorts

Linear regression models trained in sporadic PD are applied to genetic PD cohorts (GBA, LRRK2, SNCA) and healthy controls to validate generalizability.
Requirements



## Requirements
### R Version

    R >= 4.0.0

### Required Packages

```{r}
# Data manipulation
tidyverse
dplyr
tidyr
purrr
stringr
readxl
readr

# Statistics
DESeq2
rstatix
broom

# Machine Learning
caret
ranger
randomForest
glmnet

# Clustering
cluster
fpc
clValid
clusterSim
factoextra
dendextend

# Visualization
ggplot2
ggpubr
patchwork
plotly
pheatmap
RColorBrewer
circlize
networkD3

# Gene Annotation
biomaRt
msigdbr
clusterProfiler
org.Hs.eg.db
gprofiler2

# Multivariate Analysis
FactoMineR
e1071

# Parallel Processing
foreach
doParallel

# Tables
DT
kableExtra

```



