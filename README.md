# README: Neuronal Subnetwork Clustering and Statistical Analysis

## Protocol to detect neuronal subnetworks via clustering and analyze nested data using statistical modeling
Jui-Yen Huang, Gautam Chauhan, Pei-Ying Chen, Esen Tuna, and Hui-Chen Lu (2025)  

This repository contains MATLAB and R code for clustering neuronal activity data and performing statistical analysis.

## Installation

### MATLAB
- **Timing**: < 30 min
- **Requirements**: MATLAB R2023b or later with:
  - Optimization Toolbox
  - Parallel Computing Toolbox
  - Statistics and Machine Learning Toolbox
- **Steps**:
  1. Install MATLAB if not already available.
  2. Ensure required toolboxes are installed.

### R and RStudio
- **Timing**: < 30 min
- **Steps**:
  1. Install R from [CRAN](https://cran.r-project.org/).
  2. Install RStudio from [RStudio](https://www.rstudio.com/products/rstudio/).
  3. Install R packages:
     ```R
     install.packages(c("tidyverse", "janitor", "car", "sjPlot", "tictoc", "skimr", "here", "svglite", "readxl", "gridExtra", "tinytex", "extrafont", "ggpubr", "ggsci"))
     ```

## Hardware Requirements
- **CPU**: Intel® Core™ i5
- **RAM**: 32 GB
- **OS**: Windows 10, MacOS, Linux
- **Software**: MATLAB, RStudio

## Usage

### MATLAB: Clustering Analysis
1. **Setup Paths**:
   - Add custom function paths:
     ```matlab
     addpath 'Github\SubNetwork_search_multilevel2025\Clustering\Function'
     addpath 'Github\SubNetwork_search_multilevel2025\Clustering\Function_CommDetec'
     ```
   - Dependencies for Louvain community detection available from [Brain Connectivity Toolbox](https://sites.google.com/site/bctnet/).

2. **Activate Parallel Computing**:
   ```matlab
   numCores = feature('numcores');
   wish_coresToUse = 8;
   maxCoresToUse = min(wish_coresToUse, numCores);
   if isempty(gcp('nocreate'))
       fprintf('Starting parallel pool with %d cores...\n', maxCoresToUse);
       parpool('local', maxCoresToUse);
   else
       currentPool = gcp('nocreate');
       fprintf('Parallel pool already running with %d workers.\n', currentPool.NumWorkers);
   end
   ```

3. **Prepare Input Data**:
   - Use binarized neuronal activity matrix (rows: neurons, columns: time frames).
   - Example:
     ```matlab
     input_path = 'Github\SubNetwork_search_multilevel2025\Clustering_Example';
     time_per_frame = 0.065;
     detected_events = readmatrix(fullfile(input_path,'detected_events.xlsx')) ~= 0;
     start_time = readmatrix(fullfile(input_path,'start_time.csv'));
     end_time = readmatrix(fullfile(input_path,'end_time.csv'));
     Race = detected_events(:, floor(start_time/time_per_frame):floor(end_time/time_per_frame)-1);
     ```

4. **Run Clustering**:
   - Use `compute_clustering_loop.mlx` to loop through session folders (e.g., folders with "C2-1").
   - Supported algorithms: k-means, Louvain (uniform), DBSCAN.
   - Example:
     ```matlab
     run_Kmean(input_path, Race, 'CovM', 100, 5000);
     run_CommDetect_uniform(input_path, Race, 'CovM', 500, 5000);
     run_DBSCAN(input_path, Race, 'CovM', 5000);
     ```

5. **Plot Results**:
   - Plot k-means clusters:
     ```matlab
     plot_Kmean_clusters(fullfile(input_path, 'Output_Kmean_CovM/all.mat'), false, true);
     ```

### R: Statistical Analysis
1. **Setup**:
   - Download `Analysis_Develop_WT` project folder from GitHub.
   - Open `Run_Develop_Kmean.R` in RStudio.
   - Set working directory:
     ```R
     setwd("your/project/folder/path/Analysis_Develop_WT")
     ```

2. **Load Packages and Data**:
   ```R
   library(readxl)
   library(janitor)
   library(dplyr)
   library(car)
   library(sjPlot)
   library(tictoc)
   library(skimr)
   library(here)
   source(here("src", "Develop_Function.R"))
   df <- read_xlsx(here("data", "Kmean_CovM_SummaryDev.xlsx"), sheet = "Sheet1")
   df <- clean_names(df)
   ```

3. **Statistical Modeling**:
   - Define response and covariate variables:
     ```R
     rv_list <- names(df)[8:17]
     cv_list <- c("sex", "f.layer", "f.age")
     ```
   - Filter data (e.g., exclude location 4):
     ```R
     df_lite <- df %>%
       filter(location != 4) %>%
       mutate(f.age = factor(age), f.layer = factor(layer)) %>%
       select(all_of(cv_list), all_of(rv_list), "subject_id", "location")
     ```
   - Run three-way mixed-effects models and ANOVA:
     ```R
     mod_list <- lapply(rv_list, lmer_three_way)
     mod_list_select <- lapply(rv_list, model_select)
     anova_list_select <- lapply(mod_list_select, Anova, type = 3, test.statistic = "F")
     ```

4. **Bootstrap Analysis**:
   - For linear models:
     ```R
     lm_boot_list <- lapply(mod_list_lm, lm_bootstrap, R = 10000)
     ```
   - For mixed-effects models:
     ```R
     cl <- makeCluster(10)
     registerDoParallel(cl)
     case_boot_list <- lapply(mod_list_lmer, case_bootstrap, b1 = 1000, b2 = 10)
     stopCluster(cl)
     ```

5. **Visualize Results**:
   - Use `Visualizing_result.Rmd` to generate plots and reports.
   - Example plotting:
     ```R
     plot_raw_boot_facet_layer("total_cell_number", y_breaks = seq(0, 400, 100), y_limits = c(0, 400))
     ggsave("total_cell_number.svg", plot = last_plot(), path = here("output_plot"), width = 3.5, height = 2, units = "in")
     ```
   - Plot multiple variables:
     ```R
     variables_to_plot <- c('total_cell_number','n_cls_before_stat','silhs_mean_before_stat','n_cls')
     y_axis_params <- list(
       total_cell_number = list(y_breaks = seq(0, 400, by = 100), y_limits = c(0, 400)),
       n_cls_before_stat = list(y_breaks = seq(0, 100, by = 5), y_limits = c(0, 25)),
       silhs_mean_before_stat = list(y_breaks = seq(0, 100, by = 0.25), y_limits = c(0, 1)),
       n_cls = list(y_breaks = seq(0, 100, by = 5), y_limits = c(0, 15))
     )
     ```

6. **Generate Report**:
   - Knit `Visualizing_result.Rmd` to HTML using RStudio’s “Knit” button.
   - Includes animal/image counts, model summaries, and statistical reports.

## Expected Outcomes
- **MATLAB**: Outputs include `all.mat` for k-means, Louvain, and DBSCAN, with figures like `Event_cluster`, `Neuronal_Cluster`, and `Neuronal_Cluster_colored`.
- **R**: HTML report with raw data plots, bootstrapped estimates, and statistical summaries (main and interaction effects).

## Limitations
- Test similarity measures and algorithms on a data subset first.
- DBSCAN is faster; k-means and Louvain are resource-intensive.
- Statistical analysis assumes specific covariates (sex, layer, age); adapt R code for your dataset.
- Consult a statistician for complex model selection.
- Explore [R Markdown documentation](https://rmarkdown.rstudio.com/) for advanced formatting.

## Troubleshooting
- **MATLAB**: Ensure correct paths and sufficient computational resources. Check for missing input files (e.g., `detected_events.xlsx`, `start_time.csv`).
- **R**: Verify working directory and package installations. Check for missing `Develop_Function.R` or data files.
- Refer to paper for detailed troubleshooting steps.

## Contact
For issues, please open a GitHub issue.