# Refresh the R session to clear the environment
# This ensures that any previous objects, functions, or loaded packages are removed.
freshr::freshr()

# Set the working directory to the desired folder
# This is the folder where your project files are located.

# Verify the working directory
getwd() # Should print the path to your working directory

# Load necessary libraries
# These libraries provide various functions and tools needed for data manipulation, analysis, and visualization.
library(readxl)    # For reading Excel files
library(janitor)   # For cleaning data
library(dplyr)     # For data manipulation
library(car)       # For regression diagnostics
library(sjPlot)    # For data visualization and tabulation
library(tictoc)    # For timing code execution
library(skimr)     # For summarizing data
library(here)


source(here("src","Develop_Function.R"))
# import data from Excel file
df <-read_xlsx(here("data","Kmean_CovM_SummaryDev.xlsx"), sheet = "Sheet1")


# Alternatively, construct excel file path
excel_file_path <- file.path(getwd(),"data","Kmean_CovM_SummaryDev.xlsx")
df <- read_xlsx(excel_file_path, sheet = "Sheet1")


# Clean and check column names
df <- clean_names(df)
column_names_all = names(df)
column_names_all

# Define response variables (RV) list
# This will be the variables be used for running analysis
rv_list <- column_names_all[c(8:17)]
rv_list 

# Define covariate variables list
cv_list <- c("sex", "layer","age")
cv_list


# Check unique values of 'location'
table(df$subject_id)
table(df$location)
#' very few location 4, remove from analysis to decrease the complexity of data fitting

# overview of the variables
skimr::skim(df, c(all_of(cv_list), all_of(rv_list)))


# Update cv_list with factor variables
cv_list[2] = 'f.layer'
cv_list [3] = 'f.age'
cv_list


# Create a lite version including rv_list and cv_list and information required for lmer 
df_lite = df %>%
  filter(location != 4) %>%
  mutate (f.age = factor(age), f.layer = factor(layer)) %>%
  select(all_of(cv_list), all_of(rv_list),"subject_id", "location")

  
# Check if all outcome variables of interest are included in the data
summary(rv_list %in% names(df_lite))
summary(cv_list %in% names(df_lite))

# Overview of the variables
simple_test <-skimr::skim(df_lite, c(all_of(cv_list), all_of(rv_list)))
simple_test


# Density estimates of response variables
lapply(rv_list, function(x) with(df_lite, densityPlot(get(x), xlab = x)))

# Check distribution of variable as need
lapply(rv_list, function(x) with(df_lite, is.numeric(get(x))))

# Three-Way analysis
# Apply to all outcome variables in rv_list

mod_list <- lapply(rv_list, function(x) lmer_three_way(x))
names(mod_list) <- rv_list
tab_model(mod_list)

# Check which variables throw warning message
sapply(rv_list, function(x) mod_list[[x]]@optinfo$conv$lme4$messages)

# model selection ---------------------------------------------------------
mod_list_select = lapply(rv_list, model_select)
names(mod_list_select) = rv_list
tab_model(mod_list_select)

sapply(rv_list, function(x) mod_list_select[[x]]@optinfo$conv$lme4$messages)

# Wald test for linear mixed-effects models 
# Type III Wald F tests with Kenward-Roger degrees of freedom
anova_list_select = lapply(mod_list_select, Anova, type = 3, 
                           test.statistic = "F")
anova_list_select

# After determining the correct model for fitting, employ bootstrap---------
mod_class = sapply(mod_list_select, function(x) class(x)[1])
mod_list_lm = mod_list_select[which(mod_class == "lm")]
mod_list_lmer = mod_list_select[which(mod_class == "lmerMod")]

# lm bootstrap
# default: R = 10000
lm_boot_list = lapply(mod_list_lm, lm_bootstrap, R = 10000)

# lmer bootstrap
Sys.time()
tic("case bootstrap")
set.seed(47408)
cl = makeCluster(10)
registerDoParallel(cl)
# default: b1 = 625, b2 = 16 --> B = 10000
case_boot_list = lapply(mod_list_lmer, case_bootstrap, b1 = 1000, b2 = 10)
stopCluster(cl)
toc()

# Save results
save.image(file = here("output","Kmean_CovM_SummaryDev.RData"))
