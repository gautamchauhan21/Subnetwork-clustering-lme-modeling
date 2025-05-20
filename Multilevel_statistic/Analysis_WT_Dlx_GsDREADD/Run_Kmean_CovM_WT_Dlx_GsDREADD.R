# Refresh the R session to clear the environment
# This ensures that any previous objects, functions, or loaded packages are removed.
freshr::freshr()

# Verify the working directory
print(getwd()) # Should print the path to your working directory

# Load necessary libraries
# These libraries provide various functions and tools needed for data manipulation, analysis, and visualization.
library(readxl)
library(janitor)
library(dplyr)
library(car)
library(sjPlot)
library(tictoc)
library(skimr)   
library(stringr)
library(here)

# Load source file-custom functions
source(here("src", "WT_Dlx_GsDREADD_Function.R"))

# Import data from Excel file
df <- read_xlsx(here("data","Kmean_CovM_SummaryDlxGsDREADD.xlsx"), sheet = "Sheet1")

# Clean and check column names
df = clean_names(df)
column_names_all = names(df)
column_names_all

# rename the label
df <- df %>% rename(p15_cno = p15cno)

# Define covariate variables list
cv_list= c("sex", "layer", "p15_cno")
cv_list

# Define response variables (RV) list
# This will be the variables be used for running analysis
rv_list = column_names_all[c(10:19)]

# Check unique values of 'location'
table(df$location)

# Check if all variables are included in the data
summary(rv_list %in% names(df))
summary(cv_list %in% names(df))

# Overview of the variables
skimr::skim(df, c(all_of(cv_list), all_of(rv_list)))
# the one row is empty, need to remove it when generating the df_lite

# Check unique values of animal_id
table(df$subject_id)

# create a lite version and convert the list of response variables to numeric
df_lite = df %>%
  filter(!is.na(n_cls_before_stat))%>%
  select(all_of(cv_list), all_of(rv_list),subject_id,location) %>%
  mutate(f.layer = factor(layer))

# Update cv_list with factor variables
cv_list[2] = "f.layer"

# Density estimates of response variables
lapply(rv_list, function(x) with(df_lite, densityPlot(get(x), xlab = x)))

# Check distribution of RV
summary(df_lite$silhs_mean_before_stat)
summary(df_lite$n_cls)

# Overview of the variables and save the result
simple_test <-skimr::skim(df_lite, c(all_of(cv_list), all_of(rv_list)))

# Three-Way analysis
# Apply to all outcome variables of interest
mod_list = lapply(rv_list, function(x) lmer_three_way(x))
names(mod_list) = rv_list
tab_model(mod_list)

# Check which variables throw the warning message
sapply(rv_list, function(x) mod_list[[x]]@optinfo$conv$lme4$messages)


# Model selection ---------------------------------------------------------
mod_list_select = lapply(rv_list, model_select)
names(mod_list_select) = rv_list
tab_model(mod_list_select)

anova_list_select = lapply(mod_list_select, Anova, type = 3, 
                           test.statistic = "F")
anova_list_select


# After determining the correct model for fitting, employ bootstrap---------
mod_class = sapply(mod_list_select, function(x) class(x)[1])
mod_list_lm = mod_list_select[which(mod_class == "lm")]
mod_list_lmer = mod_list_select[which(mod_class == "lmerMod")]


# lm bootstrap
tic("lm Bootstrap Execution") 
lm_boot_list <- lapply(mod_list_lm, lm_bootstrap, R = 10000)
toc()


# lmer bootstrap
Sys.time()
tic("case bootstrap")
set.seed(47408)
cl = makeCluster(16)
registerDoParallel(cl)
# default: b1 = 625, b2 = 16 --> B = 10000
case_boot_list = lapply(mod_list_lmer, case_bootstrap, b1 = 625, b2 = 16)
stopCluster(cl)
toc()

# Save output ---------------------------------------------------
save.image(file = here("output","Kmean_CovM_SummaryWT_Dlx_GsDREADD.RData"))
