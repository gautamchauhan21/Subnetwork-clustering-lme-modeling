# Load necessary libraries
library(dplyr)
library(pbkrtest)      # Provides tools for parametric bootstrap and Kenward-Roger approximation
library(lme4)          # Fits linear and generalized linear mixed-effects models
library(lmeresampler)  # Provides functions for resampling mixed models
library(emmeans)       # Estimated marginal means, aka least-squares means
library(parameters)
library(bayestestR)
library(ggplot2)       # Creating Plot
library(ggbeeswarm)    # Jittering points within categories
library(knitr)         # A general-purpose tool for dynamic report generation in R

## Define function
# Extract lm bootstrapped coefficients ------------------------------------

lm_boot_extract = function(x) {
  
  # extract simpleboot object
  mod = lm_boot_list[[x]]
  
  # extract bootstrapped coefs
  boot_ls = lapply(1:length(mod$boot.list), 
                   function(x) mod$boot.list[[x]]$coef)
  boot_rep = as.data.frame(do.call(rbind, boot_ls))
  
  # keep complete cases
  boot_rep = boot_rep[complete.cases(boot_rep), ]
  
  # assign attributes and class to the replicates
  attr(boot_rep, "original_model") = mod$orig.lm
  class(boot_rep) = append(class(boot_rep), "bootstrap_model")
  class(boot_rep) = append(class(boot_rep), "see_bootstrap_model")
  
  # return replicates
  return(boot_rep)
}


# extract the replicates and assign attributes and class ------------------

lmer_boot_extract = function(x) {
  # extract lmeresamp objects
  mod = case_boot_list[[x]]
  
  # save replicates
  boot_rep = mod$replicates
  
  # keep complete cases
  boot_rep = boot_rep[complete.cases(boot_rep), ]
  
  # assign attributes and class to the replicates
  attr(boot_rep, "original_model") = mod$model
  class(boot_rep) = append(class(boot_rep), "bootstrap_model")
  class(boot_rep) = append(class(boot_rep), "see_bootstrap_model")
  
  # return replicates
  return(boot_rep)
}



# Report main effect ------------------------------------------------------
emm_var_main = function(var_name) {
  emm = emmeans(mod_boot, var_name)
  emm_contr = contrast(emm, method = "consec")
  emm_param = as.data.frame(model_parameters(
    emm_contr, centrality = c("median", "mean"), test = "pd"))
  emm_param$pval = pd_to_p(emm_param$pd)
  return(emm_param)
}


report_main = function(x, var_name = x) {
  if (x %in% names(case_boot_list)) {
    mod_boot = lmer_boot_extract(x)
  } else {
    mod_boot = lm_boot_extract(x)
  }
  
  mod_boot <<- mod_boot
  
  #emm_list = lapply(cv_list[c(1:2, 4)], emm_var_main)
  #names(emm_list) = cv_list[c(1:2,4)]
  
  emm_list = lapply(cv_list, emm_var_main)
  names(emm_list) = cv_list
  
  
  anova = anova_list_select[[x]]
  
  table1 = kable(anova,caption = "ANOVA table before bootstrapping", digits = 3)
  print(table1)
  
  table2 = kable(emm_list$sex,caption = "Post-hoc comparison with bootstrapping output", digits = 3)
  print(table2)
  
  table3 = kable(emm_list$f.layer, digits = 3)
  print(table3)
  
  table4 = kable(emm_list$p15_cno, digits = 3)
  print(table4)
  
}


# Report interaction
report_all_interaction_simple = function(x, var_int1 = "p15_cno", var_int2 = "sex", var_int3 = "f.layer") {
  # extract bootstrapped coef by model class
  if (x %in% names(case_boot_list)) {
    mod_boot = lmer_boot_extract(x)
  } else {
    mod_boot = lm_boot_extract(x)
  }
  
  f = paste(var_int1, "*", var_int2, "*", var_int3)
  emm = emmeans(mod_boot, eval(bquote(~ .(f))))
  emm_contr = contrast(emm, method = "consec")
  emm_param = as.data.frame(model_parameters(
    emm_contr, centrality = c("median", "mean"), test = "pd"))
  emm_param$pval = pd_to_p(emm_param$pd)
  
  table = kable(emm_param, caption = "Post-hoc comparison with bootstrapping output", digits = 3)
  print(table)
}



# Function for exacting data for plotting figure
boot_model_extract_raw = function(x) {
  # extract replicates
  if (x %in% names(case_boot_list)) {
    boot_rep = lmer_boot_extract(x)
  } else {
    boot_rep = lm_boot_extract(x)
  }
  
  return(boot_rep)
}


#'
#' This function takes the name of a response variable, extracts bootstrap
#' replicates from a fitted model, and generates a faceted plot showing
#' individual subject data overlaid with marginal effect estimates and their
#' confidence intervals for different cortical layers.
#'
#' @param x A character string specifying the name of the response variable
#'   present in the `df_lite` data frame and used in the fitted model. This
#'   variable is assumed to be a column in `df_lite`.
#' @param y_breaks A numeric vector specifying the breaks for the y-axis. If
#'   `NULL` (default), breaks are automatically determined.
#' @param y_limits A numeric vector of length 2 specifying the limits of the
#'   y-axis. If `NULL` (default), limits are automatically determined.
#' @param legend_position Specifies the position of the legend. It can be a
#'   numeric vector of length 2 (x and y coordinates within the plot area,
#'   ranging from 0 to 1) or a character string indicating a standard
#'   position ("top", "bottom", "left", "right"). Defaults to `c(0.4, 0.85)`.
#' @param fig_title A character string specifying the title of the plot. If
#'   not provided or set to `NULL` or `""`, no title will be added. Defaults
#'   to the value of `x`.
#'
#' @details
#' This function assumes the existence of a data frame named `df_lite` and a
#' fitted model from which bootstrap replicates can be extracted using a
#' function called `boot_model_extract_raw()`. The `df_lite` data frame is
#' expected to have columns named "subject_id", "location", "sex", and a
#' column with the name specified by the `x` argument (the response variable).
#' It also expects a factor column named "f.layer" and a factor column named
#' "p15_cno".
#'
#' The plot displays individual subject data as faint lines, colored by sex,
#' with a combined identifier of subject ID, location, and sex used for grouping.
#' Marginal effects (estimated using `emmip` from the `emmeans` package) and
#' their confidence intervals are shown as points with error bars. These marginal
#' effects represent the interaction between "p15_cno" (levels "-" and "+") and
#' "f.layer". The plot is faceted by the "f.layer" variable, with custom labels
#' "Layer 2/3" and "Layer 4".
#'
#' @return A `ggplot` object representing the generated plot.
#'
#' @import ggplot2
#' @import dplyr
#' @import emmeans
#' @import tidyr
#' @import beeswarm
#'
#' @examples
#' \dontrun{
#' # Assuming 'my_response' is a column in 'df_lite' and
#' # 'my_model_boot' is a fitted model with bootstrap replicates.
#'
#' # Example usage with default settings:
#' plot_raw_boot_facet_layer("my_response")
#'
#' # Example with custom y-axis breaks and limits:
#' plot_raw_boot_facet_layer("my_response",
#'                           y_breaks = seq(0, 0.8, by = 0.1),
#'                           y_limits = c(0.2, 0.7))
#'
#' # Example with a legend at the bottom and a custom title:
#' plot_raw_boot_facet_layer("my_response",
#'                           legend_position = "bottom",
#'                           fig_title = "Effect of Treatment by Layer")
#'
#' # Example with automatic y-axis and no title:
#' plot_raw_boot_facet_layer("my_response",
#'                           y_breaks = NULL,
#'                           y_limits = NULL,
#'                           fig_title = "")
#' }
#'

plot_raw_boot_facet_layer = function(x,
                                      y_breaks = NULL, 
                                      y_limits = NULL, 
                                      legend_position = c(0.4, 0.85),
                                      fig_title = x) {
  # Define label
  layer.label = c("Layer 2/3", "Layer 4")
  names(layer.label) = c("2", "4")
  
  sex.label = c("Female", "Male")
  names(sex.label) = c("F", "M")
  
  
  # Ensure the response variable exists in df_lite
  if (!x %in% names(df_lite)) {
    stop(paste("The variable '", x, "' is not found in the 'df_lite' data frame."))
  }
  
  
  # Extract replicates
  mod_boot = boot_model_extract_raw(x)
  
  df_lite$yvar = df_lite[, x, drop = TRUE]
  
  # Create a new variable to denote the combination of subject id, layer, and
  # location as an unit for plotting individual lines
  df_lite = df_lite %>%
    unite(col = "set", c("subject_id", "location", "sex"), remove = FALSE)
  
  # Calculate marginal effect point estimates and confidence intervals
  # p15_cno' and 'f.layer' are relevant predictors in the model
  emms = emmip(mod_boot, ~ p15_cno | f.layer, CIs = TRUE, plotit = FALSE)
  
 p<- emms %>%
    ggplot(aes(x = factor(p15_cno, levels = c("neg", "CNO")), y = yvar)) +
    geom_line(data = df_lite, 
              aes(group = set, color = sex), alpha = 0.5, linewidth = 0.5) + 
    geom_pointrange(aes(ymin = LCL, ymax = UCL),linewidth = 1, size = 0.2) + 
    geom_line(aes(group = f.layer), linewidth = 1) +
    scale_x_discrete(labels = c("-", "+")) +
    facet_wrap(~ f.layer, labeller = labeller(f.layer = layer.label), 
               strip.position = "bottom") +
    labs(x = "", y = "") +
    theme(strip.text = element_text(size = 10),
          strip.background = element_blank(),
          strip.placement = ("outside"),
          legend.text = element_text(color = "black", size  = "10"),
          legend.position = legend_position,
          legend.title = element_blank(),
          panel.spacing = unit(0, "lines"),
          axis.text.x = element_text(color = "black", size = "10"),
          axis.line.x.bottom = element_line (linewidth = 1/2),
          axis.text.y = element_text(color = "black", size = "10"),
          axis.line.y = element_line(linewidth = 1/2), 
          axis.ticks.x = element_blank(),
          plot.title = element_text(size=12, hjust = 0.5),
          plot.background = element_blank(),
          plot.margin = unit(c(0,0.5,1.5,0), 'cm'))
 
 # Conditionally set y-axis breaks and limits
 if (!is.null(y_breaks) && !is.null(y_limits)) {
   plot = plot + scale_y_continuous(breaks = y_breaks, limits = y_limits)
 } else {
   plot = plot + scale_y_continuous() # Use automatic breaks and limits
 }
 
 # Conditionally add the title
   if (!is.null(fig_title) && !identical(fig_title, "")) {
     p = p + ggtitle(fig_title)
   }
   
  return(p)
}



#' Plot Raw Bootstrap Results with Facets for Sex
plot_raw_boot_facet_sex = function(x,
                                   y_breaks = NULL,
                                   y_limits = NULL,
                                   legend_position = "none",
                                   fig_title = x) {
  # Define labels for the layer and sex facets
  layer.label = c("Layer 2/3", "Layer 4")
  names(layer.label) = c("2", "4")
  
  sex.label = c("Female", "Male")
  names(sex.label) = c("F", "M")
  
  # Extract bootstrap replicates using the assumed function
  mod_boot = boot_model_extract_raw(x)
  
  # Ensure the response variable exists in df_lite
  if (!x %in% names(df_lite)) {
    stop(paste("The variable '", x, "' is not found in the 'df_lite' data frame."))
  }
  df_lite$yvar = df_lite[, x, drop = TRUE]
  
  # Create a new variable to denote the combination of subject id, layer, and
  # location as a unit for plotting individual lines
  df_lite = df_lite %>%
    unite(col = "set", c("subject_id", "location", "f.layer"), remove = FALSE)
  
  # Calculate marginal effect point estimates and confidence intervals
  # Assuming 'p15_cno' and 'sex' are relevant predictors in the model
  emms = emmip(mod_boot, ~ p15_cno | sex, CIs = TRUE, plotit = FALSE)
  
  p <- emms %>%
    ggplot(aes(x = factor(p15_cno, levels = c("neg", "CNO")), y = yvar)) +
    # Add individual subject lines
    geom_line(data = df_lite,
              aes(group = set, color = f.layer), alpha = 0.5, linewidth = 0.5) +
    # Add marginal effect point estimates with confidence intervals
    geom_pointrange(aes(ymin = LCL, ymax = UCL), linewidth = 1, size = 0.2) +
    # Add lines connecting the marginal effect estimates within each facet
    geom_line(aes(group = sex), linewidth = 1) +
    # Customize the x-axis labels
    scale_x_discrete(labels = c("-", "+")) +
    # Facet the plot by 'sex' with custom labels and positioning
    facet_wrap(~ sex, labeller = labeller(sex = sex.label),
               strip.position = "bottom") +
    # Set axis labels
    labs(x = "", y = "") +
    # Customize the theme
    theme(strip.text = element_text(size = 10),
          strip.background = element_blank(),
          strip.placement = ("outside"),
          legend.text = element_text(color = "black", size  = "10"),
          legend.position = legend_position,
          legend.title = element_blank(),
          panel.spacing = unit(0, "lines"),
          axis.text.x = element_text(color = "black", size = "10"),
          axis.line.x.bottom = element_line (linewidth = 1/2),
          axis.text.y = element_text(color = "black", size = "10"),
          axis.line.y = element_line(linewidth = 1/2),
          axis.ticks.x = element_blank(),
          plot.title = element_text(size = 12, hjust = 0.5),
          plot.background = element_blank(),
          plot.margin = unit(c(0,0.5,1.5,0), 'cm'))
 
   # Conditionally set y-axis breaks and limits
  if (!is.null(y_breaks) && !is.null(y_limits)) {
    plot = plot + scale_y_continuous(breaks = y_breaks, limits = y_limits)
  } else {
    plot = plot + scale_y_continuous() # Use automatic breaks and limits
  }
  
  # Conditionally add the title
  if (!is.null(fig_title) && !identical(fig_title, "")) {
    p = p + ggtitle(fig_title)
  }
  
  return(p)
}



#' Plot Raw Bootstrap Results 

plot_raw_boot = function(x,
                                  y_breaks = NULL,
                                  y_limits = NULL,
                                  legend_position = "none",
                                  fig_title = x) {
  # Extract bootstrap replicates using the assumed function
  mod_boot = boot_model_extract_raw(x)
  
  # Ensure the response variable exists in df_lite
  if (!x %in% names(df_lite)) {
    stop(paste("The variable '", x, "' is not found in the 'df_lite' data frame."))
  }
  df_lite$yvar = df_lite[, x, drop = TRUE]
  
  # Create a new variable to denote the combination of subject id, layer, and
  # location as a unit for plotting individual lines, color by sex directly
  df_lite = df_lite %>%
    unite(col = "set", c("subject_id", "location", "f.layer"), remove = FALSE)
  
  # Calculate marginal effect point estimates and confidence intervals,
  # aggregated across layers (no faceting)
  emms = emmip(mod_boot, ~ p15_cno, CIs = TRUE, plotit = FALSE)
  
  p <- emms %>%
    ggplot(aes(x = factor(p15_cno, levels = c("neg", "CNO")), y = yvar)) +
    # Add individual subject lines, colored by sex
    geom_line(data = df_lite,
              aes(group = set, color = sex), alpha = 0.5, linewidth = 0.5) +
    # Add marginal effect point estimates with confidence intervals
    geom_pointrange(aes(ymin = LCL, ymax = UCL), linewidth = 1, size = 0.2) +
    # Add lines connecting the marginal effect estimates
    geom_line(aes(group = 1), linewidth = 1) + # Group by 1 for a single line
    # Customize the x-axis labels
    scale_x_discrete(labels = c("-", "+")) +
    # Set axis labels
    labs(x = "", y = "") +
    # Customize the theme
    theme(strip.text = element_text(size = 10),
          strip.background = element_blank(),
          strip.placement = ("outside"),
          legend.text = element_text(color = "black", size  = "10"),
          legend.position = legend_position,
          legend.title = element_blank(),
          panel.spacing = unit(0, "lines"),
          axis.text.x = element_text(color = "black", size = "10"),
          axis.line.x.bottom = element_line (linewidth = 1/2),
          axis.text.y = element_text(color = "black", size = "10"),
          axis.line.y = element_line(linewidth = 1/2),
          axis.ticks.x = element_blank(),
          plot.title = element_text(size = 12, hjust = 0.5),
          plot.background = element_blank(),
          plot.margin = unit(c(0,0.5,1.5,0), 'cm'))
  
  # Conditionally set y-axis breaks and limits
  if (!is.null(y_breaks) && !is.null(y_limits)) {
    plot = plot + scale_y_continuous(breaks = y_breaks, limits = y_limits)
  } else {
    plot = plot + scale_y_continuous() # Use automatic breaks and limits
  }
   
  # Conditionally add the title
  if (!is.null(fig_title) && !identical(fig_title, "")) {
    p = p + ggtitle(fig_title)
  }
  
  return(p)
}








