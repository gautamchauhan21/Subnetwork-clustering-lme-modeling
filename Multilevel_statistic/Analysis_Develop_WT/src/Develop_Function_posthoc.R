
# ---- Load and install required libraries ----

# List of required packages
required_packages <- c(
  "dplyr",         # Data manipulation
  "pbkrtest",      # Parametric bootstrap and Kenward-Roger approximation
  "lme4",          # Linear and generalized linear mixed-effects models
  "lmeresampler",  # Resampling methods for mixed models
  "emmeans",       # Estimated marginal means (least-squares means)
  "parameters",    # Model parameter extraction and summary
  "bayestestR",    # Bayesian model diagnostics and summaries
  "ggplot2",       # Data visualization
  "ggbeeswarm",    # Quasirandom plots (e.g., beeswarm-style jittering)
  "knitr"          # Report generation
)

# Install any packages that are not already installed
installed_packages <- rownames(installed.packages())
for (pkg in required_packages) {
  if (!(pkg %in% installed_packages)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

# Load all packages
lapply(required_packages, library, character.only = TRUE)

# Function to extract replicates and assign attributes and class
boot_model_extract = function(x) {
  # Extract lmeresamp objects
  mod = case_boot_list[[x]]
  #mod = residual_boot_list[[x]]
 
   # Save replicates
  mod_rep = mod$replicates
  
  # Keep complete cases
  mod_rep = mod_rep[complete.cases(mod_rep), ]
  
  # Assign attributes and class to the replicates
  attr(mod_rep, "original_model") = mod$model
  class(mod_rep) = append(class(mod_rep), "bootstrap_model")
  class(mod_rep) = append(class(mod_rep), "see_bootstrap_model")
  
  # Return replicates
  return(mod_rep)
}


# Function to extract replicates and assign attributes and class
# Returns two data frames for plotting: (1) boot sample and (2) original data
boot_model_extract_raw = function(x) {
  # Extract lmeresamp objects
  mod = case_boot_list[[x]]
 
  # Save replicates
  mod_rep = mod$replicates
  
  # Save original data
  mod_data = mod$model@frame
  
  # Keep complete cases
  mod_rep = mod_rep[complete.cases(mod_rep), ]
  
  # Assign attributes and class to the replicates
  attr(mod_rep, "original_model") = mod$model
  class(mod_rep) = append(class(mod_rep), "bootstrap_model")
  class(mod_rep) = append(class(mod_rep), "see_bootstrap_model")
  
  # Return replicates
  return(list(mod_rep, mod_data))
}


## custom contrasts --------------------------------------------------------
pw_emm = emmeans:::revpairwise.emmc(1:5)
as.matrix(names(pw_emm))
pw_emm_contr = pw_emm[, c(1, 3, 6, 10, 2, 9)]
names(pw_emm_contr) = c("age13 - age11", "age15 - age13", "age18 - age15", 
                        "age21 - age18", "age15 - age11", "age21 - age15")


## Report main effects of sex, layer and age to compare ANOVA table of before bootstrapping -------------

report_main_effect = function(x, var_name = x) {
  
  # Extract replicates
  mod_boot = boot_model_extract(x)
  
  # ANOVA of selected model-before bootstrapping
  anova = anova_list_select[[x]]
  
  # emmeans: age
  emm_age = emmeans(mod_boot, "f.age")
  #emm_contr_age = contrast(emm_age, method = pw_emm_contr)
  emm_contr_age = contrast(emm_age, method = "pairwise")
  emm_param_age = as.data.frame(model_parameters(
    emm_contr_age, centrality = c("median", "mean"), test = "pd"))
  emm_param_age$pval = pd_to_p(emm_param_age$pd)
  
  # emmeans: sex
  emm_sex = emmeans(mod_boot, "sex")
  emm_contr_sex = contrast(emm_sex, method = "pairwise")
  emm_param_sex = as.data.frame(model_parameters(
    emm_contr_sex, centrality = c("median", "mean"), test = "pd"))
  emm_param_sex$pval = pd_to_p(emm_param_sex$pd)
  
  # emmeans: layer
  emm_layer = emmeans(mod_boot, "f.layer")
  emm_contr_layer = contrast(emm_layer, method = "pairwise")
  emm_param_layer = as.data.frame(model_parameters(
    emm_contr_layer, centrality = c("median", "mean"), test = "pd"))
  emm_param_layer$pval = pd_to_p(emm_param_layer$pd)
  
  table1 = kable(anova,caption = "ANOVA table before bootstrapping", digits = 3)
  print(table1)
  
  table2 = kable(emm_param_age, caption ="Post-hoc comparison with bootstrapping output", digits = 3)
  print(table2)
  
  table3 = kable(emm_param_sex ,digits = 3)
  print(table3)
  
  table4 = kable(emm_param_layer, digits = 3)
  print(table4)
  
}



## Report age-sex-layer -----the most comprehensive comparison--------------
report_inter_age_sex_layer = function(x, var_name = x) {
  
  # Extract replicates
  mod_boot = boot_model_extract(x)
  
  # Simple comparisons
  emm = emmeans(mod_boot, ~ sex* f.layer*f.age)
  
  # age
  emm_contr_age = contrast(emm, method = "pairwise", simple = "f.age")
  #emm_contr_age = contrast(emm, method = pw_emm_contr, simple = "f.age")
  emm_param_age = as.data.frame(model_parameters(
    emm_contr_age, centrality = c("median", "mean"), test = "pd"))
  emm_param_age$pval = pd_to_p(emm_param_age$pd)
  emm_param_age$Median.1 = NULL
  
  # sex
  emm_contr_sex = contrast(emm, method = "consec", simple = "sex")
  emm_param_sex = as.data.frame(model_parameters(
    emm_contr_sex, centrality = c("median", "mean"), test = "pd"))
  emm_param_sex$pval = pd_to_p(emm_param_sex$pd)
  emm_param_sex$Median.1 = NULL
  
  # layer
  emm_contr_layer = contrast(emm, method = "consec", simple = "f.layer")
  emm_param_layer = as.data.frame(model_parameters(
    emm_contr_layer, centrality = c("median", "mean"), test = "pd"))
  emm_param_layer$pval = pd_to_p(emm_param_layer$pd)
  emm_param_layer$Median.1 = NULL
  
  # print formatted tables in sequence
  table1 = kable (emm_param_age,
                  caption = "Post-hoc comparison with bootstrapping output", digits = 3)
  print(table1)
  table2 = kable (emm_param_sex,digits = 3)
  print(table2)
  table3 = kable (emm_param_layer,digits = 3)
  print(table3)
}


##' Plot raw data and bootstrapping results, faceted by cortical layer and sex.
#'
#' This function generates a plot showing raw data points overlaid with
#' bootstrapping results (e.g., from a linear mixed-effects model) for a
#' specified variable, faceted by cortical layer and with different colors
#' representing sex. It provides options for customizing the y-axis breaks
#' and limits, legend position, and plot title. If `y_breaks` or `y_limits`
#' are not provided (i.e., remain as `NULL`), the y-axis scale will be
#' automatically determined by `ggplot2`. By default, the plot title will
#' be the name of the variable being plotted, but this can be overridden
#' with a custom title or removed entirely.
#'
#' @param var_name A string specifying the name of the variable to plot. This
#'   should correspond to a column in the data frame used by the
#'   `boot_model_extract_raw` function. The variable name is also used as the
#'   default plot title.
#' @param y_breaks A numeric vector specifying the breaks for the y-axis. If
#'   `NULL` (default), breaks are automatically determined.
#' @param y_limits A numeric vector of length 2 specifying the limits for the
#'   y-axis (e.g., `c(0, 1)`). If `NULL` (default), limits are automatically
#'   determined.
#' @param legend_position A numeric vector of length 2 specifying the x and y
#'   coordinates for the legend's position within the plot (values between 0
#'   and 1). Defaults to `c(0.85, 0.85)`. Alternatively, common legend
#'   positions like "top", "bottom", "left", "right" can also be used as
#'   strings.
#' @param fig_title A string specifying the title of the plot. Defaults to the
#'   value of `var_name`. To remove the title, set this argument to `NULL` or
#'   an empty string (`""`).
#'
#' @return A `ggplot` object representing the generated plot.
#'
#' @examples
#' \dontrun{
#' # Assuming you have a data frame and a function 'boot_model_extract_raw'
#' # defined elsewhere.
#'
#' # Plot 'expression_level' with default title and automatic y-axis
#' plot1 <- plot_raw_boot_facet_layer("expression_level")
#' print(plot1)
#'
#' # Plot 'firing_rate' with a custom title and specific y-axis breaks
#' plot2 <- plot_raw_boot_facet_layer("firing_rate",
#'                                     y_breaks = c(0, 0.5, 1.0),
#'                                     fig_title = "Neuronal Firing Rate")
#' print(plot2)
#'
#' # Plot 'protein_concentration' with specified y-limits and legend position
#' plot3 <- plot_raw_boot_facet_layer("protein_concentration",
#'                                     y_limits = c(0, 2),
#'                                     legend_position = "bottom")
#' print(plot3)
#'
#' # Plot 'gene_activity' with no title and automatic y-axis
#' plot4 <- plot_raw_boot_facet_layer("gene_activity",
#'                                     fig_title = "")
#' print(plot4)
#' }

plot_raw_boot_facet_layer = function(var_name,
                                     y_breaks = NULL,
                                     y_limits = NULL,
                                     legend_position = c(0.85, 0.85),
                                     fig_title = var_name) {
  
  # Define label for cortical layer
  layer.label = c("Layer 2/3", "Layer 4")
  names(layer.label) = c("2", "4")
  
  # Extract replicates using pre-defined function
  mod_boot = boot_model_extract_raw(var_name)[[1]]
  
  # Prepare data for plotting
  mod_data = boot_model_extract_raw(var_name)[[2]]
 
  mod_data$yvar = mod_data[, var_name]
  mod_data$xvar = mod_data$f.age
  mod_data$tvar = mod_data$sex
  
  # Start the plot
  plot = emmip(mod_boot, sex ~ f.age | f.layer, CIs = TRUE,
               CIarg = list(lwd = 1, alpha = 1),
               dodge = 0.5, linewidth = 0.5, fatten = 1) +
    geom_beeswarm(data = mod_data, dodge.width = 0.5, alpha = 0.7,
                  shape = 1, size = 1) +
    scale_x_discrete(limits = c("11", "13", "15", "18", "21"), labels = c("P11", "P13", "P15", "P18", "P21")) +
    facet_wrap(~ f.layer, labeller = labeller(f.layer = layer.label), strip.position = "bottom") +
    labs(x = "", y = "", color = "Sex") +
    theme(strip.text = element_text(size = 10),
          strip.background = element_blank(),
          strip.placement = ("outside"),
          legend.text = element_text(color = "black", size  = "10"),
          legend.position = legend_position,
          legend.title = element_blank(),
          axis.text.x = element_text(color = "black", size = "10"),
          axis.line.x.bottom = element_line (linewidth = 1/2),
          axis.text.y = element_text(color = "black", size = "10"),
          axis.line.y = element_line(linewidth = 1/2),
          axis.ticks.x = element_blank(),
          plot.title = element_text(size = 12, hjust = 0.5), # centers the title
          plot.background = element_blank(),
          plot.margin = unit(c(0, 0.5, 0, 0), 'cm')) # top, right, bottom, left
  
  # Conditionally set y-axis breaks and limits
  if (!is.null(y_breaks) && !is.null(y_limits)) {
    plot = plot + scale_y_continuous(breaks = y_breaks, limits = y_limits)
  } else {
    plot = plot + scale_y_continuous() # Use automatic breaks and limits
  }
  
  # Conditionally add the title
  if (!is.null(fig_title) && !identical(fig_title, "")) {
    plot = plot + ggtitle(fig_title)
  }
  
  return(plot)
}



##plot raw data and bootstrapping result -Facet by sex ------------------
# the code here act similar to the one prior to Facet by cortical layer
plot_raw_boot_facet_sex = function(var_name,
                                     y_breaks = NULL,
                                     y_limits = NULL,
                                     legend_position = c(0.4, 0.85),
                                     fig_title = var_name) {
  
  # Define label for cortical layer
  sex.label = c("Female", "Male")
  names(sex.label) = c("F", "M")
  
  # Extract replicates using pre-defined function
  mod_boot = boot_model_extract_raw(var_name)[[1]]
  
  # Prepare data for plotting
  mod_data = boot_model_extract_raw(var_name)[[2]]
  
  mod_data$yvar = mod_data[, var_name]
  mod_data$xvar = mod_data$f.age
  mod_data$tvar = mod_data$f.layer
  
  # Start the plot
  plot = emmip(mod_boot, f.layer ~ f.age | sex, CIs = TRUE,
               CIarg = list(lwd = 1, alpha = 1),
               dodge = 0.5, linewidth = 0.5, fatten = 1) +
    geom_beeswarm(data = mod_data, dodge.width = 0.5, alpha = 0.7,
                  shape = 1, size = 1) +
    scale_x_discrete(limits = c("11", "13", "15", "18", "21"), labels = c("P11", "P13", "P15", "P18", "P21")) +
    facet_wrap(~ sex, labeller = labeller(sex = sex.label), strip.position = "bottom") +
    labs(x = "", y = "", color = "Layer") +
    scale_color_jama()+
    theme(strip.text = element_text(size = 10),
          strip.background = element_blank(),
          strip.placement = ("outside"),
          legend.text = element_text(color = "black", size  = "10"),
          legend.position = legend_position,
          legend.title = element_blank(),
          axis.text.x = element_text(color = "black", size = "10"),
          axis.line.x.bottom = element_line (linewidth = 1/2),
          axis.text.y = element_text(color = "black", size = "10"),
          axis.line.y = element_line(linewidth = 1/2),
          axis.ticks.x = element_blank(),
          plot.title = element_text(size = 12, hjust = 0.5), # centers the title
          plot.background = element_blank(),
          plot.margin = unit(c(0, 0.5, 0, 0), 'cm'))
          
  # Conditionally set y-axis breaks and limits
  if (!is.null(y_breaks) && !is.null(y_limits)) {
    plot = plot + scale_y_continuous(breaks = y_breaks, limits = y_limits)
  } else {
    plot = plot + scale_y_continuous() # Use automatic breaks and limits
  }
  
  # Conditionally add the title
  if (!is.null(fig_title) && !identical(fig_title, "")) {
    plot = plot + ggtitle(fig_title)
  }
  
  return(plot)
}



plot_raw_boot_merge_sex = function(var_name,
                                   y_breaks = NULL,
                                   y_limits = NULL,
                                   legend_position = c(0.85, 0.85),
                                   fig_title = var_name) {
 
  # Extract replicates using pre-defined function
  mod_boot = boot_model_extract_raw(var_name)[[1]]
  
  # Prepare data for plotting
  mod_data = boot_model_extract_raw(var_name)[[2]]
  
  mod_data$yvar = mod_data[, var_name]
  mod_data$xvar = mod_data$f.age
  mod_data$tvar = mod_data$f.layer
  
  # Start the plot
  plot = emmip(mod_boot, f.layer ~ f.age, CIs = TRUE,
               CIarg = list(lwd = 1, alpha = 1),
               dodge = 0.5) +
    geom_beeswarm(data = mod_data, dodge.width = 0.5, alpha = 0.7,
                  shape = 1, size = 1) +
    scale_x_discrete(limits = c("11", "13", "15", "18", "21"), labels = c("P11", "P13", "P15", "P18", "P21")) +
    labs(x = "", y = "", color = "Layer") +
    scale_color_jama() +
    theme(legend.text = element_text(color = "black", size  = "10"),
          legend.position = legend_position,
          legend.title = element_blank(),
          axis.text.x = element_text(color = "black", size = "10"),
          axis.line.x.bottom = element_line (linewidth = 1/2),
          axis.text.y = element_text(color = "black", size = "10"),
          axis.line.y = element_line(linewidth = 1/2),
          axis.ticks.x = element_blank(),
          plot.title = element_text(size = 12, hjust = 0.5), # centers the title
          plot.background = element_blank(),
          plot.margin = unit(c(0, 0.5, 0, 0), 'cm'))# top, right, bottom, left
  
  # Conditionally set y-axis breaks and limits
  if (!is.null(y_breaks) && !is.null(y_limits)) {
    plot = plot + scale_y_continuous(breaks = y_breaks, limits = y_limits)
  } else {
    plot = plot + scale_y_continuous() # Use automatic breaks and limits
  }
          
  
  # Conditionally add the title
  if (!is.null(fig_title) && !identical(fig_title, "")) {
    plot = plot + ggtitle(fig_title)
  }
  
  return(plot)
}



