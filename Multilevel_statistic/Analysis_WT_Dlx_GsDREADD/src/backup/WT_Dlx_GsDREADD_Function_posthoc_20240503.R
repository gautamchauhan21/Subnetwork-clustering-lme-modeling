###for report----------------------------------------------------
library(dplyr)
library(emmeans)
library(parameters)
library(bayestestR)
library(ggplot2)
library(knitr)

# extract lm bootstrapped coefficients ------------------------------------

##defince function
# extract lm bootstrapped coefficients ------------------------------------

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


boot_model_extract_raw = function(x) {
  # extract replicates
  if (x %in% names(case_boot_list)) {
    boot_rep = lmer_boot_extract(x)
  } else {
    boot_rep = lm_boot_extract(x)
  }
  
  return(boot_rep)
}


plot_raw_boot_facet_layer = function(x, var_name = x,title = title,
                                   y_breaks = c(1:6 * 0.2), 
                                   y_limits = c(0, 1), 
                                   legend_position = "none") {
  # define label
  layer.label = c("Layer 2/3", "Layer 4")
  names(layer.label) = c("2", "4")
  
  sex.label = c("Female", "Male")
  names(sex.label) = c("F", "M")
  
  # extract replicates
  mod_boot = boot_model_extract_raw(x)
  
  # use original data set instead
  # https://tibble.tidyverse.org/reference/subsetting.html
  # the subsetting behavior is somewhat different in tibbles lol
  df_lite$yvar = df_lite[, x, drop = TRUE]
  
  # create a new variable to denote the combination of subject id, layer, and
  # location as an unit
  df_lite = df_lite %>%
    unite(col = "set", c("subject_id", "location", "sex"), remove = FALSE)
  
  # marginal effect point estimates and CI 
  emms = emmip(mod_boot, ~ p15_cno | f.layer, CIs = TRUE, plotit = FALSE)
  
  emms %>%
    ggplot(aes(x = factor(p15_cno, levels = c("neg", "CNO")), y = yvar)) +
    geom_line(data = df_lite, 
              aes(group = set, color = sex), alpha = 0.5, linewidth = 0.5) + 
    geom_pointrange(aes(ymin = LCL, ymax = UCL),linewidth = 1, size = 0.2) + 
    geom_line(aes(group = f.layer), linewidth = 1) +
    scale_x_discrete(labels = c("-", "+")) +
    scale_y_continuous(breaks = y_breaks, limits = y_limits) +
    facet_wrap(~ f.layer, labeller = labeller(f.layer = layer.label), 
               strip.position = "bottom") +
    ggtitle(title) +
    labs(x = "", y = "") +
    theme(legend.position = legend_position,
          legend.title = element_blank(),
          legend.text = element_text(color = "black", size  = "10"),
          panel.spacing = unit(0, "lines"),
          strip.background = element_blank(),
          strip.placement = ("outside"),
          strip.text = element_text(size = 10),
          axis.text.x = element_text(color = "black", size = "10"),
          axis.line.x.bottom = element_line (linewidth = 1/2),
          axis.text.y = element_text(color = "black", size = "10"),
          axis.line.y = element_line(linewidth = 1/2), 
          axis.ticks.x = element_blank(),
          plot.background = element_blank(),
          plot.margin = unit(c(0,0.5,1.5,0), 'cm')) ->p
      return(p)
}


plot_raw_boot = function(x, var_name = x,title = title,
                                     y_breaks = c(1:6 * 0.2), 
                                     y_limits = c(0, 1), 
                                     legend_position = "none") {
  # define label
  layer.label = c("Layer 2/3", "Layer 4")
  names(layer.label) = c("2", "4")
  
  sex.label = c("Female", "Male")
  names(sex.label) = c("F", "M")
  
  # extract replicates
  mod_boot = boot_model_extract_raw(x)
  
  # use original data set instead
  # https://tibble.tidyverse.org/reference/subsetting.html
  # the subsetting behavior is somewhat different in tibbles lol
  df_lite$yvar = df_lite[, x, drop = TRUE]
  
  # create a new variable to denote the combination of subject id, layer, and
  # location as an unit
  df_lite = df_lite %>%
    unite(col = "set", c("subject_id", "location", "sex"), remove = FALSE)
  
  # marginal effect point estimates and CI 
  emms = emmip(mod_boot, ~ p15_cno | f.layer, CIs = TRUE, plotit = FALSE)
  
  emms %>%
    ggplot(aes(x = factor(p15_cno, levels = c("neg", "CNO")), y = yvar)) +
    geom_line(data = df_lite, 
              aes(group = set, color = sex), alpha = 0.5, linewidth = 0.5) + 
    geom_pointrange(aes(ymin = LCL, ymax = UCL),linewidth = 1, size = 0.2) + 
    geom_line(aes(group = f.layer), linewidth = 1) +
    scale_x_discrete(labels = c("-", "+")) +
    scale_y_continuous(breaks = y_breaks, limits = y_limits) +
    facet_wrap(f.layer ~ sex,ncol=4, labeller = labeller(f.layer = layer.label, sex = sex.label), 
               strip.position = "bottom") +
    ggtitle(title) +
    labs(x = "", y = "") +
    theme(legend.position = legend_position,
          legend.title = element_blank(),
          legend.text = element_text(color = "black", size  = "10"),
          panel.spacing = unit(0, "lines"),
          strip.background = element_blank(),
          strip.placement = ("outside"),
          strip.text = element_text(size = 10),
          axis.text.x = element_text(color = "black", size = "10"),
          axis.line.x.bottom = element_line (linewidth = 1/2),
          axis.text.y = element_text(color = "black", size = "10"),
          axis.line.y = element_line(linewidth = 1/2), 
          axis.ticks.x = element_blank(),
          plot.background = element_blank(),
          plot.margin = unit(c(0,0.5,1.5,0), 'cm')) ->p
  return(p)
}



plot_raw_boot_facet_sex = function(x, var_name = x,title = title,
                                   y_breaks = c(1:6 * 0.2), 
                                   y_limits = c(0, 1), 
                                   legend_position = "none") {
  # define label
  layer.label = c("Layer 2/3", "Layer 4")
  names(layer.label) = c("2", "4")
  
  sex.label = c("Female", "Male")
  names(sex.label) = c("F", "M")
  
  # extract replicates
  mod_boot = boot_model_extract_raw(x)
  
  # use original data set instead
  # https://tibble.tidyverse.org/reference/subsetting.html
  # the subsetting behavior is somewhat different in tibbles lol
  df_lite$yvar = df_lite[, x, drop = TRUE]
  
  # create a new variable to denote the combination of subject id, layer, and
  # location as an unit
  df_lite = df_lite %>%
    unite(col = "set", c("subject_id", "f.layer", "location"), remove = FALSE)
  
  # marginal effect point estimates and CI 
  emms = emmip(mod_boot, ~ p15_cno | sex, CIs = TRUE, plotit = FALSE)
  
  emms %>%
    ggplot(aes(x = factor(p15_cno, levels = c("neg", "CNO")), y = yvar)) +
    geom_line(data = df_lite, 
              aes(group = set, color = f.layer), alpha = 0.5, linewidth = 0.5) + 
    geom_pointrange(aes(ymin = LCL, ymax = UCL),linewidth = 1, size = 0.2) + 
    geom_line(aes(group = sex), linewidth = 1) +
    scale_x_discrete(labels = c("-", "+")) +
    scale_fill_discrete(labels = c("L2/3", "L4")) +
    scale_y_continuous(breaks = y_breaks, limits = y_limits) +
    facet_wrap(~ sex, labeller = labeller(sex = sex.label), 
               strip.position = "bottom") +
    ggtitle(title) +
    labs(x = "", y = "") +
    theme(legend.position = legend_position,
          legend.title = element_blank(),
          legend.text = element_text(color = "black", size  = "10"),
          panel.spacing = unit(0, "lines"),
          strip.background = element_blank(),
          #strip.background = element_rect(color = "black", linewidth = 1),
          strip.placement = ("outside"),
          strip.text = element_text(size = 10),
          axis.text.x = element_text(color = "black", size = "10"),
          axis.line.x.bottom = element_line (linewidth = 1/2),
          axis.text.y = element_text(color = "black", size = "10"),
          axis.line.y = element_line(linewidth = 1/2), 
          axis.ticks.x = element_blank(),
          plot.background = element_blank(),
          plot.margin = unit(c(0,0.5,1.5,0), 'cm')) ->p
  return(p)
}




# report main effect ------------------------------------------------------
emm_var_main = function(var) {
  emm = emmeans(mod_boot, var)
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
  
  emm_list = lapply(cv_list[c(1:2, 4)], emm_var_main)
  names(emm_list) = cv_list[c(1:2,4)]
 
  anova = anova_list_select[[x]]
  
  table1 = kable(anova,caption = "ANOVA table before bootstrapping", digits = 3)
  print(table1)
  
  table2 = kable(emm_list$sex,caption = "Main effect after bootstrapping", digits = 3)
  print(table2)
  
  table3 = kable(emm_list$f.layer, digits = 3)
  print(table3)
  
  table4 = kable(emm_list$p15_cno, digits = 3)
  print(table4)
  
}


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

