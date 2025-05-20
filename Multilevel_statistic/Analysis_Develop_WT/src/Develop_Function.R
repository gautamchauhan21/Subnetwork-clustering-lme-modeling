# ---- Load and install required libraries ----

# List of required packages
required_packages <- c(
  "lme4",          # For fitting linear mixed-effects models
  "pbkrtest",      # For Kenward-Roger approximation for mixed models
  "simpleboot",    # For simple bootstrap methods
  "lmeresampler",  # For resampling methods for mixed models
  "foreach",       # For looping constructs for parallel execution
  "doParallel"     # For parallel backend to the foreach
)

# Install any packages that are not already installed
installed_packages <- rownames(installed.packages())
for (pkg in required_packages) {
  if (!(pkg %in% installed_packages)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

# Load all required packages
lapply(required_packages, library, character.only = TRUE)


# Define MLM with three-way interaction --------------------------------------------

lmer_three_way = function(x, data = df_lite) {
  # Fit the model with three-way interaction
  mod = lmer(paste(x, "~ f.age * sex * f.layer + (1 | subject_id) + (1 | subject_id:location)"), data = data, REML = TRUE)
  return(mod)
}


# Compare lmer models --------------------------------------------------

lmer_model_comparison = function() {
  # Check if models are identical and perform Kenward-Roger comparison
  if (identical(mod0, mod_select) == TRUE) {
    krmc = KRmodcomp(mod7, mod_select) # compare to the most complex model
    print(krmc)
    if (krmc$stats$p.value < 0.05) {
      return(mod7)
    } else {
      return(mod_select)
    }
  } else if (identical(mod7, mod_select) == TRUE) {
    krmc = KRmodcomp(mod_select, mod0) # compare to the least complex model
    print(krmc)
    if (krmc$stats$p.value < 0.05) {
      return(mod_select)
    } else {
      return(mod0)
    }
  } else {
    krmc1 = KRmodcomp(mod7, mod_select) # compare to the most complex model
    krmc2 = KRmodcomp(mod_select, mod0) # compare to the least complex model
    print(krmc1)
    print(krmc2)
    if (krmc1$stats$p.value > 0.05 & krmc2$stats$p.value < 0.05) {
      return(mod_select)
    } else if (krmc1$stats$p.value < 0.05 & krmc2$stats$p.value < 0.05) {
      return(mod7)
    } else if (krmc1$stats$p.value > 0.05 & krmc2$stats$p.value > 0.05) {
      return(mod0)
    } else {
      krmc = KRmodcomp(mod7, mod0)
      if (krmc$stats$p.value < 0.05) {
        return(mod7)
      } else {
        return(mod0)
      }
    }
  }
}


# Compare lm models----------------------------------------------------

lm_model_comparison = function() {
  if (identical(mod0, mod_select) == TRUE) {
    aov = anova(mod7, mod_select) # compare to the most complex model
    print(aov)
    if (aov$`Pr(>F)`[2] < 0.05) {
      return(mod7)
    } else {
      return(mod_select)
    }
  } else if (identical(mod7, mod_select) == TRUE) {
    aov = anova(mod_select, mod0) # compare to the least complex model
    print(aov)
    if (aov$`Pr(>F)`[2] < 0.05) {
      return(mod_select)
    } else {
      return(mod0)
    }
  } else {
    aov1 = anova(mod7, mod_select) # compare to the most complex model
    aov2 = anova(mod_select, mod0) # compare to the least complex model
    print(aov1)
    print(aov2)
    if (aov1$`Pr(>F)`[2] > 0.05 & aov2$`Pr(>F)`[2] < 0.05) {
      return(mod_select)
    } else if (aov1$`Pr(>F)`[2] < 0.05 & aov2$`Pr(>F)`[2] < 0.05) {
      return(mod7)
    } else if (aov1$`Pr(>F)`[2] > 0.05 & aov2$`Pr(>F)`[2] > 0.05) {
      return(mod0)
    } else {
      aov = anova(mod7, mod0)
      if (aov$`Pr(>F)`[2] < 0.05) {
        return(mod7)
      } else {
        return(mod0)
      }
    }
  }
}



# Model selection process---------------------------------------------------------

model_select = function(x, data = df_lite) {
  
  # Fit main effects model
  mod0 = lmer(paste(x, "~ f.age + sex + f.layer + (1 | subject_id) + (1 | subject_id:location)"), data = data, REML = TRUE)
  
  if (isSingular(mod0) == TRUE) {
    mod0 = lmer(paste(x, "~ f.age + sex + f.layer + (1 | subject_id)"), 
                data = data, REML = TRUE)
  }
  
  if (isSingular(mod0) == TRUE) {
    f = as.formula(paste(x, "~ f.age + sex + f.layer"))
    mod0 = eval(bquote(lm(.(f), data = data)))
  }
  
  # Add one-way interaction
  mod1 = update(mod0, . ~ . + f.layer:sex)
  mod2 = update(mod0, . ~ . + f.age:sex)
  mod3 = update(mod0, . ~ . + f.age:f.layer)
  
  # Add two way interaction
  mod4 = update(mod1, . ~ . + f.age:sex)
  mod5 = update(mod1, . ~ . + f.age:f.layer)
  mod6 = update(mod2, . ~ . + f.age:f.layer)
  
  # Add three-way interaction
  mod7 = update(mod4, . ~ . + f.age:f.layer + f.age:sex:f.layer)
  
  mod_list = list(mod0, mod1, mod2, mod3, mod4, mod5, mod6, mod7)
  
  if (class(mod0) == "lmerMod") {
    # anova (for lmer models)
    mod_compare = anova(mod0, mod1, mod2, mod3, mod4, mod5, mod6, mod7)
  } else {
    # AIC (for lm models)
    mod_compare = AIC(mod0, mod1, mod2, mod3, mod4, mod5, mod6, mod7)
  }
  
  # Select model with lowest AIC
  row_nb = which(mod_compare$AIC == min(mod_compare$AIC))
  mod_select = mod_list[row_nb][[1]]
  
  mod_select <<- mod_select
  mod0 <<- mod0
  mod7 <<- mod7
  
  # Model comparisons between selected, simplest, and most complex models
  if (class(mod_select) == "lmerMod") {
    lmer_model_comparison()
  } else {
    lm_model_comparison()
  }
}


# Case bootstrap method ----------------------------------------------------------

case_bootstrap = function(mod, b1 = 625 , b2 = 16) {
  tryCatch(
    {
      if (length(attributes(mod@flist)$names) == 1) {
        case_boot = foreach(B = rep(b1, b2), 
                            .combine = combine_lmeresamp,
                            .packages = c("lmeresampler", "lme4")) %dopar% {
                              bootstrap(mod, .f = fixef, 
                                        resample = c(TRUE, FALSE), 
                                        type = "case", B = B)
                            }
      } else {
        case_boot = foreach(B = rep(b1, b2), 
                            .combine = combine_lmeresamp,
                            .packages = c("lmeresampler", "lme4")) %dopar% {
                              bootstrap(mod, .f = fixef, 
                                        resample = c(TRUE, FALSE, FALSE), 
                                        type = "case", B = B)
                            }
      }
      return(case_boot)
    }, 
    error = function(e) {
      message("An error occurred.")
      print(e)
    }
  )
}


# Simple bootstrap method for lm models --------------------------------------------------------

lm_bootstrap = function(mod, R = 10000) {
  lm_boot = lm.boot(mod, R = R, rows = TRUE)
  return(lm_boot)
}


