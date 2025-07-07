# ---- Load and install required libraries ----

# List of required packages
required_packages <- c(
  "lme4",          # For fitting linear mixed-effects models
  "pbkrtest",      # For Kenward-Roger approximation for mixed models
  "simpleboot",    # For simple bootstrap methods
  "lmeresampler",  # For resampling methods for mixed models
  "foreach",       # For looping constructs for parallel execution
  "doParallel",    # For parallel backend to the foreach
  "doRNG"          # For parallel Reproducibility
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
    # anova (for lmer models): comparees models using likelihood ratio tests
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
#' This function performs a case bootstrap for linear mixed-effects models
#' using parallel processing via `foreach` and `lmeresampler`.
#'
#' The bootstrapping strategy is "case bootstrap", meaning it resamples entire
#' cases (rows of data), which implicitly accounts for the hierarchical structure.
#' It automatically adjusts the `resample` argument based on the number of
#' grouping factors in the mixed-effects model.
#'
#' **Important Note for Parallel Processing:**
#' This function relies on `foreach` and `doParallel` for parallel execution.
#' To ensure the function runs in parallel and produces reproducible results,
#' you *must* set up and register a parallel backend (e.g., using `makeCluster()`
#' and `registerDoParallel()`) *before* calling this function.
#' After all parallel computations are complete, remember to `stopCluster()`.
#'
#' @param mod An `lmerMod` object, typically created by `lme4::lmer()`.
#' @param b1 An integer specifying the number of bootstrap replicates to be performed
#'           by each parallel worker. The total number of replicates will be `b1 * b2`.
#'           Defaults to 625.
#' @param b2 An integer specifying the number of parallel workers (or batches) over
#'           which the total bootstrap replicates will be distributed.
#'           Defaults to 16.
#'
#' @return A `lmeresamp.bootstrap` object containing the results of the case bootstrap.
#'         This object can be further analyzed using methods provided by `lmeresampler`
#'         (e.g., `summary()`, `confint()`).
#'
#' @details
#' The function uses `lmeresampler::bootstrap()` with `type = "case"`.
#'
#' - If the model has only one grouping factor (e.g., `(1|Subject)`), `resample = c(TRUE, FALSE)`
#'   is used, indicating resampling of subjects but not residuals.
#' - If the model has multiple grouping factors (e.g., `(1|Subject) + (1|Item)`),
#'   `resample = c(TRUE, FALSE, FALSE)` is used, indicating resampling of subjects
#'   but not residuals or lower-level effects.
#'
#' **Reproducibility:** `registerDoRNG(123)` is called inside the function to ensure
#' that the parallel random number streams are reproducible. This means that
#' calling the function multiple times with the same inputs will yield identical
#' bootstrap results, provided the parallel setup is consistent. The choice of
#' `123` as the seed is arbitrary; any integer would work.
#'
#'
case_bootstrap = function(mod, b1 = 625 , b2 = 16) {
 
   registerDoRNG(123) # Choose any integer seed you like
  
   tryCatch(
    {
      if (length(attributes(mod@flist)$names) == 1) {
        case_boot = foreach(B_iter = rep(b1, b2), 
                            .combine = combine_lmeresamp,
                            .packages = c("lmeresampler", "lme4")) %dopar% {
                              bootstrap(mod, .f = fixef, 
                                        resample = c(TRUE, FALSE), 
                                        type = "case", B = B_iter)
                            }
      } else {
        case_boot = foreach(B_iter = rep(b1, b2), 
                            .combine = combine_lmeresamp,
                            .packages = c("lmeresampler", "lme4")) %dopar% {
                              bootstrap(mod, .f = fixef, 
                                        resample = c(TRUE, FALSE, FALSE), 
                                        type = "case", B = B_iter)
                            }
      }
      # Return the combined bootstrap results
      return(case_boot)
    },
    
    # Error handling block: catches any errors during execution
    error = function(e) {
      message("An error occurred.")
      print(e)
    }
  )
}


# Simple bootstrap method for lm models --------------------------------------------------------
#'
#' @description
#' This function performs a non-parametric bootstrap for a linear model (LM)
#' using the `simpleboot::lm.boot` function. It resamples rows of the original
#' dataset with replacement to generate bootstrap samples and then refits the
#' linear model for each sample.
#'
#' The function sets a fixed random seed internally to ensure that the
#' bootstrap results are fully reproducible.
#'
#' @param mod An `lm` object, typically created by `stats::lm()`.
#' @param R An integer specifying the number of bootstrap replicates to perform.
#'          Defaults to 10000.
#'
#' @return An object of class `lm.boot`, which is a list containing the
#'         bootstrap results (e.g., coefficients, residuals, fitted values)
#'         and can be further analyzed using `summary()` or `plot()` methods
#'         from the `simpleboot` package.
#'
#' @details
#' The bootstrap method used here is a "case bootstrap" as it resamples entire
#' rows (cases) from the original data. This is suitable for standard linear
#' models where independence of observations is assumed.
#'
#' **Reproducibility:** A fixed seed (`set.seed(10)`) is set at the beginning
#' of the function. This guarantees that calling `lm_bootstrap` multiple times
#' with the same input model and `R` value will produce identical bootstrap
#' results. The choice of `10` as the seed is arbitrary; any integer would work
#' for reproducibility.
#' 
#' 
lm_bootstrap = function(mod, R = 10000) {
  # Set the seed for reproducible random number generation within this function.
  # This ensures that calling this function multiple times with the same inputs
  # will always yield the exact same bootstrap results.
  set.seed(10)
  # Perform the linear model bootstrap using simpleboot::lm.boot.
  # 'rows = TRUE' indicates that entire rows (cases) of the data are resampled.
  lm_boot = lm.boot(mod, R = R, rows = TRUE)
  # Return the results of the bootstrap.
  return(lm_boot)
}


