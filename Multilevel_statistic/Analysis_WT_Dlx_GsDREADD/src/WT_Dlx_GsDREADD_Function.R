# Load necessary libraries
library(lme4)         # For fitting linear mixed-effects models
library(pbkrtest)     # For Kenward-Roger approximation for mixed models
library(simpleboot)   # For simple bootstrap methods
library(lmeresampler) # For resampling methods for mixed models
library(foreach)      # For looping constructs for parallel execution
library(doParallel)  # For parallel backend to the foreach 

#define MLM w/ three-way interaction --------------------------------------------

lmer_three_way = function(x, data = df_lite) {
  mod = lmer(paste(x, "~ p15_cno * sex * f.layer + (1 | subject_id) + (1 | subject_id:location)"), data = data, REML = TRUE)
  return(mod)
}


# Compare lmer models -------------------------------------------------------

lmer_model_comparison = function() {
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


# Compare lm models-------------------------------------------------------------------

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


# Model selection process--------------------------------------------------------------------------

model_select = function(x, data = df_lite) {
  
  # main effects only
  mod0 = lmer(paste(x, "~ p15_cno + sex + f.layer + (1 | subject_id) + (1 | subject_id:location)"), data = data, REML = TRUE)
  
  if (isSingular(mod0) == TRUE) {
    mod0 = lmer(paste(x, "~ p15_cno + sex + f.layer + (1 | subject_id)"), 
                data = data, REML = TRUE)
  }
  
  if (isSingular(mod0) == TRUE) {
    f = as.formula(paste(x, "~ p15_cno + sex + f.layer"))
    # mod0 = lm(paste(x, "~ p15_cno + sex + f.layer"), data = data)
    mod0 = eval(bquote(lm(.(f), data = data)))
  }
  
  # add one-way interaction
  mod1 = update(mod0, . ~ . + f.layer:sex)
  mod2 = update(mod0, . ~ . + p15_cno:sex)
  mod3 = update(mod0, . ~ . + p15_cno:f.layer)
  
  # add two way interaction
  mod4 = update(mod1, . ~ . + p15_cno:sex)
  mod5 = update(mod1, . ~ . + p15_cno:f.layer)
  mod6 = update(mod2, . ~ . + p15_cno:f.layer)
  
  # add three-way interaction
  mod7 = update(mod4, . ~ . + p15_cno:f.layer + p15_cno:sex:f.layer)
  
  mod_list = list(mod0, mod1, mod2, mod3, mod4, mod5, mod6, mod7)
  
  if (class(mod0) == "lmerMod") {
    # anova (for lmer models)
    mod_compare = anova(mod0, mod1, mod2, mod3, mod4, mod5, mod6, mod7)
  } else {
    # AIC (for lm models)
    mod_compare = AIC(mod0, mod1, mod2, mod3, mod4, mod5, mod6, mod7)
  }
  
  # find the model with lowest AIC
  row_nb = which(mod_compare$AIC == min(mod_compare$AIC))
  mod_select = mod_list[row_nb][[1]]
  
  mod_select <<- mod_select
  mod0 <<- mod0
  mod7 <<- mod7
  
  # model comparisons between selected, simplest, and most complex models
  if (class(mod_select) == "lmerMod") {
    lmer_model_comparison()
  } else {
    lm_model_comparison()
  }
}


# Case bootstrap ----------------------------------------------------------

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


# Simple bootstrap --------------------------------------------------------

lm_bootstrap = function(mod, R = 10000) {
  lm_boot = lm.boot(mod, R = R, rows = TRUE)
  return(lm_boot)
}

