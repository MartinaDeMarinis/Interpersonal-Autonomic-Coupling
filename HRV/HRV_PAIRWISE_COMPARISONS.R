# =====================================================================
# LIBRARIES AND SETUP
# =====================================================================
cat("\014")
rm(list = ls())

library(lmerTest)
library(lme4)
library(pbkrtest)
library(dplyr)
library(purrr)
library(ggplot2)

set.seed(42)

# =====================================================================
# DATA IMPORT AND PREPROCESSING
# =====================================================================
script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(script_dir)
file_path <- file.path(script_dir, "DATA_HRV.csv")

my_data <- read.csv(file_path, header = TRUE) %>%
  filter(FREQUENCY >= 0.15 & FREQUENCY <= 0.40) %>% #selecting HF band [0.15-0.40] Hz
  mutate(
    ID = as.factor(ID),
    logDC = log(as.numeric(DC) + 1),
    # Centering Frequency
    FREQUENCY = as.numeric(FREQUENCY) - mean(as.numeric(FREQUENCY)),
    CONDITION = factor(CONDITION, 
                       levels = c("4", "5", "3", "1", "2"),
                       labels = c("Single Resting", "Paired Resting", "Paired Watching", "Following", "Leading")),
    GENDERS = factor(GENDERS, 
                    levels = c("0", "1", "2", "3"),
                    labels = c("f->m", "m->f", "m->m", "f->f"))
  )

# Sum contrasts for ANOVA-like interpretation
contrasts(my_data$CONDITION) <- contr.sum
contrasts(my_data$GENDERS) <- contr.sum
contrasts(my_data$ID) <- contr.sum 

# Define Global Model
ctrl <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
mod_full <- lmer(logDC ~ CONDITION + GENDERS + FREQUENCY + (1 | ID), data = my_data, control = ctrl)

# =====================================================================
# FUNCTION FOR PAIRWISE COMPARISONS (Parametric Bootstrap)
# =====================================================================
# This function automates the creation of reduced models for PBtest
run_pairwise_pb <- function(target_var, data_orig, full_model, nsim = 1000) {
  
  message(paste("\n>>> STARTING PAIRWISE COMPARISONS FOR:", target_var))
  
  levels_list <- levels(data_orig[[target_var]])
  pairs <- combn(levels_list, 2, simplify = FALSE)
  
  test_single_pair <- function(pair) {
    lv1 <- pair[1]; lv2 <- pair[2]
    message("Testing: ", lv1, " vs ", lv2)
    
    # Merge the two levels being compared into one category
    dat <- data_orig
    dat$merged_var <- as.character(dat[[target_var]])
    dat$merged_var[dat$merged_var %in% c(lv1, lv2)] <- "merged_level"
    dat$merged_var <- factor(dat$merged_var)
    
    # Update formula dynamically by replacing target_var with merged_var
    base_formula_str <- as.character(formula(full_model))[3]
    reduced_formula_str <- gsub(target_var, "merged_var", base_formula_str)
    reduced_formula <- as.formula(paste("logDC ~", reduced_formula_str))
    
    tryCatch({
      mod_reduced <- lmer(reduced_formula, data = dat, control = lmerControl(optimizer = "bobyqa"))
      pb <- PBmodcomp(full_model, mod_reduced, nsim = nsim)
      
      tibble(
        Comparison = paste(lv1, "vs", lv2),
        Pval_raw = pb$test["PBtest", "p.value"],
        Stat = pb$test["PBtest", "stat"],
        Df = pb$LRTstat["df"]
      )
    }, error = function(e) {
      warning("Error in comparison ", lv1, " vs ", lv2, ": ", e$message)
      return(tibble(Comparison = paste(lv1, "vs", lv2), Pval_raw = NA, Stat = NA, Df = NA))
    })
  }
  
  # Run all pairs and apply False Discovery Rate (FDR) correction
  results <- map_df(pairs, test_single_pair) %>%
    mutate(Pval_FDR = p.adjust(Pval_raw, method = "fdr"))
  
  return(results)
}

# =====================================================================
# EXECUTION: STATISTICAL ANALYSIS
# =====================================================================

# 1. Condition Comparisons
results_cond <- run_pairwise_pb("CONDITION", my_data, mod_full)
print(results_cond)

# 2. Genders Comparisons
results_genders <- run_pairwise_pb("GENDERS", my_data, mod_full)
print(results_genders)

write.csv(results_cond, file.path(script_dir, "results_HRV_condition.csv"), row.names = FALSE)
write.csv(results_cond, file.path(script_dir, "results_HRV_genders.csv"), row.names = FALSE)
