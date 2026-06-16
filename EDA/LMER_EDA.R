# =====================================================================
# LINEAR MIXED-EFFECTS MODELS - EDA
# =====================================================================

# Setup
cat("\014")
rm(list = ls())

# Load required packages
library(lmerTest)
library(lme4)
library(pbkrtest)
library(interactions)
library(ggplot2)
library(emmeans)

set.seed(42)

# =====================================================================
# DATA IMPORT AND PREPROCESSING
# =====================================================================

# Get file path relative to script location
script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
file_path <- file.path(script_dir, "DATA_EDA.csv")

# Import and filter data (EDAsymp band: 0.05-0.25 Hz)
my_data <- read.csv(file_path, header = TRUE)
my_data <- subset(my_data, FREQUENCY >= 0.05 & FREQUENCY <=0.25)

# Prepare variables
my_data$ID <- as.factor(my_data$ID)
my_data$IPCC <- log(as.numeric(my_data$DC) + 1)
my_data$FREQUENCY <- as.numeric(my_data$FREQUENCY) - mean(my_data$FREQUENCY)
my_data$CONDITION <- factor(my_data$CONDITION, 
                            levels = c("4", "5", "3", "1", "2"),
                            labels = c("Single\nResting", "Paired\nResting", 
                                       "Paired\nWatching", "Following", "Leading"))
my_data$GENDERS <- factor(my_data$GENDERS, 
                         levels = c("0", "1", "2", "3"),
                         labels = c("f->m", "m->f", "m->m", "f->f"))

# Set contrasts
contrasts(my_data$CONDITION) <- contr.sum
contrasts(my_data$GENDERS) <- contr.sum
contrasts(my_data$ID) <- contr.sum

# Setting controls for LMER model
ctrl <- lmerControl(
  optimizer = "bobyqa",
  optCtrl = list(maxfun = 2e5)
)

N<-1000

# =====================================================================
# MODEL 1
# =====================================================================

message("=== MODEL 1 ===\n")

mod_1 <- lmer(IPCC ~ CONDITION + GENDERS + FREQUENCY + (1 | ID), 
              data = my_data, control = ctrl)

# Function for parametric bootstrap model comparison
test_effect <- function(full_model, term, nsim = N) {
  null_model <- update(full_model, as.formula(paste(". ~ . -", term)), control = ctrl)
  pb_results <- PBmodcomp(full_model, null_model, nsim = nsim)
  message(sprintf("\n--- %s ---", term))
  print(summary(pb_results), digits = 4)
  return(pb_results)
}

# Test fixed effects
test_effect(mod_1, "CONDITION")
test_effect(mod_1, "GENDERS")
test_effect(mod_1, "FREQUENCY")



# =====================================================================
# GRAPH - CONDITION PAIRWISE COMPARISONS
# =====================================================================

message("=== CONDITION PAIRWISE COMPARISONS ===\n")

# Uploading CONDITION pairwise comparisons results
file_path <- file.path(script_dir, "results_EDA.csv")
results_CONDITION_PAIRWISE <- read.csv(file_path, header = TRUE)
print(results_CONDITION_PAIRWISE)


ggplot(my_data, aes(x = CONDITION, y = IPCC)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0, size = 1.2) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  labs(
    title = "IPCC in EDA across Conditions",
    x = "CONDITION",
    y = "IPCC"
  ) +
  scale_x_discrete(labels = function(x) gsub(" ", "\n", x)) +
  theme_minimal() +
  theme(
    panel.border = element_rect(colour = "grey60", fill = NA, size = 1),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b = 10)),
    axis.title = element_text(face = "bold", size = 14),
    axis.title.x = element_text(vjust = -0.5, margin = margin(t = 0)),
    axis.title.y = element_text(angle = 90, vjust = 1, hjust = 0, margin = margin(r = 20)),
    axis.text = element_text(face = "bold", size = 12),
    legend.position = "none"
  )

# =====================================================================
# MODEL 2
# =====================================================================

message("\n=== MODEL 2 ===\n")

# Filter and prepare data for Model 2
my_data_m2 <- subset(my_data, CONDITION %in% c("Paired\nWatching", "Following", "Leading"))
my_data_m2$CONDITION <- droplevels(my_data_m2$CONDITION)

# Center arousal and valence
my_data_m2$AROUSAL <- as.numeric(my_data_m2$AROUSAL) - mean(my_data_m2$AROUSAL)
my_data_m2$VALENCE <- as.numeric(my_data_m2$VALENCE) - mean(my_data_m2$VALENCE)

# Reset contrasts for filtered data
contrasts(my_data_m2$CONDITION) <- contr.sum
contrasts(my_data_m2$GENDERS) <- contr.sum
contrasts(my_data_m2$ID) <- contr.sum

# Fit model with interaction terms
mod_2 <- lmer(IPCC ~ AROUSAL*CONDITION + VALENCE*CONDITION + GENDERS + FREQUENCY + (1 | ID), 
              data = my_data_m2, control = ctrl)

# Test interactions terms
test_effect(mod_2, "AROUSAL:CONDITION")
test_effect(mod_2, "VALENCE:CONDITION")

# Test single predictor terms
test_effect(mod_2, "AROUSAL*CONDITION")
test_effect(mod_2, "VALENCE*CONDITION")
test_effect(mod_2, "GENDERS")
test_effect(mod_2, "FREQUENCY")

# =====================================================================
# SLOPES ESTIMATION
# =====================================================================
#Estimating the effect of Arousal and Valence on IPCC in three experimental 
#conditions involving emotional stimuLi exposure: Paired Watching, Following, Leading.

# We obtained the slopes using emtrends, while the significance was tested 
#with PBtest, coehrently with the rest of the analysis

emm_options(pbkrtest.limit = 4000, lmerTest.limit = 4000)

#AROUSAL
message("\n Arousal*Condition Slopes \n")
arousal_slopes <- summary(emtrends(mod_2, ~ CONDITION, var = "AROUSAL"))[, c("CONDITION", "AROUSAL.trend")]
print(arousal_slopes)

# =====================================================================
# TESTING SLOPES
# =====================================================================

# Creating conditional interaction variables to isolate the effects of 
# Arousal and Valence within each specific level of the CONDITION factor. 
# This dummy-coding approach allows us to test the significance of the 
# slopes for each condition independently by comparing full and reduced 
# models using the PBmodcomp function.

# Interaction Effects Isolation
my_data_m2$AROUSAL_PW <- my_data_m2$AROUSAL*(my_data_m2$CONDITION == "Paired\nWatching")
my_data_m2$VALENCE_PW<-my_data_m2$VALENCE*(my_data_m2$CONDITION=="Paired\nWatching")
my_data_m2$AROUSAL_F<-my_data_m2$AROUSAL*(my_data_m2$CONDITION=="Following")
my_data_m2$VALENCE_F<-my_data_m2$VALENCE*(my_data_m2$CONDITION=="Following")
my_data_m2$AROUSAL_L<-my_data_m2$AROUSAL*(my_data_m2$CONDITION=="Leading")
my_data_m2$VALENCE_L<-my_data_m2$VALENCE*(my_data_m2$CONDITION=="Leading")

#AROUSAL*PAIRED_WATCHING
message("\n Testing Slope : AROUSAL*PAIRED WATCHING \n")
m_f<-lmer(IPCC ~ CONDITION + AROUSAL + VALENCE + AROUSAL_F + AROUSAL_L + VALENCE_F + VALENCE_L + GENDERS + FREQUENCY + (1|ID), data=my_data_m2, control = ctrl)
m_r<-lmer(IPCC ~ CONDITION + VALENCE + AROUSAL_F + AROUSAL_L + VALENCE_F + VALENCE_L + GENDERS + FREQUENCY + (1|ID), data=my_data_m2, control = ctrl)
results_A_PW<-PBmodcomp(m_f, m_r, nsim = N)
print(results_A_PW)

#AROUSAL*FOLLOWING
message("\n Testing Slope : AROUSAL*FOLLOWING \n")
m_f<-lmer(IPCC ~ CONDITION + AROUSAL + VALENCE + AROUSAL_PW + AROUSAL_L + VALENCE_PW + VALENCE_L + GENDERS + FREQUENCY + (1|ID), data=my_data_m2, control = ctrl)
m_r<-lmer(IPCC ~ CONDITION + VALENCE + AROUSAL_PW + AROUSAL_L + VALENCE_PW + VALENCE_L + GENDERS + FREQUENCY + (1|ID), data=my_data_m2, control = ctrl)
results_A_F<-PBmodcomp(m_f, m_r, nsim = N)
print(results_A_F)

#AROUSAL*LEADING
message("\n Testing Slope : AROUSAL*LEADING \n")
m_f<-lmer(IPCC ~ CONDITION + AROUSAL + VALENCE + AROUSAL_F + AROUSAL_PW + VALENCE_F + VALENCE_PW + GENDERS + FREQUENCY + (1|ID), data=my_data_m2, control = ctrl)
m_r<-lmer(IPCC ~ CONDITION + VALENCE + AROUSAL_F + AROUSAL_PW + VALENCE_F + VALENCE_PW + GENDERS + FREQUENCY + (1|ID), data=my_data_m2, control = ctrl)
results_A_L<-PBmodcomp(m_f, m_r, nsim = N)
print(results_A_L)

# =====================================================================
# GRAPHS - AROUSAL
# =====================================================================
interact_plot(
  mod_2,
  pred =AROUSAL,
  modx = CONDITION,
  interval = TRUE,
  plot.points = FALSE,
  y.label = "IPCC",
  x.label = "AROUSAL",
  legend.main = "CONDITION",
  modx.labels = c("Paired\nWatching", "Following", "Leading"),
  facet.modx = TRUE
) +
  labs(title = "Variation of IPCC in EDA with AROUSAL") +
  scale_color_manual(values = c("#2A9D8F", "#264653", "#E9C46A")) +
  scale_fill_manual(values = c("#2A9D8F", "#264653", "#E9C46A")) +
  theme_minimal() +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 9),
    plot.title = element_text(face = "bold", size = 10, hjust = 0.5),
    panel.border = element_rect(color = "grey60", fill = NA, linewidth = 0.1),
    panel.spacing = unit(0.5, "lines"),
    axis.title = element_text(face = "bold", size = 8),
    axis.text = element_text(face = "bold", size = 8)
  )



