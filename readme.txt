This repository contains the datasets and R scripts required to reproduce the results presented in the paper. 
The analysis is divided into two main physiological modalities: EDA and HRV.

Directory Structure
The repository is organized into two main folders: /EDA and /HRV. Each folder contains:

1. Data Tables
DATA_EDA / DATA_HRV: CSV files containing the calculated IPCC and its predictors.

2. Statistical Analysis (Scripts)
LMER_EDA / LMER_HRV: R scripts for the Linear Mixed-Effects Models (Model 1 and Model 2).

Note on Performance: These scripts are optimized for fast compilation. They utilize pre-calculated pairwise comparison results to avoid lengthy computation times.

EDA_PAIRWISE_COMPARISONS / HRV_PAIRWISE_COMPARISONS: Original scripts used to generate the post-hoc pairwise comparisons for Experimental Conditions and Genders. Running these scripts from scratch may require significant computational time.

3. Pre-saved Results
For faster reproduction of the LMER models, we provide the following pre-saved outputs:

results_EDA: Post-hoc results for EDA Condition effects.

results_HRV_condition: Post-hoc results for HRV Condition effects.

results_HRV_genders: Post-hoc results for EDA Gender effects.

Instructions for Use
Fast Replication: To verify the main LMER results and generate the summary tables, run LMER_EDA.R or LMER_HRV.R. These will automatically load the pre-saved results.

Full Computation: To re-calculate the bootstrap-based pairwise comparisons, refer to the _PAIRWISE_COMPARISONS scripts. Note that this process is computationally intensive.
