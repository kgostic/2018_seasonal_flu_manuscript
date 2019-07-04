This file contains data and code used to analyze AZ influenza surveillance data, searching for signatures of childhood imprinting protection in seasonal influenza age distributions. All analyses were performed using a multinomial model.

####################
FOLDERS
raw-data/ - raw data files sent by Shane Brady
processed-data/ - processed data files output by R code
-------------------


####################
R Code
00-Inputs_multinomial.R - Imports and processes raw data. Also creates multinomial model inputs, including imprintig probabilities, and age group-specific indicators.

0func-likelihood.R - Defines the likelihood.


01-multinomial_model_comparison.R - Fits multinomials and performs model comparison.
Plots outputs - Individual model fits and parameter estimates, overall fits to data, and AIC score comparison

02-calculate-Likelihood-profiles.R - Calculates likelihood profile confidence intervals.

03-format-model-comp-into-latex-table.R - Formats all outputs into tables for easy reporting. (Manuscript Table 2)

04-data-summary.R - Summarizes case counts by year. (Manuscript Table 1)


-------------------
