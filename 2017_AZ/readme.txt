This file contains data and code used to analyze AZ influenza surveillance data, searching for signatures of childhood imprinting protection in seasonal influenza age distributions. All analyses were performed using a multinomial model. A previous version of the analysis using GAMMs can be found in the Archive/ sub directory.


####################
FOLDERS
Archive/ - old files and analysis not used in publication
raw-data/ - raw data files sent by Deborah Wentworth
processed-data/ - processed data files output by R code
-------------------


####################
R Code
00-Inputs_multinomial.R - Imports and processes raw data. Also creates multinomial model inputs, including imprintig probabilities, and age group-specific indicators.


01-multinomial_model_comparison.R - Fits multinomials and performs model comparison.
Plots outputs - Individual model fits and parameter estimates, overall fits to data, and AIC score comparison
-------------------
