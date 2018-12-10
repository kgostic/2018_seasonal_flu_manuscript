This file contains data and code used to analyze INSIGHT outpatient data, searching for signatures of childhood imprinting protection in seasonal influenza age distributions. All analyses were performed using a multinomial model. A previous version of the analysis using GAMMs can be found in the Archive/ sub directory.


####################
FOLDERS
Archive/ - old files and analysis not used in publication
raw-data/ - raw data files sent by Deborah Wentworth
processed-data/ - processed data files output by R code
-------------------


####################
R Code
00-Import_FLU002-for-multinomial.R - Imports and processes raw data. Also creates multinomial model inputs, including imprintig probabilities, fractions of each age group vaccinated, treated for antivirals, and with underlying conditions.

0func-multinomial_likelihood.R - Defined likelihood function

01-Fit-multinomial-to-INSIGHT002.R - Fits multinomials and performs model comparison.

NOTE - Main text analysis fits a single age curve to all cases in the data. Files wtih ssn-specific flag run an alternate version of the analysis, in which we only analyze data from seasons in which both H1N1 and H3N2 cocirculated, and in which we fit a unique age curve to all cases observed in a given season. This approach has some ability to account for season-specific differences in age distribution, but suffers because we drop large numbers of cases from the data, and also because H1N1 and H3N2 circulation is not well balanced, even in years when both subtypes circulate.
-------------------
