setwd('~/Dropbox/R/2018_seasonal_flu/')

## OUPTUT A PDF
setwd('~/Dropbox/R/2018_seasonal_flu/')
plot1 = 'figures/Akaike_wts.pdf'

## Plot Akaike weight barplots

## Load fits to AZ and INSIGHT data
load('2017_AZ/processed-data/AZ_model_fits.RData')
load('2017_INSIGHT/processed-data/INSIGHT_fitted_models.RData')

## Make sure both weights vectors are in the same order
AZ.imp.weights
INSIGHT.wts = imp.type.weights; INSIGHT.wts
## Looks good, N, S, G, A
## Setup color vector
cols = rev(c('purple', 'dodgerblue', 'limegreen', 'goldenrod'))

## Format into a matrix
pdf(plot1, width = 2, height = 3)
par(mar = c(6,3.5,1,.5), mgp = c(2, 1, 0))
pmat = matrix(c(rev(AZ.imp.weights), rev(INSIGHT.wts)), ncol = 2, dimnames = list(c('A', 'G', 'S', 'N'),  c('AZ', 'INSIGHT')))

## Barplot
xx=barplot(pmat, col = cols, border = cols, beside = FALSE, ylab = 'Akaike weight', cex.names = .9)
legend(.3, -.3, legend = c('NA subtype', 'HA subtype', 'HA group', 'no imprinting'), fill = rev(cols), bty = 'n', xpd = NA, cex = .7, border = NA)
dev.off()
