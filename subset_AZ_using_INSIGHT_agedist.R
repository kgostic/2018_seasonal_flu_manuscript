## Compare age dists tested
## Create supplementary figure with distributions of all tested cases from INSIGHT

## Start in main seasonal flu folder
setwd('~/Dropbox/R/2018_seasonal_flu/')
## Load processed INSIGHT data
dat.002 = read.csv(file = '2017_INSIGHT/processed-data/INSIGHT002_processed.csv')

## OUTPUTS
outfile1 = 'figures/all-tested-cases.pdf'

## Extract age distribution of all tested cases
tested = hist(x = dat.002$age, breaks = 18:91)$counts
confirmed = hist(x = subset(dat.002, flutype %in% c(1,2))$age, breaks = 18:91)$counts


## Sample 2000 H1N1 caess and 2000 H3N2 cases according to the weights specified in the INSIGHT age dist
setwd('2017_AZ/')
source('00-Inputs_multinomial.R')
setwd('../')
raw.dat$yr = as.numeric( sub(pattern = "(\\d{4})\\d{2}", replacement="\\1", x = raw.dat$season) )+1
raw.dat$age = raw.dat$yr - raw.dat$birthyear


## Add sampling weights defined by INSIGHT age distribution to raw AZ data
raw.dat$weights = 0 # Initialize. No probability of sampling ages other than 18:90
wts = tested/sum(tested)
for(aa in 18:90){ # For possible age groups, fill in probabilities given by tested
  raw.dat$weights[which(raw.dat$age == aa)] = wts[aa-17]
}

# CHeck
plot(raw.dat$age, raw.dat$weights)  




### Sample AZ data using age distribution observed in INSIGHT
valid.ages = subset(raw.dat, age %in% 18:90)
samprows = sample(1:nrow(valid.ages), prob = valid.ages$weights, size = 4000, replace = FALSE)
AZsubsample = valid.ages[samprows, ]


### Plot observed distributions of H1N1 and H3N2 cases from the subsample
fullH1dist =  hist(x = subset(raw.dat, subtype == 'H1')$age, breaks = 0:100)$counts
fullH3dist =  hist(x = subset(raw.dat, subtype == 'H3')$age, breaks = 0:100)$counts

H1dist = hist(x = subset(AZsubsample, subtype == 'H1')$age, breaks = 18:91)$counts
H3dist = hist(x = subset(AZsubsample, subtype == 'H3')$age, breaks = 18:91)$counts



par(mfrow = c(2,1))
plot(0:99, fullH1dist, col = 'dodgerblue', xlab = 'age', ylab = 'count', main = 'full')
points(0:99, fullH3dist, col = 'firebrick1')

plot(18:90, H1dist, col = 'dodgerblue', xlab = 'age', ylab = 'count', main = 'subset')
points(18:90, H3dist, col = 'firebrick1')


AZ_all_confirmed = hist(raw.dat$age, breaks = 0:100)$counts

pdf(outfile1)
par(mfrow = c(2,1), mar = c(4,4,2,1))
# Counts
plot(0:99, c(rep(0, 18), tested, rep(0, 9)), xlab = 'age', ylab = 'count', pch = 16, ylim = c(0, 700))
points(0:99, c(rep(0, 18), confirmed, rep(0,9)), pch = 4, col = 'blue')
points(0:99, AZ_all_confirmed, col = 'green3')
legend('topright', legend = c('INSIGHT - all tested cases', 'INSIGHT - all confirmed influenza A', 'Arizona - all confirmed influenza A'), pch = c(16, 4, 1), col = c('black', 'blue', 'green3'))
mtext(text = 'A', side = 3, line = .5, at = -5, font = 2)

plot(0:99, c(rep(0, 18), tested/sum(tested), rep(0, 9)), xlab = 'age', ylab = 'frequency', ylim = c(0, .04), pch = 16)
points(0:99, c(rep(0, 18), confirmed/sum(confirmed), rep(0,9)), pch = 4, col = 'blue')
points(0:99, AZ_all_confirmed[1:100]/sum(AZ_all_confirmed[1:100]), col = 'green3')
mtext(text = 'B', side = 3, line = .5, at = -5, font = 2)
dev.off()
