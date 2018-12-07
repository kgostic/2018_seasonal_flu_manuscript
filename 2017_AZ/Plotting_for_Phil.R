## dat is a data frame with the following columns.

#####   INSIGHT data:
# "age" - integer from 18 to 90
# "anyvac" - binary indicating if subject was vaccinated prior to the current season.codes.NH
# "anydx" - binary indicating if subject had underling conditions
# "anyav" - binary indicating if the subject was treated with antivirals
# "season" - {"NH.09.10", "NH.10.11", "SH.11", "SH.10", "NH.11.12", "NH.12.13", "SH.12", "SH.13", "NH.13.14", "NH.14.15", "SH.14"   , "SH.15", "NH.15.16", "SH.16", "NH.16.17", "SH.17"}   
# "country" - {"Denmark", "Spain","Germany","Estonia","USA","Belgium","Portugal","Poland","Austria","UK","Australia", "Thailand",  "Argentina", "Chile", "Greece", "Peru", "Japan"}
# "flutype" - {1 = H1N1, 2 = H3N2, 3 = B, 4 = negative} 




##### AZ data:
# season - {199394, 199495, 200203, 200304, 200405, 200506, 200607, 200708, 200809, 200910, 201011, 201112, 201213, 201314, 201415}
# age - integer from 0 to 97
# subtype - {H1, H3}


## Code below was written for INSIGHT data. 
## I've commented out the country-specific parts since your data is all from the US.

## Plot observed age distributions by country, season and country-season
# ccs = unique(dat$country)
ssns = unique(dat$season)
## Write a function to tabluate ages from 18 to 90 quickly
agerange = 18:90
age.tab = function(agevec){sapply(agerange, FUN = function(aa) sum(agevec == aa))}


#pdf('INSIGHT_agedat_confirmed.pdf', height = 9)
## Overall
valid1 = subset(dat,flutype == "1") # Extract H1N1 cases
valid2 = subset(dat,flutype == "2") # Extract H3N2 cases
# Plot H1N1 cases
plot(agerange, age.tab(valid1$age)/nrow(valid1), xlab = 'age', ylab = 'frequency', xlim = c(18, 90), ylim = c(0, .1), main = paste( 'Overall \nH1 n = ', nrow(valid1), '\nH3 n = ', nrow(valid2)), col = 'dodgerblue') 
# Plot H3N2 cases in same panel
points(agerange, age.tab(valid2$age)/nrow(valid2), col = 'firebrick')

# Repeat subset and plot procedure for each country and each season
par(mfrow = c(4,3))
# for(cc in ccs){
#   valid1 = subset(dat, country == cc & flutype == "1")
#   valid2 = subset(dat, country == cc & flutype == "2")
#   plot(agerange, age.tab(valid1$age)/nrow(valid1), xlab = 'age', ylab = 'frequency', xlim = c(18, 90), ylim = c(0, .1), main = paste(cc, '\nH1 n = ', nrow(valid1), '\nH3 n = ', nrow(valid2)), col = 'dodgerblue') 
#   points(agerange, age.tab(valid2$age)/nrow(valid2), col = 'firebrick')
# }

par(mfrow = c(4,3))
for(ss in ssns){
  valid1 = subset(dat, season == ss & flutype == "1")
  valid2 = subset(dat, season == ss & flutype == "2")
  plot(agerange, age.tab(valid1$age)/nrow(valid1), xlab = 'age', ylab = 'frequency', xlim = c(18, 90), ylim = c(0, .1), main = paste(ss, '\nH1 n = ', nrow(valid1), '\nH3 n = ', nrow(valid2)), col = 'dodgerblue') 
  points(agerange, age.tab(valid2$age)/nrow(valid2), col = 'firebrick')
}
#dev.off()





## For the INSIGHT data, plot all tested cases. Denominatior data is not available for AZ data.
## Overall
# Plot H1N1 cases
plot(agerange, age.tab(dat.002)/nrow(dat.002), xlab = 'age', ylab = 'frequency', xlim = c(18, 90), ylim = c(0, .1), main = paste( 'Overall \nn = ', nrow(dat.002))) 
# Plot H3N2 cases in same panel


# Repeat subset and plot procedure for each country and each season
# par(mfrow = c(4,3))
# for(cc in ccs){
#   valid1 = subset(dat.002, country == cc)
#   plot(agerange, age.tab(valid1$age)/nrow(valid1), xlab = 'age', ylab = 'frequency', xlim = c(18, 90), ylim = c(0, .1), main = paste(cc, '\nn = ', nrow(valid1)))
# }

par(mfrow = c(4,3))
for(ss in ssns){
  valid1 = subset(dat.002, season == ss)
  plot(agerange, age.tab(valid1$age)/nrow(valid1), xlab = 'age', ylab = 'frequency', xlim = c(18, 90), ylim = c(0, .1), main = paste(ss, '\nn = ', nrow(valid1)))
}