## Get US weights
source('~/Dropbox/R/ImmuneAgeStructure/Infection.age.structure.AB_vax.R')
coverage = as.numeric(read.csv('~/Dropbox/R/ImmuneAgeStructure/USA_vax_assumpitons.csv', skip = 1, header = FALSE, stringsAsFactors = FALSE))

vax.assumptions = coverage; names(vax.assumptions) = 1918:2022
efficacy = 0.7

           


weights = get.type.weights.AB.vax(years.out = 1993:2015, Countries.out = 'USA', vax.vector = vax.assumptions*efficacy, type = 5)

weights.master.1 = weights[[1]]
weights.master.2 = weights[[2]]
weights.master.3 = weights[[3]]
weights.master.naiive = weights[[4]]

setwd('~/Dropbox/R/2017_seasonal_flu/')


