## Cleaned up cocirculation import

## Import master spreadsheet
setwd('~/Dropbox/R/ImmuneAgeStructure/')
cocirculation = read.csv('CocirculationData.csv', header = TRUE)
years = as.character(seq(1997, 2015))
source('get.country.data.R')

get.cocirculation.ref = function(Country){
## Goal: matrix with years across columns, H1, H2, H3, G1, G2 down rows for 2015:1918

# 1. Load template from .csv
#      Template includes all years from 1976:1918, which are fixed across countries
template = as.data.frame(read.csv('Cocirculation_template.csv', header = T))
rownames(template) = template[,1]; template = template[,-1]
colnames(template) = 2015:1901

# 2. Fill in 1996:1977 data from Thompson paper
#     
Thompson.data = read.csv('Cocirculation_77-97.csv', header = FALSE, skip = 1, col.names = c('Year', 'Group1', 'Group2', 'Source', 'X'), row.names = as.character(1977:2015))[,1:3]

template['Group1', as.character(1996:1977)] = Thompson.data[as.character(1996:1977), 'Group1']
template['H1', as.character(1996:1977)] = Thompson.data[as.character(1996:1977), 'Group1']
template['H2', as.character(1996:1977)] = 0
template['Group2', as.character(1996:1977)] = Thompson.data[as.character(1996:1977), 'Group2']
template['H3', as.character(1996:1977)] = Thompson.data[as.character(1996:1977), 'Group2']


# 3. Fill in country-specific data from flu net for 2015:1997 where availabl
  #get.country.data is a Function that formats data to be input into template

## Format country -specific data
country.data.1997_2015 = get.country.data(country = Country)
## Quality checks
#    - remove data from years with < 50 specimens, and place NAs throughout the column
if(any(country.data.1997_2015['A', is.na(country.data.1997_2015['A', ]) == FALSE] < 50)) country.data.1997_2015[ , which(country.data.1997_2015['A', ] < 50)] = rep(NA, 4)

## Fill in template
template['H1', as.character(2015:1997)] = country.data.1997_2015['H1', as.character(2015:1997)]
template['H3', as.character(2015:1997)] = country.data.1997_2015['H3', as.character(2015:1997)]
template['H2', as.character(2015:1997)] = 0
template['Group1', as.character(2015:1997)] = colSums(template[c('H1', 'H2'), as.character(2015:1997)])
template['Group2', as.character(2015:1997)] = colSums(template[c('H3'), as.character(2015:1997)])


# 4. Fill in general data where not available
#    Load fallback data
#    fallback.data = Thompson.data[1997:2015]

#    Figure out which years are still NA, and replace
if(any(is.na(template['Group1', ]))){
  proxy.years = as.character(2015:1997)[which(is.na(template['Group1', as.character(2015:1997)]))]
  template['Group1', proxy.years] = Thompson.data[proxy.years, 'Group1']
  template['H1', proxy.years] = Thompson.data[proxy.years, 'Group1'] #H1 was the only G1 virus circulating
  template['H2', proxy.years] = 0                                    #H2 was not circulating
  template['Group2', proxy.years] = Thompson.data[proxy.years, 'Group2']
  template['H3', proxy.years] = Thompson.data[proxy.years, 'Group2'] #H3 was the only G2 virus circulating
}


template
}

