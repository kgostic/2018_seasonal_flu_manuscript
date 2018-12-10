library(tidyverse)
setwd('~/Dropbox/R/2017_seasonal_flu/')

## By year (season)
calc_expected_dem_AG <- function(dem_data,
                                 age_groups,
                                 seasons){
    #######################
    #
    # Given demography data and age groups, 
    # this function calculates the proportion
    # of the population in each age group
    #
    #######################
    
    filtered_dem_data = dem_data %>% filter(Year %in% seasons)
    all_fracs = data.frame(Year=seasons)
    
    for (i in 1:nrow(age_groups)) {
        lb = age_groups[i, 'LB']
        ub = age_groups[i, 'UB']
        colname = paste(lb, ub, sep='-')
        filtered_dem_data %>% 
            group_by(Year) %>% 
            filter(Age <= ub & Age >=lb) %>% 
            summarize(temp=sum(pop_frac)) -> fracs
        colnames(fracs) = c('Year', colname)
        all_fracs = merge(all_fracs, fracs, by='Year')
    }

    rownames(all_fracs) = all_fracs$Year
    all_fracs = all_fracs %>% select(-Year) 
    return(all_fracs / rowSums(all_fracs))
}


calc_inc_AG <- function(inc_data,
                        age_groups,
                        seasons){
    #######################
    #
    # Given incidence data and age groups, 
    # this function calculates the number of
    # cases in each age group for H1 and H3
    #
    #######################
    
    filtered_inc_data = inc_data %>% filter(Season %in% seasons)
    inc_h1 = inc_h3 = data.frame(Year=seasons)
    
    for (i in 1:nrow(age_groups)){
        lb = age_groups[i, 'LB']
        ub = age_groups[i, 'UB']
        colname = paste(lb, ub, sep='-')
        
        filtered_inc_data %>% 
            group_by(Season) %>% 
            filter(Age <= ub & Age >=lb) %>% 
            summarize(total=sum(subtype == 'H1')) -> totals_h1
        
        filtered_inc_data %>% 
            group_by(Season) %>% 
            filter(Age <= ub & Age >=lb) %>% 
            summarize(total=sum(subtype == 'H3')) -> totals_h3
        
        colnames(totals_h1) = colnames(totals_h3) = c('Year', colname)
        inc_h1 = merge(inc_h1, totals_h1, by='Year')
        inc_h3 = merge(inc_h3, totals_h3, by='Year')
    } 
    rownames(inc_h1) = inc_h1$Year
    rownames(inc_h3) = inc_h3$Year
    inc_h1 = inc_h1 %>% select(-Year)
    inc_h3 = inc_h3 %>% select(-Year)
    return(list(h1_inc=inc_h1, h3_inc=inc_h3))
}

# Output folder for results
result_folder="Compare_Data/" 

# Demography data by age, expects columns: Age, Year, pop_frac (the fraction of the total population that consists of people of a given age in a given year)
dem_data_raw = read.csv("../Data/Census_processed/AZ_2000-2017_raw.csv")

# Incidence data by age, expects columns: Age, H1 (number of cases), H3 (number of cases), Season
incidence_data_raw = read.csv("../Data/Processed_data/AZ_seasonal_linelist.csv", check.names = FALSE, stringsAsFactors = FALSE)
names(incidence_data_raw) = c('ssn', 'birthyear', 'subtype', 'year1', 'Season', 'Age')

# Setting age groups
age_groups <- rbind(c(0, 4),
                    c(5, 9),
                    c(10, 14),
                    c(15, 19),
                    c(20, 24),
                    c(25, 29),
                    c(30, 34),
                    c(35, 39),
                    c(40, 44),
                    c(45, 49),
                    c(50, 54),
                    c(55, 59),
                    c(60, 64),
                    c(65, 69),
                    c(70, 74),
                    c(75, 79),
                    c(80, 100))
seasons = c(2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015)
colnames(age_groups) <- c("LB","UB")
group_names <- paste(age_groups[,'LB'], age_groups[, 'UB'], sep='-')
rownames(age_groups) <- group_names


# Refactors demography data into fractions of population in each age group in each season
demo_group = calc_expected_dem_AG(dem_data_raw, age_groups, seasons)

# Kind of hacky thing where we assume that the population distribution in the pandemic 
# is the same as the population distribution in 2009
inc_data = calc_inc_AG(incidence_data_raw, age_groups, seasons)

demo_group = demo_group[order(row.names(demo_group)), ]
inc_data$h1_inc = inc_data$h1_inc[order(row.names(inc_data$h1_inc)), ]
inc_data$h3_inc = inc_data$h3_inc[order(row.names(inc_data$h3_inc)), ]






## Repeat for INSIGHT data
seasons = 2009:2017
dem_INSIGHT_raw = read.csv('../Data/Census_processed/INSIGHT_2009-2017_raw.csv', stringsAsFactors = FALSE)
dem_INSIGHT_raw$Country[which(dem_INSIGHT_raw$Country == 'UnitedStates')] = 'USA'
dem_INSIGHT_raw$Country[which(dem_INSIGHT_raw$Country == 'UnitedKingdom')] = 'UK'
countries = unique(dem_INSIGHT_raw$Country)[-1]

age_groups <- rbind(c(20, 24),
                    c(25, 29),
                    c(30, 34),
                    c(35, 39),
                    c(40, 44),
                    c(45, 49),
                    c(50, 54),
                    c(55, 59),
                    c(60, 64),
                    c(65, 69),
                    c(70, 74),
                    c(75, 79),
                    c(80, 90))
colnames(age_groups) <- c("LB","UB")
group_names <- paste(age_groups[,'LB'], age_groups[, 'UB'], sep='-')
rownames(age_groups) <- group_names


## By year (season)
calc_expected_dem_AG <- function(dem_data,
                                 age_groups,
                                 seasons,
                                 country){
  #######################
  #
  # Given demography data and age groups, 
  # this function calculates the proportion
  # of the population in each age group
  #
  #######################
  
  filtered_dem_data = dem_data %>% filter(Year %in% seasons)
  all_fracs = data.frame(Year=seasons)
  
  for (i in 1:nrow(age_groups)) {
    lb = age_groups[i, 'LB']
    ub = age_groups[i, 'UB']
    colname = paste(lb, ub, sep='-')
    filtered_dem_data %>% 
      group_by(Year) %>% 
      filter(Age <= ub & Age >=lb) %>% 
      summarize(temp=sum(pop_frac)) -> fracs
    colnames(fracs) = c('Year', colname)
    all_fracs = merge(all_fracs, fracs, by='Year')
  }
  
  rownames(all_fracs) = paste(country, all_fracs$Year, sep = '_')
  all_fracs = all_fracs %>% select(-Year) 
  return(all_fracs / rowSums(all_fracs))
}



dem_by_country = function(country){
  calc_expected_dem_AG(dem_data = subset(dem_INSIGHT_raw, Country == country), age_groups, seasons, country)
}
dem_tab = dem_by_country(countries[1])
for(ii in 2:length(countries)){
  dem_tab = rbind(dem_tab, dem_by_country(countries[ii]))
}


## Load INSIGHT incidence data
load('../Data/Processed_data/INSIGHT002_inc_tables.RData')
## H1.master and H3.master loaded
ages = as.numeric(colnames(H1.master))
group = c(0,0, rep(1:13, each = 5), rep(13, 6))


## Summarize by country and season into age groups
H1.summary = sapply(1:13, FUN = function(gg){rowSums(H1.master[,group == gg])})
H3.summary = sapply(1:13, FUN = function(gg){rowSums(H3.master[,group == gg])})
colnames(H1.summary) = colnames(H3.summary) = paste(age_groups[,1], age_groups[,2], sep = "-")




## Align row names for demography and incidence

## Get universal rownames from country sesason names from incidnece data
universal_rn = gsub(pattern = '(\\w+_)\\w\\w.\\d?\\d?.?(\\d\\d)', replacement = '\\120\\2', x = rownames(H1.summary))
## Copy dem_tab rows into corresponding positions
dem_tab = dem_tab[universal_rn, ]
## Overwrite universal row names with seasons
rownames(dem_tab) = rownames(H1.summary)


## Aggregate by country
# extract country names from summary rows
summary.rns = gsub(pattern = '(\\w+)_\\wH.\\d\\d.?\\d?\\d?', replacement = '\\1', x = rownames(H1.master))
H1.country = t(sapply(countries, FUN = function(cc) colSums(H1.summary[which(summary.rns == cc), ]))) # Take colsums across all rows belonging to that country
H3.country = t(sapply(countries, FUN = function(cc) colSums(H3.summary[which(summary.rns == cc), ])))
# Take mean of demographic age distribution across all seasons. Weight seasons proportional to toatl cases observed
weighted = dem_tab*rowSums(H1.summary)
dem_country_H1 = t(sapply(countries, FUN = function(cc) colSums(weighted[which(summary.rns == cc), ])))
weighted = dem_tab*rowSums(H3.summary)
dem_country_H3 = t(sapply(countries, FUN = function(cc) colSums(weighted[which(summary.rns == cc), ])))




## Aggregate by season
# extract season names from summary rows
summary.ssns = gsub(pattern = '\\w+_(\\wH.\\d\\d.?\\d?\\d?)', replacement = '\\1', x = rownames(H1.master))
ssns = unique(summary.ssns)
H1.season = t(sapply(ssns, FUN = function(ss) colSums(H1.summary[which(summary.ssns == ss), ]))) # Take colsums across all rows belonging to that season
H3.season = t(sapply(ssns, FUN = function(ss) colSums(H3.summary[which(summary.ssns == ss), ])))
# Take mean of demographic age distribution across all seasons. Weight seasons proportional to toatl cases observed
weighted = dem_tab*rowSums(H1.summary)
dem_season_H1 = t(sapply(ssns, FUN = function(ss) colSums(weighted[which(summary.ssns == ss), ])))
weighted = dem_tab*rowSums(H3.summary)
dem_season_H3 = t(sapply(ssns, FUN = function(ss) colSums(weighted[which(summary.ssns == ss), ])))





## Remove rows from country-season summary where fewer than 90 cases observed
useH1 = which(rowSums(H1.summary)>90)
useH3 = which(rowSums(H3.summary)>90)
H1.summary = H1.summary[useH1, ]
H3.summary = H3.summary[useH3, ]
dem_tabH1 = dem_tab[useH1, ]
dem_tabH3 = dem_tab[useH3, ]
rm(dem_tab)
