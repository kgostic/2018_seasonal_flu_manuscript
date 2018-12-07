library(tidyverse)

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
            summarize(total=sum(H1)) -> totals_h1
        
        filtered_inc_data %>% 
            group_by(Season) %>% 
            filter(Age <= ub & Age >=lb) %>% 
            summarize(total=sum(H3)) -> totals_h3
        
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
result_folder="../results/MESA_eligible_separate/" 

# Demography data by age, expects columns: Age, Year, pop_frac (the fraction of the total population that consists of people of a given age in a given year)
dem_data_raw = read.csv("../../data/processed_data/MESA_demography_eligible_core.csv", check.names = FALSE)

# Incidence data by age, expects columns: Age, H1 (number of cases), H3 (number of cases), Season
incidence_data_raw = read.csv("../../POMP/separate_vaccinated_class_pandemic/filtered_eligible_observed_by_age.csv", check.names = FALSE, stringsAsFactors = FALSE)

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
seasons = c(2008, 2009, 2009.5, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018)
colnames(age_groups) <- c("LB","UB")
group_names <- paste(age_groups[,'LB'], age_groups[, 'UB'], sep='-')
rownames(age_groups) <- group_names

# Refactors demography data into fractions of population in each age group in each season
demo_group = calc_expected_dem_AG(dem_data_raw, age_groups, seasons)

# Kind of hacky thing where we assume that the population distribution in the pandemic 
# is the same as the population distribution in 2009
demo_group = rbind(demo_group, '2009.5'=demo_group['2009', ])
inc_data = calc_inc_AG(incidence_data_raw, age_groups, seasons)

demo_group = demo_group[order(row.names(demo_group)), ]
inc_data$h1_inc = inc_Dat[order(row.names(inc_data$h1_inc)), ]
inc_data$h3_inc = inc_Dat[order(row.names(inc_data$h3_inc)), ]