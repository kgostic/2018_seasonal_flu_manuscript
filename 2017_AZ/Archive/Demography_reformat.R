## Import census data
# Source: https://www.census.gov/population/international/data/idb/region.php?N=%20Results%20&T=10&A=separate&RT=0&Y=2013&R=-1&C=ID
USA.raw = read.csv('census_data_USA_1980_2020.csv') # No data available before 1980


## Reformat into a matrix with a single year on each row and a single age from 0:100 on each column
years = sort(unique(USA.raw$Year))
raw.formatted.data = matrix(NA, nrow = length(years), ncol = 89, dimnames = list(years, c(as.character(0:84), '85-89', '90-94', '95-99', '100+')))
# Fill in formatted data matrix
for(yy in 1:length(years)){
  raw.formatted.data[yy, ] = subset(USA.raw, Year == years[yy])$Both.Sexes.Population
}

# Disaggregate older ages
aggregate.ages = formatted.data[, 86:88]
formatted.data = cbind(raw.formatted.data[,1:85], matrix(apply(aggregate.ages, 2, function(x) rep(x, 5)), nrow = nrow(formatted.data))/5, raw.formatted.data[,89])
colnames(formatted.data) = 0:100

# Save the file
write.csv(x = formatted.data, file = 'USA_demography_1980_2020_formatted.csv')
