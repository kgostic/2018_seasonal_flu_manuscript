## Format results into a text table

rm(list = ls())
load('processed-data/AZ_CIs.RData')
load('processed-data/AZ_model_fits.RData')
load('processed-data/AZ_profiles.RData')

outfile1 = 'processed-data/AZ_results_table.csv'
outfile2 = 'processed-data/AIC_AZ.csv'
outfile3 = '../tables/Formatted_AZ_results.csv'

# Initialize matris
results = matrix(NA, nrow = 4, ncol = length(lk.AG$par)*3, dimnames = list(c('A', 'AN', 'AS', 'AG'),  paste(rep(names(lk.AG$par), each = 3), c('best', 'low', 'high'), sep = '_') ))

# Set column indices of best esimtates, low and high CIs
best = seq(1, 42, by = 3)
low = best+1
high = best+2

# Model A results
results[1,best[-c(1,2)]] = lk.A$par
results[1,low[-c(1,2)]] = A.CIs[1,]
results[1,high[-c(1,2)]] = A.CIs[2,]
# Model AN results
results[2,best] = lk.AN$par
results[2,low] = AN.CIs[1,]
results[2,high] = AN.CIs[2,]
# Model AS results
results[3,best] = lk.AS$par
results[3,low] = AS.CIs[1,]
results[3,high] = AS.CIs[2,]
# Model AG results
results[4,best] = lk.AG$par
results[4,low] = AG.CIs[1,]
results[4,high] = AG.CIs[2,]


results = t(round(results, 2))

write.csv(results, file = outfile1)
write.csv(del.AIC, file = outfile2)

## Write function to fomat table S1
fmt = function(rn){
  paste(results[rn,], " (", results[rn+1, ],"-", results[rn+2, ], ")", sep = "")
}

formatted = t(sapply(X = best, FUN = fmt))
colnames(formatted) = colnames(results)
rownames(formatted) = names(lk.AG$par)
write.csv(formatted, file = outfile3)

