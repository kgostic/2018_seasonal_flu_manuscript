
rm(list = ls())
load('Hm_H1_chain1_1900steps_2017-08-18.RData')


# Trace plots
par(mfrow = c(3, 2))
for(ii in 1:6){
  plot(1:ncol(ests), ests[ii, ], main = rownames(ests)[ii], type = 'l', col = 'blue')
}

# :(