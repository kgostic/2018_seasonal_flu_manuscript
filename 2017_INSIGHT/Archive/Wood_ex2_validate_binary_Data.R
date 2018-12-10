n <- 100; m <- 1
x <- runif(n)
lp <- 3*x-1
mu <- binomial()$linkinv(lp)
y <- rbinom(1:n,m,mu)
par(mfrow=c(2,2))
plot(glm(y/m ~ x,family=binomial,weights=rep(m,n)))


## Write colde to take a glm fitted to binary data. Simulate new binary data given model results. Refit the model using the resulting data and extract residuals from the resulting fits.
fit = glm(y ~ x, family = binomial)
probs = fitted(fit)
sim.dat = sapply(probs, FUN = function(pp) rbinom(n = 1, size = 1, prob = pp))
fit2 = glm(sim.dat ~ x, family = binomial)
plot(fit2)
rsd = resid(fit2)
plot(sort(rsd),(1:length(rsd)-.5)/length(rsd))


## Turn the above code into a for loop that repeats the simulation 100 times and stores the sorted residuals into a matrix
fit0 = fit
one.run = function(){
  probs = fitted(fit0)
  sim.dat = sapply(probs, FUN = function(pp) rbinom(n = 1, size = 1, prob = pp))
  fit0 = glm(sim.dat ~ x, family = binomial)
  rsd = resid(fit0)
  return(sort(rsd))
}

# Replicate 1000 residual simulations
raw.mat = replicate(1000, expr = one.run())
# Sort each column (quantile)
sorted.mat = apply(raw.mat, MARGIN = 1, FUN = sort)
plot(sorted.mat[500,],(1:length(rsd)-.5)/length(rsd), type = 'l')
lines(sorted.mat[50,], (1:length(rsd)-.5)/length(rsd), lty = 2)
lines(sorted.mat[950,], (1:length(rsd)-.5)/length(rsd), lty = 2)
for(ii in c(1:49, 951:1000)){
  lines(sorted.mat[ii,], (1:length(rsd)-.5)/length(rsd), lty = 3)
}

