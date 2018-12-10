months = c('Aug', 'Sept', 'Oct', 'Nov', 'Dec', 'Jan', 'Feb', 'March', 'April', 'May', 'June', 'July')

years = 5

dev.off()
plot.new()
plot.window(c(0, 12*years), c(0, 1))
axis(1, at = 1:(12*years), labels = rep(months, years), las = 2)
sapply(X = 12*((1:years)-1)+6, function(s) lines(x = c(s,s), y = c(-.1, .1), col = 'blue', cex = 3))
sapply(X = 12*((1:years)-1)+3, function(s) polygon(x = c(s, s, s+6, s+6), y = c(0, .1, .1, 0)+.05, col = 'darkslategray3', border = NA))

