##### Plot a figure to describe hypotheses and modeling appropach
rm(list = ls())
source('00-Import_FLU002_-for-multinomial.R')

## OUTPUTS
outfile = 'figures/Fig1_hypotheses.pdf'



{
pdf(outfile)
## Set the layout
layout(mat = matrix(c(0,1,1, 9,3,6, 2,4,7, 10, 5, 8)+1, nrow = 4,  byrow = T), heights = c(1.5, 1,1,1))
#layout.show(10)



## Define a function to make transparent colors
tns = function(col.in, aa = .5){
  cc = col2rgb(col.in)
  new.col = rgb(red = cc[1], green = cc[2], blue = cc[3], 255*aa, maxColorValue = 255)
  new.col
}


######################################
## A. Plot imprinting reconstructions
######################################
par(mar = c(6,4,4,1))
xx = barplot(rbind(proH1.master[85,], proH2.master[85,], proH3.master[85,]), col = c('dodgerblue', 'lightskyblue', 'firebrick'), border = NA, space = 0, axes = T, xlab = '', ylab = '', main = 'Population-level imprinting')
axis(side = 1, line = 0, at = xx[seq(1, 75, by = 5)], labels = NA)
mtext(text = 'age', side = 1, line = 1.75, cex = .7)
mtext(text = 'prob imprinting to subtype', side = 2, line = 1.75, cex = .7)
legend(x = 0, y = 1.22, legend = c('H1N1', 'H2N2', 'H3N2'), col = c('dodgerblue', 'lightskyblue', 'firebrick'), bty = 'n', ncol = 3, xpd = NA, pch = 15)

## A. label
mtext(text = 'A', side = 3, line = 2, at = -12, font = 2)



######################################
## B. Plot Hypothesis table
######################################
par(mar = c(5,1,3,2))
plot.new()
plot.window(xlim = c(1.4, 8.75), ylim = c(0, 5.9))
abline(h = c(1,2,3))
segments(x0 = c(4,5,7,8), y0 = rep(0,6), y1 = rep(5, 6), lty = 3) # Dashed vertical lines
segments(x0 = 6, y0 = 0, y1 = 5) # Solid vertical line
## Add seasonal challenge boxes
polygon(x = c(3,3,6,6), y = c(6,7,7,6)-1, col = tns('dodgerblue'))
text(x = 4.5, y = 5.5, 'H1N1 challenge', font = 2)
polygon(x = c(3,3,6,6)+3, y = c(6,7,7,6)-1, col = tns('firebrick'))
text(x = 4.5+3, y = 5.5, 'H3N2 challenge', font = 2)

## Add Imprinting status row
#text(x = 2.25, y = 4, 'Imprinting\nhypothesis', font = 2)
text(x = 3.5, y = 4, 'H1N1\nimprinted', font = 1)
text(x = 4.5, y = 4, 'H2N2\nimprinted', font = 1)
text(x = 5.5, y = 4, 'H3N2\nimprinted', font = 1)
text(x = 6.5, y = 4, 'H1N1\nimprinted', font = 1)
text(x = 7.5, y = 4, 'H2N2\nimprinted', font = 1)
text(x = 8.5, y = 4, 'H3N2\nimprinted', font = 1)


## Add HA group row
text(x = 2, y = 2.5, 'HA group level', font = 2)
polygon(x = c(3,4,4,3), y = c(2,2,3,3), border = NA, col = tns('forestgreen')); text(x = 3.5, y = 2.5, 'protected') 
polygon(x = c(3,4,4,3)+1, y = c(2,2,3,3), border = NA, col = tns('forestgreen')); text(x = 4.5, y = 2.5, 'protected')
text(x = 5.5, y = 2.5, '-')
text(x = 6.5, y = 2.5, '-')
text(x = 7.5, y = 2.5, '-')
polygon(x = c(3,4,4,3)+5, y = c(2,2,3,3), border = NA, col = tns('forestgreen')); text(x = 8.5, y = 2.5, 'protected')

## Add HA subtype row
text(x = 2, y = 1.5, 'HA subtype level', font = 2)
polygon(x = c(3,4,4,3), y = c(2,2,3,3)-1, border = NA, col = tns('forestgreen')); text(x = 3.5, y = 1.5, 'protected')
text(x = 4.5, y = 1.5, '-')
text(x = 5.5, y = 1.5, '-')
text(x = 6.5, y = 1.5, '-')
text(x = 7.5, y = 1.5, '-')
polygon(x = c(3,4,4,3)+5, y = c(2,2,3,3)-1, border = NA, col = tns('forestgreen')); text(x = 8.5, y = 1.5, 'protected')

## Add NA subtype row
text(x = 2, y = 0.5, 'NA subtype level', font = 2)
polygon(x = c(3,4,4,3), y = c(2,2,3,3)-2, border = NA, col = tns('forestgreen')); text(x = 3.5, y = 0.5, 'protected')
text(x = 4.5, y = 0.5, '-')
text(x = 5.5, y = 0.5, '-')
text(x = 6.5, y = 0.5, '-')
polygon(x = c(3,4,4,3)+4, y = c(2,2,3,3)-2, border = NA, col = tns('forestgreen')); text(x = 7.5, y = 0.5, 'protected')
polygon(x = c(3,4,4,3)+5, y = c(2,2,3,3)-2, border = NA, col = tns('forestgreen')); text(x = 8.5, y = 0.5, 'protected')

## B. label
mtext(text = 'B', side = 3, line = 1, at = 1.5, font = 2)






######################################
## C. Age-specific predictions
######################################
## These are made up, just for illustration
par(mar = c(2,2.5,2,1)+1)
age.vec = rep(c(3.3, 3.5, 3.7, 3.4, 2.9, 2.4, 2.3, 1.9, 1.8, 1.7), c(rep(7,9), 10))
age.vec = age.vec/sum(age.vec)
plot(18:90, age.vec, xlab = '', ylab = '', main = '', ylim = c(0, .025))
mtext(text = 'age', side = 1, line = 1.9, cex = .7)
mtext(text = 'fraction of cases', side = 2, line = 1.9, cex = .7)

## C. label
mtext(text = 'C', side = 3, line = 1, at = 1.5, font = 2)









######################################
## D-F. Imprinting predictions
######################################
rownames(prog1.master)[85] ## Use Thailand NH.16.17 as the example
plot(18:90, proH1.master[85,], col = 'dodgerblue', xlab = '', ylab = '', main = '')
points(18:90, proH3.master[85,], col = 'firebrick')
mtext(text = 'age', side = 1, line = 1.9, cex = .7)
mtext(text = 'P(imprint prot.)', side = 2, line = 1.9, cex = .7)
mtext(text = 'HA subtype-level', side = 3, line = .2, cex = .7, font = 2)

## Legend
legend(x = -5, y = 1.55+.2, legend = c('H1N1 protection', 'H3N2 protection'), col = c('dodgerblue', 'firebrick'), bty = 'n', ncol = 3, xpd = NA, pch = 15)

## Section header
polygon(x = c(-5, -5, 100, 100), y = c(1.3, 1.6, 1.6, 1.3)+.4, col = 'gray80', xpd = NA, border = NA)
text(45, 1.65+.2, 'Imprinting effects', xpd = NA, font = 2, cex = 1.2)

## D. label
mtext(text = 'D', side = 3, line = 1, at = 1.5, font = 2)


## HA subtype-level imprinting
plot(18:90, prog1.master[85,], col = 'dodgerblue', xlab = '', ylab = '', main = '')
points(18:90, prog2.master[85,], col = 'firebrick')
mtext(text = 'age', side = 1, line = 1.9, cex = .7)
mtext(text = 'P(imprint prot.)', side = 2, line = 1.9, cex = .7)
mtext(text = 'HA group-level', side = 3, line = .2, cex = .7, font = 2)
## E. label
mtext(text = 'E', side = 3, line = 1, at = 1.5, font = 2)


## NA subtype-level imprinting
plot(18:90, proN1.master[85,], col = 'dodgerblue', xlab = '', ylab = '', main = '')
points(18:90, proN2.master[85,], col = 'firebrick')
mtext(text = 'age', side = 1, line = 1.9, cex = .7)
mtext(text = 'P(imprint prot.)', side = 2, line = 1.9, cex = .7)
mtext(text = 'NA subtype-level', side = 3, line = .2, cex = .7, font = 2)
## F. label
mtext(text = 'F', side = 3, line = 1, at = 1.5, font = 2)





######################################
## G-I. Predicted risk
######################################
## HA subtype prediction
## Write a function to predict overall risk, given assumed Hm value (relative risk of infection given imprinting protection)
rr = function(Hm, pro.vec){
  pred = (Hm*pro.vec+(1-pro.vec))*age.vec
  pred/sum(pred)
}
Hm = .75 # Set the assumed relative risk of infection given imprinting protection
plot(18:90, rr(Hm, proH1.master[85,]), col = 'dodgerblue', xlab = '', ylab = '', main = '', ylim = c(0,.025))
points(18:90, rr(Hm, proH3.master[85,]), col = 'firebrick')
mtext(text = 'age', side = 1, line = 1.9, cex = .7)
mtext(text = 'fraction of cases', side = 2, line = 1.9, cex = .7)
mtext(text = 'Age + HA subtype-level imprinting', side = 3, line = .2, cex = .7, font = 2)


## Legend
legend(x = 20, y = .0436, legend = c('H1N1', 'H3N2'), col = c('dodgerblue', 'firebrick'), bty = 'n', ncol = 3, xpd = NA, pch = 15)

## Section header
polygon(x = c(-5, -5, 100, 100), y = c(.0424, .05, .05, .0424), col = 'gray80', xpd = NA, border = NA)
text(49, .046, 'Hypothetical model predictions', xpd = NA, font = 2, cex = 1.2)

## G. label
mtext(text = 'G', side = 3, line = 1, at = 1.5, font = 2)

## HA subtype-level imprinting
plot(18:90, rr(Hm, prog1.master[85,]), col = 'dodgerblue', xlab = '', ylab = '', main = '', ylim = c(0,.025))
points(18:90, rr(Hm, prog2.master[85,]), col = 'firebrick')
mtext(text = 'age', side = 1, line = 1.9, cex = .7)
mtext(text = 'fraction of cases', side = 2, line = 1.9, cex = .7)
mtext(text = 'Age + HA group-level imprinting', side = 3, line = .2, cex = .7, font = 2)
## H. label
mtext(text = 'H', side = 3, line = 1, at = 1.5, font = 2)

## NA subtype-level imprinting
plot(18:90, rr(Hm, proN1.master[85,]), col = 'dodgerblue', xlab = '', ylab = '', main = '', ylim = c(0,.025))
points(18:90, rr(Hm, proN2.master[85,]), col = 'firebrick')
mtext(text = 'age', side = 1, line = 1.9, cex = .7)
mtext(text = 'fraction of cases', side = 2, line = 1.9, cex = .7)
mtext(text = 'Age + NA subtype-level imprinting', side = 3, line = .2, cex = .7, font = 2)
## I. label
mtext(text = 'I', side = 3, line = 1, at = 1.5, font = 2)


######################################
##  Age risk section header
######################################
plot.new()
plot.window(xlim = c(min(xx), max(xx)), ylim = c(0,1))
## Section header
polygon(x = c(-15, -15, 83, 83), y = c(1.3, 1.6, 1.6, 1.3)+.4, col = 'gray80', xpd = NA, border = NA)
text(33, 1.65+.2, 'Hypothetical age effects', xpd = NA, font = 2, cex = 1.2)
#text(87, 1.86, '+', cex = 2, xpd = NA)
#text(200, 1.86, '\u27a1', cex = 2, xpd = NA)
#plot(200, 1.86, pch=-0x279EL, xpd = NA)

dev.off()
}
