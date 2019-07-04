##### Plot a figure to describe hypotheses and modeling appropach
rm(list = ls())
setwd('~/Dropbox/R/2018_seasonal_flu/2017_AZ/')
source('00-Inputs_multinomial.R')


## OUTPUTS
outfile = '../figures/Fig1_hypotheses.pdf'




pdf(outfile)
## Set the layout
layout(matrix(c(0,0,0,1,1,1, 9,9,3,3,6,6, 2,2,4,4,7,7, 10,10,5,5,8,8)+1, nrow = 4,  byrow = T), heights = c(1.5, 1,1,1))
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
h1_imprinted = (proH1.master[13,])
h2_imprinted = (proH2.master[13,])
h3_imprinted = (proH3.master[13,]) 
names(h1_imprinted)=names(h2_imprinted)=names(h3_imprinted)=2015-(as.numeric(colnames(proH1.master)))
naive = 1-h1_imprinted-h2_imprinted-h3_imprinted
xx = barplot(rbind(h1_imprinted, h2_imprinted, h3_imprinted, naive), col = c('dodgerblue', 'lightskyblue', 'firebrick', 'gray'), border = NA, space = 0, axes = T, xlab = '', ylab = '', main = '', xaxt = 'n')
axis(side = 1, line = 0, at = xx[seq(1, 96, by = 5)]-.5, labels = xx[seq(1, 96, by = 5)]-.5)
mtext(text = 'age in 2015', side = 1, line = 1.75, cex = .7)
mtext(text = 'prob imprinting to subtype', side = 2, line = 1.75, cex = .7)
legend(x = 15, y = 1.4, legend = c('H3N2', 'H2N2', 'H1N1', 'naive'), col = c('firebrick', 'lightskyblue', 'dodgerblue', 'gray'),  ncol = 4, xpd = NA, pch = 15)

## A. label
mtext(text = 'A', side = 3, line = 2, at = -12, font = 2)



######################################
## B. Plot Hypothesis table
######################################
par(mar = c(5,1,3,1))
plot.new()
plot.window(xlim = c(1.4, 6), ylim = c(0, 4))
abline(h = c(1,2,3))
segments(x0 = c(4,5,7,8), y0 = rep(0,7), y1 = rep(5, 6), lty = 3) # Dashed vertical lines
segments(x0 = 6, y0 = 0, y1 = 5) # Solid vertical line

## Add Imprinting status row
y.imp = 3.5
text(x = 2, y = y.imp, ('Child imprinted to'), font = 2, xpd = NA)
text(x = 3.5, y = y.imp, 'H1N1', font = 1)
text(x = 4.5, y = y.imp, 'H2N2', font = 1)
text(x = 5.5, y = y.imp, 'H3N2', font = 1)

## Add HA group row
text(x = 2, y = 2.5, 'HA group level', font = 2)
polygon(x = c(3,4,4,3), y = c(2,2,3,3), border = NA, col = tns('dodgerblue')); text(x = 3.5, y = 2.5, 'protected\nH1N1') 
polygon(x = c(3,4,4,3)+1, y = c(2,2,3,3), border = NA, col = tns('dodgerblue')); text(x = 4.5, y = 2.5, 'protected\nH1N1') 
polygon(x = c(3,4,4,3)+2, y = c(2,2,3,3), border = NA, col = tns('red')); text(x = 5.5, y = 2.5, 'protected\nH3N2') 


## Add HA subtype row
text(x = 2, y = 1.5, 'HA subtype level', font = 2)
polygon(x = c(3,4,4,3), y = c(2,2,3,3)-1, border = NA, col = tns('dodgerblue')); text(x = 3.5, y = 1.5, 'protected\nH1N1') 
text(x = 4.5, y = 1.5, 'no\nprotection') # No protection this row
polygon(x = c(3,4,4,3)+2, y = c(2,2,3,3)-1, border = NA, col = tns('red')); text(x = 5.5, y = 1.5, 'protected\nH3N2') 


## Add NA subtype row
text(x = 2, y = 0.5, 'NA subtype level', font = 2)
polygon(x = c(3,4,4,3), y = c(2,2,3,3)-2, border = NA, col = tns('dodgerblue')); text(x = 3.5, y = .5, 'protected\nH1N1') 
polygon(x = c(3,4,4,3)+1, y = c(2,2,3,3)-2, border = NA, col = tns('red')); text(x = 4.5, y = .5, 'protected\nH3N2') 
polygon(x = c(3,4,4,3)+2, y = c(2,2,3,3)-2, border = NA, col = tns('red')); text(x = 5.5, y = .5, 'protected\nH3N2') 


## B. label
mtext(text = 'B', side = 3, line = 1, at = 1.5, font = 2)





######################################
## C. Age-specific predictions
######################################
## These are made up, just for illustration
par(mar = c(2,2.5,1,1)+1)
age.vec = rep(c(4.1, 4.05, 3.9, 3.7, 3.5, 2.8, 2.6, 2.4, 2.3, 2.2, 1.9, 1.8, 1.7), c(5,6,rep(7,10), 17))
age.vec = age.vec/sum(age.vec)
plot(0:97, age.vec, xlab = '', ylab = '', main = '', ylim = c(0, .02))
mtext(text = 'age', side = 1, line = 1.9, cex = .7)
mtext(text = 'fraction of cases', side = 2, line = 1.9, cex = .7)

## C. label
mtext(text = 'C', side = 3, line = 1, at = 1.5, font = 2)









######################################
## D-F. Imprinting predictions
######################################
xx = -as.numeric(colnames(proH1.master)) # Convert from age to birth year for x axis
plot(xx, 1-h1_imprinted, col = 'dodgerblue', xlab = '', ylab = '', main = '', xaxt = 'n')
ax.lab = seq(xx[6], xx[96], by = 15)
axis(1, at = ax.lab, labels = -ax.lab)
points(xx, 1-h3_imprinted, col = 'firebrick')
mtext(text = 'birth year', side = 1, line = 1.9, cex = .7)
mtext(text = 'P(unprotected)', side = 2, line = 1.9, cex = .7)
mtext(text = 'HA subtype-level', side = 3, line = .2, cex = .7, font = 2)
## D. label
mtext(text = 'D', side = 3, line = 1, at = -2015, font = 2)
legend(x = -2010, y = 1.6, legend = c('H1N1 risk', 'H3N2 risk'), col = c('dodgerblue', 'firebrick'), bty = 'n', ncol = 3, xpd = NA, pch = 15)


## HA group-level imprinting
plot(xx, 1-(h1_imprinted+h2_imprinted), col = 'dodgerblue', xlab = '', ylab = '', main = '', xaxt = 'n')
axis(1, at = ax.lab, labels = -ax.lab)
points(xx, 1-h3_imprinted, col = 'firebrick')
mtext(text = 'birth year', side = 1, line = 1.9, cex = .7)
mtext(text = 'P(unprotected)', side = 2, line = 1.9, cex = .7)
mtext(text = 'HA group-level', side = 3, line = .2, cex = .7, font = 2)
## E. label
mtext(text = 'E', side = 3, line = 1, at = -2015, font = 2)


## NA subtype-level imprinting
plot(xx, 1-h1_imprinted, col = 'dodgerblue', xlab = '', ylab = '', main = '', xaxt = 'n')
axis(1, at = ax.lab, labels = -ax.lab)
points(xx, 1-(h2_imprinted+h3_imprinted), col = 'firebrick')
mtext(text = 'birth year', side = 1, line = 1.9, cex = .7)
mtext(text = 'P(unprotected)', side = 2, line = 1.9, cex = .7)
mtext(text = 'NA subtype-level', side = 3, line = .2, cex = .7, font = 2)
## F. label
mtext(text = 'F', side = 3, line = 1, at = -2015, font = 2)





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
plot(0:97, rr(Hm, h1_imprinted), col = 'dodgerblue', xlab = '', ylab = '', main = '', ylim = c(0,.02))
points(0:97, rr(Hm, h3_imprinted), col = 'firebrick')
mtext(text = 'age', side = 1, line = 1.9, cex = .7)
mtext(text = 'fraction of cases', side = 2, line = 1.9, cex = .7)
legend(x = 10, y = .0315, legend = c('H1N1 risk', 'H3N2 risk'), col = c('dodgerblue', 'firebrick'), bty = 'n', ncol = 3, xpd = NA, pch = 15)
## G. label
mtext(text = 'G', side = 3, line = 1, at = 1.5, font = 2)

## HA subtype-level imprinting
plot(0:97, rr(Hm,(h1_imprinted+h2_imprinted)), col = 'dodgerblue', xlab = '', ylab = '', main = '', ylim = c(0,.02))
points(0:97, rr(Hm, h3_imprinted), col = 'firebrick')
mtext(text = 'age', side = 1, line = 1.9, cex = .7)
mtext(text = 'fraction of cases', side = 2, line = 1.9, cex = .7)
mtext(text = 'H', side = 3, line = 1, at = 1.5, font = 2)

## NA subtype-level imprinting
plot(0:97, rr(Hm, h1_imprinted), col = 'dodgerblue', xlab = '', ylab = '', main = '', ylim = c(0,.02))
points(0:97, rr(Hm, (h2_imprinted+h3_imprinted)), col = 'firebrick')
mtext(text = 'age', side = 1, line = 1.9, cex = .7)
mtext(text = 'fraction of cases', side = 2, line = 1.9, cex = .7)
## I. label
mtext(text = 'I', side = 3, line = 1, at = 1.5, font = 2)


dev.off()

