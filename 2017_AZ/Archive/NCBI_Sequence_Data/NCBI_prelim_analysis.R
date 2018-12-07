## Import all H1N1 sequences
setwd('~/Dropbox/R/2017_seasonal_flu/NCBI_Sequence_Data/')
H1.seq = read.table('~/Dropbox/R/2017_seasonal_flu/NCBI_Sequence_Data/H1N1_combined_AGEINFO.txt', header = TRUE, sep = '\t', stringsAsFactors = FALSE)

valid = H1.seq$Year >80
H1.seq$Year[valid] = H1.seq$Year[valid] + 1900
valid = H1.seq$Year < 20
H1.seq$Year[valid] = H1.seq$Year[valid] + 2000

## Remove duplicate accession
uID = unique(H1.seq$Accession)
keep = check = 1
for(ii in 1:length(uID)){
  valid = which(H1.seq$Accession == uID[ii])
  ages = unique(H1.seq$Age[valid])
  if(length(ages) == 1){
  keep[ii] = valid[1]
  }else{
    check[ii] = uID[ii]
  }
}
keep = keep[-1]; check = check[-1]; check ## If check is not empty, go back and check those accession numbers
H1.seq = H1.seq[keep,]

table(H1.seq$Year)

H1.seq$BirthYear = H1.seq$Year - H1.seq$Age
yrs = sort(unique(H1.seq$Year))

pdf('H1N1_NCBI_dat.pdf')
par(mfrow = c(2,2))
for(ii in 1:length(yrs)){
  dat = subset(H1.seq, Year == yrs[ii])
  hist(x = dat$BirthYear, breaks = 1900:2017, main = paste(yrs[ii], ', N = ', nrow(dat)), xlab = 'Birth Year')
}
hist(x = H1.seq$BirthYear, breaks = 1900:2017, main = paste('All years, N = ', nrow(H1.seq)), xlab = 'Birth Year')
dev.off()







## Import all H3N2 sequences
setwd('~/Dropbox/R/2017_seasonal_flu/NCBI_Sequence_Data/')
H3.seq = read.table('~/Dropbox/R/2017_seasonal_flu/NCBI_Sequence_Data/H3N2_combined_ALLAGES.txt', header = TRUE, sep = ',', stringsAsFactors = FALSE)

table(H3.seq$Year)

H3.seq$BirthYear = H3.seq$Year - H3.seq$Age
yrs = sort(unique(H3.seq$Year))

## Remove duplicates
uID = unique(H3.seq$Acession)
keep = 1
for(ii in 1:length(uID)){
  valid = which(H3.seq$Acession == uID[ii])
  keep[ii] = valid[1]
}
H3.seq = H3.seq[keep,]

pdf('H3N2_NCBI_dat.pdf')
par(mfrow = c(2,2))
for(ii in 1:length(yrs)){
  dat = subset(H3.seq, Year == yrs[ii])
  hist(x = dat$BirthYear, breaks = 1870:2017, main = paste(yrs[ii], ', N = ', nrow(dat)), xlab = 'Birth Year')
}
hist(x = H3.seq$BirthYear, breaks = 1870:2017, main = paste('All years, N = ', nrow(H1.seq)), xlab = 'Birth Year')
dev.off()


## Save data
write.csv(x = H1.seq, file = 'NCBI_H1_Data.csv')
write.csv(x = H3.seq, file = 'NCBI_H3_Data.csv')
