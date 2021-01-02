setwd('F:/R/MARS')
library(earth)
a <- read.table("30_edger.txt",header=T, row.names = 1, sep = "\t" )   
fit <- earth(P30~., a,trace=4)  #Fit the model
summary(fit)   #summarize the fit
evimp(fit)   #summarize the importance of input variables
cat(format(fit))
