rm(list=ls())
#install.packages("tidyverse")
suppressPackageStartupMessages(library(tidyverse))
setwd("/Users/akshitha/Documents/GitHub/AdvanceStatistics-6310/project/pythondataprocessor")
myfile1=read.csv("finalpatients.csv",header = TRUE, sep=",", check.names=FALSE)
myfile2= read.csv("finalgenevalues.csv",header= TRUE, sep=",",check.names=FALSE)

age<- myfile1$AGE

# Performing linear regression for Age against all the genes.

Pvalues <- vector()
age<- myfile1$AGE

for ( i in 1:nrow(myfile2)){
  #print(i)
  mydata <- as.numeric(myfile2[i,4:93])
  #print(length(mydata))
  mylm<- lm(mydata ~ age)
  Pvalues[i] <- anova(mylm)$"Pr(>F)"[1]
  
}
Pvalues
hist(Pvalues,main = "Linear Regression for Genes  Vs Age",breaks =30,col = "red", xlab="P-values")
AdjPval_Lme <- sum(p.adjust(Pvalues,method = "fdr") <= 0.10)
AdjPval_Lme

# Performing T.test for Gender vs expression levels
Tanalysis<-c()
y=myfile1$SEX
for ( i in 1:nrow(myfile2)){
  x=as.numeric(myfile2[i,4:93])
  #print(i)
  Tanalysis[i]= t.test(x~y)$p.value
}
hist(Tanalysis,main = "T test for Genes  Vs Gender",breaks =30,col = "blue", xlab="P-values")
AdjPval_Gender <- sum(p.adjust(Tanalysis,method = "fdr") <= 0.1)
AdjPval_Gender

# Performing T.test for Drinking status across all the genes.
TanalysisVsDrink<-c()
y=myfile1$drinkingstatus
for ( i in 1:nrow(myfile2)){
  x=as.numeric(myfile2[i,4:93])
  #print(i)
  TanalysisVsDrink[i]= t.test(x~y)$p.value
}
hist(TanalysisVsDrink,main = "T test for Genes  Vs Drinking status",breaks =30,col = "Orange", xlab="P-values")
AdjPval_D<- sum(p.adjust(TanalysisVsDrink,method = "fdr") <= 0.1)
AdjPval_D

# Performing T.test for Smoking status across all the genes.
TanalysisVsSmoke<-c()
y=myfile1$smokingstatus
for ( i in 1:nrow(myfile2)){
  x=as.numeric(myfile2[i,4:93])
  #print(i)
  TanalysisVsSmoke[i]= t.test(x~y)$p.value
}
hist(TanalysisVsSmoke,main = "T test for Genes Vs Smoking status",breaks =30,col = "magenta", xlab="P-values")
AdjPval_S<- sum(p.adjust(TanalysisVsSmoke,method = "fdr") <= 0.1)
AdjPval_S
