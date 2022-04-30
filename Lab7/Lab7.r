rm(list=ls())
setwd("/Users/akshitha/Documents/Github/Advanced-Statistics-for-Bioinformatics/lab7")
myfile = read.table("prePostPhylum.txt",header =TRUE, sep="\t")
myfile
numCols = ncol(myfile)
#numCols
ColClasses = c(rep("character",4),rep("numeric",numCols-4))
#ColClasses
PhylaData = myfile[,5:10]
#PhylaData
cageData =myfile$cage
#cageData
genotypeData = myfile$genotype
#genotypeData
TimeData = myfile$time
#TimeData
PCOA = princomp(PhylaData)

components =summary(PCOA)
components
PCAscores=PCOA$scores
PCA1 =PCAscores[,1]

PCA2= PCAscores[,2]
# Question 2
# (2)	Graph PCA1 vs. PCA2.  Make three versions of the graph.  One colored by genotype,
# one colored by cage and one colored by timepoint (pre-vs-post)

#install.packages('ggfortify')
library(ggfortify)
# PCA 1 vs PCA2 with respect to genotype
genotype_plot = autoplot(PCOA,data = myfile, colour = 'genotype', main ="PCA 1 Vs PCA 2 by Genotype")
genotype_plot
# uncomment the below statement if the graphs does'nt work and recompile the graph code.
#dev.off()
# PCA 1 vs PCA2 with respect to cage
cage_plot = autoplot(PCOA,
                 data = myfile,
                 colour = 'cage',main ="PCA 1 Vs PCA 2 by cage")
cage_plot
timepoint_plot = autoplot(PCOA,
                          data = myfile,
                          colour = 'time', main ="PCA 1 Vs PCA 2 by Time")
timepoint_plot

# Question 3
#(3)	Fill in the following table for p-values testing the null hypothesis for PCA 1 and and 2.  
# For cage, use a way one-ANOVA.  For genotype and timepoint (“pre” vs “post”) use a t-test
# Which variable seems to be most associated with the first PCA axis?  
#Which variable is most associated with the second PCA axis?  Does cage seem to be having 
#an effect on these data?

#T.Test by genotype
#genotype and PCA1

PCA1_genotype = t.test(PCA1~genotypeData)
PCA1_genotype
PCA2_genotype = t.test(PCA2~genotypeData)
PCA2_genotype

# t.test for time point
PCA1_timepoint = t.test(PCA1~TimeData)
PCA1_timepoint
PCA2_timepoint = t.test(PCA2~TimeData)
PCA2_timepoint

# Anova for cage
PCA1_cage = aov(PCA1~cageData, data = myfile)
PCA1_cage
summary(PCA1_cage)
PCA2_cage =aov(PCA2~cageData, data = myfile)
PCA2_cage
summary(PCA2_cage)
# Question 4 A
postTimeData = myfile[myfile$time=="POST",]
postTimeData
PostphylaData= postTimeData[,5:10]
PostphylaData
PostcageData = postTimeData$cage
PostgenotypeData = postTimeData$genotype
fileDataframe = data.frame(PostphylaData,PostcageData,PostgenotypeData)
fileDataframe
# Linear Model
lmPvals = vector()
lmePvals =vector()
adj_Pvalues_lme = vector()
GlsPvals = vector()
rhoPvals= vector()
adj_Pvalues_gls = vector()
#CorrelationCoeff= vector()
phylaName = colnames(fileDataframe)[1:6]
phylaName
colors=c("blue","red","green","yellow","orange","cyan")
library(nlme)
for (i in 1:ncol(fileDataframe[1:6])){
  # print(fileDataframe[,i])
  # Linear Model
  mylm = lm(fileDataframe[,i]~fileDataframe$PostcageData)
  lmPvals[i]<- anova(mylm)$"Pr(>F)"[1] 
  
  # Mixed model

  PostphylaData = as.numeric(fileDataframe[,i])
  mylme= lme(PostphylaData~PostgenotypeData, random = ~1| PostcageData, data=fileDataframe)
  lmePvals[PostphylaData] =unclass(summary(mylme))$tTable[2,5]
  adj_Pvalues_lme = p.adjust(lmePvals[PostphylaData],method ="fdr")
  
  # gls model
  myGls= gls(PostphylaData~PostgenotypeData, correlation = corCompSymm(form=~1|PostcageData), method ="REML", data= fileDataframe)
  GlsPvals[PostphylaData] = anova(myGls)$ "p-value"[2]
  rhoPvals[PostphylaData]= coef(myGls$modelStruct[1]$corStruct,unconstrained=FALSE)[[1]]
  
  adj_Pvalues_gls = p.adjust(GlsPvals[PostphylaData],method ="fdr")
  #par(mfrow = c(1,2)) 
layout.matrix <- matrix(c(4,3,2,1,0,3), nrow=3, ncol=2)
plot(fileDataframe[,i] ~ fileDataframe$PostcageData, main = paste(c("P-value:",lmPvals[i],phylaName[i]),sep="\n"),sub= paste("Rho_value=",rhoPvals[PostphylaData]),las=2,cex.axis= 0.6, xlab = "", ylab = "Relative Abundance",col=colors[i])
stripchart(fileDataframe[,i] ~ fileDataframe$PostcageData, data = fileDataframe,vertical = TRUE, pch = 19, add=TRUE)

}

# Except for the phylum- Actinobacteria, all other phylas exhibits cage effect.
# by getting the pvalues from above,

lmPvals 
#lmePvals
#adj_Pvalues_lme 
#GlsPvals 
#rhoPvals
#adj_Pvalues_gls 

#par(mfrow = c(1,1))
# Question 4B
hist(lmPvals,main = "Linear Model-Phyla P-values Vs Cage",breaks =30,col = "red", xlab="P-values")
hist(lmePvals,main = "Mixed Linear Model-Phyla P-values Vs Cage",breaks =30,col = "green",xlab="P-values")
hist(GlsPvals,main = "Gls-Mixed Linear Model-Phyla vs Cage",breaks =30,col = "Orange",xlab="P-values")

# with 10% FDR
Gls_FDR <-sum(GlsPvals <= 0.1) 
lme_FDR <- sum(lmePvals <= 0.1)
#Gls_FDR
AdjPval_Gls <- sum(p.adjust(GlsPvals,method = "fdr")<= 0.1)
AdjPval_Gls
AdjPval_Lme <- sum(p.adjust(lmePvals,method = "fdr") <= 0.1)
AdjPval_Lme
#
#From the above graph we can get all the rho values are as follows
# Tenericutes = 0.6075008, Verrucomicrobia = 0.64679876, Bacteroidetes = 0.5806221, Actinobacteria= -0.04015915
# Firmicutes = 0.56426547, Proteobacteria = 0.4274917
# Except for Actinobacteria, remaining all seem to look same. Based on the rho value of Actinobacteria, we can say that there is no effect of genotype + cage.


