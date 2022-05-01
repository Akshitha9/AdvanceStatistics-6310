rm(list=ls())
setwd("/Lab6/longitdunalRNASeqData")
myfile= read.table("nc101_scaff_dataCounts.txt", sep="\t", header=TRUE, row.names =1)
#myT
myfile <- myfile[ apply( myfile,1, median)> 5,]
#myT
myfileNorm <- myfile
#myTNorm
for ( i in 1:ncol(myfile))
{
  colSum = sum(myfile[,i])
  myfileNorm[,i] =myfileNorm[,i]/colSum
}
myfileNorm
# Question 1A
Pvalues_categories <- vector()
for ( i in 1:nrow(myfileNorm)){
  day2 <- as.numeric(myfileNorm[i,1:3])
  days_86 <- as.numeric(myfileNorm[i,4:6])
  days_128 <- as.numeric(myfileNorm[i,7:11])
  mydata <- c(day2,days_86,days_128)
  
  genes <- c(rep("day2",length(day2)),rep("days_86",length(days_86)),rep("days_128",length(days_128)))
  
  #genes <- factor(genes)
  mylm <- lm(mydata ~ genes, x = TRUE)
  Pvalues_categories[i] <- anova(mylm)$"Pr(>F)"[1]
}
#Pvalues
hist(Pvalues_categories,main="Anova P values for DAY 2, DAY 86 and DAY 128", xlab="P-values", ylab=" Frequency",col="magenta")
adjPvalsCategorial <- sum(p.adjust(Pvalues_categories, method = "BH")<= 0.05)
adjPvalsCategorial

# Question 1B
Pvalues_time <- vector()
for(i in 1:nrow(myfileNorm)){
  day2 <- as.numeric(myfileNorm[i,1:3])
  day_86 <- as.numeric(myfileNorm[i,4:6])
  day_128 <- as.numeric(myfileNorm[i,7:11])
  mydata <- c(day2,day_86,day_128)
  Time <- c(rep(2,length(day2)),rep(86,length(day_86)),rep(128,length(day_128)))
  mylm_time <- lm(mydata ~ Time, x= TRUE)
  Pvalues_time[i] = anova(mylm_time)$"Pr(>F)"[1]
}
#Pvalues_time
hist(Pvalues_time,main="Regression P values as a function of Time", breaks=50, xlab="P-values",ylab="Frequency", col="cyan")
adjPvalsTime <- sum(p.adjust(Pvalues_time, method = "BH")<= 0.05)
adjPvalsTime

#Question 1C

onewayPvalue <- vector()
linearPvalues <- vector()
diff_model_pvalues <- vector()
index <- vector()
for (i in 1:nrow(myfileNorm)) {
  index[i] <- i
  day3 <- as.numeric(myfileNorm[i,1:3])
  week12 <- as.numeric(myfileNorm[i,4:6])
  week18 <- as.numeric(myfileNorm[i,7:11])
  weeksOld<- c(day3,week12,week18)
  genes <- factor(c(rep("day3",length(day3)),rep("week12",length(week12)),rep("week18",length(week18))))
  Time <- c(rep(2,length(day3)),rep(86,length(week12)),rep(128,length(week18)))
  OnewayAnova <- lm(weeksOld ~ genes)
  Linearreg<- lm(weeksOld ~ Time)
  onewayPvalue[i] <- anova(OnewayAnova)$"Pr(>F)"[1]
  linearPvalues[i] <- anova(Linearreg)$"Pr(>F)"[1]
  OnewayResiduals <- sum(residuals(OnewayAnova)^2)
  LinearRegResiduals <- sum(residuals(Linearreg)^2)
  diff <- ((LinearRegResiduals - OnewayResiduals)/2)/(OnewayResiduals/8)
  diff_model_pvalues[i] <- pf(diff,2,8,lower.tail = FALSE)
}
hist(diff_model_pvalues, main = "Histogram of Difference model",breaks = 40, xlab="P-values",ylab="Frequency",col="orange")
adjPvalsDiff<- sum(p.adjust(diff_model_pvalues, method = "BH")<= 0.05)
adjPvalsDiff

# Question 1D
categories <- factor( c( rep("day3",3),rep("week12",3),rep("week20",5)  ))
weeksOld <- c(day3,week12,week18)
Time <- c(rep(2,length(day3)),rep(86,length(week12)),rep(128,length(week18)))
myFrame <- data.frame( index, onewayPvalue,linearPvalues,diff_model_pvalues)
Anova_model <- myFrame[ order(myFrame$onewayPvalue), ] 
boxplot( as.numeric( myfileNorm[ Anova_model$index[1],]) ~ categories, main="Anova Model", xlab="Category",ylab="Relative Abundance" ,col="blue") 

Regression_model <- myFrame[ order(myFrame$linearPvalues), ] 
boxplot( as.numeric( myfileNorm[ Regression_model$index[1],]) ~ Time,main="Regression Model",xlab="Days",ylab="Relative Abundance",col="green" ) 
abline(lm(weeksOld ~ Time), col="red")

Difference_Model<- myFrame[ order(myFrame$diff_model_pvalues), ] 
boxplot( as.numeric( myfileNorm[ Difference_Model$index[1],]) ~ categories , main="Difference Model", xlab="Category",ylab="Relative Abundance",col="yellow") 

# Question 1E
# I think, the three parameter model is more appropriate for this 
# dataset than the 2 parameter model as the P value for two parameter model is low. 
# The 3 parameter model is much better eventhough losing 2 degrees of freedom. 
# Also, the slope does not fit the data properly in the two 
# parameter model.

