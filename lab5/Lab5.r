
# question 1
#install.packages("ggplot2) in the console
library(ggplot2)
library(scales)
myfile<- read.table("/Users/akshitha/Desktop/Final_sem/Advancedstatistics/lab5/cancerRisk.txt", header = TRUE,sep = "\t")
myfile
data <- ggplot(myfile, aes(x = CumulativeCellDivisions,y = Lifetime_cancer_risk, label = myfile$Cancer_type)) +geom_point() + scale_x_log10(labels = trans_format("log10", math_format(10^.x))) + scale_y_log10(labels = trans_format("log10", math_format(10^.y)))
data
data + geom_text(check_overlap = TRUE, size = 1.5) 
#Question 1B
myLm<-lm(Lifetime_cancer_risk~CumulativeCellDivisions,data=myfile)
myLm
par(mfrow=c(2,2))
plot(myLm, which=2,col="blue")
#devAskNewPage(ask = FALSE)
abline(myLm,col="red")
#Question 1C
summary(myLm)
# the p-value is 0.002028 which is significant and the  R squared for the model is 0.2839.
# Question 1D
# Yes, from te above QQplot we can observe that all the data points are perfectly aligned with 
# the fitted line by which we can say that the variance and the residual errors is constant by making
# the model reasonable.

# Question 2

caseControlData<- read.table("/Users/akshitha/Desktop/Final_sem/Advancedstatistics/lab5/caseControlData.txt",header = TRUE,sep = "\t")
#caseControlData
BMI<-read.table("/Users/akshitha/Desktop/Final_sem/Advancedstatistics/lab5/BMI_Data.txt",header=TRUE, sep="\t")
BMI
datafilter<- vector() # removing case and control from the sample id's
for (i in 1:nrow(caseControlData)){
  datafilter[i]<- sub("case","",caseControlData$sample[i])
  datafilter[i]<- sub("control","",datafilter[i])

#datafilter
#removing the extraneous information from the suffix of the filtered data
 datafilter[i]<- strsplit(datafilter[i],"_")[[1]][1]
}
caseControlData$sample<- datafilter
caseControlData$sample
names(caseControlData)[names(caseControlData)=="sample"] <- "studyid"

# merging both the tables
newtable <- merge(BMI,caseControlData,by.x="studyid")
#newtable<- na.omit(newtable)
newtable
# linear model
p_values= vector()
for (x in 1: ncol(newtable-2))
{
  y<-x+2
  OTU<- newtable[,2]
  bodyMassIndex<-newtable[,y]
  mylm<-lm(bodyMassIndex~OTU)

  An<-anova(mylm)
  #print(An)
  p_values[x]<- anova(mylm)$"Pr(>F)"
  }
p_values
hist(p_values, main= "Histogram of P_values for each OTU",col="blue")
# The p-values shows a uniform distribution thus indicating the microbial community appears to be influencing the body weight in this cohort.
#For false discory rate
adjustedp_values = sum(p.adjust(p_values,method="BH")<=0.10)
adjustedp_values
# with the correction of P-values at 10% false discovery rate there were no significant associations found.

# there are no significant associations found at a 10% false discovery rate.