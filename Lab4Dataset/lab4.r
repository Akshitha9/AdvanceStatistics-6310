#Question1
setwd("/Users/akshitha/Desktop/Final_sem/Advancedstatistics/Lab4Dataset")
mytable<- read.table("nc101_scaff_dataCounts.txt",header=TRUE,row.names=1)
head(mytable)

#Question2
#install.packages("ggplot2)
library(ggplot2)	
library(scales)
countsplot <- ggplot(mytable, aes(x = D2_01,y = D2_02)) +geom_point() + scale_x_log10(labels = trans_format("log10", math_format(10^.x)),oob = scales::squish_infinite) + scale_y_log10(labels = trans_format("log10", math_format(10^.y)),oob = scales::squish_infinite)
countsplot
#They appear to have similar patterns of expression.
#Question 3

plot( log10(apply(mytable,1,mean)), log10(apply(mytable,1,var)),main="Mean vs Variance")
lines(log10(apply(mytable,1,mean)), log10(apply(mytable,1,mean)), col="red")

#Question 4
D2_01_NC101_00003<-mytable$D2_01[1]
#D2_01_NC101_00003
D2_01_notNC101_00003<-sum(mytable$D2_01)-D2_01_NC101_00003
#D2_01_notNC101_00003
D2_02_NC101_00003<-mytable$D2_02[1]
#D2_02_NC101_00003
D2_02_notNC101_00003<-sum(mytable$D2_02)-D2_02_NC101_00003
#D2_02_notNC101_00003
#Convert to matrix
contingencyTable1<-matrix(c(D2_01_NC101_00003,D2_01_notNC101_00003,D2_02_NC101_00003,D2_02_notNC101_00003),nrow = 2)
colnames(contingencyTable1) <- c("Sequences in D2_01","Sequences in D2_01")
rownames(contingencyTable1)<- c("Assigned to NC101_00003","Not assigned to NC101_00003")
contingencyTable1
FisherTest1<-(fisher.test(contingencyTable1))
FisherTest1
Pvalue <- FisherTest1$p.value
Pvalue

# Question 5

#For all the genes
pValues1<-vector()
nrow(mytable)
for (i in 1:nrow(mytable)){
  
D2_01withgene<-mytable$D2_01[i]

D2_01_withoutgene<-sum(mytable$D2_01)-D2_01withgene

D2_02withgene<-mytable$D2_02[i]

D2_02_withoutgene<-sum(mytable$D2_02)-D2_02withgene

contingencyTable2<-matrix(c(D2_01withgene,D2_01_withoutgene,D2_02withgene,D2_02_withoutgene),nrow = 2)
#colnames(contingencyTable) <- c("Sequences in D2_01","Sequences in D2_01")
#rownames(contingencyTable)<- c("Assigned to gene","Not assigned to gene")
FisherTest2<-(fisher.test(contingencyTable2))
P<- FisherTest2$p.value
pValues1[i]=P
}
#length(pValues1)
hist(pValues1,breaks = 20, main= "Fisher test P_values across all genes", col="green")
# The plot displays non uniform distribution. We can expect that from the intuition that some genes depend on other genes for expression.
# under a uniform distribution, we would expect that the probability of the expression of all the genes are same. But in reality, the genes
# will not have same frequency of expression. In this case, there are some genes which are highly expressed,
# some genes are less expressed while others are moderately expressed.
#If we remove the low abundance genes from the given table:
mytable1 = mytable[(mytable$D2_01 +mytable$D2_02>50),]
#head(mynewtable)
newpValues = vector()
for (i in 1:nrow(mytable1))
{
  D2_01withgene<-mytable1$D2_01[i]
  
  D2_01_withoutgene<-sum(mytable1$D2_01)-D2_01withgene
  
  D2_02withgene<-mytable1$D2_02[i]
  
  D2_02_withoutgene<-sum(mytable1$D2_02)-D2_02withgene
  
  contingencyTable3<-matrix(c(D2_01withgene,D2_01_withoutgene,D2_02withgene,D2_02_withoutgene),nrow = 2)
  #colnames(contingencyTable) <- c("Sequences in D2_01","Sequences in D2_01")
  #rownames(contingencyTable)<- c("Assigned to gene","Not assigned to gene")
  FisherTest3<-(fisher.test(contingencyTable3))
  P<- FisherTest3$p.value
  newpValues[i]=P
}

#tail(newpValues)
hist(newpValues, breaks =20)
# there is a reduction in the frequency of the p-values.

# Question 6
#Add 1 to every value in the table ( with something like myT = myT + 1 ).  This is called adding a pseudo-count.  Now consider the first gene (NC101_00003 ) again.  From the first experiment, calculate 

#expected frequency = p = 
  # Assigned to NC101_00003 in D2_01)/total # of sequences in D2_01)
    
    #Now use poisson.test to assign a p-value for the null hypothesis that value of p derived from D2_01 could have produced the number of reads observed for this gene in D2_02 .
  
mytable2= mytable+1
mytable2
Assigned_NC101_00003_D2_01=mytable2$D2_01[1] # no. assigned to the gene NC101_00003_D2_01 im D2_01
Assigned_NC101_00003_D2_01

TotalSequences = sum(mytable2$D2_01)
TotalSequences

ExpectedFreq = Assigned_NC101_00003_D2_01/TotalSequences
ExpectedFreq
#For D2_02
Assigned_NC101_00003_D2_02 =mytable2$D2_02[1]
Assigned_NC101_00003_D2_02

TotalnumberOfSeq_D2_02 = sum(mytable2$D2_02)
TotalnumberOfSeq_D2_02
PT<-poisson.test(Assigned_NC101_00003_D2_02 ,TotalnumberOfSeq_D2_02,ExpectedFreq)
PT$p.value

# Question 7
poissonPvalues = vector()
mytable3= mytable+1
mytable3
for (i in 1:nrow(mytable3)) 
{
Assignedgenes=mytable3$D2_01[i] # no. assigned to the gene NC101_00003_D2_01 im D2_01
Assignedgenes

TotalSequences = sum(mytable3$D2_01)
TotalSequences

ExpectedFreq = Assignedgenes/TotalSequences
ExpectedFreq
PT<-poisson.test(Assignedgenes ,TotalnumberOfSeq_D2_02,ExpectedFreq)
poissonPvalues[i]<-PT$p.value
}
head(poissonPvalues)
#fisherVspoisson<- ggplot(mytable, aes(x = pValues1,y = poissonPvalues)) +geom_point() + scale_x_log10(labels = trans_format("log10", math_format(10^.x)),oob = scales::squish_infinite) + scale_y_log10(labels = trans_format("log10", math_format(10^.y)),oob = scales::squish_infinite)
#fisherVspoisson
#length(poissonPvalues)
plot(log10(pValues1),log10(poissonPvalues),main="Fisher test Vs Poisson test", xlab="Fisher p-values",ylab="Poisson p-values")

