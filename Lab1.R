#Question 1
loadedDie<-c(1,2,3,4,5,6)
probilityOfDie<- c(rep(0.1,5),0.5)
chanceOfDie <- loadedDie* probilityOfDie
chanceOfDie
meanValue<- 0
for (i in 1:length(chanceOfDie))
{
  meanValue <- meanValue + chanceOfDie[i]
}
meanValue
# For Calculating Variance
VarianceValue<-0

for (i in 1:length(loadedDie))
{
  
  #VarianceVector[i]<- ((loadedDie[i]-meanValue)^2)*probilityOfDie[i]
  VarianceValue <- VarianceValue+((loadedDie[i]-meanValue)^2)*probilityOfDie[i]
}
VarianceValue  

#Question 2
rollLoadedDie <-function(Rolls)
{myRolls= sample(loadedDie, size =Rolls, replace =TRUE,prob =c(0.1,0.1,0.1,0.1,0.1,0.5))
return (myRolls)}
output =rollLoadedDie(100000)
output

#Question 3
hist(output,main="Large number of rolls for a Loaded Die",xlab= "Faces of the Die", ylab="Frequency of the Rolls",col="Red")
"No, the distribution is not uniform as the outcome of a loaded die is not equally likely."

# Question 4
trailSizes<- c(5,10,15,20,25,30,40,50,100,200,300,400,500,1000,2000,3000,4000,5000,10000,20000,30000,100000)
means<-vector(mode="double",length = length(trailSizes))
variances <- vector(mode="double", length =length(trailSizes))
for (i in 1:length(trailSizes))
{
  rolls <- vector(length=trailSizes[i], mode ="double")
  for (j in 1:trailSizes[i])
  { rolls[j]<-sample(1:6, size =trailSizes, replace =TRUE,prob =c(0.1,0.1,0.1,0.1,0.1,0.5)) }
  means[i]<-mean(rolls);
  variances[i]<-var(rolls)
}
plot(log10(trailSizes),means)
lines(log10(trailSizes),rep(meanValue,length(trailSizes)))

plot(log10(trailSizes),variances)
lines(log10(trailSizes),rep(VarianceValue,length(trailSizes)))
"The convergence of the expected values for the mean and variance appears to be between 5000 and 10,000, however it took 10,000 rolls to get the convergence to the expected values."
