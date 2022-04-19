# Problem 1
#Question 1
rolls<- c(2,3,2,6,3,5,6,2,6,6,2,6,6,2,3,6,6,6,5,6,6,5,6,6,6,6,6,4,6,3,3,3,6,6,5,6,6)

posteriorProbabilities<- function(rolls)
{
  #posteriorList <- vector(mode ="list",length= length(rolls))
  posterior<- vector();
  titleStr<- ""
  FairDie<-c(1/6,1/6,1/6,1/6,1/6,1/6)
  LoadedDie<-  c(1/10,1/10,1/10,1/10,1/10,5/10)
  FairDieProb =0.99 #probibility that it is fair
  loadedDieProb = 0.01 # probability that it is unfair.
  
  for (i in 1:length(rolls)){
    
    posterior[i] <- loadedDieProb
    denominator <- (FairDie[rolls[i]]*FairDieProb) +(LoadedDie[rolls[i]]*loadedDieProb)
    #print(denom)
    loadedDieProb = (LoadedDie[rolls[i]]*loadedDieProb)/denominator
    FairDieProb = (FairDie[rolls[i]]*FairDieProb)/denominator
    titleStr<- paste(titleStr,rolls[i],sep ="")
    plot(1:i,posterior, main = paste("Posterior probability of the Loaded Die\n",titleStr),xlim=c(1,length(rolls)+1),ylim=c(0,1),xlab="Rolls",ylab= "Loaded die probability")
    Sys.sleep(1)
  }
}
posteriorProbabilities(rolls)

#Question 2

# It takes approximately 23-25 rolls inorder to be 0.99% confident that it is a loaded die.

# Problem 2
#Question 2A
#pwd= person with disease
#pwod= person without disease
#ppwd= positive that a person has a disease
#ppwod = positive that a person does not have a disease
#npwd= negative that a person have a disease
#npwod= negative that a person does not have a disease
# Question 2 A
pwd= 0.001
pwod=0.9999
liklihood_has_disease = c(0.91,0.09) # ppwd=0.91; npwd= 0.09
liklihood_no_disease = c(0.16,0.84) # ppwod =0.16; npwod = 0.84
#cost_test =1$
requiredPosterior=0.99999
numOftests=1:30
numOfSim=1000000
AvgPosteriorVals <- c(length=length(numOftests))
estimation1 <- c(length=length(numOftests))
likelihoodData<-function(TypeOfliklihood,times)
{
  v <- vector(mode = "integer", length = times)
  for (i in 1:times) 
  { 
    if(runif(1) <= TypeOfliklihood[1])
    {
      v[i] <- 1
    }
    else 
    {
      v[i] <- 2
    }
  }
  return(v)
}

#likelihoodData(liklihood_has_disease,30)
for (i in numOftests) 
{
  posteriorvalues <- c(length=numOfSim)
  for (j in 1:numOfSim) 
  { 
    prior <- c(0.001,0.9999)
    
    df <- likelihoodData(liklihood_has_disease,i)
    
    for (k in 1:length(df)) 
    { 
      denominator <- prior[1]*liklihood_has_disease[df[i]]+prior[2]*liklihood_no_disease[df[k]]
      
      prior[1] = prior[1] * liklihood_has_disease[df[k]] / denominator
      prior[2] = prior[2] * liklihood_no_disease[df[k]] / denominator
    }
    posteriorvalues[j] = prior[1]
  }
  AvgPosteriorVals[i] <- mean(posteriorvalues)
  #print(AvgPosteriorVals)
  estimation1[i] <- sum(posteriorvalues >= requiredPosterior)/numOfSim 
}
estimation1
plot(numOftests,estimation1, main = "Tests performed for a diseased person", xlab="Number of tests", ylab= "Estimation for accuracy")
# from the distribution, we can say that the hospital requires atleast 19 tests to say that a person is diseased.
#Question 2B
pwd= 0.001
pwod=0.9999
liklihood_has_disease = c(0.91,0.09) # ppwd=0.91; npwd= 0.09
liklihood_no_disease = c(0.16,0.84) # ppwod =0.16; npwod = 0.84
#cost_test =1$
requiredPosterior=0.99999
numOftests=1:30
numOfSim=1000000
AvgPosteriorVals <- c(length=length(numOftests))
estimation2 <- c(length=length(numOftests))
likelihoodData<-function(TypeOfliklihood,times)
{
  v <- vector(mode = "integer", length = times)
  for (i in 1:times) 
  { 
    if(runif(1) <= TypeOfliklihood[1])
    {
      v[i] <- 1
    }
    else 
    {
      v[i] <- 2
    }
  }
  return(v)
}

for (i in numOftests) 
{
  posteriorvalues <- c(length=numOfSim)
  for (j in 1:numOfSim) 
  { 
    prior <- c(0.001,0.9999)
    
    df <- likelihoodData(liklihood_no_disease,i)
    
    for (k in 1:length(df)) 
    { 
      denominator <- prior[2]*liklihood_no_disease[df[i]]+prior[1]*liklihood_has_disease[df[k]]
      
      prior[1] = prior[1] * liklihood_has_disease[df[k]] / denominator
      prior[2] = prior[2] * liklihood_no_disease[df[k]] / denominator
    }
    posteriorvalues[j] = prior[2]
  }
  AvgPosteriorVals[i] <- mean(posteriorvalues)
  #print(AvgPosteriorVals)
  estimation2[i] <- sum(posteriorvalues >= requiredPosterior)/numOfSim 
}
estimation2
plot(numOftests,estimation2,main = "Tests performed for a undiseased person", xlab="Number of tests", ylab= "Estimation for accuracy") 
# from the distribution, we can say that the hospital requires atleast 10 tests to say that a person is not diseased.
#Question 2C
#For the Hospital to perform tests with the required posterior of 0.9999 for a million people with the cost of 1$ for each tests, it would cost around 10 - 12 million dollars per year based on what the patient is tested.  
