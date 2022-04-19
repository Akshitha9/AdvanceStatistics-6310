# 1A)  Plot the prior graph for a situation for a coin where the prior belief
# for p(head) is represented by the following R code :
# dexp(x, rate =5) / 0.9932621 
# for values of 0 <= x <= 1 and 0 otherwise.  
# (We choose the denominator to make the Integral between 0 and 1 sum to 1).
x= seq(0,1,0.001)
plot(dexp(x, rate =5) / 0.9932621,main="Prior Graph", xlab="probability of heads", ylab="prior for p(heads)")
# (1B)  Calculate the posterior graph with both the Metropolis algorithm and
#grid approximation for a case with 14 heads and 10 tails (where x = prob(head)) .
# Show the two methods roughly agree.  Compare these to a plot with a posterior for 
# new data of 14 heads and 10 tails with a prior with beta(40,40).
# (So for the observation of 14 heads and 10 tails you will end up with a graph with
#three plots superimposed: 
#(i) the Metropolis algorithm 
piOld =0.5 #let us choose the initial value of pi old as 0.5
numberOfIterations= 50000
posteriorDistr=vector()
for (i in 1:numberOfIterations)
{
  # prior heads=9 and prior tails =9 (let's choose)
  # new data heads= 14 and tails=10
  
  pOld = (dexp(piOld,rate=5)/0.9932621)* dbinom(14,24,piOld)
  
  piNew =piOld +rnorm(1,0,sd=0.01);
  if (piNew >1)
    piNew=1;
  if (piNew <0)
    piNew=0;
  pNew = (dexp(piNew,rate=5)/0.9932621)*dbinom(14,24,piNew)
  
  ratio = pNew/pOld

  if( ratio >1 || ratio>=runif(1))
    piOld =piNew;
  posteriorDistr[i] =piOld;
  
  if (i %% 100 == 0)
  {
    myHist = hist(posteriorDistr,breaks =200, plot=FALSE)
    plot(myHist$mids, myHist$counts/i, main= paste ("iterations",i))
    dbetasum =sum(dbeta(myHist$mids, 14+40,10+40))
    lines(myHist$mids,dbeta(myHist$mids,14+40,10+40)/dbetasum, col="red")
    #Sys.sleep(.1)
    }
}
#(ii) grid approximation  
numBreaks = 1000;
posteriorDistr = vector (length =numBreaks)
xVals <- seq(0,1,1/numBreaks);

i=1;
sum = 0;
for( x in xVals)
{
  # our prior with 9 heads and 9 tails
  # our new data with 14 heads and 10 tails
  posteriorDistr[i] =(dexp(x,rate=5)/0.9932621)* dbinom(14,24,x)
  sum=sum +posteriorDistr[i];
  i = i+1;
}
plot(posteriorDistr/sum, main="Grid approximation")
lines(dbeta(xVals, 40+14,40+10)/sum(dbeta(xVals, 40+14,40+10)),col ="green")

#(iii) exact analytical solution from a beta(40,40) prior make the plots different
#colors so you can visualize them…)
x = seq(0,1,0.001)
plot(dbeta(x,14+40,10+40),col="blue")

# the above three methods gave similar posterior plots.

#(1C)Repeat the above calculation but for a case of 583 heads and 417 tails. 
#(You may need to adjust your model step parameters to try and get the grid 
#and Metropolis graphs to match up).  How do the three posterior curves relate
#to each other now?  Why does this plot look different than the plot in (1B)?
# (i) Metropolis Algorithm
piOld =0.5 #let us choose the initial value of pi old as 0.5
numberOfIterations= 50000
posteriorDistr=vector()
for (i in 1:numberOfIterations)
{
  # prior heads=9 and prior tails =9 (let's choose)
  # new data heads= 14 and tails=10
  
  pOld = (dexp(piOld,rate=5)/0.9932621)* dbinom(583,1000,piOld)
  
  piNew =piOld +rnorm(1,0,sd=0.01);
  if (piNew >1)
    piNew=1;
  if (piNew <0)
    piNew=0;
  pNew = (dexp(piNew,rate=5)/0.9932621)*dbinom(583,1000,piNew)
  
  ratio = pNew/pOld
  
  if( ratio >1 || ratio>=runif(1))
    piOld =piNew;
  posteriorDistr[i] =piOld;
  
  if (i %% 100 == 0)
  {
    myHist = hist(posteriorDistr,breaks =200, plot=FALSE)
    plot(myHist$mids, myHist$counts/i, main= paste ("iterations",i))
    dbetasum =sum(dbeta(myHist$mids, 583+40,417+40))
    lines(myHist$mids,dbeta(myHist$mids,583+40,417+40)/dbetasum, col="red")
    #Sys.sleep(.1)
  }
}
#(ii) grid approximation  
numBreaks = 1000;
posteriorDistr = vector (length =numBreaks)
xVals <- seq(0,1,1/numBreaks);
#heads=583; tails=417;
i=1;
sum = 0;
for( x in xVals)
{
  # our prior with 9 heads and 9 tails
  # our new data with 14 heads and 10 tails
  posteriorDistr[i] =dexp(x,rate=5)/0.9932621* dbinom(583,1000,x)
  sum=sum +posteriorDistr[i];
  i = i+1;
}
plot(posteriorDistr/sum, main="Grid approximation for 583 heads and 417 tails")
lines(dbeta(xVals, 40+583, 40+417)/sum(dbeta(xVals, 40+583, 40+417)),col ="green")

#(iii) exact analytical solution from a beta(40,40) prior make the plots different
#colors so you can visualize them…)
x = seq(0,1,0.001)
plot(dbeta(x,583+40,417+40),col="blue", main="Exact analytical for 583 heads and 417 tails") 

# The three posterior curves show similar probability curves. By comparing the plots between 1b and 1c, as the no of flips increased, the relative probability converges the target distribution.

