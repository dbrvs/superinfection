
## author = amalia magaret
## date started = 2016 Sep 08
## PI = Daniel Reeves
## objective = run EM algorithm for determining number of true
##   infections from participants who got n samples each

## to use
rm(list=ls())
source('superinfection.R')

## here are the parameters you need to choose
lambda <- 1.4  ## true Poisson parameter for number of infections
alpha <- 2.2  ## exponential abundance parameter
n <- 5  ## let there be n samples per person

## make simulated data
rinf <- rpois(n=80000,lambda=lambda)
## make a table of true number of infections
## need to fill in zeros when rinf not contiguous
##  so add 1 case for every value bwtn 1 and max in the vector, 
##    and then subtract a count of 1 for each item in the table
maxr <- max(rinf)
Nr <- table(c(rinf[rinf>0],1:maxr))-rep(1,maxr)
truenum <- 1:maxr

## vectrue is an ordered vector of numbers of strains, 
##    one for each person, with true number of strains for each
##  or could just sort rinf ##@ all(vectrue==sort(rinf[rinf>0]))
vectrue <- rep(truenum,Nr)

## now implement measurement error via sampling
## and for each person, simulate which strain detection for each sampling
vecobs <- rep(1,Nr[1])  ## will build on this vector
for (j in 2:maxr) { ## first loop over each true number of infections
  ## can skip j=1 because can only sample one if one true infection 
  ## also skip if no simulated persons with that number of strains
  if (Nr[j] > 0) {
    psamp <- exp(-alpha*(1:j))/sum(exp(-alpha*(1:j)))
    for (inds in 1:n) { ## then loop over samples for each person
      detect <- rmultinom(Nr[j],size=1,psamp)
      ## dim(detect) ## will be length(psamp) rows by Nr[j] columns
      ## rbind(round(rowMeans(detect),4),round(psamp,4))
      ## add them up, to see which ones got detected
      if (inds==1) {keepdet <- detect}  else {keepdet <- keepdet + detect}
    }  ## at end, if above zero, detected that strain at least once
    everdet <- 1*(keepdet > 0)  
    howmanydet <- colSums(everdet)  ## how many unique strains per person
    vecobs <- c(vecobs,howmanydet)  ## append observed numbers of strains
  }
}

## make a matrix of true by simulated numbers of strains
tabb <- table(vectrue,vecobs)
probsest <- tabb / rowSums(tabb)

## run to maximize lambda, see if get right one
Nr0 <- table(c(vecobs,1:n))-rep(1,n)
tryit <- EM_superinfection_anym(n=n,vecn=Nr0,alpha=alpha)
