
## author = amalia magaret
## date started = 2017 Jan 29
## PI = Daniel Reeves
## objective = run EM algorithm for determining number of true
##   infections from participants who got n samples each
## this time specify a different alpha from simulated
## pick 3 samples per each of 500 persons

## to use
rm(list=ls())
source('superinfection.R')

## here are the parameters you need to choose
##   lambda = true Poisson parameter for number of infections
##   alpha = exponential abundance parameter
##   n = samples per person

## make range of other conditions
## start with range of superinfection rates and transform into lambdas
## from psuper <- 1 - lamcurr/(exp(lamcurr)-1)
## use interpolation, first show relation, then pull out desired lambda
lamck <- seq(.005,1,.005)
pck <- 1 - lamck/(exp(lamck)-1)
superuse <- c(.01,.02,.05,.1,.2,.3)  ## superinfection rates
lamuse <- approx(x=pck,y=lamck,xout=superuse)$y
nuse <- 3  
npers <- 80000
nsimulations <- 2000 ## number of times to simulate and estimate lambda

## run evenness conditions first, where alpha = 0
alpha <- 0
alphatry <- alpha + c(-1,-.5,-.25,-0.1,0,0.1,0.25,.5,1)
alphause <- alphatry[alphatry >= 0]
## loop over simulations to estimate lambda over lamuse and alphause
for (ll in 1:length(lamuse)) {
  for (aa in 1:length(alphause)) {
    filename <- paste("results/lamest_sup",superuse[ll],"_alpha",alpha,"_ause",alphause[aa],".txt",sep="")
    for (nsim in 1:nsimulations) {
    ## make simulated data, use first 500 persons infected
    rinf <- rpois(n=npers,lambda=lamuse[ll])
    ## you will get 500 b/c start with 80000 and lowest superuse is 1% so 800
    rpos <- rinf[rinf>0][1:500]
    ## make a table of true number of infections
    ## need to fill in zeros when rpos not contiguous
    ##  so add 1 case for every value bwtn 1 and max in the vector, 
    ##    and then subtract a count of 1 for each item in the table
    maxr <- max(rpos)  ## note that even for highest lam
    Nr <- table(c(rpos,1:maxr))-rep(1,maxr)
    truenum <- 1:maxr
    
    ## vectrue is an ordered vector of numbers of strains, 
    ##    one for each person, with true number of strains for each
    ##  or could just sort rinf ##@ all(vectrue==sort(rpos))
    vectrue <- rep(truenum,Nr)
    
    ## now implement measurement error via sampling
    ## and for each person, simulate which strain detection for each sampling
    vecobs <- rep(1,Nr[1])  ## will build on this vector
    if (maxr > 1) {  ## but don't need to if only one true infection will observe ones
    for (j in 2:maxr) { ## first loop over each true number of infections
      ## can skip j=1 because can only sample one if one true infection 
      ## also skip if no simulated persons with that number of strains
      if (Nr[j] > 0) {
        psamp <- exp(-alpha*(1:j))/sum(exp(-alpha*(1:j)))
        for (inds in 1:nuse) { ## then loop over samples for each person
          detect <- rmultinom(Nr[j],size=1,psamp)
          ## dim(detect) ## will be length(psamp) rows by Nr[j] columns
          ## rbind(round(rowMeans(detect),4),round(psamp,4))
          ## add them up, to see which ones got detected
          if (inds==1) {keepdet <- detect}  else {keepdet <- keepdet + detect}
        }  ## at end, if above zero, detected that strain at least once
        everdet <- 1*(keepdet > 0)  
        howmanydet <- colSums(everdet)  ## how many unique strains per person
        vecobs <- c(vecobs,howmanydet)  ## append observed numbers of strains
        ## can make a matrix of true by simulated numbers of strains
        ##@ tabb <- table(vectrue,vecobs)
        ##@ probsest <- tabb / rowSums(tabb)
        ## can show that probsest matches matprobs from EM_superinfection_anym.R
      }  ## end over whether Nr[j] > 0
    } ## end over j from 2 to maxr
    } ## end over whether observed any superinf (maxr > 1)  
    ## run to maximize lambda, see if get right one
    ## make table of vecobs that includes all possible values to nuse
    Nr0 <- table(c(vecobs,1:nuse))-rep(1,nuse)
    anss <- EM_superinfection_anym(n=nuse,vecn=Nr0,alpha=alphause[aa])
    lamest <- anss$lamcurr
    converged <- anss$conv
    print(lamest)
    write.table(x=lamest,file=filename,append=T, row.names=F,col.names=F)
    } ## end loop over nsimulations
    ##@    lammat[ll,nn] <- lamest (old version with one iteration)
  } ## end loop over number of samples to take
}  ## end loop over true lambda


