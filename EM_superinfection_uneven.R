
## author = amalia magaret
## date started = 2016 Apr 29
## PI = Daniel Reeves
## objective = build EM algorithm for determining number of true
##   infections from participants who got m samples each
## this version adds allowance that the strains 
##   are not evenly present, that their prevalence or at
##   least the probability of detection decreases exponentially
##   with setting alpha = 2

## to use
## source('z://Schiffer/Reeves/EM_Superinfection_uneven.R')

## inputs
## n = number of participants
## m = number of samples each: require 2 for now
## n2 = number of participants with 2 sequences identified over the two samples

n <- 400
m <- 2
n2 <- 32
alpha <- 2
  
EMsuperinfun <- function(n,m,n2,alpha) {
  ################## initialize ############################
  ## start with grid search for lambda that satisfies E(K)
  ## set reasonable range for potential values of lambda
  ## these endpoints give range between 1 and 2 infections pp
  lammin <- .01  ;  lammax <- 1.6
  lamtry <- seq(lammin,lammax,by=0.001)
  expectK <- (n+n2)/n  
  explam <- function(lam) { lam/(1-exp(-lam)) }
  ## see which values of lam give us expectK close to observed
  lamtest <- explam(lam=lamtry)
  lamcurr <- lamtry[abs(lamtest-expectK)==min(abs(lamtest-expectK))]
  
  ## need two inital values of lambda since loop depends on
  ##   the distance between successive estimates from M step
  lam <- lamcurr + .5
  ## max k should be infty but just need large enough
  ##   so that its probability gets small, # infections wont be huge
  maxk <- 150
  
  ## build vectors and matrices for help in computing probs
  ## value of k corresponds to row number, k0 to column number
  seqm <- 1:m ; matm <- rep(1,maxk)%*%t(seqm)
  seqk <- 1:maxk ; matk <- seqk%*%t(rep(1,m))
    
  #################### start loop ###########################
  ntries <- 0
  ## do loop no more than 50 times and stop if estimates close
  while (abs(lam-lamcurr) > .001 & ntries < 50) {
    ntries <- ntries + 1
    #################### E step ##############################
    ## first time use initialized lambda, o.w. from end of M step
    lam <- lamcurr  
    ## this is the old numerator of the conditional probability
    ##@ tobuild <- (1/seqk * lam^seqk / factorial(seqk))%*%t(rep(1,m)) *
    ##@   (matk-1)^(matm-1)
    ## this is the new numerator using the exponential decay
    matek <- exp(-matk*alpha)
    ek.sq.then.sum <- apply(X=matek^2,MARGIN=2,FUN=cumsum)
    ek.sum.then.sq <- apply(X=matek,MARGIN=2,FUN=cumsum)^2
    tobuild <- 1/(ek.sum.then.sq) * ek.sq.then.sum^(2-matm) * 
      (ek.sum.then.sq - ek.sq.then.sum)^(matm-1) *
      (lam^matk / factorial(matk))
    
    ## zero out tobuild where k < k0 if need be?
    condprobs <- tobuild/(rep(1,maxk)%*%t(colSums(tobuild)))
    
    
    #################### M step ##############################
    ## to solve this, make left side equal to right side by grid search
    ## the left side of the equation uses the condition probs
    ##    and the right side has to do with lambda
    addup <- colSums(matk*condprobs)
    lside <- sum(c(n-n2,n2)*addup)
    rside <- function(lam,n) { n*lam*exp(lam)/(exp(lam)-1) }
    lamtest <- rside(lam=lamtry,n=n)
    ## once you compute the right side of the equation at all these values
    ##    then pick the one that makes left and ride sides closest
    lamcurr <- lamtry[abs(lamtest-lside)==min(abs(lamtest-lside))]
    cat("ntries=",ntries,'\n')
    print(round(condprobs[1:4,],2)) ; print(lamcurr)
  }  ## end while loop
  ## compute expected infection pp and probability of superinf
  expectK <- lamcurr/(1-exp(-lamcurr))
  ## p(super) = p(K>=2) = 1 - p(K=1)
  psuper <- 1 - lamcurr/(exp(lamcurr)-1)
  list(lamcurr=lamcurr,expectK=expectK,psuper=psuper)  
} ## end of function
  
ans <- EMsuperinfun(n=400,m=2,n2=32,alpha=2)   
ans  
  
  