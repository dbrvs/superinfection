## author = amalia magaret
## date started = 2016 Apr 29
## PI = Daniel Reeves
## uses EM algorithm to infer the best fit ZTP parameter (thereby determining the average richness or superinfection probability) from participants in a trial with a virus like HSV-2 that only allows a single strain to be detected in a single sample. The strain abundance is modeled with an exponential function, which could be changed.

## uses helper functions
source('mk_matind.R') #creates all possible combinations of true strains resulting in the number of strains observed
source('mk_matposs.R') #gives the number of ways to split up n samples to get the the sizes of the numbers of strains observed

## main code for HSV superinfection
HSV_superinfection <- function(n,vecn,alpha) {
  
  ## inputs
  # n = number of samples per person, must be >2
  # vecn = vector of ppts with 1, 2, ... n observed strains
  # alpha = superinfection parameter (if alpha=0 strains are evenly abundant)
  
  #check length of vecn
  if (length(vecn) > n) {print("error: more observed infections than samples")
    break}
  
  N <- sum(vecn) #total number of trial participants
  
  lammin <- .01  ;  lammax <- 4 #boundaries for ZTP parameter

  lamtry <- seq(lammin,lammax,by=0.001) #all lambdas to try
  
  expectR <- (N+vecn[2])/N  #
  
  explam <- function(lam) {lam/(1-exp(-lam))} #avg Richness from ZTP
    
  ## see which values of lam give us expectK close to observed
  lamtest <- explam(lam=lamtry)
  lamcurr <- lamtry[abs(lamtest-expectR)==min(abs(lamtest-expectR))]
  
  ## need two inital values of lambda since loop depends on the distance between successive estimates from M step
  lam <- lamcurr + .5
  ## max r should be infty but just need large enough so that its probability gets small, infections wont be huge
  maxr <- 12
  seqr <- 1:maxr ; matr <- seqr%*%t(rep(1,n))
  
  #################### start loop ###########################
  ntries <- 0
  ## do loop no more than 50 times and stop if estimates close
  while (abs(lam-lamcurr) > .001 & ntries < 50) {
    ntries <- ntries + 1
    #################### E step ##############################
    ## first time use initialized lambda, o.w. from end of M step
    lam <- lamcurr  
    ## this is the new numerator using the exponential decay
    ##   but if alpha = 0 simplifies to evenness assumption
    
    ## start with the conditional probabilities of r0 given r and n true
    ## make all possibilies of true and observed numbers of strains
    matprobs <- matrix(0,nrow=maxr,ncol=n)
    dimnames(matprobs) <- list(paste('r',1:maxr,sep=''),paste('r0',1:n,sep=''))
    
    rr <- 1 
    for (indr in 1:maxr) {  ## loop over number of possible true strains
      ## peq is the probability of each strain from a single sample
      seqrloop <- 1:indr 
      peq <- exp(-seqrloop*alpha) / sum(exp(-seqrloop*alpha))
      ## r0 must be no greater than n samples, no greater than r true strains 
      for (indr0 in 1:min(indr,n)) {  ## loop over number of observed strains
        ## for each iteration of this loop, 
        ##   the conditional probabilities p(r0|r,n) should add to 1
        matind <- mkmatind(r=indr,r0=indr0) ; di <- max(1,dim(matind)[1])
        matposs <- mkmatposs(r0=indr0,n=n) ; dp <- max(1,dim(matposs)[1])
        ## now need to make one row for every combination or matind and matposs
        ## combination will be all possible ways to observe indr0 strains
        ## bigind gives which strains observed
        bigind <- matrix(as.vector(t(matind)),nrow=di*dp,ncol=indr0,byrow=T)
        ## bigposs says how many of each of those named in bigind
        if (dp == 1 ) {   ## placeholder only need if r0=1
          bigposs <- matrix(matposs,nrow=di,ncol=indr0,byrow=T) 
        }
        if (dp > 1 ) {
          bigposs <- matrix(matposs[1,],nrow=di,ncol=indr0,byrow=T)
          for (pp in 2:dp) {  addon <- matrix(matposs[pp,],nrow=di,ncol=indr0,byrow=T)
             bigposs <- rbind(bigposs,addon) } }
        ## the product of each row of pmat is the probability of that combination
        ##   of strains and numbers of strains making up the observed number
        ## these first part is the probabilities in order
        matpeq <- matrix(peq[bigind],nrow=di*dp,ncol=indr0,byrow=F)
        ## next exponentiate them for the number of times drawn
        pmat <- (matpeq)^bigposs 
        ## and compute probability for this combo as the product
        jointprobs <- apply(pmat,1,prod)
        matprobs[indr,indr0] <- sum(jointprobs)
      } ## end loop over r0
    }  ## end loop over r
    
    ## now finishing up E step
    ## recall matprobs is f(r0|r,n)
    ## so tosum is f(r0,r|n) = f(r0|r,n) * f(r)
    tosum <- matprobs * (lam^matr / factorial(matr))    
    ## and condprobs is f(r|r0,n) = f(r0,r|n) / f(r0|n)
    ##   where f(r0|n) = sum over r of f(r0,r|n) using colSums
    condprobs <- tosum/(rep(1,maxr)%*%t(colSums(tosum)))

        
    #################### M step ##############################
    ## to solve this, make left side equal to right side by grid search
    ## the left side of the equation uses the condition probs
    ##    and the right side has to do with lambda
    addup <- colSums(matr*condprobs)
    lside <- sum(vecn*addup)
    rside <- function(lam,N) { N*lam*exp(lam)/(exp(lam)-1) }
    lamtest <- rside(lam=lamtry,N=N)
    ## once you compute the right side of the equation at all these values
    ##    then pick the one that makes left and ride sides closest
    lamcurr <- lamtry[abs(lamtest-lside)==min(abs(lamtest-lside))]
    #cat("ntries=",ntries,'\n')
    #print(round(condprobs[1:4,],2)) ; print(lamcurr)
  }  ## end while loop
  ## compute expected infection pp and probability of superinf
  expectr <- lamcurr/(1-exp(-lamcurr))
  ## p(super) = p(r>=2) = 1 - p(r=1)
  psuper <- 1 - lamcurr/(exp(lamcurr)-1)
  list(lamcurr=lamcurr,expectr=expectr,psuper=psuper,condprobs=condprobs)  
} ## end of function

