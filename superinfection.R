## this script builds the EM algorithm for determining number of true infections from participants who got m samples each, including uneven strain prevalence baesd on exponential ##
# author = amalia magaret
# contributor = dbr
# date started = 2016 Apr 29

## requires 2 helper files
source('mk_matind.R') # enumerate possible combinations that admit Robs
source('mk_matposs.R') # enumerate ways to split up n samples to get the the sizes of the numbers of strains observed

explam <- function(lam) { lam/(1-exp(-lam)) } #average of ero-truncated poisson distribution

#function that computes the EM algorithm
EM_superinfection_anym <- function(n,vecn,alpha) {
 
  ################## inputs ############################
  ## n = number of samples per person
  ## vecn = number of persons with 1, 2, ... n observed infections N(Robs)
  ## alpha = exponential parameter for evenness
 
  #check that vecn makes sense
  if (length(vecn) > n) {print("error: more observed infections than samples")
    break}
  
  N <- sum(vecn) #total number of samples in trial
  
  ################## initialize ############################
  lammin <- .01  ;  lammax <- 4  #start with unrealistically wide range of lambdas
  lamtry <- seq(lammin,lammax,by=0.001) #make a huge list of lambdas
  expectR <- sum((1:length(vecn))*vecn)/N #use naive MLE at first order 
  lamtest <- explam(lam=lamtry) #use TP distribution
  lamcurr <- lamtry[abs(lamtest-expectR)==min(abs(lamtest-expectR))] #best lambda given <R>
  lam <- lamcurr + 0.5   # 2nd estimate to compute initial error 
  maxr <- 12 #we are interested in relatively small richnesses in this paper
  seqr <- 1:maxr ; matr <- seqr%*%t(rep(1,n))
  
  #################### start loop ###########################
  ntries <- 0
  ## do loop no more than 50 times and stop if estimates close
  while (abs(lam-lamcurr) > .001 & ntries < 50) {
    ntries <- ntries + 1
    #################### E step ##############################
    lam <- lamcurr  #update lam with new value from M step
    
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
    cat("ntries=",ntries,'\n')
    #print(round(condprobs[1:4,],2)) ; print(lamcurr)
  }  ## end while loop
  ## compute expected infection pp and probability of superinf
  expectr <- lamcurr/(1-exp(-lamcurr))
  psuper <- 1 - lamcurr/(exp(lamcurr)-1)
  list(lamcurr=lamcurr,expectr=expectr,psuper=psuper,condprobs=condprobs)  
} ## end of function
  
