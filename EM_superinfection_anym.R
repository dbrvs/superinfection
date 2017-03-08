
## author = amalia magaret
## date started = 2016 Apr 29
## PI = Daniel Reeves
## objective = build EM algorithm for determining number of true
##   infections from participants who got m samples each
## this version adds allowance that the strains 
##   are not evenly present, that their prevalence or at
##   least the probability of detection decreases exponentially
##   with setting alpha = 2
## it also allows n = number of samples drawn to be be treated generally 
##   rather than specifically
##   using n as an input rather than a separate formula for each 

## to use
## source('z://Schiffer/Reeves/EM_Superinfection_anym.R')
## this one gives all possible combinations of possible true strains 
##   to get the number of strains observed
source('mk_matind.R')
## this one gives the number of ways to split up n samples
##   to get the the sizes of the numbers of strains observed
source('mk_matposs.R')

## inputs
## vecn = number of participants with each number of observed specimens
##    so vecn = c(N1, N2, ..., Nn)  where can observe up to n
## n = number of samples each: require 2 for now
## N2 = number of participants with 2 sequences identified over the two samples

##@ vecn <- c(400-32,32) ; n <- 2
##@ vecn <- c(400-32,30,2) ; n <- 3
##@ alpha <- 2
  
EM_superinfection_anym <- function(n,vecn,alpha) {
  ################## inputs ############################
  ## n = number of samples per person, each with single winner, 
  ##    will be 2 or 3 or larger
  ## vecn = number of persons with 1, 2, ... n observed infections
  ## alpha = exponential parameter
  if (length(vecn) > n) {print("error: more observed infections than samples")
    break}
  N <- sum(vecn)
  ## N2 <- vecn[2] ; if (n>2) {N3 <- vecn[3]}
  ################## initialize ############################
  ## start with grid search for lambda that satisfies E(K)
  ## set reasonable range for potential values of lambda
  ## these endpoints give range between 1 and 2 infections pp
  ##@ lammin <- .01  ;  lammax <- 1.6
  ## but use wider range so can explore extreme cases
  lammin <- .01  ;  lammax <- 4
  lamtry <- seq(lammin,lammax,by=0.001)
  expectR <- (N+vecn[2])/N  
  explam <- function(lam) { lam/(1-exp(-lam)) }
  ## see which values of lam give us expectK close to observed
  lamtest <- explam(lam=lamtry)
  lamcurr <- lamtry[abs(lamtest-expectR)==min(abs(lamtest-expectR))]
  
  ## need two inital values of lambda since loop depends on
  ##   the distance between successive estimates from M step
  lam <- lamcurr + .5
  ## max r should be infty but just need large enough
  ##   so that its probability gets small, # infections wont be huge
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
    cat("ntries=",ntries,'\n')
    print(round(condprobs[1:4,],2)) ; print(lamcurr)
  }  ## end while loop
  ## compute expected infection pp and probability of superinf
  expectr <- lamcurr/(1-exp(-lamcurr))
  ## p(super) = p(r>=2) = 1 - p(r=1)
  psuper <- 1 - lamcurr/(exp(lamcurr)-1)
  list(lamcurr=lamcurr,expectr=expectr,psuper=psuper,condprobs=condprobs)  
} ## end of function
  
## do an example from HCV paper
 seqq <- seq(0,3,.25)
 pest <- seqq
 for (aa in 1:length(seqq)) 
 {
 pest[aa] <- EM_superinfection_anym(n=2,vecn=c(150,6),alpha=seqq[aa])$psuper   
 }
 
 ## do an example from Redd survey
 N2=c(1,3,1,13,1,2,2,3,2,7,12,1,3,2,10) 
 N=c(13,78,16,58,8,14,46,145,147,149,56,7,44,130,145)
 Redd_data=cbind(N-N2,N2)

 pest <- N1

 for (i in 1:length(N1)) 
 {
 pest[i] <- EM_superinfection_anym(n=2,vecn=Redd_data[i,],alpha=2)$psuper 
 }
