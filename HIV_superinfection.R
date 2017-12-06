## author = amalia magaret
## date started = 2017 July 13
## PI = Daniel Reeves
## objective = build EM algorithm for determining number of true
##   infections from participants who got one draw containing m reads

## to use
library(permute)

HIVsuperinfection <- function(matn,alpha) {
  ################## inputs ############################
  ## for each person, get vecn = c(n_1, n_2, ..., n_R^0)  
  ## so dataset consists of one row per person and 
  ##   N columns where N is the max over persons of n_R^0
  ## alpha = exponential parameter
  robs <- apply(1*(matn>0),1,sum)
  maxr <- max(robs)
  N <- dim(matn)[1]
  ################## initialize ############################
  ## start with grid search for lambda that satisfies E(K)
  ## set reasonable range for potential values of lambda
  ## these endpoints give range between 1 and 2 infections pp
  ##@ lammin <- .01  ;  lammax <- 1.6
  ## but use wider range so can explore extreme cases
  lammin <- .0001  ;  lammax <- 4
  lamtry <- seq(lammin,lammax,by=0.0001)
  expectR <- (sum(matn>0))/N  
  explam <- function(lam) { lam/(1-exp(-lam)) }
  ## see which values of lam give us expectK close to observed
  lamtest <- explam(lam=lamtry)
  lamcurr <- lamtry[abs(lamtest-expectR)==min(abs(lamtest-expectR))]
  laminit <- lamcurr  ## store this one, the raw pre-alg estimate
  
  ## store lambda based on raw proportion with superinfection too
  lamck <- seq(.0001,1,.005)
  pck <- 1 - lamck/(exp(lamck)-1)
  obssuper <- mean(rowSums(matn>0)>1)
  lamraw <- approx(x=pck,y=lamck,xout=obssuper)$y
    
  #print(lamraw)
  #print(laminit)
    
  ## need two inital values of lambda since loop depends on
  ##   the distance between successive estimates from M step
  lam <- lamcurr + .5
  ## max r should be infty but just need large enough
  ##   so that its probability gets small, # infections wont be huge
  maxr <- 6 ## the probability of having more than 6 true infections 
  ##   given the range of simulated lambdas is remote  --
  ##   we simulated 80000 people 10X and nobody had more than 6
  seqr <- 1:maxr ; matr <- rep(1,N)%*%t(seqr)
  dimnames(matr) <- list(paste('id',1:N,sep=''),paste('R',seqr,sep=''))
  
  #################### start loop ###########################
  ntries <- 0
  ## do loop no more than 50 times and stop if estimates close
  while (abs(lam-lamcurr) > .0005 & ntries < 150) {
    ntries <- ntries + 1
    #################### E step ##############################
    ## first time use initialized lambda, o.w. from end of M step
    lam <- lamcurr  
    ## this is the new numerator using the exponential decay
    ##   but if alpha = 0 simplifies to evenness assumption
    
    ## start with the conditional probabilities of r0 given r and n true
    ## make all possibilies of true and observed numbers of strains
    matprobs <- matrix(0,nrow=N,ncol=maxr)
    dimnames(matprobs) <- list(paste('id',1:N,sep=''),paste('R',seqr,sep=''))
    
    for (indn in 1:N) {  ## loop over persons
      vecn <- matn[indn,][matn[indn,]>0] ; m <- sum(vecn)
      r0 <- sum(vecn>0)
      ## now need to loop over possible numbers of true infections
      for (indr in r0:maxr) {  ## loop over number of possible true strains
        ## peq is the probability of each strain from a single sample
        seqrloop <- 1:indr 
        peq <- exp(-seqrloop*alpha) / sum(exp(-seqrloop*alpha))
        ## all possible orders of 1:R0 things
        matords <- rbind(1:r0,allPerms(r0,control=how()))
        dp <- max(1,dim(matords)[1])
        matposs <- matrix(t(vecn[matords]),ncol=r0,nrow=dp,byrow=F)
        ## combn gives the possible combinations of the R true infections
        ##    that give R0 observed infections
        matind <- t(combn(indr,r0))
        di <- max(1,dim(matind)[1])
        ## now need to make one row for every combination or matind and matposs
        ## combination will be all possible ways to observe r0 strains
        ## bigind gives which strains observed
        bigind <- matrix(as.vector(t(matind)),nrow=di*dp,ncol=r0,byrow=T)
        ## bigposs says how many of each of those named in bigind
        if (dp == 1 ) {   ## placeholder only need if r0=1
          bigposs <- matrix(matposs,nrow=di,ncol=r0,byrow=T) 
        }
        if (dp > 1 ) {
          bigposs <- matrix(matposs[1,],nrow=di,ncol=r0,byrow=T)
          for (pp in 2:dp) {  addon <- matrix(matposs[pp,],nrow=di,ncol=r0,byrow=T)
             bigposs <- rbind(bigposs,addon) } }
        ## the product of each row of pmat is the probability of that combination
        ##   of strains and numbers of strains making up the observed number
        ## these first part is the probabilities in order
        matpeq <- matrix(peq[bigind],nrow=di*dp,ncol=r0,byrow=F)
        ## choose only handles two gps at a time, so used a product of choose fncs
        ## eg (n choose n1,n2,n3) = (n choose n1,(n2+n3)) * ((n2+n3) choose (n2,n3))
        ## vector leftt below is the numerator for the subsequent choose functions
        ##   after each element in vecn is sequentially chosen
        leftt <- rev(cumsum(rev(vecn)))
        choosef <- prod(choose(leftt,vecn))
        ## next exponentiate them for the number of times drawn
        pmat <- (matpeq)^bigposs 
        ## and compute probability for this combo as the product
        jointprobs <- apply(pmat,1,prod)
        matprobs[indn,indr] <- choosef*sum(jointprobs)
      } ## end loop over indr
    }  ## end loop over indn
    
    ## now finishing up E step
    ## recall matprobs is f(r0|r,n)
    ## so tosum is f(r0,r|n) = f(r0|r,n) * f(r)
    tosum <- matprobs * (lam^matr / factorial(matr))    
    ## and condprobs is f(r|r0,n) = f(r0,r|n) / f(r0|n)
    ##   where f(r0|n) = sum over r of f(r0,r|n) using colSums
    condprobs <- tosum/(rowSums(tosum)%*%t(rep(1,maxr)))

        
    #################### M step ##############################
    ## to solve this, make left side equal to right side by grid search
    ## the left side of the equation uses the condition probs
    ##    and the right side has to do with lambda
    lside <- sum(matr*condprobs)
    rside <- function(lam,N) { N*lam*exp(lam)/(exp(lam)-1) }
    lamtest <- rside(lam=lamtry,N=N)
    ## once you compute the right side of the equation at all these values
    ##    then pick the one that makes left and ride sides closest
    lamcurr <- lamtry[abs(lamtest-lside)==min(abs(lamtest-lside))]
    #cat("ntries=",ntries,'\n')
    #print(lamcurr)
  }  ## end while loop
  ## compute expected infection pp and probability of superinf
  expectr <- lamcurr/(1-exp(-lamcurr))
  ## p(super) = p(r>=2) = 1 - p(r=1)
  psuper <- 1 - lamcurr/(exp(lamcurr)-1)
  list(lamcurr=lamcurr,laminit=laminit,lamraw=lamraw,expectr=expectr,
       psuper=psuper,condprobs=condprobs,conv=(ntries<150))  
} ## end of function
  

