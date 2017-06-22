
## author = amalia magaret + daniel reeves
## date started = 2017 Jan 29
## PI = Daniel Reeves
## objective = simulate data to check the naiive estimator

## to use
rm(list=ls())

## here are the parameters you need to choose
##   lambda = true Poisson parameter for number of infections
##   alpha = exponential abundance parameter
##   n = samples per person

#lam_list = c(0.0201,0.0403 ,0.1017,0.2071,0.4308,0.6755)
#Ps = c(0.01001633,0.02001466,0.04998824,0.09997835,0.19998191,0.30001105)

#parameters to pick!
alpha=2
lamuse=0.7
nN=2000
nuse=2  

#calculate a few things
npers=round(nN/nuse)
avgR = lamuse/(1-exp(-lamuse))

nsimulations <- 200 ## number of times to simulate and estimate lambda
n_est=0
for (nsim in 1:nsimulations) {
    ## make simulated data, use first 500 persons infected
    rinf <- rpois(80000,lamuse)
    ## you will get 500 b/c start with 80000 and lowest superuse is 1% so 800
    rpos <- rinf[rinf>0][1:npers]
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
    
	naive = sum(Nr0/sum(Nr0)*c(1:length(Nr0)))
		
	n_est = n_est + (naive-avgR)/avgR
		
    } ## end loop over nsimulations

print('average percent error in nsimulations')
print(n_est/nsimulations*100)
#return(n_est/nsimulations)

