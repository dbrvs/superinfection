## short script for HIV examples (or any other many strain per sample virus)
# DBR 12/6/17

source('HIV_superinfection.R') #import the computational method

#generate some example data, could be commented out to use other data

lam_sim=1.5 #simulated ZTP parameter
alpha_sim=2  #simulated superinfection parameter
n=200 #number of samples per person, the number of sequences
maxr=20 #maximum richness possible
N=10 #number of trial ppts, number of rows in data

#make a sample trial population from ZTP
trial <- rpois(10^3,lam_sim) #generate poisson numbers, many more than number of actual ppts in trial
true_trial = trial[trial>0][1:N] #only interested in ppts with more than 1 strain detected 
 
maxRobs <- max(true_trial) #the maximum observed richness

counts <- table(c(1:maxRobs,true_trial))-1 ; #the counts of number of participants found with n strains

# make data for all ppts who are not multiply infected
data <- matrix(c(n,rep(0,(maxRobs-1))),nrow=counts[1],ncol=maxRobs,byrow=T)

#loop through all multi-infected ppts do the sampling from multinomial
for (indr0 in 2:maxRobs) { 
  if (counts[indr0] > 0) {
    
    psamp <- exp(-alpha_sim*(1:indr0))/sum(exp(-alpha_sim*(1:indr0))) #the exponential sampling probability for each strain beyond 1
    
    ppt_counts = t(rmultinom(counts[indr0],size=n,prob=psamp))
    
    ppt_counts <- cbind(ppt_counts,matrix(0,nrow=counts[indr0],ncol=(maxRobs-indr0))) #make sure number of rows is right
    
    data <- rbind(data,ppt_counts) #add that ppts counts to the total counts
  }
}

#now run the inference on the simulated data using the EM algorithm
alpha_guess=alpha_sim
EM_sol <- HIVsuperinfection(data,alpha_guess)

cat("true avg richness=",lam_sim/(1-exp(-lam_sim)),'\n')
cat("naive avg richness=",sum(counts*c(1:length(counts))/N),'\n')
cat("EM inferred avg richness=",EM_sol$expectr,'\n')

#EM_sol$psuper #can also look at prevalence of superinfection
