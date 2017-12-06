## short script for HSV examples, easily modifiable
# DBR 12/6/17

source('HSV_superinfection.R') #import the computational method

data = c(471,18) #number with n=1, n=2, and so on, right now set to values from example from HSV-2 Johnston et al. Plos Med 2017

n=length(data) #number of samples
N=sum(data) #number of participants

naive_Ps = 1-data[1]/N #naive superinfection prob
naive_avgR = sum(data/N*c(1:n)) #naive avg Richness

alpha=0 #superinfection parameter, must be>0 and alpha=0 means even strains, and most parsimonious estimate without extra data

EM_sol = HSV_superinfection(n,data,alpha)

## the values we want!
EM_Ps = EM_sol$psuper #EM inference of superinfection prob
EM_avgR = EM_sol$expectr #EM inference of avg Richness