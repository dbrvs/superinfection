## examples
source('EM_superinfection_anym.R')

## do an example from HCV paper
data = c(150,4,2,0,0)
naive = sum(data/sum(data)*c(1:length(data)))

 al_list <- seq(0,3,.25)
 est_Ps <- seqq
 est_R <- seqq
 for (aa in 1:length(al_list)) 
 {
 est_Ps[aa] <- EM_superinfection_anym(length(data),data,alpha=al_list[aa])$psuper  
 est_R[aa] <- EM_superinfection_anym(length(data),data,alpha=al_list[aa])$expectr  
 }

## do examples with n=2 from Redd survey
N2=c(3. ,2. ,1,3 ,13,7 ,1 )
N =c(196,126,7,78,58,36,16)

 N1=N-N2
 
 Redd_data=cbind(N1,N2)

 est_Ps <- N
 est_R <- N
 
 for (i in 1:length(N)) 
 { 	
 	est_Ps[i] = EM_superinfection_anym(n=2,vecn=c(N1[i],N2[i]),alpha=2)$psuper
 }
