
## author = Amalia Magaret
## date started = 6 Sep 2016 
##  enumerate ways to split up n samples into r0 groups
##   ie if there are n=3 samples and you observe r0=2 strains,
##   then could get 1 of first strain and 2 of second
##      or 2 of first and 1 of the second strain
##   this one does assumes the strains are already specified in order
## use mk_matind to enumerate possible ways to get r0=2 observed
##      out of r true strains
library(permute)

mkmatposs <- function(r0,n) {
  matposs <- matrix(NA,nrow=1,ncol=r0)
  if (r0==1) {
      matposs <- matrix(rbind(matposs,n))  }  
  if (r0==2) {
    for (i in 1:(n-r0+1)){ 
      matposs <- rbind(matposs,c(i,n-i)) }  }  
  if (r0==3) {
    for (i in 1:(n-r0+1)){ 
      for (j in 1:(n-r0+2-i)) { 
        matposs <- rbind(matposs,c(i,j,n-(i+j))) } }  }  
  if (r0==4) {
    for (i in 1:(n-r0+1)){ 
      for (j in 1:(n-r0+2-i)) { 
        for (l in 1:(n-r0+3-(i+j))) { 
          matposs <- rbind(matposs,c(i,j,l,n-(i+j+l))) } } }  }  
  if (r0==5) {
    for (i in 1:(n-r0+1)){ 
      for (j in 1:(n-r0+2-i)) { 
        for (l in 1:(n-r0+3-(i+j))) { 
          for (p in 1:(n-r0+4-(i+j+l))) { 
            matposs <- rbind(matposs,c(i,j,l,p,n-(i+j+l+p))) } } } }  }
  if (r0==6) {
    for (i in 1:(n-r0+1)){ 
      for (j in 1:(n-r0+2-i)) { 
        for (l in 1:(n-r0+3-(i+j))) { 
          for (p in 1:(n-r0+4-(i+j+l))) { 
            for (q in 1:(n-r0+5-(i+j+l+p))) { 
              matposs <- rbind(matposs,c(i,j,l,p,q,n-(i+j+l+p+q))) } } } } }  }
  if (r0==7) {
    for (i in 1:(n-r0+1)){ 
      for (j in 1:(n-r0+2-i)) { 
        for (l in 1:(n-r0+3-(i+j))) { 
          for (p in 1:(n-r0+4-(i+j+l))) { 
            for (q in 1:(n-r0+5-(i+j+l+p))) { 
              for (s in 1:(n-r0+6-(i+j+l+p+q))) { 
                matposs <- rbind(matposs,c(i,j,l,p,q,s,n-(i+j+l+p+q+s))) } } } } } }  } 
  if (r0==8) {
    for (i in 1:(n-r0+1)){ 
      for (j in 1:(n-r0+2-i)) { 
        for (l in 1:(n-r0+3-(i+j))) { 
          for (p in 1:(n-r0+4-(i+j+l))) { 
            for (q in 1:(n-r0+5-(i+j+l+p))) { 
              for (s in 1:(n-r0+6-(i+j+l+p+q))) { 
                for (t in 1:(n-r0+7-(i+j+l+p+q+s))) { 
                  matposs <- rbind(matposs,c(i,j,l,p,q,s,t,n-(i+j+l+p+q+s+t))) } } } } } } }  } 
  if (r0==9) {
    for (i in 1:(n-r0+1)){ 
      for (j in 1:(n-r0+2-i)) { 
        for (l in 1:(n-r0+3-(i+j))) { 
          for (p in 1:(n-r0+4-(i+j+l))) { 
            for (q in 1:(n-r0+5-(i+j+l+p))) { 
              for (s in 1:(n-r0+6-(i+j+l+p+q))) { 
                for (t in 1:(n-r0+7-(i+j+l+p+q+s))) { 
                  for (u in 1:(n-r0+8-(i+j+l+p+q+s+t))) { 
                    matposs <- rbind(matposs,c(i,j,l,p,q,s,t,u,n-(i+j+l+p+q+s+t+u))) } } } } } } } }  } 
  if (r0==10) {
    for (i in 1:(n-r0+1)){ 
      for (j in 1:(n-r0+2-i)) { 
        for (l in 1:(n-r0+3-(i+j))) { 
          for (p in 1:(n-r0+4-(i+j+l))) { 
            for (q in 1:(n-r0+5-(i+j+l+p))) { 
              for (s in 1:(n-r0+6-(i+j+l+p+q))) { 
                for (t in 1:(n-r0+7-(i+j+l+p+q+s))) { 
                  for (u in 1:(n-r0+8-(i+j+l+p+q+s+t))) { 
                    for (v in 1:(n-r0+9-(i+j+l+p+q+s+t+u))) { 
                      matposs <- rbind(matposs,c(i,j,l,p,q,s,t,u,v,n-(i+j+l+p+q+s+t+u+v))) } } } } } } } } }  } 
  dd <- dim(matposs)[1]
  ## remove first row matrix
  matposs <- matposs[2:dd,]
  ## make sure it is still a matrix, if only one row
  if (length(matposs)==r0) {matposs <- matrix(matposs,nrow=1,ncol=r0)}
  
  ## now duplicate rows that have ambiguous ordering
  for (rn in 1:(dd-1)) {
    roww <- matposs[rn,] 
    ## already know sum of roww is n and length is r0
    ## repeat each row by number of unique orderings
    nreps <- factorial(n) / prod(factorial(roww))
    if (rn==1) { newposs <- matrix(matposs[rn,],nrow=nreps,ncol=r0,byrow=T) }  else {
      newposs <- rbind(newposs,matrix(matposs[rn,],nrow=nreps,ncol=r0,byrow=T))
    }
  }
  return(newposs)
}
  
