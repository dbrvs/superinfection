
## author = Amalia Magaret
## date started = 17 Aout 2016 (Gabonese Independence Day)
## objective = give up on recursive loops and just make
## separate functions for each r
## this one makes up all possible combinations of possible strains
##   out of r true strains so that you have observed r0
##   ie if you have r=3 true strains and you observe r0=2
##   then you could have observed strains 1, 2 or 1, 3 or 2, 3
## aha think about if order is important

mkmatind <- function(r,r0) {
  matind <- matrix(0,nrow=choose(r,r0),ncol=r0)
  rownum <- 1
  if (r0==1) {
    for (i in 1:r) { 
        matind[rownum,] <- i
        rownum <- rownum + 1 }  }
  if (r0==2) {
    for (i in 1:(r-1)){ 
      for (j in (i+1):r) { 
        matind[rownum,] <- c(i,j)
        rownum <- rownum + 1 } }  }
  if (r0==3) {
    for (i in 1:(r-2)){ 
      for (j in (i+1):(r-1)) { 
        for (l in (j+1):r) { 
          matind[rownum,] <- c(i,j,l)
          rownum <- rownum + 1 } } } }
  if (r0==4) {
    for (i in 1:(r-3)){ 
      for (j in (i+1):(r-2)) { 
        for (l in (j+1):(r-1)) { 
          for (p in (l+1):r) { 
            matind[rownum,] <- c(i,j,l,p)
            rownum <- rownum + 1 } } } } }
  if (r0==5) {
    for (i in 1:(r-4)){ 
      for (j in (i+1):(r-3)) { 
        for (l in (j+1):(r-2)) { 
          for (p in (l+1):(r-1)) { 
            for (q in (p+1):r) { 
              matind[rownum,] <- c(i,j,l,p,q)
              rownum <- rownum + 1 } } } } } }
  if (r0==6) {
    for (i in 1:(r-5)){ 
      for (j in (i+1):(r-4)) { 
        for (l in (j+1):(r-3)) { 
          for (p in (l+1):(r-2)) { 
            for (q in (p+1):(r-1)) { 
              for (s in (q+1):r) { 
                matind[rownum,] <- c(i,j,l,p,q,s)
                rownum <- rownum + 1 } } } } } } }
  if (r0==7) {
    for (i in 1:(r-6)){ 
      for (j in (i+1):(r-5)) { 
        for (l in (j+1):(r-4)) { 
          for (p in (l+1):(r-3)) { 
            for (q in (p+1):(r-2)) { 
              for (s in (q+1):(r-1)) { 
                for (t in (s+1):r) { 
                  matind[rownum,] <- c(i,j,l,p,q,s,t)
                  rownum <- rownum + 1 } } } } } } } }
  if (r0==8) {
    for (i in 1:(r-7)){ 
      for (j in (i+1):(r-6)) { 
        for (l in (j+1):(r-5)) { 
          for (p in (l+1):(r-4)) { 
            for (q in (p+1):(r-3)) { 
              for (s in (q+1):(r-2)) { 
                for (t in (s+1):(r-1)) { 
                  for (u in (t+1):r) { 
                    matind[rownum,] <- c(i,j,l,p,q,s,t,u)
                    rownum <- rownum + 1 } } } } } } } } }
  if (r0==9) {
    for (i in 1:(r-8)){ 
      for (j in (i+1):(r-7)) { 
        for (l in (j+1):(r-6)) { 
          for (p in (l+1):(r-5)) { 
            for (q in (p+1):(r-4)) { 
              for (s in (q+1):(r-3)) { 
                for (t in (s+1):(r-2)) { 
                  for (u in (t+1):(r-1)) { 
                    for (v in (u+1):r) { 
                      matind[rownum,] <- c(i,j,l,p,q,s,t,u,v)
                      rownum <- rownum + 1 } } } } } } } } } }
  if (r0==10) {
    for (i in 1:(r-9)){ 
      for (j in (i+1):(r-8)) { 
        for (l in (j+1):(r-7)) { 
          for (p in (l+1):(r-6)) { 
            for (q in (p+1):(r-5)) { 
              for (s in (q+1):(r-4)) { 
                for (t in (s+1):(r-3)) { 
                  for (u in (t+1):(r-2)) { 
                    for (v in (u+1):(r-1)) { 
                      for (w in (v+1):r) { 
                        matind[rownum,] <- c(i,j,l,p,q,s,t,u,v,w)
                        rownum <- rownum + 1 } } } } } } } } } } }
  return(matind)
}


