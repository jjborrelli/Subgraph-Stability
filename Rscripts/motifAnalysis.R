
## Create matrices for each of the 13 possible 3 species configurations
# "s" denotes single links only

s1<-matrix(c(0,1,0,-1,0,1,0,-1,0),nrow=3,ncol=3)
s2<-matrix(c(0,1,1,-1,0,1,-1,-1,0),nrow=3,ncol=3)
s3<-matrix(c(0,1,-1,-1,0,1,1,-1,0),nrow=3,ncol=3)
s4<-matrix(c(0,0,1,0,0,1,-1,-1,0),nrow=3,ncol=3)
s5<-matrix(c(0,1,1,-1,0,0,-1,0,0),nrow=3,ncol=3)

# "d" denotes that double links are included
d1<-matrix(c(0,1,1,1,0,1,-1,-1,0),nrow=3,ncol=3)
d2<-matrix(c(0,1,1,-1,0,1,-1,1,0),nrow=3,ncol=3)
d3<-matrix(c(0,1,1,1,0,0,-1,0,0),nrow=3,ncol=3)
#d4<-matrix(c(0,1,1,0,0,0,0,1,0),nrow=3,ncol=3)
d5<-matrix(c(0,1,1,-1,0,1,1,-1,0),nrow=3,ncol=3)
d6<-matrix(c(0,1,1,1,0,1,1,1,0),nrow=3,ncol=3)
d7<-matrix(c(0,1,1,1,0,1,1,-1,0),nrow=3,ncol=3)
d8<-matrix(c(0,1,1,1,0,0,1,0,0),nrow=3,ncol=3)

# Create a list of all 13 matrices
mot.lst <- list(s1, s2, s3, s4, s5, d1, d2, d3, d4, d5, d6, d7, d8)

eigen_unif <- function(m, params, self = -1){
  # For when I want to use uniform distribution
  # Params is dataframe of mean and standard deviation for relative impact of prey on pred 
  ## and pred on prey
  ev <- c()
  for(i in 1:nrow(m)){
    for (j in 1:nrow(m)){
      if(m[i, j] == 1){
        m[i, j] <- runif(1, params$pred1, params$pred2)
        m[j, i] <- runif(1, params$prey1, params$prey2)
      }
    }
  }
  diag(m) <- self
  ev <- max(Re(eigen(m)$values))
  return(ev)
}


ran.unif <- function(motmat){
  newmat <- apply(motmat, c(1,2), function(x){
    if(x==1){runif(1, 0, 10)}else if(x==-1){runif(1, -1, 0)} else{0}
  })
  diag(newmat) <- runif(3, -1, 0)
  return(newmat)
}

maxRE <- function(rmat){
  lam.max <- max(Re(eigen(rmat)$values))
  return(lam.max)
}

eig.analysis <- function(n, matrices){
  cols <- length(matrices)
  rows <- n
  eigenMATRIX <- matrix(nrow = rows, ncol = cols)
  for(i in 1:n){
    ranmat <- lapply(matrices, ran.unif)
    eigs <- sapply(ranmat, maxRE)
    eigenMATRIX[i,] <- eigs
  }
  return(eigenMATRIX)
}

system.time(
  mot.stab<- eig.analysis(10000, mot.lst)
)

boxplot(mot.stab)
abline(h=0)
ranmat <- lapply(mot.lst, ran.unif)



sapply(ranmat, maxRE)
