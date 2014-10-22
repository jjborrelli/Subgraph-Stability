## Create matrices for each of the 13 possible 3 species configurations
# "s" denotes single links only

s1<-matrix(c(0,1,0,-1,0,1,0,-1,0),nrow=3,ncol=3)
s2<-matrix(c(0,1,1,-1,0,1,-1,-1,0),nrow=3,ncol=3)
s3<-matrix(c(0,1,-1,-1,0,1,1,-1,0),nrow=3,ncol=3)
s4<-matrix(c(0,1,1,-1,0,0,-1,0,0),nrow=3,ncol=3)
s5<-matrix(c(0,0,1,0,0,1,-1,-1,0),nrow=3,ncol=3)

# "d" denotes that double links are included
d1<-matrix(c(0,1,1,-1,0,1,-1,1,0),nrow=3,ncol=3)
d2<-matrix(c(0,1,1,1,0,1,-1,-1,0),nrow=3,ncol=3)
d3<-matrix(c(0,1,1,1,0,0,-1,0,0),nrow=3,ncol=3)
d4<-matrix(c(0,1,1,-1,0,0,1,0,0),nrow=3,ncol=3)
d5<-matrix(c(0,1,1,-1,0,1,1,-1,0),nrow=3,ncol=3)
d6<-matrix(c(0,1,1,1,0,1,1,1,0),nrow=3,ncol=3)
d7<-matrix(c(0,1,1,1,0,1,1,-1,0),nrow=3,ncol=3)
d8<-matrix(c(0,1,1,1,0,0,1,0,0),nrow=3,ncol=3)

# Create a list of all 13 matrices
mot.lst <- list(s1, s2, s3, s4, s5, d1, d2, d3, d4, d5, d6, d7, d8)
names(mot.lst) <- c("s1", "s2", "s3", "s4", "s5", "d1", "d2", "d3", "d4", "d5", "d6", "d7", "d8")

ran.unif <- function(motmat){
  newmat <- apply(motmat, c(1,2), function(x){
    if(x==1){runif(1, 0, 10)}else if(x==-1){runif(1, -1, 0)} else{0}
  })
  diag(newmat) <- -1
  return(newmat)
}

maxRE <- function(rmat){
  lam.max <- max(Re(eigen(rmat)$values))
  return(lam.max)
}

eig_analysis <- function(matrices){
  ranmat <- lapply(matrices, ran.unif)
  eigs <- t(sapply(ranmat, maxRE))
  return(eigs)
}

niche.model <- function(S,C){
  new.mat<-matrix(0,nrow=S,ncol=S)
  ci<-vector()
  niche<-runif(S,0,1)
  r<-rbeta(S,1,((1/(2*C))-1))*niche
  
  for(i in 1:S){
    ci[i]<-runif(1,r[i]/2,niche[i])
  }
  
  r[which(niche==min(niche))]<-.00000001
  
  for(i in 1:S){
    
    for(j in 1:S){
      if(niche[j]>(ci[i]-(.5*r[i])) && niche[j]<(ci[i]+.5*r[i])){
        new.mat[j,i]<-1
      }
    }
  }
  
  new.mat<-new.mat[,order(apply(new.mat,2,sum))]
  return(new.mat)
}

conversion <- function(tm){
  for(i in 1:nrow(tm)){
    for(j in 1:ncol(tm)){
      if(tm[i,j] == 1){tm[j,i] <- -1}
    }
  }
  return(tm)
}



n = 10
nm <- list()
for(i in 1:n){
  nm[[i]] <- niche.model(25, .1)
}

nmc <- lapply(nm, conversion)
nmg <- lapply(nm, graph.adjacency)

mots <- motif_counter(nmg, webs = 1:n)
motN <- as.matrix(mots[,2:14])
system.time(
nme <- replicate(1000, eig_analysis(nmc))
)

qss1 <- apply(nme, 2, function(x){sum(x<0)})
qss2 <- apply(nme, 2, function(x){sum(x>=0)})
m <- cbind(motN[,1:5], double = rowSums(motN[,6:13]))

mod1 <- glm(cbind(qss1, qss2)~m, family = "binomial")
summary(mod1)
qss

library(foreach)
library(doParallel)

cl <- makeCluster(3)
registerDoParallel(cl)
#.export = c("eig_analysis", "maxRE", "ran.unif")
system.time(
inPAR <- foreach(i = 1:100, .combine = "rbind") %dopar% {eig_analysis(nmc)}
)