
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

# This shows that some of my motif matrices are WRONG
# NEEDS TO BE FIXED!!!!
testA2 <- lapply(mot.lst, graph.adjacency)
motif_counter(testA2, webs = names(mot.lst))

## Functions to get eigenvalues of randomly sampled matrices
ran.unif <- function(motmat){
  newmat <- apply(motmat, c(1,2), function(x){
    if(x==1){runif(1, 0, 5)}else if(x==-1){runif(1, -1, 0)} else{0}
  })
  diag(newmat) <- -1
  return(newmat)
}

### using a lognormal
ran.lnorm <- function(motmat){
  newmat <- apply(motmat, c(1,2), function(x){
    if(x==1){rlnorm(1, -1, 1)}else if(x==-1){-rlnorm(1, -5, 1)} else{0}
  })
  diag(newmat) <- -rlnorm(1, -5, 1)
  return(newmat)
}

### Calculate largest eigenvalue (the real part)
maxRE <- function(rmat){
  lam.max <- max(Re(eigen(rmat)$values))
  return(lam.max)
}

### Wrap previous two functions together
eig.analysis <- function(n, matrices, mode = "unif"){
  cols <- length(matrices)
  rows <- n
  eigenMATRIX <- matrix(0, nrow = rows, ncol = cols)
  for(i in 1:n){
    if(mode == "unif"){
      ranmat <- lapply(matrices, ran.unif)
    }else if(mode == "lnorm"){
      ranmat <- lapply(matrices, ran.lnorm)
    }
    eigs <- sapply(ranmat, maxRE)
    eigenMATRIX[i,] <- eigs
  }
  return(eigenMATRIX)
}

## Run with uniform
n <- 5000
system.time(
  mot.stab<- eig.analysis(n, mot.lst, mode = "unif")
)

colnames(mot.stab) <- names(mot.lst)

## Run with lognormal
n <- 10000
system.time(
  mot.stab.l<- eig.analysis(n, mot.lst, mode = "lnorm")
)

colnames(mot.stab.l) <- names(mot.lst)

require(reshape2)
m <- melt(mot.stab)
m.l <- melt(mot.stab.l)
require(ggplot2)

ggplot(m, aes(x = Var2, y = value)) + geom_boxplot() + geom_hline(aes(yintercept = 0))
ggplot(m.l, aes(x = Var2, y = value)) + geom_boxplot() + geom_hline(aes(yintercept = 0))

# Unif
mot.qss <- apply(mot.stab, 2, function(x){sum(x<0)/n})
sorted <- sort(mot.qss, decreasing = T)
sort.qss <- data.frame(sorted, names = names(sorted))
#Lnorm
mot.qss.l <- apply(mot.stab.l, 2, function(x){sum(x<0)/n})
sorted.l <- sort(mot.qss.l, decreasing = T)
sort.qss.l <- data.frame(sorted.l, names = names(sorted))

# Read in and reshape subgraph frequency data
motcount <- read.csv("C:/Users/borre_000/Desktop/GitHub/Ecological-Networks/FoodWebs/Tables/zscore_both.csv", row.names = 1)
df.freq <- data.frame(motcount)

mfreq <- melt(df.freq[,names(sorted)])

g <- ggplot(mfreq, aes(x = variable, y = value)) + geom_boxplot() 
g <- g + geom_line(data = sort.qss, aes(x = 1:13, y = sorted), size = 1.5)
g <- g + geom_point(data = sort.qss, aes(x = 1:13, y = sorted), size = 4, col = "blue")
g <- g + geom_line(data = sort.qss.l, aes(x = 1:13, y = sorted.l), size = 1.5)
g <- g + geom_point(data = sort.qss.l, aes(x = 1:13, y = sorted.l), size = 4, col = "darkred")
g <- g + geom_hline(aes(yintercept = 0), lty = 2, col = "red")
g + xlab("Subgraph") + ylab("Frequency")


#######

motif.df <- read.table("C:/Users/borre_000/Desktop/GitHub/Ecological-Networks/FoodWebs/Tables/motifCOUNTS.csv", header = T, sep = ",", row.names = 1)
sub.counts <- motif.df[,2:14]
row.names(sub.counts) <- motif.df[,1]

permutes_prob <- function (mat, iter = 100, FUN, ...){
  cS <- colSums(mat)/nrow(mat)
  rS <- rowSums(mat)/ncol(mat)
  pmat <- matrix(0, nrow = nrow(mat), ncol = ncol(mat))
  for (i in 1:nrow(mat)) {
    for (j in 1:ncol(mat)) {
      pmat[i, j] <- sum((rS[i] + cS[j])/2)
    }
  }
  res <- list()
  for (q in 1:iter) {
    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat))
    for (i in 1:nrow(mat)) {
      for (j in 1:ncol(mat)) {
        mat2[i, j] <- rbinom(1, 1, prob = pmat[i, j])
      }
    }
    res[[q]] <- FUN(list(graph.adjacency(mat2)), ...)
  }
  return(res)
}


require(data.table)
require(doSNOW)
require(foreach)

clus <- makeCluster(3, "SOCK")
clusterExport(cl = clus, list = c("rbindlist", "permutes_prob", "motif_counter"))
registerDoSNOW(cl)

system.time(
pmat <- parLapply(cl, web.matrices, function(x){rbindlist(permutes_prob(x, iter = 1000, FUN = motif_counter, webs = 1))})
)
#45 min for 1000 iterations

stopCluster(clus)

nullmean <- t(sapply(pmat, function(x){apply(x, 2, mean)}))[,2:14]
nullsd <- t(sapply(pmat, function(x){apply(x, 2, sd)}))[,2:14]
nullmot <- (sub.counts - nullmean)/nullsd
nullmot[is.na(nullmot)] <- 0

nullmot <- apply(nullmot, 2, function(x){x/sqrt(sum(x^2))})

boxplot(nullmot[,names(sorted)])

mfreq.rcp <- melt(nullmot[,names(sorted)])

p1 <- ggplot(mfreq.rcp, aes(x = Var2, y = value)) + geom_boxplot() 
p1 <- p1 + geom_line(data = sort.qss, aes(x = 1:13, y = sorted), size = 1.5)
p1 <- p1 + geom_point(data = sort.qss, aes(x = 1:13, y = sorted), size = 4, col = "blue")
p1 <- p1 + geom_line(data = sort.qss.l, aes(x = 1:13, y = sorted.l), size = 1.5)
p1 <- p1 + geom_point(data = sort.qss.l, aes(x = 1:13, y = sorted.l), size = 4, col = "darkred")
p1 <- p1 + geom_hline(aes(yintercept = 0), lty = 2, col = "red")
p1 + xlab("Subgraph") + ylab("Frequency")


# Subgraph models

conversion <- function(tm){
  for(i in 1:nrow(tm)){
    for(j in 1:ncol(tm)){
      if(tm[i,j] == 1){tm[j,i] <- -1}
    }
  }
  return(tm)
}

nm <- list()
for(i in 1:10){
  nm[[i]] <- niche.model(50, .15)
}

nmL <- lapply(nm, graph.adjacency)
mcN <- motif_counter(nmL, webs = 1:10)
motN <- as.matrix(mcN[,2:14])
motN2 <- matrix(c(motN[,1:5], colSums(motN[,6:13])), nrow = 10)

nm2 <- lapply(nm, conversion)
mcEIG <- eig.analysis(1000, nm2, mode = "unif")
qss <- apply(mcEIG, 2, function(x){sum(x<0)/1000})
qss
m <- apply(mcEIG, 2, min)


summary(betareg(qss~motN2))

setwd("C:/Users/borre_000/Desktop/GitHub/Subgraph-Stability/")
#save.image("subgraphSTABILITY.Rdata")
load("subgraphSTABILITY.Rdata")
