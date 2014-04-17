
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
d4<-matrix(c(0,1,1,-1,0,0,1,0,0),nrow=3,ncol=3)
d5<-matrix(c(0,1,1,-1,0,1,1,-1,0),nrow=3,ncol=3)
d6<-matrix(c(0,1,1,1,0,1,1,1,0),nrow=3,ncol=3)
d7<-matrix(c(0,1,1,1,0,1,1,-1,0),nrow=3,ncol=3)
d8<-matrix(c(0,1,1,1,0,0,1,0,0),nrow=3,ncol=3)

# Create a list of all 13 matrices
mot.lst <- list(s1, s2, s3, s4, s5, d1, d2, d3, d4, d5, d6, d7, d8)
names(mot.lst) <- c("s1", "s2", "s3", "s4", "s5", "d1", "d2", "d3", "d4", "d5", "d6", "d7", "d8")

## Functions to get eigenvalues of randomly sampled matrices
ran.unif <- function(motmat){
  newmat <- apply(motmat, c(1,2), function(x){
    if(x==1){runif(1, 0, 10)}else if(x==-1){runif(1, -1, 0)} else{0}
  })
  diag(newmat) <- runif(1, -1, 0)
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
motcount <- read.csv("~/Desktop/GitHub/Ecological-Networks/FoodWebs/Tables/zscore_both.csv", row.names = 1)
df.freq <- data.frame(motcount)

mfreq <- melt(df.freq[,names(sorted)])

g <- ggplot(mfreq, aes(x = variable, y = value)) + geom_boxplot() 
g <- g + geom_line(data = sort.qss, aes(x = 1:13, y = sorted), size = 1.5)
g <- g + geom_point(data = sort.qss, aes(x = 1:13, y = sorted), size = 4, col = "blue")
g <- g + geom_line(data = sort.qss.l, aes(x = 1:13, y = sorted.l), size = 1.5)
g <- g + geom_point(data = sort.qss.l, aes(x = 1:13, y = sorted.l), size = 4, col = "darkred")
g <- g + geom_hline(aes(yintercept = 0), lty = 2, col = "red")
g + xlab("Subgraph") + ylab("Frequency")


## After loading food web data

conversion <- function(tm){
  for(i in 1:nrow(tm)){
    for(j in 1:ncol(tm)){
      if(tm[i,j] == 1){tm[j,i] <- -1}
    }
  }
  return(tm)
}

webmats <- lapply(web.matrices, conversion)

reps = 1000
system.time(
  fw.stab.u <- eig.analysis(reps, webmats, mode = "unif")
)

fwQSSu <- apply(fw.stab.u, 2, function(x){sum(x<0)/reps})

system.time(
  fw.stab.l <- eig.analysis(reps, webmats, mode = "lnorm")
)

fwQSSl <- apply(fw.stab.l, 2, function(x){sum(x<0)/reps})
