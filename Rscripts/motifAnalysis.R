
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

## Run
n <- 10000
system.time(
  mot.stab<- eig.analysis(n, mot.lst)
)

colnames(mot.stab) <- names(mot.lst)

require(reshape2)
m <- melt(mot.stab)
require(ggplot2)
ggplot(m, aes(x = Var2, y = value)) + geom_boxplot() + geom_hline(aes(yintercept = 0))


mot.qss <- apply(mot.stab, 2, function(x){sum(x<0)/n})
sorted <- sort(mot.qss, decreasing = T)
sort.qss <- data.frame(sorted, names = names(sorted))

# Read in and reshape subgraph frequency data
motcount <- read.csv("~/Desktop/GitHub/Ecological-Networks/FoodWebs/Tables/zscore_both.csv", row.names = 1)
df.freq <- data.frame(motcount)

mfreq <- melt(df.freq[,names(sorted)])

g <- ggplot(mfreq, aes(x = variable, y = value)) + geom_boxplot() 
g <- g + geom_line(data = sort.qss, aes(x = 1:13, y = sorted), size = 1.5)
g <- g + geom_point(data = sort.qss, aes(x = 1:13, y = sorted), size = 4, col = "blue")
g <- g + geom_hline(aes(yintercept = 0), lty = 2, col = "red")
g + xlab("Subgraph") + ylab("Frequency")

