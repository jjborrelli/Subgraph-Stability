
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
d3<-matrix(c(0,1,-1,1,0,0,1,0,0),nrow=3,ncol=3)
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
    if(x==1){runif(1, 0, 10)}else if(x==-1){runif(1, -1, 0)} else{0}
  })
  diag(newmat) <- runif(3, -1, 0)
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
n <- 10000
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
motcount <- read.csv("C:/Users/borre_000/Desktop/GitHub/Ecological-Networks/FoodWebs/Tables/zscore_both2.csv", row.names = 1)
df.freq <- data.frame(motcount)

mfreq <- melt(df.freq[,names(sorted)])
ord <- c("s5", "s4", "s1", "s2", "d3", "d1", "d4", "d2", "d8", "d5", "s3", "d7", "d6")

g <- ggplot(mfreq, aes(x = variable, y = value)) + geom_boxplot() 
#g <- g + geom_line(data = sort.qss, aes(x = 1:13, y = sorted), size = 1.5)
g <- g + geom_point(data = sort.qss, aes(x = 1:13, y = sorted), size = 4, col = "blue")
#g <- g + geom_line(data = sort.qss.l, aes(x = 1:13, y = sorted.l), size = 1.5)
#g <- g + geom_point(data = sort.qss.l, aes(x = 1:13, y = sorted.l), size = 4, col = "darkred")
g <- g + geom_hline(aes(yintercept = 0), lty = 2, col = "red")
g + xlab("Subgraph") + ylab(expression(over("Log(Frequency)","Max(Log(Frequency))")))

cor.test(apply(z.both[,names(sorted)], 2, function(x){mean(x)}), sorted)

temp.df <- data.frame(matrix(rep(mot.qss.l, 50), ncol = 13, byrow = T))
colnames(temp.df) <- names(mot.lst)
mstab <- melt(temp.df)

dat1 <- data.frame(sub = names(sorted), qss = sorted, freq = apply(sub.counts[,names(sorted)], 2, mean))
fit <- glm(round(cbind(10000*dat1$qss, 10000-10000*dat1$qss))~dat1$freq, family = "binomial")
plot(dat1$freq, dat1$qss)
points(fit$fitted.values~dat1$freq, pch = 20, col = "red")


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
pmat <- parLapply(cl, web.matrices, function(x){rbindlist(permutes_prob(x, iter = 200, FUN = motif_counter, webs = 1))})
)
#45 min for 1000 iterations

stopCluster(clus)

nullmean2 <- t(sapply(pmat, function(x){apply(x, 2, mean)}))[,2:14]
nullsd2 <- t(sapply(pmat, function(x){apply(x, 2, sd)}))[,2:14]
nullmot2 <- (sub.counts - nullmean2)/nullsd2
nullmot2[is.na(nullmot2)] <- 0

nullmot2a <- t(apply(nullmot2, 1, function(x){x/sqrt(sum(x^2))}))

boxplot(nullmot2a[,names(sorted)])

mfreq.rcp <- melt(nullmot2a[,names(sorted)])

p1 <- ggplot(mfreq.rcp, aes(x = Var2, y = value)) + geom_boxplot() 
#p1 <- p1 + geom_line(data = sort.qss, aes(x = 1:13, y = sorted), size = 1.5)
#p1 <- p1 + geom_point(data = sort.qss, aes(x = 1:13, y = sorted), size = 4, col = "blue")
#p1 <- p1 + geom_line(data = sort.qss.l, aes(x = 1:13, y = sorted.l), size = 1.5)
#p1 <- p1 + geom_point(data = sort.qss.l, aes(x = 1:13, y = sorted.l), size = 4, col = "darkred")
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



###### MY NULL MODEL ROWS AND COLUMN SUMS 

permutes_rc <- function(mat, iter = 100){
  
  pattern1 <- matrix(c(0,1,1,0), nrow = 2, ncol = 2)
  pattern2 <- matrix(c(1,0,0,1), nrow = 2, ncol = 2)
  count <- 0
  
  mat.list <- list()
  
  while(count < iter){
    srow <- sample(1:nrow(mat), 2)
    scol <- sample(1:ncol(mat), 2)
    
    test <- mat[srow, scol]
    
    if(sum(test == pattern1) == 4){
      count <- count + 1
      mat[srow, scol] <- pattern2
      mat.list[[count]] <- mat
      
      next
    } else if(sum(test == pattern2) == 4){
      count <- count + 1
      mat[srow, scol] <- pattern1
      mat.list[[count]] <- mat
      
      next
    } else {next}
  }
  
  matrices <- lapply(mat.list, as.matrix)
  return(permuted.matrices = matrices)
}

require(igraph)
pmot <- list()

system.time(
for(i in 1:length(web.matrices)){
  p <- permutes_rc(web.matrices[[i]], 1000)
  g <- lapply(p, graph.adjacency)
  pmot[[i]] <- motif_counter(g, 1:1000)
  cat(i, names(web.matrices)[i], "\n")
}
)

mus <- t(sapply(pmot, function(x){colMeans(x[,2:14])}))
sig <- t(sapply(pmot, function(x){apply(x[,2:14], 2, sd)}))
z <- (sub.counts - mus)/sig
zmat <- as.matrix(z)
zmat[is.nan(zmat)] <- 0
profile <- apply(zmat, 2, function(x){x/sqrt(rowSums(zmat^2))})

require(ggplot2)
require(reshape2)

norm.pro <- melt(profile[,names(sorted)])

p1 <- ggplot(norm.pro, aes(x = Var2, y = value)) + geom_boxplot() 
p1 <- p1 + geom_line(data = sort.qss, aes(x = 1:13, y = sorted), size = 1.5)
p1 <- p1 + geom_point(data = sort.qss, aes(x = 1:13, y = sorted), size = 4, col = "blue")
p1 <- p1 + geom_line(data = sort.qss.l, aes(x = 1:13, y = sorted.l), size = 1.5)
p1 <- p1 + geom_point(data = sort.qss.l, aes(x = 1:13, y = sorted.l), size = 4, col = "darkred")
p1 <- p1 + geom_hline(aes(yintercept = 0), lty = 2, col = "red")
p1 + xlab("Subgraph") + ylab("Frequency")

cor.test(sorted, apply(profile[,names(sorted)], 2, median))

d <- data.frame(s = names(sorted), q = sorted)
ggplot(d, aes(x = s, y = q)) + geom_point(size = 3) + 
  scale_x_discrete(limits = c("s1", "s4", "s5", "s2", "d4", "d3",
                              "s3", "d2", "d1", "d5", "d7", "d8", "d6")) +
  xlab("Subgraph") + ylab("Quasi Sign-Stability")

setwd("C:/Users/borre_000/Desktop/GitHub/Subgraph-Stability/")
#save.image("subgraphSTABILITY2.Rdata")
load("subgraphSTABILITY.Rdata")
