
## Create matrices for each of the 13 possible 3 species configurations
# "s" denotes single links only

s1<-matrix(c(0,1,0,-1,0,1,0,-1,0),nrow=3,ncol=3)
s2<-matrix(c(0,1,1,-1,0,1,-1,-1,0),nrow=3,ncol=3)
s3<-matrix(c(0,1,-1,-1,0,1,1,-1,0),nrow=3,ncol=3)
s4<-matrix(c(0,0,1,0,0,1,-1,-1,0),nrow=3,ncol=3)
s5<-matrix(c(0,1,1,-1,0,0,-1,0,0),nrow=3,ncol=3)

# "d" denotes that double links are included
d1<-matrix(c(0,1,1,1,0,1,0,0,0),nrow=3,ncol=3)
d2<-matrix(c(0,1,1,0,0,1,0,1,0),nrow=3,ncol=3)
d3<-matrix(c(0,1,1,1,0,0,0,0,0),nrow=3,ncol=3)
d4<-matrix(c(0,1,1,0,0,0,0,1,0),nrow=3,ncol=3)
d5<-matrix(c(0,1,1,0,0,1,1,0,0),nrow=3,ncol=3)
d6<-matrix(c(0,1,1,1,0,1,1,1,0),nrow=3,ncol=3)
d7<-matrix(c(0,1,1,1,0,1,1,0,0),nrow=3,ncol=3)
d8<-matrix(c(0,1,1,1,0,0,1,0,0),nrow=3,ncol=3)

# Create a list of all 13 matrices
mot.lst <- list(s1, s2, s3, s4, s5, d1, d2, d3, d4, d5, d6, d7, d8)


