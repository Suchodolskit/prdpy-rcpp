#path <- "~/Desktop/untitled folder 2/prdpy2-master/"
path <- "C:\\Users\\Tomasz\\Desktop\\PRDPY-project2\\"

Sys.setenv("PKG_CXXFLAGS"="-std=c++14")

library("fpc")
library("mclust")
library("dendextend")

Rcpp::sourceCpp(paste(path,'spectral_aux.cpp',sep=""))

Laplacian_eigen <- function(G, k)
{
  ev <- eigen(diag(rowSums(G))-G)
  k2 <- ev$vectors[,order(ev$values)][,1:k+1]
}

spectral_clustering <- function(X,k,M)
{
  k6<-Laplacian_eigen(Mnn_graph(Mnn(X,M)),k)
  kmeans(Laplacian_eigen(Mnn_graph(Mnn(X,M)),k),k)
}


DistributedCircle <- function(n,R,delta)
{
  angle <- seq(0 , 2*pi, (2*pi/(n-1)))
  
  r <- runif(length(angle), min = R-delta, max = R+delta)
  
  x <- r *sin(angle)
  y <- r *cos(angle)
  m <- matrix(c(x,y),nrow = length(x))
  colnames(m) <- c("x","y")
  return(m)
}

circlesClust <- function(n)
{
  c1 <- DistributedCircle(n,5,10)
  c2 <- DistributedCircle(n,100,10)
  c3 <- DistributedCircle(n,200,10)
  c4 <- DistributedCircle(n,300,10)
  return(rbind(c1,c2,c3,c4))
}

spiral <- function(n,R,delta, beginAngle = 0, endAngle = 2*pi)
{
  angle <- seq(beginAngle , endAngle, ((endAngle-beginAngle)/(n-1)))
  
  r <- runif(length(angle), min = R-delta, max = R+delta) + seq(10,100,length.out = length(angle))
  
  x <- r *sin(angle)
  y <- r *cos(angle)
  m <- matrix(c(x,y),nrow = length(x))
  colnames(m) <- c("x","y")
  return(m)
}

twoSpirals <- function(n)
{
  s1 <- spiral(n,5,1.5)
  s2 <- spiral(n,5,1.5,pi,3*pi)
  return(rbind(s1,s2))
}

clust <- function(x=0, y=0, n=100, sd=1)
{
  xx <- rnorm(n, x, sd)
  yy <- rnorm(n, y, sd)
  m <- matrix(c(xx,yy),nrow = length(xx))
  colnames(m) <- c("x","y")
  return(m)
}

fiveClusts <- function(n)
{
  c1 <- clust(0  ,0  ,n,1)
  c2 <- clust(10 ,10 ,n,1)
  c3 <- clust(10 ,-10,n,1)
  c4 <- clust(-10,10 ,n,1)
  c5 <- clust(-10,-10,n,1)
  return(rbind(c1,c2,c3,c4,c5))
}


plotManyClusts <- function(cl, n, color)
{
  plot(cl[,1],cl[,2],col = rep(color, each = (length(cl)/2/n)),
  main = "Zbiory skupieñ", xlab = "x", ylab = "y")
}


clust1 <- fiveClusts(80)
plotManyClusts(clust1,5,c("red","green","blue","black","magenta"))
indexes1 <- rep(c(1,2,3,4,5),each = 80)

clust2 <- twoSpirals(200)
plotManyClusts(clust2, 2, c("red","blue"))
indexes2 <- rep(c(1,2),each = 200)

clust3 <- circlesClust(100)
plotManyClusts(clust3,4,c("red","green","blue","black"))
indexes3 <- rep(c(1,2,3,4),each = 100)

adjRandIndex <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(adjRandIndex) <- c("5 skupieñ", "2 spirale", "4 ko³a")
FMIndex <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(FMIndex) <- c("5 skupieñ", "2 spirale", "4 ko³a")


for(val in 1:100)
{
  c1 <- spectral_clustering(clust1,5,val)[[1]]
  c2 <- spectral_clustering(clust2,5,val)[[1]]
  c3 <- spectral_clustering(clust3,5,val)[[1]]
  
  adjRandIndex[val,1] =adjustedRandIndex(c1, indexes1)
  adjRandIndex[val,2] =adjustedRandIndex(c2, indexes2)
  adjRandIndex[val,3] =adjustedRandIndex(c3, indexes3)
  
  FMIndex[val,1] =FM_index(c1, indexes1)
  FMIndex[val,2] =FM_index(c2, indexes2)
  FMIndex[val,3] =FM_index(c3, indexes3)
  
  print(val)
}







lastwd <- getwd()
setwd(paste(path,"\\my",sep=""))


write.csv(data.frame(clust1),"set1.data")
write.csv(data.frame(clust2),"set2.data")
write.csv(data.frame(clust3),"set3.data")

write.csv(indexes1,"set1.labels0")
write.csv(indexes2,"set2.labels0")
write.csv(indexes3,"set3.labels0")

setwd(lastwd)



