Sys.setenv("PKG_CXXFLAGS"="-std=c++14")

#path <- "~/Desktop/untitled folder 2/prdpy2-master/"
path <- "C:\\Users\\Tomasz\\Desktop\\PRDPY-project2\\"
Rcpp::sourceCpp(paste(path,'spectral_aux.cpp',sep=""))

#install.packages("kernlab")
#install.packages("fpc")
#install.packages("mclust")
#install.packages("dendextend")
#install.packages("stringi")
#install.packages("genie")
#install.packages("cluster")

library("kernlab")
library("fpc")
library("mclust")
library("dendextend")
library("stringi")
library("genie")
library("cluster")



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


ReadAllFiles <- function(path2)
{
  lastwd <- getwd()
  setwd(paste(path2,'pd2-zbiory-benchmarkowe',sep=""))
  
  file_list <- list.files(recursive=TRUE)
  
  tmp <-read.csv(file_list[1],sep = "\t")
  k <-matrix(unlist(tmp), ncol = length(tmp)) 
  
  objects <- list()
  
  k <- stri_subset_regex(file_list,".data.gz", negate=FALSE)
  
  
  for (val in k)
  {
    k2 <- stri_sub(val,to = stri_locate_first_fixed(val,".")[,1]-1)
    matrix1 <- read.csv(paste(k2[1],".data.gz", sep=""),sep = "\t")
    nc = ncol(matrix1)
    matrix1 <- mapply(matrix1, FUN=as.numeric)
    colnames(matrix1) <-NULL
    matrix2 <- matrix(unlist(matrix1), ncol = nc)
    labels <- read.csv(paste(k2[1],".labels0.gz", sep=""),sep = "\t")
    objects  <- append(objects,list(k2, matrix2, labels))
  }
  setwd(lastwd)
  return (objects)
}


adjRandIndex <- data.frame(matrix(ncol = 11, nrow = 0))
colnames(adjRandIndex) <- c("Spec_Clust", "ward.D", "ward.D2", "single", "complete", "mcquitty", "average", "median", "centroid", "Genie_useVpTree_F", "cluster_pam")
FMIndex <- data.frame(matrix(ncol = 11, nrow = 0))
colnames(FMIndex) <- c("Spec_Clust", "ward.D", "ward.D2", "single", "complete", "mcquitty", "average", "median", "centroid", "Genie_useVpTree_F", "cluster_pam")

data <- ReadAllFiles(path)

cl <- list()

for(val in 1:(length(data)/3))
{
  ma = max(data[[val*3]])
  mi = min(data[[val*3]])
  
  cl[[1]] <- spectral_clustering(data[[val*3-1]],(ma-mi+1),3)[[1]]
  cl[[2]] <- cutree(hclust(dist(data[[val*3-1]]),method = "ward.D")   , k = mi:ma)[,(ma-mi+1)]
  cl[[3]] <- cutree(hclust(dist(data[[val*3-1]]),method = "ward.D2")  , k = mi:ma)[,(ma-mi+1)]
  cl[[4]] <- cutree(hclust(dist(data[[val*3-1]]),method = "single")   , k = mi:ma)[,(ma-mi+1)]
  cl[[5]] <- cutree(hclust(dist(data[[val*3-1]]),method = "complete") , k = mi:ma)[,(ma-mi+1)]
  cl[[6]] <- cutree(hclust(dist(data[[val*3-1]]),method = "average")  , k = mi:ma)[,(ma-mi+1)]
  cl[[7]] <- cutree(hclust(dist(data[[val*3-1]]),method = "mcquitty") , k = mi:ma)[,(ma-mi+1)]
  cl[[8]] <- cutree(hclust(dist(data[[val*3-1]]),method = "median")   , k = mi:ma)[,(ma-mi+1)]
  cl[[9]] <- cutree(hclust(dist(data[[val*3-1]]),method = "centroid") , k = mi:ma)[,(ma-mi+1)]
  cl[[10]] <- cutree(hclust2(dist(data[[val*3-1]]),useVpTree = FALSE)  , k = mi:ma)[,(ma-mi+1)]
  if(val!=20&&val!=21)
  {
    cl[[11]] <- pam(data[[val*3-1]],(ma-mi+1))[[3]]
  }

  
  adjRandIndex[val,1] =adjustedRandIndex( cl[[1]],unlist(data[[val*3]]))
  adjRandIndex[val,2] =adjustedRandIndex( cl[[2]],unlist(data[[val*3]]))
  adjRandIndex[val,3] =adjustedRandIndex( cl[[3]],unlist(data[[val*3]]))
  adjRandIndex[val,4] =adjustedRandIndex( cl[[4]],unlist(data[[val*3]]))
  adjRandIndex[val,5] =adjustedRandIndex( cl[[5]],unlist(data[[val*3]]))
  adjRandIndex[val,6] =adjustedRandIndex( cl[[6]],unlist(data[[val*3]]))
  adjRandIndex[val,7] =adjustedRandIndex( cl[[7]],unlist(data[[val*3]]))
  adjRandIndex[val,8] =adjustedRandIndex( cl[[8]],unlist(data[[val*3]]))
  adjRandIndex[val,9] =adjustedRandIndex( cl[[9]],unlist(data[[val*3]]))
  adjRandIndex[val,10]=adjustedRandIndex( cl[[10]],unlist(data[[val*3]]))
  adjRandIndex[val,11]=adjustedRandIndex( cl[[11]],unlist(data[[val*3]]))
  
  FMIndex[val,1] =FM_index( cl[[1]],unlist(data[[val*3]]))
  FMIndex[val,2] =FM_index( cl[[2]],unlist(data[[val*3]]))
  FMIndex[val,3] =FM_index( cl[[3]],unlist(data[[val*3]]))
  FMIndex[val,4] =FM_index( cl[[4]],unlist(data[[val*3]]))
  FMIndex[val,5] =FM_index( cl[[5]],unlist(data[[val*3]]))
  FMIndex[val,6] =FM_index( cl[[6]],unlist(data[[val*3]]))
  FMIndex[val,7] =FM_index( cl[[7]],unlist(data[[val*3]]))
  FMIndex[val,8] =FM_index( cl[[8]],unlist(data[[val*3]]))
  FMIndex[val,9] =FM_index( cl[[9]],unlist(data[[val*3]]))
  FMIndex[val,10]=FM_index( cl[[10]],unlist(data[[val*3]]))
  FMIndex[val,11]=FM_index( cl[[11]],unlist(data[[val*3]]))
}

ReadMyFiles <- function(path2)
{
  lastwd <- getwd()
  setwd(paste(path2,'my',sep=""))
  
  file_list <- list.files(recursive=TRUE)
  
  tmp <-read.csv(file_list[1],sep = "\t")
  k <-matrix(unlist(tmp), ncol = length(tmp)) 
  
  objects <- list()
  
  k <- stri_subset_regex(file_list,".data", negate=FALSE)
  
  
  for (val in k)
  {
    k2 <- stri_sub(val,to = stri_locate_first_fixed(val,".")[,1]-1)
    matrix1 <- read.csv(paste(k2[1],".data", sep=""),sep = ",")
    print(matrix1)
    nc = ncol(matrix1)
    matrix1 <- mapply(matrix1, FUN=as.numeric)
    colnames(matrix1) <-NULL
    matrix2 <- matrix(unlist(matrix1), ncol = nc)
    labels <- read.csv(paste(k2[1],".labels0", sep=""),sep = ",")
    objects  <- append(objects,list(k2, matrix2, labels))
  }
  
  setwd(lastwd)
  return (objects)
}

data <- ReadMyFiles(path)

data[[3]]$X <- NULL
data[[6]]$X <- NULL
data[[9]]$X <- NULL

colnames(data[[3]]) <- NULL
colnames(data[[6]]) <- NULL
colnames(data[[9]]) <- NULL

data[[2]] <-data[[2]][,2:3]
data[[5]] <-data[[5]][,2:3]
data[[8]] <-data[[8]][,2:3]


for(val in 1:(length(data)/3))
{
  ma = max(data[[val*3]])
  mi = min(data[[val*3]])
  
  
  
  print(val*3-1)
  print(ma)
  print(mi)
  cl[[1]] <- spectral_clustering(data[[val*3-1]],(ma-mi+1),3)[[1]]
  cl[[2]] <- cutree(hclust(dist(data[[val*3-1]]),method = "ward.D")   , k = mi:ma)[,(ma-mi+1)]
  cl[[3]] <- cutree(hclust(dist(data[[val*3-1]]),method = "ward.D2")  , k = mi:ma)[,(ma-mi+1)]
  cl[[4]] <- cutree(hclust(dist(data[[val*3-1]]),method = "single")   , k = mi:ma)[,(ma-mi+1)]
  cl[[5]] <- cutree(hclust(dist(data[[val*3-1]]),method = "complete") , k = mi:ma)[,(ma-mi+1)]
  cl[[6]] <- cutree(hclust(dist(data[[val*3-1]]),method = "average")  , k = mi:ma)[,(ma-mi+1)]
  cl[[7]] <- cutree(hclust(dist(data[[val*3-1]]),method = "mcquitty") , k = mi:ma)[,(ma-mi+1)]
  cl[[8]] <- cutree(hclust(dist(data[[val*3-1]]),method = "median")   , k = mi:ma)[,(ma-mi+1)]
  cl[[9]] <- cutree(hclust(dist(data[[val*3-1]]),method = "centroid") , k = mi:ma)[,(ma-mi+1)]
  cl[[10]] <- cutree(hclust2(dist(data[[val*3-1]]),useVpTree = FALSE)  , k = mi:ma)[,(ma-mi+1)]
  cl[[11]] <- pam(data[[val*3-1]],(ma-mi+1))[[3]]

  
  adjRandIndex[val+43,1] =adjustedRandIndex( cl[[1]],unlist(data[[val*3]]))
  adjRandIndex[val+43,2] =adjustedRandIndex( cl[[2]],unlist(data[[val*3]]))
  adjRandIndex[val+43,3] =adjustedRandIndex( cl[[3]],unlist(data[[val*3]]))
  adjRandIndex[val+43,4] =adjustedRandIndex( cl[[4]],unlist(data[[val*3]]))
  adjRandIndex[val+43,5] =adjustedRandIndex( cl[[5]],unlist(data[[val*3]]))
  adjRandIndex[val+43,6] =adjustedRandIndex( cl[[6]],unlist(data[[val*3]]))
  adjRandIndex[val+43,7] =adjustedRandIndex( cl[[7]],unlist(data[[val*3]]))
  adjRandIndex[val+43,8] =adjustedRandIndex( cl[[8]],unlist(data[[val*3]]))
  adjRandIndex[val+43,9] =adjustedRandIndex( cl[[9]],unlist(data[[val*3]]))
  adjRandIndex[val+43,10]=adjustedRandIndex( cl[[10]],unlist(data[[val*3]]))
  adjRandIndex[val+43,11]=adjustedRandIndex( cl[[11]],unlist(data[[val*3]]))
  
  FMIndex[val+43,1] =FM_index( cl[[1]],unlist(data[[val*3]]))
  FMIndex[val+43,2] =FM_index( cl[[2]],unlist(data[[val*3]]))
  FMIndex[val+43,3] =FM_index( cl[[3]],unlist(data[[val*3]]))
  FMIndex[val+43,4] =FM_index( cl[[4]],unlist(data[[val*3]]))
  FMIndex[val+43,5] =FM_index( cl[[5]],unlist(data[[val*3]]))
  FMIndex[val+43,6] =FM_index( cl[[6]],unlist(data[[val*3]]))
  FMIndex[val+43,7] =FM_index( cl[[7]],unlist(data[[val*3]]))
  FMIndex[val+43,8] =FM_index( cl[[8]],unlist(data[[val*3]]))
  FMIndex[val+43,9] =FM_index( cl[[9]],unlist(data[[val*3]]))
  FMIndex[val+43,10]=FM_index( cl[[10]],unlist(data[[val*3]]))
  FMIndex[val+43,11]=FM_index( cl[[11]],unlist(data[[val*3]]))
}

adjRandIndex[1,12] <- "tmp"
colnames(adjRandIndex)[[12]] <- "Names"
adjRandIndex[44,12] <- data[[1]]
adjRandIndex[45,12] <- data[[4]]
adjRandIndex[46,12] <- data[[7]]





setwd("C:/Users/Tomasz/Desktop/PRDPY-project2")
write.csv(adjRandIndex,"adjRandIndex.data")
write.csv(FMIndex,"FMIndex.data")

adjRandIndex$X <- NULL

data <- ReadAllFiles(path)
for(val in 1:(length(data)/3))
{
  adjRandIndex[val,12] <- data[3*val-2]
}



adjRandIndex <- read.csv("adjRandIndex.data")
FMIndex <- read.csv("FMIndex.data")
adjRandIndex$X <- NULL
FMIndex$X <- NULL
