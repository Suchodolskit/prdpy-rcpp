Sys.setenv("PKG_CXXFLAGS"="-std=c++14")
Rcpp::sourceCpp('C:\\Users\\Tomasz\\Desktop\\PRDPY - project 2\\Mnn.Cpp')
library("stringi")



v <-rep(c(1,2,3,10,11,12,20,21,22),4)
dim(v) <- c(9,4)




Laplacian_eigen <- function(G, k)
{
  diag(rowSums(G))-G;
  ev <- eigen(G)
  k2 <- ev$vectors[,order(ev$values)][,1:k+1]
}

spectral_clustering <- function(X,k,M)
{
  kmeans(Laplacian_eigen(Mnn_graph(Mnn(X,M)),k),k)
}


ReadAllFiles <- function()
{
  setwd("C:\\Users\\Tomasz\\Desktop\\PRDPY - project 2\\pd2-zbiory-benchmarkowe")
  
  file_list <- list.files(recursive=TRUE)
  
  
  tmp <-read.csv(file_list[1],sep = "\t")
  k <-matrix(unlist(tmp), ncol = length(tmp)) 
  
  objects <- list()
  
  k <- stri_subset_regex(file_list,".data.gz", negate=FALSE)

  
  for (val in k)
  {
    k2 <- stri_sub(val,to = stri_locate_first_fixed(val,".")[,1]-1)
    matrix <- read.csv(paste(k2[1],".data.gz", sep=""),sep = "\t")
    colnames(matrix) <-NULL
    matrix2 <- matrix(unlist(matrix), ncol = length(matrix))
    labels <- read.csv(paste(k2[1],".labels0.gz", sep=""),sep = "\t")
    objects  <- append(objects,list(k2, matrix2, labels))
  }
  return (objects)
}

p <- spectral_clustering(ddd[[2]],max(ddd[[3]]),2)

ddd[2][2]
head(ddd,5)

c2 <- spectral_clustering(ddd[[2]],2,2)
adjustedRandIndex(c2[[1]],unlist(ddd[[3]]))












