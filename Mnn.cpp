#include <Rcpp.h>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <functional>
#include <array>
using namespace Rcpp;
using namespace std;


NumericVector EuclidesDistanceRowWectors(NumericMatrix X, IntegerVector row1, IntegerVector row2)
{
  NumericVector v = NumericVector(0);
  double sumpow2 = 0;
  for(int i=0;i<X.ncol();i++)
  {
    sumpow2+=(X(row1[0],i)-X(row2[0],i))*(X(row1[0],i)-X(row2[0],i));
  }
  return NumericVector::create(sqrt(sumpow2));
}


// [[Rcpp::export]]
IntegerMatrix Mnn(NumericMatrix X, IntegerVector M) {
  IntegerMatrix ret = IntegerMatrix(X.nrow(),M[0]);
  
  static NumericMatrix distances = NumericMatrix(X.nrow(),X.nrow());
  
 
  
  for(int i=0;i<X.nrow();i++)
  {
    for(int j=0;j<X.nrow();j++)
    {
      distances(i,j) = EuclidesDistanceRowWectors(X,IntegerVector::create(i),IntegerVector::create(j))(0);
    }
  }
  
  
  static int k=0;
  

  for(k=0;k<X.nrow();k++)
  {
    
    IntegerVector indexes = IntegerVector(X.nrow());
    
    for(int i=0;i<X.nrow();i++)
    {
      indexes[i]=i;
    }
    
    sort(indexes.begin(),indexes.end(),[](const int& x,const int& y){return distances(x,k)<distances(y,k);});
    
    int i=0;
    int j=0;
    while(j<M[0])
    {
      if(indexes[i]==k)i++;
      else
      {
        ret(k,j) = indexes[i];
        j++; 
        i++;
      }
      
    }
    //cout<<indexes<<endl;
  }
  return ret;
}


void DFSrec(IntegerMatrix Nmatrix,int vertex, int* searched, int partNr, int& searchedcount)
{
  //cout<<"searched: "<<vertex<<endl;
  searchedcount++;
  searched[vertex] = partNr;
  for(int i=0;i<Nmatrix.ncol();i++)
  {
    //cout<<vertex<<" "<<i<<" "<<Nmatrix(vertex,i)<<" "<<searched[i]<<""<<endl;
    if(Nmatrix(vertex,i)==1 && searched[i]==0)
    {
      DFSrec(Nmatrix,i,searched,partNr,searchedcount); 
    }    
  }
}

// [[Rcpp::export]]
IntegerMatrix Mnn_graph(IntegerMatrix S)
{
  IntegerMatrix Nmatrix = IntegerMatrix(S.nrow(),S.nrow());
  for(int i=0;i<S.nrow();i++)
  {
    for(int j=0;j<S.ncol();j++)
    {
      Nmatrix(i,S(i,j))=1;
    }
  }
  
  
  int* searched = new int[S.nrow()];
  
  for(int i=0;i<S.nrow();i++)
  {
    searched[i]=0;
  }
  
  int k=0;
  int num=0;
  while(k<S.nrow())
  {
    num++;
    for(int i=0;i<S.nrow();i++)
    {
      if(searched[i]==0)
      {
        DFSrec(Nmatrix,i,searched,num,k);
        break;
      }
    }

    //cout<<"k: "<<k<<endl;
  }
  
  //wiêcej ni¿ 1 sk³adowa-³¹czenie naiwne
  if(num>1)
  {
    int actnum=2;
    for(int i=1;i<S.nrow();i++)
    {
      if(searched[i]==actnum)
      {
        Nmatrix(0,i)=1;
        Nmatrix(i,0)=1;
        actnum++;
      }
      
    }
  }
  
  //for(int i=0;i<S.nrow();i++)
  //{
  //  cout<<searched[i]<<", ";
  //}
  //cout<<endl;
  
  delete[] searched;
  return Nmatrix;
}





//Sys.setenv("PKG_CXXFLAGS"="-std=c++14")
//Rcpp::sourceCpp('C:\\Users\\Tomasz\\Desktop\\PRDPY - project 2\\Mnn.Cpp')
//v <-rep(1:10,4)
//dim(v) <- c(10,4)

