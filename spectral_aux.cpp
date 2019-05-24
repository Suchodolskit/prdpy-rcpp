#include <Rcpp.h>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <functional>
#include <array>
using namespace Rcpp;
using namespace std;

class fun
{
  public:
    NumericVector distances;
    int k;
    fun(NumericVector d, int k)
    {
      distances = d;
      this->k = k;
    }
    int operator()(const int& x,const int& y)
    {
      return distances(x,k)<distances(y,k);
    }
};

NumericVector EuclidesDistanceRowWectors(NumericMatrix X, int row1, int row2)
{
  NumericVector v = NumericVector(0);
  double sumpow2 = 0;
  for(int i=0;i<X.ncol();i++)
  {
    sumpow2+=(X(row1,i)-X(row2,i))*(X(row1,i)-X(row2,i));
  }
  return NumericVector::create(sqrt(sumpow2));
}


// [[Rcpp::export]]
IntegerMatrix Mnn(NumericMatrix X, IntegerVector M) {
  IntegerMatrix ret = IntegerMatrix(X.nrow(),M[0]);
  
  NumericMatrix distances = NumericMatrix(X.nrow(),X.nrow());

  for(int i=0;i<X.nrow();i++)
  {
    for(int j=0;j<X.nrow();j++)
    {
      distances(i,j) = EuclidesDistanceRowWectors(X,i,j)(0);
    }
  }
  
  int k=0;
  
  fun f = fun(distances,k);

  for(k=0;k<X.nrow();k++)
  {
    f.k=k;
    
    IntegerVector indexes = IntegerVector(X.nrow());
    
    for(int i=0;i<X.nrow();i++)
    {
      indexes[i]=i;
    }
    
    sort(indexes.begin(),indexes.end(),f);
    
    
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
  }
  return ret;
}


void DFSrec(IntegerMatrix Nmatrix,int vertex, int* searched, int partNr, int& searchedcount)
{
  searchedcount++;
  searched[vertex] = partNr;
  for(int i=0;i<Nmatrix.ncol();i++)
  {
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
      Nmatrix(S(i,j),i)=1;
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
  }
  
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
  
  delete[] searched;
  return Nmatrix;
}


