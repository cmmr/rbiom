// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>
using namespace Rcpp;
using namespace RcppParallel;

#include <cmath>
#include <cstring>
#include <algorithm>


// Ennumerate the different algorithms
#define BRAYCURTIS_W 1
#define BRAYCURTIS_U 2
#define EUCLIDEAN_W  3
#define EUCLIDEAN_U  4
#define MANHATTAN_W  5
#define MANHATTAN_U  6
#define JACCARD_W    7
#define JACCARD_U    8




//======================================================
// Bray-Curtis Weighted
// sum(abs(x-y))/sum(x+y)
//======================================================
template <typename InputIterator1, typename InputIterator2>
inline double BrayCurtis_W(InputIterator1 begin1, InputIterator1 end1, InputIterator2 begin2) {
  
  // set iterators to beginning of ranges
  InputIterator1 it1 = begin1;
  InputIterator2 it2 = begin2;
  
  double diffs = 0;
  double sums  = 0;
  
  // for each input item
  while (it1 != end1) {
    
    // take the value and increment the iterator
    double x = *it1++;
    double y = *it2++;
    
    // accumulate if appropriate
    sums += x + y;
    if (x > y) diffs += x - y;
    if (y > x) diffs += y - x;
  }
  
  // value to return
  double rval = diffs / sums;
  return rval;  
}


//======================================================
// Bray-Curtis Unweighted
// (sum(x>0)+sum(y>0)-2*sum(x&y))/(sum(x>0)+sum(y>0))
//======================================================
template <typename InputIterator1, typename InputIterator2>
inline double BrayCurtis_U(InputIterator1 begin1, InputIterator1 end1, InputIterator2 begin2) {
  
  
  // set iterators to beginning of ranges
  InputIterator1 it1 = begin1;
  InputIterator2 it2 = begin2;
  
  double  x_nonzero = 0;
  double  y_nonzero = 0;
  double xy_nonzero = 0;
  
  // for each input item
  while (it1 != end1) {
    
    // take the value and increment the iterator
    double x = *it1++;
    double y = *it2++;
    
    // accumulate if appropriate
    if (x > 0)           x_nonzero++;
    if (y > 0)           y_nonzero++;
    if (x > 0 && y > 0) xy_nonzero++;
  }
  
  // value to return
  double rval = (x_nonzero + y_nonzero - 2 * xy_nonzero) / (x_nonzero + y_nonzero);
  
  return rval;  
}


//======================================================
// Euclidean Weighted
// sqrt(sum((x-y)^2))
//======================================================
template <typename InputIterator1, typename InputIterator2>
inline double Euclidean_W(InputIterator1 begin1, InputIterator1 end1, InputIterator2 begin2) {
  
  // set iterators to beginning of ranges
  InputIterator1 it1 = begin1;
  InputIterator2 it2 = begin2;
  
  // value to return
  double rval = 0;
  
  // for each input item
  while (it1 != end1) {
    
    // take the value and increment the iterator
    double x = *it1++;
    double y = *it2++;
    
    // accumulate if appropriate
    rval += (x - y) * (x - y);
  }
  
  rval = sqrt(rval);
  
  return rval;  
}


//======================================================
// Euclidean Unweighted
// sqrt(sum(x>0)+sum(y>0)-2*sum(x&y))
//======================================================
template <typename InputIterator1, typename InputIterator2>
inline double Euclidean_U(InputIterator1 begin1, InputIterator1 end1, InputIterator2 begin2) {
  
  
  // set iterators to beginning of ranges
  InputIterator1 it1 = begin1;
  InputIterator2 it2 = begin2;
  
  double  x_nonzero = 0;
  double  y_nonzero = 0;
  double xy_nonzero = 0;
  
  // for each input item
  while (it1 != end1) {
    
    // take the value and increment the iterator
    double x = *it1++;
    double y = *it2++;
    
    // accumulate if appropriate
    if (x > 0)           x_nonzero++;
    if (y > 0)           y_nonzero++;
    if (x > 0 && y > 0) xy_nonzero++;
  }
  
  // value to return
  double rval = sqrt(x_nonzero + y_nonzero - 2 * xy_nonzero);
  
  return rval;  
}


//======================================================
// Manhattan Weighted
// sum(abs(x-y))
//======================================================
template <typename InputIterator1, typename InputIterator2>
inline double Manhattan_W(InputIterator1 begin1, InputIterator1 end1, InputIterator2 begin2) {
  
  // value to return
  double rval = 0;
  
  // set iterators to beginning of ranges
  InputIterator1 it1 = begin1;
  InputIterator2 it2 = begin2;
  
  // for each input item
  while (it1 != end1) {
    
    // take the value and increment the iterator
    double x = *it1++;
    double y = *it2++;
    
    // accumulate if appropriate
    if (x > y) rval += x - y;
    if (y > x) rval += y - x;
  }
  return rval;  
}


//======================================================
// Manhattan Unweighted
// sum(x>0)+sum(y>0)-2*sum(x&y)
//======================================================
template <typename InputIterator1, typename InputIterator2>
inline double Manhattan_U(InputIterator1 begin1, InputIterator1 end1, InputIterator2 begin2) {
  
  
  // set iterators to beginning of ranges
  InputIterator1 it1 = begin1;
  InputIterator2 it2 = begin2;
  
  double  x_nonzero = 0;
  double  y_nonzero = 0;
  double xy_nonzero = 0;
  
  // for each input item
  while (it1 != end1) {
    
    // take the value and increment the iterator
    double x = *it1++;
    double y = *it2++;
    
    // accumulate if appropriate
    if (x > 0)           x_nonzero++;
    if (y > 0)           y_nonzero++;
    if (x > 0 && y > 0) xy_nonzero++;
  }
  
  // value to return
  double rval = x_nonzero + y_nonzero - 2 * xy_nonzero;
  
  return rval;  
}


//======================================================
// Jaccard Weighted
// 2 * BrayCurtis_W / (1 + BrayCurtis_W)
//======================================================
template <typename InputIterator1, typename InputIterator2>
inline double Jaccard_W(InputIterator1 begin1, InputIterator1 end1, InputIterator2 begin2) {
  
  double bray = BrayCurtis_W(begin1, end1, begin2);
  
  // value to return
  double rval = 2 * bray / (1 + bray);
  return rval;
}


//======================================================
// Jaccard Unweighted
// 2 * BrayCurtis_U / (1 + BrayCurtis_U)
//======================================================
template <typename InputIterator1, typename InputIterator2>
inline double Jaccard_U(InputIterator1 begin1, InputIterator1 end1, InputIterator2 begin2) {
  
  double bray = BrayCurtis_U(begin1, end1, begin2);
  
  // value to return
  double rval = 2 * bray / (1 + bray);
  return rval;
}



//======================================================
// Entry point when new threads are spun up
//======================================================
struct BetaDivWorker : public Worker {
  
  const RMatrix<double> mat;
        RMatrix<double> rmat;
  const int             method;
  
  BetaDivWorker(const NumericMatrix mat, NumericMatrix rmat, int method)
    : mat(mat), rmat(rmat), method(method) {}
  
  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; i++) {
      for (std::size_t j = 0; j < i; j++) {
        
        // rows we will operate on
        RMatrix<double>::Row row1 = mat.row(i);
        RMatrix<double>::Row row2 = mat.row(j);
        
        // calculate distance and write to output matrix
        switch (method) {
          case BRAYCURTIS_W: { rmat(i,j) = BrayCurtis_W( row1.begin(), row1.end(), row2.begin() ); break; }
          case BRAYCURTIS_U: { rmat(i,j) = BrayCurtis_U( row1.begin(), row1.end(), row2.begin() ); break; }
          case EUCLIDEAN_W:  { rmat(i,j) = Euclidean_W(  row1.begin(), row1.end(), row2.begin() ); break; }
          case EUCLIDEAN_U:  { rmat(i,j) = Euclidean_U(  row1.begin(), row1.end(), row2.begin() ); break; }
          case MANHATTAN_W:  { rmat(i,j) = Manhattan_W(  row1.begin(), row1.end(), row2.begin() ); break; }
          case MANHATTAN_U:  { rmat(i,j) = Manhattan_U(  row1.begin(), row1.end(), row2.begin() ); break; }
          case JACCARD_W:    { rmat(i,j) = Jaccard_W(    row1.begin(), row1.end(), row2.begin() ); break; }
          case JACCARD_U:    { rmat(i,j) = Jaccard_U(    row1.begin(), row1.end(), row2.begin() ); break; }
        }
      }
    }
  }
};



//======================================================
// Interface to R
//======================================================
// [[Rcpp::export]]
NumericMatrix par_beta_div(NumericMatrix mat, const char* method, bool weighted) {
  
  
  // Map method and weighted to individual functions to call
  int metric = 0;
  if (( weighted) && (strcmp(method, "bray-curtis") == 0)) { metric = BRAYCURTIS_W; }
  if ((!weighted) && (strcmp(method, "bray-curtis") == 0)) { metric = BRAYCURTIS_U; }
  if (( weighted) && (strcmp(method, "euclidean")   == 0)) { metric = EUCLIDEAN_W;  }
  if ((!weighted) && (strcmp(method, "euclidean")   == 0)) { metric = EUCLIDEAN_U;  }
  if (( weighted) && (strcmp(method, "manhattan")   == 0)) { metric = MANHATTAN_W;  }
  if ((!weighted) && (strcmp(method, "manhattan")   == 0)) { metric = MANHATTAN_U;  }
  if (( weighted) && (strcmp(method, "jaccard")     == 0)) { metric = JACCARD_W;    }
  if ((!weighted) && (strcmp(method, "jaccard")     == 0)) { metric = JACCARD_U;    }
  
  
  // allocate the matrix we will return
  NumericMatrix rmat(mat.nrow(), mat.nrow());
  
  // create the worker and call it with parallelFor
  if (metric > 0) {
    BetaDivWorker betaDivWorker(mat, rmat, metric);
    parallelFor(0, mat.nrow(), betaDivWorker);
  }
  
  return rmat;
}

