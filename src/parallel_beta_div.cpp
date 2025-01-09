#include <Rcpp.h>
#include <RcppThread.h>
using namespace Rcpp;


//======================================================
// Bray-Curtis Weighted
// sum(abs(x-y))/sum(x+y)
//======================================================

// [[Rcpp::export]]
NumericVector bdiv_braycurtis_w(
    NumericVector counts, int n_otus,
    IntegerVector pairs,  int n_pairs,
    int nThreads ) {

  if (nThreads < 0)
    nThreads = std::thread::hardware_concurrency();

  // allocate the vector we will return
  NumericVector distances(n_pairs);

  RcppThread::parallelFor(0, n_pairs, [&] (int pair) {

    // find the offsets in the counts matrix
    int sample_1 = pairs[(pair * 2) + 0] * n_otus;
    int sample_2 = pairs[(pair * 2) + 1] * n_otus;

    double diffs = 0;
    double sums  = 0;

    for (int otu = 0; otu < n_otus; otu++) {

      // abundance of this OTU in the two samples
      double x = counts[sample_1 + otu];
      double y = counts[sample_2 + otu];

      // accumulate if appropriate
      sums += x + y;
      if (x > y) diffs += x - y;
      if (y > x) diffs += y - x;
    }

    // value to return
    distances[pair] = diffs / sums;

  }, nThreads);

  return distances;
}



//======================================================
// Bray-Curtis Unweighted
// (sum(x>0)+sum(y>0)-2*sum(x&y))/(sum(x>0)+sum(y>0))
//======================================================

// [[Rcpp::export]]
NumericVector bdiv_braycurtis_u(
    NumericVector counts, int n_otus,
    IntegerVector pairs,  int n_pairs,
    int nThreads ) {

  if (nThreads < 0)
    nThreads = std::thread::hardware_concurrency();

  // allocate the vector we will return
  NumericVector distances(n_pairs);

  RcppThread::parallelFor(0, n_pairs, [&] (int pair) {

    // find the offsets in the counts matrix
    int sample_1 = pairs[(pair * 2) + 0] * n_otus;
    int sample_2 = pairs[(pair * 2) + 1] * n_otus;

    double  x_nz = 0;
    double  y_nz = 0;
    double xy_nz = 0;

    for (int otu = 0; otu < n_otus; otu++) {

      // abundance of this OTU in the two samples
      double x = counts[sample_1 + otu];
      double y = counts[sample_2 + otu];

      // accumulate if appropriate
      if (x > 0)           x_nz++;
      if (y > 0)           y_nz++;
      if (x > 0 && y > 0) xy_nz++;
    }

    // value to return
    distances[pair] = (x_nz + y_nz - 2 * xy_nz) / (x_nz + y_nz);

  }, nThreads);

  return distances;
}


//======================================================
// Euclidean Weighted
// sqrt(sum((x-y)^2))
//======================================================

// [[Rcpp::export]]
NumericVector bdiv_euclidean_w(
    NumericVector counts, int n_otus,
    IntegerVector pairs,  int n_pairs,
    int nThreads ) {

  if (nThreads < 0)
    nThreads = std::thread::hardware_concurrency();

  // allocate the vector we will return
  NumericVector distances(n_pairs);

  RcppThread::parallelFor(0, n_pairs, [&] (int pair) {

    // find the offsets in the counts matrix
    int sample_1 = pairs[(pair * 2) + 0] * n_otus;
    int sample_2 = pairs[(pair * 2) + 1] * n_otus;

    double result = 0;

    for (int otu = 0; otu < n_otus; otu++) {

      // abundance of this OTU in the two samples
      double x = counts[sample_1 + otu];
      double y = counts[sample_2 + otu];

      // accumulate if appropriate
      result += (x - y) * (x - y);
    }

    // value to return
    distances[pair] = sqrt(result);

  }, nThreads);

  return distances;
}


//======================================================
// Euclidean Unweighted
// sqrt(sum(x>0)+sum(y>0)-2*sum(x&y))
//======================================================

// [[Rcpp::export]]
NumericVector bdiv_euclidean_u(
    NumericVector counts, int n_otus,
    IntegerVector pairs,  int n_pairs,
    int nThreads ) {

  if (nThreads < 0)
    nThreads = std::thread::hardware_concurrency();

  // allocate the vector we will return
  NumericVector distances(n_pairs);

  RcppThread::parallelFor(0, n_pairs, [&] (int pair) {

    // find the offsets in the counts matrix
    int sample_1 = pairs[(pair * 2) + 0] * n_otus;
    int sample_2 = pairs[(pair * 2) + 1] * n_otus;

    double  x_nz = 0;
    double  y_nz = 0;
    double xy_nz = 0;

    for (int otu = 0; otu < n_otus; otu++) {

      // abundance of this OTU in the two samples
      double x = counts[sample_1 + otu];
      double y = counts[sample_2 + otu];

      // accumulate if appropriate
      if (x > 0)           x_nz++;
      if (y > 0)           y_nz++;
      if (x > 0 && y > 0) xy_nz++;
    }

    // value to return
    distances[pair] = sqrt(x_nz + y_nz - 2 * xy_nz);

  }, nThreads);

  return distances;
}


//======================================================
// Manhattan Weighted
// sum(abs(x-y))
//======================================================

// [[Rcpp::export]]
NumericVector bdiv_manhattan_w(
    NumericVector counts, int n_otus,
    IntegerVector pairs,  int n_pairs,
    int nThreads ) {

  if (nThreads < 0)
    nThreads = std::thread::hardware_concurrency();

  // allocate the vector we will return
  NumericVector distances(n_pairs);

  RcppThread::parallelFor(0, n_pairs, [&] (int pair) {

    // find the offsets in the counts matrix
    int sample_1 = pairs[(pair * 2) + 0] * n_otus;
    int sample_2 = pairs[(pair * 2) + 1] * n_otus;

    double result = 0;

    for (int otu = 0; otu < n_otus; otu++) {

      // abundance of this OTU in the two samples
      double x = counts[sample_1 + otu];
      double y = counts[sample_2 + otu];

      // accumulate if appropriate
      if (x > y) result += x - y;
      if (y > x) result += y - x;
    }

    // value to return
    distances[pair] = result;

  }, nThreads);

  return distances;
}


//======================================================
// Manhattan Unweighted
// sum(x>0)+sum(y>0)-2*sum(x&y)
//======================================================

// [[Rcpp::export]]
NumericVector bdiv_manhattan_u(
    NumericVector counts, int n_otus,
    IntegerVector pairs,  int n_pairs,
    int nThreads ) {

  if (nThreads < 0)
    nThreads = std::thread::hardware_concurrency();

  // allocate the vector we will return
  NumericVector distances(n_pairs);

  RcppThread::parallelFor(0, n_pairs, [&] (int pair) {

    // find the offsets in the counts matrix
    int sample_1 = pairs[(pair * 2) + 0] * n_otus;
    int sample_2 = pairs[(pair * 2) + 1] * n_otus;

    double  x_nz = 0;
    double  y_nz = 0;
    double xy_nz = 0;

    for (int otu = 0; otu < n_otus; otu++) {

      // abundance of this OTU in the two samples
      double x = counts[sample_1 + otu];
      double y = counts[sample_2 + otu];

      // accumulate if appropriate
      if (x > 0)           x_nz++;
      if (y > 0)           y_nz++;
      if (x > 0 && y > 0) xy_nz++;
    }

    // value to return
    distances[pair] = x_nz + y_nz - 2 * xy_nz;

  }, nThreads);

  return distances;
}


//======================================================
// Jaccard Weighted
// 2 * BrayCurtis_W / (1 + BrayCurtis_W)
//======================================================

// [[Rcpp::export]]
NumericVector bdiv_jaccard_w(
    NumericVector counts, int n_otus,
    IntegerVector pairs,  int n_pairs,
    int nThreads ) {

  NumericVector distances = bdiv_braycurtis_w(counts, n_otus, pairs, n_pairs, nThreads);
  distances = 2 * distances / (1 + distances);

  return distances;
}


//======================================================
// Jaccard Unweighted
// 2 * BrayCurtis_U / (1 + BrayCurtis_U)
//======================================================

// [[Rcpp::export]]
NumericVector bdiv_jaccard_u(
    NumericVector counts, int n_otus,
    IntegerVector pairs,  int n_pairs,
    int nThreads ) {

  NumericVector distances = bdiv_braycurtis_u(counts, n_otus, pairs, n_pairs, nThreads);
  distances = 2 * distances / (1 + distances);

  return distances;
}



/*

 =================================
 -- OLD CODE USING RcppParallel --
 =================================

 --> Converted to RcppThreads on 1/8/2025 to avoid UBSAN errors on CRAN checks.


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

    const RMatrix<double> counts;
    const RMatrix<int>    pairs;
    const int             algorithm;
    RVector<double> distances;

    BetaDivWorker(const NumericMatrix counts, const IntegerMatrix pairs, int algorithm, NumericVector distances)
      : counts(counts), pairs(pairs), algorithm(algorithm), distances(distances) {}

    void operator()(std::size_t begin, std::size_t end) {

      // Rcpp::checkUserInterrupt();

      for (std::size_t i = begin; i < end; i++) {

        // rows we will operate on
        RMatrix<double>::Row row1 = counts.row(pairs(i,0));
        RMatrix<double>::Row row2 = counts.row(pairs(i,1));

        // calculate distance and write to output matrix
        switch (algorithm) {
        case BRAYCURTIS_W: { distances[i] = BrayCurtis_W( row1.begin(), row1.end(), row2.begin() ); break; }
        case BRAYCURTIS_U: { distances[i] = BrayCurtis_U( row1.begin(), row1.end(), row2.begin() ); break; }
        case EUCLIDEAN_W:  { distances[i] = Euclidean_W(  row1.begin(), row1.end(), row2.begin() ); break; }
        case EUCLIDEAN_U:  { distances[i] = Euclidean_U(  row1.begin(), row1.end(), row2.begin() ); break; }
        case MANHATTAN_W:  { distances[i] = Manhattan_W(  row1.begin(), row1.end(), row2.begin() ); break; }
        case MANHATTAN_U:  { distances[i] = Manhattan_U(  row1.begin(), row1.end(), row2.begin() ); break; }
        case JACCARD_W:    { distances[i] = Jaccard_W(    row1.begin(), row1.end(), row2.begin() ); break; }
        case JACCARD_U:    { distances[i] = Jaccard_U(    row1.begin(), row1.end(), row2.begin() ); break; }
        }

      }
    }
  };



  //======================================================
  // Interface to R
  //======================================================
  // [[Rcpp::export]]
  NumericVector par_beta_div(NumericMatrix counts, IntegerMatrix pairs, const char* bdiv, bool weighted) {


    // Map bdiv method and weighted to individual functions to call
    int algorithm = 0;
    if (( weighted) && (strcmp(bdiv, "bray-curtis") == 0)) { algorithm = BRAYCURTIS_W; }
    if ((!weighted) && (strcmp(bdiv, "bray-curtis") == 0)) { algorithm = BRAYCURTIS_U; }
    if (( weighted) && (strcmp(bdiv, "euclidean")   == 0)) { algorithm = EUCLIDEAN_W;  }
    if ((!weighted) && (strcmp(bdiv, "euclidean")   == 0)) { algorithm = EUCLIDEAN_U;  }
    if (( weighted) && (strcmp(bdiv, "manhattan")   == 0)) { algorithm = MANHATTAN_W;  }
    if ((!weighted) && (strcmp(bdiv, "manhattan")   == 0)) { algorithm = MANHATTAN_U;  }
    if (( weighted) && (strcmp(bdiv, "jaccard")     == 0)) { algorithm = JACCARD_W;    }
    if ((!weighted) && (strcmp(bdiv, "jaccard")     == 0)) { algorithm = JACCARD_U;    }


    // allocate the vector we will return
    NumericVector distances(pairs.nrow());

    // create the worker and call it with parallelFor
    if (algorithm > 0) {
      BetaDivWorker betaDivWorker(counts, pairs, algorithm, distances);
      parallelFor(0, pairs.nrow(), betaDivWorker);
    }

    return distances;
  }

*/

