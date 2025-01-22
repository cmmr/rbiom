#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // calloc, free
#include <math.h>   // sqrt
#include "pair.h"

// Detect if pthread is available.
#if defined __has_include
#  if __has_include (<pthread.h>)
#    include <pthread.h>
#    include <string.h> // memcpy
#    define HAVE_PTHREAD
#  endif
#endif


typedef struct {
  dbl_ptr_pair_t *col_pairs;
  int             n_otus;
  int             n_pairs;
  int             n_threads;
  int             thread_i;
  double         *result;
} bdiv_t;



//======================================================
// Bray-Curtis Weighted
// sum(abs(x-y))/sum(x+y)
//======================================================

void *bdiv_braycurtis_w(void *arg) {
  
  bdiv_t setup = *((bdiv_t *) arg);
  
  dbl_ptr_pair_t *col_pairs = setup.col_pairs;
  int             n_otus    = setup.n_otus;
  int             n_pairs   = setup.n_pairs;
  int             n_threads = setup.n_threads;
  int             thread_i  = setup.thread_i;
  double         *result    = setup.result;
  
  for (int pair = thread_i; pair < n_pairs; pair += n_threads) {
    
    // pointers to each sample's column in otu_mtx
    double *sample_1_otus = col_pairs[pair].A;
    double *sample_2_otus = col_pairs[pair].B;
    
    double diffs = 0;
    double sums  = 0;
    
    for (int otu = 0; otu < n_otus; otu++) {
      
      // abundance of this OTU in the two samples
      double x = sample_1_otus[otu];
      double y = sample_2_otus[otu];
      
      // accumulate if appropriate
      sums += x + y;
      if (x > y) diffs += x - y;
      if (y > x) diffs += y - x;
    }
    
    // value to return
    result[pair] = diffs / sums;
  }
  
  return NULL;
}



//======================================================
// Bray-Curtis Unweighted
// (sum(x>0)+sum(y>0)-2*sum(x&y))/(sum(x>0)+sum(y>0))
//======================================================

void *bdiv_braycurtis_u(void *arg) {
  
  bdiv_t setup = *((bdiv_t *) arg);
  
  dbl_ptr_pair_t *col_pairs = setup.col_pairs;
  int             n_otus    = setup.n_otus;
  int             n_pairs   = setup.n_pairs;
  int             n_threads = setup.n_threads;
  int             thread_i  = setup.thread_i;
  double         *result    = setup.result;
  
  for (int pair = thread_i; pair < n_pairs; pair += n_threads) {
    
    // pointers to each sample's column in otu_mtx
    double *sample_1_otus = col_pairs[pair].A;
    double *sample_2_otus = col_pairs[pair].B;

    double  x_nz = 0;
    double  y_nz = 0;
    double xy_nz = 0;

    for (int otu = 0; otu < n_otus; otu++) {
      
      // abundance of this OTU in the two samples
      double x = sample_1_otus[otu];
      double y = sample_2_otus[otu];
      
      // accumulate if appropriate
      if (x > 0)           x_nz++;
      if (y > 0)           y_nz++;
      if (x > 0 && y > 0) xy_nz++;
    }

    // value to return
    result[pair] = (x_nz + y_nz - 2 * xy_nz) / (x_nz + y_nz);
  }
  
  return NULL;
}


//======================================================
// Euclidean Weighted
// sqrt(sum((x-y)^2))
//======================================================

void *bdiv_euclidean_w(void *arg) {
  
  bdiv_t setup = *((bdiv_t *) arg);
  
  dbl_ptr_pair_t *col_pairs = setup.col_pairs;
  int             n_otus    = setup.n_otus;
  int             n_pairs   = setup.n_pairs;
  int             n_threads = setup.n_threads;
  int             thread_i  = setup.thread_i;
  double         *result    = setup.result;
  
  for (int pair = thread_i; pair < n_pairs; pair += n_threads) {
    
    // pointers to each sample's column in otu_mtx
    double *sample_1_otus = col_pairs[pair].A;
    double *sample_2_otus = col_pairs[pair].B;

    double distance = 0;

    for (int otu = 0; otu < n_otus; otu++) {
      
      // abundance of this OTU in the two samples
      double x = sample_1_otus[otu];
      double y = sample_2_otus[otu];

      // accumulate if appropriate
      distance += (x - y) * (x - y);
    }

    // value to return
    result[pair] = sqrt(distance);
  }
  
  return NULL;
}


//======================================================
// Euclidean Unweighted
// sqrt(sum(x>0)+sum(y>0)-2*sum(x&y))
//======================================================

void *bdiv_euclidean_u(void *arg) {
  
  bdiv_t setup = *((bdiv_t *) arg);
  
  dbl_ptr_pair_t *col_pairs = setup.col_pairs;
  int             n_otus    = setup.n_otus;
  int             n_pairs   = setup.n_pairs;
  int             n_threads = setup.n_threads;
  int             thread_i  = setup.thread_i;
  double         *result    = setup.result;
  
  for (int pair = thread_i; pair < n_pairs; pair += n_threads) {
    
    // pointers to each sample's column in otu_mtx
    double *sample_1_otus = col_pairs[pair].A;
    double *sample_2_otus = col_pairs[pair].B;

    double  x_nz = 0;
    double  y_nz = 0;
    double xy_nz = 0;

    for (int otu = 0; otu < n_otus; otu++) {
      
      // abundance of this OTU in the two samples
      double x = sample_1_otus[otu];
      double y = sample_2_otus[otu];

      // accumulate if appropriate
      if (x > 0)           x_nz++;
      if (y > 0)           y_nz++;
      if (x > 0 && y > 0) xy_nz++;
    }

    // value to return
    result[pair] = sqrt(x_nz + y_nz - 2 * xy_nz);
  }
  
  return NULL;
}


//======================================================
// Manhattan Weighted
// sum(abs(x-y))
//======================================================

void *bdiv_manhattan_w(void *arg) {
  
  bdiv_t setup = *((bdiv_t *) arg);
  
  dbl_ptr_pair_t *col_pairs = setup.col_pairs;
  int             n_otus    = setup.n_otus;
  int             n_pairs   = setup.n_pairs;
  int             n_threads = setup.n_threads;
  int             thread_i  = setup.thread_i;
  double         *result    = setup.result;
  
  for (int pair = thread_i; pair < n_pairs; pair += n_threads) {
    
    // pointers to each sample's column in otu_mtx
    double *sample_1_otus = col_pairs[pair].A;
    double *sample_2_otus = col_pairs[pair].B;

    double distance = 0;

    for (int otu = 0; otu < n_otus; otu++) {
      
      // abundance of this OTU in the two samples
      double x = sample_1_otus[otu];
      double y = sample_2_otus[otu];

      // accumulate if appropriate
      if (x > y) distance += x - y;
      if (y > x) distance += y - x;
    }

    // value to return
    result[pair] = distance;
  }
  
  return NULL;
}


//======================================================
// Manhattan Unweighted
// sum(x>0)+sum(y>0)-2*sum(x&y)
//======================================================

void *bdiv_manhattan_u(void *arg) {
  
  bdiv_t setup = *((bdiv_t *) arg);
  
  dbl_ptr_pair_t *col_pairs = setup.col_pairs;
  int             n_otus    = setup.n_otus;
  int             n_pairs   = setup.n_pairs;
  int             n_threads = setup.n_threads;
  int             thread_i  = setup.thread_i;
  double         *result    = setup.result;
  
  for (int pair = thread_i; pair < n_pairs; pair += n_threads) {
    
    // pointers to each sample's column in otu_mtx
    double *sample_1_otus = col_pairs[pair].A;
    double *sample_2_otus = col_pairs[pair].B;

    double  x_nz = 0;
    double  y_nz = 0;
    double xy_nz = 0;

    for (int otu = 0; otu < n_otus; otu++) {
      
      // abundance of this OTU in the two samples
      double x = sample_1_otus[otu];
      double y = sample_2_otus[otu];

      // accumulate if appropriate
      if (x > 0)           x_nz++;
      if (y > 0)           y_nz++;
      if (x > 0 && y > 0) xy_nz++;
    }

    // value to return
    result[pair] = x_nz + y_nz - 2 * xy_nz;
  }
  
  return NULL;
}


//======================================================
// Jaccard Weighted
// 2 * BrayCurtis_W / (1 + BrayCurtis_W)
//======================================================

void *bdiv_jaccard_w(void *arg) {

  bdiv_braycurtis_w(arg);
  
  bdiv_t setup = *((bdiv_t *) arg);
  
  int     n_pairs   = setup.n_pairs;
  int     n_threads = setup.n_threads;
  int     thread_i  = setup.thread_i;
  double *result    = setup.result;
  
  for (int pair = thread_i; pair < n_pairs; pair += n_threads) {
    
    double distance = result[pair];
    
    result[pair] = 2 * distance / (1 + distance);
  }
  
  return NULL;
}


//======================================================
// Jaccard Unweighted
// 2 * BrayCurtis_U / (1 + BrayCurtis_U)
//======================================================

void *bdiv_jaccard_u(void *arg) {
  
  bdiv_braycurtis_u(arg);
  
  bdiv_t setup = *((bdiv_t *) arg);
  
  int     n_pairs   = setup.n_pairs;
  int     n_threads = setup.n_threads;
  int     thread_i  = setup.thread_i;
  double *result    = setup.result;
  
  for (int pair = thread_i; pair < n_pairs; pair += n_threads) {
    
    double distance = result[pair];
    
    result[pair] = 2 * distance / (1 + distance);
  }
  
  return NULL;
}



//======================================================
// R interface. Dispatches threads to bdiv algorithms.
//======================================================
SEXP C_beta_div(
    SEXP sexp_otu_mtx,  SEXP sexp_pair_mtx,  
    SEXP sexp_algoritm, SEXP sexp_n_threads ) {
  
  double     *otu_mtx   = REAL(sexp_otu_mtx);
  int_pair_t *pairs     = (int_pair_t *)INTEGER(sexp_pair_mtx);
  int         algoritm  = asInteger(sexp_algoritm);
  int         n_threads = asInteger(sexp_n_threads);
  
  int         n_otus    = nrows(sexp_otu_mtx);
  int         n_pairs   = ncols(sexp_pair_mtx);
  
  
  // allocate the return vector
  SEXP sexp_result = PROTECT(allocVector(REALSXP, n_pairs));
  double *result   = REAL(sexp_result);
  
  // function to run
  void * (*bdiv_func)(void *) = NULL;
  switch (algoritm) {
    case 1:  bdiv_func = bdiv_braycurtis_w; break;
    case 2:  bdiv_func = bdiv_braycurtis_u; break;
    case 3:  bdiv_func = bdiv_euclidean_w;  break;
    case 4:  bdiv_func = bdiv_euclidean_u;  break;
    case 5:  bdiv_func = bdiv_manhattan_w;  break;
    case 6:  bdiv_func = bdiv_manhattan_u;  break;
    case 7:  bdiv_func = bdiv_jaccard_w;    break;
    case 8:  bdiv_func = bdiv_jaccard_u;    break;
    default: error("Invalid bdiv metric."); return R_NilValue;
  }
  
  // pre-compute pointers to each sample's column in otu_mtx
  dbl_ptr_pair_t *col_pairs = calloc(n_pairs, sizeof(dbl_ptr_pair_t));
  for (int pair = 0; pair < n_pairs; pair++) {
    col_pairs[pair].A = &otu_mtx[(pairs[pair].A - 1) * n_otus];
    col_pairs[pair].B = &otu_mtx[(pairs[pair].B - 1) * n_otus];
  }
  
  
  // common values for all the threads
  bdiv_t setup;
  setup.col_pairs = col_pairs;
  setup.n_otus    = n_otus;
  setup.n_pairs   = n_pairs;
  setup.n_threads = n_threads;
  setup.result    = result;
  
  
  // Run WITH multithreading
  #ifdef HAVE_PTHREAD
    if (n_threads > 1) {
      
      // threads and their individual input arguments
      pthread_t *tids = calloc(n_threads, sizeof(pthread_t));
      bdiv_t    *args = calloc(n_threads, sizeof(bdiv_t));
      
      for (int i = 0; i < n_threads; i++) {
        memcpy(&args[i], &setup, sizeof(bdiv_t));
        args[i].thread_i = i;
        pthread_create(&tids[i], NULL, bdiv_func, &args[i]);
      }
      
      for (int i = 0; i < n_threads; i++) {
        pthread_join(tids[i], NULL);
      }
      
      free(tids);
      free(args);
      free(col_pairs);
      
      UNPROTECT(1);
      
      return sexp_result;
    }
  #endif
  
  
  // Run WITHOUT multithreading
  setup.n_threads = 1;
  setup.thread_i  = 0;
  bdiv_func(&setup);
  
  free(col_pairs);
  
  UNPROTECT(1);
  
  return sexp_result;
}

