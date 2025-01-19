#include <R.h>
#include <Rinternals.h>
#include <math.h> // log, pow

// Detect if pthread is available.
#if defined __has_include
#  if __has_include (<pthread.h>)
#    include <pthread.h>
#    include <stdlib.h> // calloc, free
#    include <string.h> // memcpy
#    define HAVE_PTHREAD
#  endif
#endif

#define SHANNON     1
#define CHAO1       2
#define SIMPSONS    4


typedef struct {
  double *otu_mtx;
  int     algorithms;
  int     n_otus;
  int     n_samples;
  int     n_threads;
  int     thread_i;
  double *result;
} adiv_t;


void *adiv_worker(void *arg) {
  
  adiv_t setup = *((adiv_t *) arg);
  
  double *otu_mtx    = setup.otu_mtx;
  int     algorithms = setup.algorithms;
  int     n_otus     = setup.n_otus;
  int     n_samples  = setup.n_samples;
  int     n_threads  = setup.n_threads;
  int     thread_i   = setup.thread_i;
  double *result     = setup.result;
  
  double *res_depth       = result + n_samples * 0;
  double *res_nnz         = result + n_samples * 1;
  double *res_shannon     = result + n_samples * 2;
  double *res_chao1       = result + n_samples * 3;
  double *res_simpson     = result + n_samples * 4;
  double *res_inv_simpson = result + n_samples * 5;
  
  int compute_shannon   = algorithms & SHANNON;
  int compute_chao1     = algorithms & CHAO1;
  int compute_non_chao1 = algorithms ^ CHAO1;
  int compute_simpsons  = algorithms & SIMPSONS;
  
  
  
  //======================================================
  // Calculate measures of alpha diversity
  // depth <- sum(v)
  // nnz   <- sum(v > 0)
  // Chao1 <- n + (sum(v == 1) ** 2) / (2 * sum(v == 2))
  // Shannon    <- -sum(p * log(p))
  // Simpson    <- 1 - sum(p ** 2)
  // InvSimpson <- 1 / sum(p ** 2)
  //======================================================
  
  for (int sample = thread_i; sample < n_samples; sample += n_threads) {
    
    double *values = otu_mtx + (sample * n_otus);
    
    double depth      = 0;
    double nnz        = 0;
    double chao1_ones = 0;
    double chao1_twos = 0;
    
    for (int otu = 0; otu < n_otus; otu++) {
      
      double value = values[otu];
      if (value == 0) continue;
      
      nnz++;
      depth += value;
      
      if (compute_chao1) {
        if      (value <= 1) chao1_ones++;
        else if (value <= 2) chao1_twos++;
      }
    }
    
    res_depth[sample] = depth;
    res_nnz[sample]   = nnz;
    
    if (compute_chao1) {
      res_chao1[sample] = nnz + (pow(chao1_ones, 2) / (2 * chao1_twos));
    }
    
    
    if (compute_non_chao1) {
    
      double shannon_sum = 0;
      double simpson_sum = 0;
      
      for (int otu = 0; otu < n_otus; otu++) {
        
        double value = values[otu];
        if (value == 0) continue;
        
        double p = values[otu] / depth;
        
        if (compute_shannon)  shannon_sum += p * log(p);
        if (compute_simpsons) simpson_sum += pow(p, 2);
      }
      
      if (compute_shannon) {
        res_shannon[sample] = -1 * shannon_sum;
      }
      
      if (compute_simpsons) {
        res_simpson[sample]     = 1 - simpson_sum;
        res_inv_simpson[sample] = 1 / simpson_sum;
      }
    }
    
  }
  
  
  
  return NULL;
}



//======================================================
// R interface. Dispatches threads on compute methods.
//======================================================
SEXP C_alpha_div(SEXP sexp_otu_mtx, SEXP sexp_algoritms, SEXP sexp_n_threads) {
  
  double *otu_mtx    = REAL(sexp_otu_mtx);
  int     n_otus     = nrows(sexp_otu_mtx);
  int     n_samples  = ncols(sexp_otu_mtx);
  int     algorithms = asInteger(sexp_algoritms);
  int     n_threads  = asInteger(sexp_n_threads);
  
  // allocate the return vector
  SEXP sexp_result = PROTECT(allocVector(REALSXP, n_samples * 6));
  double *result   = REAL(sexp_result);
  
  // common values for all the threads
  adiv_t setup;
  setup.otu_mtx    = otu_mtx;
  setup.algorithms = algorithms;
  setup.n_otus     = n_otus;
  setup.n_samples  = n_samples;
  setup.n_threads  = n_threads;
  setup.result     = result;
  
  
  // Run WITH multithreading
  #ifdef HAVE_PTHREAD
    if (n_threads > 1) {
      
      // threads and their individual input arguments
      pthread_t *tids = calloc(n_threads, sizeof(pthread_t));
      adiv_t    *args = calloc(n_threads, sizeof(adiv_t));
      
      for (int i = 0; i < n_threads; i++) {
        memcpy(&args[i], &setup, sizeof(adiv_t));
        args[i].thread_i = i;
        pthread_create(&tids[i], NULL, adiv_worker, &args[i]);
      }
      
      for (int i = 0; i < n_threads; i++)
        pthread_join(tids[i], NULL);
      
      free(tids);
      free(args);
      
      UNPROTECT(1);
      
      return sexp_result;
    }
  #endif
  
  
  // Run WITHOUT multithreading
  setup.n_threads = 1;
  setup.thread_i  = 0;
  adiv_worker(&setup);
  
  UNPROTECT(1);
  
  return sexp_result;
}
