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


typedef struct {
  double *otu_mtx;
  double *depths;
  int     n_otus;
  int     n_samples;
  int     n_threads;
  int     thread_i;
  double *result;
} adiv_t;



//======================================================
// Calculate Shannon Diversity Index
// Shannon <- -sum(p * log(p))
//======================================================

void *adiv_shannon(void *arg) {
  
  adiv_t setup = *((adiv_t *) arg);
  
  double *otu_mtx   = setup.otu_mtx;
  double *depths    = setup.depths;
  int     n_otus    = setup.n_otus;
  int     n_samples = setup.n_samples;
  int     n_threads = setup.n_threads;
  int     thread_i  = setup.thread_i;
  double *result    = setup.result;
  
  for (int sample = thread_i; sample < n_samples; sample += n_threads) {
    
    double *values = otu_mtx + (sample * n_otus);
    int     depth  = depths[sample];
    
    double sum = 0;
    for (int otu = 0; otu < n_otus; otu++) {
      double p = values[otu] / depth;
      if (p > 0) sum += p * log(p);
    }
    
    double shannon = -1 * sum;
    result[sample] = shannon;
  }
  
  return NULL;
}


//======================================================
// Calculate Chao1 Index
// v     <- ceiling(v)
// Chao1 <- n + (sum(v == 1) ** 2) / (2 * sum(v == 2))
//======================================================

void *adiv_chao1(void *arg) {
  
  adiv_t setup = *((adiv_t *) arg);
  
  double *otu_mtx   = setup.otu_mtx;
  int     n_otus    = setup.n_otus;
  int     n_samples = setup.n_samples;
  int     n_threads = setup.n_threads;
  int     thread_i  = setup.thread_i;
  double *result    = setup.result;
  
  for (int sample = thread_i; sample < n_samples; sample += n_threads) {
    
    double *values = otu_mtx + (sample * n_otus);
    
    int nz   = 0;
    int ones = 0;
    int twos = 0;
    for (int otu = 0; otu < n_otus; otu++) {
      
      double value = values[otu];
      if (value == 0) continue;
      
      nz++;
      if      (value <= 1) ones++;
      else if (value <= 2) twos++;
    }
    
    double chao1 = nz + (pow(ones, 2) / (2 * twos));
    result[sample] = chao1;
  }
  
  return NULL;
}


//======================================================
// Calculate Simpson Index
// Simpson <- 1 - sum(p ** 2)
//======================================================

void *adiv_simpson(void *arg) {
  
  adiv_t setup = *((adiv_t *) arg);
  
  double *otu_mtx   = setup.otu_mtx;
  double *depths    = setup.depths;
  int     n_otus    = setup.n_otus;
  int     n_samples = setup.n_samples;
  int     n_threads = setup.n_threads;
  int     thread_i  = setup.thread_i;
  double *result    = setup.result;
  
  for (int sample = thread_i; sample < n_samples; sample += n_threads) {
    
    double *values = otu_mtx + (sample * n_otus);
    int     depth  = depths[sample];
    
    double sum = 0;
    for (int otu = 0; otu < n_otus; otu++) {
      double p = values[otu] / depth;
      if (p > 0) sum += pow(p, 2);
    }
    
    double simpson = 1 - sum;
    result[sample] = simpson;
  }
  
  return NULL;
}


//======================================================
// Calculate Inverse Simpson Index
// InvSimpson <- 1 / sum(p ** 2)
//======================================================

void *adiv_inv_simpson(void *arg) {
  
  adiv_simpson(arg);
  
  
  adiv_t setup = *((adiv_t *) arg);
  
  int     n_samples = setup.n_samples;
  int     n_threads = setup.n_threads;
  int     thread_i  = setup.thread_i;
  double *result    = setup.result;
  
  for (int sample = thread_i; sample < n_samples; sample += n_threads) {
    
    double simpson     = result[sample];
    double inv_simpson = 1 / (1 - simpson);
    
    result[sample] = inv_simpson;
  }
  
  return NULL;
}


//======================================================
// R interface. Dispatches threads on compute methods.
//======================================================
SEXP C_alpha_div(
    SEXP sexp_otu_mtx,  SEXP sexp_depths, 
    SEXP sexp_algoritm, SEXP sexp_n_threads ) {
  
  double *otu_mtx   = REAL(sexp_otu_mtx);
  double *depths    = REAL(sexp_depths);
  int     algoritm  = asInteger(sexp_algoritm);
  int     n_threads = asInteger(sexp_n_threads);
  
  int     n_otus    = nrows(sexp_otu_mtx);
  int     n_samples = ncols(sexp_otu_mtx);
  
  
  // allocate the return vector
  SEXP sexp_result = PROTECT(allocVector(REALSXP, n_samples));
  double *result   = REAL(sexp_result);
  
  // function to run
  void * (*adiv_func)(void *) = NULL;
  switch (algoritm) {
  case 1:  adiv_func = adiv_shannon;     break;
  case 2:  adiv_func = adiv_chao1;       break;
  case 3:  adiv_func = adiv_simpson;     break;
  case 4:  adiv_func = adiv_inv_simpson; break;
  default: error("Invalid adiv metric."); return R_NilValue;
  }
  
  // common values for all the threads
  adiv_t setup;
  setup.otu_mtx   = otu_mtx;
  setup.depths    = depths;
  setup.n_otus    = n_otus;
  setup.n_samples = n_samples;
  setup.n_threads = n_threads;
  setup.result    = result;
  
  
  // Run WITH multithreading
  #ifdef HAVE_PTHREAD
    if (n_threads > 1) {
      
      // threads and their individual input arguments
      pthread_t *tids = calloc(n_threads, sizeof(pthread_t));
      adiv_t    *args = calloc(n_threads, sizeof(adiv_t));
      
      for (int i = 0; i < n_threads; i++) {
        memcpy(&args[i], &setup, sizeof(adiv_t));
        args[i].thread_i = i;
        pthread_create(&tids[i], NULL, adiv_func, &args[i]);
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
  adiv_func(&setup);
  
  UNPROTECT(1);
  
  return sexp_result;
}
