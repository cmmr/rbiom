#include <R.h>
#include <Rinternals.h>
#include <string.h> // memset, memcpy

// Detect if pthread is available.
#if defined __has_include
#  if __has_include (<pthread.h>)
#    include <pthread.h>
#    include <stdlib.h> // calloc, free
#    define HAVE_PTHREAD
#  endif
#endif


typedef struct {
  int *values;
  int *n_otus;
  int *depths;
  int  target;
  int  n_samples;
  int *rand_ints;
  int  n_threads;
  int *result;
  int  remainder;
} rarefy_t;




//======================================================
// Limits itself to samples based on modulo.
//======================================================
void *rarefy_worker(void *arg) {
  
  rarefy_t setup = *((rarefy_t *) arg);
  
  int *values    = setup.values;
  int *n_otus    = setup.n_otus;
  int *depths    = setup.depths;
  int  target    = setup.target;
  int  n_samples = setup.n_samples;
  int *rand_ints = setup.rand_ints;
  int  n_threads = setup.n_threads;
  int *result    = setup.result;
  int  remainder = setup.remainder;
  
  //======================================================
  // Loop over each sample
  //======================================================
  
  int start, end = 0;
  
  for (int sample_i = 0; sample_i < n_samples; sample_i++) {
    
    int n     = n_otus[sample_i];
    int depth = depths[sample_i];
    
    start = end;
    end   = start + n;
    
    if (sample_i % n_threads != remainder) continue;
    
    
    //======================================================
    // Insufficient sequences - leave as all zeroes.
    //======================================================
    
    if (depth < target) {
      continue;
    }
    
    
    //======================================================
    // Already rarefied - copy input to output
    //======================================================
    
    if (depth == target) {
      memcpy(result + start, values + start, n * sizeof(int));
      continue;
    }
    
    
    //======================================================
    // Knuth algorithm for removing (target - depth) seqs.
    //======================================================
    
    int value = values[start], pos = start, kept = 0;
    
    for (int tried = 0; tried < depth && kept < target; ++tried) {
      
      int not_tried  = depth - tried;
      int still_need = target - kept;
      
      if (rand_ints[tried] % not_tried < still_need) {
        result[pos]++; // retain this observation
        kept++;
      }
      
      if (--value == 0) {
        pos++; // move to next OTU
        value = values[pos];
      }
    }
    
  }
  
  
  return NULL;
}



//======================================================
// R interface. Assigns samples to worker threads.
//======================================================
SEXP C_rarefy(
    SEXP sexp_values,    SEXP sexp_n_otus, 
    SEXP sexp_depths,    SEXP sexp_target, 
    SEXP sexp_rand_ints, SEXP sexp_n_threads) {
  
  int n_threads = asInteger(sexp_n_threads);
  int n_values  = length(sexp_values);
  
  
  // allocate the return vector
  SEXP sexp_result = PROTECT(allocVector(INTSXP, n_values));
  int *result      = INTEGER(sexp_result);
  memset(result, 0, n_values * sizeof(int));
  
  // common values for all the threads
  rarefy_t setup;
  setup.values    = INTEGER(sexp_values);
  setup.n_otus    = INTEGER(sexp_n_otus);
  setup.depths    = INTEGER(sexp_depths);
  setup.target    = asInteger(sexp_target);
  setup.rand_ints = INTEGER(sexp_rand_ints);
  setup.n_samples = length(sexp_n_otus);
  setup.n_threads = n_threads;
  setup.result    = result;
  
  
  // Run WITH multithreading
  #ifdef HAVE_PTHREAD
    if (n_threads > 1) {
      
      // threads and their individual input arguments
      pthread_t *tids = calloc(n_threads, sizeof(pthread_t));
      rarefy_t  *args = calloc(n_threads, sizeof(rarefy_t));
      
      for (int i = 0; i < n_threads; i++) {
        memcpy(args + i, &setup, sizeof(rarefy_t));
        args[i].remainder = i;
        pthread_create(&tids[i], NULL, rarefy_worker, &args[i]);
      }
      
      for (int i = 0; i < n_threads; i++) {
        pthread_join(tids[i], NULL);
      }
      
      free(tids);
      free(args);
      
      UNPROTECT(1);
      
      return sexp_result;
    }
  #endif
  
  
  // Run WITHOUT multithreading
  setup.n_threads = 1;
  setup.remainder = 0;
  rarefy_worker(&setup);
  
  UNPROTECT(1);
  
  return sexp_result;
}

