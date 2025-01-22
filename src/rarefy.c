#include <R.h>
#include <Rinternals.h>
#include <string.h> // memset, memcpy
#include "get.h"

// Detect if pthread is available.
#if defined __has_include
#  if __has_include (<pthread.h>)
#    include <pthread.h>
#    include <stdlib.h> // calloc, free
#    define HAVE_PTHREAD
#  endif
#endif


typedef struct {
  int *samples;
  int *values;
  int  n_samples;
  int *indices;
  int  n_indices;
  int  target;
  int *rand_ints;
  int  n_threads;
  int  thread_i;
  int *result;
} rarefy_t;




//======================================================
// Limits itself to samples based on modulo.
//======================================================
void *rarefy_worker(void *arg) {
  
  rarefy_t setup = *((rarefy_t *) arg);
  
  int *samples   = setup.samples;
  int *values    = setup.values;
  int  n_samples = setup.n_samples;
  int *indices   = setup.indices;
  int  n_indices = setup.n_indices;
  int  target    = setup.target;
  int *rand_ints = setup.rand_ints;
  int  n_threads = setup.n_threads;
  int  thread_i  = setup.thread_i;
  int *result    = setup.result;
  
  //======================================================
  // `indices` orders the values by sample, then otu.
  // Reduces scanning; consistent seeded randomness.
  //======================================================
  
  int start, end = 0;
  
  for (int sample = thread_i; sample < n_samples; sample += n_threads) {
    
    //======================================================
    // Seek to the beginning of this sample's values
    //======================================================
    
    for (start = end; start < n_indices; start++) {
      int idx = indices[start] - 1;
      if (samples[idx] - 1 == sample) break;
    }
    
    if (start == n_indices) continue; // sample has no values
    
    
    
    //======================================================
    // Seek to end, summing depth along the way.
    //======================================================
    
    int depth = 0;
    
    for (end = start; end < n_indices; end++) {
      int idx = indices[end] - 1;
      if (samples[idx] - 1 != sample) break;
      depth += values[idx];
    }
    
    
    
    //======================================================
    // Insufficient sequences - leave as all zeroes.
    //======================================================
    
    if (depth < target) continue;
    
    
    //======================================================
    // Already rarefied - copy input to output
    //======================================================
    
    if (depth == target) {
      for (int i = start; i < end; i++) {
        int idx = indices[i] - 1;
        result[idx] = values[idx];
      }
      continue;
    }
    
    
    //======================================================
    // Knuth algorithm for choosing target seqs from depth.
    //======================================================
    
    int tried = 0, kept = 0;
    
    for (int i = start; i < end && kept < target; i++) {
      
      int idx  = indices[i] - 1;
      int seqs = values[idx];
      
      for (int seq = 0; seq < seqs && kept < target; seq++) {
        
        int not_tried  = depth - tried;
        int still_need = target - kept;
        
        if (rand_ints[tried] % not_tried < still_need) {
          result[idx]++; // retain this observation
          kept++;
        }
        
        tried++;
      }
    }
    
  }
  
  
  return NULL;
}



//======================================================
// R interface. Assigns samples to worker threads.
//======================================================
SEXP C_rarefy(
    SEXP sexp_otu_slam_mtx, SEXP sexp_indices, 
    SEXP sexp_target,       SEXP sexp_rand_ints, 
    SEXP sexp_n_threads ) {
  
  int *samples      = INTEGER(  get(sexp_otu_slam_mtx, "j"));
  int *values       = INTEGER(  get(sexp_otu_slam_mtx, "v"));
  int  n_samples    = asInteger(get(sexp_otu_slam_mtx, "ncol"));
  
  int *indices      = INTEGER(sexp_indices);
  int  n_indices    = length(sexp_indices);
  
  int  target       = asInteger(sexp_target);
  int *rand_ints    = INTEGER(sexp_rand_ints);
  int  n_threads    = asInteger(sexp_n_threads);
  
  
  // allocate the return vector
  SEXP sexp_result = PROTECT(allocVector(INTSXP, n_indices));
  int *result      = INTEGER(sexp_result);
  memset(result, 0, n_indices * sizeof(int));
  
  // common values for all the threads
  rarefy_t setup;
  setup.samples   = samples;
  setup.values    = values;
  setup.n_samples = n_samples;
  setup.indices   = indices;
  setup.n_indices = n_indices;
  setup.target    = target;
  setup.rand_ints = rand_ints;
  setup.n_threads = n_threads;
  setup.result    = result;
  
  
  // Run WITH multithreading
  #ifdef HAVE_PTHREAD
    if (n_threads > 1) {
      
      // threads and their individual input arguments
      pthread_t *tids = calloc(n_threads, sizeof(pthread_t));
      rarefy_t  *args = calloc(n_threads, sizeof(rarefy_t));
      
      if (tids == NULL || args == NULL) {
        free(tids); free(args);
        error("Unable to allocate memory for multithreaded rarefaction.");
        return R_NilValue;
      }
      
      for (int i = 0; i < n_threads; i++) {
        memcpy(args + i, &setup, sizeof(rarefy_t));
        args[i].thread_i = i;
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
  setup.thread_i  = 0;
  rarefy_worker(&setup);
  
  UNPROTECT(1);
  
  return sexp_result;
}

