#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // calloc, free
#include <string.h> // strcmp, memset, memcpy

// Detect if pthread is available.
#if defined __has_include
#  if __has_include (<pthread.h>)
#    include <pthread.h>
#    define HAVE_PTHREAD
#  endif
#endif


typedef struct {
  int    *otus;
  int    *samples;
  double *values;
  int     n_otus;
  int     n_samples;
  int    *pair_mtx;
  int    *edge_mtx;
  double *edge_lengths;
  int     weighted;
  int     n_values;
  int     n_pairs;
  int     n_edges;
  double *sample_depths;
  int    *child_at;
  int    *edge_to_leaves;
  double *sample_edge_wt;
  int     n_threads;
  int     thread_i;
  double *result;
} unifrac_t;


//======================================================
// Trace the path from each leaf to the root node.
//======================================================

void *calc_edge_to_leaves(void *arg) {
  
  unifrac_t setup = *((unifrac_t *) arg);
  
  int  n_otus         = setup.n_otus;
  int  n_edges        = setup.n_edges;
  int *child_at       = setup.child_at;
  int *edge_mtx       = setup.edge_mtx;
  int *edge_to_leaves = setup.edge_to_leaves;
  int  n_threads      = setup.n_threads;
  int  thread_i       = setup.thread_i;
  
  // IntegerMatrix edge2leaves = IntegerMatrix(nTreeEdges, nOTUs);
  // RcppThread::parallelFor(1, nOTUs + 1, [&] (int leaf) {
  //   int node = leaf;
  //   do {
  //     edge2leaves(childAt[node], leaf - 1) = 1;
  //     node = treeEdge(childAt[node], 0);
  //   } while (childAt[node] != -1);
  // }, nThreads);
  
  for (int leaf = thread_i; leaf < n_otus; leaf += n_threads) {
    
    int node = leaf;
    int edge = child_at[node];
    do {
      edge_to_leaves[edge * n_otus + leaf] = 1;
      node = edge_mtx[0 * n_edges + edge] - 1;
      edge = child_at[node];
    } while (edge != -1);
    
  }
  
  return NULL;
}


//======================================================
// Map samples to their branch weights.
//======================================================

void *calc_sample_edge_wt(void *arg) {
  
  unifrac_t setup = *((unifrac_t *) arg);
  
  int    *otus           = setup.otus;
  int    *samples        = setup.samples;
  double *values         = setup.values;
  double *sample_depths  = setup.sample_depths;
  double *edge_lengths   = setup.edge_lengths;
  int    *edge_to_leaves = setup.edge_to_leaves;
  double *sample_edge_wt = setup.sample_edge_wt;
  int     n_edges        = setup.n_edges;
  int     n_otus         = setup.n_otus;
  int     n_samples      = setup.n_samples;
  int     n_values       = setup.n_values;
  int     weighted       = setup.weighted;
  int     n_threads      = setup.n_threads;
  int     thread_i       = setup.thread_i;
  
  if (weighted) {
  
    for (int i = thread_i; i < n_values; i += n_threads) {
      
      int    otu          = otus[i]    - 1;
      int    sample       = samples[i] - 1;
      double value        = values[i];
      double sample_depth = sample_depths[sample];
      
      double wt = 0;
      for (int edge = 0; edge < n_edges; edge++) {
        if (edge_to_leaves[edge * n_otus + otu] == 0) continue;
        wt = edge_lengths[edge] * (value / sample_depth);
        sample_edge_wt[edge * n_samples + sample] += wt;
      }
    }
    
  } else {
    
    for (int i = thread_i; i < n_values; i += n_threads) {
      
      int otu    = otus[i]    - 1;
      int sample = samples[i] - 1;
    
      for (int edge = 0; edge < n_edges; edge++) {
        if (edge_to_leaves[edge * n_otus + otu] == 0) continue;
        sample_edge_wt[edge * n_samples + sample] = edge_lengths[edge];
      }
    }
    
  }
  
  return NULL;
}


//======================================================
// Compute pairwise distances.
//======================================================

void *calc_result(void *arg) {
  
  unifrac_t setup = *((unifrac_t *) arg);
  
  int     weighted       = setup.weighted;
  int     thread_i       = setup.thread_i;
  int     n_threads      = setup.n_threads;
  int     n_pairs        = setup.n_pairs;
  int    *pair_mtx       = setup.pair_mtx;
  int     n_edges        = setup.n_edges;
  double *sample_edge_wt = setup.sample_edge_wt;
  int     n_samples      = setup.n_samples;
  double *result         = setup.result;
  
  if (weighted) {
    
    for (int pair = thread_i; pair < n_pairs; pair += n_threads) {
      
      int sample1 = pair_mtx[0 * n_pairs + pair] - 1;
      int sample2 = pair_mtx[1 * n_pairs + pair] - 1;
      
      double amt_diff = 0;
      for (int edge = 0; edge < n_edges; edge++) {
        double x = sample_edge_wt[edge * n_samples + sample1];
        double y = sample_edge_wt[edge * n_samples + sample2];
        amt_diff += (x >= y) ? (x - y) : (y - x);
      }
      result[pair] = amt_diff;
    }
    
  } else {
    
    for (int pair = thread_i; pair < n_pairs; pair += n_threads) {
      
      int sample1 = pair_mtx[0 * n_pairs + pair] - 1;
      int sample2 = pair_mtx[1 * n_pairs + pair] - 1;
      
      double amt_diff = 0, amt_total = 0;
      for (int edge = 0; edge < n_edges; edge++) {
        double x = sample_edge_wt[edge * n_samples + sample1];
        double y = sample_edge_wt[edge * n_samples + sample2];
        amt_diff  += (x && y) ? (0) : (x + y);
        amt_total += (x && y) ? (x) : (x + y);
      }
      result[pair] = amt_diff / amt_total;
    }
    
  }
    
  
  return NULL;
}



//======================================================
// Extract list elements by name.
//======================================================
SEXP get(SEXP sexp_list, const char *str) {
  
  SEXP names = getAttrib(sexp_list, R_NamesSymbol);
  
  for (int i = 0; i < length(sexp_list); i++) {
    if (strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
      return VECTOR_ELT(sexp_list, i);
    }
  }
    
  return R_NilValue;
}


int any_null (void *a, void *b, void *c, void *d) {
  return(a == NULL || b == NULL || c == NULL || d == NULL);
}
void *free_4 (void *a, void *b, void *c, void *d) {
  free(a); free(b); free(c); free(d); return NULL;
}




//======================================================
// R interface. Dispatches threads on compute methods.
//======================================================
SEXP C_unifrac(
    SEXP sexp_otu_slam_mtx, SEXP sexp_phylo_tree,
    SEXP sexp_pair_mtx, SEXP sexp_weighted, SEXP sexp_n_threads) {
  
  int    *otus         = INTEGER(  get(sexp_otu_slam_mtx, "i"));
  int    *samples      = INTEGER(  get(sexp_otu_slam_mtx, "j"));
  double *values       = REAL(     get(sexp_otu_slam_mtx, "v"));
  int     n_values     = length(   get(sexp_otu_slam_mtx, "v"));
  int     n_otus       = asInteger(get(sexp_otu_slam_mtx, "nrow"));
  int     n_samples    = asInteger(get(sexp_otu_slam_mtx, "ncol"));
  
  int    *edge_mtx     = INTEGER(  get(sexp_phylo_tree, "edge"));
  int     n_edges      = nrows(    get(sexp_phylo_tree, "edge"));
  double *edge_lengths = REAL(     get(sexp_phylo_tree, "edge.length"));
  
  int    *pair_mtx     = INTEGER(  sexp_pair_mtx);
  int     n_pairs      = nrows(    sexp_pair_mtx);
  int     weighted     = asLogical(sexp_weighted);
  int     n_threads    = asInteger(sexp_n_threads);
  
  
  // allocate the return vector
  SEXP sexp_result = PROTECT(allocVector(REALSXP, n_pairs));
  double *result   = REAL(sexp_result);
  
  // intermediary values
  double *sample_depths  = calloc(n_samples,           sizeof(double));
  int    *child_at       = calloc(n_edges + 1,         sizeof(int));
  int    *edge_to_leaves = calloc(n_edges * n_otus,    sizeof(int));
  double *sample_edge_wt = calloc(n_edges * n_samples, sizeof(double));
  
  if (any_null(sample_depths, child_at, edge_to_leaves, sample_edge_wt)) {
    free_4(sample_depths, child_at, edge_to_leaves, sample_edge_wt);
    error("Unable to allocate memory for UniFrac calculation.");
    return R_NilValue;
  }
  
  
  // common values for all the threads
  unifrac_t setup;
  setup.otus           = otus;
  setup.samples        = samples;
  setup.values         = values;
  setup.n_otus         = n_otus;
  setup.n_samples      = n_samples;
  setup.pair_mtx       = pair_mtx;
  setup.edge_mtx       = edge_mtx;
  setup.edge_lengths   = edge_lengths;
  setup.weighted       = weighted;
  setup.n_threads      = n_threads;
  setup.n_values       = n_values;
  setup.n_pairs        = n_pairs;
  setup.n_edges        = n_edges;
  setup.sample_depths  = sample_depths;
  setup.child_at       = child_at;
  setup.edge_to_leaves = edge_to_leaves;
  setup.sample_edge_wt = sample_edge_wt;
  setup.result         = result;
  
  
  
  
  //======================================================
  // Count the number of sequences in each sample.
  //======================================================
  
  memset(sample_depths, 0, n_samples * sizeof(double));
  for (int i = 0; i < n_values; i++) {
    sample_depths[samples[i] - 1] += values[i];
  }
  
  
  
  //======================================================
  // Index where a given node is the child. Root node = -1
  //======================================================
  
  // IntegerVector childAt = IntegerVector(max(treeEdge) + 1, -1);
  // for (int i = 0; i < nTreeEdges; i++) {
  //   childAt[treeEdge(i, 1)] = i;
  // }
  
  for (int edge = 0; edge <= n_edges; edge++) { child_at[edge] = -1; }
  for (int edge = 0; edge <  n_edges; edge++) { 
    int node = edge_mtx[1 * n_edges + edge] - 1;
    child_at[node] = edge;
  }
  
  
  // Run WITH multithreading
  #ifdef HAVE_PTHREAD
    if (n_threads > 1) {
      
      // threads and their individual input arguments
      pthread_t *tids = calloc(n_threads, sizeof(pthread_t));
      unifrac_t *args = calloc(n_threads, sizeof(unifrac_t));
      
      if (tids == NULL || args == NULL) {
        free_4(sample_depths, child_at, edge_to_leaves, sample_edge_wt);
        free(tids); free(args);
        error("Unable to allocate memory for multithreaded UniFrac calculation.");
        return R_NilValue;
      }
      
      int i, n = n_threads;
      
      for (i = 0; i < n; i++) memcpy(args + i, &setup, sizeof(unifrac_t));
      for (i = 0; i < n; i++) args[i].thread_i = i;
      
      for (i = 0; i < n; i++) pthread_create(&tids[i], NULL, calc_edge_to_leaves, &args[i]);
      for (i = 0; i < n; i++) pthread_join(   tids[i], NULL);
      
      for (i = 0; i < n; i++) pthread_create(&tids[i], NULL, calc_sample_edge_wt, &args[i]);
      for (i = 0; i < n; i++) pthread_join(   tids[i], NULL);
      
      for (i = 0; i < n; i++) pthread_create(&tids[i], NULL, calc_result, &args[i]);
      for (i = 0; i < n; i++) pthread_join(   tids[i], NULL);
      
      
      free_4(sample_depths, child_at, edge_to_leaves, sample_edge_wt);
      free(tids); free(args);
      
      UNPROTECT(1);
      
      return sexp_result;
    }
  #endif
  
  
  // Run WITHOUT multithreading
  setup.n_threads = 1;
  setup.thread_i  = 0;
  calc_edge_to_leaves(&setup);
  calc_sample_edge_wt(&setup);
  calc_result(&setup);
  
  free_4(sample_depths, child_at, edge_to_leaves, sample_edge_wt);
  
  UNPROTECT(1);
  
  return sexp_result;
}


