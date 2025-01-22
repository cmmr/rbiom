#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // calloc, free
#include <string.h> // memcpy
#include "get.h"
#include "pair.h"

// Detect if pthread is available.
#if defined __has_include
#  if __has_include (<pthread.h>)
#    include <pthread.h>
#    define HAVE_PTHREAD
#  endif
#endif


typedef struct {
  int    edge;
  int    parent;
  double length;
} node_t;


typedef struct {
  double      *otu_mtx;
  int          n_otus;
  int          n_samples;
  int          n_edges;
  double      *edge_lengths;
  int_pair_t  *pairs;
  int          n_pairs;
  int          weighted;
  int          normalized;
  double      *weight_mtx;
  node_t      *nodes;
  double      *sample_norms;
  int          n_threads;
  int          thread_i;
  double      *result;
} unifrac_t;



//======================================================
// Find the tree's edge weights for each sample.
//======================================================

void *unifrac_weight_mtx(void *arg) {

  unifrac_t setup = *((unifrac_t *) arg);
  
  double      *otu_mtx      = setup.otu_mtx;
  int          n_otus       = setup.n_otus;
  int          n_samples    = setup.n_samples;
  int          n_edges      = setup.n_edges;
  int          weighted     = setup.weighted;
  double      *weight_mtx   = setup.weight_mtx;
  node_t      *nodes        = setup.nodes;
  double      *sample_norms = setup.sample_norms;
  int          n_threads    = setup.n_threads;
  int          thread_i     = setup.thread_i;
  
  for (int sample = thread_i; sample < n_samples; sample += n_threads) {
    
    // initialize calloc()ed memory to zeroes
    for (int edge = 0; edge < n_edges; edge++) {
      weight_mtx[sample * n_edges + edge] = 0;
    }
    
    double sample_depth  = 0; // for relative abundance calculation
    sample_norms[sample] = 0; // denominator of normalization step
    
    if (weighted) {            
      for (int otu = 0; otu < n_otus; otu++) {
        sample_depth += otu_mtx[sample * n_otus + otu];
      }
    }
    
    for (int otu = 0; otu < n_otus; otu++) {
      
      double abundance = otu_mtx[sample * n_otus + otu];
      if (abundance == 0) continue; // OTU not present in sample
      
      int node = otu; // start at OTU tip/leaf in tree
      while (node > -1) { // traverse until we hit the tree's root
        
        int    edge   = nodes[node].edge;
        double length = nodes[node].length;
        int    idx    = sample * n_edges + edge;
        
        if (weighted) {
          
          // relative abundance, weighted by branch length
          double relative_abundance = abundance / sample_depth;
          double weighted_abundance = length * relative_abundance;
          weight_mtx[idx]      += weighted_abundance;
          sample_norms[sample] += weighted_abundance;
          
        } else {
          
          if (weight_mtx[idx]) break; // already traversed
          
          weight_mtx[idx] = 1;
        }
        
        node = nodes[node].parent; // proceed on up the tree
      }
      
    }
  }
  
  
  return NULL;
}



//======================================================
// Second pass, now with weight_mtx.
//======================================================

void *unifrac_result(void *arg) {
  
  unifrac_t setup = *((unifrac_t *) arg);
  
  int          n_edges      = setup.n_edges;
  double      *edge_lengths = setup.edge_lengths;
  int_pair_t  *pairs        = setup.pairs;
  int          n_pairs      = setup.n_pairs;
  int          weighted     = setup.weighted;
  int          normalized   = setup.normalized;
  double      *weight_mtx   = setup.weight_mtx;
  double      *sample_norms = setup.sample_norms;
  int          n_threads    = setup.n_threads;
  int          thread_i     = setup.thread_i;
  double      *result       = setup.result;
  
  
  //======================================================
  // WEIGHTED UniFrac
  //======================================================
  if (weighted) {
    
    for (int pair = thread_i; pair < n_pairs; pair += n_threads) {
      
      int sample_1 = pairs[pair].A - 1;
      int sample_2 = pairs[pair].B - 1;
      
      double sum = 0;
      
      for (int edge = 0; edge < n_edges; edge++) {
        
        double weight_1 = weight_mtx[sample_1 * n_edges + edge];
        double weight_2 = weight_mtx[sample_2 * n_edges + edge];
        
        if (weight_1 > weight_2) { sum += weight_1 - weight_2; }
        else                     { sum += weight_2 - weight_1; }
        
      }
      
      if (normalized) {
        sum /= sample_norms[sample_1] + sample_norms[sample_2];
      }
      
      result[pair] = sum;
    }
    
  }
  
  
  //======================================================
  // UNWEIGHTED UniFrac
  //======================================================
  else {
    
    for (int pair = thread_i; pair < n_pairs; pair += n_threads) {
      
      int sample_1 = pairs[pair].A - 1;
      int sample_2 = pairs[pair].B - 1;
      
      double distinct = 0, shared = 0;
      
      for (int edge = 0; edge < n_edges; edge++) {
        
        double length   = edge_lengths[edge];
        double weight_1 = weight_mtx[sample_1 * n_edges + edge];
        double weight_2 = weight_mtx[sample_2 * n_edges + edge];
        
        if      (weight_1 && weight_2) { shared   += length; }
        else if (weight_1 || weight_2) { distinct += length; }
      }
      
      
      result[pair] = distinct / (distinct + shared);
    }
    
  }
  
  return NULL;
}



//======================================================
// R interface. Dispatches threads on compute methods.
//======================================================
SEXP C_unifrac(
    SEXP sexp_otu_mtx,  SEXP sexp_phylo_tree, SEXP sexp_pair_mtx, 
    SEXP sexp_weighted, SEXP sexp_normalized, SEXP sexp_n_threads ) {
  
  double     *otu_mtx      = REAL( sexp_otu_mtx);
  int         n_otus       = nrows(sexp_otu_mtx);
  int         n_samples    = ncols(sexp_otu_mtx);
  
  int        *edge_mtx     = INTEGER(get(sexp_phylo_tree, "edge"));
  int         n_edges      = nrows(  get(sexp_phylo_tree, "edge"));
  double     *edge_lengths = REAL(   get(sexp_phylo_tree, "edge.length"));
  
  int_pair_t *pairs        = (int_pair_t *)INTEGER(sexp_pair_mtx);
  int         n_pairs      = ncols(sexp_pair_mtx);
  
  int         weighted     = asInteger(sexp_weighted);
  int         normalized   = asInteger(sexp_normalized);
  int         n_threads    = asInteger(sexp_n_threads);
  
  
  // allocate the return vector
  SEXP sexp_result = PROTECT(allocVector(REALSXP, n_pairs));
  double *result   = REAL(sexp_result);
  
  
  // intermediary values
  double  *weight_mtx   = calloc(n_samples * n_edges, sizeof(double));
  double  *sample_norms = calloc(n_samples,           sizeof(double));
  node_t  *nodes        = calloc(n_edges,             sizeof(node_t));
  
  if (weight_mtx == NULL || sample_norms == NULL || nodes == NULL) {
    free(weight_mtx); free(sample_norms); free(nodes);
    error("Unable to allocate memory for UniFrac calculation.");
    return R_NilValue;
  }
  
  
  // sort edge data by child node
  for (int edge = 0; edge < n_edges; edge++) {
    
    int    parent = edge_mtx[0 * n_edges + edge] - 2;
    int    child  = edge_mtx[1 * n_edges + edge] - 1;
    double length = edge_lengths[edge];
    
    if (child  > n_otus) child--;
    if (parent < n_otus) parent = -1;
    
    nodes[child].edge   = edge;
    nodes[child].parent = parent;
    nodes[child].length = length;
  }
  
  
  // common values for all the threads
  unifrac_t setup;
  
  setup.otu_mtx      = otu_mtx;
  setup.n_otus       = n_otus;
  setup.n_samples    = n_samples;
  setup.n_edges      = n_edges;
  setup.edge_lengths = edge_lengths;
  setup.pairs        = pairs;
  setup.n_pairs      = n_pairs;
  setup.weighted     = weighted;
  setup.normalized   = normalized;
  setup.weight_mtx   = weight_mtx;
  setup.nodes        = nodes;
  setup.sample_norms = sample_norms;
  setup.n_threads    = n_threads;
  setup.result       = result;
  
  
  
  // Run WITH multithreading
  #ifdef HAVE_PTHREAD
    if (n_threads > 1) {
      
      // threads and their individual input arguments
      pthread_t *tids = calloc(n_threads, sizeof(pthread_t));
      unifrac_t *args = calloc(n_threads, sizeof(unifrac_t));
      
      if (tids == NULL || args == NULL) {
        free(weight_mtx); free(sample_norms); free(nodes); free(tids); free(args);
        error("Unable to allocate memory for multithreaded UniFrac calculation.");
        return R_NilValue;
      }
      
      int i, n = n_threads;
      
      for (i = 0; i < n; i++) memcpy(args + i, &setup, sizeof(unifrac_t));
      for (i = 0; i < n; i++) args[i].thread_i = i;
      
      for (i = 0; i < n; i++) pthread_create(&tids[i], NULL, unifrac_weight_mtx, &args[i]);
      for (i = 0; i < n; i++) pthread_join(   tids[i], NULL);
      
      for (i = 0; i < n; i++) pthread_create(&tids[i], NULL, unifrac_result, &args[i]);
      for (i = 0; i < n; i++) pthread_join(   tids[i], NULL);
      
      
      free(weight_mtx); free(sample_norms); free(nodes);
      free(tids); free(args);
      
      UNPROTECT(1);
      
      return sexp_result;
    }
  #endif
  
  
  // Run WITHOUT multithreading
  setup.n_threads = 1;
  setup.thread_i  = 0;
  unifrac_weight_mtx(&setup);
  unifrac_result(&setup);
  
  free(weight_mtx); free(sample_norms); free(nodes);
  
  UNPROTECT(1);
  
  return sexp_result;
}


/*
 
#====================================
# Example output from phyloseq
#====================================

# library(phyloseq)
# data(esophagus)
#
# > UniFrac(esophagus, weighted=FALSE)
#           B         C
# C 0.5175550          
# D 0.5182284 0.5422394
#
# > UniFrac(esophagus, weighted=TRUE, normalized=FALSE)
#           B         C
# C 0.1050480          
# D 0.1401124 0.1422409
#
# > UniFrac(esophagus, weighted=TRUE, normalized=TRUE)
#           B         C
# C 0.2035424          
# D 0.2603371 0.2477016



#====================================
# Simple Implementation in R
#====================================

unifrac_r <- function (biom) {
  
  biom <- rbiom::as_rbiom(biom)
  
  otu_mtx   <- as.matrix(biom$counts)
  n_samples <- biom$n_samples
  n_otus    <- biom$n_otus
  edges     <- biom$tree$edge
  lengths   <- biom$tree$edge.length
  n_edges   <- nrow(edges)
  
  
  m <- array(0, c(n_samples, n_edges))
  
  for (sample in seq_len(n_samples)) {
    
    sample_depth <- 0
    for (otu in seq_len(n_otus)) {
      abundance    <- otu_mtx[otu,sample]
      sample_depth <- sample_depth + abundance
    }
    
    for (otu in seq_len(n_otus)) {
      
      abundance <- otu_mtx[otu,sample]
      if (abundance == 0) next
      
      abundance <- abundance / sample_depth
      
      edge <- which(edges[,2] == otu)
      while (length(edge) == 1) {
        
        len <- lengths[edge]
        m[sample,edge] <- m[sample,edge] + (len * abundance)
        
        parent <- edges[edge,1]
        edge   <- which(edges[,2] == parent)
      }
        
    }
    
  }
  
  pairs   <- t(combn(seq_len(n_samples), 2))
  n_pairs <- nrow(pairs)
  
  print('Unweighted:')
  for (p in seq_len(n_pairs)) {
    i <- m[pairs[p,1],]
    j <- m[pairs[p,2],]
    print(sum(lengths[xor(i,j)]) / sum(lengths[i | j]))
  }
  
  print('\n')
  print('Weighted Non-Normalized:')
  for (p in seq_len(n_pairs)) {
    i <- m[pairs[p,1],]
    j <- m[pairs[p,2],]
    print(sum(abs(i - j)))
  }
  
  print('\n')
  print('Weighted Normalized:')
  for (p in seq_len(n_pairs)) {
    i <- m[pairs[p,1],]
    j <- m[pairs[p,2],]
    print(sum(abs(i - j)) / sum(i + j))
  }
        
}

*/
 
