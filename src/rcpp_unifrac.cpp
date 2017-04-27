#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;

  

// [[Rcpp::export]]
NumericVector rcpp_unifrac(int iJob, int nJobs, List sparseMatrix, List tree, bool weighted) {
  
  //======================================================
  // Assumptions:
  //   - sparseMatrix sorted by j (samples), then i (taxa)
  //   - tip.labels are in same order as sparseMatrix$j
  //======================================================

  
  //======================================================
  // Convert biom and tree to C++ objects
  //======================================================
  
  IntegerVector mtxSample    = as<IntegerVector>(sparseMatrix["j"]);
  IntegerVector mtxTaxa      = as<IntegerVector>(sparseMatrix["i"]);
  NumericVector mtxAbundance = as<NumericVector>(sparseMatrix["v"]);
  
  int mtxSize  = mtxAbundance.size();
  int nOTUs    = as<int>(sparseMatrix["nrow"]);
  int nSamples = as<int>(sparseMatrix["ncol"]);
  
  List            dimnames    = as<List>(sparseMatrix["dimnames"]);
  CharacterVector sampleNames = as<CharacterVector>(dimnames[1]);
  
  IntegerMatrix   treeEdge    = as<IntegerMatrix>(tree["edge"]);
  NumericVector   treeLengths = as<NumericVector>(tree["edge.length"]);
  int             nTreeEdges  = treeEdge.nrow();
  
  
  //======================================================
  // Index the starting offset of each sample
  //======================================================
  IntegerVector sPtr(nSamples + 1);
  for (int i = 0; i < mtxSize; i++) {
    sPtr[mtxSample[i]] = i + 1;
  }
  
  
  //======================================================
  // Count the number of sequences in each sample
  //======================================================
  NumericVector nSeq(nSamples + 1);
  for (int i = 0; i < mtxSize; i++) {
    nSeq[mtxSample[i]] += mtxAbundance[i];
  }
  
  
  //======================================================
  // R style dist object, initialized to zero
  //======================================================
  int nIndices = (nSamples * nSamples - nSamples) / 2;
  NumericVector results(nIndices);
  results.attr("class")  = "dist";
  results.attr("Size")   = nSamples;
  results.attr("Upper")  = false;
  results.attr("Diag")   = false;
  results.attr("Labels") = sampleNames;
  
  
  
  //======================================================
  // Map edges/branches to the leaves that are under them.
  //======================================================
  
  LogicalMatrix edge2leaves = LogicalMatrix(nTreeEdges, nOTUs);
  NumericVector otuDepths   = NumericVector(nOTUs);
  
  
  // Index where a given node is the child. Root node == -1
  IntegerVector childAt = IntegerVector(max(treeEdge) + 1, -1);
  for (int i = 0; i < nTreeEdges; i++) {
    childAt[treeEdge(i, 1)] = i;
  }
  
  // Trace the path from each leaf to the root node
  for (int leaf = 1; leaf <= nOTUs; leaf++) {
    int node = leaf;
    do {
      edge2leaves(childAt[node], leaf - 1) = true;
      node = treeEdge(childAt[node], 0);
    } while (childAt[node] != -1);
  }
  
  
  //======================================================
  // Map samples to their branch weights
  //======================================================
  
  NumericMatrix sampleEdgeWt = NumericMatrix(nSamples, nTreeEdges);
  
  if (weighted) {
    
    double wt = 0.0;
    
    for (int i = 0; i < mtxSize; i++) {
      for (int edge = 0; edge < nTreeEdges; edge++) {
        if (edge2leaves(edge, mtxTaxa[i] - 1) == false) continue;
        
        wt = treeLengths[edge] * (mtxAbundance[i] / nSeq[mtxSample[i]]);
        sampleEdgeWt(mtxSample[i] - 1, edge) += wt;
      }
    }
  } else {
    
    for (int i = 0; i < mtxSize; i++) {
      for (int edge = 0; edge < nTreeEdges; edge++) {
        if (edge2leaves(edge, mtxTaxa[i] - 1) == false) continue;
        
        sampleEdgeWt(mtxSample[i] - 1, edge) = treeLengths[edge];
      }
    }
  }
  
  
  
  // Find the starting and ending indices for this job
  int startIdx = floor(nIndices * ((iJob - 1) / static_cast<double>(nJobs)));
  int lastIdx  = floor(nIndices * (iJob / static_cast<double>(nJobs))) - 1;
  int nIters   = lastIdx - startIdx + 1;
  
  
  // Which two samples are we comparing first?
  int i  = startIdx + 1;
  int n  = nSamples;
  int s1 = ceil(.5 * (-1 * pow(-8 * (i - 1) + 4 * pow(n, 2) - 4 * n - 7, .5) + 2 * n - 1) - 1) + 1;
  int s2 = n - (s1 * (n - 1 - s1) + (s1 * (s1 + 1)) / 2) + i;
  
  // Convert from 1-based to 0-based indexing
  i--;
  s1--;
  s2--;
  
  // This gets re-incremented on loop entry
  s2--;
  
    
  double x = 0, y = 0, amt_diff = 0, amt_total = 0;
  
  if (weighted) {
    
    while (nIters-- > 0) { 
      if (++s2 >= nSamples) { s1++; s2 = s1 + 1; }
      
      amt_diff = 0;
      for (int edge = 0; edge < nTreeEdges; edge++) {
        x = sampleEdgeWt(s1, edge);
        y = sampleEdgeWt(s2, edge);
        amt_diff += (x >= y) ? (x - y) : (y - x);
      }
      
      results[i++] = amt_diff;
    }
    
  } else {
    
    while (nIters-- > 0) {
      if (++s2 >= nSamples) { s1++; s2 = s1 + 1; }
      
      amt_diff = 0; amt_total = 0;
      for (int edge = 0; edge < nTreeEdges; edge++) {
        x = sampleEdgeWt(s1, edge);
        y = sampleEdgeWt(s2, edge);
        amt_diff  += (x && y) ? (0) : (x + y);
        amt_total += (x && y) ? (x) : (x + y);
      }
      
      results[i++] = amt_diff / amt_total;
    }
    
  }
  
  return results;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

// /*** R
// 
// library(slam)
// library(rbiom)
// 
// infile        <- system.file("extdata", "hmp50.biom", package = "rbiom")
// biom          <- select(read.biom(infile), 1:3)
// 
// biom$counts   <- biom$counts[as.character(biom$phylogeny$tip.label),]
// 
// ord           <- order(biom$counts$j, biom$counts$i)
// biom$counts$i <- biom$counts$i[ord]
// biom$counts$j <- biom$counts$j[ord]
// biom$counts$v <- biom$counts$v[ord]
// 
// rcpp_unifrac(1L, 1L, biom$counts, biom$phylogeny, TRUE)
// 
// */


//============================================
// QIIME Weighted Unifrac --------------------
//         HMP01            HMP02
// HMP02   0.571402944323
// HMP03   0.726275816051   0.343551123533
// 
// 
// QIIME Unweighted Unifrac ------------------
//         HMP01            HMP02
// HMP02   0.478055388912
// HMP03   0.536688384399   0.340015248802
//============================================
