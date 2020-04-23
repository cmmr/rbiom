// [[Rcpp::depends(RcppParallel)]]
#include <math.h>
#include <stdio.h>
#include <Rcpp.h>
#include <RcppParallel.h>
using namespace Rcpp;
using namespace RcppParallel;



//======================================================
// Parallel Task 1: Find the leaves under each node
//======================================================

struct TraverseTree : public Worker
{
  const RVector<int> childAt;
  const RMatrix<int> treeEdge;
        RMatrix<int> edge2leaves;
  
  TraverseTree(
    const IntegerVector childAt, 
    const IntegerMatrix treeEdge, 
          IntegerMatrix edge2leaves
  ) : 
    childAt(childAt), 
    treeEdge(treeEdge),
    edge2leaves(edge2leaves)
  {}
  
  void operator()(std::size_t begin_in, std::size_t end_in) {
    
    int begin = static_cast<int>(begin_in);
    int end   = static_cast<int>(end_in);
    
    for (int leaf = begin; leaf < end; leaf++) {
      int node = leaf;
      do {
        edge2leaves(childAt[node], leaf - 1) = 1;
        node = treeEdge(childAt[node], 0);
      } while (childAt[node] != -1);
    }
  }
};



//======================================================
// Parallel Task 2: Assign Weights to each Tree Edge
//======================================================

struct WeighEdges : public Worker
{
  const RMatrix<int>    edge2leaves;
  const RVector<int>    weighted;
  const RVector<int>    mtxSample;
  const RVector<int>    mtxTaxa;
  const RVector<double> mtxAbundance;
  const RVector<double> treeLengths;
  const RVector<double> nSeq;
        RMatrix<double> sampleEdgeWt;
  
  WeighEdges(
    const IntegerMatrix edge2leaves,
    const IntegerVector weighted,
    const IntegerVector mtxSample,
    const IntegerVector mtxTaxa,
    const NumericVector mtxAbundance,
    const NumericVector treeLengths,
    const NumericVector nSeq,
          NumericMatrix sampleEdgeWt
  ) :
    edge2leaves(edge2leaves),
    weighted(weighted),
    mtxSample(mtxSample),
    mtxTaxa(mtxTaxa),
    mtxAbundance(mtxAbundance),
    treeLengths(treeLengths),
    nSeq(nSeq),
    sampleEdgeWt(sampleEdgeWt)
  {}
  
  void operator()(std::size_t begin_in, std::size_t end_in) {
    
    int begin = static_cast<int>(begin_in);
    int end   = static_cast<int>(end_in);
    
    int nTreeEdges = sampleEdgeWt.ncol();
    
    
    if (weighted[0] == 1) {
      
      double wt = 0.0;
      
      for (int i = begin; i < end; i++) {
        for (int edge = 0; edge < nTreeEdges; edge++) {
          if (edge2leaves(edge, mtxTaxa[i] - 1) == 0) continue;
          
          wt = treeLengths[edge] * (mtxAbundance[i] / nSeq[mtxSample[i]]);
          sampleEdgeWt(mtxSample[i] - 1, edge) += wt;
        }
      }
    } else {
      
      for (int i = begin; i < end; i++) {
        for (int edge = 0; edge < nTreeEdges; edge++) {
          if (edge2leaves(edge, mtxTaxa[i] - 1) == 0) continue;
          
          sampleEdgeWt(mtxSample[i] - 1, edge) = treeLengths[edge];
        }
      }
    }
    
  }
};



//======================================================
// Parallel Task 3: Pairwise sample distances
//======================================================

struct PairwiseDist : public Worker
{
  const RMatrix<double> sampleEdgeWt;
  const RVector<int>    weighted;
        RVector<double> results;
  
  PairwiseDist(
    const NumericMatrix sampleEdgeWt,
    const IntegerVector weighted,
          NumericVector results
  ) :
    sampleEdgeWt(sampleEdgeWt),
    weighted(weighted),
    results(results)
  {}
  
  void operator()(std::size_t begin_in, std::size_t end_in) {
    
    int begin = static_cast<int>(begin_in);
    int end   = static_cast<int>(end_in);
    
    int nIters     = end - begin;
    int nSamples   = sampleEdgeWt.nrow();
    int nTreeEdges = sampleEdgeWt.ncol();
    
    
    // Which two samples are we comparing first?
    int i  = begin + 1;
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
    
    if (weighted[0] == 1) {
      
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
    
    
  }
};





// [[Rcpp::export]]
NumericVector par_unifrac(List sparseMatrix, List tree, IntegerVector weighted) {
  
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
  // Index where a given node is the child. Root node = -1
  //======================================================
  
  IntegerVector childAt = IntegerVector(max(treeEdge) + 1, -1);
  for (int i = 0; i < nTreeEdges; i++) {
    childAt[treeEdge(i, 1)] = i;
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
  // Trace the path from each leaf to the root node
  //======================================================
  
  IntegerMatrix edge2leaves = IntegerMatrix(nTreeEdges, nOTUs);
  TraverseTree traverseTree(
      childAt, 
      treeEdge, 
      edge2leaves);
  parallelFor(1, nOTUs + 1, traverseTree, 10000);
  
  
  //======================================================
  // Map samples to their branch weights
  //======================================================
  
  NumericMatrix sampleEdgeWt = NumericMatrix(nSamples, nTreeEdges);
  WeighEdges weighEdges(
      edge2leaves,
      weighted,
      mtxSample,
      mtxTaxa,
      mtxAbundance,
      treeLengths,
      nSeq,
      sampleEdgeWt);
  parallelFor(0, mtxSize, weighEdges, 10000);
  
  
  //======================================================
  // Compute Pairwise Distances
  //======================================================
  PairwiseDist pairwiseDist(
      sampleEdgeWt,
      weighted,
      results);
  parallelFor(0, results.size(), pairwiseDist, 10000);
  
  
  return results;
}

