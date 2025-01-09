#include <Rcpp.h>
#include <RcppThread.h>
using namespace Rcpp;



// [[Rcpp::export]]
NumericVector parallel_unifrac(
    List sparseMatrix, IntegerMatrix pairs, List tree,
    bool weighted, int nThreads ) {


  //======================================================
  // Assumptions:
  //   - sparseMatrix sorted by j (samples), then i (taxa)
  //   - tip.labels are in same order as sparseMatrix$j
  //======================================================

  //======================================================
  // Convert slam matrix and phylo tree to C++ objects.
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
  // Default to using all available CPU cores
  //======================================================

  if (nThreads < 0)
    nThreads = std::thread::hardware_concurrency();


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
  // Trace the path from each leaf to the root node
  //======================================================

  IntegerMatrix edge2leaves = IntegerMatrix(nTreeEdges, nOTUs);
  RcppThread::parallelFor(1, nOTUs + 1, [&] (int leaf) {
    int node = leaf;
    do {
      edge2leaves(childAt[node], leaf - 1) = 1;
      node = treeEdge(childAt[node], 0);
    } while (childAt[node] != -1);
  }, nThreads);


  //======================================================
  // Map samples to their branch weights
  //======================================================

  NumericMatrix sampleEdgeWt = NumericMatrix(nSamples, nTreeEdges);

  if (weighted) {

    RcppThread::parallelFor(0, mtxSize, [&] (int i) {
      double wt = 0;
      for (int edge = 0; edge < nTreeEdges; edge++) {
        if (edge2leaves(edge, mtxTaxa[i] - 1) == 0) continue;
        wt = treeLengths[edge] * (mtxAbundance[i] / nSeq[mtxSample[i]]);
        sampleEdgeWt(mtxSample[i] - 1, edge) += wt;
      }
    }, nThreads);

  } else {

    RcppThread::parallelFor(0, mtxSize, [&] (int i) {
      for (int edge = 0; edge < nTreeEdges; edge++) {
        if (edge2leaves(edge, mtxTaxa[i] - 1) == 0) continue;
        sampleEdgeWt(mtxSample[i] - 1, edge) = treeLengths[edge];
      }
    }, nThreads);

  }


  //======================================================
  // Compute Pairwise Distances
  //======================================================

  NumericVector distances(pairs.nrow());

  if (weighted) {

    RcppThread::parallelFor(0, distances.size(), [&] (int i) {
      double x = 0, y = 0, amt_diff = 0;
      for (int edge = 0; edge < nTreeEdges; edge++) {
        x = sampleEdgeWt(pairs(i,0), edge);
        y = sampleEdgeWt(pairs(i,1), edge);
        amt_diff += (x >= y) ? (x - y) : (y - x);
      }
      distances[i] = amt_diff;
    }, nThreads);

  } else {

    RcppThread::parallelFor(0, distances.size(), [&] (int i) {
      double x = 0, y = 0, amt_diff = 0, amt_total = 0;
      for (int edge = 0; edge < nTreeEdges; edge++) {
        x = sampleEdgeWt(pairs(i,0), edge);
        y = sampleEdgeWt(pairs(i,1), edge);
        amt_diff  += (x && y) ? (0) : (x + y);
        amt_total += (x && y) ? (x) : (x + y);
      }
      distances[i] = amt_diff / amt_total;
    }, nThreads);

  }


  return distances;
}


/*

 =================================
 -- OLD CODE USING RcppParallel --
 =================================

 --> Converted to RcppThreads on 1/8/2025 to avoid UBSAN errors on CRAN checks.


  // [[Rcpp::depends(RcppParallel)]]
  #include <cmath>
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

      // Rcpp::checkUserInterrupt();

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

      // Rcpp::checkUserInterrupt();

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
    const RMatrix<int>    pairs;
    const RVector<int>    weighted;
    RVector<double> distances;

    PairwiseDist(
      const NumericMatrix sampleEdgeWt,
      const IntegerMatrix pairs,
      const IntegerVector weighted,
      NumericVector distances
    ) :
      sampleEdgeWt(sampleEdgeWt),
      pairs(pairs),
      weighted(weighted),
      distances(distances)
    {}

    void operator()(std::size_t begin, std::size_t end) {

      // Rcpp::checkUserInterrupt();

      int nTreeEdges = sampleEdgeWt.ncol();

      double x = 0, y = 0, amt_diff = 0, amt_total = 0;


      if (weighted[0] == 1) {

        for (std::size_t i = begin; i < end; i++) {

          amt_diff = 0;
          for (int edge = 0; edge < nTreeEdges; edge++) {
            x = sampleEdgeWt(pairs(i,0), edge);
            y = sampleEdgeWt(pairs(i,1), edge);
            amt_diff += (x >= y) ? (x - y) : (y - x);
          }
          distances[i] = amt_diff;

        }

      } else {

        for (std::size_t i = begin; i < end; i++) {

          amt_diff = 0; amt_total = 0;
          for (int edge = 0; edge < nTreeEdges; edge++) {
            x = sampleEdgeWt(pairs(i,0), edge);
            y = sampleEdgeWt(pairs(i,1), edge);
            amt_diff  += (x && y) ? (0) : (x + y);
            amt_total += (x && y) ? (x) : (x + y);
          }
          distances[i] = amt_diff / amt_total;
        }

      }

    }
  };





  // [[Rcpp::export]]
  NumericVector par_unifrac(List sparseMatrix, IntegerMatrix pairs, List tree, IntegerVector weighted) {

    //======================================================
    // Assumptions:
    //   - sparseMatrix sorted by j (samples), then i (taxa)
    //   - tip.labels are in same order as sparseMatrix$j
    //======================================================

    //======================================================
    // Convert slam matrix and phylo tree to C++ objects.
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
    // Trace the path from each leaf to the root node
    //======================================================

    IntegerMatrix edge2leaves = IntegerMatrix(nTreeEdges, nOTUs);
    TraverseTree traverseTree(childAt, treeEdge, edge2leaves);
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

    NumericVector distances(pairs.nrow());
    PairwiseDist pairwiseDist(sampleEdgeWt, pairs, weighted, distances);
    parallelFor(0, distances.size(), pairwiseDist, 10000);


    return distances;
  }

*/



