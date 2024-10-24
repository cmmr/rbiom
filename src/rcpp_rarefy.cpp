#include <Rcpp.h>
using namespace Rcpp;

#include <stdint.h> // uint64_t uint32_t


// [[Rcpp::export]]
List rcpp_rarefy(List sparseMatrix, int depth, int seed = 0) {
  
  
  //======================================================
  // Don't change the original object (passed by ref)
  //======================================================
  
  List returnMatrix = clone(sparseMatrix);
  
  
  //======================================================
  // Convert slam matrix to C++ objects
  //======================================================
  
  IntegerVector   mtxSample     = as<IntegerVector>(returnMatrix("j"));
  IntegerVector   mtxTaxa       = as<IntegerVector>(returnMatrix("i"));
  IntegerVector   mtxAbundance  = as<IntegerVector>(returnMatrix("v"));
  List            dimnames      = as<List>(returnMatrix["dimnames"]);
  CharacterVector taxaNames     = as<CharacterVector>(dimnames[0]);
  CharacterVector sampleNames   = as<CharacterVector>(dimnames[1]);
  
  
  //======================================================
  // Loop over each sample
  //======================================================
  
  int mtxSize   = mtxAbundance.size();
  int first_idx = min(mtxSample);
  int last_idx  = max(mtxSample);

  for (int sample_idx = first_idx; sample_idx <= last_idx; sample_idx++) {
    
    // if (sample_idx % 1000 == 0) {
    //   Rcpp::checkUserInterrupt();
    // }


    //======================================================
    // Count the number of sequences for this sample
    //======================================================

    int nSeqs = 0;
    for (int i = 0; i < mtxSize; i++) {
      if (mtxSample[i] == sample_idx) {
        nSeqs = nSeqs + mtxAbundance[i];
      }
    }


    //======================================================
    // Already rarefied - finished
    //======================================================

    if (nSeqs == depth) {
      continue;
    }


    //======================================================
    // Insufficient sequences; discard all abundances
    //======================================================

    if (nSeqs < depth) {
      for (int i = 0; i < mtxSize; i++) {
        if (mtxSample[i] == sample_idx) {
          mtxAbundance[i] = 0;
        }
      }
      continue;
    }


    //======================================================
    // Create a vector with all of this sample's sequences
    //======================================================

    IntegerVector taxa = IntegerVector(nSeqs);
    int n = 0;
    for (int i = 0; i < mtxSize; i++) {
      if (mtxSample[i] == sample_idx) {
        for (int j = 0; j < mtxAbundance[i]; j++) {
          taxa[n] = mtxTaxa[i];
          n++;
        }
        mtxAbundance[i] = 0;
      }
    }
    
    
    //======================================================
    // Middle Square Weyl Sequence PRNG
    //======================================================

    uint64_t prng_x = 0;
    uint64_t prng_w = 0;
    uint64_t prng_s = 0xb5ad4eceda1ce2a9;
    // UINT32_MAX = 0xFFFFFFFF = 4294967295;
    
    prng_x = ((uint64_t)nSeqs * (uint64_t)seed) % 4294967295;


    //======================================================
    // Randomly pick 'depth' number of sequences to keep
    //======================================================

    IntegerVector kept = IntegerVector(depth);
    
    for (int i = 0; i < depth; i++) {
      
      prng_x *= prng_x; 
      prng_x += (prng_w += prng_s); 
      prng_x = (prng_x>>32) | (prng_x<<32);
      
      int   j = ((uint32_t)(prng_x)) % (nSeqs - i);
      kept[i] = taxa[j];
      taxa[j] = taxa(nSeqs - i - 1);
    }

    IntegerVector tabulated = IntegerVector(mtxSize + 1);
    for (int i = 0; i < depth; i++) {
      tabulated(kept[i])++;
    }

    for (int i = 0; i < mtxSize; i++) {
      if (mtxSample[i] == sample_idx) {
        mtxAbundance[i] = tabulated(mtxTaxa[i]);
      }
    }


  }
  
  
  //======================================================
  // Drop any zero-abundance entries
  //======================================================
  
  LogicalVector nonZero = LogicalVector(mtxSize);
  for (int i = 0; i < mtxSize; i++) {
    nonZero[i] = mtxAbundance[i] > 0;
  }

  IntegerVector returnSample    = mtxSample[nonZero];
  IntegerVector returnTaxa      = mtxTaxa[nonZero];
  IntegerVector returnAbundance = mtxAbundance[nonZero];
  
  mtxSize = returnAbundance.size();

  int nInitSample = sampleNames.size();
  int nInitTaxa   = taxaNames.size();


  //======================================================
  // Drop zero-sum Samples and re-number 1:n
  //======================================================

  LogicalVector nonZeroSample = LogicalVector(nInitSample);
  nonZeroSample.fill(false);

  for (int i = 0; i < mtxSize; i++) {
    nonZeroSample[returnSample[i] - 1] = true;
  }

  int nextSampleIdx = 1;

  IntegerVector newSampleIdx = IntegerVector(nInitSample);
  for (int i = 0; i < nInitSample; i++) {
    if (nonZeroSample[i]) {
      newSampleIdx[i] = nextSampleIdx;
      nextSampleIdx++;
    }
  }

  for (int i = 0; i < mtxSize; i++) {
    returnSample[i] = newSampleIdx[returnSample[i] - 1];
  }

  CharacterVector returnSampleNames = sampleNames[nonZeroSample];



  //======================================================
  // Drop zero-sum Taxa and re-number 1:n
  //======================================================

  LogicalVector nonZeroTaxa = LogicalVector(nInitTaxa);
  nonZeroTaxa.fill(false);

  for (int i = 0; i < mtxSize; i++) {
    nonZeroTaxa[returnTaxa[i] - 1] = true;
  }

  int nextTaxaIdx = 1;
  IntegerVector newTaxaIdx = IntegerVector(nInitTaxa);
  for (int i = 0; i < nInitTaxa; i++) {
    if (nonZeroTaxa[i]) {
      newTaxaIdx[i] = nextTaxaIdx;
      nextTaxaIdx++;
    }
  }

  for (int i = 0; i < mtxSize; i++) {
    returnTaxa[i] = newTaxaIdx[returnTaxa[i] - 1];
  }

  CharacterVector returnTaxaNames = taxaNames[nonZeroTaxa];

  
  
  //======================================================
  // Update and return the new slam matrix
  //======================================================
  
  returnMatrix("j")        = returnSample;
  returnMatrix("i")        = returnTaxa;
  returnMatrix("v")        = returnAbundance;
  returnMatrix["nrow"]     = (int)(max(returnTaxa));
  returnMatrix["ncol"]     = (int)(max(returnSample));
  returnMatrix("dimnames") = List::create(returnTaxaNames, returnSampleNames);
  
  return returnMatrix;
}

