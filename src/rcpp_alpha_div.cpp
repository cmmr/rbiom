#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::export]]
List rcpp_alpha_div(List sparseMatrix) {
  
  
  //======================================================
  // Convert slam matrix to C++ objects
  //======================================================
  
  IntegerVector SampleIdx  = as<IntegerVector>(sparseMatrix("j"));
  IntegerVector TaxaIdx    = as<IntegerVector>(sparseMatrix("i"));
  NumericVector Abundances = as<NumericVector>(sparseMatrix("v"));
  
  List            dimnames    = as<List>(sparseMatrix["dimnames"]);
  CharacterVector SampleNames = as<CharacterVector>(dimnames[1]);
  
  int nSamples = SampleNames.size();
  int nEntries = Abundances.size();
  
  
  //======================================================
  // Initialize the output vectors
  //======================================================
  
  NumericVector FinalDepth      = NumericVector(nSamples);
  NumericVector FinalOTUs       = NumericVector(nSamples);
  NumericVector FinalShannon    = NumericVector(nSamples);
  NumericVector FinalChao1      = NumericVector(nSamples);
  NumericVector FinalSimpson    = NumericVector(nSamples);
  NumericVector FinalInvSimpson = NumericVector(nSamples);
  
  
  //======================================================
  // First Pass = Tabulate Depth and OTUs per Sample
  //======================================================
  
  for (int i = 0; i < nEntries; i++) {
    FinalDepth[SampleIdx[i] - 1] += Abundances[i];
    FinalOTUs[SampleIdx[i] - 1]++;
  }
  
  
  //======================================================
  // Loop over each sample
  //======================================================

  for (int s = 1; s <= nSamples; s++) {
    
    // if (s % 1000 == 0) {
    //   Rcpp::checkUserInterrupt();
    // }
    
    
    int   nOTUs = (int)(FinalOTUs[s - 1]);
    float Depth = (float)(FinalDepth[s - 1]);
    
    NumericVector x = NumericVector(nOTUs); // Observed Abundances
    NumericVector p = NumericVector(nOTUs); // Probability of OTUs
    
    
    
    //======================================================
    // Find all this sample's entries => Populate x and p
    //======================================================
    
    int n = 0;
    for (int i = 0; i < nEntries; i++) {
      if (SampleIdx[i] == s) {
        x[n] = Abundances[i];
        p[n] = Abundances[i] / Depth;
        n++;
      }
    }
    
    
    //======================================================
    // Calculate Shannon Diversity Index
    // Shannon <- -sum(p * log(p))
    //======================================================
    
    float total = 0;
    for (int i = 0; i < nOTUs; i++) {
      total += p[i] * log(p[i]);
    }
    FinalShannon[s - 1] = -1 * total;
    
    
    //======================================================
    // Calculate Chao1 Index
    // x     <- ceiling(x)
    // Chao1 <- nOTUs + (sum(x == 1) ** 2) / (2 * sum(x == 2))
    //======================================================
    
    float ones = 0;
    float twos = 0;
    NumericVector c = ceiling(x);
    for (int i = 0; i < nOTUs; i++) {
      if (c[i] == 1) ones++;
      if (c[i] == 2) twos++;
    }
    FinalChao1[s - 1] = nOTUs + (pow(ones, 2) / (2 * twos));
    
    
    //======================================================
    // Calculate Simpson Index and Inverse Simpson Index
    // Simpson    <- 1 - sum(p ** 2)
    // InvSimpson <- 1 / sum(p ** 2)
    //======================================================
    total = 0;
    for (int i = 0; i < nOTUs; i++) {
      total += pow(p[i], 2);
    }
    FinalSimpson[s - 1]    = 1 - total;
    FinalInvSimpson[s - 1] = 1 / total;
  
  }
  
  
  //======================================================
  // Create and return a data frame with the results
  //======================================================
  
  DataFrame df = DataFrame::create(
    Named(".sample")    = SampleNames,
    Named("Depth")      = FinalDepth,
    Named("OTUs")       = FinalOTUs,
    Named("Shannon")    = FinalShannon,
    Named("Chao1")      = FinalChao1,
    Named("Simpson")    = FinalSimpson,
    Named("InvSimpson") = FinalInvSimpson );
  
  return df;
}

