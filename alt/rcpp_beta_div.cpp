#include <Rcpp.h>
#include <math.h>
#include <string.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector rcpp_beta_div(int iJob, int nJobs, List sparseMatrix, const char* method, bool weighted) {
  
  IntegerVector samp = as<IntegerVector>(sparseMatrix["j"]);
  IntegerVector taxa = as<IntegerVector>(sparseMatrix["i"]);
  NumericVector amt  = as<NumericVector>(sparseMatrix["v"]);
  
  List dimnames               = as<List>(sparseMatrix["dimnames"]);
  CharacterVector sampleNames = as<CharacterVector>(dimnames[1]);
  
  // Ennumerate the different algorithms
  int fn = 0;
  if (( weighted) && (strcmp(method, "bray-curtis") == 0)) { fn = 1; }
  if ((!weighted) && (strcmp(method, "bray-curtis") == 0)) { fn = 2; }
  if (( weighted) && (strcmp(method, "euclidean")   == 0)) { fn = 3; }
  if ((!weighted) && (strcmp(method, "euclidean")   == 0)) { fn = 4; }
  if (( weighted) && (strcmp(method, "manhattan")   == 0)) { fn = 5; }
  if ((!weighted) && (strcmp(method, "manhattan")   == 0)) { fn = 6; }
  if (( weighted) && (strcmp(method, "jaccard")     == 0)) { fn = 1; } // jaccard = 2*bray/(1+bray)
  if ((!weighted) && (strcmp(method, "jaccard")     == 0)) { fn = 2; } // jaccard = 2*bray/(1+bray)
  
  int n = as<int>(sparseMatrix["ncol"]); // Number of distinct samples
  
  // Index the starting location of each sample in the vectors
  // Assumes these vectors are sorted by j (samples), then i (taxa)
  IntegerVector sPtr(n + 1);
  for (int i = 0; i < samp.size(); i++) {
    sPtr[samp[i]] = i + 1;
  }
  
  // R style dist object, initialized to zero
  int nIndices = (n * n - n) / 2;
  NumericVector results(nIndices);
  results.attr("class")  = "dist";
  results.attr("Size")   = n;
  results.attr("Upper")  = false;
  results.attr("Diag")   = false;
  results.attr("Labels") = sampleNames;
  
  
  // Find the starting and ending indices for this job
  int startIdx = floor(nIndices * ((iJob - 1) / static_cast<double>(nJobs)));
  int lastIdx  = floor(nIndices * (iJob / static_cast<double>(nJobs))) - 1;
  int nIters   = lastIdx - startIdx + 1;
  
  
  // Which two samples are we comparing first?
  int i  = startIdx + 1;
  int s1 = ceil(.5 * (-1 * pow(-8 * (i - 1) + 4 * pow(n, 2) - 4 * n - 7, .5) + 2 * n - 1) - 1) + 1;
  int s2 = n - (s1 * (n - 1 - s1) + (s1 * (s1 + 1)) / 2) + i;
  
  // Convert from 1-based to 0-based indexing
  i--;
  s1--;
  s2--;
  
  // These get re-incremented on loop entry
  i--;
  s2--;
  
  while (nIters-- > 0) {
    
    i++;
    if (++s2 >= n) {
      s1++;
      s2 = s1 + 1;
    }
    
    // What are these samples' first/last positions in the taxa and amt vectors?
    int s1_first_idx = sPtr[s1];
    int s1_last_idx  = sPtr[s1 + 1] - 1;
    int s2_first_idx = sPtr[s2];
    int s2_last_idx  = sPtr[s2 + 1] - 1;
    
    // The number of taxa found for sample 1 and 2
    int s1_nTaxa = s1_last_idx - s1_first_idx + 1;
    int s2_nTaxa = s2_last_idx - s2_first_idx + 1;
    if (s1_nTaxa == 0 && s2_nTaxa == 0) { continue; }
    
    double amt_sum  = 0.0; // amt_sum  = sum(x+y)
    double amt_diff = 0.0; // amt_diff = sum(abs(x-y))
    int common_nTaxa = 0;  // common_nTaxa = sum(x&y)
    
    
    switch (fn) {
    
      //======================================================
      // Bray-Curtis Weighted
      // sum(abs(x-y))/sum(x+y)
      //======================================================
    
      case 1: {
        
        int y = s2_first_idx;
        for (int x = s1_first_idx; x <= s1_last_idx; x++) {
          
          amt_sum += amt[x];
          
          while (taxa[y] < taxa[x] && y <= s2_last_idx) {
            amt_sum  += amt[y];
            amt_diff += amt[y];
            y++;
          }
          
          if (taxa[y] == taxa[x]) {
            amt_sum  += amt[y];
            amt_diff += fabs(amt[x] - amt[y]);
            y++;
          } else {
            amt_diff += amt[x];
          }
        }
        for (; y <= s2_last_idx; y++) {
          amt_sum  += amt[y];
          amt_diff += amt[y];
        }
        
        double result = amt_diff / amt_sum;
        results[i] = result;
        break;
      }
      
      
      //======================================================
      // Bray-Curtis Unweighted
      // (sum(x>0)+sum(y>0)-2*sum(x&y))/(sum(x>0)+sum(y>0))
      //======================================================
      
      case 2: {
        
        int y = s2_first_idx;
        for (int x = s1_first_idx; x <= s1_last_idx; x++) {
          while (taxa[y] < taxa[x] && y <= s2_last_idx) { y++;                 }
          if (y > s2_last_idx)                          { break;               }
          if (taxa[y] == taxa[x])                       { common_nTaxa++; y++; }
        }
        
        double result = (s1_nTaxa + s2_nTaxa - 2.0 * common_nTaxa) / (s1_nTaxa + s2_nTaxa);
        results[i] = result;
        break;
      }
      
      
      //======================================================
      // Euclidean Weighted
      // sqrt(sum((x-y)^2))
      //======================================================
      
      case 3: {
        
        int y = s2_first_idx;
        for (int x = s1_first_idx; x <= s1_last_idx; x++) {
          
          while (taxa[y] < taxa[x] && y <= s2_last_idx) {
            amt_diff += pow(amt[y], 2);
            y++;
          }
          
          if (taxa[y] == taxa[x]) {
            amt_diff += pow(amt[x] - amt[y], 2);
            y++;
          } else {
            amt_diff += pow(amt[x], 2);
          }
        }
        for (; y <= s2_last_idx; y++) {
          amt_diff += pow(amt[y], 2);
        }
        
        double result = sqrt(amt_diff);
        results[i] = result;
        break;
      }
      
      
      //======================================================
      // Euclidean Unweighted
      // sqrt(sum(x>0)+sum(y>0)-2*sum(x&y))
      //======================================================
      
      case 4: {
        
        int y = s2_first_idx;
        for (int x = s1_first_idx; x <= s1_last_idx; x++) {
          while (taxa[y] < taxa[x] && y <= s2_last_idx) { y++;                 }
          if (y > s2_last_idx)                          { break;               }
          if (taxa[y] == taxa[x])                       { common_nTaxa++; y++; }
        }
        
        double result = sqrt(s1_nTaxa + s2_nTaxa - 2 * common_nTaxa);
        results[i] = result;
        break;
      }
      
      
      //======================================================
      // Manhattan Weighted
      // sum(abs(x-y))
      //======================================================
      
      case 5: {
        
        int y = s2_first_idx;
        for (int x = s1_first_idx; x <= s1_last_idx; x++) {
          
          while (taxa[y] < taxa[x] && y <= s2_last_idx) {
            amt_diff += amt[y];
            y++;
          }
          
          if (taxa[y] == taxa[x]) {
            amt_diff += fabs(amt[x] - amt[y]);
            y++;
          } else {
            amt_diff += amt[x];
          }
        }
        for (; y <= s2_last_idx; y++) {
          amt_diff += amt[y];
        }
        
        double result = amt_diff;
        results[i] = result;
        break;
      }
      
      
      //======================================================
      // Manhattan Unweighted
      // sum(x>0)+sum(y>0)-2*sum(x&y)
      //======================================================
      
      case 6: {
        
        int y = s2_first_idx;
        for (int x = s1_first_idx; x <= s1_last_idx; x++) {
          while (taxa[y] < taxa[x] && y <= s2_last_idx) { y++;                 }
          if (y > s2_last_idx)                          { break;               }
          if (taxa[y] == taxa[x])                       { common_nTaxa++; y++; }
        }
        
        double result = s1_nTaxa + s2_nTaxa - 2 * common_nTaxa;
        results[i] = result;
        break;
      }
      
      
      default: {
        results[i] = -1.0;
        break;
      }
    
    }
    
  }
  
  
  //======================================================
  // Convert Bray-Curtis to Jaccard
  // 2 * bray / (1 + bray)
  //======================================================
  
  if (strcmp(method, "jaccard") == 0) {
    for (int i = 0; i < results.size(); i++) {
      results[i] = 2 * results[i] / (1 + results[i]);
    }
  }
  
  
  return results;
}



// /*** R
//   library(rbiom)
//   library(slam)
//   
//   infile <- system.file("extdata", "hmp50.biom", package = "rbiom")
//   biom   <- read.biom(infile)
//   
//   otus   <- biom$counts
//   ord    <- order(otus$j, otus$i)
//   otus$i <- otus$i[ord]
//   otus$j <- otus$j[ord]
//   otus$v <- otus$v[ord]
//   t.mtx  <- t(as.matrix(otus))
//   
//   g_wdm <- distance(1L, 1L, otus, "bray-curtis", TRUE)
//   g_udm <- distance(1L, 1L, otus, "bray-curtis", FALSE)
//   h_wdm <- vegan::vegdist(t.mtx, method="bray", binary=FALSE)
//   h_udm <- vegan::vegdist(t.mtx, method="bray", binary=TRUE)
//   identical(as.vector(g_wdm), as.vector(h_wdm))
//   identical(as.vector(g_udm), as.vector(h_udm))
//   
//   g_wdm <- distance(1L, 1L, otus, "manhattan", TRUE)
//   g_udm <- distance(1L, 1L, otus, "manhattan", FALSE)
//   h_wdm <- vegan::vegdist(t.mtx, method="manhattan", binary=FALSE)
//   h_udm <- vegan::vegdist(t.mtx, method="manhattan", binary=TRUE)
//   identical(as.vector(g_wdm), as.vector(h_wdm))
//   identical(as.vector(g_udm), as.vector(h_udm))
//   
//   g_wdm <- distance(1L, 1L, otus, "euclidean", TRUE)
//   g_udm <- distance(1L, 1L, otus, "euclidean", FALSE)
//   h_wdm <- vegan::vegdist(t.mtx, method="euclidean", binary=FALSE)
//   h_udm <- vegan::vegdist(t.mtx, method="euclidean", binary=TRUE)
//   identical(as.vector(g_wdm), as.vector(h_wdm))
//   identical(as.vector(g_udm), as.vector(h_udm))
//   
//   g_wdm <- distance(1L, 1L, otus, "jaccard", TRUE)
//   g_udm <- distance(1L, 1L, otus, "jaccard", FALSE)
//   h_wdm <- vegan::vegdist(t.mtx, method="jaccard", binary=FALSE)
//   h_udm <- vegan::vegdist(t.mtx, method="jaccard", binary=TRUE)
//   identical(as.vector(g_wdm), as.vector(h_wdm))
//   identical(as.vector(g_udm), as.vector(h_udm))
//   
//   rbenchmark::benchmark(
//     columns = c("test", "replications", "elapsed", "relative"),
//     distance(1L, 1L, otus, "braycurtis", TRUE),
//     distance(1L, 1L, otus, "braycurtis", FALSE),
//     distance(1L, 1L, otus, "manhattan",  TRUE),
//     distance(1L, 1L, otus, "manhattan",  FALSE),
//     distance(1L, 1L, otus, "euclidean",  TRUE),
//     distance(1L, 1L, otus, "euclidean",  FALSE),
//     distance(1L, 1L, otus, "jaccard",    TRUE),
//     distance(1L, 1L, otus, "jaccard",    FALSE),
//     vegan::vegdist(t.mtx, method="bray",      binary=FALSE),
//     vegan::vegdist(t.mtx, method="bray",      binary=TRUE),
//     vegan::vegdist(t.mtx, method="manhattan", binary=FALSE),
//     vegan::vegdist(t.mtx, method="manhattan", binary=TRUE),
//     vegan::vegdist(t.mtx, method="euclidean", binary=FALSE),
//     vegan::vegdist(t.mtx, method="euclidean", binary=TRUE),
//     vegan::vegdist(t.mtx, method="jaccard",   binary=FALSE),
//     vegan::vegdist(t.mtx, method="jaccard",   binary=TRUE) )
// */
