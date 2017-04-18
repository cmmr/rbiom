// #include <Rcpp.h>
// #include <math.h>
// #include <stdio.h>
// using namespace Rcpp;
// 
// // braycurtis_unweighted
// // (sum(x>0)+sum(y>0)-2*sum(x&y))/(sum(x>0)+sum(y>0))
// 
// 
// // [[Rcpp::export]]
// NumericVector braycurtis_unweighted(int mod, int rem, List sparseMatrix) {
//   
//   IntegerVector samp = as<IntegerVector>(sparseMatrix["j"]);
//   IntegerVector taxa = as<IntegerVector>(sparseMatrix["i"]);
//   //double amt = as<double>(sparseMatrix["v"]);
//   
//   List dimnames               = as<List>(sparseMatrix["dimnames"]);
//   CharacterVector sampleNames = as<CharacterVector>(dimnames[1]);
//   
//   int n = as<int>(sparseMatrix["ncol"]); // Number of distinct samples
//   int c = (n * n - n) / 2; // Number of unique combinations of 2 samples
//   
//   // Index the starting location of each sample in the vectors
//   // Assumes these vectors are sorted by j (samples), then i (taxa)
//   IntegerVector sPtr(n + 1);
//   for (int i = 0; i < samp.size(); i++) {
//     sPtr[samp[i]] = i + 1;
//   }
//   
//   // R style dist object, initialized to zero
//   NumericVector results(c);
//   results.attr("class")  = "dist";
//   results.attr("Size")   = n;
//   results.attr("Upper")  = false;
//   results.attr("Diag")   = false;
//   results.attr("Labels") = sampleNames;
//   
//   for (int i = 1; i <= c; i++) {
// 
//     if (i % mod != rem) { continue; }
// 
//     // Which two samples are we comparing?
//     int s1 = ceil(.5 * (-1 * pow(-8 * (i - 1) + 4 * pow(n, 2) - 4 * n - 7, .5) + 2 * n - 1) - 1) + 1;
//     int s2 = n - (s1 * (n - 1 - s1) + (s1 * (s1 + 1)) / 2) + i;
// 
//     // Convert from 1-based to 0-based indexing
//     s1--;
//     s2--;
// 
//     // What are these samples' first/last positions in the taxa and amt vectors?
//     int s1_first_idx = sPtr[s1];
//     int s1_last_idx  = sPtr[s1 + 1] - 1;
//     int s2_first_idx = sPtr[s2];
//     int s2_last_idx  = sPtr[s2 + 1] - 1;
// 
//     // The number of taxa found for sample 1 and 2
//     int s1_nTaxa = s1_last_idx - s1_first_idx + 1;
//     int s2_nTaxa = s2_last_idx - s2_first_idx + 1;
//     if (s1_nTaxa == 0 && s2_nTaxa == 0) {
//       results[i - 1] = 0.0;
//       continue;
//     }
// 
//     // The number of taxa common to samples 1 and 2
//     int common_nTaxa = 0;
//     int y = s2_first_idx;
//     for (int x = s1_first_idx; x <= s1_last_idx; x++) {
//       while (taxa[y] < taxa[x] && y <= s2_last_idx) { y++;                 }
//       if (y > s2_last_idx)                          { break;               }
//       if (taxa[y] == taxa[x])                       { common_nTaxa++; y++; }
//     }
// 
//     double result = (s1_nTaxa + s2_nTaxa - 2.0 * common_nTaxa) / (s1_nTaxa + s2_nTaxa);
//     results[i - 1] = result;
//     //printf("%i = [%i, %i] = %.3f\n", i, s1, s2, result);
//   }
// 
//   return results;
// }
// 
// 
// // You can include R code blocks in C++ files processed with sourceCpp
// // (useful for testing and development). The R code will be automatically 
// // run after the compilation.
// //
// 
// /*** R
// library(rbiom)
// library(slam)
// 
// infile <- system.file("extdata", "hmp50.biom", package = "rbiom")
// biom <- read.biom(infile)
// #biom <- select(biom, 1:3)
// 
// otus <- biom$counts
// ord  <- order(otus$j, otus$i)
// otus$i <- otus$i[ord]
// otus$j <- otus$j[ord]
// otus$v <- otus$v[ord]
// 
// g1 <- braycurtis_unweighted(2, 0, otus)
// g2 <- braycurtis_unweighted(2, 1, otus)
// g <- g1 + g2
// h <- vegan::vegdist(t(as.matrix(otus)), binary=TRUE)
// identical(as.vector(g), as.vector(h))
// */
