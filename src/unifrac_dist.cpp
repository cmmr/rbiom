// #include <Rcpp.h>
// #include <math.h>
// using namespace Rcpp;
// 
// 
// // [[Rcpp::export]]
// NumericVector timesTwo(int mod, int rem, List sparseMatrix) {
// 
//   IntegerVector samp = as<IntegerVector>(sparseMatrix["j"]);
//   IntegerVector taxa = as<IntegerVector>(sparseMatrix["i"]);
//   NumericVector amt  = as<NumericVector>(sparseMatrix["v"]);
// 
//   List dimnames               = as<List>(sparseMatrix["dimnames"]);
//   CharacterVector sampleNames = as<CharacterVector>(dimnames[1]);
// 
//   int n = as<int>(sparseMatrix["ncol"]); // Number of distinct samples
// 
//   // Index the starting location of each sample in the vectors
//   // Assumes these vectors are sorted by j (samples), then i (taxa)
//   IntegerVector sPtr(n + 1);
//   for (int i = 0; i < samp.size(); i++) {
//     sPtr[samp[i]] = i + 1;
//   }
// 
//   // R style dist object, initialized to zero
//   NumericVector results((n * n - n) / 2);
//   results.attr("class")  = "dist";
//   results.attr("Size")   = n;
//   results.attr("Upper")  = false;
//   results.attr("Diag")   = false;
//   results.attr("Labels") = sampleNames;
// 
//   int i = -1;
//   for (int s1 = 0; s1 < n; s1++) {
//     for (int s2 = s1 + 1; s2 < n; s2++) {
// 
//       if (++i % mod != rem) { continue; }
// 
//       // What are these samples' first/last positions in the taxa and amt vectors?
//       int s1_first_idx = sPtr[s1];
//       int s1_last_idx  = sPtr[s1 + 1] - 1;
//       int s2_first_idx = sPtr[s2];
//       int s2_last_idx  = sPtr[s2 + 1] - 1;
// 
// 
//       x  <- sample2branchwts[[sampleNames[[rowIdx]]]]
//       y  <- sample2branchwts[[sampleNames[[colIdx]]]]
// 
//       if (weighted == FALSE) {
//         res <- sum(x[!names(x) %in% names(y)], y[!names(y) %in% names(x)])
//         res <- res / sum(x, y[!names(y) %in% names(x)])
// 
//       } else {
// 
//         z <- intersect(names(x), names(y))
//         z <- abs(x[z] - y[z])
//         x <- x[!names(x) %in% names(z)]
//         y <- y[!names(y) %in% names(z)]
// 
//         res <- sum(x, y, z)
//         res <- res / sum(otuDepths * otus[,rowIdx], otuDepths * otus[,colIdx])
//       }
// 
//       double result = sqrt(amt_diff);
//       results[i] = result;
//     }
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
// timesTwo(42)
// */
