#include <Rcpp.h>
#include <stdio.h>      /* printf */
#include <stdlib.h>     /* atof */
#include <string.h>     /* strncpy */
using namespace Rcpp;


List rcpp_read_tree(const char*);
void readtree2(
    const char     *tree, 
    unsigned int    x1, 
    unsigned int    x2, 
    unsigned int    parent,
    unsigned int    eIdx, 
    unsigned int    nIdx, 
    unsigned int    lIdx, 
    NumericMatrix   edge, 
    NumericVector   eLen,
    CharacterVector nLab, 
    CharacterVector lLab);
char* extractname(
    const char   *tree, 
    unsigned int x1, 
    unsigned int x2);



// [[Rcpp::export]]
List rcpp_read_tree(const char* tree) {
  
  // Start and End positions of the newick string
  unsigned int x1 = 0; 
  unsigned int x2 = strlen(tree) - 1;
  
  // Determine how many nodes are in this tree
  unsigned int nNodes = 0;
  for (unsigned int i = x1; i <= x2; i++) {
    
    // Ignore special characters inside single quotes
    if (tree[i] == '\'') {
      do { i++; } while (tree[i] != '\'' && i <= x2);
      continue;
    }
    
    // Only consider the first tree in the file
    if (tree[i] == ';') {
      x2 = i - 1;
      break;
    }
    
    if (tree[i] == ',') nNodes++;
  }
  unsigned int nLeafs = nNodes + 1;
  unsigned int n = nNodes + nLeafs;
  
  
  // Convert the result into a phylo-compatible data structure for R
  NumericMatrix   retEdges(n - 1, 2);
  NumericVector   retEdgeLengths(n - 1);
  CharacterVector retLeafLabels(nLeafs);
  CharacterVector retNodeLabels(nNodes);
  
  // Start recursing at the highest level of parentheses; i.e. the whole tree
  readtree2(tree, x1, x2, 0, 0, 0, 0, retEdges, retEdgeLengths, retNodeLabels, retLeafLabels);
  
  
  List ret = List::create(
      Named("edge")        = retEdges,
      Named("Nnode")       = nNodes,
      Named("tip.label")   = retLeafLabels,
      Named("edge.length") = retEdgeLengths,
      Named("node.label")  = retNodeLabels
  );
  ret.attr("class") = "phylo";
  ret.attr("order") = "cladewise";
  
  return ret;
}








void readtree2(
    const char     *tree, 
    unsigned int    x1, 
    unsigned int    x2, 
    unsigned int    parent,
    unsigned int    eIdx, 
    unsigned int    nIdx, 
    unsigned int    lIdx, 
    NumericMatrix   edge, 
    NumericVector   eLen,
    CharacterVector nLab, 
    CharacterVector lLab) {
  
  unsigned int i;
  
  // Trim off whitespace from beginning and end of string section
  while ((tree[x1] == ' ' || tree[x1] == '\t') && x1 <= x2) x1++;
  while ((tree[x2] == ' ' || tree[x2] == '\t') && x1 <= x2) x2--;
  
  // Read backwards, extracting name and length if present
  for (i = x2; i >= x1; i--) {
    
    // Ignore special characters inside single quotes
    if (tree[i] == '\'') {
      do { i--; } while (tree[i] != '\'' && i >= x1);
      continue;
    }
    
    // Text after the colon is the branch length
    if (tree[i] == ':') {
      
      if (eIdx > 0) {
        char* strLength = new char[x2 - i + 1];
        strncpy(strLength, tree + i + 1, x2 - i);
        strLength[x2 - i]  = '\0';
        eLen[eIdx - 1] = atof(strLength);
      }
      
      x2 = i - 1;
    }
    
    // Text after the end-paren is the node name
    if (tree[i] == ')') {
      
      if (i < x2)
        nLab[nIdx] = extractname(tree, i + 1, x2);
      
      if (eIdx > 0) {  
        edge(eIdx - 1, 0) = parent + lLab.size() + 1;
        edge(eIdx - 1, 1) = nIdx   + lLab.size() + 1;
      }
      
      // Peel off parentheses
      x1++;
      x2 = i - 1;
      
      break;
    }
  }
  
  
  // No parens means we're at a leaf
  if (i <= x1) {
    
    if (x1 < x2)
      lLab[lIdx] = extractname(tree, x1, x2);
    
    if (eIdx > 0) {
      edge(eIdx - 1, 0) = parent + lLab.size() + 1;
      edge(eIdx - 1, 1) = lIdx + 1;
    }
    
    return;
  }
  
  // Find the comma for the current block (pivot)
  // Also track how many nodes are in the first half (n)
  unsigned int level = 0;
  unsigned int pivot = 0;
  unsigned int n     = 0;
  for (i = x1; i <= x2; i++) {
    
    // Ignore special characters inside single quotes
    if (tree[i] == '\'') {
      do { i++; } while (tree[i] != '\'' && i <= x1);
      continue;
    }
    
    if (tree[i] == '(') {
      level++;
    } else if (tree[i] == ')') {
      level--;
    }
    
    if (tree[i] == ',') {
      if (level == 0) {
        pivot = i;
        break;
      }
      n++;
    }
  }
  unsigned int nl = n + n + 1; // nodes and leafs in first half
  
  
  // Recurse into the two subtrees
  //        tree  x1         x2  parent
  readtree2(tree, x1, pivot - 1, nIdx, eIdx + 1,      nIdx + 1,     lIdx + 0,     edge, eLen, nLab, lLab);
  readtree2(tree, pivot + 1, x2, nIdx, eIdx + 1 + nl, nIdx + 1 + n, lIdx + 1 + n, edge, eLen, nLab, lLab);
  
  return;
}



char* extractname(const char *tree, unsigned int x1, unsigned int x2) {
  
  bool quoted = tree[x1] == '\'' && tree[x2] == '\'';

  // Quoted Name ==> Strip off quote marks
  if (quoted) {
    x1++;
    x2--;
  }
  
  char* nodeName = new char[x2 - x1 + 2];
  strncpy(nodeName, tree + x1, x2 - x1 + 1);
  nodeName[x2 - x1 + 1] = '\0';
  
  // Unquoted Name ==> Replace underscores with spaces
  if (!quoted) {
    for (unsigned int j = 0; j <= x2 - x1; j++) {
      if (nodeName[j] == '_') nodeName[j] = ' ';
    }
  }
      
  return nodeName;
}


