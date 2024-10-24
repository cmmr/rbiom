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
    unsigned int    *eIdx, 
    unsigned int    *nIdx, 
    unsigned int    *lIdx, 
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
  unsigned int nLeafs = 1;
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
    
    if (tree[i] == '(') nNodes++;
    if (tree[i] == ',') nLeafs++;
  }
  unsigned int nEdges = nNodes + nLeafs - 1;
  
  
  // Convert the result into a phylo-compatible data structure for R
  NumericMatrix   retEdges(nEdges, 2);
  NumericVector   retEdgeLengths(nEdges);
  CharacterVector retNodeLabels(nNodes);
  CharacterVector retLeafLabels(nLeafs);
  
  // Track next open indices in output vectors using "global" counters
  unsigned int eIdx = 0;
  unsigned int nIdx = 0;
  unsigned int lIdx = 0;
  
  
  // Start recursing at the highest level of parentheses; i.e. the whole tree
  readtree2(tree, x1, x2, 0, &eIdx, &nIdx, &lIdx, retEdges, retEdgeLengths, retNodeLabels, retLeafLabels);
  
  
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
    unsigned int    *eIdx, 
    unsigned int    *nIdx, 
    unsigned int    *lIdx, 
    NumericMatrix   edge, 
    NumericVector   eLen,
    CharacterVector nLab, 
    CharacterVector lLab) {
  
  
  // if ((*eIdx) % 1000 == 0) {
  //   Rcpp::checkUserInterrupt();
  // }
  
  
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
      
      if ((*eIdx) > 0) {
        char* strLength = new char[x2 - i + 1];
        strncpy(strLength, tree + i + 1, x2 - i);
        strLength[x2 - i]  = '\0';
        eLen[(*eIdx) - 1] = atof(strLength);
      }
      
      x2 = i - 1;
    }
    
    // Text after the end-paren is the node name
    if (tree[i] == ')') {
      
      if (i < x2)
        nLab[(*nIdx)] = extractname(tree, i + 1, x2);
      
      if ((*eIdx) > 0) {  
        edge((*eIdx) - 1, 0) = parent  + lLab.size();
        edge((*eIdx) - 1, 1) = (*nIdx) + lLab.size() + 1;
      }
      (*nIdx)++;
      (*eIdx)++;
      
      // Peel off parentheses
      x1++;
      x2 = i - 1;
      
      break;
    }
  }
  
  
  // No parens means we're at a leaf
  if (i <= x1) {
    
    if (x1 <= x2)
      lLab[(*lIdx)] = extractname(tree, x1, x2);
    
    if (*eIdx > 0) {
      edge((*eIdx) - 1, 0) = parent + lLab.size();
      edge((*eIdx) - 1, 1) = (*lIdx) + 1;
    }
    (*lIdx)++;
    (*eIdx)++;
    
    return;
  }
  
  
  // Recurse into each of the subtrees
  parent = (*nIdx);
  unsigned int level  = 0;
  for (i = x1; i <= x2; i++) {
    
    // Ignore special characters inside single quotes
    if (tree[i] == '\'') {
      do { i++; } while (tree[i] != '\'' && i <= x2);
      continue;
    }
    
    // Find the other end of the current subclade
    if (tree[i] == '(') {
      level++;
    } else if (tree[i] == ')') {
      level--;
    } else if (tree[i] == ',' && level == 0) {
      readtree2(tree, x1, i - 1, parent, eIdx, nIdx, lIdx, edge, eLen, nLab, lLab);
      x1 = i + 1;
    }
    
  }
  readtree2(tree, x1, x2, parent, eIdx, nIdx, lIdx, edge, eLen, nLab, lLab);
  
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


