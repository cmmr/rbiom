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



// [[Rcpp::export]]
List rcpp_read_tree(const char* tree) {
  
  // Start and End positions of the newick string
  unsigned int x1 = 0; 
  unsigned int x2 = strlen(tree) - 1;
  
  // Determine how many nodes are in this tree
  unsigned int nNodes = 0;
  for (unsigned int i = x1; i <= x2; i++) {
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
  
  // Read backwards, extracting name and length if present
  for (unsigned int i = x2; i >= x1; i--) {
    
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
      char* nodeName = new char[x2 - i + 1];
      strncpy(nodeName, tree + i + 1, x2 - i);
      nodeName[x2 - i] = '\0';
      
      nLab[nIdx] = nodeName;
      
      if (eIdx > 0) {  
        edge(eIdx - 1, 0) = parent + lLab.size() + 1;
        edge(eIdx - 1, 1) = nIdx   + lLab.size() + 1;
      }
      
      // Peel off parentheses
      x1++;
      x2 = i - 1;
      
      break;
    }
    
    // No parens means we're at a leaf
    if (i == x1) {
      
      char* leafName = new char[x2 - i + 2];
      strncpy(leafName, tree + i, x2 - i + 1);
      leafName[x2 - i + 1]  = '\0';
      
      lLab[lIdx] = leafName;
      
      if (eIdx > 0) {
        edge(eIdx - 1, 0) = parent + lLab.size() + 1;
        edge(eIdx - 1, 1) = lIdx + 1;
      }
      
      return;
    }
  }
  
  
  // Find the comma for the current block (pivot)
  // Also track how many nodes are in the first half (n)
  unsigned int level = 0;
  unsigned int pivot = 0;
  unsigned int n     = 0;
  for (unsigned int i = x1; i <= x2; i++) {
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


/*** R
x <- readtree("((walrus:10,seal:12)sea:8,((monkey:100,tiger:50)jungle:20,ferret:20)land:2)all:1;")
*/




