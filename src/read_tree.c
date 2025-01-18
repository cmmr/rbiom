#include <R.h>
#include <Rinternals.h>
#include <stdio.h>      // printf
#include <stdlib.h>     // strtod
#include <string.h>     // strlen, strncpy
#include <stdbool.h>    // bool


SEXP _read_tree(SEXP sexp_tree);
void readtree2(
    const char   *tree, 
    unsigned int  x1, 
    unsigned int  x2, 
    unsigned int  parent,
    unsigned int *eIdx, 
    unsigned int *nIdx, 
    unsigned int *lIdx, 
    int          *edge_mtx, 
    double       *eLen,
    unsigned int  nEdges,
    unsigned int  nLeafs,
    SEXP          nLab, 
    SEXP          lLab);
char* extractname(
    const char   *tree, 
    unsigned int  x1, 
    unsigned int  x2);



SEXP C_read_tree(SEXP sexp_tree) {
  
  const char *tree = CHAR(asChar(sexp_tree));
  
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
  SEXP retEdges       = PROTECT(allocMatrix(INTSXP,  nEdges, 2));
  SEXP retEdgeLengths = PROTECT(allocVector(REALSXP, nEdges));
  SEXP retNodeLabels  = PROTECT(allocVector(STRSXP,  nNodes));
  SEXP retLeafLabels  = PROTECT(allocVector(STRSXP,  nLeafs));
  
  // Track next open indices in output vectors using "global" counters
  unsigned int eIdx = 0;
  unsigned int nIdx = 0;
  unsigned int lIdx = 0;
  
  int    *edge_mtx = INTEGER(retEdges);
  double *eLen     = REAL(retEdgeLengths);
  
  
  // Start recursing at the highest level of parentheses; i.e. the whole tree
  readtree2(tree, x1, x2, 0, &eIdx, &nIdx, &lIdx, edge_mtx, eLen, nEdges, nLeafs, retNodeLabels, retLeafLabels);
  
  SEXP result = PROTECT(allocVector(VECSXP, 5));
  SET_VECTOR_ELT(result, 0, retEdges);
  SET_VECTOR_ELT(result, 1, ScalarInteger(nNodes));
  SET_VECTOR_ELT(result, 2, retLeafLabels);
  SET_VECTOR_ELT(result, 3, retEdgeLengths);
  SET_VECTOR_ELT(result, 4, retNodeLabels);
  
  
  UNPROTECT(5);
  
  return result;
}








void readtree2(
    const char   *tree, 
    unsigned int  x1, 
    unsigned int  x2, 
    unsigned int  parent,
    unsigned int *eIdx, 
    unsigned int *nIdx, 
    unsigned int *lIdx, 
    int          *edge_mtx, 
    double       *eLen,
    unsigned int  nEdges,
    unsigned int  nLeafs,
    SEXP          nLab, 
    SEXP          lLab) {
  
  
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
        char *junk_ptr;
        eLen[(*eIdx) - 1] = strtod(tree + i + 1, &junk_ptr);
      }
      
      x2 = i - 1;
    }
    
    // Text after the end-paren is the node name
    if (tree[i] == ')') {
      
      if (i < x2)
        SET_STRING_ELT(nLab, (*nIdx), mkChar(extractname(tree, i + 1, x2)));
      
      if ((*eIdx) > 0) {
        int eRow = (*eIdx) - 1;
        edge_mtx[eRow + 0]      = parent  + nLeafs;
        edge_mtx[eRow + nEdges] = (*nIdx) + nLeafs + 1;
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
      SET_STRING_ELT(lLab, (*lIdx), mkChar(extractname(tree, x1, x2)));
    
    if (*eIdx > 0) {
      int eRow = (*eIdx) - 1;
      edge_mtx[eRow + 0]      = parent + nLeafs;
      edge_mtx[eRow + nEdges] = (*lIdx) + 1;
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
      readtree2(tree, x1, i - 1, parent, eIdx, nIdx, lIdx, edge_mtx, eLen, nEdges, nLeafs, nLab, lLab);
      x1 = i + 1;
    }
    
  }
  readtree2(tree, x1, x2, parent, eIdx, nIdx, lIdx, edge_mtx, eLen, nEdges, nLeafs, nLab, lLab);
  
  return;
}



char* extractname(const char *tree, unsigned int x1, unsigned int x2) {
  
  bool quoted = tree[x1] == '\'' && tree[x2] == '\'';
  
  // Quoted Name ==> Strip off quote marks
  if (quoted) {
    x1++;
    x2--;
  }
  
  char* nodeName = (char*) malloc(x2 - x1 + 2);
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


