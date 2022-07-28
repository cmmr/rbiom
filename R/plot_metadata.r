
#________________________________________________________
# Subset a data.frame of metadata to meet the
# user's criteria for colors, facets, etc.
#________________________________________________________

subset_by_params <- function (df, params) {
  
  groupvars <- c()
  
  #________________________________________________________
  # Drop NA and unused values in columns of interest
  #________________________________________________________
  for (i in c("xval", "color", "shape", "facet", "pattern", "label", "order")) {
    
    colname <- params[[paste0(i, ".by")]] %||% next
    colvals <- params[[paste0(i, "s")]]
    
    colname <- sub("^[!=]=", "", as.vector(colname))
    
    # 'facet' might specify multiple columns, e.g. c(".rank", "Body Site")
    for (cname in colname) {
      
      if (!hasName(df, cname)) next
      groupvars %<>% c(cname)
      
      if (is.character(df[[cname]]))
        df[[cname]] %<>% as.factor()
      
      if (is.factor(df[[cname]]) && !is.null(colvals)) {
        
        cvals <- colvals                          # c('Saliva', 'Stool')
        if (is.list(cvals)) {
          alt_cname <- sub("^\\.(.)", "\\U\\1", cname, perl = TRUE) # match '.rank'= or 'Rank'=
          cvals     <- cvals[[cname]] %||%        # list('Body Site' = c('Saliva', 'Stool'))
            cvals[[alt_cname]] %||%               # list('Rank' = c('Phylum', 'Genus'))
            cvals[[which(cname == colname)]] %||% # list(c('Phylum', 'Genus'), c('Saliva', 'Stool'))
            cvals[[1]]                            # list(c('Saliva', 'Stool'))
          remove("alt_cname")
        }
        
        new_levels <- names(cvals) %||% unname(cvals)
        if (all(new_levels %in% levels(df[[cname]])))
          df[[cname]] %<>% factor(levels = new_levels)
      }
      
      df <- df[!is.na(df[[cname]]), , drop = FALSE]
      
    }
  }
  
  attr(df, 'groupvars') <- groupvars
  
  return (df)
}





#________________________________________________________
# Subset a data.frame, removing NAs and 
# refactoring according to arguments.
#________________________________________________________

metadata_filters <- function (ggdata, biom = NULL, params = NULL) {
  
  params %||=% attr(ggdata, 'params', exact = TRUE)
  biom   %||=% params[['biom']]
  
  
  # Re-attach after subsetting
  attribs <- attributes(ggdata)
  ggdata  <- subset_by_params(ggdata, params) # also sets groupvars attr
  
  
  #________________________________________________________
  # Drop columns not referenced above
  #________________________________________________________
  groupvars       <- attr(ggdata, 'groupvars', exact = TRUE)
  non_ggdata_cols <- setdiff(colnames(ggdata), colnames(metadata(biom)))
  ggdata <- ggdata[, unique(c(non_ggdata_cols, groupvars)), drop=FALSE]
  
  remove("non_ggdata_cols")
  
  
  #________________________________________________________
  # Add a group col to ensure groups are always separated
  #________________________________________________________
  if (length(groupvars) == 0) {
    ggdata[[".group"]] <- "all"
    
  } else {
    df <- ggdata[, unique(groupvars), drop=FALSE]
    ggdata[[".group"]] <- unname(apply(df, 1L, paste, collapse="|"))
    
    remove("df")
  }
  
  
  for (i in setdiff(names(attribs), names(attributes(data.frame()))))
    attr(ggdata, i) <- attribs[[i]]
  
  
  return (ggdata)
}




#________________________________________________________
# Handle configuration of ggplot2 layers with
# respect to metadata columns.
#________________________________________________________

metadata_layers <- function (layers, ...) {
  
  dots   <- list(...)
  ggdata <- attr(layers, 'data',   exact = TRUE)
  params <- attr(layers, 'params', exact = TRUE)
  
  
  #________________________________________________________
  # Colors, shapes, and patterns
  #________________________________________________________
  for (layer in c("pattern", "color", "shape")) {
    
    colname <- params[[paste0(layer, ".by")]] %||% next
    usrvals <- params[[paste0(layer, "s")]]
    
    colvals <- levels(ggdata[[colname]])
    palette <- get_palette(layer, length(colvals))
    
    # browser()
    
    if (is.character(usrvals)) {
      if (is.null(names(usrvals))) {
        if (!all(colvals %in% usrvals))
          palette <- usrvals # colors = c('#00B9EB', '#ED5F16')
      } else {
        palette <- unname(usrvals[colvals])  # colors = c(a = '#00B9EB', b = '#ED5F16')
      }
    }
    
    
    # palette %<>% setNames(colvals)
    setLayer(values = palette)
    
    if (layer == "color")
      setLayer(layer = "fill", values = palette)
  }
  
  
  return (layers)
}



