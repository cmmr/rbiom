

#________________________________________________________
#' Convert all *.by params into long form.
#'
#' Examples values for color.by:
#' 
#' "Body Site"
#' c("Body Site" = "okabe")
#' list("Body Site")
#' list("Body Site" = "okabe")
#' list("Body Site" = list(values = "okabe", name = "Swab Location"))
#' list("Body Site" = list(values = c("red", "green", "blue"), name = "Swab Location"))
#' list("Body Site" = list(values = c(Saliva = "red", Stool = "green"), name = "Swab Location"))
#' 
#' 
#' c("Body Site", "Sex")
#' c("Body Site" = "okabe", "Sex")
#' c("Body Site" = "okabe", "Sex" = "muted")
#' list("Body Site", "Sex")
#' 
metadata_params <- function (params, contraints = list()) {
  
  
  biom   <- params[['biom']]
  md     <- if (is(biom, "BIOM")) metadata(biom) else biom
  new_md <- md
  
  by_params <- c(
    "x", "color.by", "shape.by", "facet.by", 
    "pattern.by", "label.by", "order.by", "stat.by", "limit.by" )
  group_params <- c(
    "x", "color.by", "shape.by", "facet.by", 
    "pattern.by", "stat.by" )
  
  
  
  get_col_type <- function (col_name) {
    
    col_name <- sub("^([!=]=|\\-)", "", col_name)
    
    if (!hasName(md, col_name))    return (NULL)
    if (is.factor(md[[col_name]])) return ("cat")
    
    if (is.character(md[[col_name]])) {
      md[[col_name]]     <<- factor(md[[col_name]])
      new_md[[col_name]] <<- factor(new_md[[col_name]])
      return ("cat")
    }
    
    return ("num")
  }
  
  
  
  #________________________________________________________
  # Tracks which metadata columns are used.
  #________________________________________________________
  md_params <- intersect(by_params, names(params))
  md_cols   <- c() # Metadata columns needed in ggdata.
  md_cmps   <- c() # For beta div, == or != prefixed cols.
  md_revs   <- c() # For order.by, reverse if '-' prefix.
  group_by  <- c() # ggdata cols to use for ggplot group.
  
  
  
  #________________________________________________________
  # Standardize the data structure for metadata specs.
  #________________________________________________________
  for (param in md_params) {
    
    results <- list()
    
    for (i in seq_along(params[[param]])) {
      
      result   <- list()
      col_name <- names(params[[param]][i])
      col_spec <- params[[param]][[i]]
      col_type <- NULL
      
      
      if (identical(col_name, "") || is_null(col_name)) {
        
        col_name <- col_spec
        col_type <- get_col_type(col_name)
        
      } else {
        
        col_type <- get_col_type(col_name)
        
        if (is_list(col_spec))           { result             <- col_spec
        } else if (is_palette(col_spec)) { result[['colors']] <- col_spec
        } else if (col_type == "num")    { result[['range']]  <- col_spec
        } else                           { result[['values']] <- col_spec }
      }
      
      
      stopifnot(is_scalar_character(col_name))
      
      
      #________________________________________________________
      # Reverse sort any columns prefixed with a '-'.
      #________________________________________________________
      if (identical(param, 'order.by')) {
        prefixed <- startsWith(col_name, "-")
        col_name <- sub("^\\-", "", col_name)
        md_revs[[col_name]] <- prefixed
        remove("prefixed")
      }
      
      
      #________________________________________________________
      # For bdiv_boxplot, hande within/between/all comparisons.
      #________________________________________________________
      if (isTRUE(params[['.cmp']])) {
        prefix   <- substr(col_name, 1, 2)
        col_name <- sub("^[!=]=", "", col_name)
        if (!hasName(md_cmps, col_name)) md_cmps[[col_name]] <- ""
        if (identical(prefix, "=="))     md_cmps[[col_name]] <- "=="
        if (identical(prefix, "!="))     md_cmps[[col_name]] <- "!="
        remove("prefix")
      }
      
      
      
      if (hasName(result, 'range')) {
        
        if (!identical(col_type, "num"))
          stop ("Can't apply `range` to non-numeric '", col_name, "' column.")
        
        if (!is.numeric(result[['range']]))
          stop ("Argument for '", col_name, "' `range` must be numeric.")
        
        
        result[['range']] <- range(result[['range']])
        new_md[which(md[[col_name]] < result[['range']][[1]]), col_name] <- NA
        new_md[which(md[[col_name]] > result[['range']][[2]]), col_name] <- NA
        
        
      } else if (identical(col_type, "cat") && hasName(result, 'values')) {
        
        values <- result[['values']]
        
        if (!identical(col_type, "cat"))
          stop ("Can't filter by `values` on non-categorical '", col_name, "' column.")
        
        if (!is.character(values) && !identical(param, "shape.by"))
          stop ("Values for '", col_name, "' must be a character vector.")
        
        
        if (is_null(names(values))) { # c("Saliva", "Stool")
          values <- intersect(values, levels(md[[col_name]]))
          result[['values']] <- NULL
          
        } else { # c("Saliva" = "blue", "Stool" = "orange")
          values <- intersect(names(values), levels(md[[col_name]]))
          result[['values']] <- result[['values']][values]
        }
        
        new_md[[col_name]] %<>% factor(levels=values)
        remove("values")
      }
      
      
      
      if (param %in% group_params)   group_by %<>% c(col_name)
      if (hasName(new_md, col_name)) md_cols  %<>% c(col_name)
      
      results[[col_name]] <- result
    }
    
    
    #________________________________________________________
    # Enforce validation constraints.
    #________________________________________________________
    validate_arg(
      vals     = names(results), 
      biom     = new_md, 
      arg      = param, 
      mode     = 'meta', 
      n        = contraints[[param]][['n']], 
      col_type = contraints[[param]][['col_type']], 
      default  = contraints[[param]][['default']] )
    
    
    params[[param]] <- results
  }
  
  
  
  #________________________________________________________
  # Which columns to concatenate and use for ggplot group.
  #________________________________________________________
  params[['.group.by']] <- unique(group_by)
  
  
  #________________________________________________________
  # Other column and column attribute tracking.
  #________________________________________________________
  params[['.md.cols']]  <- unique(md_cols)
  params[['.md.cmps']]  <- md_cmps
  
  
  
  #________________________________________________________
  # Drop rows with NA's from the metadata table.
  #________________________________________________________
  new_md <- new_md[,unique(md_cols),drop=FALSE]
  new_md <- new_md[complete.cases(new_md),,drop=FALSE]
  
  
  
  #________________________________________________________
  # Update all the *.by params to ensure subsets match.
  #________________________________________________________
  for (param in md_params) {
    
    for (i in seq_along(params[[param]])) {
      
      col_name <- names(params[[param]][i])
      result   <- params[[param]][[i]]
      values   <- result[['values']]
      col_type <- get_col_type(col_name)
      
      
      if (identical(col_type, "cat")) {
        
        col_vals <- levels(new_md[[col_name]])
        n        <- length(col_vals)
        
        
        # Will need additional colors/shapes/patterns for bdiv.
        if (isTRUE(params[['.cmp']]) && n > 1)
          n <- as.integer(local({
            cmp <- md_cmps[[col_name]]
            if (identical(cmp, "==")) return (n)
            if (identical(cmp, "!=")) return ((n * (n - 1) / 2))
            if (identical(cmp, ""))   return ((n * (n - 1) / 2) + n)
          }))
        
        
        if (is_null(values))
          values <- switch(param,
            "color.by"   = get_n_colors(n, result[['colors']]),
            "shape.by"   = get_n_shapes(n),
            "pattern.by" = get_n_patterns(n) )
        
        
        if (!is_null(names(values)))
          values <- values[col_vals]
        
        
      } else if (identical(col_type, "num")) {
        if (hasName(result, 'colors')) values <- result[['colors']]
        if (is_palette(values))        values %<>% get_palette()
      }
      
      result[['colors']]   <- NULL
      result[['values']]   <- values
      params[[param]][[i]] <- result
    }
  }
  
  
  
  #________________________________________________________
  # Simplify all *.by params except color/shape/pattern.
  #________________________________________________________
  for (param in md_params) {
    
    if (length(params[[param]]) == 0) {
      params[param] <- list(NULL)
      
    } else if (!param %in% c("color.by", "shape.by", "pattern.by")) {
      params[[param]] %<>% names()
    }
  }
  
  
  
  #________________________________________________________
  # Re-order the samples according to metadata.
  #________________________________________________________
  for (i in rev(params[['order.by']]))
    new_md <- new_md[order(new_md[[i]], decreasing = md_revs[[i]]),,drop=FALSE]
  
  
  
  if (is(biom, 'BIOM')) {
    biom             <- select(biom, rownames(new_md))
    metadata(biom)   <- new_md
    params[['biom']] <- biom
    
  } else {
    params[['biom']] <- new_md
  }
  
  
  return (params)
}




#________________________________________________________
# Add a group col to ensure groups are always separated.
#________________________________________________________
metadata_group <- function (ggdata, params = NULL) {
  
  params %||=% attr(ggdata, 'params', exact = TRUE)
  groupvars <- params[['.group.by']]
  
  if (length(groupvars) == 0) {
    ggdata[[".group"]] <- "all"
    
  } else {
    df <- ggdata[,groupvars,drop=FALSE]
    ggdata[[".group"]] <- unname(apply(df, 1L, paste, collapse="|"))
  }
  
  
  return (ggdata)
}



# #________________________________________________________
# # Subset a data.frame of metadata to meet the
# # user's criteria for colors, facets, etc.
# #________________________________________________________
# 
# subset_by_params <- function (df, params) {
#   
#   groupvars <- c()
#   
#   #________________________________________________________
#   # Drop NA and unused values in columns of interest
#   #________________________________________________________
#   for (i in c("xval", "color", "shape", "facet", "pattern", "label", "order", "stat")) {
#     
#     colname <- params[[paste0(i, ".by")]] %||% next
#     colvals <- params[[paste0(i, "s")]]
#     
#     colname <- sub("^[!=]=", "", as.vector(colname))
#     
#     # 'facet' might specify multiple columns, e.g. c(".rank", "Body Site")
#     for (cname in colname) {
#       
#       if (!hasName(df, cname)) next
#       groupvars %<>% c(cname)
#       
#       if (is.character(df[[cname]]))
#         df[[cname]] %<>% as.factor()
#       
#       if (is.factor(df[[cname]]) && !is_null(colvals)) {
#         
#         cvals <- colvals                          # c('Saliva', 'Stool')
#         if (is_list(cvals)) {
#           alt_cname <- sub("^\\.(.)", "\\U\\1", cname, perl = TRUE) # match '.rank'= or 'Rank'=
#           cvals     <- cvals[[cname]] %||%        # list('Body Site' = c('Saliva', 'Stool'))
#             cvals[[alt_cname]] %||%               # list('Rank' = c('Phylum', 'Genus'))
#             cvals[[which(cname == colname)]] %||% # list(c('Phylum', 'Genus'), c('Saliva', 'Stool'))
#             cvals[[1]]                            # list(c('Saliva', 'Stool'))
#           remove("alt_cname")
#         }
#         
#         new_levels <- names(cvals) %||% unname(cvals)
#         if (all(new_levels %in% levels(df[[cname]])))
#           df[[cname]] %<>% factor(levels = new_levels)
#       }
#       
#       df <- df[!is.na(df[[cname]]), , drop = FALSE]
#       
#     }
#   }
#   
#   attr(df, 'groupvars') <- groupvars
#   
#   return (df)
# }





# #________________________________________________________
# # Subset a data.frame, removing NAs and 
# # refactoring according to arguments.
# #________________________________________________________
# 
# metadata_filters <- function (ggdata, biom = NULL, params = NULL) {
#   
#   params %||=% attr(ggdata, 'params', exact = TRUE)
#   biom   %||=% params[['biom']]
#   
#   
#   # Re-attach after subsetting
#   attribs <- attributes(ggdata)
#   ggdata  <- subset_by_params(ggdata, params) # also sets groupvars attr
#   
#   
#   #________________________________________________________
#   # Drop columns not referenced above
#   #________________________________________________________
#   groupvars       <- attr(ggdata, 'groupvars', exact = TRUE)
#   non_ggdata_cols <- setdiff(colnames(ggdata), colnames(metadata(biom)))
#   ggdata <- ggdata[, unique(c(non_ggdata_cols, groupvars)), drop=FALSE]
#   
#   remove("non_ggdata_cols")
#   
#   
#   #________________________________________________________
#   # Add a group col to ensure groups are always separated
#   #________________________________________________________
#   if (length(groupvars) == 0) {
#     ggdata[[".group"]] <- "all"
#     
#   } else {
#     df <- ggdata[, unique(groupvars), drop=FALSE]
#     ggdata[[".group"]] <- unname(apply(df, 1L, paste, collapse="|"))
#     
#     remove("df")
#   }
#   
#   
#   for (i in setdiff(names(attribs), names(attributes(data.frame()))))
#     attr(ggdata, i) <- attribs[[i]]
#   
#   
#   return (ggdata)
# }




# #________________________________________________________
# # Handle configuration of ggplot2 layers with
# # respect to metadata columns.
# #________________________________________________________
# 
# metadata_layers <- function (layers, ...) {
#   
#   dots   <- list(...)
#   ggdata <- attr(layers, 'data',   exact = TRUE)
#   params <- attr(layers, 'params', exact = TRUE)
#   
#   
#   #________________________________________________________
#   # Colors, shapes, and patterns
#   #________________________________________________________
#   for (layer in c("pattern", "color", "shape")) {
#     
#     colname <- params[[paste0(layer, ".by")]] %||% next
#     usrvals <- params[[paste0(layer, "s")]]
#     
#     colvals <- levels(ggdata[[colname]])
#     n       <- length(colvals)
#     palette <- get_palette(layer, n)
#     
#     
#     if (!is_null(usrvals)) {
#       
#       if (length(usrvals) == 1 && usrvals %in% names(PALETTES)) {
#         palette <- color_palette(pal = usrvals, n = n) # colors = 'vibrant'
#         
#       } else if (is_null(names(usrvals))) {
#         if (!all(colvals %in% usrvals))
#           palette <- usrvals # colors = c('#00B9EB', '#ED5F16')
#         
#       } else {
#         palette <- unname(usrvals[colvals])  # colors = c(a = '#00B9EB', b = '#ED5F16')
#       }
#     }
#     
#     
#     # palette %<>% setNames(colvals)
#     setLayer(values = palette)
#     
#     if (layer == "color")
#       setLayer(layer = "fill", values = palette)
#   }
#   
#   
#   return (layers)
# }



