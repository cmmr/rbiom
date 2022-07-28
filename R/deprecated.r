
# For consistency's sake, use underscores in function names instead of periods.

alpha.div        <- function (...) { .Deprecated("adiv_table");    adiv_table(...)       }
as.percent       <- function (...) { .Deprecated("as_percent");    as_percent(...)       }
has.phylogeny    <- function (...) { .Deprecated("has_phylogeny"); has_phylogeny(...)    }
has.sequences    <- function (...) { .Deprecated("has_sequences"); has_sequences(...)    }
init.cache       <- function (...) { .Deprecated("init_cache");    init_cache(...)       }
is.rarefied      <- function (...) { .Deprecated("is_rarefied");   is_rarefied(...)      }
read.biom        <- function (...) { .Deprecated("read_biom");     read_biom(...)        }
read.fasta       <- function (...) { .Deprecated("read_fasta");    read_fasta(...)       }
read.tree        <- function (...) { .Deprecated("read_tree");     read_tree(...)        }
sample.names     <- function (...) { .Deprecated("sample_names");  sample_names(...)     }
`sample.names<-` <- function (...) { .Deprecated("sample_names");  `sample_names<-`(...) }
sample.sums      <- function (...) { .Deprecated("sample_sums");   sample_sums(...)      }
stats.table      <- function (...) { .Deprecated("stats_table");   stats_table(...)      }
taxa.means       <- function (...) { .Deprecated("taxa_means");    taxa_means(...)       }
taxa.names       <- function (...) { .Deprecated("taxa_names");    taxa_names(...)       }
`taxa.names<-`   <- function (...) { .Deprecated("taxa_names");    `taxa_names<-`(...)   }
taxa.ranks       <- function (...) { .Deprecated("taxa_ranks");    taxa_ranks(...)       }
`taxa.ranks<-`   <- function (...) { .Deprecated("taxa_ranks");    `taxa_ranks<-`(...)   }
taxa.rollup      <- function (...) { .Deprecated("taxa_rollup");   taxa_rollup(...)      }
taxa.sums        <- function (...) { .Deprecated("taxa_sums");     taxa_sums(...)        }
top.taxa         <- function (...) { .Deprecated("top_taxa");      top_taxa(...)         }
write.biom       <- function (...) { .Deprecated("write_biom");    write_biom(...)       }
write.fasta      <- function (...) { .Deprecated("write_fasta");   write_fasta(...)      }
write.tree       <- function (...) { .Deprecated("write_tree");    write_tree(...)       }
write.xlsx       <- function (...) { .Deprecated("write_xlsx");    write_xlsx(...)       }




# Depending on its arguments, beta.div() returned either a dist object or a 
# data.frame. Now it's split into two functions - bdiv_dist() which always 
# returns a dist, and bdiv_table() which always returns a data.frame.

beta.div <- function (
  biom, method, weighted = TRUE, tree = NULL, long = FALSE,
  md = FALSE, safe = FALSE, stat.by = NULL, seed = 0, perms = 999 ) {
  
  msg <- paste0(
    "The rbiom function beta.div() is deprecated.\n",
    "Use bdiv_%s() for generating a %s instead." )
  
  if (isTRUE(long) || !isFALSE(md)) {
    .Deprecated("bdiv_table", msg = sprintf(msg, "table", "data.frame"))
    bdiv_table(biom, method, weighted, tree, md, safe, stat.by, seed, perms)
    
  } else {
    .Deprecated("bdiv_dist", msg = sprintf(msg, "dist", "distance matrix"))
    bdiv_dist(biom, method, weighted, tree, stat.by, seed, perms)
  }
}
  
