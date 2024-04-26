#' Human Microbiome Project - demo dataset (n = 50)
#'
#' @format An rbiom object with 50 samples. 
#'         Includes metadata, taxonomy, phylogeny, and sequences.
#' \itemize{
#'   \item{Sex - }{Male or Female}
#'   \item{Body Site - }{Anterior nares, Buccal mucosa, Mid vagina, Saliva, or Stool}
#'   \item{Age - }{21 - 40}
#'   \item{BMI - }{19 - 32}
#' }
#' @source \url{https://portal.hmpdacc.org}
"hmp50"


#' Global Enteric Multicenter Study (n = 1,006)
#'
#' @format An rbiom object with 1,006 samples. 
#'         Includes metadata and taxonomy.
#' \itemize{
#'   \item{Diarrhea - }{Case or Control}
#'   \item{Sex - }{Male or Female}
#'   \item{Age (months) - }{0 - 58}
#'   \item{Country - }{Bangladesh, Gambia, Kenya, or Mali}
#' }
#' @source \url{https://doi.org/10.1186/gb-2014-15-6-r76} and \url{https://doi.org/10.1093/nar/gkx1027}
"gems"


#' Longitudinal Stool Samples from Infants (n = 2,684)
#'
#' @format An rbiom object with 2,684 samples. 
#'         Includes metadata and taxonomy.
#' \itemize{
#'   \item{Subject ID - }{ID1, ID2, ..., ID12}
#'   \item{Sex - }{Male or Female}
#'   \item{Age (days) - }{1 - 266}
#'   \item{Child's diet - }{"Breast milk", "Breast milk and formula", or "Formula"}
#'   \item{Sample collection - }{"Frozen upon collection" or "Stored in alcohol"}
#'   \item{Antibiotic exposure - }{Yes or No}
#'   \item{Antifungal exposure - }{Yes or No}
#'   \item{Delivery mode - }{Cesarean or Vaginal}
#'   \item{Solid food introduced (Age) - }{116 - 247}
#' }
#' @source \url{https://www.nature.com/articles/s41467-018-04641-7} and \url{https://doi.org/10.1038/s41467-017-01973-8}
"babies"
