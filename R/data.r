#' Human Microbiome Project - demo dataset (n = 50)
#'
#' @format An rbiom object with 50 samples. 
#'         Includes metadata, taxonomy, phylogeny, and sequences.
#' \describe{
#'   \item{Sex - }{Male or Female}
#'   \item{Body Site - }{Anterior nares, Buccal mucosa, Mid vagina, Saliva, or Stool}
#'   \item{Age - }{21 - 40}
#'   \item{BMI - }{19 - 32}
#' }
#' @source \doi{10.1101/gr.096651.109}
#' @family Built-In Datasets
#' @usage hmp50
#' @name hmp50
#' @examples
#' hmp50
#' hmp50$metadata
#' @export
delayedAssign(
  x     = 'hmp50', 
  value = local({
    f <- system.file(package = 'rbiom', 'extdata', 'hmp50.bz2')
    if (nzchar(f)) rbiom::as_rbiom(f) else NULL
  })
)


#' Global Enteric Multicenter Study (n = 1,006)
#'
#' @format An rbiom object with 1,006 samples. 
#'         Includes metadata and taxonomy.
#' \describe{
#'   \item{diarrhea - }{Case or Control}
#'   \item{age - }{0 - 4.8 (years old)}
#'   \item{country - }{Bangladesh, Gambia, Kenya, or Mali}
#' }
#' @source \doi{10.1186/gb-2014-15-6-r76} and \doi{10.1093/nar/gkx1027}
#' @family Built-In Datasets
#' @usage gems
#' @name gems
#' @examples
#' gems
#' table(gems$metadata$country)
#' @export
delayedAssign(
  x     = 'gems', 
  value = local({
    f <- system.file(package = 'rbiom', 'extdata', 'gems.bz2')
    if (nzchar(f)) rbiom::as_rbiom(f) else NULL
  })
)


#' Longitudinal Stool Samples from Infants (n = 2,684)
#'
#' @format An rbiom object with 2,684 samples. 
#'         Includes metadata and taxonomy.
#' \describe{
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
#' @source \url{https://www.nature.com/articles/s41467-018-04641-7} and \doi{10.1038/s41467-017-01973-8}
#' @family Built-In Datasets
#' @usage babies
#' @name babies
#' @examples
#' babies
#' head(babies$metadata$Age)
#' @export
delayedAssign(
  x     = 'babies', 
  value = local({
    f <- system.file(package = 'rbiom', 'extdata', 'babies.bz2')
    if (nzchar(f)) rbiom::as_rbiom(f) else NULL
  })
)
