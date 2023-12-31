#' a S4 class to operate AA mutation
#'
#' @slot mut AA mutations
#' @slot site the sites of mutations
#'
#' @export
#'
setClass("BiologyAAmut",
  slots = c(
    mut = "character",
    site = "numeric"
  ),
  prototype = list(
    mut = c(),
    site = c()
  )
)

#' create BiologyAAmut object
#'
#' @param mut vector of mutations
#'
#' @return BiologyAAmut object
#' @export
#'
#' @examples
#' mut <- BiologyAAmut(c("D123G", "D6F", "C878C", "D123G", "C878T"))
#'
#' mut
#'
#' sort(mut)
#'
#' unique(sort(mut))
#'
#' unique(sort(mut), bysite = TRUE)
#'
BiologyAAmut <- function(mut) {
  aa_alphabeta <- paste(Biostrings::AA_ALPHABET[1:28], collapse = "") # nolint
  match_mtx <- stringr::str_match(
    mut,
    stringr::str_glue("^[{aa_alphabeta}](\\d+)[{aa_alphabeta}]$")
  )
  nomatch <- mut[is.na(match_mtx[, 1])]
  nomatch_text <- stringr::str_c(nomatch, collapse = ", ") # nolint

  if (length(nomatch) > 0) {
    stop(stringr::str_glue(
      "please input mutations as character vector like D123G: {nomatch_text}"
    ))
  }

  site <- as.numeric(match_mtx[, 2])

  new("BiologyAAmut",
    mut = mut,
    site = site
  )
}


# show
#' @export
setMethod("show", "BiologyAAmut", function(object) {
  cat(is(object)[[1]], "\n")
  cat(" @mut", "(", length(object@mut), "): ", sep = "")
  cat(object@mut)
})


# sort
#' @export
setMethod("sort", "BiologyAAmut", function(x, decreasing = FALSE) {
  ord <- order(x@site, decreasing = decreasing)
  res <- BiologyAAmut(x@mut[ord])

  return(res)
})

# unique
#' @export
setMethod("unique", "BiologyAAmut", function(x, bysite = FALSE) {
  if (bysite == TRUE) {
    mask <- !duplicated(x@site)
    res <- BiologyAAmut(x@mut[mask])
  } else {
    res <- BiologyAAmut(unique(x@mut))
  }

  return(res)
})


# select mutations
setGeneric("select_mut", function(x, start, end) standardGeneric("select_mut"))
#' @export
setMethod("select_mut", "BiologyAAmut", function(x, start, end) {
  mask <- x@site %in% start:end
  res <- BiologyAAmut(x@mut[mask])

  return(res)
})
