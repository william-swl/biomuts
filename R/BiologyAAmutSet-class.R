#' a S4 class to operate AA mutation set
#'
#' @slot muts list of AA mutations
#'
#' @include BiologyAAmut-class.R
#'
#' @export
#'
setClass("BiologyAAmutSet",
  slots = c(
    muts = "list"
  ),
  prototype = list(
    muts = list()
  )
)



#' create BiologyAAmutSet object
#'
#' @param mut_list list of AA mutations
#'
#' @return `BiologyAAmutSet` object
#' @export
#'
#' @examples
#' mut_list <- list(
#'   mut1 = c("D123G", "D6F", "C878C", "C878T"),
#'   mut2 = c("D123G", "D6F"),
#'   mut3 = c("C878C", "D123G", "C878T")
#' )
#'
#' muts <- BiologyAAmutSet(mut_list)
#'
BiologyAAmutSet <- function(mut_list) {
  trans_class <- function(x) {
    if (is(x, "character")) {
      res <- BiologyAAmut(x)
    } else if (is(x, "BiologyAAmut")) {
      res <- x
    }
    return(res)
  }

  res <- purrr::map(mut_list, trans_class)

  if (is.null(names(mut_list))) {
    names(res) <- seq_along(res)
  }

  new("BiologyAAmutSet",
    muts = res
  )
}


# show
#' @export
setMethod("show", "BiologyAAmutSet", function(object) {
  cat(is(object)[[1]], "\n")
  for (n in names(object@muts)) {
    n_text <- stringr::str_trunc(n, 9, side = "right", ellipsis = "...")
    n_text <- stringr::str_pad(
      n_text, max(nchar(names(object@muts))),
      side = "right", pad = " "
    )
    mut_vec <- object@muts[[n]]@mut
    mut_text <- paste(mut_vec, collapse = ",")
    cat(n_text, "|(", length(mut_vec), ") ", mut_text, "\n", sep = "")
  }
})

# names
setMethod("names", "BiologyAAmutSet", function(x) {
  return(names(x@muts))
})

setMethod("names<-", "BiologyAAmutSet", function(x, value) {
  names(x@muts) <- value
  return(x)
})




# subsettable
#' @export
setMethod("[",
  signature = signature(
    x = "BiologyAAmutSet",
    i = "numeric", j = "missing"
  ),
  function(x, i, j = "missing") BiologyAAmutSet(x@muts[i])
)


#' @export
setMethod("[[", "BiologyAAmutSet", function(x) x@muts[[i]]
)


# sort
#' @export
setMethod("sort", "BiologyAAmutSet", function(x, decreasing = FALSE) {
  res <- lapply(x@muts, function(y) sort(y, decreasing = decreasing))
  res <- BiologyAAmutSet(res)
  return(res)
})


# unique
#' @export
setMethod("unique", "BiologyAAmutSet", function(x, method = "term") {
  if (method == "mut") {
    res <- lapply(x@muts, function(y) unique(y))
  } else if (method == "site") {
    res <- lapply(x@muts, function(y) unique(y, bysite = TRUE))
  } else if (method == "term") {
    res <- x@muts[!duplicated(x@muts)]
  }

  res <- BiologyAAmutSet(res)
  return(res)
})


# select mutations
#' @export
setMethod("select_mut", "BiologyAAmutSet", function(x, start, end) {
  res <- lapply(x@muts, function(y) select_mut(y, start, end))
  res <- BiologyAAmutSet(res)

  return(res)
})
