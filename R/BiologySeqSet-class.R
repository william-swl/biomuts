#' a S4 class to operate DNA together with AA sequence set
#'
#' @slot DNA DNAStringSet.
#' @slot AA AAStringSet.
#'
#' @include BiologySeq-class.R
#'
#' @export
#'
setClass("BiologySeqSet",
  slots = c(
    DNA = "DNAStringSet",
    AA = "AAStringSet"
  ),
  prototype = list(
    DNA = Biostrings::DNAStringSet(),
    AA = Biostrings::AAStringSet()
  )
)

#' create BiologySeqSet object
#'
#' @param DNA character vector of DNA letters
#' @param corr_gaps see `?DNA_gaps_corr`
#'
#' @return BiologySeqSet object
#'

#' @export
#'
#' @examples
#'
#' bss <- BiologySeqSet(c("ATGAAA---", "ATG---AAA"))
#'
#' DNA(bss)
#'
#' AA(bss)
#'
BiologySeqSet <- function(DNA, corr_gaps = FALSE) {
  no3 <- which(Biostrings::nchar(DNA) %% 3 != 0)
  if (length(no3) > 0) {
    no3text <- stringr::str_c(no3, collapse = ",") # nolint
    stop(stringr::str_glue(
      "check sequence {no3text}: nchar(@DNA) is not multiples of 3!"
    ))
  }

  if (corr_gaps == TRUE) {
    DNA <- purrr::map_chr(DNA, DNA_gaps_corr)
  }
  AA <- purrr::map_chr(DNA, dna2aa)
  if (is.null(names(DNA))) {
    names(DNA) <- seq_along(DNA)
    names(AA) <- seq_along(AA)
  }

  new("BiologySeqSet",
    DNA = Biostrings::DNAStringSet(DNA),
    AA = Biostrings::AAStringSet(AA)
  )
}

# valid
setValidity("BiologySeqSet", function(object) {
  no3 <- which(Biostrings::nchar(object@DNA) %% 3 != 0)

  if (length(no3) > 0) {
    no3text <- stringr::str_c(no3, collapse = ",")
    stringr::str_glue(
      "check sequence {no3text}: nchar(@DNA) is not multiples of 3!"
    )
  } else {
    TRUE
  }
})

# show
setMethod("show", "BiologySeqSet", function(object) {
  cat(is(object)[[1]], "\n @DNA: ")
  show(object@DNA)
  cat(" @AA:  ")
  show(object@AA)
})

# subsettable
#' @export
setMethod("[", "BiologySeqSet", function(x, i) BiologySeqSet(x@DNA[i]))
#' @export
setMethod("[[", "BiologySeqSet", function(x, i) BiologySeq(x@DNA[[i]]))


# names
setMethod("names", "BiologySeqSet", function(x) {
  return(names(x@DNA))
})

setMethod("names<-", "BiologySeqSet", function(x, value) {
  names(x@DNA) <- value
  names(x@AA) <- value
  return(x)
})

#' @export
setMethod("length", "BiologySeqSet", function(x) length(x@DNA))

#' @export
setMethod("DNA", "BiologySeqSet", function(x) x@DNA)

#' @export
setMethod("AA", "BiologySeqSet", function(x) x@AA)

#' @export
setMethod("DNA<-", "BiologySeqSet", function(x, value) {
  x <- BiologySeqSet(DNA = value)
  return(x)
})


# DNAsite_by_AA
#' @export
setMethod("DNAsite_by_AA", "BiologySeqSet", function(x, start, end) {
  start_nt <- 3 * (start - 1) + 1
  end_nt <- 3 * (end - 1) + 3
  res <- XVector::subseq(x@DNA, start = start_nt, end = end_nt) %>%
    as.data.frame() %>%
    dplyr::rename(seq_nt = x) %>%
    dplyr::mutate(start = start_nt, end = end_nt, .before = 1)

  return(res)
})

# AAsite_by_DNA
#' @export
setMethod("AAsite_by_DNA", "BiologySeqSet", function(x, start, end) {
  start_aa <- floor((start - 1) / 3) + 1
  end_aa <- floor((end - 1) / 3) + 1
  res <- XVector::subseq(x@AA, start = start_aa, end = end_aa) %>%
    as.data.frame() %>%
    dplyr::rename(seq_aa = 1) %>%
    dplyr::mutate(start = start_aa, end = end_aa, .before = 1)

  return(res)
})
