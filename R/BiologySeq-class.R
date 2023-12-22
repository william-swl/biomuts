#' a S4 class to operate DNA together with AA sequence
#'
#' @slot DNA DNAString.
#' @slot AA AAString.
#'
#' @include utils.R
#'
#' @export
#'
setClass("BiologySeq",
  slots = c(
    DNA = "DNAString",
    AA = "AAString"
  ),
  prototype = list(
    DNA = Biostrings::DNAString(),
    AA = Biostrings::AAString()
  )
)

#' create BiologySeq object
#'
#' @param DNA character of DNA letters
#' @param corr_gaps see `?DNA_gaps_corr`
#'
#' @return BiologySeq object
#' @export
#'
#' @examples
#'
#' bs <- BiologySeq("ATGAAA---CCCTTTGGG---")
#'
#' DNA(bs)
#'
#' AA(bs)
#'
BiologySeq <- function(DNA, corr_gaps = FALSE) {
  if (Biostrings::nchar(DNA) %% 3 != 0) {
    stop("nchar(@DNA) is not multiples of 3!")
  }
  if (corr_gaps == TRUE) {
    DNA <- DNA_gaps_corr(DNA)
  }
  AA <- dna2aa(DNA)
  new("BiologySeq",
    DNA = Biostrings::DNAString(DNA),
    AA = Biostrings::AAString(AA)
  )
}

# valid
setValidity("BiologySeq", function(object) {
  if (Biostrings::nchar(object@DNA) %% 3 != 0) {
    "nchar(@DNA) is not multiples of 3!"
  } else {
    TRUE
  }
})

# show
setMethod("show", "BiologySeq", function(object) {
  cat(is(object)[[1]], "\n",
    " @DNA: ", toString(object@DNA), "\n",
    " @AA:  ", toString(object@AA), "\n",
    sep = ""
  )
})


# get slot
setGeneric("DNA", function(x) standardGeneric("DNA"))
#' @export
setMethod("DNA", "BiologySeq", function(x) x@DNA)

setGeneric("AA", function(x) standardGeneric("AA"))
#' @export
setMethod("AA", "BiologySeq", function(x) x@AA)


# set slot
setGeneric("DNA<-", function(x, value) standardGeneric("DNA<-"))
#' @export
setMethod("DNA<-", "BiologySeq", function(x, value) {
  x <- BiologySeq(value)
  return(x)
})


# DNAsite_by_AA
setGeneric("DNAsite_by_AA",
  function(x, start, end) standardGeneric("DNAsite_by_AA"),
  signature = "x"
)
#' @export
setMethod("DNAsite_by_AA", "BiologySeq", function(x, start, end) {
  start_nt <- 3 * (start - 1) + 1
  end_nt <- 3 * (end - 1) + 3
  res <- XVector::subseq(x@DNA, start = start_nt, end = end_nt) %>%
    toString() %>%
    as.data.frame() %>%
    dplyr::rename(seq_nt = 1) %>%
    dplyr::mutate(start = start_nt, end = end_nt, .before = 1)

  return(res)
})


# AAsite_by_DNA
setGeneric("AAsite_by_DNA",
  function(x, start, end) standardGeneric("AAsite_by_DNA"),
  signature = "x"
)
#' @export
setMethod("AAsite_by_DNA", "BiologySeq", function(x, start, end) {
  start_aa <- floor((start - 1) / 3) + 1
  end_aa <- floor((end - 1) / 3) + 1
  res <- XVector::subseq(x@AA, start = start_aa, end = end_aa) %>%
    toString() %>%
    as.data.frame() %>%
    dplyr::rename(seq_aa = 1) %>%
    dplyr::mutate(start = start_aa, end = end_aa, .before = 1)

  return(res)
})
