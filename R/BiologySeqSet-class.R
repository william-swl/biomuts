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

# names
setMethod("names", "BiologySeqSet", function(x) {
  return(names(x@DNA))
})

setMethod("names<-", "BiologySeqSet", function(x, value) {
  names(x@DNA) <- value
  names(x@AA) <- value
  return(x)
})
