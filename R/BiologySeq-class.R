#' a S4 class to operate DNA together with AA sequence
#'
#' @slot DNA DNAString.
#' @slot AA AAString.
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
