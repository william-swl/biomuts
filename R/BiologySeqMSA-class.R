#' a S4 class to operate aligned DNA together with AA sequence set
#'
#' @slot DNA DNAStringSet.
#' @slot AA AAStringSet.
#' @slot params MSA params
#' @slot consDNA consensus DNAString
#' @slot consDNAfreq DNA consensus frequency
#' @slot consAA consensus AAString
#' @slot consAAfreq AA consensus frequence
#'
#' @include BiologySeqSet-class.R
#' @export
#'
setClass("BiologySeqMSA",
  contains = "BiologySeqSet",
  slots = c(
    params = "list",
    consDNA = "DNAString",
    consDNAfreq = "numeric",
    consAA = "AAString",
    consAAfreq = "numeric"
  ),
  prototype = list(
    params = list(),
    consDNA = Biostrings::DNAString(),
    consDNAfreq = c(),
    consAA = Biostrings::AAString(),
    consAAfreq = c()
  )
)

#' create BiologySeqMSA object
#'
#' @param x BiologySeqSet
#' @param corr_gaps correct gaps codon frame after MSA, `TRUE` as default.
#' see `?DNA_gaps_corr`
#'
#' @return BiologySeqMSA object
#' @export
#'
#' @examples
#'
#' bss <- BiologySeqSet(c("ATGAAA---", "ATG---AAA"))
#'
#' BiologySeqMSA(bss)
#'
BiologySeqMSA <- function(x, corr_gaps = TRUE) {
  DNA <- DNA(x)
  aln <- msa::msa(DNA, type = "dna", order = "input")
  params <- S4Vectors::params(aln)
  alnDNA <- Biostrings::DNAStringSet(aln)

  if (corr_gaps == TRUE) {
    alnDNA <- purrr::map_chr(alnDNA, DNA_gaps_corr)
  }

  # consensus
  alnbs <- BiologySeqSet(alnDNA)

  # DNA consensus
  codes <- c("A", "C", "T", "G", "-")
  consmtx <- Biostrings::consensusMatrix(DNA(alnbs))[codes, ]
  consmtx <- round(consmtx / length(DNA(alnbs)), 3)
  cons_dna_freq <- apply(consmtx, 2, max)
  cons_bool <- consmtx == matrix(rep(cons_dna_freq, length(codes)),
    nrow = length(codes), byrow = TRUE
  )
  cons_dna <- Biostrings::DNAString(
    paste(rep("-", dim(consmtx)[2]), collapse = "")
  )

  for (i in seq_along(codes)) {
    cons_dna[cons_bool[i, ], ] <- codes[i]
  }

  # AA consensus
  consmtx <- Biostrings::consensusMatrix(AA(alnbs))
  codes <- rownames(consmtx)
  consmtx <- round(consmtx / length(AA(alnbs)), 3)
  cons_aa_freq <- apply(consmtx, 2, max)
  cons_bool <- consmtx == matrix(rep(cons_aa_freq, length(codes)),
    nrow = length(codes), byrow = TRUE
  )
  cons_aa <- Biostrings::AAString(
    paste(rep("-", dim(consmtx)[2]), collapse = "")
  )

  for (i in seq_along(codes)) {
    cons_aa[cons_bool[i, ], ] <- codes[i]
  }


  new("BiologySeqMSA",
    alnbs,
    params = params,
    consDNA = cons_dna,
    consDNAfreq = cons_dna_freq,
    consAA = cons_aa,
    consAAfreq = cons_aa_freq
  )
}

# show
setMethod("show", "BiologySeqMSA", function(object) {
  cat(is(object)[[1]], "\n")
  cat(" @consDNA: ")
  show(object@consDNA)
  cat(" @consAA: ")
  show(object@consAA)
  cat(" @DNA: ")
  show(object@DNA)
  cat(" @AA:  ")
  show(object@AA)
})



setGeneric("consDNA", function(x) standardGeneric("consDNA"))
#' @export
setMethod("consDNA", "BiologySeqMSA", function(x) x@consDNA)

setGeneric("consDNAfreq", function(x) standardGeneric("consDNAfreq"))
#' @export
setMethod("consDNAfreq", "BiologySeqMSA", function(x) x@consDNAfreq)


setGeneric("consAA", function(x) standardGeneric("consAA"))
#' @export
setMethod("consAA", "BiologySeqMSA", function(x) x@consAA)

setGeneric("consAAfreq", function(x) standardGeneric("consAAfreq"))
#' @export
setMethod("consAAfreq", "BiologySeqMSA", function(x) x@consAAfreq)
