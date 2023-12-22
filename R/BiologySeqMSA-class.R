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
    consSeq = "BiologySeq",
    consDNA = "DNAString",
    consDNAfreq = "numeric",
    consAA = "AAString",
    consAAfreq = "numeric"
  ),
  prototype = list(
    params = list(),
    consSeq = BiologySeq(""),
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
#' @param method msa method, `ClustalOmega|Muscle|ClustalW`.
#' @param order sequences order, `input` to keep input order, `aligned` to use
#' aligned order
#' @param skip_align skip align or not
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
BiologySeqMSA <- function(x, corr_gaps = TRUE, skip_align = FALSE,
                          order = "input", method = "Muscle") {
  DNA <- DNA(x)

  if (skip_align == FALSE) {
    aln <- msa::msa(DNA, type = "dna", order = order, method = method)
    params <- S4Vectors::params(aln)
    params <- c(list("method" = aln@version), params)

    alnDNA <- Biostrings::DNAStringSet(aln)
  } else {
    alnDNA <- DNA
    if (any(nchar(alnDNA) != max(nchar(alnDNA)))) {
      stop("the nchar of sequences should be identical if align skipped")
    }
    params <- list("method" = "skip align")
  }

  if (corr_gaps == TRUE) {
    alnDNA <- purrr::map_chr(alnDNA, DNA_gaps_corr)
  }

  # consensus
  alnbs <- BiologySeqSet(alnDNA)

  # DNA consensus
  codes <- c("-", "A", "G", "C", "T")
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

  cons_seq <- BiologySeq(cons_dna)

  if (cons_seq@AA != cons_aa) {
    warning("AA(consSeq) is not same as consAA")
  }


  new("BiologySeqMSA",
    alnbs,
    params = params,
    consSeq = cons_seq,
    consDNA = cons_dna,
    consDNAfreq = cons_dna_freq,
    consAA = cons_aa,
    consAAfreq = cons_aa_freq
  )
}

# show
setMethod("show", "BiologySeqMSA", function(object) {
  cat(is(object)[[1]], "\n")
  cat(" @consSeq: ", is(object@consSeq)[[1]], "\n")
  cat("  @DNA: ", toString(object@consSeq@DNA), "\n")
  cat("  @AA: ", toString(object@consSeq@AA), "\n")
  cat(" @DNA: ")
  show(object@DNA)
  cat(" @AA:  ")
  show(object@AA)
})

# subsettable
#' @export
setMethod("[[", "BiologySeqMSA", function(x, i) BiologySeq(x@DNA[[i]]))

setGeneric("consSeq", function(x) standardGeneric("consSeq"))
#' @export
setMethod("consSeq", "BiologySeqMSA", function(x) x@consSeq)


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

setGeneric("aln_params", function(x) standardGeneric("aln_params"))
#' @export
setMethod("aln_params", "BiologySeqMSA", function(x) x@params)



setGeneric(
  "call_AAmutSet",
  function(querys, ...) standardGeneric("call_AAmutSet")
)

#' call AA mutations from `BiologySeqMSA` object
#'
#' @param querys BiologySeqMSA object
#' @param ref `consensus`, or the name of sequence
#'
#' @return `BiologyAAmutSet` object
#' @export
setMethod(
  "call_AAmutSet", "BiologySeqMSA",
  function(querys, ref = "consensus", rm_ref = TRUE) {
    if (ref == "consensus") {
      ref_ob <- querys@consSeq # nolint
    } else {
      ref_ob <- querys[[ref]]
      if (rm_ref == TRUE) {
        querys <- querys[-which(names(querys) == ref)]
      }
    }

    res <- purrr::map(
      querys,
      ~ call_mut(query = .x, ref = ref_ob, ignore_silence = TRUE)$AA
    )
    res <- BiologyAAmutSet(res)

    return(res)
  }
)
