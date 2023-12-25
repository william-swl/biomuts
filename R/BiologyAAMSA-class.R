#' a S4 class to operate aligned AA sequence set
#'
#' @slot AA AAStringSet.
#' @slot params MSA params
#' @slot consAA consensus AAString
#' @slot consAAfreq AA consensus frequence
#' @slot AAnumbering numbering
#'
#' @include BiologySeqMSA-class.R
#'
#' @export
#'
setClass("BiologyAAMSA",
  slots = c(
    AA = "AAStringSet",
    params = "list",
    consAA = "AAString",
    consAAfreq = "numeric",
    AAnumbering = "character"
  ),
  prototype = list(
    AA = Biostrings::AAStringSet(),
    params = list(),
    consAA = Biostrings::AAString(),
    consAAfreq = c(),
    AAnumbering = c()
  )
)

setValidity("BiologyAAMSA", function(object) {
  consAA_length <- Biostrings::nchar(object@consAA)

  if (any(duplicated(object@AAnumbering))) {
    "the numbering label must be unique!"
  } else if (length(object@AAnumbering) != consAA_length) {
    "the numbering label length must be identical with sequence AA!"
  } else {
    TRUE
  }
})



#' create BiologyAAMSA object
#'
#' @param x AAStringSet
#' @param skip_align skip align or not
#' @param order sequences order, `input` to keep input order, `aligned` to use
#' aligned order
#' @param method msa method, `ClustalOmega|Muscle|ClustalW`.
#'
#' @return BiologyAAMSA
#' @export
#'
#' @examples
#' AA <- c("MQVNPTE", "MQVTE", "MQVTV")
#'
#' aln <- BiologyAAMSA(AA)
#'
#' AA(aln)
#'
#' consAA(aln)
#'
#' aln[[2]]
#'
BiologyAAMSA <- function(x, skip_align = FALSE,
                         order = "input", method = "Muscle") {
  if (is.null(names(x))) {
    names(x) <- seq_along(x)
  }

  if (skip_align == FALSE) {
    aln <- msa::msa(x, type = "protein", order = order, method = method)
    params <- S4Vectors::params(aln)
    params <- c(list("method" = aln@version), params)

    alnAA <- Biostrings::AAStringSet(aln)
  } else {
    alnAA <- x
    if (any(Biostrings::nchar(alnAA) != max(Biostrings::nchar(alnAA)))) {
      stop("the nchar of sequences should be identical if align skipped")
    }
    params <- list("method" = "skip align")
  }

  # AA consensus
  consmtx <- Biostrings::consensusMatrix(alnAA)
  codes <- rownames(consmtx)
  consmtx <- round(consmtx / length(alnAA), 3)
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

  consAA_length <- Biostrings::nchar(cons_aa)

  consAA_numbering <- as.character(seq_len(consAA_length))
  names(consAA_numbering) <- as.character(seq_len(consAA_length))

  new("BiologyAAMSA",
    AA = alnAA,
    params = params,
    consAA = cons_aa,
    consAAfreq = cons_aa_freq,
    AAnumbering = consAA_numbering
  )
}


# show
setMethod("show", "BiologyAAMSA", function(object) {
  cat(is(object)[[1]], "\n")
  cat(" @consAA: ", toString(object@consAA), "\n")
  cat(" @AA: \n")
  show(object@AA)
  cat(" @AAnumbering:  ")
  show(object@AAnumbering)
})

# names
setMethod("names", "BiologyAAMSA", function(x) {
  return(names(x@AA))
})

# subsettable
#' @export
setMethod("[", "BiologyAAMSA", function(x, i) x@AA[i])
#' @export
setMethod("[[", "BiologyAAMSA", function(x, i) Biostrings::AAString(x@AA[[i]]))

#' @export
setMethod("AA", "BiologyAAMSA", function(x) x@AA)

#' @export
setMethod("consAA", "BiologyAAMSA", function(x) x@consAA)

#' @export
setMethod("consAAfreq", "BiologyAAMSA", function(x) x@consAAfreq)

#' @export
setMethod("aln_params", "BiologyAAMSA", function(x) x@params)

#' @export
setMethod("length", "BiologyAAMSA", function(x) length(x@AA))

#' @export
setMethod("AAnumbering", "BiologyAAMSA", function(x) x@AAnumbering)

#' @export
setMethod(
  "AAnumbering<-", "BiologyAAMSA",
  function(x, value) {
    value <- as.character(value)
    consAA_length <- Biostrings::nchar(x@consAA)
    if (is.null(names(value))) {
      names(value) <- seq_len(consAA_length)
    }
    x@AAnumbering <- value
    validObject(x)

    return(x)
  }
)


#' call AA mutations from `BiologyAAMSA` object
#'
#' @param querys BiologyAAMSA object
#' @param ref `consensus`, or the name of sequence
#'
#' @return `BiologyAAmutSet` object
#' @export
setMethod(
  "call_AAmutSet", "BiologyAAMSA",
  function(querys, ref = "consensus", rm_ref = TRUE) {
    AAnumbering <- querys@AAnumbering

    if (ref == "consensus") {
      ref_ob <- querys@consAA # nolint
    } else {
      ref_ob <- querys[[ref]]
      if (rm_ref == TRUE) {
        querys <- querys[-which(names(querys) == ref)]
      }
    }

    res <- purrr::map(
      querys,
      ~ call_AAmut(query = .x, ref = ref_ob)
    )
    res <- BiologyAAmutSet(res)
    numbering(res) <- AAnumbering

    return(res)
  }
)
