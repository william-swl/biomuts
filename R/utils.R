#' correct the gaps location to fit codon frames
#'
#' @param x character of DNA letters
#' @param dismiss_no3gaps ignore gap widths that are not multiples of 3
#'
#' @return character of DNA letters
#' @export
#'
#' @examples DNA_gaps_corr("AT---G")
#'
DNA_gaps_corr <- function(x, dismiss_no3gaps = FALSE) {
  if (!is(x, "DNAString")) {
    x <- Biostrings::DNAString(x)
  }

  gaps <- IRanges::reduce(Biostrings::matchPattern("-", x))

  if (any(IRanges::width(gaps) %% 3 != 0) && dismiss_no3gaps == FALSE) {
    stop("Some of the gaps are not multiples of 3!")
  }

  for (i in seq_along(gaps)) {
    if (IRanges::width(gaps) %% 3 == 0) {
      start <- IRanges::start(gaps)[i]
      end <- IRanges::end(gaps)[i]
      if (start %% 3 == 2) {
        temp_letter <- x[start - 1]
        x[start - 1] <- "-"
        x[end] <- temp_letter
      }
      if (start %% 3 == 0) {
        temp_letter <- x[end + 1]
        x[end + 1] <- "-"
        x[start] <- temp_letter
      }
    }
  }

  return(toString(x))
}


#' translate DNA to AA (gaps allowed)
#'
#' @param x character of DNA letters, or Biostrings::DNAString object
#'
#' @return character of AA letters
#' @export
#'
#' @examples dna2aa("ATG---AAA")
dna2aa <- function(x) {
  if (!is(x, "DNAString")) {
    x <- Biostrings::DNAString(x)
  }
  gaps <- IRanges::reduce(Biostrings::matchPattern("-", x))

  if (any(IRanges::width(gaps) %% 3 != 0)) {
    stop("Some of the gaps are not multiples of 3!")
  }
  if (any(IRanges::start(gaps) %% 3 != 1)) {
    stop("Some of the gaps are out of frame!")
  }

  start_aa <- floor((IRanges::start(gaps) - 1) / 3) + 1
  end_aa <- floor((IRanges::end(gaps) - 1) / 3) + 1

  x <- Biostrings::DNAString(gsub("---", "TAA", x))

  AA <- Biostrings::translate(x, no.init.codon = TRUE)

  for (i in seq_along(start_aa)) {
    start <- start_aa[i]
    end <- end_aa[i]
    AA[start:end] <- "-"
  }
  return(toString(AA))
}
