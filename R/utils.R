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


#' call DNA and AA mutations from two BiologySeq objects
#'
#' @param query BiologySeq
#' @param ref BiologySeq
#' @param ignore_silence ignore AA silence mutations or not
#'
#' @return list
#' @export
#'
#' @examples
#'
#' ref <- BiologySeq("CTTCCCTTC")
#' query <- BiologySeq("ATGCCCTTT")
#' call_mut(query, ref)
#'
call_mut <- function(query, ref, ignore_silence = FALSE) {
  if (!is(query, "BiologySeq")) {
    stop("query must be BiologySeq object!")
  }
  if (!is(ref, "BiologySeq")) {
    stop("ref must be BiologySeq object!")
  }

  queryDNA <- DNA(query)
  refDNA <- DNA(ref)
  queryAA <- AA(query)
  refAA <- AA(ref)

  mut_idx <- which(as.matrix(queryDNA) != as.matrix(refDNA))

  site_aa <- floor((mut_idx - 1) / 3) + 1

  site_counts <- table(site_aa)

  site_aa <- unique(site_aa)

  # mut str
  ref_vec <- as.character(as.matrix(refDNA[mut_idx]))
  query_vec <- as.character(as.matrix(queryDNA[mut_idx]))
  mut_dna <- stringr::str_c(ref_vec, mut_idx, query_vec, sep = "")
  names(mut_dna) <- mut_idx



  ref_vec <- as.character(as.matrix(refAA[site_aa]))
  query_vec <- as.character(as.matrix(queryAA[site_aa]))
  if (ignore_silence == TRUE) {
    mask <- query_vec != ref_vec
    ref_vec <- ref_vec[mask]
    query_vec <- query_vec[mask]
    site_aa <- site_aa[mask]
  }
  mut_aa <- stringr::str_c(ref_vec, site_aa, query_vec, sep = "")
  names(mut_aa) <- site_aa

  # mut table
  site2mutTable <- function(x) {
    site_nt <- 3 * (x - 1) + 1

    res <- c(
      "codon_start" = site_nt,
      "codon_changes" = unname(site_counts[toString(x)]),
      "ref_nt" = toString(refDNA[(site_nt):(site_nt + 2)]),
      "query_nt" = toString(queryDNA[site_nt:(site_nt + 2)]),
      "site_aa" = x,
      "ref_aa" = toString(refAA[x]),
      "query_aa" = toString(queryAA[x])
    )

    res["valid"] <- res["query_aa"] == dna2aa(res["query_nt"])
    res["silence"] <- res["ref_aa"] == res["query_aa"]

    return(res)
  }

  mut_tb <- purrr::map_dfr(site_aa, site2mutTable)

  res <- list(
    DNA = mut_dna,
    AA = mut_aa,
    table = mut_tb
  )

  return(res)
}


#' call AA mutations from two AAString objects
#'
#' @param query Biostrings::AAString
#' @param ref Biostrings::AAString
#'
#' @return vector of AA mutations
#' @export
#'
#' @examples
#'
#' call_AAmut(
#'   Biostrings::AAString("MQVNPTE"),
#'   Biostrings::AAString("MQVCTTE")
#' )
call_AAmut <- function(query, ref) {
  if (!is(query, "AAString")) {
    stop("query must be AAString object!")
  }
  if (!is(ref, "AAString")) {
    stop("ref must be AAString object!")
  }

  site_aa <- which(as.matrix(query) != as.matrix(ref))

  # mut str
  ref_vec <- as.character(as.matrix(ref[site_aa]))
  query_vec <- as.character(as.matrix(query[site_aa]))
  mut_aa <- stringr::str_c(ref_vec, site_aa, query_vec, sep = "")

  return(mut_aa)
}


#' count AA mutations from `BiologyAAmutSet` object
#'
#' @param muts `BiologyAAmutSet` object
#'
#' @return tibble
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
#' count_muts(muts)
#'
count_muts <- function(muts) {
  if (!is(muts, "BiologyAAmutSet")) {
    stop("muts must be BiologyAAmut object!")
  }
  mut_list <- purrr::map(muts@muts, function(x) x@mut)
  mut_vec <- sort(BiologyAAmut(unique(unlist(mut_list))))@mut
  count_tb <- mut_list %>% purrr::map_dfr(~ mut_vec %in% .x)
  res <- count_tb %>% dplyr::mutate(mut_aa = mut_vec, .before = 1)

  return(res)
}


#' amino acid features
#'
#' @param x aa
#' @param raw return feature raw index instead of scaled index
#'
#' @return data.frame
#' @export
#'
#' @examples aa_info("C")
aa_info <- function(x, raw = FALSE) {
  if (raw == TRUE) {
    res <- aa_feature_raw[x, ]
  } else {
    res <- aa_feature[x, ]
  }

  return(res)
}


#' compare feature index between AAs
#'
#' @param query query aa
#' @param ref ref aa
#'
#' @return data.frame
#' @export
#'
#' @examples compare_aa(c("A", "T"), "C")
compare_aa <- function(query, ref, raw = FALSE) {
  if (length(ref) == 1) {
    ref <- rep(ref, length(query))
  } else if (length(ref) != length(query)) {
    stop("length(ref) must be 1 or length(query)")
  }

  if (raw == TRUE) {
    res <- aa_feature_raw[ref, ] - aa_feature_raw[query, ]
  } else {
    res <- aa_feature[ref, ] - aa_feature[query, ]
  }

  rownames(res) <- stringr::str_glue("{query}>{ref}")
  return(res)
}
