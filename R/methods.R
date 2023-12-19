# get slot
setGeneric("DNA", function(x) standardGeneric("DNA"))
setMethod("DNA", "BiologySeq", function(x) x@DNA)
setMethod("DNA", "BiologySeqSet", function(x) x@DNA)


setGeneric("AA", function(x) standardGeneric("AA"))
setMethod("AA", "BiologySeq", function(x) x@AA)
setMethod("AA", "BiologySeqSet", function(x) x@AA)


# set slot
setGeneric("DNA<-", function(x, value) standardGeneric("DNA<-"))

setMethod("DNA<-", "BiologySeq", function(x, value) {
  x <- BiologySeq(value)
  return (x)
})

setMethod("DNA<-", "BiologySeqSet", function(x, value) {
  x <- BiologySeqSet(DNA=value)
  return (x)
})

# DNAsite_by_AA
setGeneric("DNAsite_by_AA",
           function(x, start, end) standardGeneric("DNAsite_by_AA"),
           signature='x')

setMethod("DNAsite_by_AA", "BiologySeq", function(x, start, end) {
  start_nt <- 3 * (start - 1) + 1
  end_nt <- 3 * (end - 1) + 3
  res <- XVector::subseq(x@DNA, start=start_nt, end=end_nt) %>%
    toString %>% as.data.frame %>%
    dplyr::rename(seq_nt=1) %>%
    dplyr::mutate(start=start_nt, end=end_nt, .before=1)

  return(res)
})

setMethod("DNAsite_by_AA", "BiologySeqSet", function(x, start, end) {
  start_nt <- 3 * (start - 1) + 1
  end_nt <- 3 * (end - 1) + 3
  res <- XVector::subseq(x@DNA, start=start_nt, end=end_nt) %>%
    as.data.frame %>%
    dplyr::rename(seq_nt=x) %>%
    dplyr::mutate(start=start_nt, end=end_nt, .before=1)

  return(res)
})



# AAsite_by_DNA
setGeneric("AAsite_by_DNA",
           function(x, start, end) standardGeneric("AAsite_by_DNA"),
           signature='x')

setMethod("AAsite_by_DNA", "BiologySeq", function(x, start, end) {
  start_aa <- floor((start - 1) / 3) + 1
  end_aa <- floor((end - 1) / 3) + 1
  res <- XVector::subseq(x@AA, start=start_aa, end=end_aa) %>%
    toString %>% as.data.frame %>%
    dplyr::rename(seq_aa=1) %>%
    dplyr::mutate(start=start_aa, end=end_aa, .before=1)

  return(res)
})

setMethod("AAsite_by_DNA", "BiologySeqSet", function(x, start, end) {
  start_aa <- floor((start - 1) / 3) + 1
  end_aa <- floor((end - 1) / 3) + 1
  res <- XVector::subseq(x@AA, start=start_aa, end=end_aa) %>%
    as.data.frame %>%
    dplyr::rename(seq_aa=1) %>%
    dplyr::mutate(start=start_aa, end=end_aa, .before=1)

  return(res)
})



