
<!-- README.md is generated from README.Rmd. Please edit that file -->

# biomuts

<!-- badges: start -->
<!-- badges: end -->

## installation

You can install the development version of `biomuts` like so:

``` r
devtools::install_github("william-swl/biomuts")
```

## S4 classes

### BiologySeq

- a S4 class to operate DNA together with AA sequence

``` r
bs <- BiologySeq("ATGAAA---CCCTTTGGG---")

DNA(bs)
#> 21-letter DNAString object
#> seq: ATGAAA---CCCTTTGGG---

AA(bs)
#> 7-letter AAString object
#> seq: MK-PFG-

DNA(bs) <- "ATGAAACCCTTTGGG"

bs
#> BiologySeq
#>  @DNA: ATGAAACCCTTTGGG
#>  @AA:  MKPFG

DNA(bs)[1:6] <- "A"

bs
#> BiologySeq
#>  @DNA: AAAAAACCCTTTGGG
#>  @AA:  KKPFG

# get DNA site by AA site
DNAsite_by_AA(bs, start = 2, end = 3)
#>   start end seq_nt
#> 1     4   9 AAACCC

# get AA site by DNA site
AAsite_by_DNA(bs, start = 1, end = 9)
#>   start end seq_aa
#> 1     1   3    KKP
```

### BiologySeqSet

- a S4 class to operate DNA together with AA sequence set

``` r
bss <- BiologySeqSet(c("ATGAAA---", "ATG---AAA"))

DNA(bss)
#> DNAStringSet object of length 2:
#>     width seq                                               names               
#> [1]     9 ATGAAA---                                         1
#> [2]     9 ATG---AAA                                         2

AA(bss)
#> AAStringSet object of length 2:
#>     width seq                                               names               
#> [1]     3 MK-                                               1
#> [2]     3 M-K                                               2

DNA(bss) <- c("ATGAAACCC", "ATGCCCAAA")

bss
#> BiologySeqSet 
#>  @DNA: DNAStringSet object of length 2:
#>     width seq                                               names               
#> [1]     9 ATGAAACCC                                         1
#> [2]     9 ATGCCCAAA                                         2
#>  @AA:  AAStringSet object of length 2:
#>     width seq                                               names               
#> [1]     3 MKP                                               1
#> [2]     3 MPK                                               2

# get DNA site by AA site
DNAsite_by_AA(bss, start = 2, end = 3)
#>   start end seq_nt
#> 1     4   9 AAACCC
#> 2     4   9 CCCAAA

# get AA site by DNA site
AAsite_by_DNA(bss, start = 3, end = 9)
#>   start end seq_aa
#> 1     1   3    MKP
#> 2     1   3    MPK
```

### BiologySeqMSA

- a S4 class to operate aligned DNA together with AA sequence set
- please set `method="ClustalOmega"` if [Muscle is
  crashed](https://stackoverflow.com/questions/76663781/using-msa-package-in-r-and-it-is-crashing)

``` r
bss <- BiologySeqSet(
  c("ATGCAGGTAAACCCTACTGAG", "ATGCAGGTTACTGAG", "ATGCAGGTAACTGTG")
)
alnbs <- BiologySeqMSA(bss, method = "Muscle")

alnbs
#> BiologySeqMSA 
#>  @consDNA: 21-letter DNAString object
#> seq: ATGCAGGTA------ACTGAG
#>  @consAA: 7-letter AAString object
#> seq: MQV--TE
#>  @DNA: DNAStringSet object of length 3:
#>     width seq                                               names               
#> [1]    21 ATGCAGGTAAACCCTACTGAG                             1
#> [2]    21 ATGCAGGTT------ACTGAG                             2
#> [3]    21 ATGCAGGTA------ACTGTG                             3
#>  @AA:  AAStringSet object of length 3:
#>     width seq                                               names               
#> [1]     7 MQVNPTE                                           1
#> [2]     7 MQV--TE                                           2
#> [3]     7 MQV--TV                                           3

consDNA(alnbs)
#> 21-letter DNAString object
#> seq: ATGCAGGTA------ACTGAG

consDNAfreq(alnbs)
#>  [1] 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 0.667 0.667 0.667 0.667
#> [13] 0.667 0.667 0.667 1.000 1.000 1.000 1.000 0.667 1.000

consAA(alnbs)
#> 7-letter AAString object
#> seq: MQV--TE

consAAfreq(alnbs)
#> [1] 1.000 1.000 1.000 0.667 0.667 1.000 0.667

aln_params(alnbs)
#> $method
#> [1] "MUSCLE 3.8.31  "
#> 
#> $gapOpening
#> [1] 400
#> 
#> $gapExtension
#> [1] 0
#> 
#> $maxiters
#> [1] 16
#> 
#> $verbose
#> [1] FALSE
#> 
#> $inputSeqIsFileFlag
#> [1] FALSE
#> 
#> $le
#> [1] FALSE
#> 
#> $sp
#> [1] FALSE
#> 
#> $sv
#> [1] FALSE
#> 
#> $spn
#> [1] TRUE
#> 
#> $brenner
#> [1] FALSE
#> 
#> $diags
#> [1] FALSE
#> 
#> $diags1
#> [1] FALSE
#> 
#> $diags2
#> [1] FALSE
#> 
#> $dimer
#> [1] FALSE
#> 
#> $noanchors
#> [1] FALSE
#> 
#> $nocore
#> [1] FALSE
#> 
#> $profile
#> [1] FALSE
#> 
#> $refine
#> [1] FALSE
#> 
#> $refinew
#> [1] FALSE
#> 
#> $spscore
#> [1] FALSE
#> 
#> $version
#> [1] FALSE
```

- if the frequency of letters are the same, show with the priority of
  `TCGA-`

``` r
bss <- BiologySeqSet(c("ACC", "ATC", "GTA", "GCA"))

alnbs <- BiologySeqMSA(bss)
#> using Gonnet
```

## utils

- correct the gaps location to fit codon frames

``` r
DNA_gaps_corr("AT---G")
#> [1] "ATG---"
```

- translate DNA to AA (gaps allowed)

``` r
dna2aa("ATG---AAA")
#> [1] "M-K"
```
