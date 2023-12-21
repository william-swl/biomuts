
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

bss[[1]]
#> BiologySeq
#>  @DNA: ATGAAACCC
#>  @AA:  MKP

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
alnbs <- BiologySeqMSA(bss)
#> using Gonnet

alnbs
#> BiologySeqMSA 
#>  @consSeq:  BiologySeq 
#>   @DNA:  ATGCAGGTAACT------GAG 
#>   @AA:  MQVT--E 
#>  @DNA: DNAStringSet object of length 3:
#>     width seq                                               names               
#> [1]    21 ATGCAGGTAAACCCTACTGAG                             1
#> [2]    21 ATGCAGGTTACT------GAG                             2
#> [3]    21 ATGCAGGTAACT------GTG                             3
#>  @AA:  AAStringSet object of length 3:
#>     width seq                                               names               
#> [1]     7 MQVNPTE                                           1
#> [2]     7 MQVT--E                                           2
#> [3]     7 MQVT--V                                           3

alnbs[[1]]
#> BiologySeq
#>  @DNA: ATGCAGGTAAACCCTACTGAG
#>  @AA:  MQVNPTE

consDNA(alnbs)
#> 21-letter DNAString object
#> seq: ATGCAGGTAACT------GAG

consDNAfreq(alnbs)
#>  [1] 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 0.667 1.000 0.667 0.667
#> [13] 0.667 0.667 0.667 0.667 0.667 0.667 1.000 0.667 1.000

consAA(alnbs)
#> 7-letter AAString object
#> seq: MQVT--E

consAAfreq(alnbs)
#> [1] 1.000 1.000 1.000 0.667 0.667 0.667 0.667

aln_params(alnbs)[1:5]
#> $method
#> [1] "ClustalOmega 1.2.0"
#> 
#> $gapOpening
#> [1] "default"
#> 
#> $gapExtension
#> [1] "default"
#> 
#> $maxiters
#> [1] 0
#> 
#> $verbose
#> [1] FALSE
```

- if the frequency of letters are the same, show with the priority of
  `TCGA-`

``` r
bss <- BiologySeqSet(c("ACC", "ATC", "GTA", "GCA"))

alnbs <- BiologySeqMSA(bss)
#> using Gonnet
```

### BiologyAAmut

- a S4 class to operate AA mutation

``` r
mut <- BiologyAAmut(c("D123G", "D6F", "C878C", "D123G", "C878T"))

mut
#> BiologyAAmut 
#>  @mut(5): D123G D6F C878C D123G C878T

# sort(mut)

# unique(sort(mut))

# unique(sort(mut), bysite = TRUE)

select_mut(mut, 1, 150)
#> BiologyAAmut 
#>  @mut(3): D123G D6F D123G
```

### BiologyAAmutSet

``` r
mut_list <- list(
  mut1 = c("D123G", "D6F", "C878C", "C878T"),
  mut2 = c("D123G", "D6F"),
  mut3 = c("C878C", "D123G", "C878T")
)

muts <- BiologyAAmutSet(mut_list)

muts[1:2]
#> BiologyAAmutSet 
#> mut1|(4) D123G,D6F,C878C,C878T
#> mut2|(2) D123G,D6F

muts[[3]]
#> BiologyAAmut 
#>  @mut(3): C878C D123G C878T

select_mut(muts, 1, 130)
#> BiologyAAmutSet 
#> mut1|(2) D123G,D6F
#> mut2|(2) D123G,D6F
#> mut3|(1) D123G

# sort(muts)

# unique(muts, method = "mut")

# unique(muts, method = "site")

muts2 <- BiologyAAmutSet(list("D12C", "D12C"))

muts2
#> BiologyAAmutSet 
#> 1|(1) D12C
#> 2|(1) D12C

# unique(muts2, method = "term")
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

- call DNA and AA mutations from two `BiologySeq` objects

``` r
call_mut(BiologySeq("ATGCCCTTT"), BiologySeq("CTTCCCTTC"))
#> $DNA
#>     1     3     9 
#> "C1A" "T3G" "C9T" 
#> 
#> $AA
#>     1     3 
#> "L1M" "F3F" 
#> 
#> $table
#> # A tibble: 2 × 9
#>   codon_start codon_changes ref_nt query_nt site_aa ref_aa query_aa valid
#>   <chr>       <chr>         <chr>  <chr>    <chr>   <chr>  <chr>    <chr>
#> 1 1           2             CTT    ATG      1       L      M        TRUE 
#> 2 7           1             TTC    TTT      3       F      F        TRUE 
#> # ℹ 1 more variable: silence <chr>
```

- call DNA and AA mutations from two `BiologySeq` objects

``` r
mut_list <- list(
  mut1 = c("D123G", "D6F", "C878C", "C878T"),
  mut2 = c("D123G", "D6F"),
  mut3 = c("C878C", "D123G", "C878T")
)

muts <- BiologyAAmutSet(mut_list)

count_muts(muts)
#> # A tibble: 4 × 4
#>   mut_aa mut1  mut2  mut3 
#>   <chr>  <lgl> <lgl> <lgl>
#> 1 D6F    TRUE  TRUE  FALSE
#> 2 D123G  TRUE  TRUE  TRUE 
#> 3 C878C  TRUE  FALSE TRUE 
#> 4 C878T  TRUE  FALSE TRUE
```
