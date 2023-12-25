
<!-- README.md is generated from README.Rmd. Please edit that file -->

# biomuts

<!-- badges: start -->
<!-- badges: end -->

## installation

You can install the development version of `biomuts` like so:

``` r
devtools::install_github("william-swl/biomuts")
```

And use `biomuts` like so:

    library(biomuts)
    library(msa)

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

alnbs
#> BiologySeqMSA 
#>  @consSeq:  BiologySeq 
#>   @DNA:  ATGCAGGTA------ACTGAG 
#>   @AA:  MQV--TE 
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
#>  @AAnumbering:    1   2   3   4   5   6   7 
#> "1" "2" "3" "4" "5" "6" "7"

alnbs[[1]]
#> BiologySeq
#>  @DNA: ATGCAGGTAAACCCTACTGAG
#>  @AA:  MQVNPTE

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

aln_params(alnbs)[1:5]
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
```

- if the frequency of letters are the same, show with the priority of
  `TCGA-`

``` r
bss <- BiologySeqSet(c("ACC", "ATC", "GTA", "GCA"))

alnbs <- BiologySeqMSA(bss, method = "ClustalOmega")
#> using Gonnet
```

- call AA mutations from `BiologySeqMSA` object

``` r
bss <- BiologySeqSet(
  c("ATGCAGGTAAACCCTACTGAG", "ATGCAGGTTACTGAG", "ATGCAGGTAACTGTG")
)
alnbs <- BiologySeqMSA(bss)

call_AAmutSet(alnbs)
#> BiologyAAmutSet 
#> @muts
#> 1|(2) -4N,-5P
#> 2|(0) 
#> 3|(1) E7V
#> @numbering
#>   1   2   3   4   5   6   7 
#> "1" "2" "3" "4" "5" "6" "7"

call_AAmutSet(alnbs, ref = 3)
#> BiologyAAmutSet 
#> @muts
#> 1|(3) -4N,-5P,V7E
#> 2|(1) V7E
#> @numbering
#>   1   2   3   4   5   6   7 
#> "1" "2" "3" "4" "5" "6" "7"
```

### BiologyAAMSA

- a S4 class to operate aligned AA sequence set
- please set `method="ClustalOmega"` if [Muscle is
  crashed](https://stackoverflow.com/questions/76663781/using-msa-package-in-r-and-it-is-crashing)

``` r
AA <- c("MQVNPTE", "MQVTE", "MQVTV")

aln <- BiologyAAMSA(AA)

AA(aln)
#> AAStringSet object of length 3:
#>     width seq                                               names               
#> [1]     7 MQVNPTE                                           1
#> [2]     7 MQV--TE                                           2
#> [3]     7 MQV--TV                                           3

consAA(aln)
#> 7-letter AAString object
#> seq: MQV--TE

consAAfreq(aln)
#> [1] 1.000 1.000 1.000 0.667 0.667 1.000 0.667

aln[[2]]
#> 7-letter AAString object
#> seq: MQV--TE

call_AAmutSet(aln)
#> BiologyAAmutSet 
#> @muts
#> 1|(2) -4N,-5P
#> 2|(0) 
#> 3|(1) E7V
#> @numbering
#>   1   2   3   4   5   6   7 
#> "1" "2" "3" "4" "5" "6" "7"
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
#> @muts
#> mut1|(4) D123G,D6F,C878C,C878T
#> mut2|(2) D123G,D6F
#> @numbering
#> NULL

muts[[3]]
#> BiologyAAmut 
#>  @mut(3): C878C D123G C878T

select_mut(muts, 1, 130)
#> BiologyAAmutSet 
#> @muts
#> mut1|(2) D123G,D6F
#> mut2|(2) D123G,D6F
#> mut3|(1) D123G
#> @numbering
#> NULL

numbering(muts) <- c("123" = "site123", "6" = "site6", "87" = "site87")

number_muts(muts)[[1]]
#> [1] "D[site123]G" "D[site6]F"   "C878C"       "C878T"

# sort(muts)

# unique(muts, method = "mut")

# unique(muts, method = "site")

muts2 <- BiologyAAmutSet(list("D12C", "D12C"))

muts2
#> BiologyAAmutSet 
#> @muts
#> 1|(1) D12C
#> 2|(1) D12C
#> @numbering
#> NULL

# unique(muts2, method = "term")
```

## utils

- correct the gaps location to fit codon frames

``` r
DNA_gaps_corr("AT---G")
#> [1] "ATG---"
```

- translate DNA to AA (gap allowed)

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

- call AA mutations from two `Biostrings::AAString` objects

``` r
call_AAmut(Biostrings::AAString("MQVNPTE"), Biostrings::AAString("MQVCTTE"))
#> [1] "C4N" "T5P"
```

- count AA mutations from `BiologyAAmutSet` objects

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

- amino acid features

``` r
aa_info("C")
#>   hydrophilicity polarity charge volume
#> C          -0.79    -0.85      0  -0.39
```

- compare feature index between AAs

``` r
compare_aa(c("A", "T"), "C")
#>     hydrophilicity polarity charge volume
#> A>C          -0.05    -0.64      0   0.23
#> T>C          -0.10    -0.76      0  -0.10
```

- number the mutations

``` r
number_mut(c("A12T", "C18P"), c("12" = "12", "18" = "18"))
#> [1] "A[12]T" "C[18]P"
```
