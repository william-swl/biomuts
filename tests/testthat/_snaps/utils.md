# dna2aa

    Code
      call_mut(BiologySeq("ATGCCCTTT"), BiologySeq("CTTCCCTTC"))
    Output
      $DNA
          1     3     9 
      "C1A" "T3G" "C9T" 
      
      $AA
          1     3 
      "L1M" "F3F" 
      
      $table
      # A tibble: 2 x 9
        codon_start codon_changes ref_nt query_nt site_aa ref_aa query_aa valid
        <chr>       <chr>         <chr>  <chr>    <chr>   <chr>  <chr>    <chr>
      1 1           2             CTT    ATG      1       L      M        TRUE 
      2 7           1             TTC    TTT      3       F      F        TRUE 
      # i 1 more variable: silence <chr>
      

# call_mut

    Code
      call_mut(query, ref)
    Output
      $DNA
          1     3     9 
      "C1A" "T3G" "C9T" 
      
      $AA
          1     3 
      "L1M" "F3F" 
      
      $table
      # A tibble: 2 x 9
        codon_start codon_changes ref_nt query_nt site_aa ref_aa query_aa valid
        <chr>       <chr>         <chr>  <chr>    <chr>   <chr>  <chr>    <chr>
      1 1           2             CTT    ATG      1       L      M        TRUE 
      2 7           1             TTC    TTT      3       F      F        TRUE 
      # i 1 more variable: silence <chr>
      

# call_AAmutSet

    Code
      call_AAmutSet(alnbs)
    Output
      BiologyAAmutSet 
      1|(2) -4N,-5P
      2|(0) 
      3|(1) E7V

---

    Code
      call_AAmutSet(alnbs, ref = "3")
    Output
      BiologyAAmutSet 
      1|(3) -4N,-5P,V7E
      2|(1) V7E

# count_muts

    Code
      count_muts(muts)
    Output
      # A tibble: 4 x 4
        mut_aa mut1  mut2  mut3 
        <chr>  <lgl> <lgl> <lgl>
      1 D6F    TRUE  TRUE  FALSE
      2 D123G  TRUE  TRUE  TRUE 
      3 C878C  TRUE  FALSE TRUE 
      4 C878T  TRUE  FALSE TRUE 

