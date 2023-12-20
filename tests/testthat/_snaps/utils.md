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
      

