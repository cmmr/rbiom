# adiv_boxplot_1

    Code
      p[["stats"]]
    Output
      $pairwise
          n         Test         Group1        Group2   .p.val   .adj.p
      1  20 Mann-Whitney Anterior nares Buccal mucosa 9.71e-01 9.71e-01
      2  20 Mann-Whitney Anterior nares    Mid vagina 1.50e-03 1.88e-03
      3  20 Mann-Whitney Anterior nares        Saliva 1.08e-05 4.33e-05
      4  19 Mann-Whitney Anterior nares         Stool 2.17e-05 4.33e-05
      5  20 Mann-Whitney  Buccal mucosa    Mid vagina 1.50e-03 1.88e-03
      6  20 Mann-Whitney  Buccal mucosa        Saliva 2.17e-05 4.33e-05
      7  19 Mann-Whitney  Buccal mucosa         Stool 1.45e-03 1.88e-03
      8  20 Mann-Whitney     Mid vagina        Saliva 1.08e-05 4.33e-05
      9  19 Mann-Whitney     Mid vagina         Stool 2.17e-05 4.33e-05
      10 19 Mann-Whitney         Saliva         Stool 2.99e-03 3.32e-03
      

# adiv_boxplot_2

    Code
      p[["stats"]]
    Output
      $groupwise
           Sex .metric  n           Test .p.val .adj.p
      1 Female    <NA> 60 Kruskal-Wallis 0.0131 0.0263
      2   Male    <NA> 38 Kruskal-Wallis 0.0740 0.0740
      

# adiv_boxplot_3

    Code
      p[["stats"]]
    Output
      $groupwise
             Body Site  n         Test .p.val .adj.p
      1 Anterior nares 10 Mann-Whitney  1.000      1
      2  Buccal mucosa 10 Mann-Whitney  0.310      1
      3     Mid vagina 10         <NA>  1.000      1
      4         Saliva 10 Mann-Whitney  0.548      1
      5          Stool  9 Mann-Whitney  0.905      1
      

