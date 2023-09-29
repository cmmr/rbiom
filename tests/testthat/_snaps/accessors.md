# biom_info

    Code
      biom_info(hmp50)
    Output
      $id
      [1] "Human Microbiome Project - 50 Sample Demo"
      
      $type
      [1] "OTU table"
      
      $format
      [1] "1.0.0"
      
      $format_url
      [1] "http://biom-format.org"
      
      $generated_by
      [1] "rbiom 1.0.2.9026"
      
      $date
      [1] "2021-07-01T18:48:45Z"
      
      $matrix_type
      [1] "sparse"
      
      $matrix_element_type
      [1] "int"
      
      $shape
      [1] 490  50
      
      $comment
      [1] ""
      

# sample_sums

    Code
      sample_sums(hmp5)
    Output
      HMP01 HMP02 HMP03 HMP04 HMP05 
       1660  1371  1353  1895  3939 

---

    Code
      sample_sums(hmp5, long = TRUE)
    Output
        Sample Reads
      1  HMP01  1660
      2  HMP02  1371
      3  HMP03  1353
      4  HMP04  1895
      5  HMP05  3939

---

    Code
      sample_sums(hmp5, long = TRUE, md = TRUE)
    Output
        Sample Reads Age BMI     Body Site    Sex
      1  HMP01  1660  22  20 Buccal mucosa Female
      2  HMP02  1371  24  23 Buccal mucosa   Male
      3  HMP03  1353  28  26        Saliva   Male
      4  HMP04  1895  25  23        Saliva   Male
      5  HMP05  3939  27  24 Buccal mucosa Female

---

    Code
      sample_sums(min5, long = TRUE, md = TRUE)
    Output
        Sample Reads
      1  HMP01  1660
      2  HMP02  1371
      3  HMP03  1353
      4  HMP04  1895
      5  HMP05  3939

# rarefaction_level

    Code
      rarefaction_level(hmp5)
    Output
      NULL

---

    Code
      rarefaction_level(rare5)
    Output
      [1] 1353

---

    Code
      rarefaction_level(as_percent(rare5))
    Output
      [1] 1

# has_tree

    Code
      has_tree(hmp5)
    Output
      [1] TRUE

---

    Code
      has_tree(min5)
    Output
      [1] FALSE

# has_sequences

    Code
      has_sequences(hmp5)
    Output
      [1] TRUE

---

    Code
      has_sequences(min5)
    Output
      [1] FALSE

# is_rarefied

    Code
      is_rarefied(hmp5)
    Output
      [1] FALSE

---

    Code
      is_rarefied(rare5)
    Output
      [1] TRUE

# n_samples

    Code
      n_samples(hmp5)
    Output
      [1] 5

# n_otus

    Code
      n_otus(hmp5)
    Output
      [1] 132

# n_ranks

    Code
      n_ranks(hmp5)
    Output
      [1] 6

---

    Code
      n_ranks(min5)
    Output
      NULL

