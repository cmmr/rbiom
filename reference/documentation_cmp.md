# documentation_cmp

documentation_cmp

## Metadata Comparisons

Prefix metadata fields with `==` or `!=` to limit comparisons to within
or between groups, respectively. For example, `stat.by = '==Sex'` will
run calculations only for intra-group comparisons, returning "Male" and
"Female", but NOT "Female vs Male". Similarly, setting
`stat.by = '!=Body Site'` will only show the inter-group comparisons,
such as "Saliva vs Stool", "Anterior nares vs Buccal mucosa", and so on.

The same effect can be achieved by using the `within` and `between`
parameters. `stat.by = '==Sex'` is equivalent to
`stat.by = 'Sex', within = 'Sex'`.
