# What is go_enrich?

This script is used to calculate enriched gene ontology (GO) terms from a set of input genes and a gene universe using the `topGO` and `AnnotationDbi` packages.

# Usage

The following parameters can be passed to the script. Default values are in [square brackets].

 * `--gene_list`: Path to file containing input genes. Each gene has to be separated by a newline.
 * `--gene_universe`: Path to file containing genes to test enrichment against. Each gene has to be separated by a newline.
 * `--outdir`: Directory to write output to.
 * `--species`: Species from which the data is derived. Currently, only human or mouse are supported. [human]
 * `--gene_ids`: Gene ids of genes in gene_list, e.g. ALIAS, SYMBOL, ENSEMBL, ENTREZID. For a complete list, see the [AnnotationDbi vignette](https://bioconductor.org/packages/release/bioc/vignettes/AnnotationDbi/inst/doc/IntroToAnnotationPackages.pdf). [SYMBOL]
 * `--node_size`: An integer larger or equal to 1, used to prune the GO hierarchy from the terms which have less than `node_size` annotated genes (after the [true path rule](https://homes.di.unimi.it/~valentini/papers/vale.TPR.hier.revised.pdf) is applied). Values between 5 and 10 are recommended. [5]
 * `--pvalue_cutoff`: A numeric value between 0 and 1, used to remove non-signifcant GO terms from the result table. [0.05]
 * `--dag`: If TRUE, a directed acyclic graph of the top 5 enriched GO terms is created. [TRUE]
 * `--dag_format`: output format for DAG: svg or pdf. [svg]

# Example

```
Rscript go_enrich.R  \
  --gene_list example_input/gene_list.txt \
  --gene_universe example_input/gene_universe.txt \  
  --outdir example_output \ 
  --species human \
  --gene_ids SYMBOL \
  --node_size 10 \
  --pvalue_cutoff 0.05 \
  --dag TRUE \
  --dag_format svg
```

# Output

See `example_output`