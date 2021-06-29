#!/usr/bin/env Rscript

################################################
# Load packages
################################################

packages <- c(
  'topGO',
  'org.Hs.eg.db',
  'org.Mm.eg.db',
  'dplyr',
  'optparse',
  'tidyverse',
  'ggplot2',
  'Rgraphviz'
)
loaded <- lapply(X = packages, function(x) suppressMessages(suppressWarnings(require(x, character.only = TRUE))))

################################################
# Set up options and check input
################################################

script_name <- 'go_enrich'
start_time <- Sys.time()
message(paste0("Starting ", script_name, " at ", start_time))

# Parse parameters
option_list = list(
  make_option(c("--gene_list"), type="character", 
              help="Path to file containing input genes. Each gene has to be separated by a newline.", metavar="character"),
  make_option(c("--gene_universe"), type="character", 
              help="Path to file containing genes to test enrichment against. Each gene has to be separated by a newline.", metavar="character"),
  make_option(c("--outdir"), type="character",
              help="Directory to write output to.", metavar="character"),
  make_option(c("--species"), type="character", default="human", 
              help="Species from which the data is derived. Currently, only human or mouse are supported. [default = %default]", metavar="character"),
  make_option(c("--gene_ids"), type="character", default="SYMBOL", 
              help="Gene ids of genes in gene_list, e.g. ALIAS, SYMBOL, ENSEMBL, ENTREZID. For a complete list, see the AnnotationDbi vignette [default = %default]", metavar="character"),
  make_option(c("--node_size"), type="numeric", default=5, 
              help="An integer larger or equal to 1, used to prune the GO hierarchy from the terms which have less than node_size annotated genes (after the true path rule is applied). Values between 5 and 10 are recommended. [default = %default]", metavar="character"),
  make_option(c("--pvalue_cutoff"), type="numeric", default=0.05, 
              help="A numeric value between 0 and 1, used to remove non-signifcant GO terms from the result table. [default = %default]", metavar="character"),
  make_option(c("--dag"), type="logical", default=TRUE, 
              help="If TRUE, a directed acyclic graph of the top 5 enriched GO terms is created. [default = %default]", metavar="character"),
  make_option(c("--dag_format"), type="character", default="svg", 
              help="Output format for DAG: svg or pdf. [default = %default]", metavar="character")
); 
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Used for testing
dotest <- FALSE
if(dotest){
  setwd("/home/rstudio/02_code/go_enrichment/go_enrich")
  opt$gene_list <- 'example_input/gene_list.txt'
  opt$gene_universe <- 'example_input/gene_universe.txt'
  opt$outdir <- 'example_output'
  opt$gene_ids <- "SYMBOL"
  opt$species <- "human"
  opt$pvalue_cutoff <- 0.05
  opt$node_size <- 10
  opt$dag <- TRUE
  opt$dag_format <- "svg"
}

# Check gene list
if (is.null(opt$gene_list)) { # Check if file was provided
  stop(paste0(script_name, '::No input gene list given.'))
} else if(!file.exists(opt$gene_list)){ # Check if file exists
  stop(paste0(script_name, '::Input file does not exist. Check input: ', opt$gene_list))
}

# Check gene universe
if (is.null(opt$gene_universe)) { # Check if file was provided
  stop(paste0(script_name, '::No input gene universe given.'))
} else if(!file.exists(opt$gene_universe)){ # Check if file exists
  stop(paste0(script_name, '::Input file does not exist. Check input: ', opt$gene_universe))
}

# Check output directory
if (!file.exists(opt$outdir)){ # Check if directory exists
  stop(paste0(script_name, '::Output directory does not exist. Check input: ', opt$outdir))
}

# Select gene database to use, and check species
if(tolower(opt$species) == 'human') {
  annotationdbi_obj <- org.Hs.eg.db
} else if(tolower(opt$species) == 'mouse') {
  annotationdbi_obj <- org.Mm.eg.db
} else{
  stop(paste0(script_name, '::Invalid species: ', opt$species, '. Current version supports human or mouse.'))
}

# Check p-value cutoff
if (opt$pvalue_cutoff < 0 | opt$pvalue_cutoff > 1) {
  stop(paste0(script_name, '::P-value cutoff must be >= 0 or <= 1. Check input: ', opt$pvalue_cutoff))
}

# Check node_size
if (opt$node_size <= 0) {
  stop(paste0(script_name, '::Node size must be > 0. Check input: ', opt$node_size))
}

# Check gene ids
if(!opt$gene_ids %in% keytypes(annotationdbi_obj)) {
  stop(paste0(script_name, '::Invalid Gene id. Check input: ', opt$gene_ids, '. Must be one of: ', paste(keytypes(annotationdbi_obj), collapse=", "), '.' ))
}

# Check DAG boolean
if (!opt$dag %in% c(TRUE, FALSE)) { # Check if dag format is correct
  stop(paste0(script_name, '::DAG boolean must be "TRUE" or "FALSE". Check input: ', opt$dag))
}

# Check DAG format
if (!opt$dag_format %in% c("svg", "pdf")) { # Check if dag format is correct
  stop(paste0(script_name, '::DAG format must be "svg" or "pdf". Check input: ', opt$dag_format))
}

################################################
# GO enrichment
################################################

# Read gene universe
gene_universe <- read.table(opt$gene_universe, stringsAsFactors=FALSE)[,1]

# Create table containing GO IDs per genes
gene_universe_annotated <- AnnotationDbi::select(
  x = annotationdbi_obj,
  keys = gene_universe,
  columns = c(opt$gene_ids, 'GO'),
  keytype = opt$gene_ids
)

# Filter out duplicates
gene_universe_annotated <- unique(gene_universe_annotated[,c(1, 2, 4)])

# Create gene to go mapping
gene_2_go <- unstack(unique(gene_universe_annotated[,c(2, 1)])) # unique because different evidence codes create duplicate entries

# Read gene list  
gene_list_raw <- read.table(opt$gene_list, stringsAsFactors=FALSE)[,1]

# Remove genes that don't have GO annotation
gene_list <- gene_list_raw[gene_list_raw %in% gene_universe_annotated[,1]]
message(paste0(script_name, "::Genes removed without GO annotation: ", sum(!gene_list %in% gene_universe_annotated[,1]), " out of ", length(gene_list_raw)))

# Define all genes vector. All genes start with a 0 value
allGenes <- rep(0, length(gene_universe))
names(allGenes) <- gene_universe

# Genes that are in the gene list get a 1 value (allGenes has to be a numeric vector)
for (g in gene_list) {
  allGenes[names(allGenes) == g] <- 1
}

# Loop over each ontology
for (ont in c("BP", "CC", "MF")) {
  
  # Create topGOdata object
  GOdata <- new("topGOdata",
                ontology=ont,
                allGenes=allGenes,
                annot=annFUN.gene2GO,
                gene2GO=gene_2_go,
                geneSel=function(x) return(x == 1), # Select genes in gene_list
                nodeSize=opt$node_size) # pruning GO tree

  # Run enrichment test
  result <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")

  no_of_go_terms_scored <- length(result@score)
  
  # Print some statistics
  #message(paste0(script_name, "::Number of genes from the gene universe that can be used in the analysis: ", sum(GOdata@feasible)), " out of ", length(allGenes))
  #message(paste0(script_name, "::Number of genes from the gene list that can be used in the analysis: ", numSigGenes(GOdata)), " out of ", length(gene_list))
  #message(paste0(script_name, "::Number of GO terms scored: ", no_of_go_terms_scored))
  
  # Create results table of all enriched GO terms
  result_table_full <- GenTable(GOdata, classicFisher=result, topNodes=no_of_go_terms_scored, numChar=1000)
  result_table_full$classicFisher[result_table_full$classicFisher %in% c("<1e-30", "< 1e-30")] <- "1e-30"
  result_table_full$classicFisher <- as.numeric(result_table_full$classicFisher)
  
  # Remove non-significant GO terms
  result_table_full <- dplyr::filter(result_table_full, classicFisher < opt$pvalue_cutoff)
  
  # Retrieve all genes per GO term from GOdata object
  allGO <- genesInTerm(GOdata)
  
  # Rename columns
  result_table_full <- result_table_full %>% dplyr::rename("nAnnotated" = "Annotated", 
                                      "nSignificant" = "Significant")
  
  # Add genes for nAnnotated
  result_table_full$Annotated <- unlist(lapply(result_table_full$GO.ID, function(x) paste0(allGO[x][[1]], collapse=", ")))
  
  # Add genes for nSignficant
  result_table_full$Significant <- unlist(lapply(result_table_full$Annotated, function(x) paste0(strsplit(x, ", ")[[1]][strsplit(x, ", ")[[1]] %in% gene_list], collapse=", ")))

  # Rename column
  result_table_full <- result_table_full %>% dplyr::rename("Pvalue" = "classicFisher")
  
  # Add Rank
  result_table_full$Rank <- 1:nrow(result_table_full)

  # Rearrange and remove "Expected" column (no documentation on what exactly this is)
  result_table_full <- result_table_full %>% dplyr::select("Rank", "GO.ID", "Term", "Pvalue", "nAnnotated", "nSignificant", "Annotated", "Significant")

  # Save results table
  write.table(result_table_full, sprintf("%s/go_enrich_%s.tsv", opt$outdir, GOdata@ontology), quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
  
  # Select top 10 GO terms
  result_table_top10 <- result_table_full[1:10,]

  # Format new column as combination of GO ID and Term
  result_table_top10$Term <- paste0(result_table_top10$Term, " (", result_table_top10$GO.ID, ")")
  result_table_top10$Term <- factor(result_table_top10$Term, levels=rev(result_table_top10$Term))
  
  # Set xlab for plotting
  if (GOdata@ontology == "BP") {
    xtitle <- "Biological Process"
  } else if (GOdata@ontology == "MF") {
    xtitle <- "Molecular Function"
  } else if (GOdata@ontology == "CC") {
    xtitle <- "Molecular Function"
  }
  
  plt <- ggplot(result_table_top10, aes(x=Term, y=-log10(Pvalue))) +
    stat_summary(geom = "bar", fun = mean, position = "dodge") +
    scale_y_continuous(breaks = round(seq(0, -log10(1e-30), by = 2), 1)) +
    xlab(xtitle) +
    ylab("Enrichment (-log10 Pval)") +
    ggtitle(sprintf("Go Enrichment %s", GOdata@ontology)) +
    theme_bw() +
    guides(colour=guide_legend(override.aes=list(size=2.5))) +
    coord_flip() 
  ggsave(sprintf("%s/go_enrich_%s_top10_barplot.png",  opt$outdir, GOdata@ontology), height=4, width=10)
  
  # Select top 30 GO terms
  # If there are less than 30 significant GO terms, dont select 30 (would introduce NAs...)
  n_top_gos <- min(30, nrow(result_table_full))
  result_table_top30 <- result_table_full[1:n_top_gos,]
  
  # Format new column as combination of GO ID and Term
  result_table_top30$Term <- paste0(result_table_top30$Term, " (", result_table_top30$GO.ID, ")")
  result_table_top30$Term <- factor(result_table_top30$Term, levels=rev(result_table_top30$Term))
  
  # Calculate GeneRatio
  result_table_top30$GeneRatio <- round(result_table_top30$nSignificant / result_table_top30$nAnnotated, 2)
  
  # Calculate GeneCount
  result_table_top30$GeneCount <- result_table_top30$nSignificant
  result_table_top30$Term <- factor(result_table_top30$Term, levels = rev(result_table_top30$Term)) # fixes order
  
  # Create plot
  plt <- ggplot(result_table_top30, aes(x=Term, 
                                        y=GeneRatio, 
                                        size=GeneCount, 
                                        fill=-log10(Pvalue))) +
    geom_point(shape = 21) +
    scale_fill_continuous(low = 'blue', high = 'red') +
    ylim(0, 1) + # consistently 0 to 1 Gene Ratio
    labs(
      x = "",
      y = "GeneRatio",
      fill = "Enrichment (-log10 p-value)",
      title = sprintf("GO enrichment %s", GOdata@ontology), 
      subtitle = sprintf("Top %s terms ordered by p-value", n_top_gos)
    ) +
    theme_bw() +
    coord_flip() + 
    theme(legend.text = element_text(size = 10)) # Legend text size
  ggsave(sprintf("%s/go_enrich_%s_top30_dotplot.png", opt$outdir, GOdata@ontology), height=8, width=10)
    
  if (opt$dag) {
    # Create and save DAG
    if (opt$dag_format == "svg") {
      svg(sprintf("%s/go_enrich_%s_top5_DAG.svg", opt$outdir, GOdata@ontology), height=20, width=20)
    } else {
      pdf(sprintf("%s/go_enrich_%s_top5_DAG.pdf", opt$outdir, GOdata@ontology), height=20, width=20)
    }
    showSigOfNodes(GOdata, score(result), firstSigNodes = 5, useInfo = "all")
    dev.off()
  }
}

end_time <- Sys.time()
message(paste0("Finished ", script_name, " at ", end_time))
