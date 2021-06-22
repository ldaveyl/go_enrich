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
  make_option(c("--outdir"), type="character", default=".",
              help="Directory to write output to [default = go_enich]", metavar="character"),
  make_option(c("--species"), type="character", default='human', 
              help="Species from which the data is derived. Currently, only human or mouse are supported [default = %default]", metavar="character"),
  make_option(c("--gene_ids"), type="character", default='SYMBOL', 
              help="Gene ids, e.g. ALIAS, SYMBOL, ENSEMBL, ENTREZID. For a complete list, see the AnnotationDBI vignette [default = %default]", metavar="character"),
  make_option(c("--node_size"), type="numeric", default='5', 
              help="An integer larger or equal to 1, used to prune the GO hierarchy from the terms which have less than nodeSize annotated genes. Values between 5 and 10 are recommended. [default = %default]", metavar="character"),
  make_option(c("--pvalue_cutoff"), type="numeric", default='0.05', 
              help="A numeric value between 0 and 1, used to remove non-signifcant GO terms from the result table. [default = %default]", metavar="character")
); 
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# TESTPARAM
dotest <- FALSE
if(dotest){
  opt$gene_list <- '~/Documents/contribution-day/20210621_contribution_day/go_enrich/gene_list.txt'
  opt$gene_universe <- '~/Documents/contribution-day/20210621_contribution_day/go_enrich/gene_universe.txt'
  opt$outdir <- '~/Documents/contribution-day/20210621_contribution_day/go_enrich/test_output'
  opt$node_size <- 5
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
  stop(paste0(script_name_name, '::Input file does not exist. Check input: ', opt$gene_universe))
}

# Check output directory
if (!file.exists(opt$outdir)){ # Check if file exists
  stop(paste0(script_name_name, '::Output directory does not exist. Check input: ', opt$outdir))
}

# Select gene database to use, and check species
if(tolower(opt$species) == 'human') {
  annotationdbi_obj <- org.Hs.eg.db
} else if(tolower(opt$species) == 'mouse') {
  annotationdbi_obj <- org.Mm.eg.db
} else{
  stop(paste0(script_name, '::Invalid species: ', opt$species, '. Current version supports human or mouse.'))
}

################################################
# GO enrichment
################################################

# Read gene universe
gene_universe <- read.table(opt$gene_universe, stringsAsFactors=FALSE)[,1]

# Create table containing GO IDs per genes
gene_universe_annotated <- AnnotationDbi::select(
  x = org.Hs.eg.db,
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

# Define all genes vector. All genes start with a FALSE value
allGenes <- rep(FALSE, length(gene_universe))
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
  result_table_full <- dplyr::filter(result_table_full, classicFisher < 0.05)
  
  # Save results table
  write.table(result_table_full, sprintf("%s/go_enrich_%s.tsv", opt$outdir, GOdata@ontology), quote=FALSE, row.names=TRUE, col.names=TRUE, sep="\t")
  
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
  
  plt <- ggplot(result_table_top10, aes(x=Term, y=-log10(classicFisher))) +
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
  result_table_top30$GeneRatio <- round(result_table_top30$Significant / result_table_top30$Annotated, 2)
  
  # Calculate GeneCount
  result_table_top30$GeneCount <- result_table_top30$Significant
  result_table_top30$Term <- factor(result_table_top30$Term, levels = rev(result_table_top30$Term)) # fixes order
  
  # Create plot
  plt <- ggplot(result_table_top30, aes(x=Term, 
                                        y=GeneRatio, 
                                        size=GeneCount, 
                                        fill=-log10(classicFisher))) +
    geom_point(shape = 21) +
    scale_fill_continuous(low = 'blue', high = 'red') +
    ylim(0, 1) + # consistently 0 to 1 Gene Ratio
    labs(
      x = "",
      y = "Gene Ratio",
      fill = "Enrichment (-log10 p-value)",
      title = sprintf("GO enrichment %s", GOdata@ontology), 
      subtitle = sprintf("Top %s terms ordered by p-value", n_top_gos)
    ) +
    theme_bw() +
    coord_flip() + 
    theme(legend.text = element_text(size = 10)) # Legend text size
  ggsave(sprintf("%s/go_enrich_%s_top30_dotplot.png", opt$outdir, GOdata@ontology), height=8, width=10)
  
  # Create and save DAG
  svg(sprintf("%s/go_enrich_%s_top5_DAG.svg", opt$outdir, GOdata@ontology), height=20, width=20)
  showSigOfNodes(GOdata, score(result), firstSigNodes = 5, useInfo = "all")
  dev.off()
}

end_time <- Sys.time()
message(paste0("Finished ", script_name, " at ", end_time))













  
  
