#!/usr/bin/env Rscript
# setwd("/cashew/users/sur/exp/fraserv/2020/today")
library(tidyverse)
library(broom)
library(argparser)

process_arguments <- function(){
  p <- arg_parser(paste("Associate genes with features"))
  
  # Positional arguments
  p <- add_argument(p, "genes",
                    help = paste("Table of genes by genome. Must have species column"),
                    type = "character")
  p <- add_argument(p, "feats",
                    help = paste("Table of feats by genome. Must have Genome column"),
                    type = "character")
  p <- add_argument(p, "prefix",
                    help = paste("Prefix for output files"),
                    type = "character")
  
  # Optional arguments
  p <- add_argument(p, "--outdir",
                     help = paste("Output directory"),
                     type = "character",
                     default = "output/")
  p <- add_argument(p, "--n_pcs",
                    help = paste("Number of PCs to use to control for relatedness"),
                    default = 10,
                    type = "numeric")
  p <- add_argument(p, "--min_genomes_test",
                    help = paste("Minimum number of genomes in which a gene needs",
                                 "to be present for being tested"),
                    default = 100,
                    type = "numeric")
                     
  # Read arguments
  cat("Processing arguments...\n")
  args <- parse_args(p)
  
  # Process arguments
  
  return(args)
}

#' Title
#'
#' @param dat 
#' @param gene_names 
#' @param feat_names 
#' @param meta_names 
#' @param trans_fun 
#'
#' @return
#' @export
#'
#' @examples
phyloprof_lm <- function(dat, gene_names, feat_names, meta_names = NULL, trans_fun = FALSE){
  
  if(!is.null(trans_fun)){
    cat("Transforming...\n")
    trans_fun <- match.fun(trans_fun)
    dat <- trans_fun(dat)
  }
  
  cat("Fitting models...\n")
  expand_grid(gene_name = gene_names, feat_name = feat_names) %>%
    pmap_dfr(function(gene_name, feat_name, dat = dat, meta_names = NULL){
      # gene_name <- "K00001"
      # feat_name <- "Number.of.duplications.in.the.genome"
      
      d <- dat[,c(gene_name, feat_name, meta_names)]
      colnames(d) <- c("gene", "feat", meta_names)
      # f1 <- paste(feat_name, "~", gene_name)
      f1 <- formula(feat ~ .)
      tidy(lm(f1, data = d)) %>%
        mutate(q.value = p.adjust(p.value, 'fdr'),
               feat = feat_name,
               gene = gene_name) %>%
        filter(term == "gene") %>%
        select(-term) 
    }, dat = dat, meta_names = meta_names)
}

#' Title
#'
#' @param genes 
#' @param gene_summaries 
#' @param prop_genomes 
#' @param n_pcs 
#'
#' @return
#' @export
#'
#' @examples
get_gene_pcs <- function(genes, gene_summaries, prop_genomes = 0.5, n_pcs = 10){
  # n_pcs <- 10
  
  # Select genes for PCA
  pca_genes <- gene_summaries$gene[gene_summaries$n_genomes > prop_genomes * nrow(genes) ]
  pca_genes <- intersect(colnames(genes), pca_genes)
  
  # Make PCA
  genome_ids <- genes$Genome
  pca_mat <- genes %>%
    select(all_of(pca_genes)) %>%
    as.matrix()
  # pca_mat
  # genes_pca <- prcomp(pca_mat)
  genes_pca <- gmodels::fast.prcomp(pca_mat, center = TRUE, scale. = TRUE)
  # genes_pca <- gmodels::fast.prcomp(pca_mat)
  
  # Crate output
  genes_pca$x[, 1:n_pcs] %>%
    as_tibble() %>%
    bind_cols(Genome = genome_ids)
  # genes_pca
  # summary(genes_pca)
  # biplot(genes_pca)
  # heatmap(pca_mat)
  # genes_svd <- svd(pca_mat, nu = 10)
}

args <- process_arguments()
# args <- list(genes = "/cashew/shared_data/grins/2020-05-12.eggnot_annots/tabs/KEGG_KOs.txt.gz",
#              feats = "/cashew/users/nivina/2020-03-24.Streptomyces_GRINS_detection/GRINS_detected_in_genomes_and_BGCs.txt",
#              outdir = "output",
#              prefix = "KEGG_KOs",
#              n_pcs = 10,
#              min_genomes_test = 100)
# args <- list(genes = "/cashew/shared_data/grins/2020-05-12.eggnot_annots/tabs/BiGG_reactions.txt.gz",
#              feats = "/cashew/users/nivina/2020-03-24.Streptomyces_GRINS_detection/GRINS_detected_in_genomes_and_BGCs.txt",
#              outdir = "output",
#              prefix = "BiGG_reactions",
#              n_pcs = 10,
#              min_genomes_test = 100)


# Prepare output dir
if(!dir.exists(args$outdir)){
  dir.create(args$outdir)
}

# Read genes
genes <- read_tsv(args$genes,
                  col_types = cols(species = col_character(),
                                   .default = col_number()))
genes <- genes %>%
  mutate(Genome = str_remove(species, "_genomic$")) %>%
  select(-species) %>%
  replace(is.na(.), 0)

# genes <- genes %>%
#   select(Genome, starts_with("K00"))

# Read features
feats <- read_tsv(args$feats,
                  col_types = cols(Genome = col_character(),
                                   .default = col_number()))
feats <- feats %>%
  rename_all(~str_replace_all(., " ", "."))

# Define constant feats
meta_names <- c("Genome.length", "Number.of.contigs.in.the.genome")

# Filter genes
genes <- genes %>%
  filter(Genome %in% feats$Genome)
gene_summaries <- genes %>%
  pivot_longer(-Genome, names_to = "gene", values_to = "count") %>%
  group_by(gene) %>%
  summarise(n_genomes = sum(count > 0),
            count_sd = sd(count),
            count_iqr = IQR(count),
            min_count = min(count),
            max_count = max(count)) %>%
  ungroup()
gene_summaries %>%
  write_tsv(path = file.path(args$outdir, paste0(args$prefix, "_gene.summaries.txt")))
# hist(gene_summaries$n_genomes)

# Remove constant
genes <- genes %>%
  select(!all_of(gene_summaries$gene[gene_summaries$count_sd == 0]))

# Log transform
feats <- feats %>%
  mutate_at(vars(-Genome), function(x) log10(x + 1))
genes <- genes %>%
  mutate_at(vars(-Genome), function(x) log10(x + 1))
  
if(args$n_pcs > 0){
  cat("Calculating gene PCSs\n")
  gene_pcs <- get_gene_pcs(genes = genes, gene_summaries = gene_summaries,
                           prop_genomes = 0.5, n_pcs = args$n_pcs)
  feats <- feats %>%
    left_join(gene_pcs, by = "Genome")
  
  meta_names <- c(meta_names, paste0("PC", 1:args$n_pcs))
}

# filter genes before test
genes_not_tested <- gene_summaries$gene[gene_summaries$n_genomes < args$min_genomes_test]
genes_not_tested <- intersect(genes_not_tested, colnames(genes))
genes <- genes %>%
  select(!all_of(genes_not_tested))

# Combine all data
dat <- feats %>%
  left_join(genes, by = "Genome") %>%
  select(-Genome)
# dat

# Create variable name vectors
gene_names <- setdiff(colnames(genes), c("Genome", meta_names))
feat_names <- setdiff(colnames(feats), c("Genome", meta_names))

# Linear models
Res <- phyloprof_lm(dat = dat,
                    gene_names = gene_names,
                    feat_names = feat_names,
                    meta_names = meta_names,
                    trans_fun = NULL)
# Res %>%
#   arrange(q.value, p.value) %>%
#   filter(feat == "Number.of.GRINS.detected.in.the.genome") %>%
#   filter(q.value <= 0.05) %>%
#   select(gene) %>%
#   print(n = 100)
Res %>%
  write_tsv(path = file.path(args$outdir, paste0(args$prefix, "_lms.txt")))

# # Log models
# Res <-phyloprof_lm(dat = dat,
#                    gene_names = gene_names,
#                    feat_names = feat_names,
#                    meta_names = meta_names,
#                    trans_fun = function(x) log10(x + 1))
# Res %>%
#   write_tsv(path = file.path(args$outdir, paste0(args$prefix, "_log.lms.txt")))




