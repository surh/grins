#!/usr/bin/env Rscript
# setwd("/cashew/users/sur/exp/fraserv/2020/today3")
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
        mutate(feat = feat_name,
               gene = gene_name) %>%
        filter(term == "gene") %>%
        select(-term)
    }, dat = dat, meta_names = meta_names)
}


args <- process_arguments()
# args <- list(genes = "/cashew/shared_data/grins/2020-05-12.eggnot_annots/tabs/KEGG_KOs.txt.gz",
#              feats = "/cashew/users/nivina/2020-03-24.Streptomyces_GRINS_detection/GRINS_detected_in_genomes_and_BGCs.txt",
#              outdir = "output",
#              prefix = "KEGG_KOs")


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
#   select(Genome, starts_with("K000"))

# Read features
feats <- read_tsv(args$feats,
                  col_types = cols(Genome = col_character(),
                                   .default = col_number()))
feats <- feats %>%
  rename_all(~str_replace_all(., " ", "."))

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

# Remove constant
genes <- genes %>%
  select(!all_of(gene_summaries$gene[gene_summaries$count_sd == 0]))

# Combine all data
dat <- feats %>%
  left_join(genes, by = "Genome") %>%
  select(-Genome)
dat

# Create variable name vectors
meta_names <- c("Genome.length", "Number.of.contigs.in.the.genome")
gene_names <- setdiff(colnames(genes), c("Genome", meta_names))
feat_names <- setdiff(colnames(feats), c("Genome", meta_names))

# Linear models
Res <- phyloprof_lm(dat = dat,
                    gene_names = gene_names,
                    feat_names = feat_names,
                    meta_names = meta_names,
                    trans_fun = NULL)
Res %>%
  write_tsv(path = file.path(args$outdir, paste0(args$prefix, "_lms.txt")))

# Log models
Res <-phyloprof_lm(dat = dat,
                   gene_names = gene_names,
                   feat_names = feat_names,
                   meta_names = meta_names,
                   trans_fun = function(x) log10(x + 1))
Res %>%
  write_tsv(path = file.path(args$outdir, paste0(args$prefix, "_log.lms.txt")))




