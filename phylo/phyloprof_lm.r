#!/usr/bin/env Rscript
# setwd("/cashew/users/sur/exp/fraserv/2020/today2")
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

# Remove constant
genes <- genes %>%
  select(!all_of(gene_summaries$gene[gene_summaries$count_sd == 0]))

# Combine all data
dat <- feats %>%
  left_join(genes, by = "Genome")

# Create variable name vectors
gene_names <- setdiff(colnames(genes), "Genome")
feat_names <- setdiff(colnames(feats), "Genome")

# Fit
models <- expand_grid(gene_name = gene_names, feat_name = feat_names)
cat("Fitting models\n")
Res <- models %>%
  pmap_dfr(function(gene_name, feat_name, dat = dat){
    d <- dat[,c(gene_name, feat_name)]
    colnames(d) <- c("gene", "feat")
    # f1 <- paste(feat_name, "~", gene_name)
    f1 <- formula(feat ~ gene)
    tidy(lm(f1, data = d)) %>%
      mutate(feat = feat_name,
             gene = gene_name) %>%
      filter(term != "(Intercept)") %>%
      select(-term)
  }, dat = dat)

cat("Fitting log models\n")
Res.log <- models %>%
  pmap_dfr(function(gene_name, feat_name, dat = dat){
    d <- dat[,c(gene_name, feat_name)]
    colnames(d) <- c("gene", "feat")
    # f1 <- paste("log10(",feat_name,"+ 1 )", "~", gene_name)
    f1 <- formula(log(feat + 1) ~ gene)
    tidy(lm(f1, data = d)) %>%
      mutate(feat = feat_name,
             gene = gene_name) %>%
      filter(term != "(Intercept)") %>%
      select(-term)
    }, dat = dat)

# Write results
Res %>%
  write_tsv(path = file.path(args$outdir, paste0(args$prefix, "_lms.txt")))

Res.log %>%
  write_tsv(path = file.path(args$outdir, paste0(args$prefix, "_log.lms.txt")))

gene_summaries %>%
  write_tsv(path = file.path(args$outdir, paste0(args$prefix, "_gene.summaries.txt")))

