#!/usr/bin/env Rscript

# (C) Copyright 2019 Sur Herrera Paredes

library(tidyverse)
library(argparser)

#' Plot statistics from checkm
#'
#' @param x Checkm results table as tibble or data frame
#' @param completeness threshold
#' @param contamination theshold
#' @param outdir 
#' @param prefix 
#' @param breaks1 breaks for contamination and completeness
#' @param breaks2 breaks for strain heterogeneity
#'
#' @return
#' @export
#'
#' @examples
plot_checkm_results <- function(x, completeness = 98, contamination = 2,
                                outdir = "./", prefix = NULL,
                                breaks1 = c(0, 80, 90,95,96,97,98,99,100),
                                breaks2 =  c(0, 1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)){
  # Plot main statistics using custom breakpoints
  # breaks <- c(0, 80, 90,95,96,97,98,99,100)
  
  # Completeness
  p1 <- ggplot(tibble(Completeness = cut(Res$Completeness,
                                         breaks1,
                                         include.lowest = TRUE)),
               aes(x = Completeness)) +
    geom_histogram(stat = "count") +
    scale_y_log10() +
    theme_classic()
  # p1
  filename <- paste0(outdir, "/", prefix, "completeness_histogram.svg")
  ggsave(filename, p1, width = 6, height = 4)
  
  # Contamination
  p1 <- ggplot(tibble(Contamination = cut(Res$Contamination,
                                          100 - breaks1,
                                          include.lowest = TRUE, right = FALSE)),
               aes(x = Contamination)) +
    geom_histogram(stat = "count") +
    scale_y_log10() +
    theme_classic()
  # p1
  filename <- paste0(outdir, "/", prefix, "contamination_histogram.svg")
  ggsave(filename, p1, width = 6, height = 4)
  # Res %>% filter(Contamination > 100)
  
  # Strain heterogeneity
  # breaks <- c(0, 1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)
  p1 <- ggplot(tibble(Strain.heterogeneity = cut(Res$`Strain heterogeneity`,
                                                 breaks2,
                                                 include.lowest = TRUE)),
               aes(x = Strain.heterogeneity)) +
    geom_histogram(stat = "count") +
    scale_y_log10() +
    theme_classic()
  # p1
  filename <- paste0(outdir, "/", prefix, "heterogeneity_histogram.svg")
  ggsave(filename, p1, width = 6, height = 4)
  
  # Relationship between the three
  p1 <- ggplot(Res, aes(x = Completeness, y = Contamination)) +
    geom_point(alpha = 0.2) +
    # geom_point(aes(col = `Strain heterogeneity`)) +
    # scale_colour_continuous(trans='log2') +
    scale_y_log10() +
    # scale_x_sqrt() +
    geom_hline(yintercept = contamination, col = "magenta") +
    geom_vline(xintercept = completeness, col = "magenta") +
    theme_classic()
  # p1
  filename <- paste0(outdir, "/", prefix, "completeness_vs_contamination.svg")
  ggsave(filename, p1, width = 6, height = 6)
}

#' Process argumemtns
#'
#' @return
#' @export
#'
#' @examples
process_arguments <- function(){
  p <- arg_parser(paste0("Reads an output table from checkM, plots the resutlts, ",
                         "and selects genomes to keep based on thresholds"))
  
  p <- add_argument(p, "file",
                    help = paste0("CheckM results tab-delimited file"),
                    type = "character")
  p <- add_argument(p, "--outdir",
                    help = paste0("Directory to write (and create if needed) to write ",
                                  "the output. Must not exist."),
                    default = "output/",
                    type = "character")
  p <- add_argument(p, "--completeness",
                    help = paste0("Minumum completeness for a genome to pass the threshold."),
                    default = 98,
                    type = "integer")
  p <- add_argument(p, "--contamination",
                    help = paste0("Maximum contamination for a genome to pass the threshold"),
                    type = "integer",
                    default = 2)
  p <- add_argument(p, "--heterogeneity",
                    help = paste0("Maximim heterogeneity for a genome to pass the threshold"),
                    type = "integer",
                    default = 0)
  
  
  # Read arguments
  cat("Processing arguments...\n")
  args <- parse_args(p) 
  
  # Process arguments
  #cat("\t", args$heterogeneity, "\n")
  
  return(args)
}

#############################
# args <- list(file = "checkm_results.txt",
#              outdir = "output/",
#              completeness = 98,
#              contamination = 2,
#              heterogeneity = NULL)
args <- process_arguments()
# print(args)

# Read data and prepater output directory
cat("Reading data...\n")
Res <- read_tsv(file = args$file, col_types = 'ccnnnnnnnnnnnn')
if(dir.exists(args$outdir)){
  stop("ERROR: output directory already exists")
}else{
  dir.create(args$outdir)
}

# Plot overall distribution
cat("Plotting results...\n")
plot_checkm_results(x = Res, completeness = args$completeness,
                    contamination = args$contamination,
                    outdir = args$outdir, prefix = "all_")

# Filter
cat("Filtering...\n")
n_genomes <- nrow(Res)
Res <- Res %>% filter(Completeness >= args$completeness & Contamination <= args$contamination)
if(!is.na(args$heterogeneity)){
  Res <- Res %>% filter(`Strain heterogeneity` <= args$heterogeneity)
}
n_passed <- nrow(Res)
n_discarded <- n_genomes - n_passed

# Plot final selection
cat("Plotting final selection...\n")
plot_checkm_results(x = Res, completeness = args$completeness,
                    contamination = args$contamination,
                    outdir = args$outdir, prefix = "chosen_")

# Write results
cat("Writing results...\n")
filename <- paste0(args$outdir, "/chosen_checkm_results.txt")
write_tsv(Res, filename)

cat(n_genomes, " genomes analyzed by checkM.\n")
cat(n_discarded, " genomes discarded.\n")
cat(n_passed, " genomes passed the thresholds.\n")
