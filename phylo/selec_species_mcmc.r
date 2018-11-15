library(tidyverse)
library(ape)

# setwd("/godot/users/sur/exp/fraserv/2018/today3/")
mat_file <- "outfile.txt"

dd.mat <- read.table(mat_file, header = TRUE, row.names = 1, sep = "\t")
#dd.mat %>% head

dd.mat <- as.matrix(dd.mat)


# MCMC

seed <- 12345
N <- 100
s <- 50
min_s <- 3
iter <- 10000
burnin <- 10
keep <- 2

set.seed(seed)

all_leaves <- colnames(dd.mat)
# Initialize
curr_leaves <- sample(length(all_leaves), size = N, replace = FALSE)
curr_remain_leaves <- all_leaves[-curr_leaves]
curr_leaves <- all_leaves[curr_leaves]
curr_leaves
curr_tree <- as.phylo(hclust(as.dist(dd.mat[curr_leaves, curr_leaves])))
curr_dist <- sum(curr_tree$edge.length)

max_leaves <- curr_leaves
max_dist <- curr_dist
Dist <- c(curr_dist)
Leaves <- list(curr_leaves)
for(i in 1:(iter-1)){
  # Linear simulated annealing
  s <- ceiling(s * (iter-i) / iter)
  s <- max(s,min_s)
  
  # Mutate
  mutated <- sample(N, size = s, replace = FALSE)
  new_remain_leaves <- c(curr_leaves[mutated], curr_remain_leaves)
  new_leaves <- curr_leaves[-mutated]
  mutation <- sample(length(new_remain_leaves), size = s, replace = FALSE)
  new_leaves <- c(new_leaves, new_remain_leaves[mutation])
  new_remain_leaves <- new_remain_leaves[-mutation]
  
  # Check
  new_tree <- as.phylo(hclust(as.dist(dd.mat[new_leaves, new_leaves])))
  new_dist <- sum(new_tree$edge.length)
  
  if(new_dist > max_dist){
    max_dist <- new_dist
    max_leaves <- new_leaves
    cat(max_leaves, "\n")
    cat("---------------\n")
  }
  
  alpha <- new_dist / curr_dist
  u <- runif(1)
  if(u <= alpha){
    curr_tree <- new_tree
    curr_dist <- new_dist
    
    curr_leaves <- new_leaves
    curr_remain_leaves <- new_remain_leaves
  }
  
  Dist <- c(Dist, curr_dist)
  Leaves <- c(Leaves, list(curr_leaves))
}

write_tsv(as.tibble(max_leaves), "selected_leaves_mcmc.txt")
max_leaves
max_dist

png("mcmc_dist.png", width = 1000, height = 800)
plot(1:iter, Dist)
dev.off()
