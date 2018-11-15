library(tidyverse)
library(ape)
library(ggtree)
library(treeio)
# setwd("/godot/users/sur/exp/fraserv/2018/today3/")

# Read data
mat_file <- "outfile.txt"
dd.mat <- read.table(mat_file, header = TRUE, row.names = 1, sep = "\t")

# Convert matrix to nj tree
dd.mat <- as.matrix(dd.mat)
dd.nj <- nj(dd)
plot(dd.nj)


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
curr_tree <- keep.tip(dd.nj, tip = leaves)
curr_dist <- sum(curr_tree$edge.length)

max_leaves <- curr_leaves
max_dist <- curr_dist
Dist <- c(curr_dist)
Leaves <- list(curr_leaves)
for(i in 2:iter){
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
  new_tree <- keep.tip(dd.nj, tip = new_leaves)
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

dd.gg <- as.treedata(dd.nj)
p1 <- ggtree(groupOTU(dd.gg, .node = max_leaves),
             aes(color=group, size = group)) +
  scale_size_manual(values = c(0.2,1)) +
  geom_tiplab()
p1
ggsave("max_tree.png", width = 6, height = 18)

