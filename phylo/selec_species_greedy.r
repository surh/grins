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

# Greedy
N <- 100
leaves <- NULL

# Find indices
ind <- arrayInd(which(dd.mat == max(dd.mat)), dim(dd.mat))         
ind <- ind[ ind[ ,1 ] < ind[ ,2 ], , drop = FALSE]
# cat(ind,"\n")

# Find Leaves
new_leaves <- colnames(dd.mat)[ind[1,]]
cat(new_leaves, "\n")

leaves <- c(leaves, new_leaves)
cat(leaves, "\n")
remaining_leaves <- setdiff(colnames(dd.mat), leaves)

# Update tree
curr_tree <-  keep.tip(dd.nj, tip = leaves)
cat(sum(curr_tree$edge.length), "\n")

while(length(leaves) < N){
  cat("===============\n")
  prev_dist <- sum(curr_tree$edge.length)
  chosen_leave <- NULL
  for(l in remaining_leaves){
    temp_leaves <- c(leaves, l)
    temp_tree <- keep.tip(dd.nj, tip = temp_leaves)
    temp_dist <- sum(temp_tree$edge.length)
    
    if(temp_dist > prev_dist){
      chosen_leave <- l
      prev_dist <- temp_dist
    }
  }
  
  # Find Leaves
  new_leaves <- chosen_leave
  cat(new_leaves, "\n")
  
  leaves <- c(leaves, new_leaves)
  cat(leaves, "\n")
  remaining_leaves <- setdiff(colnames(dd.mat), leaves)
  
  # Update tree
  curr_tree <- keep.tip(dd.nj, tip = leaves)
  cat(sum(curr_tree$edge.length), "\n")
  
  cat("-------------\n")
}

write_tsv(as.tibble(leaves), "selected_leaves_greedy.txt")

# plot(curr_tree)

dd.gg <- as.treedata(curr_tree)
p1 <- ggtree(groupOTU(dd.gg, .node = leaves),
             aes(color=group, size = group)) +
  scale_size_manual(values = c(0.2,1)) +
  geom_tiplab()
p1
ggsave("max_tree_greedy.png", width = 6, height = 18)