library(tidyverse)
library(ape)

# setwd("/godot/users/sur/exp/fraserv/2018/today3/")
mat_file <- "outfile.txt"

dd.mat <- read.table(mat_file, header = TRUE, row.names = 1, sep = "\t")
#dd.mat %>% head

dd.mat <- as.matrix(dd.mat)

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
curr_tree <- as.phylo(hclust(as.dist(dd.mat[leaves, leaves])))
cat(sum(curr_tree$edge.length), "\n")

while(length(leaves) < N){
  cat("===============\n")
  prev_dist <- sum(curr_tree$edge.length)
  chosen_leave <- NULL
  for(l in remaining_leaves){
    temp_leaves <- c(leaves, l)
    temp_tree <- as.phylo(hclust(as.dist(dd.mat[temp_leaves, temp_leaves])))
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
  curr_tree <- as.phylo(hclust(as.dist(dd.mat[leaves, leaves])))
  cat(sum(curr_tree$edge.length), "\n")
  
  cat("-------------\n")
}

write_tsv(as.tibble(leaves), "selected_leaves_greedy.txt")

dd.phylo <- ape::as.phylo.hclust(dd.clus)
plot(dd.phylo)
plot(curr_tree)
