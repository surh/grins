# Copyright (C) 2018 Sur Herrera Paredes
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

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
dd.nj <- nj(as.dist(dd.mat))
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

dd.gg <- as.treedata(dd.nj) %>% 
  full_join(tibble(label=dd.nj$tip.label) %>%
              mutate(Selected = label %in% leaves),
            by = "label")
dd.gg <- groupOTU(dd.gg, .node = leaves)

p1 <- ggtree(dd.gg,
             aes(color=group, size = group)) +
  scale_size_manual(values = c(0.2,1)) +
  geom_tiplab()
p1
ggsave("max_tree_greedy.png", p1, width = 6, height = 18)
