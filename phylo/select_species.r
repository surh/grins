library(tidyverse)
library(ape)

# set.seed(12345)
# 
# mat <- cbind(matrix(rnorm(20, mean = 0), nrow = 10),
#              matrix(rnorm(20, mean = 1), nrow = 10),
#              matrix(rnorm(20, mean = 2), nrow = 10),
#              matrix(rnorm(20, mean = 3), nrow = 10),
#              matrix(rnorm(20, mean = 4), nrow = 10))
# row.names(mat) <- LETTERS[1:10]
# dd <- dist(mat)
# dd.clus <- hclust(dd)
# plot(dd.clus)
# dd.phylo <- ape::as.phylo.hclust(dd.clus)
# dd.dendro <- as.dendrogram(dd.clus)

setwd("/godot/users/sur/exp/fraserv/2018/today3/")
mat_file <- "outfile.txt"

dd.mat <- read.table(mat_file, header = TRUE, row.names = 1, sep = "\t")
dd.mat %>% head

dd.mat <- as.matrix(dd.mat)

##########################

N <- 20
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

rm(curr_tree, ind, temp_tree, chosen_leave, l,
   leaves, mat_file, N, new_leaves, prev_dist,
   remaining_leaves, temp_dist, temp_leaves)
############## Use metropolis hastings ###########

# set.seed(12345)
# 
# mat <- cbind(matrix(rnorm(20, mean = 0), nrow = 10),
#              matrix(rnorm(20, mean = 1), nrow = 10),
#              matrix(rnorm(20, mean = 2), nrow = 10),
#              matrix(rnorm(20, mean = 3), nrow = 10),
#              matrix(rnorm(20, mean = 4), nrow = 10))
# row.names(mat) <- LETTERS[1:10]
# dd <- dist(mat)
# dd.mat <- as.matrix(dd)

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

plot(1:iter, Dist)

rm(curr_tree, Leaves, new_tree,
   all_leaves, alpha, burnin, curr_dist, curr_leaves,
   curr_remain_leaves, Dist, i, iter, keep, max_dist,
   max_leaves, min_s, mutated, mutation, N,
   new_dist, new_leaves, new_remain_leaves, s,seed, u)

# 
# 
# dat <- Leaves[9001:10000]
# dat <- apply(sapply(dat, sort),2,paste0,collapse="")
# sort(table(dat))
# max(Dist)
# 
# Dist[ 9000 + which(dat == "BCDF") ]
# Dist[ 9000 + which(dat == "AHIJ") ]
# Dist[ 9000 + which(dat == "BDFH") ]
# 
# Dist[ 9000 + which(dat == "ABDEFGI") ]
# 
# 
# 
