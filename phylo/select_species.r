library(tidyverse)
library(ape)


set.seed(12345)

mat <- cbind(matrix(rnorm(20, mean = 0), nrow = 10),
             matrix(rnorm(20, mean = 1), nrow = 10),
             matrix(rnorm(20, mean = 2), nrow = 10),
             matrix(rnorm(20, mean = 3), nrow = 10),
             matrix(rnorm(20, mean = 4), nrow = 10))
row.names(mat) <- LETTERS[1:10]
dd <- dist(mat)
# dd.clus <- hclust(dd)
# plot(dd.clus)
# dd.phylo <- ape::as.phylo.hclust(dd.clus)
# dd.dendro <- as.dendrogram(dd.clus)


dd.mat <- as.matrix(dd)

N <- 4
leaves <- NULL

# Find indices
ind <- arrayInd(which(dd.mat == max(dd.mat)), dim(dd.mat))         
ind <- ind[ ind[ ,1 ] < ind[ ,2 ], , drop = FALSE]
cat(ind,"\n")

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

write_tsv(as.matrix(leaves), "selected_leaves_gridy")

dd.phylo <- ape::as.phylo.hclust(dd.clus)
plot(dd.phylo)
plot(curr_tree)

rm(list = ls())
############## Use metropolis hastings ###########

set.seed(12345)

mat <- cbind(matrix(rnorm(20, mean = 0), nrow = 10),
             matrix(rnorm(20, mean = 1), nrow = 10),
             matrix(rnorm(20, mean = 2), nrow = 10),
             matrix(rnorm(20, mean = 3), nrow = 10),
             matrix(rnorm(20, mean = 4), nrow = 10))
row.names(mat) <- LETTERS[1:10]
dd <- dist(mat)
dd.mat <- as.matrix(dd)

N <- 7
s <- 2
iter <- 10000
burnin <- 10
keep <- 2

set.seed(12345)

all_leaves <- colnames(dd.mat)
# Initialize
curr_leaves <- sample(length(all_leaves), size = N, replace = FALSE)
curr_remain_leaves <- all_leaves[-curr_leaves]
curr_leaves <- all_leaves[curr_leaves]
curr_leaves
curr_tree <- as.phylo(hclust(as.dist(dd.mat[curr_leaves, curr_leaves])))
curr_dist <- sum(curr_tree$edge.length)

Dist <- c(curr_dist)
Leaves <- list(curr_leaves)
for(i in 1:(iter-1)){
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

plot(1:iter, Dist)



dat <- Leaves[9001:10000]
dat <- apply(sapply(dat, sort),2,paste0,collapse="")
sort(table(dat))
max(Dist)

Dist[ 9000 + which(dat == "BCDF") ]
Dist[ 9000 + which(dat == "AHIJ") ]
Dist[ 9000 + which(dat == "BDFH") ]

Dist[ 9000 + which(dat == "ABDEFGI") ]



