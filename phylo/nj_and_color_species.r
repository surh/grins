library(tidyverse)
library(ape)
library(ggtree)
library(treeio)

set.seed(12345)

mat <- cbind(matrix(rnorm(20, mean = 0), nrow = 10),
             matrix(rnorm(20, mean = 1), nrow = 10),
             matrix(rnorm(20, mean = 2), nrow = 10),
             matrix(rnorm(20, mean = 3), nrow = 10),
             matrix(rnorm(20, mean = 4), nrow = 10))
row.names(mat) <- LETTERS[1:10]
dd <- dist(mat)

dd.nj <- nj(dd)
plot(dd.nj)

leaves <- c('C','B','D','H')

dd.gg <- as.treedata(dd.nj)

# dat <- tibble(label = dd.gg@phylo$tip.label) %>%
#   mutate(Selected = label %in% leaves)
# dat

p1 <- ggtree(groupOTU(dd.gg, .node = leaves),
             aes(color=group, size = group)) +
  scale_size_manual(values = c(0.2,1)) +
  geom_tiplab()
p1
