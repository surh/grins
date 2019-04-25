library(ggplot2)
library(ggtree)

# Read tree
tree <- read.tree(file = "RAxML_bipartitions.example")
tree

# Make ultrametric. Requires ape (installed automatically when ggtree installed)
tree <- ape::chronopl(tree, lambda = 1)

# Plot
mytree.plot <- ggtree(tree, ladderize = TRUE, layout = "rectangular") + 
  geom_tiplab() +
  geom_nodelab()
mytree.plot

# Create some fake meta data
dat <- data.frame(GRINS = rep(c("GRINS", "nonGRINS"), each = 5),
                  cluster = rep(c("pks1", "pks2", "pks3", "pks4", "pks5"), length.out = 10),
                  stringsAsFactors = FALSE)
row.names(dat) <- tree$tip.label
dat

# Use facet plot to plot different variables
# Need to create one data frame per value
# First add GRINS column
d.temp <- data.frame(id = row.names(dat), GRINS = dat$GRINS, stringsAsFactors = FALSE)
# Create color colum
d.temp$fill <- d.temp$GRINS
d.temp
p1 <- facet_plot(mytree.plot, panel = "GRINS", data=d.temp,
                 geom=geom_tile, mapping = aes(x = 1, fill=GRINS))
p1

# Second add a second pannel
d.temp <- data.frame(id = row.names(dat), cluster = dat$cluster)
p1 <- facet_plot(p1, panel = "cluster", data=d.temp,
                 geom=geom_tile, mapping = aes(x = 1, fill=cluster))
p1
