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
                  GCSkew = runif(10, -1, 1),
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

# Add a second panel with cluster annotation
d.temp <- data.frame(id = row.names(dat), cluster = dat$cluster)
p1 <- facet_plot(p1, panel = "cluster", data=d.temp,
                 geom=geom_tile, mapping = aes(x = 1, fill=cluster))
p1

# Add a third panel with cluster label
d.temp <- data.frame(id = row.names(dat), cluster = dat$cluster)
p1 <- facet_plot(p1, panel = "cluster.label", data=d.temp,
                 geom=geom_text, mapping = aes(x = 1, label=cluster))
p1

# Add a fourth panel with GC skew
d.temp <- data.frame(id = row.names(dat), GCSkew = dat$GCSkew)
p1 <- facet_plot(p1, panel = "GCSkew", data=d.temp,
                 geom=geom_segment, mapping = aes(x = 0, xend = GCSkew, y = y, yend = y), col = "red", stat = "identity")
p1

