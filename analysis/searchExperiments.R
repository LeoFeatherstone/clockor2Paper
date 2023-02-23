# toying with code to compare all subtrees
library(ape)
library(devtools)
install_github("YuLab-SMU/ggtree")
library(ggtree)

# trying clade and stem data
trees <- read.tree("./sims/flcStem.nwk.tree")
load("./flcStemClade.RData")

t <- trees[[1]]
c <- clades[[1]]

sub <- subtrees(t)

# start out plot animation
plot.phylo(t, show.tip.label = F)
edge.color(t, groups = sub[[20]]$tip.label)         

# looking into edge labels 24-08-22
library(ape)
t <- rtree(8, rooted = T)
plot.phylo(t, root.edge = T, show.tip.label=F)
nodelabels(cex=2)
tiplabels(cex=2)
