### Script that moves through the tree and assess support for local clocks
require(ape)
require(phytools)
require(animation)
tree <- rtree(20, rooted = T)
plot.phylo(tree, show.tip.label=T)
nodelabels(cex=2)


tiplabels(cex=1.5)
tree$tip.label

minCladeSize = 3

# now, plot subclades of tree
getNodes <- function(tree, minCladeSize){
    n <- length(tree$tip.label)
    nodes <- (n+1):(2*n-1)
    tips <- sapply(nodes, function (x) getDescendants(tree, x))
    nodes <- nodes[which(sapply(tips, function(x) length(x)) > minCladeSize)]

    nTips <- sapply(tips, function (x) length(x))
    nTips <- nTips[which(sapply(tips, function(x) length(x)) > minCladeSize)]
    nodes <- nodes[order(nTips, decreasing=T)] # now ordered for colouring adn nesting
    return(nodes)
}

# set num Clocks = 2, print trees. Later, enter functionality to nest colours
# expecing combn nodes
# does new colour need to be done
colVec <- function(tree, selectedNodes){
    edgeCol <- rep(1, length=(2*length(tree$tip.label)-1))
    for (i in 1:length(selectedNodes)){
      edgeCol[which(tree$edge[,2] %in% getDescendants(tree, selectedNodes[i]))] <- i+1
    }
    #edgeCol[which(tree$edge[,2] %in% getDescendants(tree, selectedNode))] <- 'red'
    return(edgeCol)
}

plot(1:10, 1:10, pch=16, cex=3, col=1:8)

tree <- read.tree("./sims/flcStem.nwk.tree")[[1]]
#nodes <- getNodes(tree, minCladeSize = 3)

# attempting to plot now with generalised functions
# Need to include likelihood calculation now

showAlgorithm <- function(nClocks, tree, minCladeSize){
  nodes <- getNodes(tree, minCladeSize = minCladeSize)
  combos <- combn(nodes, m=nClocks, simplify=T)

  for (i in 1:dim(combos)[2]){
    col <- colVec(tree, selectedNodes = as.vector(combos[,i]))
    # test to make sure nesting avoids clusters with too many tips
    if (any(table(col)<minCladeSize)) {next}
    # Plotting part. Give titles for EG
    par(mfrow=c(1,2))
    date <- as.numeric(gsub(tree$tip.label, pattern = ".+@", replacement = ''))
    height <- NELSI::allnode.times(tree, tipsonly=T, keeproot=F)

    plot(x=date, y=height, pch=21, bg=col[-which(tree$edge[,2] %in% tree$edge[,1])])    
    plot.phylo(tree, show.tip.label = F, edge.color = col)
    
  }
  print("PLOT TEST DONE")
}

#####################################################
################## TEST TIME WITH PLOTTING ##########
#####################################################
t1 <- Sys.time()
showAlgorithm(tree, nClocks=3, minCladeSize = 20)
t2 <- Sys.time()

#####################################################
################## TEST TIME WITHOUT PLOTTING #######
#####################################################
showAlgorithm <- function(nClocks, tree, minCladeSize){
  nodes <- getNodes(tree, minCladeSize = minCladeSize)
  combos <- combn(nodes, m=nClocks, simplify=T)

  for (i in 1:dim(combos)[2]){
    col <- colVec(tree, selectedNodes = as.vector(combos[,i]))
    # test to make sure nesting avoids clusters with too many tips
    if (any(table(col)<minCladeSize)) {next}
    # Plotting part. Give titles for EG
    #par(mfrow=c(1,2))
    date <- as.numeric(gsub(tree$tip.label, pattern = ".+@", replacement = ''))
    height <- NELSI::allnode.times(tree, tipsonly=T, keeproot=F)

    #plot(x=date, y=height, pch=21, bg=col[-which(tree$edge[,2] %in% tree$edge[,1])])    
    #plot.phylo(tree, show.tip.label = F, edge.color = col)
    
  }
  print("NO-PLOT TEST DONE")
}

t3 <- Sys.time()
showAlgorithm(tree, nClocks=3, minCladeSize = 20)
t4 <- Sys.time()

print(c("WITH PLOT: ", t2-t1, "WITHOUT PLOT: ", t4-t3))

animation::saveGIF(showAlgorithm(tree, nClocks=4, minCladeSize = 20), movie.name = "algoEg.gif", interval=1)

animation::saveGIF(showAlgorithm(tree, nClocks=2, minCladeSize = 10), movie.name = "algoEg2Clocks.gif", interval=1)