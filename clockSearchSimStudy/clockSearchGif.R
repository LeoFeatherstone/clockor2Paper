### Script that moves through the tree and assess support for local clocks
require(ape)
require(phytools)
require(animation)

# Nodes 
getNodes <- function(tree, minCladeSize){
  n <- length(tree$tip.label)
  nodes <- (n+1):(2*n-1)
  tips <- sapply(nodes, function (x) getDescendants(tree, x))
  nodes <- nodes[which(sapply(tips, function(x) length(x)) >= minCladeSize)]
  
  nTips <- sapply(tips, function (x) length(x))
  nTips <- nTips[which(sapply(tips, function(x) length(x)) > minCladeSize)]
  nodes <- nodes[order(nTips, decreasing=T)] # now ordered for col and nesting
  return(nodes)
}

# Add colours to groups
colScheme = RColorBrewer::brewer.pal(3,name = "Set1")[2:1]
colVec <- function(tree, selectedNodes){
  edgeCol <- rep(colScheme[1], length=(2*length(tree$tip.label)-1))
  edgeCol[
    which(tree$edge[,2] %in% getDescendants(tree, selectedNodes[1]))
  ] <- colScheme[2]
#  for (i in 2:length(selectedNodes)){
#    edgeCol[
#      which(tree$edge[,2] %in% getDescendants(tree, selectedNodes[i]))
#      ] <- colScheme[i]
#  }
  return(edgeCol)
}

# Will need further updates to handle n clocks
showAlgorithm2Clocks <- function(nClocks, tree, minCladeSize){
  nodes <- getNodes(tree, minCladeSize = minCladeSize)
  combos <- combn(nodes, m=nClocks-1, simplify=T)
  
  descTips <- adephylo::listTips(tree)
  
  for (i in 1:dim(combos)[2]){
    col <- colVec(tree, selectedNodes = as.vector(combos[,i]))
    # make sure nesting avoids clusters with too few tips
    if (any(table(col)<minCladeSize)) {next}
    # Plotting part.
    par(mfrow=c(1,2))
    date <- as.Date(lubridate::date_decimal(
      as.numeric(gsub(tree$tip.label, pattern = "t.+_", replacement = '')) + 2000
      ))
    
    height <- NELSI::allnode.times(tree, tipsonly=T, keeproot=F)
    #bg=col[-which(tree$edge[,2] %in% tree$edge[,1])]
   plotCols <- rep(colScheme[1], length(tree$tip.label))
   plotCols[getDescendants(tree, as.vector(combos[,i]))] <- colScheme[2]

 #col[sort(tree$edge[,2][-which(tree$edge[,2] %in% tree$edge[,1])])],
    plot(x=date, y=jitter(height), pch=21, bg=plotCols,
         ylab = "Height (subs/site)", xlab = "Date")    
    plot.phylo(tree, show.tip.label = F, edge.color = col)
    
  }

}
## Matching up tip lab node number??

load("./simStudyData.Rdata")
tree <- stem[[40]]
tree$tip.label = gsub(tree$tip.label, pattern = "_Group.+", replacement = "")
tree$edge.length <- sapply(tree$edge.length, function(x) x * runif(0.75, 1.25, n=1))

## Plotting
animation::saveGIF(
  showAlgorithm2Clocks(
    tree, 
    nClocks=2, 
    minCladeSize = 50
  ), 
  movie.name = "clockSearchEg2Clocks.gif", 
  interval=1
)

