### Simulating trees to test branching regresion
library(ape)
library(devtools)
#install_github("sebastianduchene/NELSI")
library(NELSI)
library(TreeSim)
library(geiger)

## function to get clades with 3-50 tips. I.e. the min possible to regress up to half
getRandClade <- function(tr, minSize, maxSize){
    # num nodes=2*ntips-1
    minNode <- length(tr$tip.label)+1
    maxNode <- 2*length(tr$tip.label)-1
    # search for clade with between 3 and 50 descendants
    found <- FALSE
    while (!found){
        node <- sample((minNode:maxNode), size=1)
        nDesc <- length(geiger::tips(tr, node))
        if (nDesc >= minSize & nDesc <= maxSize){
            clade <- geiger::tips(tr, node)
            found <- TRUE
        }
    }
    return(clade)
}

## Now including Sebastian's NELSI function to simulate local rates
simulate.FLC <- function(tree, params = list(clade.list, stem.clade.indicator, background.rate, local.rates)){
    clade.list = params$clade.list
    stem.clade.indicator = params$stem.clade.indicator
    background.rate = params$background.rate
    local.rates = params$local.rates
        
  check.tip.labels <- function(x){
    if(is.numeric(x)){
      return(x)
    }else if(is.character(x)){
      return(match(x, tree$tip.label))
    }
  }
  
  clade.list <- lapply(clade.list, function(x) check.tip.labels(x))
  if(length(clade.list) != length(local.rates)) stop('The length of local.rates and clade.list must be thesame')
  check.monophyly <- sapply(clade.list, function(x) is.monophyletic(tree, x))
  if(any(!check.monophyly)) stop('At least one of the clades defined is not monophyletic. Please check with is.monophyletic')
  data.matrix <- get.tree.data.matrix(tree)
  data.matrix[, 'branch.rate'] <- background.rate # Note that here all branches are populated and then local clocks are populated below
  
  for(i in 1:length(clade.list)){
    # Steps below need to be done for every local clock
    mrca.node <- get.mrca(tree, clade.list[[i]])
    all.descendants <- get.descending.nodes.branches(tree, mrca.node)
    clade.branches <- all.descendants$descending.branches
    stem.branch <- data.matrix[data.matrix[, 'daughter.node'] == mrca.node, 'branch.index']
    if(stem.clade.indicator[[i]][1]){ # if it applies to the stem
      data.matrix[stem.branch, 'branch.rate'] = local.rates[[i]]
    }
    if(stem.clade.indicator[[i]][2]){ # if it applies to the clade
      data.matrix[clade.branches, 'branch.rate'] = local.rates[[i]]
    }
  }
  
    data.matrix[, 'length.subst'] <- data.matrix[, 'length.time'] * data.matrix[, 'branch.rate']
    tree$edge.length <- data.matrix[, 'length.subst']
    res <- list(tree, data.matrix)
    names(res) <- c('phylogram', 'tree.data.matrix')
    class(res) <- 'ratesim'
    return(res)
}


## Simulating 100 test trees. Using BD for heterochronous trees
set.seed(1234)
trees <- lapply(1:100, function(x) sim.bdsky.stt(n=100, lambdasky=2.5, deathsky=1, sampprobsky=1, timesky=0, timestop=0))
# Add time to tip names
for (i in 1:length(trees)){
	trees[[i]][[1]]$tip.label <- paste0(trees[[i]][[1]]$tip.label, '@', diag(vcv.phylo(trees[[i]][[1]])))
}

## Getting clades between 3 and 50 tips for a second clock
clades <- lapply(trees, function(x) getRandClade(x[[1]], minSize=10, maxSize=50))
# Check dist of clade sizes
#hist(unlist(lapply(clades, function(x) length(x))))

## Generate simple FLC trees with bg rate of 1e-3, and local rate of 2e-3
newTrees <- list()
for (i in 1:length(trees)){
    newTrees[[i]] <- simulate.FLC(trees[[i]][[1]], list(clade.list=clades[i], stem.clade.indicator=list(c(T, T)), background.rate=1, local.rates=2))$phylogram
}
class(newTrees) <- 'multiPhylo'
write.tree(newTrees, file=paste0(getwd(), '/sims/', 'flcStemClade.nwk.tree'))
