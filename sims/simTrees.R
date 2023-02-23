# script smulates trees for test data
library(ape)
library(devtools)
#install_github("sebastianduchene/NELSI")
library(NELSI)
library(TreeSim)
library(geiger)
library(phangorn)

## function to get clades within a specified size. I.e. the min possible to half of the num tips.
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


## Simulating 100 test trees. Using BD for heterochronous trees. Height added to 2000.
## Using a Unif[100, 1000] for number of tips
trees <- list()
class(trees) <- "multiPhylo"

set.seed(1234)
for (i in 1:100) {
  trees[i] <- sim.bdsky.stt(n=ceiling(runif(100, 1000, n = 1)),
    lambdasky=2.5,
    deathsky=1,
    sampprobsky=1,
    timesky=0,
    timestop=0)
}

# Add time to tip names
for (i in seq_along(trees)){
  trees[[i]]$tip.label <- paste0(trees[[i]]$tip.label,
    '_',
    as.Date(
      date_decimal(
        diag(vcv.phylo(trees[[i]])) + 2000)))
}

## Getting clades between 25 and floor(nTips / 2) tips for a second clock
clades <- lapply(trees, 
  function(x) getRandClade(
    x, 
    minSize=25, 
    maxSize = floor(length(x$tip.label) / 2)))

for (i in seq_along(clades)){
  match <- which(trees[[i]]$tip.label %in% clades[[i]])
  clades[[i]] <- paste0(clades[[i]], "_C1")
  trees[[i]]$tip.label[match] <- paste0(trees[[i]]$tip.label[match], "_C1")
  trees[[i]]$tip.label[-match] <- paste0(trees[[i]]$tip.label[-match], "_C0")
}
# TODO: Modify tip labs to add group name

save(trees, file='baseTrees.RData')
# Check dist of clade sizes
#hist(unlist(lapply(clades, function(x) length(x))))

## Generate simple FLC trees with bg rate of 1e-3, and local rate of 5e-3
newStemClade <- list()
newStem <- list()
newClade <- list()

for (i in 1:length(trees)){
    newStemClade[[i]] <- simulate.FLC(trees[[i]], list(clade.list=clades[i], stem.clade.indicator=list(c(T, T)), background.rate=0.001, local.rates=0.005))$phylogram
    newStem[[i]] <- simulate.FLC(trees[[i]], list(clade.list=clades[i], stem.clade.indicator=list(c(T, F)), background.rate=0.001, local.rates=0.005))$phylogram
    newClade[[i]] <- simulate.FLC(trees[[i]], list(clade.list=clades[i], stem.clade.indicator=list(c(F, T)), background.rate=0.001, local.rates=0.005))$phylogram
    # NB, the notional 1,1 element in ths treatment matrix is the regular trees themselves
    print(i)
}

# saving trees
for (i in 1:length(trees)){
  write.tree(newStemClade[[i]], file=paste0(getwd(), '/sims/', "stemClade", i, ".nwk"))
  write.tree(newStem[[i]], file=paste0(getwd(), '/sims/', "stem", i, ".nwk"))
  write.tree(newClade[[i]], file=paste0(getwd(), '/sims/', "clade", i, ".nwk"))
}

class(trees) <- 'multiPhylo'
write.tree(trees, file=paste0(getwd(), '/sims/', 'baseTrees.nwk.tree')) # doesn't want to work right now

class(newStem) <- 'multiPhylo'
write.tree(newStem, file=paste0(getwd(), '/sims/', 'flcStem.nwk.tree'))

class(newClade) <- 'multiPhylo'
write.tree(newClade, file=paste0(getwd(), '/sims/', 'flcClade.nwk.tree'))

class(newStemClade) <- 'multiPhylo'
write.tree(newStemClade, file=paste0(getwd(), '/sims/', 'flcStemClade.nwk.tree'))
