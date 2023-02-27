###############################################S######################################
#####################################################################################
######################### Simulation Study for ClockSearch() ########################
#####################################################################################
#####################################################################################
library(ape)
library(TreeSim)
#devtools::install_github("sebastianduchene/NELSI")
library(NELSI)

setwd(paste0(getwd(), "/clockSearchSimStudy/"))

## Get random clade within a specified size. I.e. the min possible to half of the num tips.
## Returns vector of tips comprising clade
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



## Simulating 100 test trees. Using BD for heterochronous trees. Height added to 2000.
## Using a Unif[100, 1000] for number of tips
trees <- list()
class(trees) <- "multiPhylo"

set.seed(1234)
for (i in 1:100) {

  trees[i] <- sim.bdsky.stt(n=250,
    lambdasky=2.5,
    deathsky=1,
    sampprobsky=1,
    timesky=0,
    timestop=0)

}

# Add time to tip names
for (i in seq_along(trees)){

  trees[[i]]$tip.label <- paste0(
    trees[[i]]$tip.label,
    '_',
    diag(vcv.phylo(trees[[i]]))
    )

}

## Getting clades between between 125 and 375 tips (25%-75% of tips)
clades <- lapply(trees, 
  function(x) getRandClade(
    x, 
    minSize = 50, 
    maxSize = 150
    )
  )

for (i in seq_along(clades)){

  match <- which(trees[[i]]$tip.label %in% clades[[i]])
  clades[[i]] <- paste0(clades[[i]], "_Group1")
  trees[[i]]$tip.label[match] <- paste0(trees[[i]]$tip.label[match], "_Group1")
  trees[[i]]$tip.label[-match] <- paste0(trees[[i]]$tip.label[-match], "_Group2")

}


## Generate simple FLC trees with global rate of 1e-3 and local rate of 5e-3
stemClade <- list()
stem <- list()
class(stemClade) <- "multiPhylo"
class(stem) <- "multiPhylo"

for (i in 1:length(trees)){

    stemClade[[i]] <- NELSI::simulate.FLC(
      trees[[i]], 
      list(
        clade.list=clades[i], 
        stem.clade.indicator=list(c(T, T)), 
        background.rate=0.001, 
        local.rates=0.005))$phylogram

    stem[[i]] <- NELSI::simulate.FLC(
      trees[[i]], 
      list(
        clade.list=clades[i], 
        stem.clade.indicator=list(c(T, F)), 
        background.rate=0.001, 
        local.rates=0.005))$phylogram
}

dataStemClade <- data.frame()
dataStem <- data.frame()
# wrapper function to use in testing clockSearch()
testClockSearch <- function(maxClocks, testTrees, baseTrees, trueClades) {
  data <- data.frame()

  for (i in seq_along(baseTrees)){
    cmd <- paste0("'", write.tree(testTrees[[i]]), "' ", 1, ' ', maxClocks)
    op <- system(paste('node ./clockSearchWrapper.js', cmd), intern = TRUE)
    
    inferredNumClocks <- as.numeric(op)
    grp <- read.table("./tmpTips.txt", header = FALSE)

    pcMatch <- max(
      (length(which(grp$V1 %in% trueClades[[i]])) / length(grp$V1)),
      (length(
        which(grp$V1 %in% baseTrees[[i]]$tip.label[-(which(baseTrees[[i]]$tip.label %in% trueClades[[i]]))])) / length(grp$V1))
    )

    data <- rbind(
      data,
      c(maxClocks, op, pcMatch)
    )
  }

  colnames(data) <- c("maxClocks", "numGroups", "pcMatch")
  system("rm ./tmpTips.txt")

  return(data)
}

testMaxClocks <- c(2:4)

dataStemClade <- data.frame()
dataStem <- data.frame()

for (c in testMaxClocks) {
  dataStemClade <- rbind(
    dataStemClade,
    testClockSearch(
      maxClocks = c,
      testTrees = stemClade,
      baseTrees = trees,
      trueClades = clades
    )
  )
  print(paste(c, "Done"))
}

for (c in testMaxClocks) {
  dataStem <- rbind(
    dataStem,
    testClockSearch(
      maxClocks = c,
      testTrees = stem,
      baseTrees = trees,
      trueClades = clades
    )
  )
  print(paste(c, "Done"))
}