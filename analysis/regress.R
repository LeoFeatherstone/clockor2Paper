# Will start simple with regressing one tree from one sim type: flcStem
require(ape)
require(NELSI)
require(tidyverse)

# function to give AIC for >= 1 clades
# let clades be a list of vectors in future. A vector for now.
getIC <- function(tree, clades){
    date <- as.numeric(gsub(tree$tip.label, pattern = ".+@", replacement = ''))
    height <- allnode.times(tree, tipsonly=T, keeproot=F)
    cladeIndex <- which(tree$tip.label %in% clades)

    k <- 2
    ksplit <- 4 # include argument for length of clades later
    n <- length(tree$tip.label)
    likWhole <- as.numeric(logLik(lm(height~date)))
    lik1 <- as.numeric(logLik(lm(height[cladeIndex]~date[cladeIndex])))
    lik2 <- as.numeric(logLik(lm(height[-cladeIndex]~date[-cladeIndex])))

    # standard AIC
    wholeAIC <- 2*k-2*likWhole # assume 2 estimated parms (k=2)
    splitAIC <- 2*k-2*lik1 + 2*k-2*lik2

    # AICc
    wholeAICc <- (2*(k^2)+2*k)/(n-k-1) + wholeAIC
    splitAICc <- (2*ksplit+(2*(ksplit^2)+2*ksplit)/(n-ksplit-1))-2*(lik1+lik1)

    # BIC
    wholeBIC <- k*log(n)-2*(likWhole)
    splitBIC <- ksplit*log(n)-2*(lik1+lik2)

    op <- c(wholeAIC, splitAIC, wholeAICc, splitAICc, wholeBIC, splitBIC)
    names(op) <- c("wholeAIC", "splitAIC", "wholeAICc", "splitAICc", "wholeBIC", "splitBIC")

    #return(list(AIC=c(wholeAIC, splitAIC), AICc=c(wholeAICc, splitAICc), BIC=c(wholeBIC, splitBIC)))
    return(op)
}

# trying clade and stem data
trees <- read.tree("../sims/flcStem.nwk.tree")
load("../flcStemClade.RData")
stemCladeData <- data.frame()
for (i in 1:length(trees)){
    stemCladeData <- rbind(stemCladeData, getIC(tree=trees[[i]], clades=clades[[i]]))
}

trees <- read.tree("../sims/flcClade.nwk.tree")
load("../flcClades.RData")
cladeData <- data.frame()
for (i in 1:length(trees)){
    cladeData <- rbind(cladeData, getIC(tree=trees[[i]], clades=clades[[i]]))
}

trees <- read.tree("../sims/flcStem.nwk.tree")
load("../flcStem.RData")
stemData <- data.frame()
for (i in 1:length(trees)){
    stemData <- rbind(stemData, getIC(tree=trees[[i]], clades=clades[[i]]))
}


c(any(stemCladeData[,1] < stemCladeData[,2]), any(stemCladeData[,3] < stemCladeData[,4]), any(stemCladeData[,5] < stemCladeData[,6]))
c(any(cladeData[,1] < cladeData[,2]), any(cladeData[,3] < cladeData[,4]), any(cladeData[,5] < cladeData[,6]))
c(any(stemData[,1] < stemData[,2]), any(stemData[,3] < stemData[,4]), any(stemData[,5] < stemData[,6]))
# ... so all support the split model as expected

# plotting
ggplot() +
    geom_smooth(data=data, aes(x=date, y=height), se=F, method='lm', col="darkslate") +
    geom_smooth(data=data, aes(x=date, y=height, group=clade, col=clade), se=F, method='lm', lty=2) +
    geom_point(data=data, aes(x=date, y=height, col=clade), shape=4, alpha=0.5, size=3) +
    theme_minimal()
