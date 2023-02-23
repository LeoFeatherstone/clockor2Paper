## Testing regressions in R first. 

library(ape)
library(NELSI)

# read and process data
t <-  read.tree("sims/flcStemClade.nwk.tree")
heights <- lapply(t, function(x) allnode.times(x, tipsonly = T, reverse = T)) # double check that these are in fact heights
times <- lapply(t, function(x) as.numeric(gsub(x$tip.label, pattern='.+@', replacement='')))



# base regressions
baseMod <- list()
for (i in 1:length(t)){
    baseMod[[i]] <- summary(lm(heights[[i]]~times[[i]]))
    plot(x=times[[i]], y=heights[[i]])
    abline(baseMod[[i]], col='dodgerblue', lty=2)
    
}
# regressions with added group
load('flcClades.RData')


splitReg <- function() {}# later convert to a function newMod(trees, clades)


dat <- data.frame(times = times[[1]], tipLabs = t[[1]]$tip.label, heights = heights[[1]])
baseMod[[1]] <- list(base=summary(lm(dat$heights~dat$times)))
datPrime <- dat[which(dat$tipLabs %in% clades[[1]]),]
dat <- dat[-which(dat$tipLabs %in% clades[[1]]),]
newMod <- list(base = summary(lm(dat$heights~dat$times)), 
                split = summary(lm(datPrime$heights~datPrime$times)))

# now plot two models

# get information criteria


## up to trying to store this in a sensible objects



for (i in 1:length(t)){

    newMod[i] <- list()
    #for (j in 1:2) {
        newMod[i][[1]] <- summary(lm(heights[[i]]~times[[i]]))
        newMod[i][[2]] <- summary(lm(heights[[i]]~times[[i]]))
    #}
    newMod[[i]] <- summary(lm(heights[[i]]~times[[i]]))
    plot(x=times[[i]], y=heights[[i]])
    abline(newMod[[i]], col='dodgerblue', lty=2)
    
}
