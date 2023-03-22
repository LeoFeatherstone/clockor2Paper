library(ape)

### Generating trees to test clockor2 bfr performance

## Trees for bfr timing, want num tips = 10^2,...,10^5
set.seed(1234)
nTips <- c(100, 500, 1000, 5000, 10000)
trees <- lapply(
    nTips, # initally wanted tree of 10^5 tips, but too large or R
    function(x) rtree(x)
)

# add noise to bl and tips names
for (i in seq_along(trees)) {
    trees[[i]]$edge.length <- trees[[i]]$edge.length*runif(0, 2, n=length(trees[[i]]$edge.length))
    dates <- as.Date(lubridate::date_decimal(diag(vcv.phylo(trees[[i]])) + 2000))
    trees[[i]]$tip.label <- paste0(trees[[i]]$tip.label, "_", dates)
}

for (i in seq_along(trees)) {
    write.tree(trees[[i]], file=paste0("bfrTimeTest", nTips[i], "tips.nwk"))
}

# record times
time <- data.frame(
    nTips =        c(100,      500,     1000,      5000,    10000),
    clockor2Time = c("0.47", "2.609", "7.819" , "235.979",   "1121.833"    ),
    tempestTime  = c("< 1",   "1.52",  "5.51",    "166.33", "820.0")
)