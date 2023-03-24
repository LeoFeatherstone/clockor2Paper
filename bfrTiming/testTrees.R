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

# on Lenovo Thinkpad with 11th Gen intel i7 processor
# TempEst v1.5.3
# chromium v111.0.5563.64
time <- data.frame(
    nTips =        c(100,      500,     1000,      5000,    10000),
    clockor2Time = c("0.47", "2.609", "7.819" , "235.979",   "1121.833"    ),
    tempestTime  = c("< 1",   "1.52",  "5.51",    "166.33", "820.0")
)

# on mac M1
# chrome v 110.5481.177
# tempest v1.5.3
time <- data.frame(
    nTips =        c(100,      500,     1000,      5000,    10000),
    clockor2Time = c("1.239",  "1.70", "4.43",    "95.26",  "422.11"),
    tempestTime  = c("<1s",    "2.40", "10.28",  "272.29", "1310.34")
)