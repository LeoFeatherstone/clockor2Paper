library(ape)
library(scales)
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
R2 <- data.frame(
    nTips =        c(100,      500,     1000,      5000,    10000),
    tempestTime  = c("< 1",   "1.52",  "5.51",    "166.33", "820.0"),
    clockor2Time = c("0.097", "2.234", "8.933" , "313.048",   ""    )
)

RMS <- data.frame(
    nTips =        c(100,      500,     1000,      5000,    10000),
    tempestTime = c("",    "",   "",    "",  ""),
    clockor2Time  = c("0.162", "0.419", "2.152",  "46.935", "")
)

# on mac M1
# chrome v 113.0.5627.126
# tempest v1.5.3
R2 <- data.frame(
    nTips =        c(100,      500,     1000,      5000,    10000),
    tempestTime  = c("0.76",    "2.40", "10.28",    "272.29",  "1310.34"),
    clockor2Time = c("0.313",  "1.370", "3.476",    "78.013",  "306.821")
)

RMS <- data.frame(
    nTips =        c(100,      500,     1000,      5000,    10000),
    tempestTime = c("0.05",    "1.5",   "5.43",    "122",  "951.00"),
    clockor2Time  = c("0.129", "0.502", "1.514",  "24.370", "94.992")
)

tab <- bind_rows("R-Squared" = R2, "Residual Mean Squared" = RMS, .id = "Objective") %>% 
    pivot_longer(cols = ends_with("Time"), names_to = "Tool", values_to = "Time") %>%
    mutate(Time = as.numeric(Time)) %>%
    pivot_wider()

print(xtable::xtable(RMS, type = "latex"), file = "rmsTime.tex")
print(xtable::xtable(R2, type = "latex"), file = "r2Time.tex")


# try plot
bind_rows("R-Squared" = R2, "Residual Mean Squared" = RMS, .id = "Objective") %>% 
    pivot_longer(cols = ends_with("Time"), names_to = "Tool", values_to = "Time") %>%
    mutate(Time = as.numeric(Time)) %>%
    ggplot() +
    geom_line(aes(x = nTips, y = Time, col = Tool, lty = Objective), size = 2) +
    scale_y_continuous(trans = "log10") +
    scale_color_manual(values = c("red", "dodgerblue"))
