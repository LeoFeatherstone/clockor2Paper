################################################################################
################################################################################
######################### Simulation Study for clockSearch() ###################
################################################################################
################################################################################
library(ape)
library(TreeSim)
# devtools::install_github("sebastianduchene/NELSI")
library(NELSI)
library(tidyverse)
library(ggplot2)
library(latex2exp)
library(rjson)

setwd("./clockSearchSimStudy/")

## Get random clade within a specified size. I.e. the min possible to half of the num tips.
## Returns vector of tips comprising clade
getRandClade <- function(tr, minSize, maxSize) {
  # num nodes=2*ntips-1
  minNode <- length(tr$tip.label) + 1
  maxNode <- 2 * length(tr$tip.label) - 1

  # search for clade with between 3 and 50 descendants
  found <- FALSE
  while (!found) {
    node <- sample((minNode:maxNode), size = 1)
    nDesc <- length(geiger::tips(tr, node))
    if (nDesc >= minSize & nDesc <= maxSize) {
      clade <- geiger::tips(tr, node)
      found <- TRUE
    }
  }
  return(clade)
}

## Simulating 100 test trees. Using BD for heterochronous trees. Height added to 2000.
trees <- list()
class(trees) <- "multiPhylo"

set.seed(1234)
for (i in 1:100) {
  trees[i] <- sim.bdsky.stt(
    n = 250,
    lambdasky = 2.5,
    deathsky = 1,
    sampprobsky = 1,
    timesky = 0,
    timestop = 0
  )
}

# Add time to tip names
for (i in seq_along(trees)) {
  trees[[i]]$tip.label <- paste0(
    trees[[i]]$tip.label,
    "_",
    diag(vcv.phylo(trees[[i]]))
  )
}

## Getting clades between between 125 and 375 tips (25%-75% of tips)
clades <- lapply(
  trees,
  function(x) {
    getRandClade(
      x,
      minSize = 50,
      maxSize = 150
    )
  }
)

for (i in seq_along(clades)) {
  match <- which(trees[[i]]$tip.label %in% clades[[i]])
  clades[[i]] <- paste0(clades[[i]], "_Group1")
  trees[[i]]$tip.label[match] <- paste0(trees[[i]]$tip.label[match], "_Group1")
  trees[[i]]$tip.label[-match] <- paste0(trees[[i]]$tip.label[-match], "_Group2")
}

## Add noise to branch lengths
for (i in seq_along(trees)) {
  trees[[i]]$edge.length <- sapply(
    trees[[i]]$edge.length,
    function(x) {
      return(x + runif(min = -0.5 * x, max = 0.5 * x, n = 1))
    }
  )
}

## Generate simple FLC trees with global rate of 1e-3 and local rate of 5e-3
stemClade <- list()
stem <- list()
class(stemClade) <- "multiPhylo"
class(stem) <- "multiPhylo"

for (i in 1:length(trees)) {
  stemClade[[i]] <- NELSI::simulate.FLC(
    trees[[i]],
    list(
      clade.list = clades[i],
      stem.clade.indicator = list(c(T, T)),
      background.rate = 0.001,
      local.rates = 0.005
    )
  )$phylogram

  stem[[i]] <- NELSI::simulate.FLC(
    trees[[i]],
    list(
      clade.list = clades[i],
      stem.clade.indicator = list(c(T, F)),
      background.rate = 0.001,
      local.rates = 0.005
    )
  )$phylogram

  # rescale trees to subs/site
  trees[[i]]$edge.length <- 0.001 / trees[[i]]$edge.length
}

# wrapper function to use in testing clockSearch()
testClockSearch <- function(maxClocks, testTrees, baseTrees, trueClades) {
  df <- data.frame()
  names <- list()
  for (i in seq_along(baseTrees)) {
    cmd <- paste0("'", write.tree(testTrees[[i]]), "' ", 50, " ", maxClocks)
    system(
      # Enter node path
      paste("/home/leo/.nvm/versions/node/v20.9.0/bin/node ./clockSearchWrapper.js", cmd),
      intern = TRUE
    )
    op <- fromJSON(file = "tmp.json")
    # handle baseClock output
    if (class(op[[1]]) != "list") {
      op <- list(op)
    }

    data <- do.call(
      rbind,
      lapply(
        op,
        function(x) {
          return(
            c(
              i, maxClocks, length(op), x[["slope"]],
              max(
                length(which(x[["tip"]] %in% trueClades[[i]])) / length(trueClades[[i]]),
                1 - (length(which(x[["tip"]] %in% trueClades[[i]])) / length(trueClades[[i]]))
              ),
              length(x[["tip"]])
            )
          )
        }
      )
    )
    df <- rbind(df, data)
    print(paste("#", i, "/", length(baseTrees)))
  }
  colnames(df) <- c("id", "testMax", "nChosen", "slope", "maxPcMatch", "size")
  system("rm ./tmp.json")

  return(df)
}

globalClock <- lapply(
  c(2:4),
  function(x) {
    testClockSearch(
      x,
      trees,
      trees,
      clades
    )
  }
)

stemData <- lapply(
  c(2:4),
  function(x) {
    testClockSearch(
      x,
      stem,
      trees,
      clades
    )
  }
)

stemCladeData <- lapply(
  c(2:4),
  function(x) {
    testClockSearch(
      x,
      stemClade,
      trees,
      clades
    )
  }
)

data <- c(globalClock, stemData, stemCladeData)
names(data) <- c(rep("global", 3), rep("stem", 3), rep("stem+clade", 3))

data <- bind_rows(
  data,
  .id = "type"
)
save(data, trees, clades, stem, stemClade, file = "simStudyData.RData")

############################################################################
###################### Results and Figures below ###########################
############################################################################
load("./simStudyData.RData")

### Inferred clock group figs
data %>%
  select(
    type, testMax, nChosen, id
  ) %>%
  distinct() %>%
  group_by(
    type,
    testMax,
    nChosen
  ) %>%
  summarise(
    n = n()
  ) %>%
  ungroup() %>%
  complete(type, testMax, nChosen = 1:4, fill=list(n=0)) %>%
  ggplot(aes(testMax, nChosen, fill = n)) +
  geom_tile() +
  geom_text(aes(label = n)) +
  scale_fill_gradient(
    low = "lightgrey",
    high = "dodgerblue",
    name = "Identified Clocks"
  ) +
  xlab("Maximum Clocks in Local Clock Search") +
  ylab("Number of Clocks Inferred") +
  facet_wrap(~type) +
  coord_fixed() +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank()
  )
ggsave(
  file = "../inferredClocks.pdf",
  width = 6,
  height = 3,
  dpi = 300,
  units = "in"
)

### Rate error fig
scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}

data %>%
  select(
    type, testMax, slope, id, nChosen
  ) %>%
  group_by(
    type, testMax, slope, id, nChosen
  ) %>%
  mutate(
    error = case_when(
      type == "global" ~ abs(0.001 - slope),
      .default = min(
      abs(0.001 - slope), abs(0.005 - slope)
      )
    )
  ) %>%
  ggplot(aes(x = log10(error))) +
  geom_histogram(fill = alpha("dodgerblue", 0.6)) +
  #scale_x_continuous(label = scientific_10) +
  xlab("Absolute error in evolutionary rate (log10)") +
  ylab("Count") +
  facet_grid(type ~ nChosen, scales = "free") +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x =element_text(angle = 60, hjust = 1)
  )

ggsave(
  file = "../clockError.pdf",
  width = 6,
  height = 6,
  dpi = 300,
  units = "in"
)

### Clade identity fig
data %>%
  select(
    type, testMax, maxPcMatch, id, nChosen
  ) %>%
  ggplot(aes(x = 100 * maxPcMatch)) +
  geom_histogram(fill = alpha("dodgerblue", 0.6)) +
  xlab("Tue Clade Match (%)") +
  ylab("Count") +
  facet_grid(type ~ nChosen, scales = "free") +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x =element_text(angle = 60, hjust = 1)
  )

ggsave(
  file = "../cladeMatch.pdf",
  width = 6,
  height = 6,
  dpi = 300,
  units = "in"
)

#### example local clock Fig
t1 <- stem[[40]]
# t1$edge.length <- sapply(diag(vcv.phylo(t1)), function(x) x * runif(0.9, 1.1, n=1))
stemHeights <- sapply(diag(vcv.phylo(stem[[40]])), function(x) x * runif(0.9, 1.1, n = 1))
dateStem <- as.Date(lubridate::date_decimal(diag(vcv.phylo(trees[[40]])) + 2000))
namesStem <- trees[[40]]$tip.label
colsStem <- case_when(
  grepl(namesStem, pattern = "Group1") ~ "2",
  grepl(namesStem, pattern = "Group2") ~ "1"
)

df1 <- as_tibble(data.frame(dateStem, stemHeights, colsStem))
rtt1 <- ggplot(df1, aes(x = dateStem, y = stemHeights, fill = colsStem)) +
  geom_point(shape = 21) +
  scale_fill_brewer(palette = "Set1", direction = -1, "Local Clock") +
  annotate("text",
    x = as.Date("2003-01-01"),
    y = 0.001,
    label = "italic(r[1]) == 10^-3~subs/site/yr",
    parse = TRUE
  ) +
  annotate("text",
    x = as.Date("2002-01-01"),
    y = 0.01,
    label = "italic(r[2]) == 10^-3~subs/site/yr",
    parse = TRUE
  ) +
  xlab("") +
  ylab("Height (subs/site)") +
  theme_minimal() +
  theme(
    legend.position = "bottom"
  )

leg <- cowplot::get_legend(rtt1)
rtt1 <- rtt1 +
  geom_smooth(se = FALSE, method = "lm", col = "black", lwd = 0.5) +
  theme(legend.position = "none")

metadata1 <- data.frame("taxa" = namesStem, "tipCol" = colsStem)
tree1 <- ggtree(stem[[40]])
tree1 <- tree1 %<+%
  metadata1 +
  geom_tippoint(aes(fill = tipCol), shape = 21) +
  scale_fill_brewer(palette = "Set1", direction = -1) +
  theme(legend.position = "none")

t2 <- stemClade[[41]]
dateStemClade <- as.Date(lubridate::date_decimal(diag(vcv.phylo(trees[[41]])) + 2000))

stemCladeHeights <- sapply(diag(vcv.phylo(stemClade[[41]])), function(x) x * runif(0.9, 1.1, n = 1))

namesStemClade <- trees[[41]]$tip.label
colsStemClade <- case_when(
  grepl(namesStemClade, pattern = "Group1") ~ "2",
  grepl(namesStemClade, pattern = "Group2") ~ "1"
)

df2 <- as_tibble(data.frame(dateStemClade, stemCladeHeights, colsStemClade))
rtt2 <- ggplot(df2, aes(x = dateStemClade, y = stemCladeHeights, fill = colsStemClade)) +
  geom_point(shape = 21) +
  geom_smooth(se = FALSE, method = "lm", col = "black", lwd = 0.5) +
  scale_fill_brewer(palette = "Set1", direction = -1, "Local Clock") +
  annotate("text",
    x = as.Date("2002-07-01"),
    y = 0.00025,
    label = "italic(r[1]) == 10^-3~subs/site/yr",
    parse = TRUE
  ) +
  annotate("text",
    x = as.Date("2001-06-01"),
    y = 0.0125,
    label = "italic(r[2]) == 5%*%10^-3~subs/site/yr",
    parse = TRUE
  ) +
  xlab("Date") +
  ylab("Height (subs/site)") +
  theme_minimal() +
  theme(legend.position = "botnonetom")

metadata2 <- data.frame("taxa" = namesStemClade, "tipCol" = colsStemClade)
tree2 <- ggtree(stemClade[[41]])
tree2 <- tree2 %<+%
  metadata2 +
  geom_tippoint(aes(fill = tipCol), shape = 21) +
  scale_fill_brewer(palette = "Set1", direction = -1) +
  theme(legend.position = "none")

p1 <- cowplot::plot_grid(rtt1, tree1, rtt2, tree2, labels = "AUTO")
p2 <- cowplot::plot_grid(p1, leg, ncol = 1, rel_heights = c(6, 1))

pdf(file = "egRTT.pdf", width = 6, height = 6, useDingbats = F)
p2
dev.off()
