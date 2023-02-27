# Clockor2:  Inferring global and local strict molecular clocks using root-to-tip regression

## Abstract

## Introduction
- Need to give distrinction of strict, global and local molecular clocks.
- RTT is still an attractive method for calibrating string molecular clocks due to its simplicity. However, as phylodynamics grows to consider larger datasets, and an increasing diversity of non-viral pathogens, strict molecular clocks with a global rate becomes an increasingly restrictive assumption. To this end we introduce clockor2, a root-to-tip regression tool that  can infer local as well as global string molecular clocks.
- Need to include specific examples of where local clocks emerge: VOCS, bacteria, host - mink
- Data also getting larger - think hill, du plessis, nordic? (will surely be many others)

- Mention usual way of fitting local clocks - beast and how many hours it can expect to take

## Methods
### Dependencies 
- Clockor2 makes extensive use of the `phylotree.js` JavaScript Library for internal representation and manipulation of trees.
- Plotting of trees uses `phylocanvas.js`, which was chosen due to its speed using WebGL plotting
- '`Plotly.js` Was used for plotting regression interface
### General model for global and local strict clocks
Given a tree and dates associated with each tip, a strict molecular clock can be expressed as a simple linear regression of tip heights against dates:
`height = rate * date + origin`,
where rate is the molecular evolutionary rate, typically in units of `subs/sites/time`.
We denote the set of all tips *`T`*, with subsets of *`T`* termed *collections* of tips *`Ci`*.  *`Ci`* are non-overlapping sets of tips whose union makes up *`T`*.  For each *`Ci`*, we can model height against dates to infer a molecular local to that group of tips. Note, we refer to *collections* of tips rather than *clades* for since collections need not represent an entire clade. For example, if a collection of tips is selected from a basal branch, inducing a clade, and another is selected from a descendent node, then the original collection no longer comprises an entire clade. Below, we give an example of the model for a tree with two local clocks.
- Insert a figure of a coloured tree for clarifications sake.
- Then describe factoring of likelihood to calcuate information criteria for the model

### Algorithm for local clock search
- Describe algorithm in pseudocode and include a gif of it working

### Best fitting reroot
-We implement functionality to select the best fitting root identically to Tempest.
- Briefly, this involves looping through each node to root at the preceeding branch.
- The optimal placement of the root along this branch is then chose using a `Golden section` algorithm, optimising the value of `R^2` 
- Note that we calculate the BFR for the the global clock case, against which all other models can be compared in a maximum parsimony manner.

## Results
- Include results about speed and accuracy of IC based model selection. Perhaps some RoC cuves?
- Include point about being able to handle trees of 10^5 tips (phylocanvas anyway). Time this once we get web-workers happening for the BFR.

## Discussion
- Clockor2 is a flexible framework that is part of the borader shift to web based applications with functionality that can eventually be converted to faster code, such as though bio-rust.
- Future Direction - a pairwise regression that eliminates pseudoreplication in RTT.
- Easier o extend than tempest

## Supplementary Material
- Derivations for AIC, AICc, and BIC for general strict clock models
