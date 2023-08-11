#testing_musse#
library(here)
library(ape)
library (diversitree)
library(phytools)
library(geiger)
tree <- read.tree('data/raw/tree.txt')
X <- read.csv("data/raw/SSE.csv", row.names = 1)
plot(tree)

plotTree(tree)
ltt(tree, show.tree=TRUE, lwd=3)

Diet <- setNames(X$Diet, rownames(X))
samplingf<-c(0.58, 0.78, 0.85)
p <- starting.point.musse(tree, k=3)
obj <- birthdeath(tree)
bd(obj)
musse <- make.musse(tree, states=Diet, k=3)
result <- find.mle(musse, x.init=p[argnames(musse)])  
result$par
