rm(list = ls())
library(ggplot2)
# source("Helper_SimulateR_1_3.R") # The location of the helper function, please change it! 
# load("PPnormalizedcount.RData")
source("C:/Users/Cory/Desktop/UROP2019/Simulation/helper_SimulateR_1_3.R") 
load("C:/Users/Cory/Desktop/UROP2019/Simulation/PPnormalizedcount.RData")

# remove non-expression genes
select = apply(PPnormalizedcount,1,function(s){sum(s>0)>0})
sum(!select) # 12068 genes removed
PPnormalizedcount = PPnormalizedcount[select,]

# generate data
n.indv <- 300
n.gene <- nrow(PPnormalizedcount)
n.cluster <- 3
perc.cluster.gene = 0.05
n.cluster.gene = floor(n.gene*perc.cluster.gene)

vars = get_vars(PPnormalizedcount)
means = get_means(PPnormalizedcount)

result <- gene.simulate(real.data = PPnormalizedcount,
                        n.indv = n.indv, 
                        n.gene = n.gene, 
                        n.cluster = n.cluster,
                        cluster.effect.sd = 3 * sqrt(vars), # questinable
                        perc.cluster.gene = perc.cluster.gene, 
                        dist.cluster = rep(1,n.cluster), 
                        n.biop = 1,
                        s = 1, 
                        seed = 1234,
                        designMat.heatmap = F,
                        betaMat.heatmap = F)
k = result[[1]]
beta = result[[2]]
x = result[[3]]
size = result[[4]]
mu = result[[5]]


k.inv <- t(apply(k,1,function(x){ qnorm((rank(x)-(3/8))/(length(x)-2*(3/8)+1))}))
library(pheatmap)

png(filename = 'heatmap.png')
pheatmap(t(k.inv)[,1:n.cluster.gene], cluster_cols=T,cluster_rows=T, main = "individuals X genes")
dev.off()
