####################
####################
####################
#Preconference "Estimating and interpreting psychological networks"
#
#Code for example
#
####################
####################
####################

##########
####Get started

###set your working directory
setwd("C:/Users/Jens Lange/OneDrive/Zusatz/Eingeladene Vortr?ge/2019/Oslo/Workshop")

#R version 3.6.1
###packages
library(qgraph)            #plotting networks; version 1.6.3
library(bootnet)           #estimating networks; version 1.2.3
library(psych)             #psych statistics; version 1.8.12
library(rio)               #import/export data; version 0.5.16
library(CliquePercolation) #clique percolation community detection; version 0.2.0 

###load data - ratings of Obama
data <- import("Data_mixed.csv")

###Items
###Scale from from 0 (not at all) to 6 (very much)
##Awe
#AS_1 = "I feel my jaw drop"
#AS_2 = "I gasp"
#AS_3 = "I feel that I am in the presence of something grand"
#AS_4 = "I sense things momentarily slow down"
#AS_5 = "I have goosebumps"
#AS_6 = "I feel challenged to mentally process what I am experiencing"
#AS_7 = "I feel my eyes widen"
#AS_8 = "I have the sense of being connected to everything"
#AS_9 = "I have chills"
#AS_10 = "I feel awed"
#AS_11 = "I feel that my sense of self is diminished"
#
##Kama Muta
#KA_1 = "I have moist eyes"
#KA_2 = "I have positive feelings"
#KA_3 = "I feel like telling someone how much I care about them"
#KA_4 = "I have difficulty speaking"
#KA_5 = "A warm feeling in the center of the chest"
#KA_6 = "I feel refreshed, energized, or exhilarated"
#KA_7 = "A lump in the throat"
#KA_8 = "I am moved"
#KA_9 = "I feel buoyant or light"
#KA_10 = "Some feeling in the center of the chest"
#KA_11 = "I feel/observed an incredible bond"
#KA_12 = "It is heartwarming"
#KA_13 = "I smile"
#KA_14 = "I shed tears"
#KA_15 = "I feel choked up"
#KA_16 = "I am touched"

##########
####Descriptive Statistics

describe(data)
#-> one missing value
#-> no obvious skewness or kurtosis

#automatic detection of reasonable kind of correlation
#here polychoric correlations
correlations <- cor_auto(data)
correlations

###plotting correlations with qgraph
cor_graph <- qgraph(correlations)

#+ using Fruchterman-Reingold
cor_graph <- qgraph(correlations, layout = "spring")

#+ change for color-blind people
cor_graph <- qgraph(correlations, layout = "spring", theme = "colorblind")

#+ change shape of nodes
cor_graph <- qgraph(correlations, layout = "spring", theme = "colorblind",
                    shape = "triangle")

#+ change edge width scale
cor_graph <- qgraph(correlations, layout = "spring", theme = "colorblind",
                    shape = "triangle", esize = 3)

#for list of features see options in qgraph help page
?qgraph

##########
####Network Estimation - Gaussian Graphical Model (Continuous Data)

###estimate regularized partial correlation network
###EBICglasso (gLASSO with EBIC model selection)
###correlations determined via cor_auto
###pairwise deletion of missing values
GGM_net <- estimateNetwork(data, default = "EBICglasso", corMethod = "cor_auto",
                           missing = "pairwise")

###plot network with qgraph
GGM_graph <- plot(GGM_net, layout = "spring", theme = "colorblind")

###declare awe and kama muta items separate groups
groups <- c(rep("Awe",11), rep("Kama Muta",16))

###plot network with groups; no legend
GGM_graph <- plot(GGM_net, layout = "spring", theme = "colorblind",
                  groups = groups, legend = FALSE)

#save layout of the graph
layout <- GGM_graph$layout

###get weights matrix
Wmat_GGM <- getWmat(GGM_graph)
Wmat_GGM

##########
###Node Centrality

#determine values
cent <- centrality(GGM_net)

#get values
cent$OutDegree
cent$Closeness

#plot values
centralityPlot(GGM_net, scale = "raw", include = c("Strength","Closeness"))

##########
###Network Stability and Difference Tests

###non-parametric bootstrap
#-> for edge stability and edge as well as centrality difference tests
#run bootstrap
# set.seed(4186)
# boot1 <- bootnet(GGM_net, statistics = c("edge","Strength","Closeness"),
#                  nboots = 1000, nCores = 2, type = "nonparametric")
# save(boot1, file = "boot_edges.RData")
load("boot_edges.RData")

###case-dropping bootstrap
#-> for centrality stability
#run bootstrap
# set.seed(4186)
# boot2 <- bootnet(GGM_net, statistics = c("Strength","Closeness"),
#                  nboots = 1000, nCores = 2, type = "case")
# save(boot2, file = "boot_centrality.RData")
load("boot_centrality.RData")

###stability
#plot edge CIs
plot(boot1, statistics = "edge", labels = TRUE, order = "sample")
#plot edge CIs + make labels visible by creating large pdf
pdf("edge_stability.pdf", height = 50)
plot(boot1, statistics = "edge", labels = TRUE, order = "sample")
dev.off()
#plot centrality stability
plot(boot2, statistics = c("Strength","Closeness"))
#check stability of centrality indices
corStability(boot2)
#-> strength can be interpreted with care

###difference tests
#plot edge weights difference test
plot(boot1, statistics = "edge", plot = "difference",
     onlyNonZero = TRUE, order = "sample")
#plot strength difference test
plot(boot1, statistics = "strength", plot = "difference",
     order = "sample")

##########
###Clique Percolation Community Detection

#run threshold function to determine k and I for a range of k and I values
threshold <- cpThreshold(W = GGM_graph, method = "weighted",
                         k.range = c(3,4,5),
                         I.range = c(seq(0.25, 0.01, by = -0.001)), 
                         threshold = c("largest.components.ratio","chi","entropy"))

#use permutation function to derive solutions more surprising than chance
# set.seed(4186)
# threshold_permute <- cpPermuteEntropy(W = GGM_graph, cpThreshold.object = threshold,
#                                       n = 100, interval = 0.95)
# save(threshold_permute, file = "threshold_permutation_mixed.Rdata")
load("threshold_permutation_mixed.Rdata") #loads previous permutation

#inspect results
threshold_permute$Confidence.Interval
threshold_permute$Extracted.Rows

#highest threshold for fewest isolated nodes obtained with k = 3 and I = .129
cp_ka_k3I.129 <- cpAlgorithm(W = GGM_graph, k = 3, 
                             method = "weighted", I = 0.129)

#inspect communities
cp_ka_k3I.129$list.of.communities.labels
#inspect shared nodes
cp_ka_k3I.129$shared.nodes.labels
#inspect isolated nodes
cp_ka_k3I.129$isolated.nodes.labels

##plot colored graph and community graph
#define which items belong to which construct, first construct awe, second kama muta. 
list.of.sets <- list(seq(from = 1, to = 11),
                     seq(from = 12, to = 27))
#colored graph
col_graph <- cpColoredGraph(GGM_graph, cp_ka_k3I.129$list.of.communities.numbers, layout=layout, 
                            list.of.sets = list.of.sets, set.palettes.size = 6,
                            theme = "colorblind")
#community graph
comm_graph <- cpCommunityGraph(cp_ka_k3I.129$list.of.communities.numbers,
                               node.size.method = "proportional",
                               max.node.size = 18, layout = "spring", theme = "colorblind")

##########
###Network Estimation - Ising Model (Binary Data)

#can be estimated via estimateNetwork
#variables automatically binarized at median
#estimate regularized logistic node-wise regression network
#define where to binarize variables
#eLASSO (LASSO with EBIC model selection)
#listwise deletion of missing values (pairwise not possible for regressions)
Ising_net <- estimateNetwork(data, default = "IsingFit",
                             missing = "listwise", rule = "OR")

#extract thresholds
thresholds <- Ising_net$intercepts

#extract weights matrix
Wmat_Ising <- getWmat(Ising_net)
Wmat_Ising

###compare GGM and Ising model
#correlating weights matrices
cor.test(Wmat_GGM[upper.tri(Wmat_GGM)], Wmat_Ising[upper.tri(Wmat_Ising)]) #.80

##########
###Further Comments

#centrality indices, stability, difference tests, and clique percolation are performed as before
