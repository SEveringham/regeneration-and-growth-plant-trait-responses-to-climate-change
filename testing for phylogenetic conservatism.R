### Testing phylogenetic signal

require(dplyr)
require(tidyr)
require(ggtree)
require(viridis)
require(phylobase)
require(phylosignal)
require(readr)
require(phytools)

# read in data

data <- read.csv("mydata.csv") %>%
        mutate(Species=as.character(Species)) %>%
        filter(Species!="Acacia hakeoidesT") %>%
        mutate(Species = ifelse(Species=="Acacia hakeoidesB", "Acacia hakeoides", Species))


#vascular plant tree
bigtree <- read.tree(file = "PhytoPhylo.tre") #reading tree
specieslist <- read.csv("specieslist.csv", header = T) #reading species list
nodes <- read_table2("nodes.csv") #reading nodes - nodes were downloaded from the S.phylomaker github page
source("https://raw.githubusercontent.com/jinyizju/S.PhyloMaker/master/R_codes%20for%20S.PhyloMaker")
resultphylo <- S.PhyloMaker(spList = specieslist, tree = bigtree, nodes = nodes, scenarios = "S3") #pruning tree
str(resultphylo) #looking at results of my tree

tiff("mytree.tiff", width = 6, height = 8, units = 'in', res = 300)
plot(resultphylo$Scenario.3,main="Scenarion Three") ## plotting tree
dev.off()

PrunedTree <- resultphylo$Scenario.3 ## this is my pruned tree, i picked scen 3 because it deals with missing spp best

Species_ <- gsub(" ", "_", data$Species)
data <- cbind(data, Species_)
data <- rename(data, "PrunedTree$tip.label"=Species_)
data <- right_join(data, as.data.frame(PrunedTree$tip.label), by = "PrunedTree$tip.label")
data <- rename(data, Species_="PrunedTree$tip.label")
species.vector <- as.vector(data$Species_)

#adds a column "phylo" to your data 
data[, "phylo"] <- species.vector
branchlengths <- compute.brlen(PrunedTree)

#List of taxa

phylo4d <- phylo4d(branchlengths, data$yi) ### in here and in the next two test rows below, 
                                            ##Change the "yi" variable to yi.x, yi.y, yi.x.x, 
                                            #etc... until you have analysed all the trait variables, 
                                              #note that the variables are:
                                            # yi= stem density
                                            # yi.x = seed mass
                                            # yi.y = seed shape
                                            # yi.x.x = seed viability
                                            # yi.y.y = seed dormancy
                                            # yi.x.x.x = germination success
                                            # yi.y.y.y = plant height
                                            # yi.x.x.x.x = biomass
                                            # yi.y.y.y.y = root:shoot ratio

teststat <- phyloSignal(phylo4d, rep =999)

teststat

## Two methods for getting signal (more in phylosig documentation)
#Method for signal 1 (K)
phylosig(branchlengths, (data$yi), method = "K", test = T)

#Method for signal 2 (lamba)
phylosig(branchlengths, (data$yi), method = "lambda", test = T)
  
#looking for local areas of strong signal
#i.e. groups with high correlation of future warming tolerance
#Run this all together to make a plot of signal
local.i <- lipaMoran(phylo4d, prox.phylo = "nNodes", as.p4d = TRUE)
points.col <- lipaMoran(phylo4d, prox.phylo = "nNodes")$p.value
save(points.col, file = "./signal")
points.col <- ifelse(points.col < 0.05, "red", "black") #red dots mean it's significant
dotplot.phylo4d(local.i, dot.col = points.col)
pdf("./Outputs/germphysignalplot.pdf", width = 15, height = 30)
