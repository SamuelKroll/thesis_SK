source('http://bioconductor.org/biocLite.R')
biocLite('phyloseq')
library(phyloseq)
library(GGally)
library(ggplot2)
rm(list=ls())

#############################################################
# Set up the working directory
#############################################################

wdir <- "/Volumes/projects/dep_psl/grp_kemen/working_directory_sam/Sequencing/Time_course/KeystoneTest/"

setwd(wdir)
outputFolder <- paste( wdir, "/Output_EdgeDependencies/", sep = "" )

# Read in the node table to get important info
NT <- "./Lists/Dec_Jan_1_Node_Table_dot.csv"
NodeTable <- read.table(file=NT,row.names=1, header = T, sep = ",", stringsAsFactors = F)

GroupsVector <- c()

#BacV5
DT <- "./Output_EdgeDependencies/BacV5_EdgeDepend_slope.txt"
BacV5_Dependencies <- read.table(file=DT, header=T, sep="\t")
DT2 <- "./Output_EdgeDependencies/BacV5_EdgeDepend_stmrsqadj.txt"
BacV5_Dependencies_stmrsqadj <- read.table(file=DT2, header=T, sep="\t")

new_colnames <- colnames(BacV5_Dependencies)
new_colnames <- c("Original", new_colnames)
new_colnames <- new_colnames[-length(new_colnames)]
new_colnames <- unlist(lapply(new_colnames, FUN=function(x){paste("BV5taxa",x[1],sep="_")}))
colnames(BacV5_Dependencies) <- new_colnames
colnames(BacV5_Dependencies_stmrsqadj) <- new_colnames
GroupsVector  <- c(GroupsVector, rep("BacV5", length(colnames(BacV5_Dependencies))))

#BacV3
DT <- "./Output_EdgeDependencies/BacV3_EdgeDepend_slope.txt"
BacV3_Dependencies <- read.table(file=DT, header=T, sep="\t")
DT2 <- "./Output_EdgeDependencies/BacV3_EdgeDepend_stmrsqadj.txt"
BacV3_Dependencies_stmrsqadj <- read.table(file=DT2, header=T, sep="\t")

new_colnames <- colnames(BacV3_Dependencies)
new_colnames <- c("Original", new_colnames)
new_colnames <- new_colnames[-length(new_colnames)]
new_colnames <- unlist(lapply(new_colnames, FUN=function(x){paste("BV3taxa",x[1],sep="_")}))
colnames(BacV3_Dependencies) <- new_colnames
colnames(BacV3_Dependencies_stmrsqadj) <- new_colnames
GroupsVector  <- c(GroupsVector, rep("BacV3", length(colnames(BacV3_Dependencies))))

#FITS2
DT <- "./Output_EdgeDependencies/FITS2_EdgeDepend_slope.txt"
FITS2_Dependencies <- read.table(file=DT, header=T, sep="\t")
DT2 <- "./Output_EdgeDependencies/FITS2_EdgeDepend_stmrsqadj.txt"
FITS2_Dependencies_stmrsqadj <- read.table(file=DT2, header=T, sep="\t")

new_colnames <- colnames(FITS2_Dependencies)
new_colnames <- c("Original", new_colnames)
new_colnames <- new_colnames[-length(new_colnames)]
new_colnames <- unlist(lapply(new_colnames, FUN=function(x){paste("FITS2taxa",x[1],sep="_")}))
colnames(FITS2_Dependencies) <- new_colnames
colnames(FITS2_Dependencies_stmrsqadj) <- new_colnames
GroupsVector  <- c(GroupsVector, rep("FITS2", length(colnames(FITS2_Dependencies))))

#Ftrad
DT <- "./Output_EdgeDependencies/Ftrad_EdgeDepend_slope.txt"
Ftrad_Dependencies <- read.table(file=DT, header=T, sep="\t")
DT2 <- "./Output_EdgeDependencies/Ftrad_EdgeDepend_stmrsqadj.txt"
Ftrad_Dependencies_stmrsqadj <- read.table(file=DT2, header=T, sep="\t")

new_colnames <- colnames(Ftrad_Dependencies)
new_colnames <- c("Original", new_colnames)
new_colnames <- new_colnames[-length(new_colnames)]
new_colnames <- unlist(lapply(new_colnames, FUN=function(x){paste("Ftradtaxa",x[1],sep="_")}))
colnames(Ftrad_Dependencies) <- new_colnames
colnames(Ftrad_Dependencies_stmrsqadj) <- new_colnames
GroupsVector  <- c(GroupsVector, rep("Ftrad", length(colnames(Ftrad_Dependencies))))

#OITS2
DT <- "./Output_EdgeDependencies/OITS2_EdgeDepend_slope.txt"
OITS2_Dependencies <- read.table(file=DT, header=T, sep="\t")
DT2 <- "./Output_EdgeDependencies/OITS2_EdgeDepend_stmrsqadj.txt"
OITS2_Dependencies_stmrsqadj <- read.table(file=DT2, header=T, sep="\t")

new_colnames <- colnames(OITS2_Dependencies)
new_colnames <- c("Original", new_colnames)
new_colnames <- new_colnames[-length(new_colnames)]
new_colnames <- unlist(lapply(new_colnames, FUN=function(x){paste("OITS2taxa",x[1],sep="_")}))
colnames(OITS2_Dependencies) <- new_colnames
colnames(OITS2_Dependencies_stmrsqadj) <- new_colnames
GroupsVector  <- c(GroupsVector, rep("OITS2", length(colnames(OITS2_Dependencies))))

# Otrad
# DT <- "./Output_EdgeDependencies/Otrad_EdgeDepend_slope.txt"
# Otrad_Dependencies <- read.table(file=DT, header=T, sep="\t")
# DT2 <- "./Output_EdgeDependencies/Otrad_EdgeDepend_stmrsqadj.txt"
# Otrad_Dependencies_stmrsqadj <- read.table(file=DT2, header=T, sep="\t")

# new_colnames <- colnames(Otrad_Dependencies)
# new_colnames <- c("Original", new_colnames)
# new_colnames <- new_colnames[-length(new_colnames)]
# new_colnames <- unlist(lapply(new_colnames, FUN=function(x){paste("Otrad",x[1],sep="_")}))
# colnames(Otrad_Dependencies) <- new_colnames
# colnames(Otrad_Dependencies_stmrsqadj) <- new_colnames
# GroupsVector  <- c(GroupsVector, rep("Otrad", length(colnames(Otrad_Dependencies))))

#Bind all Dependencies

All_Dependencies <- cbind(BacV5_Dependencies, BacV3_Dependencies, FITS2_Dependencies, Ftrad_Dependencies, OITS2_Dependencies)
All_Dependencies_stmrsqadj <- cbind(BacV5_Dependencies_stmrsqadj, BacV3_Dependencies_stmrsqadj, FITS2_Dependencies_stmrsqadj, Ftrad_Dependencies_stmrsqadj, OITS2_Dependencies_stmrsqadj)
All_Dependencies_stmrsqadj <- All_Dependencies_stmrsqadj[,-(grep("Original", colnames(All_Dependencies_stmrsqadj))[2: length(grep("Original", colnames(All_Dependencies_stmrsqadj)))])]
All_Dependencies_stmrsqadj_binary <- as.matrix((All_Dependencies_stmrsqadj<3) + 0)
LostEdges <- colSums(All_Dependencies_stmrsqadj_binary)

edge_table <- otu_table(All_Dependencies, taxa_are_rows=T)
bray_weighted_edges <- distance(edge_table, method ="euclidean")

k <- 2
pcoa <- cmdscale(bray_weighted_edges, k=k, eig=T)
points <- pcoa$points
eig <- pcoa$eig
points <- as.data.frame(points)
colnames(points) <- c("x", "y")

points$GroupsVector <- GroupsVector
points$GroupsVector[grep("Original", rownames(points))] <- "Original"
points$GroupsVector <- factor(points$GroupsVector)
# Get rid of the extra "Original" rows in points
points <- points[-(grep("Original", rownames(points))[2: length(grep("Original", rownames(points)))]),]
points$LostEdges <- LostEdges

#Put the nodetable in the same order as the points and add important vectors to points
y <- (NodeTable$name)
#NodeTable$SimpleName <- unlist(lapply(y, FUN=function(x){paste(x[2],x[3],sep="_")}))
NodeTable$WholeName <- paste(NodeTable$FromWhich, NodeTable$name, sep="_")
reorder_NodeTable <- NodeTable[match(row.names(points),NodeTable$WholeName),]

# Put the vectors that are desired into points
points$Degree <- reorder_NodeTable$Degree
points$BetweennessCentrality <- reorder_NodeTable$BetweennessCentrality
points$ClosenessCentrality <- reorder_NodeTable$ClosenessCentrality
points$EigenvectorCentrality <- reorder_NodeTable$Eigenvector
points$NeighborhoodConnectivity <- reorder_NodeTable$NeighborhoodConnectivity
points_plot <- points[-1,]

points_pairsplot <- points[-1,]
ggpairs(points_pairsplot, columns=3:8, ggplot2::aes(colour=GroupsVector, alpha=0.4))

labelColors=c("red","orange","yellow","brown","green","blue")
hc <- hclust(dist(points_plot[,4:9]), "ave")
pdf(file="test.pdf", width=15, height = 5)
plot(as.phylo(hc), type="unrooted", tip.color=labelColors[points_plot$GroupsVector])
dev.off()

pdf(file="KeystonenessByGroup.pdf", width=6, height=5)
ggplot(points, aes(x=x, y=y, color=GroupsVector))+ geom_point(alpha=.7, size=2) + theme(legend.position="none")
dev.off()
pdf(file="KeystonenessByGroup_Legend.pdf", width=6, height=5)
ggplot(points, aes(x=x, y=y, color=GroupsVector))+ geom_point(alpha=.7, size=2)
dev.off()

pdf(file="KeystonenessByDegree.pdf", width=6, height=5)
ggplot(points, aes(x=x, y=y, color=Degree))+ geom_point(alpha=.6, size=2) + theme(legend.position="none")
dev.off()
pdf(file="KeystonenessByDegree_Legend.pdf", width=6, height=5)
ggplot(points, aes(x=x, y=y, color=Degree))+ geom_point(alpha=.6, size=2)
dev.off()

pdf(file="KeystonenessByBC.pdf", width=6, height=5)
ggplot(points, aes(x=x, y=y, color=BetweennessCentrality))+ geom_point(alpha=.6, size=2) + theme(legend.position="none")
dev.off()
pdf(file="KeystonenessByBC_Legend.pdf", width=6, height=5)
ggplot(points, aes(x=x, y=y, color=BetweennessCentrality))+ geom_point(alpha=.6, size=2)
dev.off()

pdf(file="KeystonenessByCC.pdf", width=6, height=5)
ggplot(points, aes(x=x, y=y, color=ClosenessCentrality))+ geom_point(alpha=.6, size=2) + theme(legend.position="none")
dev.off()
pdf(file="KeystonenessByCC_Legend.pdf", width=6, height=5)
ggplot(points, aes(x=x, y=y, color=ClosenessCentrality))+ geom_point(alpha=.6, size=2)
dev.off()

pdf(file="KeystonenessByEC.pdf", width=6, height=5)
ggplot(points, aes(x=x, y=y, color=EigenvectorCentrality))+ geom_point(alpha=.6, size=2) + theme(legend.position="none")
dev.off()
pdf(file="KeystonenessByEC_Legend.pdf", width=6, height=5)
ggplot(points, aes(x=x, y=y, color=EigenvectorCentrality))+ geom_point(alpha=.6, size=2)
dev.off()

pdf(file="KeystonenessByNC.pdf", width=6, height=5)
ggplot(points, aes(x=x, y=y, color=NeighborhoodConnectivity))+ geom_point(alpha=.6, size=2) + theme(legend.position="none")
dev.off()
pdf(file="KeystonenessByNC_Legend.pdf", width=6, height=5)
ggplot(points, aes(x=x, y=y, color=NeighborhoodConnectivity))+ geom_point(alpha=.6, size=2)
dev.off()

eig[1]/sum(eig)
eig[2]/sum(eig)

All_Dependencies_stmrsqadj2 <- cbind(BacV5_Dependencies_stmrsqadj, BacV5_Dependencies_stmrsqadj, FITS2_Dependencies_stmrsqadj, FITS2_Dependencies_stmrsqadj, Otrad_Dependencies_stmrsqadj, Otrad_Dependencies_stmrsqadj)
All_Dependencies_stmrsqadj_binary2 <- as.matrix((All_Dependencies_stmrsqadj2>3) + 0)
edge_table_2 <- otu_table(All_Dependencies_stmrsqadj_binary2, taxa_are_rows=T)
bray_binary_edges_2 <- distance(edge_table_2, method ="euclidean")

k <- 2
pcoa_binary <- cmdscale(bray_binary_edges_2, k=k, eig=T)
points_binary <- pcoa_binary$points
eig_binary <- pcoa_binary$eig
points_binary <- as.data.frame(points_binary)
colnames(points_binary) <- c("x", "y")

points_binary$GroupsVector <- GroupsVector
points_binary$GroupsVector[grep("Original", rownames(points_binary))] <- "Original"
points_binary$GroupsVector <- factor(points_binary$GroupsVector)
# Get rid of the extra "Original" rows in points
points_binary <- points_binary[-(grep("Original", rownames(points_binary))[2: length(grep("Original", rownames(points_binary)))]),]

ggplot(points_binary, aes(x=x, y=y, color=GroupsVector))+ geom_point(alpha=.7, size=2) + theme(legend.position="none")

ggplot(points, aes(x=x, y=y, color=GroupsVector))+ geom_point(alpha=.8, size=5) 
ggplot(points_binary, aes(x=x, y=y, color="red"))+ geom_point(alpha=.8, size=5)

bray_weighted_edges_matrix <- as.matrix(bray_weighted_edges)