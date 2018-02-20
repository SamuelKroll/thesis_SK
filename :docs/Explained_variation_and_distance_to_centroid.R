#######Script cpoyrights belong to Samuel Kroll
###### for questions pls contact kroll@mpipz.mpg.de


library('grid')
library("ggplot2")
library("ape")
library('permute')
library("lattice")
library("vegan")
library("cluster")

setwd("/Volumes/projects/dep_psl/grp_kemen/working_directory_sam/Sequencing/Time_course/New_analysis/Figures/Variation/")
getwd()

######################   A script PcoA   ############################################
######## BV5 ####################################################################
pre_OTU_tab_BV5 <- read.csv("BV5_otu_table_mc2_s2n50_filtered_sub_8000_L6.txt", row.names =1, sep = "", header =T)
OTU_tab_BV5 <- t(pre_OTU_tab_BV5)
log_OTU_BV5 <- log10(OTU_tab_BV5 + 1)


map_BV5 <- read.csv("mapfile_BV5_final_trimmed_2.txt", sep = "\t", row.names =1, header = T)
ordered_Map_BV5 <- map_BV5[match(row.names(OTU_tab_BV5), row.names(map_BV5)),]
re_ordered_Map_BV5 <- ordered_Map_BV5[complete.cases(ordered_Map_BV5),]

#Calculate the PCoA via Bray-Curtis dissimilarity 
#BV5
dist_bray_BV5 <- vegdist(log_OTU_BV5, method= "bray", binary=FALSE)
dist_bray_binary_BV5 <- vegdist(log_OTU_BV5, method="bray", binary=TRUE)

cap_bray_BV5 <- capscale(dist_bray_BV5 ~ 1)
cap_bray_bin_BV5 <- capscale(dist_bray_binary_BV5 ~ 1)

################ ################ ################ ################ ################ 
################ ADONIS EXPLAINED VARIATION ################ ################ 
################ ################ ################ ################ ################ 

### including interaction#

adonis(log_OTU_BV5 ~ re_ordered_Map_BV5$Ecotype*re_ordered_Map_BV5$Experiment*re_ordered_Map_BV5$Month*re_ordered_Map_BV5$Temperatur*re_ordered_Map_BV5$Preci, permutations = 10000, method = 'bray') 
 
Variation <- paste(adonis(log_OTU_BV5 ~ re_ordered_Map_BV5$Ecotype*re_ordered_Map_BV5$Experiment*re_ordered_Map_BV5$Month*re_ordered_Map_BV5$Temperatur*re_ordered_Map_BV5$Preci, permutations = 10000, method = 'bray') )


B5 <- paste("BacV5_varitation.txt",sep="")
write.table(Variation,file=B5, sep = "")


##additive model#
adonis(log_OTU_BV5 ~ re_ordered_Map_BV5$Ecotype+re_ordered_Map_BV5$Year*re_ordered_Map_BV5$Month, permutations = 10000, method = 'bray') 

##including nested effects of month within years##
adonis (log_OTU_BV5 ~ re_ordered_Map_BV5$Ecotype+re_ordered_Map_BV5$Year+re_ordered_Map_BV5$Year:re_ordered_Map_BV5$Month, permutations = 1000, method = 'bray') 


######  Variation between treatments, Distance to centroid    ########

re_ordered_Map_BV5$Month <- factor(re_ordered_Map_BV5$Month, levels = c("Nov","Dec","Jan","Feb","Mar"))
re_ordered_Map_BV5$Timepoint <- factor(re_ordered_Map_BV5$Timepoint, levels = c("Nov-14","Nov-15","Nov-16","Dec-14","Dec-15","Dec-16","Jan-15","Jan-16","Jan-17","Feb-15","Feb-16","Feb-17","Mar-15","Mar-16","Mar-17"))
boxplot(betadisper(vegdist(log_OTU_BV5, method="bray"), re_ordered_Map_BV5$Month, type = c("median"), bias.adjust = TRUE),las=2,cex.axis = 0.7, col=(c("coral", "mediumseagreen","gold")))

BV5_vegdist <- betadisper(vegdist(log_OTU_BV5, method="bray"), re_ordered_Map_BV5$Month, type = c("median"), bias.adjust = TRUE)      
       
TukeyHSD(BV5_vegdist, ordered = TRUE,conf.level = 0.95)

outlier_values <- boxplot.stats(BV5_vegdist$Nov)$out
outlier_values

