#######Script cpoyrights belong to Samuel Kroll
###### for questions pls contact kroll@mpipz.mpg.de

setwd("/Volumes/projects/dep_psl/grp_kemen/working_directory_sam/Sequencing/Time_course/New_analysis/Figures/Alpha_div/")
library("ggplot2")
library("ape")
library("agricolae")
library("multcomp")
library("mvtnorm")
library("survival")
library("TH.data")
library("MASS")

####        BV5      ###
Time_course <- read.table("BV5_div_combined.txt", header =TRUE, sep="\t")
class(Time_course)
summary(Time_course)
means <- colMeans(Time_course)
sds <- apply(Time_course, 2, sd)

Time_data <- as.data.frame(Time_course)

boxplot(Time_course, par(cex.axis = 0.80), las = 2, main = "Bacterial alpha diversity V5", ylab="Chao1 at 8000 reads",
        col=(c("#56A5CC")))
stripchart(vertical = TRUE, Time_course, 
           method = "jitter", add = TRUE, pch = 20, col =c("black"))
        
outlier_values <- boxplot.stats(Time_course$March)$out
outlier_values
