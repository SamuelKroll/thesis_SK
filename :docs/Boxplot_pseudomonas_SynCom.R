#######Script cpoyrights belong to Samuel Kroll
###### for questions pls contact kroll@mpipz.mpg.de


setwd("/Volumes/projects/dep_psl/grp_kemen/working_directory_sam/Sequencing/Time_course/New_analysis/Figures/SynCom_core/")
library("ggplot2")
library("ape")
library("agricolae")
library("multcomp")
library("mvtnorm")
library("survival")
library("TH.data")
library("MASS")

####        BV5      ###
Time_course <- read.table("Pseudomonas_syncom.txt", header =TRUE, sep="\t")
class(Time_course)
summary(Time_course)
means <- colMeans(Time_course)
sds <- apply(Time_course, 2, sd)
hist(Time_course$total.SynCom)


Time_data <- as.data.frame(Time_course)

log_Timecourse <- log10(Time_data + 1)

write.table(log_Timecourse, file ="log_Pst.txt", sep = ",")


boxplot(log_Timecourse, par(cex.axis = 0.80), las = 2, ylab="Pst log10(CFU/cm2)",
        col=(c("#9EDCEF","#E7EA88","#E7EA88","#ABDBA7","#ABDBA7","#ABDBA7","#EABB88","#EABB88","#EABB88","#BACFEE"))) 
stripchart(vertical = TRUE, log_Timecourse, 
           method = "jitter", add = TRUE, pch = 20, col =c("black"))
        

outlier_values <- boxplot.stats(Time_course$SynCom.w.o.Methylobacterium)$out
outlier_values

ggplot(log_Timecourse) + 
  geom_boxplot(outlier.shape=NA) + #avoid plotting outliers twice
  geom_jitter(position=position_jitter(width=.1, height=0))


#Statistics with Tukey HSD

expt1 <- read.table("stats_pseudomonas_syncom_adjust.txt",header=T)
amod <- aov(Y~X,data=expt1)
summary(amod)
library("multcomp") 
tmod <- glht(amod,linfct=mcp(X="Tukey"))
summary(tmod)
TukeyHSD(amod, conf.level=0.95)

#Data input for Tukey HSD
# Y         	X
# 0.28551035	A
# 0.338524035	A
# 0.088313218	A
# 0.205930807	A
# 0.363240102	A
# 0.52173913	B
# 0.763358779	B
# 0.32546786	B
# 0.425305688	B
# 0.378071834	B
# 0.989119683	C
# 1.192718142	C
# 0.788288288	C
# 0.549176236	C
# 0.544588155	C
# 1.26705653	D
# 1.625320787	D
# 1.266108976	D
# 1.154187629	D
# 1.268498943	D
# 1.069518717	D
