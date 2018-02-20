#######Script cpoyrights belong to Samuel Kroll
###### for questions pls contact kroll@mpipz.mpg.de


#Keystone plot

library(ggplot2) 
library(RColorBrewer)
library(MASS)
setwd("/Volumes/biodata/dep_psl/grp_kemen/working_directory_sam/Sequencing/Time_course/New_analysis/Figures/Keystone//")

#rohdaten einlesen
data<-read.csv("keystoness_all_years_percentage_sample.txt", sep="\t", dec = ".")
str(data)

# calculate p-values for cut-off
# First Betweenness Centrality (Same for other parameters)

# Fit log-transformed data to normal distribution
fitlognorm <- fitdistr(log(data$Keystoness[which(data$Keystoness>0)]), "normal")

# Get the distribution parameters
ln_mean <- fitlognorm$estimate[1]
ln_sd <- fitlognorm$estimate[2]

# Minimum betweenness centrality for p = 0.05
bcmin <- exp(qnorm(0.95, mean=ln_mean, sd=ln_sd))
# Number of nodes above this level
bclen <- length(which(data$Keystoness > bcmin))
# List of nodes above this level
list <- which(data$Keystoness > bcmin)
bclist <- data$Taxonomy[list]

# Same for Closeness Centrality 

# Fit log-transformed data to normal distribution
fitlognorm_close <- fitdistr(log(data$Max_abund[which(data$Max_abund>0)]), "normal")

# Get the distribution parameters
ln_mean_close <- fitlognorm_close$estimate[1]
ln_sd_close <- fitlognorm_close$estimate[2]

# Minimum betweenness centrality for p = 0.05
bcmin_close <- exp(qnorm(0.95, mean=ln_mean_close, sd=ln_sd_close))
# Number of nodes above this level
bclen_close <- length(which(data$Max_abund > bcmin_close))
# List of nodes above this level
list_close <- which(data$Max_abund > bcmin_close)
bclist_close <- data$Taxonomy[list]


#definiere farben
#getPalette = colorRampPalette(brewer.pal(11, "RdYlBu")[c(10,2,5)])
#myColors <- getPalette(3)

expression <- ggplot(data, aes(x = Max_abund, y= Keystoness, fill= Organism, size =Degree, ymax=max(Keystoness)*1.1)) +
  theme(text = element_text(size=15)) 
z <- expression + geom_point(stat= "identity", position= position_dodge(),alpha= 0.7, shape= 21) + scale_fill_manual(values=(c("#56A5CC","#DF3A01","#62C152"))) +
  scale_size_area(max_size = 15)+
  geom_vline(xintercept = bcmin_close, linetype="dashed") +
  geom_hline(yintercept = bcmin, linetype= "dashed") +
  #geom_text(data=data,aes(x=ClosenessCentrality, y=BetweennessCentrality,label=Taxonomy),hjust = 1, size = 5,position=position_jitter()) +
  #stat_smooth(method=rlm, size=0.5, colour="black", level=0.95)+
  scale_x_continuous("Max_abund") +  
  scale_y_continuous("Keystoness")  

z + theme(panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.line.x = element_line(colour = "black"),
          axis.line.y = element_line(colour = "black"),
          axis.text = element_text(colour = "black"),
          axis.ticks = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.text.x = element_blank(),
          strip.background = element_blank(),
          legend.position ="top")


#speichern in working directory
ggsave("Keystones_vs_percentage_in_samples.pdf", useDingbats = FALSE)

