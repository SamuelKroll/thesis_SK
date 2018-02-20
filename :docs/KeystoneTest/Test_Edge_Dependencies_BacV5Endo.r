# This script will check the dependencies of edges of a predetermined network on a given set of abundances of microorganisms
library(MASS)
library(ppcor)
rm(list=ls())

#############################################################
# Set up the working directory
#############################################################

wdir <- "/Volumes/projects/dep_psl/grp_kemen/working_directory_sam/Sequencing/Time_course/New_analysis/Figures/Keystone/KeystoneTest/"

setwd(wdir)
outputFolder <- paste( wdir, "/Output_EdgeDependencies/", sep = "" )

#############################################################
# Load the , Edge Table, Test Table, OTU tables and mapfiles
#############################################################

ET <- "./Lists/Master_EdgeTable6_cleaned_header.txt"
EdgeTable <- read.table(ET, row.names=1, header=T, sep="\t")
NumEdges <- nrow(EdgeTable) #ncol

NT <- "./Lists/Node_table_year2.csv"
NodeTable <- read.table(NT, row.names=1, header=T, sep="\t")
NumNode <- nrow(NodeTable)

TN <- "./Output_Testtables/BacV3_TestTable_dec_jan_1.txt"
TestTable <- read.table(TN, row.names=1, header=T, sep="\t")
NumTests <- ncol(TestTable) - 1

BacV5 <- "./OTU_Tables/BV5_otu_table_mc2_s2n50_filtered_sub_8000.txt"
BacV3 <- "./OTU_Tables/BV3_otu_table_mc2_s2n50_filtered_sub_3000.txt"
FITS2 <- "./OTU_Tables/FITS2_otu_table_mc2_s2n50_filtered_sub_500.txt"
Ftrad <- "./OTU_Tables/Ftrad_otu_table_mc2_s2n50_filtered_sub_3000.txt"
OITS2 <- "./OTU_Tables/OITS2_otu_table_mc2_s2n50_filtered_sub_100.txt"
Otrad <- "./OTU_Tables/Otrad_otu_table_mc2_s2n50_filtered_sub_2000.txt"

BacV5_OTUtable <- read.table(BacV5,  header=T, row.names=1, sep="\t", blank.lines.skip = FALSE)
BacV3_OTUtable <- read.table(BacV3,  header=T, row.names=1, sep="\t", blank.lines.skip = FALSE)
FITS2_OTUtable <- read.table(FITS2,  header=T, row.names=1, sep="\t", blank.lines.skip = FALSE)
Ftrad_OTUtable <- read.table(Ftrad,  header=T, row.names=1, sep="\t", blank.lines.skip = FALSE)
OITS2_OTUtable <- read.table(OITS2,  header=T, row.names=1, sep="\t", blank.lines.skip = FALSE)
Otrad_OTUtable <- read.table(Otrad,  header=T, row.names=1, sep="\t", blank.lines.skip = FALSE)

BacV5 <- "./Mapfiles/mapfile_BV5.txt"
BacV3 <- "./Mapfiles/mapfile_BV3.txt"
FITS2 <- "./Mapfiles/mapfile_FITS2.txt"
Ftrad <- "./Mapfiles/mapfile_Ftrad.txt"
OITS2 <- "./Mapfiles/mapfile_OITS2.txt"
Otrad <- "./Mapfiles/mapfile_Otrad.txt"

BacV5_Mapfile <- read.table(BacV5, header=T, row.names=1, comment.char="", sep="\t")
BacV3_Mapfile <- read.table(BacV3, header=T, row.names=1, comment.char="", sep="\t")
FITS2_Mapfile <- read.table(FITS2, header=T, row.names=1, comment.char="", sep="\t")
Ftrad_Mapfile <- read.table(Ftrad, header=T, row.names=1, comment.char="", sep="\t")
OITS2_Mapfile <- read.table(OITS2, header=T, row.names=1, comment.char="", sep="\t")
Otrad_Mapfile <- read.table(Otrad, header=T, row.names=1, comment.char="", sep="\t")

#############################################################
# Sync the OTU tables and mapfiles to one another
#############################################################

BacV5_OTUtable <- BacV5_OTUtable[,-ncol(BacV5_OTUtable)]
BacV5_OTUtable <- t(BacV5_OTUtable)
reorder_BacV5_Mapfile <- BacV5_Mapfile[match(row.names(BacV5_OTUtable),row.names(BacV5_Mapfile)),]
nona_reorder_BacV5_Mapfile <- reorder_BacV5_Mapfile[complete.cases(reorder_BacV5_Mapfile),]
reorder_BacV5_OTUtable <- BacV5_OTUtable[match(row.names(nona_reorder_BacV5_Mapfile),row.names(BacV5_OTUtable)),]
nona_reorder_BacV5_OTUtable <- reorder_BacV5_OTUtable[complete.cases(reorder_BacV5_OTUtable),]

BacV3_OTUtable <- BacV3_OTUtable[,-ncol(BacV3_OTUtable)]
BacV3_OTUtable <- t(BacV3_OTUtable)
reorder_BacV3_Mapfile <- BacV3_Mapfile[match(row.names(BacV3_OTUtable),row.names(BacV3_Mapfile)),]
nona_reorder_BacV3_Mapfile <- reorder_BacV3_Mapfile[complete.cases(reorder_BacV3_Mapfile),]
reorder_BacV3_OTUtable <- BacV3_OTUtable[match(row.names(nona_reorder_BacV3_Mapfile),row.names(BacV3_OTUtable)),]
nona_reorder_BacV3_OTUtable <- reorder_BacV3_OTUtable[complete.cases(reorder_BacV3_OTUtable),]


FITS2_OTUtable <- FITS2_OTUtable[,-ncol(FITS2_OTUtable)]
FITS2_OTUtable <- t(FITS2_OTUtable)
reorder_FITS2_Mapfile <- FITS2_Mapfile[match(row.names(FITS2_OTUtable),row.names(FITS2_Mapfile)),]
nona_reorder_FITS2_Mapfile <- reorder_FITS2_Mapfile[complete.cases(reorder_FITS2_Mapfile),]
reorder_FITS2_OTUtable <- FITS2_OTUtable[match(row.names(nona_reorder_FITS2_Mapfile),row.names(FITS2_OTUtable)),]
nona_reorder_FITS2_OTUtable <- reorder_FITS2_OTUtable[complete.cases(reorder_FITS2_OTUtable),]

Ftrad_OTUtable <- Ftrad_OTUtable[,-ncol(Ftrad_OTUtable)]
Ftrad_OTUtable <- t(Ftrad_OTUtable)
reorder_Ftrad_Mapfile <- Ftrad_Mapfile[match(row.names(Ftrad_OTUtable),row.names(Ftrad_Mapfile)),]
nona_reorder_Ftrad_Mapfile <- reorder_Ftrad_Mapfile[complete.cases(reorder_Ftrad_Mapfile),]
reorder_Ftrad_OTUtable <- Ftrad_OTUtable[match(row.names(nona_reorder_Ftrad_Mapfile),row.names(Ftrad_OTUtable)),]
nona_reorder_Ftrad_OTUtable <- reorder_Ftrad_OTUtable[complete.cases(reorder_Ftrad_OTUtable),]

OITS2_OTUtable <- OITS2_OTUtable[,-ncol(OITS2_OTUtable)]
OITS2_OTUtable <- t(OITS2_OTUtable)
reorder_OITS2_Mapfile <- OITS2_Mapfile[match(row.names(OITS2_OTUtable),row.names(OITS2_Mapfile)),]
nona_reorder_OITS2_Mapfile <- reorder_OITS2_Mapfile[complete.cases(reorder_OITS2_Mapfile),]
reorder_OITS2_OTUtable <- OITS2_OTUtable[match(row.names(nona_reorder_OITS2_Mapfile),row.names(OITS2_OTUtable)),]
nona_reorder_OITS2_OTUtable <- reorder_OITS2_OTUtable[complete.cases(reorder_OITS2_OTUtable),]

Otrad_OTUtable <- Otrad_OTUtable[,-ncol(Otrad_OTUtable)]
Otrad_OTUtable <- t(Otrad_OTUtable)
reorder_Otrad_Mapfile <- Otrad_Mapfile[match(row.names(Otrad_OTUtable),row.names(Otrad_Mapfile)),]
nona_reorder_Otrad_Mapfile <- reorder_Otrad_Mapfile[complete.cases(reorder_Otrad_Mapfile),]
reorder_Otrad_OTUtable <- Otrad_OTUtable[match(row.names(nona_reorder_Otrad_Mapfile),row.names(Otrad_OTUtable)),]
nona_reorder_Otrad_OTUtable <- reorder_Otrad_OTUtable[complete.cases(reorder_Otrad_OTUtable),]

#############################################################
#### Loop through the Edge Table and do the business   ######
#############################################################
OutTable_stmrsqadj <- data.frame(matrix(ncol = ncol(TestTable), nrow = nrow(EdgeTable)))
rownames(OutTable_stmrsqadj) <- rownames(EdgeTable)
colnames(OutTable_stmrsqadj) <- colnames(TestTable)

OutTable_slope <- data.frame(matrix(ncol = ncol(TestTable), nrow = nrow(EdgeTable)))
rownames(OutTable_slope) <- rownames(EdgeTable)
colnames(OutTable_slope) <- colnames(TestTable)

OutTable_rsq <- data.frame(matrix(ncol = ncol(TestTable), nrow = nrow(EdgeTable)))
rownames(OutTable_rsq) <- rownames(EdgeTable)
colnames(OutTable_rsq) <- colnames(TestTable)

OutTable_adjrsq <- data.frame(matrix(ncol = ncol(TestTable), nrow = nrow(EdgeTable)))
rownames(OutTable_adjrsq) <- rownames(EdgeTable)
colnames(OutTable_adjrsq) <- colnames(TestTable)

OutTable_pval <- data.frame(matrix(ncol = ncol(TestTable), nrow = nrow(EdgeTable)))
rownames(OutTable_pval) <- rownames(EdgeTable)
colnames(OutTable_pval) <- colnames(TestTable)

for(i in 1:NumEdges){
	
	OutEdge_stmrsqadj <- c()
	OutEdge_slope <- c()
	OutEdge_rsq <- c()
	OutEdge_adjrsq <- c()
	OutEdge_pval <- c()
	
	# For each edge, recover the species names and groups
	Sp1_name <- as.character(EdgeTable$OrgID[i])
	#splitname <- strsplit(Sp1, "_")
	#Sp1_name <- paste(splitname[[1]][2],splitname[[1]][3], sep="_")
	Sp2_name <- as.character(EdgeTable$targetID[i])
	#splitname <- strsplit(Sp2, "_")
	#Sp2_name <- paste(splitname[[1]][2],splitname[[1]][3], sep="_")
	
	Obs_In <- as.character(EdgeTable$Obs_In[i])
	splitname <- strsplit(Obs_In, "_")
	Sp1_Group <- splitname[[1]][1]
	Sp2_Group <- splitname[[1]][2]
	
	# Now get the vectors that will be used to calculate the correlations
	if(Sp1_Group == "BV5taxa"){
		Sp1_Vector <- nona_reorder_BacV5_OTUtable[,which(colnames(nona_reorder_BacV5_OTUtable) == Sp1_name)]
		Sp1_Reps <- nona_reorder_BacV5_Mapfile$Replicate
	}
	if(Sp1_Group == "BV3taxa"){
		Sp1_Vector <- nona_reorder_BacV3_OTUtable[,which(colnames(nona_reorder_BacV3_OTUtable) == Sp1_name)]
		Sp1_Reps <- nona_reorder_BacV3_Mapfile$Replicate
	}
	if(Sp1_Group == "FITS2taxa"){
		Sp1_Vector <- nona_reorder_FITS2_OTUtable[,which(colnames(nona_reorder_FITS2_OTUtable) == Sp1_name)]
		Sp1_Reps <- nona_reorder_FITS2_Mapfile$Replicate
	}
	if(Sp1_Group == "Ftradtaxa"){
		Sp1_Vector <- nona_reorder_Ftrad_OTUtable[,which(colnames(nona_reorder_Ftrad_OTUtable) == Sp1_name)]
		Sp1_Reps <- nona_reorder_Ftrad_Mapfile$Replicate
	}
	if(Sp1_Group == "OITS2taxa"){
		Sp1_Vector <- nona_reorder_OITS2_OTUtable[,which(colnames(nona_reorder_OITS2_OTUtable) == Sp1_name)]
		Sp1_Reps <- nona_reorder_OITS2_Mapfile$Replicate
	}
	if(Sp1_Group == "Otradtaxa"){
		Sp1_Vector <- nona_reorder_Otrad_OTUtable[,which(colnames(nona_reorder_Otrad_OTUtable) == Sp1_name)]
		Sp1_Reps <- nona_reorder_Otrad_Mapfile$Replicate
	}
	if(Sp2_Group == "BV5taxa"){
		Sp2_Vector <- nona_reorder_BacV5_OTUtable[,which(colnames(nona_reorder_BacV5_OTUtable) == Sp2_name)]
		Sp2_Reps <- nona_reorder_BacV5_Mapfile$Replicate
	}
	if(Sp2_Group == "BV3taxa"){
		Sp2_Vector <- nona_reorder_BacV3_OTUtable[,which(colnames(nona_reorder_BacV3_OTUtable) == Sp2_name)]
		Sp2_Reps <- nona_reorder_BacV3_Mapfile$Replicate
	}
	if(Sp2_Group == "FITS2taxa"){
		Sp2_Vector <- nona_reorder_FITS2_OTUtable[,which(colnames(nona_reorder_FITS2_OTUtable) == Sp2_name)]
		Sp2_Reps <- nona_reorder_FITS2_Mapfile$Replicate
	}
	if(Sp2_Group == "Ftradtaxa"){
		Sp2_Vector <- nona_reorder_Ftrad_OTUtable[,which(colnames(nona_reorder_Ftrad_OTUtable) == Sp2_name)]
		Sp2_Reps <- nona_reorder_Ftrad_Mapfile$Replicate
	}
	if(Sp2_Group == "OITS2taxa"){
		Sp2_Vector <- nona_reorder_OITS2_OTUtable[,which(colnames(nona_reorder_OITS2_OTUtable) == Sp2_name)]
		Sp2_Reps <- nona_reorder_OITS2_Mapfile$Replicate
	}
  if(Sp2_Group == "Otradtaxa"){
    Sp2_Vector <- nona_reorder_Otrad_OTUtable[,which(colnames(nona_reorder_Otrad_OTUtable) == Sp2_name)]
    Sp2_Reps <- nona_reorder_Otrad_Mapfile$Replicate
  }
	
	# Get the two vectors aligned to one another by Replicates
	Sp1All <- data.frame(Sp1_Vector, Sp1_Reps)
	Sp2All <- data.frame(Sp2_Vector, Sp2_Reps)
	reorder_Sp1All <- Sp1All[match(Sp2All$Sp2_Reps, Sp1All$Sp1_Reps),]
	nona_reorder_Sp1All <- reorder_Sp1All[complete.cases(reorder_Sp1All),]
	reorder_Sp2All <- Sp2All[match(nona_reorder_Sp1All$Sp1_Reps, Sp2All$Sp2_Reps),]
	nona_reorder_Sp2All <- reorder_Sp2All[complete.cases(reorder_Sp2All),]
	
	# Get the sum to max ratio and convert the abundances to logscale
	smr1 <- sum(nona_reorder_Sp1All$Sp1_Vector)/max(nona_reorder_Sp1All$Sp1_Vector)
	smr2 <- sum(nona_reorder_Sp2All$Sp2_Vector)/max(nona_reorder_Sp2All$Sp2_Vector)
	smr_ave <- (smr1 + smr2) / 2
	nona_reorder_Sp1All$Sp1_Vector <- log10(nona_reorder_Sp1All$Sp1_Vector + 1)
	nona_reorder_Sp2All$Sp2_Vector <- log10(nona_reorder_Sp2All$Sp2_Vector + 1)
	
	# Get the statistical values for the basic correlation
	resultlm <- summary(lm(nona_reorder_Sp1All$Sp1_Vector ~ nona_reorder_Sp2All$Sp2_Vector))
    slope <- resultlm$coefficients[2,1]
    r_square <- resultlm$r.squared
    adj_r_square <- resultlm$adj.r.squared
    if(slope < 0){
    		r_square <- r_square * -1
         adj_r_square <- adj_r_square * -1
    }
    p_value <- pf(resultlm$fstatistic[1], resultlm$fstatistic[2], resultlm$fstatistic[3], lower.tail=F)
	stm_rsqadj <- abs(smr_ave * adj_r_square)
	
	# Take the first value that is wanted
	OutEdge_stmrsqadj <- c(stm_rsqadj)
	OutEdge_slope <- c(slope)
	OutEdge_rsq <- c(r_square)
	OutEdge_adjrsq <- c(adj_r_square)
	OutEdge_pval <- c(p_value)
	
	# Align the test organisms table to the main vector organisms and vice-versa
	reorder_TestTable <- TestTable[match(nona_reorder_Sp1All$Sp1_Reps, TestTable$Replicate),]
	nona_reorder_TestTable <- reorder_TestTable[complete.cases(reorder_TestTable),]
	nona_reorder2_Sp1All <- nona_reorder_Sp1All[match(nona_reorder_TestTable$Replicate, nona_reorder_Sp1All$Sp1_Reps),]
	nona2_reorder2_Sp1All <- nona_reorder2_Sp1All[complete.cases(nona_reorder2_Sp1All),]
	nona_reorder2_Sp2All <- nona_reorder_Sp2All[match(nona_reorder_TestTable$Replicate, nona_reorder_Sp2All$Sp2_Reps),]
	nona2_reorder2_Sp2All <- nona_reorder2_Sp2All[complete.cases(nona_reorder2_Sp2All),]
	
    
	# Loop through the Test organisms and calculate the partial correlations
	for(j in 1:NumTests){
		
        # Convert the conditioning vector to logscale
        log_condition <- log10(nona_reorder_TestTable[,j] + 1)
        
		resultlm <- pcor.test(nona2_reorder2_Sp1All$Sp1_Vector, nona2_reorder2_Sp2All$Sp2_Vector, log_condition, method = c("pearson"))
        r_square <- resultlm$estimate**2
        p_value <- resultlm$p.value
        adj_r_square <- 1 - ((1 - r_square) * ((length(nona2_reorder2_Sp1All$Sp1_Vector) - 1)/(length(nona2_reorder2_Sp1All$Sp1_Vector) - 2 - 1)))
        if(resultlm$estimate < 0){
        		r_square <- r_square * -1
             adj_r_square <- adj_r_square * -1
        	}
		stm_rsqadj <- abs(smr_ave * adj_r_square)
		
		# To get the partial slope, we need the lm model, which is also a partial model
		resultlm <- summary(lm(nona2_reorder2_Sp1All$Sp1_Vector ~ nona2_reorder2_Sp2All$Sp2_Vector + log_condition))
		slope <- resultlm $coefficients[2,1]
		
		OutEdge_stmrsqadj <- c(OutEdge_stmrsqadj, stm_rsqadj)
		OutEdge_slope <- c(OutEdge_slope, slope)
		OutEdge_rsq <- c(OutEdge_rsq, r_square)
		OutEdge_adjrsq <- c(OutEdge_adjrsq, adj_r_square)
		OutEdge_pval <- c(OutEdge_pval, p_value)
	}
	OutTable_stmrsqadj[i,] <- OutEdge_stmrsqadj
	OutTable_slope[i,] <- OutEdge_slope
	OutTable_rsq[i,] <- OutEdge_rsq
	OutTable_adjrsq[i,] <- OutEdge_adjrsq
	OutTable_pval[i,] <- OutEdge_pval
}


WT <- paste(outputFolder,"FITS2_EdgeDepend_stmrsqadj.txt",sep="")
write.table(OutTable_stmrsqadj, file=WT, sep="\t")

WT <- paste(outputFolder,"FITS2_EdgeDepend_slope.txt",sep="")
write.table(OutTable_slope, file=WT, sep="\t")

WT <- paste(outputFolder,"FITS2_EdgeDepend_rsq.txt",sep="")
write.table(OutTable_rsq, file=WT, sep="\t")

WT <- paste(outputFolder,"FITS2_EdgeDepend_adjrsq.txt",sep="")
write.table(OutTable_adjrsq, file=WT, sep="\t")

WT <- paste(outputFolder,"FITS2_EdgeDepend_pval.txt",sep="")
write.table(OutTable_pval, file=WT, sep="\t")

