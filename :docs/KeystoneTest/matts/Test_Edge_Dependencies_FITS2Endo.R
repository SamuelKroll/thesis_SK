# This script will check the dependencies of edges of a predetermined network on a given set of abundances of microorganisms
library(ppcor)
rm(list=ls())

#############################################################
# Set up the working directory
#############################################################

wdir <- "/Users/agler/LabStuff/Biodiversity/Project884_AtoH/MockCommunity/Improved_OPU_Pipeline/OPUs_Final/OPU_Network_B5F2O1/KeystoneTest"

setwd(wdir)
outputFolder <- paste( wdir, "/Output_EdgeDependencies/", sep = "" )

#############################################################
# Load the , Edge Table, Test Table, OPU tables and mapfiles
#############################################################

EN <- "./Lists/Master_EdgeTable_clean_TubSupp_stmrsqadj3.txt"
EdgeTable <- read.table(EN, row.names=1, header=T, sep="\t")
NumEdges <- nrow(EdgeTable)

TN <- "./Output_Testtables/FITS2Endo_TestTable.txt"
TestTable <- read.table(TN, row.names=1, header=T, sep="\t")
NumTests <- ncol(TestTable) - 1

BacV5Endo <- "./OPU_Tables/BacV5_OPU_table_s2n50_filtered_phylamin_Field_Endo_7000.txt"
BacV5Epi <- "./OPU_Tables/BacV5_OPU_table_s2n50_filtered_phylamin_Field_Epi_14000.txt"
FITS2Endo <- "./OPU_Tables/FITS2_OPU_table_s2n50_filtered_phylamin_Field_Endo_noOom_1000.txt"
FITS2Epi <- "./OPU_Tables/FITS2_OPU_table_s2n50_filtered_phylamin_Field_Epi_noOom_7000.txt"
OtradEndo <- "./OPU_Tables/Otrad_OPU_table_s2n50_filtered_phylamin_Field_Endo_noFungi_1000.txt"
OtradEpi <- "./OPU_Tables/Otrad_OPU_table_s2n50_filtered_phylamin_Field_Epi_noFungi_10000.txt"

BacV5Endo_OPUtable <- read.table( BacV5Endo,  header=T, row.names=1, sep="\t", blank.lines.skip = FALSE)
BacV5Epi_OPUtable <- read.table( BacV5Epi,  header=T, row.names=1, sep="\t", blank.lines.skip = FALSE)
FITS2Endo_OPUtable <- read.table( FITS2Endo,  header=T, row.names=1, sep="\t", blank.lines.skip = FALSE)
FITS2Epi_OPUtable <- read.table( FITS2Epi,  header=T, row.names=1, sep="\t", blank.lines.skip = FALSE)
OtradEndo_OPUtable <- read.table( OtradEndo,  header=T, row.names=1, sep="\t", blank.lines.skip = FALSE)
OtradEpi_OPUtable <- read.table( OtradEpi,  header=T, row.names=1, sep="\t", blank.lines.skip = FALSE)

BacV5Endo <- "./Mapfiles/Mapfile_new_BacV5_withdepth_withmaxtax_Field_Endo.txt"
BacV5Epi <- "./Mapfiles/Mapfile_new_BacV5_withdepth_withmaxtax_Field_Epi.txt"
FITS2Endo <- "./Mapfiles/Mapfile_new_FITS2_withdepth_withmaxtax_Field_Endo.txt"
FITS2Epi <- "./Mapfiles/Mapfile_new_FITS2_withdepth_withmaxtax_Field_Epi.txt"
OtradEndo <- "./Mapfiles/Mapfile_new_Otrad_withdepth_withmaxtax_Field_Endo.txt"
OtradEpi <- "./Mapfiles/Mapfile_new_Otrad_withdepth_withmaxtax_Field_Epi.txt"

BacV5Endo_Mapfile <- read.table(BacV5Endo, header=T, row.names=1, comment.char="", sep="\t")
BacV5Epi_Mapfile <- read.table(BacV5Epi, header=T, row.names=1, comment.char="", sep="\t")
FITS2Endo_Mapfile <- read.table(FITS2Endo, header=T, row.names=1, comment.char="", sep="\t")
FITS2Epi_Mapfile <- read.table(FITS2Epi, header=T, row.names=1, comment.char="", sep="\t")
OtradEndo_Mapfile <- read.table(OtradEndo, header=T, row.names=1, comment.char="", sep="\t")
OtradEpi_Mapfile <- read.table(OtradEpi, header=T, row.names=1, comment.char="", sep="\t")

#############################################################
# Sync the OPU tables and mapfiles to one another
#############################################################

BacV5Endo_OPUtable <- BacV5Endo_OPUtable[,-ncol(BacV5Endo_OPUtable)]
BacV5Endo_OPUtable <- t(BacV5Endo_OPUtable)
reorder_BacV5Endo_Mapfile <- BacV5Endo_Mapfile[match(row.names(BacV5Endo_OPUtable),row.names(BacV5Endo_Mapfile)),]
nona_reorder_BacV5Endo_Mapfile <- reorder_BacV5Endo_Mapfile[complete.cases(reorder_BacV5Endo_Mapfile),]
reorder_BacV5Endo_OPUtable <- BacV5Endo_OPUtable[match(row.names(nona_reorder_BacV5Endo_Mapfile),row.names(BacV5Endo_OPUtable)),]
nona_reorder_BacV5Endo_OPUtable <- reorder_BacV5Endo_OPUtable[complete.cases(reorder_BacV5Endo_OPUtable),]

BacV5Epi_OPUtable <- BacV5Epi_OPUtable[,-ncol(BacV5Epi_OPUtable)]
BacV5Epi_OPUtable <- t(BacV5Epi_OPUtable)
reorder_BacV5Epi_Mapfile <- BacV5Epi_Mapfile[match(row.names(BacV5Epi_OPUtable),row.names(BacV5Epi_Mapfile)),]
nona_reorder_BacV5Epi_Mapfile <- reorder_BacV5Epi_Mapfile[complete.cases(reorder_BacV5Epi_Mapfile),]
reorder_BacV5Epi_OPUtable <- BacV5Epi_OPUtable[match(row.names(nona_reorder_BacV5Epi_Mapfile),row.names(BacV5Epi_OPUtable)),]
nona_reorder_BacV5Epi_OPUtable <- reorder_BacV5Epi_OPUtable[complete.cases(reorder_BacV5Epi_OPUtable),]

FITS2Endo_OPUtable <- FITS2Endo_OPUtable[,-ncol(FITS2Endo_OPUtable)]
FITS2Endo_OPUtable <- t(FITS2Endo_OPUtable)
reorder_FITS2Endo_Mapfile <- FITS2Endo_Mapfile[match(row.names(FITS2Endo_OPUtable),row.names(FITS2Endo_Mapfile)),]
nona_reorder_FITS2Endo_Mapfile <- reorder_FITS2Endo_Mapfile[complete.cases(reorder_FITS2Endo_Mapfile),]
reorder_FITS2Endo_OPUtable <- FITS2Endo_OPUtable[match(row.names(nona_reorder_FITS2Endo_Mapfile),row.names(FITS2Endo_OPUtable)),]
nona_reorder_FITS2Endo_OPUtable <- reorder_FITS2Endo_OPUtable[complete.cases(reorder_FITS2Endo_OPUtable),]

FITS2Epi_OPUtable <- FITS2Epi_OPUtable[,-ncol(FITS2Epi_OPUtable)]
FITS2Epi_OPUtable <- t(FITS2Epi_OPUtable)
reorder_FITS2Epi_Mapfile <- FITS2Epi_Mapfile[match(row.names(FITS2Epi_OPUtable),row.names(FITS2Epi_Mapfile)),]
nona_reorder_FITS2Epi_Mapfile <- reorder_FITS2Epi_Mapfile[complete.cases(reorder_FITS2Epi_Mapfile),]
reorder_FITS2Epi_OPUtable <- FITS2Epi_OPUtable[match(row.names(nona_reorder_FITS2Epi_Mapfile),row.names(FITS2Epi_OPUtable)),]
nona_reorder_FITS2Epi_OPUtable <- reorder_FITS2Epi_OPUtable[complete.cases(reorder_FITS2Epi_OPUtable),]

OtradEndo_OPUtable <- OtradEndo_OPUtable[,-ncol(OtradEndo_OPUtable)]
OtradEndo_OPUtable <- t(OtradEndo_OPUtable)
reorder_OtradEndo_Mapfile <- OtradEndo_Mapfile[match(row.names(OtradEndo_OPUtable),row.names(OtradEndo_Mapfile)),]
nona_reorder_OtradEndo_Mapfile <- reorder_OtradEndo_Mapfile[complete.cases(reorder_OtradEndo_Mapfile),]
reorder_OtradEndo_OPUtable <- OtradEndo_OPUtable[match(row.names(nona_reorder_OtradEndo_Mapfile),row.names(OtradEndo_OPUtable)),]
nona_reorder_OtradEndo_OPUtable <- reorder_OtradEndo_OPUtable[complete.cases(reorder_OtradEndo_OPUtable),]

OtradEpi_OPUtable <- OtradEpi_OPUtable[,-ncol(OtradEpi_OPUtable)]
OtradEpi_OPUtable <- t(OtradEpi_OPUtable)
reorder_OtradEpi_Mapfile <- OtradEpi_Mapfile[match(row.names(OtradEpi_OPUtable),row.names(OtradEpi_Mapfile)),]
nona_reorder_OtradEpi_Mapfile <- reorder_OtradEpi_Mapfile[complete.cases(reorder_OtradEpi_Mapfile),]
reorder_OtradEpi_OPUtable <- OtradEpi_OPUtable[match(row.names(nona_reorder_OtradEpi_Mapfile),row.names(OtradEpi_OPUtable)),]
nona_reorder_OtradEpi_OPUtable <- reorder_OtradEpi_OPUtable[complete.cases(reorder_OtradEpi_OPUtable),]

#############################################################
# Loop through the Edge Table and do the business
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
	Sp1 <- as.character(EdgeTable$OrgID[i])
	splitname <- strsplit(Sp1, "_")
	Sp1_name <- paste(splitname[[1]][2],splitname[[1]][3], sep="_")
	Sp2 <- as.character(EdgeTable$targetID[i])
	splitname <- strsplit(Sp2, "_")
	Sp2_name <- paste(splitname[[1]][2],splitname[[1]][3], sep="_")
	
	Obs_In <- as.character(EdgeTable$Obs_In[i])
	splitname <- strsplit(Obs_In, "_")
	Sp1_Group <- splitname[[1]][1]
	Sp2_Group <- splitname[[1]][2]
	
	# Now get the vectors that will be used to calculate the correlations
	if(Sp1_Group == "BacV5Endo"){
		Sp1_Vector <- nona_reorder_BacV5Endo_OPUtable[,which(colnames(nona_reorder_BacV5Endo_OPUtable) == Sp1_name)]
		Sp1_Reps <- nona_reorder_BacV5Endo_Mapfile$Replicate
	}
	if(Sp1_Group == "BacV5Epi"){
		Sp1_Vector <- nona_reorder_BacV5Epi_OPUtable[,which(colnames(nona_reorder_BacV5Epi_OPUtable) == Sp1_name)]
		Sp1_Reps <- nona_reorder_BacV5Epi_Mapfile$Replicate
	}
	if(Sp1_Group == "FITS2Endo"){
		Sp1_Vector <- nona_reorder_FITS2Endo_OPUtable[,which(colnames(nona_reorder_FITS2Endo_OPUtable) == Sp1_name)]
		Sp1_Reps <- nona_reorder_FITS2Endo_Mapfile$Replicate
	}
	if(Sp1_Group == "FITS2Epi"){
		Sp1_Vector <- nona_reorder_FITS2Epi_OPUtable[,which(colnames(nona_reorder_FITS2Epi_OPUtable) == Sp1_name)]
		Sp1_Reps <- nona_reorder_FITS2Epi_Mapfile$Replicate
	}
	if(Sp1_Group == "OtradEndo"){
		Sp1_Vector <- nona_reorder_OtradEndo_OPUtable[,which(colnames(nona_reorder_OtradEndo_OPUtable) == Sp1_name)]
		Sp1_Reps <- nona_reorder_OtradEndo_Mapfile$Replicate
	}
	if(Sp1_Group == "OtradEpi"){
		Sp1_Vector <- nona_reorder_OtradEpi_OPUtable[,which(colnames(nona_reorder_OtradEpi_OPUtable) == Sp1_name)]
		Sp1_Reps <- nona_reorder_OtradEpi_Mapfile$Replicate
	}
	if(Sp2_Group == "BacV5Endo"){
		Sp2_Vector <- nona_reorder_BacV5Endo_OPUtable[,which(colnames(nona_reorder_BacV5Endo_OPUtable) == Sp2_name)]
		Sp2_Reps <- nona_reorder_BacV5Endo_Mapfile$Replicate
	}
	if(Sp2_Group == "BacV5Epi"){
		Sp2_Vector <- nona_reorder_BacV5Epi_OPUtable[,which(colnames(nona_reorder_BacV5Epi_OPUtable) == Sp2_name)]
		Sp2_Reps <- nona_reorder_BacV5Epi_Mapfile$Replicate
	}
	if(Sp2_Group == "FITS2Endo"){
		Sp2_Vector <- nona_reorder_FITS2Endo_OPUtable[,which(colnames(nona_reorder_FITS2Endo_OPUtable) == Sp2_name)]
		Sp2_Reps <- nona_reorder_FITS2Endo_Mapfile$Replicate
	}
	if(Sp2_Group == "FITS2Epi"){
		Sp2_Vector <- nona_reorder_FITS2Epi_OPUtable[,which(colnames(nona_reorder_FITS2Epi_OPUtable) == Sp2_name)]
		Sp2_Reps <- nona_reorder_FITS2Epi_Mapfile$Replicate
	}
	if(Sp2_Group == "OtradEndo"){
		Sp2_Vector <- nona_reorder_OtradEndo_OPUtable[,which(colnames(nona_reorder_OtradEndo_OPUtable) == Sp2_name)]
		Sp2_Reps <- nona_reorder_OtradEndo_Mapfile$Replicate
	}
	if(Sp2_Group == "OtradEpi"){
		Sp2_Vector <- nona_reorder_OtradEpi_OPUtable[,which(colnames(nona_reorder_OtradEpi_OPUtable) == Sp2_name)]
		Sp2_Reps <- nona_reorder_OtradEpi_Mapfile$Replicate
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

WT <- paste(outputFolder,"FITS2Endo_EdgeDepend_stmrsqadj.txt",sep="")
write.table(OutTable_stmrsqadj, file=WT, sep="\t")

WT <- paste(outputFolder,"FITS2Endo_EdgeDepend_slope.txt",sep="")
write.table(OutTable_slope, file=WT, sep="\t")

WT <- paste(outputFolder,"FITS2Endo_EdgeDepend_rsq.txt",sep="")
write.table(OutTable_rsq, file=WT, sep="\t")

WT <- paste(outputFolder,"FITS2Endo_EdgeDepend_adjrsq.txt",sep="")
write.table(OutTable_adjrsq, file=WT, sep="\t")

WT <- paste(outputFolder,"FITS2Endo_EdgeDepend_pval.txt",sep="")
write.table(OutTable_pval, file=WT, sep="\t")