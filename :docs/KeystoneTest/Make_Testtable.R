# Generating a table with data to test partial correlations when controlling for node profiles. The script will take in an OTU table and 
# the corresponding mapfile as well as a list of nodes. It will then get the data from the OTU table for the nodes of interest. This can
# be repeated for any OTU tables as necessary.

rm(list=ls())

# #############################################################
# # USER INPUT: Provide some information about what sort of table we are working on
# #############################################################
# # Amplicon region ("BacV3", "BacV5", "Ftrad", "FITS2", "Otrad", "OITS2", etc.)
# Region <- "BacV5"
# FromWhich <- paste(Region,, sep="")

#############################################################
# Set up the working directory
#############################################################

wdir <- "/Volumes/projects/dep_psl/grp_kemen/working_directory_sam/Sequencing/Time_course/KeystoneTest"

setwd(wdir)
outputFolder <- paste( wdir, "/Output_Testtables/", sep = "" )

#############################################################


#############################################################
# Load the OTU table, mapfile and Node Table
#############################################################

#MN <- "./Mapfiles/Mapfile.txt"
#Mapfile <- read.table(MN, header=T, row.names=1, comment.char="", sep="\t")

NN <- "./Lists/Dec_Jan_1_Node_Table.csv"
NodeTable <- read.table(NN, row.names=1, header=T, sep=",")
NumNodes <- nrow(NodeTable)

BacV5 <- "./OTU_Tables/BV5_OTUtable_mc2_taxa_s2n50_filtered_sub_8000_L6_dec_jan_1.txt"
BacV3 <- "./OTU_Tables/BV3_OTUtable_mc2_taxa_s2n50_filtered_sub_8000_L6_dec_jan_1.txt"
Ftrad <- "./OTU_Tables/Ftrad_OTUtable_mc2_taxa_s2n50_filtered_sub_3000_L6_dec_jan_1.txt"
FITS2 <- "./OTU_Tables/FITS2_OTUtable_mc2_taxa_s2n50_filtered_sub_3000_L6_dec_jan_1.txt"
Otrad <- "./OTU_Tables/Otrad_OTUtable_mc2_taxa_s2n50_filtered_sub_2000_L6_dec_jan_1.txt"
OITS2 <- "./OTU_Tables/OITS2_OTUtable_mc2_taxa_s2n50_filtered_sub_2000_L6_dec_jan_1.txt"


BacV5_OTUtable <- read.table(BacV5,  header=T, row.names=1, sep="\t", blank.lines.skip = FALSE)
BacV3_OTUtable <- read.table(BacV3,  header=T, row.names=1, sep="\t", blank.lines.skip = FALSE)
Ftrad_OTUtable <- read.table(Ftrad,  header=T, row.names=1, sep="\t", blank.lines.skip = FALSE)
FITS2_OTUtable <- read.table(FITS2,  header=T, row.names=1, sep="\t", blank.lines.skip = FALSE)
Otrad_OTUtable <- read.table(Otrad,  header=T, row.names=1, sep="\t", blank.lines.skip = FALSE)
OITS2_OTUtable <- read.table(OITS2,  header=T, row.names=1, sep="\t", blank.lines.skip = FALSE)

BacV5 <- "./Mapfiles/mapfile_BV5_trimmed_dec_jan_1.txt"
BacV3 <- "./Mapfiles/mapfile_BV3_trimmed_dec_jan_1.txt"
Ftrad <- "./Mapfiles/mapfile_Ftrad_trimmed_dec_jan_1.txt"
FITS2 <- "./Mapfiles/mapfile_FITS2_trimmed_dec_jan_1.txt"
Otrad <- "./Mapfiles/mapfile_Otrad_trimmed_dec_jan_1.txt"
OITS2 <- "./Mapfiles/mapfile_OITS2_trimmed_dec_jan_1.txt"

BacV5_Mapfile <- read.table(BacV5, header=T, row.names=1, comment.char="", sep="\t")
BacV3_Mapfile <- read.table(BacV3, header=T, row.names=1, comment.char="", sep="\t")
Ftrad_Mapfile <- read.table(Ftrad, header=T, row.names=1, comment.char="", sep="\t")
FITS2_Mapfile <- read.table(FITS2, header=T, row.names=1, comment.char="", sep="\t")
Otrad_Mapfile <- read.table(Otrad, header=T, row.names=1, comment.char="", sep="\t")
OITS2_Mapfile <- read.table(OITS2, header=T, row.names=1, comment.char="", sep="\t")

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

Ftrad_OTUtable <- Ftrad_OTUtable[,-ncol(Ftrad_OTUtable)]
Ftrad_OTUtable <- t(Ftrad_OTUtable)
reorder_Ftrad_Mapfile <- Ftrad_Mapfile[match(row.names(Ftrad_OTUtable),row.names(Ftrad_Mapfile)),]
nona_reorder_Ftrad_Mapfile <- reorder_Ftrad_Mapfile[complete.cases(reorder_Ftrad_Mapfile),]
reorder_Ftrad_OTUtable <- Ftrad_OTUtable[match(row.names(nona_reorder_Ftrad_Mapfile),row.names(Ftrad_OTUtable)),]
nona_reorder_Ftrad_OTUtable <- reorder_Ftrad_OTUtable[complete.cases(reorder_Ftrad_OTUtable),]

FITS2_OTUtable <- FITS2_OTUtable[,-ncol(FITS2_OTUtable)]
FITS2_OTUtable <- t(FITS2_OTUtable)
reorder_FITS2_Mapfile <- FITS2_Mapfile[match(row.names(FITS2_OTUtable),row.names(FITS2_Mapfile)),]
nona_reorder_FITS2_Mapfile <- reorder_FITS2_Mapfile[complete.cases(reorder_FITS2_Mapfile),]
reorder_FITS2_OTUtable <- FITS2_OTUtable[match(row.names(nona_reorder_FITS2_Mapfile),row.names(FITS2_OTUtable)),]
nona_reorder_FITS2_OTUtable <- reorder_FITS2_OTUtable[complete.cases(reorder_FITS2_OTUtable),]

Otrad_OTUtable <- Otrad_OTUtable[,-ncol(Otrad_OTUtable)]
Otrad_OTUtable <- t(Otrad_OTUtable)
reorder_Otrad_Mapfile <- Otrad_Mapfile[match(row.names(Otrad_OTUtable),row.names(Otrad_Mapfile)),]
nona_reorder_Otrad_Mapfile <- reorder_Otrad_Mapfile[complete.cases(reorder_Otrad_Mapfile),]
reorder_Otrad_OTUtable <- Otrad_OTUtable[match(row.names(nona_reorder_Otrad_Mapfile),row.names(Otrad_OTUtable)),]
nona_reorder_Otrad_OTUtable <- reorder_Otrad_OTUtable[complete.cases(reorder_Otrad_OTUtable),]

OITS2_OTUtable <- OITS2_OTUtable[,-ncol(OITS2_OTUtable)]
OITS2_OTUtable <- t(OITS2_OTUtable)
reorder_OITS2_Mapfile <- OITS2_Mapfile[match(row.names(OITS2_OTUtable),row.names(OITS2_Mapfile)),]
nona_reorder_OITS2_Mapfile <- reorder_OITS2_Mapfile[complete.cases(reorder_OITS2_Mapfile),]
reorder_OITS2_OTUtable <- OITS2_OTUtable[match(row.names(nona_reorder_OITS2_Mapfile),row.names(OITS2_OTUtable)),]
nona_reorder_OITS2_OTUtable <- reorder_OITS2_OTUtable[complete.cases(reorder_OITS2_OTUtable),]


#############################################################
# Collect nodes
#############################################################
#Loop through the node table and make the output table
BacV5_OutputTable <- data.frame(nona_reorder_BacV5_OTUtable)
BacV3_OutputTable <- data.frame(row.names=rownames(nona_reorder_BacV3_OTUtable))
Ftrad_OutputTable <- data.frame(row.names=rownames(nona_reorder_Ftrad_OTUtable))
FITS2_OutputTable <- data.frame(row.names=rownames(nona_reorder_FITS2_OTUtable))
Otrad_OutputTable <- data.frame(row.names=rownames(nona_reorder_Otrad_OTUtable))
OITS2_OutputTable <- data.frame(row.names=rownames(nona_reorder_OITS2_OTUtable))
j <- 1
k <- 1
l <- 1
m <- 1
n <- 1
o <- 1
for(i in 1:NumNodes){

	if(NodeTable$FromWhich[i] == "BV5taxa"){

		currentname <- as.character(NodeTable$name[i])
		# splitname <- strsplit(currentname, ";")
		# speciesname <- paste(splitname[[1]][2],splitname[[1]][3],sep="_")
		collectedcolumn <- nona_reorder_BacV5_OTUtable[,which(colnames(nona_reorder_BacV5_OTUtable) == currentname)]
		BacV5_OutputTable[, j] <- collectedcolumn
		colnames(BacV5_OutputTable)[j] <- currentname
		j <- j + 1
	}
	if(NodeTable$FromWhich[i] == "BV3taxa"){

		currentname <- as.character(NodeTable$name[i])
		# splitname <- strsplit(currentname, ";")
		# speciesname <- paste(splitname[[1]][2],splitname[[1]][3],sep="_")
		collectedcolumn <- nona_reorder_BacV3_OTUtable[,which(colnames(nona_reorder_BacV3_OTUtable) == currentname)]
		BacV3_OutputTable[, k] <- collectedcolumn
		colnames(BacV3_OutputTable)[k] <- currentname
		k <- k + 1
	}
	if(NodeTable$FromWhich[i] == "Ftradtaxa"){

		currentname <- as.character(NodeTable$name[i])
		# splitname <- strsplit(currentname, ";")
		# speciesname <- paste(splitname[[1]][2],splitname[[1]][3],sep="_")
		collectedcolumn <- nona_reorder_Ftrad_OTUtable[,which(colnames(nona_reorder_Ftrad_OTUtable) == currentname)]
		Ftrad_OutputTable[, l] <- collectedcolumn
		colnames(Ftrad_OutputTable)[l] <- currentname
		l <- l + 1
	}
	if(NodeTable$FromWhich[i] == "FITS2taxa"){

		currentname <- as.character(NodeTable$name[i])
		# splitname <- strsplit(currentname, ";")
		# speciesname <- paste(splitname[[1]][2],splitname[[1]][3],sep="_")
		collectedcolumn <- nona_reorder_FITS2_OTUtable[,which(colnames(nona_reorder_FITS2_OTUtable) == currentname)]
		FITS2_OutputTable[, m] <- collectedcolumn
		colnames(FITS2_OutputTable)[m] <- currentname
		m <- m + 1
	}
	if(NodeTable$FromWhich[i] == "Otradtaxa"){

		currentname <- as.character(NodeTable$name[i])
		# splitname <- strsplit(currentname, ";")
		# speciesname <- paste(splitname[[1]][2],splitname[[1]][3],sep="_")
		collectedcolumn <- nona_reorder_Otrad_OTUtable[,which(colnames(nona_reorder_Otrad_OTUtable) == currentname)]
		Otrad_OutputTable[, n] <- collectedcolumn
		colnames(Otrad_OutputTable)[n] <- currentname
		n <- n + 1
	}
	if(NodeTable$FromWhich[i] == "OITS2taxa"){

		currentname <- as.character(NodeTable$name[i])
		# splitname <- strsplit(currentname, ";")
		# speciesname <- paste(splitname[[1]][2],splitname[[1]][3],sep="_")
		collectedcolumn <- nona_reorder_OITS2_OTUtable[,which(colnames(nona_reorder_OITS2_OTUtable) == currentname)]
		OITS2_OutputTable[, o] <- collectedcolumn
		colnames(OITS2_OutputTable)[o] <- currentname
		o <- o + 1
	}
}

#############################################################
# Finalize Table and Print it
#############################################################

# Add in the name of the replicate to the table
BacV5_OutputTable$Replicate <- nona_reorder_BacV5_Mapfile$Replicate
BacV3_OutputTable$Replicate <- nona_reorder_BacV3_Mapfile$Replicate
Ftrad_OutputTable$Replicate <- nona_reorder_Ftrad_Mapfile$Replicate
FITS2_OutputTable$Replicate <- nona_reorder_FITS2_Mapfile$Replicate
Otrad_OutputTable$Replicate <- nona_reorder_Otrad_Mapfile$Replicate
OITS2_OutputTable$Replicate <- nona_reorder_OITS2_Mapfile$Replicate

# Print the table
WT <- paste(outputFolder,"BacV5_TestTable_dec_jan_1.txt",sep="")
write.table(BacV5_OutputTable, file=WT, sep = "\t")
WT <- paste(outputFolder,"BacV3_TestTable_dec_jan_1.txt",sep="")
write.table(BacV3_OutputTable, file=WT, sep = "\t")
WT <- paste(outputFolder,"Ftrad_TestTable_dec_jan_1.txt",sep="")
write.table(Ftrad_OutputTable, file=WT, sep = "\t")
WT <- paste(outputFolder,"FITS2_TestTable_dec_jan_1.txt",sep="")
write.table(FITS2_OutputTable, file=WT, sep = "\t")
WT <- paste(outputFolder,"Otrad_TestTable_dec_jan_1.txt",sep="")
write.table(Otrad_OutputTable, file=WT, sep = "\t")
WT <- paste(outputFolder,"OITS2_TestTable_dec_jan_1.txt",sep="")
write.table(OITS2_OutputTable, file=WT, sep = "\t")
