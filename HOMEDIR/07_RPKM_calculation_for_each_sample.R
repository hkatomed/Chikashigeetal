# Specify HOMEDIR before run.
HOMEDIR <- "XXX"

library("openxlsx")

SEQIDs <- c("ERR2003520", "ERR2003521", "ERR2003524", "ERR2003525", "ERR2003529", 
	"ERR2003535", "ERR2003536", "ERR2003539", "ERR2003540", "ERR2003543", 
	"ERR2003544", "ERR2003550", "ERR2003551", "ERR2003553", "ERR2003554")

setwd(HOMEDIR)
setwd("RPKM_2")
setwd("RData")
for(i in 1:length(SEQIDs)){
	SEQID <- SEQIDs[i]
	
	if(i == 1){
		filename <- paste(SEQID, "_all_ORFs.RData", sep = "")
		load(file = filename)
		objname <- paste(SEQID, ".alll.ORFs", sep = "")
		temp <- get(objname)
		temp.allPhase <- temp[,c(1:9, 18)]
		temp.phase1 <- temp[,c(1:9, 19)]
		temp.phase2 <- temp[,c(1:9, 20)]
		temp.phase3 <- temp[,c(1:9, 21)]
		colnames(temp.allPhase)[9+i] <- SEQID
		colnames(temp.phase1)[9+i] <- SEQID
		colnames(temp.phase2)[9+i] <- SEQID
		colnames(temp.phase3)[9+i] <- SEQID
	}

	if(i != 1){
		filename <- paste(SEQID, "_all_ORFs.RData", sep = "")
		load(file = filename)
		objname <- paste(SEQID, ".alll.ORFs", sep = "")
		temp <- get(objname)
		temp.allPhase <- cbind(temp.allPhase, temp[, 18])
		temp.phase1 <- cbind(temp.phase1, temp[, 19])
		temp.phase2 <- cbind(temp.phase2, temp[, 20])
		temp.phase3 <- cbind(temp.phase3, temp[, 21])
		colnames(temp.allPhase)[9+i] <- SEQID
		colnames(temp.phase1)[9+i] <- SEQID
		colnames(temp.phase2)[9+i] <- SEQID
		colnames(temp.phase3)[9+i] <- SEQID
	}
}
setwd("../")
dir.create("summary")
setwd("summary")
summary.allPhase <- temp.allPhase
summary.phase1 <- temp.phase1
summary.phase2 <- temp.phase2
summary.phase3 <- temp.phase3
write.xlsx(x = summary.allPhase, file = "summary_allPhase_RPKM.xlsx")
write.xlsx(x = summary.phase1, file = "summary_phase1_RPKM.xlsx")
write.xlsx(x = summary.phase2, file = "summary_phase2_RPKM.xlsx")
write.xlsx(x = summary.phase3, file = "summary_phase3_RPKM.xlsx")

setwd(HOMEDIR)
setwd("RPKM_2")
setwd("RData")
for(i in 1:length(SEQIDs)){
	SEQID <- SEQIDs[i]
	
	if(i == 1){
		filename <- paste(SEQID, "_all_ORFs.RData", sep = "")
		load(file = filename)
		objname <- paste(SEQID, ".alll.ORFs", sep = "")
		temp <- get(objname)
		temp.allPhase <- temp[,c(1:9, 14)]
		temp.phase1 <- temp[,c(1:9, 15)]
		temp.phase2 <- temp[,c(1:9, 16)]
		temp.phase3 <- temp[,c(1:9, 17)]
		colnames(temp.allPhase)[9+i] <- SEQID
		colnames(temp.phase1)[9+i] <- SEQID
		colnames(temp.phase2)[9+i] <- SEQID
		colnames(temp.phase3)[9+i] <- SEQID
	}

	if(i != 1){
		filename <- paste(SEQID, "_all_ORFs.RData", sep = "")
		load(file = filename)
		objname <- paste(SEQID, ".alll.ORFs", sep = "")
		temp <- get(objname)
		temp.allPhase <- cbind(temp.allPhase, temp[, 14])
		temp.phase1 <- cbind(temp.phase1, temp[, 15])
		temp.phase2 <- cbind(temp.phase2, temp[, 16])
		temp.phase3 <- cbind(temp.phase3, temp[, 17])
		colnames(temp.allPhase)[9+i] <- SEQID
		colnames(temp.phase1)[9+i] <- SEQID
		colnames(temp.phase2)[9+i] <- SEQID
		colnames(temp.phase3)[9+i] <- SEQID
	}
}
setwd("../")
setwd("summary")
summary.allPhase <- temp.allPhase
summary.phase1 <- temp.phase1
summary.phase2 <- temp.phase2
summary.phase3 <- temp.phase3
write.xlsx(x = summary.allPhase, file = "summary_allPhase_RPM.xlsx")
write.xlsx(x = summary.phase1, file = "summary_phase1_RPM.xlsx")
write.xlsx(x = summary.phase2, file = "summary_phase2_RPM.xlsx")
write.xlsx(x = summary.phase3, file = "summary_phase3_RPM.xlsx")

setwd(HOMEDIR)
setwd("RPKM_2")
setwd("RData")
for(i in 1:length(SEQIDs)){
	SEQID <- SEQIDs[i]
	
	if(i == 1){
		filename <- paste(SEQID, "_all_ORFs.RData", sep = "")
		load(file = filename)
		objname <- paste(SEQID, ".alll.ORFs", sep = "")
		temp <- get(objname)
		temp.allPhase <- temp[,c(1:9, 10)]
		temp.phase1 <- temp[,c(1:9, 11)]
		temp.phase2 <- temp[,c(1:9, 12)]
		temp.phase3 <- temp[,c(1:9, 13)]
		colnames(temp.allPhase)[9+i] <- SEQID
		colnames(temp.phase1)[9+i] <- SEQID
		colnames(temp.phase2)[9+i] <- SEQID
		colnames(temp.phase3)[9+i] <- SEQID
	}

	if(i != 1){
		filename <- paste(SEQID, "_all_ORFs.RData", sep = "")
		load(file = filename)
		objname <- paste(SEQID, ".alll.ORFs", sep = "")
		temp <- get(objname)
		temp.allPhase <- cbind(temp.allPhase, temp[, 10])
		temp.phase1 <- cbind(temp.phase1, temp[, 11])
		temp.phase2 <- cbind(temp.phase2, temp[, 12])
		temp.phase3 <- cbind(temp.phase3, temp[, 13])
		colnames(temp.allPhase)[9+i] <- SEQID
		colnames(temp.phase1)[9+i] <- SEQID
		colnames(temp.phase2)[9+i] <- SEQID
		colnames(temp.phase3)[9+i] <- SEQID
	}
}
setwd("../")
setwd("summary")
summary.allPhase <- temp.allPhase
summary.phase1 <- temp.phase1
summary.phase2 <- temp.phase2
summary.phase3 <- temp.phase3
write.xlsx(x = summary.allPhase, file = "summary_allPhase_reads.xlsx")
write.xlsx(x = summary.phase1, file = "summary_phase1_reads.xlsx")
write.xlsx(x = summary.phase2, file = "summary_phase2_reads.xlsx")
write.xlsx(x = summary.phase3, file = "summary_phase3_reads.xlsx")

setwd(HOMEDIR)
mapping.summary <- data.frame(sample = character(length = length(SEQIDs)), 
			total = integer(length = length(SEQIDs)), 
			mapped = integer(length = length(SEQIDs)), 
			ratio = numeric(length = length(SEQIDs)), 
			stringsAsFactors = FALSE)
for(i in 1:length(SEQIDs)){
	SEQID <- SEQIDs[i]
	mapping.summary$sample[i] <- SEQID
	setwd("totalReads_2")
	filename <- paste(SEQID, "_totalReads.txt", sep = "")
	mapping.summary$total[i] <- as.integer(readLines(con = filename))
	setwd("../mappedReads_2")
	filename <- paste(SEQID, "_mappedReads.txt", sep = "")
	mapping.summary$mapped[i] <- as.integer(readLines(con = filename))
	mapping.summary$ratio[i] <- mapping.summary$mapped[i]/mapping.summary$total[i]
	setwd("../")
}
setwd("RPKM_2")
dir.create("bwa_mapping")
setwd("bwa_mapping")
write.xlsx(x = mapping.summary, file = "mapping_summary.xlsx")

