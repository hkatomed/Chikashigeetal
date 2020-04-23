# Specify HOMEDIR before run.
HOMEDIR <- "XXX"


library("openxlsx")
setwd("HOMEDIR")
setwd("RPKM_2")

setwd("summary")
summary.allPhase <- read.xlsx(xlsxFile = "summary_allPhase_RPKM.xlsx")
summary.phase1 <- read.xlsx(xlsxFile = "summary_phase1_RPKM.xlsx")

setwd("HOMEDIR")
setwd("RPKM_2")
setwd("bwa_mapping")
mapping.summary <- read.xlsx(xlsxFile = "mapping_summary.xlsx")

summary.allPhase$WTwo3AT <- as.numeric(0)
summary.allPhase$WTwith3AT <- as.numeric(0)
summary.allPhase$gcn2wo3AT <- as.numeric(0)
summary.allPhase$gcn2with3AT <- as.numeric(0)
summary.allPhase$WTratio3AT <- as.numeric(0)
summary.allPhase$gcn2ratio3AT <- as.numeric(0)

for(i in 1:nrow(summary.allPhase)){
	summary.allPhase$WTwo3AT[i] <- 
	(summary.allPhase$ERR2003543[i] * mapping.summary$mapped[mapping.summary$sample == "ERR2003543"]+
	summary.allPhase$ERR2003544[i] * mapping.summary$mapped[mapping.summary$sample == "ERR2003544"]+
	summary.allPhase$ERR2003550[i] * mapping.summary$mapped[mapping.summary$sample == "ERR2003550"]+
	summary.allPhase$ERR2003551[i] * mapping.summary$mapped[mapping.summary$sample == "ERR2003551"]+
	summary.allPhase$ERR2003553[i] * mapping.summary$mapped[mapping.summary$sample == "ERR2003553"]+
	summary.allPhase$ERR2003554[i] * mapping.summary$mapped[mapping.summary$sample == "ERR2003554"])/
	(mapping.summary$mapped[mapping.summary$sample == "ERR2003543"]+
	mapping.summary$mapped[mapping.summary$sample == "ERR2003544"]+
	mapping.summary$mapped[mapping.summary$sample == "ERR2003550"]+
	mapping.summary$mapped[mapping.summary$sample == "ERR2003551"]+
	mapping.summary$mapped[mapping.summary$sample == "ERR2003553"]+
	mapping.summary$mapped[mapping.summary$sample == "ERR2003554"])

	summary.allPhase$WTwith3AT[i] <- 
	(summary.allPhase$ERR2003529[i] * mapping.summary$mapped[mapping.summary$sample == "ERR2003529"]+
	summary.allPhase$ERR2003535[i] * mapping.summary$mapped[mapping.summary$sample == "ERR2003535"]+
	summary.allPhase$ERR2003536[i] * mapping.summary$mapped[mapping.summary$sample == "ERR2003536"]+
	summary.allPhase$ERR2003539[i] * mapping.summary$mapped[mapping.summary$sample == "ERR2003539"]+
	summary.allPhase$ERR2003540[i] * mapping.summary$mapped[mapping.summary$sample == "ERR2003540"])/
	(mapping.summary$mapped[mapping.summary$sample == "ERR2003529"]+
	mapping.summary$mapped[mapping.summary$sample == "ERR2003535"]+
	mapping.summary$mapped[mapping.summary$sample == "ERR2003536"]+
	mapping.summary$mapped[mapping.summary$sample == "ERR2003539"]+
	mapping.summary$mapped[mapping.summary$sample == "ERR2003540"])

	summary.allPhase$gcn2wo3AT[i] <- 
	(summary.allPhase$ERR2003524[i] * mapping.summary$mapped[mapping.summary$sample == "ERR2003524"]+
	summary.allPhase$ERR2003525[i] * mapping.summary$mapped[mapping.summary$sample == "ERR2003525"])/
	(mapping.summary$mapped[mapping.summary$sample == "ERR2003524"]+
	mapping.summary$mapped[mapping.summary$sample == "ERR2003525"])

	summary.allPhase$gcn2with3AT[i] <- 
	(summary.allPhase$ERR2003520[i] * mapping.summary$mapped[mapping.summary$sample == "ERR2003520"]+
	summary.allPhase$ERR2003521[i] * mapping.summary$mapped[mapping.summary$sample == "ERR2003521"])/
	(mapping.summary$mapped[mapping.summary$sample == "ERR2003520"]+
	mapping.summary$mapped[mapping.summary$sample == "ERR2003521"])

	summary.allPhase$WTratio3AT[i] <- summary.allPhase$WTwith3AT[i]/summary.allPhase$WTwo3AT[i]
	summary.allPhase$gcn2ratio3AT[i] <- summary.allPhase$gcn2with3AT[i]/summary.allPhase$gcn2wo3AT[i]
}

summary.phase1$WTwo3AT <- as.numeric(0)
summary.phase1$WTwith3AT <- as.numeric(0)
summary.phase1$gcn2wo3AT <- as.numeric(0)
summary.phase1$gcn2with3AT <- as.numeric(0)
summary.phase1$WTratio3AT <- as.numeric(0)
summary.phase1$gcn2ratio3AT <- as.numeric(0)

for(i in 1:nrow(summary.phase1)){
	summary.phase1$WTwo3AT[i] <- 
	(summary.phase1$ERR2003543[i] * mapping.summary$mapped[mapping.summary$sample == "ERR2003543"]+
	summary.phase1$ERR2003544[i] * mapping.summary$mapped[mapping.summary$sample == "ERR2003544"]+
	summary.phase1$ERR2003550[i] * mapping.summary$mapped[mapping.summary$sample == "ERR2003550"]+
	summary.phase1$ERR2003551[i] * mapping.summary$mapped[mapping.summary$sample == "ERR2003551"]+
	summary.phase1$ERR2003553[i] * mapping.summary$mapped[mapping.summary$sample == "ERR2003553"]+
	summary.phase1$ERR2003554[i] * mapping.summary$mapped[mapping.summary$sample == "ERR2003554"])/
	(mapping.summary$mapped[mapping.summary$sample == "ERR2003543"]+
	mapping.summary$mapped[mapping.summary$sample == "ERR2003544"]+
	mapping.summary$mapped[mapping.summary$sample == "ERR2003550"]+
	mapping.summary$mapped[mapping.summary$sample == "ERR2003551"]+
	mapping.summary$mapped[mapping.summary$sample == "ERR2003553"]+
	mapping.summary$mapped[mapping.summary$sample == "ERR2003554"])

	summary.phase1$WTwith3AT[i] <- 
	(summary.phase1$ERR2003529[i] * mapping.summary$mapped[mapping.summary$sample == "ERR2003529"]+
	summary.phase1$ERR2003535[i] * mapping.summary$mapped[mapping.summary$sample == "ERR2003535"]+
	summary.phase1$ERR2003536[i] * mapping.summary$mapped[mapping.summary$sample == "ERR2003536"]+
	summary.phase1$ERR2003539[i] * mapping.summary$mapped[mapping.summary$sample == "ERR2003539"]+
	summary.phase1$ERR2003540[i] * mapping.summary$mapped[mapping.summary$sample == "ERR2003540"])/
	(mapping.summary$mapped[mapping.summary$sample == "ERR2003529"]+
	mapping.summary$mapped[mapping.summary$sample == "ERR2003535"]+
	mapping.summary$mapped[mapping.summary$sample == "ERR2003536"]+
	mapping.summary$mapped[mapping.summary$sample == "ERR2003539"]+
	mapping.summary$mapped[mapping.summary$sample == "ERR2003540"])

	summary.phase1$gcn2wo3AT[i] <- 
	(summary.phase1$ERR2003524[i] * mapping.summary$mapped[mapping.summary$sample == "ERR2003524"]+
	summary.phase1$ERR2003525[i] * mapping.summary$mapped[mapping.summary$sample == "ERR2003525"])/
	(mapping.summary$mapped[mapping.summary$sample == "ERR2003524"]+
	mapping.summary$mapped[mapping.summary$sample == "ERR2003525"])

	summary.phase1$gcn2with3AT[i] <- 
	(summary.phase1$ERR2003520[i] * mapping.summary$mapped[mapping.summary$sample == "ERR2003520"]+
	summary.phase1$ERR2003521[i] * mapping.summary$mapped[mapping.summary$sample == "ERR2003521"])/
	(mapping.summary$mapped[mapping.summary$sample == "ERR2003520"]+
	mapping.summary$mapped[mapping.summary$sample == "ERR2003521"])

	summary.phase1$WTratio3AT[i] <- summary.phase1$WTwith3AT[i]/summary.phase1$WTwo3AT[i]
	summary.phase1$gcn2ratio3AT[i] <- summary.phase1$gcn2with3AT[i]/summary.phase1$gcn2wo3AT[i]
}

summary.allPhase$WTwo3ATvsCDS <- as.numeric(0)
summary.allPhase$WTwith3ATvsCDS <- as.numeric(0)
summary.allPhase$gcn2wo3ATvsCDS <- as.numeric(0)
summary.allPhase$gcn2with3ATvsCDS <- as.numeric(0)
summary.allPhase$WTwo3AT_CDS <- as.numeric(0)	
summary.allPhase$WTwith3AT_CDS <- as.numeric(0)	
summary.allPhase$gcn2wo3AT_CDS <- as.numeric(0)	
summary.allPhase$gcn2with3AT_CDS <- as.numeric(0)
for(i in 1:nrow(summary.allPhase)){
	tid <- summary.allPhase$transcript_id[i]
	gene_id <- substr(tid, start = 1, stop = nchar(tid) - 2)
	CDS.name <- paste(gene_id, "_CDS", sep = "")
	temp <- subset(summary.allPhase, name == CDS.name)
	summary.allPhase$WTwo3ATvsCDS[i] <- summary.allPhase$WTwo3AT[i]/temp$WTwo3AT[1]
	summary.allPhase$WTwith3ATvsCDS[i] <- summary.allPhase$WTwith3AT[i]/temp$WTwith3AT[1]
	summary.allPhase$gcn2wo3ATvsCDS[i] <- summary.allPhase$gcn2wo3AT[i]/temp$gcn2wo3AT[1]
	summary.allPhase$gcn2with3ATvsCDS[i] <- summary.allPhase$gcn2with3AT[i]/temp$gcn2with3AT[1]
	summary.allPhase$WTwo3AT_CDS[i] <- temp$WTwo3AT[1]
	summary.allPhase$WTwith3AT_CDS[i] <- temp$WTwith3AT[1]
	summary.allPhase$gcn2wo3AT_CDS[i] <- temp$gcn2wo3AT[1]
	summary.allPhase$gcn2with3AT_CDS[i] <- temp$gcn2with3AT[1]
}

summary.phase1$WTwo3ATvsCDS <- as.numeric(0)
summary.phase1$WTwith3ATvsCDS <- as.numeric(0)
summary.phase1$gcn2wo3ATvsCDS <- as.numeric(0)
summary.phase1$gcn2with3ATvsCDS <- as.numeric(0)
summary.phase1$WTwo3AT_CDS <- as.numeric(0)	
summary.phase1$WTwith3AT_CDS <- as.numeric(0)	
summary.phase1$gcn2wo3AT_CDS <- as.numeric(0)	
summary.phase1$gcn2with3AT_CDS <- as.numeric(0)	
for(i in 1:nrow(summary.phase1)){
	tid <- summary.phase1$transcript_id[i]
	gene_id <- substr(tid, start = 1, stop = nchar(tid) - 2)
	CDS.name <- paste(gene_id, "_CDS", sep = "")
	temp <- subset(summary.phase1, name == CDS.name)
	summary.phase1$WTwo3ATvsCDS[i] <- summary.phase1$WTwo3AT[i]/temp$WTwo3AT[1]
	summary.phase1$WTwith3ATvsCDS[i] <- summary.phase1$WTwith3AT[i]/temp$WTwith3AT[1]
	summary.phase1$gcn2wo3ATvsCDS[i] <- summary.phase1$gcn2wo3AT[i]/temp$gcn2wo3AT[1]
	summary.phase1$gcn2with3ATvsCDS[i] <- summary.phase1$gcn2with3AT[i]/temp$gcn2with3AT[1]
	summary.phase1$WTwo3AT_CDS[i] <- temp$WTwo3AT[1]
	summary.phase1$WTwith3AT_CDS[i] <- temp$WTwith3AT[1]
	summary.phase1$gcn2wo3AT_CDS[i] <- temp$gcn2wo3AT[1]
	summary.phase1$gcn2with3AT_CDS[i] <- temp$gcn2with3AT[1]
}

setwd("HOMEDIR")
setwd("RPKM_2")
setwd("summary")
write.xlsx(x = summary.allPhase, file = "summary_allPhase_RPKM_vsCDS_v2.xlsx")
write.xlsx(x = summary.phase1, file = "summary_phase1_RPKM_vsCDS_v2.xlsx")


setwd("HOMEDIR")
setwd("RPKM_2")
setwd("summary")
summary.phase2 <- read.xlsx(xlsxFile = "summary_phase2_RPKM.xlsx")
summary.phase3 <- read.xlsx(xlsxFile = "summary_phase3_RPKM.xlsx")

setwd("HOMEDIR")
setwd("RPKM_2")
setwd("bwa_mapping")
mapping.summary <- read.xlsx(xlsxFile = "mapping_summary.xlsx")

summary.phase2$WTwo3AT <- as.numeric(0)
summary.phase2$WTwith3AT <- as.numeric(0)
summary.phase2$gcn2wo3AT <- as.numeric(0)
summary.phase2$gcn2with3AT <- as.numeric(0)
summary.phase2$WTratio3AT <- as.numeric(0)
summary.phase2$gcn2ratio3AT <- as.numeric(0)

for(i in 1:nrow(summary.phase2)){
	summary.phase2$WTwo3AT[i] <- 
	(summary.phase2$ERR2003543[i] * mapping.summary$mapped[mapping.summary$sample == "ERR2003543"]+
	summary.phase2$ERR2003544[i] * mapping.summary$mapped[mapping.summary$sample == "ERR2003544"]+
	summary.phase2$ERR2003550[i] * mapping.summary$mapped[mapping.summary$sample == "ERR2003550"]+
	summary.phase2$ERR2003551[i] * mapping.summary$mapped[mapping.summary$sample == "ERR2003551"]+
	summary.phase2$ERR2003553[i] * mapping.summary$mapped[mapping.summary$sample == "ERR2003553"]+
	summary.phase2$ERR2003554[i] * mapping.summary$mapped[mapping.summary$sample == "ERR2003554"])/
	(mapping.summary$mapped[mapping.summary$sample == "ERR2003543"]+
	mapping.summary$mapped[mapping.summary$sample == "ERR2003544"]+
	mapping.summary$mapped[mapping.summary$sample == "ERR2003550"]+
	mapping.summary$mapped[mapping.summary$sample == "ERR2003551"]+
	mapping.summary$mapped[mapping.summary$sample == "ERR2003553"]+
	mapping.summary$mapped[mapping.summary$sample == "ERR2003554"])

	summary.phase2$WTwith3AT[i] <- 
	(summary.phase2$ERR2003529[i] * mapping.summary$mapped[mapping.summary$sample == "ERR2003529"]+
	summary.phase2$ERR2003535[i] * mapping.summary$mapped[mapping.summary$sample == "ERR2003535"]+
	summary.phase2$ERR2003536[i] * mapping.summary$mapped[mapping.summary$sample == "ERR2003536"]+
	summary.phase2$ERR2003539[i] * mapping.summary$mapped[mapping.summary$sample == "ERR2003539"]+
	summary.phase2$ERR2003540[i] * mapping.summary$mapped[mapping.summary$sample == "ERR2003540"])/
	(mapping.summary$mapped[mapping.summary$sample == "ERR2003529"]+
	mapping.summary$mapped[mapping.summary$sample == "ERR2003535"]+
	mapping.summary$mapped[mapping.summary$sample == "ERR2003536"]+
	mapping.summary$mapped[mapping.summary$sample == "ERR2003539"]+
	mapping.summary$mapped[mapping.summary$sample == "ERR2003540"])

	summary.phase2$gcn2wo3AT[i] <- 
	(summary.phase2$ERR2003524[i] * mapping.summary$mapped[mapping.summary$sample == "ERR2003524"]+
	summary.phase2$ERR2003525[i] * mapping.summary$mapped[mapping.summary$sample == "ERR2003525"])/
	(mapping.summary$mapped[mapping.summary$sample == "ERR2003524"]+
	mapping.summary$mapped[mapping.summary$sample == "ERR2003525"])

	summary.phase2$gcn2with3AT[i] <- 
	(summary.phase2$ERR2003520[i] * mapping.summary$mapped[mapping.summary$sample == "ERR2003520"]+
	summary.phase2$ERR2003521[i] * mapping.summary$mapped[mapping.summary$sample == "ERR2003521"])/
	(mapping.summary$mapped[mapping.summary$sample == "ERR2003520"]+
	mapping.summary$mapped[mapping.summary$sample == "ERR2003521"])

	summary.phase2$WTratio3AT[i] <- summary.phase2$WTwith3AT[i]/summary.phase2$WTwo3AT[i]
	summary.phase2$gcn2ratio3AT[i] <- summary.phase2$gcn2with3AT[i]/summary.phase2$gcn2wo3AT[i]
}

summary.phase3$WTwo3AT <- as.numeric(0)
summary.phase3$WTwith3AT <- as.numeric(0)
summary.phase3$gcn2wo3AT <- as.numeric(0)
summary.phase3$gcn2with3AT <- as.numeric(0)
summary.phase3$WTratio3AT <- as.numeric(0)
summary.phase3$gcn2ratio3AT <- as.numeric(0)

for(i in 1:nrow(summary.phase3)){
	summary.phase3$WTwo3AT[i] <- 
	(summary.phase3$ERR2003543[i] * mapping.summary$mapped[mapping.summary$sample == "ERR2003543"]+
	summary.phase3$ERR2003544[i] * mapping.summary$mapped[mapping.summary$sample == "ERR2003544"]+
	summary.phase3$ERR2003550[i] * mapping.summary$mapped[mapping.summary$sample == "ERR2003550"]+
	summary.phase3$ERR2003551[i] * mapping.summary$mapped[mapping.summary$sample == "ERR2003551"]+
	summary.phase3$ERR2003553[i] * mapping.summary$mapped[mapping.summary$sample == "ERR2003553"]+
	summary.phase3$ERR2003554[i] * mapping.summary$mapped[mapping.summary$sample == "ERR2003554"])/
	(mapping.summary$mapped[mapping.summary$sample == "ERR2003543"]+
	mapping.summary$mapped[mapping.summary$sample == "ERR2003544"]+
	mapping.summary$mapped[mapping.summary$sample == "ERR2003550"]+
	mapping.summary$mapped[mapping.summary$sample == "ERR2003551"]+
	mapping.summary$mapped[mapping.summary$sample == "ERR2003553"]+
	mapping.summary$mapped[mapping.summary$sample == "ERR2003554"])

	summary.phase3$WTwith3AT[i] <- 
	(summary.phase3$ERR2003529[i] * mapping.summary$mapped[mapping.summary$sample == "ERR2003529"]+
	summary.phase3$ERR2003535[i] * mapping.summary$mapped[mapping.summary$sample == "ERR2003535"]+
	summary.phase3$ERR2003536[i] * mapping.summary$mapped[mapping.summary$sample == "ERR2003536"]+
	summary.phase3$ERR2003539[i] * mapping.summary$mapped[mapping.summary$sample == "ERR2003539"]+
	summary.phase3$ERR2003540[i] * mapping.summary$mapped[mapping.summary$sample == "ERR2003540"])/
	(mapping.summary$mapped[mapping.summary$sample == "ERR2003529"]+
	mapping.summary$mapped[mapping.summary$sample == "ERR2003535"]+
	mapping.summary$mapped[mapping.summary$sample == "ERR2003536"]+
	mapping.summary$mapped[mapping.summary$sample == "ERR2003539"]+
	mapping.summary$mapped[mapping.summary$sample == "ERR2003540"])

	summary.phase3$gcn2wo3AT[i] <- 
	(summary.phase3$ERR2003524[i] * mapping.summary$mapped[mapping.summary$sample == "ERR2003524"]+
	summary.phase3$ERR2003525[i] * mapping.summary$mapped[mapping.summary$sample == "ERR2003525"])/
	(mapping.summary$mapped[mapping.summary$sample == "ERR2003524"]+
	mapping.summary$mapped[mapping.summary$sample == "ERR2003525"])

	summary.phase3$gcn2with3AT[i] <- 
	(summary.phase3$ERR2003520[i] * mapping.summary$mapped[mapping.summary$sample == "ERR2003520"]+
	summary.phase3$ERR2003521[i] * mapping.summary$mapped[mapping.summary$sample == "ERR2003521"])/
	(mapping.summary$mapped[mapping.summary$sample == "ERR2003520"]+
	mapping.summary$mapped[mapping.summary$sample == "ERR2003521"])

	summary.phase3$WTratio3AT[i] <- summary.phase3$WTwith3AT[i]/summary.phase3$WTwo3AT[i]
	summary.phase3$gcn2ratio3AT[i] <- summary.phase3$gcn2with3AT[i]/summary.phase3$gcn2wo3AT[i]
}

summary.phase2$WTwo3ATvsCDS <- as.numeric(0)
summary.phase2$WTwith3ATvsCDS <- as.numeric(0)
summary.phase2$gcn2wo3ATvsCDS <- as.numeric(0)
summary.phase2$gcn2with3ATvsCDS <- as.numeric(0)
summary.phase2$WTwo3AT_CDS <- as.numeric(0)	
summary.phase2$WTwith3AT_CDS <- as.numeric(0)	
summary.phase2$gcn2wo3AT_CDS <- as.numeric(0)	
summary.phase2$gcn2with3AT_CDS <- as.numeric(0)	
for(i in 1:nrow(summary.phase2)){
	tid <- summary.phase2$transcript_id[i]
	gene_id <- substr(tid, start = 1, stop = nchar(tid) - 2)
	CDS.name <- paste(gene_id, "_CDS", sep = "")
	temp <- subset(summary.phase2, name == CDS.name)
	summary.phase2$WTwo3ATvsCDS[i] <- summary.phase2$WTwo3AT[i]/temp$WTwo3AT[1]
	summary.phase2$WTwith3ATvsCDS[i] <- summary.phase2$WTwith3AT[i]/temp$WTwith3AT[1]
	summary.phase2$gcn2wo3ATvsCDS[i] <- summary.phase2$gcn2wo3AT[i]/temp$gcn2wo3AT[1]
	summary.phase2$gcn2with3ATvsCDS[i] <- summary.phase2$gcn2with3AT[i]/temp$gcn2with3AT[1]
	summary.phase2$WTwo3AT_CDS[i] <- temp$WTwo3AT[1]
	summary.phase2$WTwith3AT_CDS[i] <- temp$WTwith3AT[1]
	summary.phase2$gcn2wo3AT_CDS[i] <- temp$gcn2wo3AT[1]
	summary.phase2$gcn2with3AT_CDS[i] <- temp$gcn2with3AT[1]
}

summary.phase3$WTwo3ATvsCDS <- as.numeric(0)
summary.phase3$WTwith3ATvsCDS <- as.numeric(0)
summary.phase3$gcn2wo3ATvsCDS <- as.numeric(0)
summary.phase3$gcn2with3ATvsCDS <- as.numeric(0)
summary.phase3$WTwo3AT_CDS <- as.numeric(0)	
summary.phase3$WTwith3AT_CDS <- as.numeric(0)	
summary.phase3$gcn2wo3AT_CDS <- as.numeric(0)	
summary.phase3$gcn2with3AT_CDS <- as.numeric(0)	
for(i in 1:nrow(summary.phase3)){
	tid <- summary.phase3$transcript_id[i]
	gene_id <- substr(tid, start = 1, stop = nchar(tid) - 2)
	CDS.name <- paste(gene_id, "_CDS", sep = "")
	temp <- subset(summary.phase3, name == CDS.name)
	summary.phase3$WTwo3ATvsCDS[i] <- summary.phase3$WTwo3AT[i]/temp$WTwo3AT[1]
	summary.phase3$WTwith3ATvsCDS[i] <- summary.phase3$WTwith3AT[i]/temp$WTwith3AT[1]
	summary.phase3$gcn2wo3ATvsCDS[i] <- summary.phase3$gcn2wo3AT[i]/temp$gcn2wo3AT[1]
	summary.phase3$gcn2with3ATvsCDS[i] <- summary.phase3$gcn2with3AT[i]/temp$gcn2with3AT[1]
	summary.phase3$WTwo3AT_CDS[i] <- temp$WTwo3AT[1]
	summary.phase3$WTwith3AT_CDS[i] <- temp$WTwith3AT[1]
	summary.phase3$gcn2wo3AT_CDS[i] <- temp$gcn2wo3AT[1]
	summary.phase3$gcn2with3AT_CDS[i] <- temp$gcn2with3AT[1]
}

setwd("HOMEDIR")
setwd("RPKM_2")
setwd("summary")
write.xlsx(x = summary.phase2, file = "summary_phase2_RPKM_vsCDS_v2.xlsx")
write.xlsx(x = summary.phase3, file = "summary_phase3_RPKM_vsCDS_v2.xlsx")

