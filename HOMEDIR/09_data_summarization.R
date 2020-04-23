## RPKM values in the output bundledDF_20200312.xlsx are used in Figs. 4AB, S5 and S6

# Specify HOMEDIR before run.
HOMEDIR <- "XXX"


library(openxlsx)
library(Biostrings)
library(beeswarm)
library(ggplot2)
library(ggseqlogo)

setwd(HOMEDIR)
setwd("RPKM_2")
setwd("summary")
sumary_phase1_RPKM_vsCDS_v2 <- read.xlsx(xlsxFile = "summary_phase1_RPKM_vsCDS_v2.xlsx")
sumary_phase2_RPKM_vsCDS_v2 <- read.xlsx(xlsxFile = "summary_phase2_RPKM_vsCDS_v2.xlsx")
sumary_phase3_RPKM_vsCDS_v2 <- read.xlsx(xlsxFile = "summary_phase3_RPKM_vsCDS_v2.xlsx")

sumary_phase1_RPKM_vsCDS_v2$phase2WTwo3AT <- sumary_phase2_RPKM_vsCDS_v2$WTwo3AT
sumary_phase1_RPKM_vsCDS_v2$phase3WTwo3AT <- sumary_phase3_RPKM_vsCDS_v2$WTwo3AT

sumary_phase1_RPKM_vsCDS_v2$phase2WTwith3AT <- sumary_phase2_RPKM_vsCDS_v2$WTwith3AT
sumary_phase1_RPKM_vsCDS_v2$phase3WTwith3AT <- sumary_phase3_RPKM_vsCDS_v2$WTwith3AT

sumary_phase1_RPKM_vsCDS_v2$phase2gcn2wo3AT <- sumary_phase2_RPKM_vsCDS_v2$gcn2wo3AT
sumary_phase1_RPKM_vsCDS_v2$phase3gcn2wo3AT <- sumary_phase3_RPKM_vsCDS_v2$gcn2wo3AT

sumary_phase1_RPKM_vsCDS_v2$phase2gcn2with3AT <- sumary_phase2_RPKM_vsCDS_v2$gcn2with3AT
sumary_phase1_RPKM_vsCDS_v2$phase3gcn2with3AT <- sumary_phase3_RPKM_vsCDS_v2$gcn2with3AT

setwd(HOMEDIR)
setwd("pombase")
load(file = "RData/splicedSeqs2.RData")
load(file = "RData/splicedSeqsM2.RData")

sumary_phase1_RPKM_vsCDS_v2$ratioWT <- sumary_phase1_RPKM_vsCDS_v2$WTwith3ATvsCDS / sumary_phase1_RPKM_vsCDS_v2$WTwo3ATvsCDS
sumary_phase1_RPKM_vsCDS_v2$ratiogcn2 <- sumary_phase1_RPKM_vsCDS_v2$gcn2with3ATvsCDS / sumary_phase1_RPKM_vsCDS_v2$gcn2wo3ATvsCDS
sumary_phase1_RPKM_vsCDS_v2$ratioWT[is.nan(sumary_phase1_RPKM_vsCDS_v2$ratioWT) == TRUE] <- NA
sumary_phase1_RPKM_vsCDS_v2$ratiogcn2[is.nan(sumary_phase1_RPKM_vsCDS_v2$ratiogcn2) == TRUE] <- NA
sumary_phase1_RPKM_vsCDS_v2$ratioWT[sumary_phase1_RPKM_vsCDS_v2$ratioWT == Inf] <- NA
sumary_phase1_RPKM_vsCDS_v2$ratiogcn2[sumary_phase1_RPKM_vsCDS_v2$ratiogcn2 == Inf] <- NA

sumary_phase1_RPKM_vsCDS_v2$ratioWo3AT <- sumary_phase1_RPKM_vsCDS_v2$gcn2wo3ATvsCDS / sumary_phase1_RPKM_vsCDS_v2$WTwo3ATvsCDS
sumary_phase1_RPKM_vsCDS_v2$ratioWith3AT <- sumary_phase1_RPKM_vsCDS_v2$gcn2with3ATvsCDS / sumary_phase1_RPKM_vsCDS_v2$WTwith3ATvsCDS
sumary_phase1_RPKM_vsCDS_v2$ratioWo3AT[is.nan(sumary_phase1_RPKM_vsCDS_v2$ratioWo3AT) == TRUE] <- NA
sumary_phase1_RPKM_vsCDS_v2$ratioWith3AT[is.nan(sumary_phase1_RPKM_vsCDS_v2$ratioWith3AT) == TRUE] <- NA
sumary_phase1_RPKM_vsCDS_v2$ratioWo3AT[sumary_phase1_RPKM_vsCDS_v2$ratioWo3AT == Inf] <- NA
sumary_phase1_RPKM_vsCDS_v2$ratioWith3AT[sumary_phase1_RPKM_vsCDS_v2$ratioWith3AT == Inf] <- NA

sumary_phase1_RPKM_vsCDS_v2$masked <- FALSE

checkMasked <- function(tid, start, end){
	Seq2 <- splicedSeqs2[[tid]][start:end]
	SeqM <- splicedSeqsM[[tid]][start:end]
	return(SeqM != Seq2)
}

sumary_phase1_RPKM_vsCDS_v2$masked <- unlist(mcMap(f = checkMasked, 
		tid = sumary_phase1_RPKM_vsCDS_v2$transcript_id, 
		start = sumary_phase1_RPKM_vsCDS_v2$start, 
		end = sumary_phase1_RPKM_vsCDS_v2$end, USE.NAMES = FALSE, mc.cores = 24))

maskedTid <- unique(subset(sumary_phase1_RPKM_vsCDS_v2, CDS == TRUE & masked == TRUE)$transcript_id)

sumary_phase1_RPKM_vsCDS_v2$masked2 <- FALSE
for(i in 1:length(maskedTid)){
	sumary_phase1_RPKM_vsCDS_v2$masked2[sumary_phase1_RPKM_vsCDS_v2$transcript_id == maskedTid[i]] <- TRUE
}

sumary_phase1_RPKM_vsCDS_v2$masked3 <- FALSE

checkNinORF <- function(tid, start, end, Nnum = 22){
	# Seq2 <- splicedSeqs2[[tid]][start:end]
	SeqM <- splicedSeqsM[[tid]][start:end]
	PATTERN <- character(length = Nnum)
	for(i in 1:Nnum) PATTERN[i] <- "N"
	PATTERN <- paste(PATTERN, collapse = "", sep = "")
	# > PATTERN
	# [1] "NNNNNNNNNNNNNNNNNNNNNN"
	if(countPattern(pattern = PATTERN, SeqM) > 0) return(TRUE)
	if(countPattern(pattern = PATTERN, SeqM) == 0) return(FALSE)
}

sumary_phase1_RPKM_vsCDS_v2$masked3 <- unlist(mcMap(f = checkNinORF, 
		tid = sumary_phase1_RPKM_vsCDS_v2$transcript_id, 
		start = sumary_phase1_RPKM_vsCDS_v2$start, 
		end = sumary_phase1_RPKM_vsCDS_v2$end, Nnum = 22, 
		USE.NAMES = FALSE, mc.cores = 24))

maskedTid <- unique(subset(sumary_phase1_RPKM_vsCDS_v2, CDS == TRUE & masked3 == TRUE)$transcript_id)

sumary_phase1_RPKM_vsCDS_v2$masked4 <- FALSE
for(i in 1:length(maskedTid)){
	sumary_phase1_RPKM_vsCDS_v2$masked4[sumary_phase1_RPKM_vsCDS_v2$transcript_id == maskedTid[i]] <- TRUE
}

unique.tid <- unique(sumary_phase1_RPKM_vsCDS_v2$transcript_id)

bundledDF <- sumary_phase1_RPKM_vsCDS_v2[0,]
for(i in 1:length(unique.tid)){
	cat(i, ",", sep = "")
	# i <- 3
	tid <- unique.tid[i]
	temp.uORF <- subset(sumary_phase1_RPKM_vsCDS_v2, transcript_id == tid & uORF == TRUE)
	temp.CDS <- subset(sumary_phase1_RPKM_vsCDS_v2, transcript_id == tid & CDS == TRUE)
	
	if(nrow(temp.uORF) > 0){
		temp.uORF$uORFoverlap <- FALSE
		endPrev <- 0
		for(o in 1:nrow(temp.uORF)){
			end <- temp.uORF$end[o]
			if(end == endPrev) temp.uORF$uORFoverlap[o] <- TRUE
			endPrev <- end
		}
		temp.uORF <- subset(temp.uORF, uORFoverlap == FALSE)
		for(o in 1:nrow(temp.uORF)){
			temp.uORF$name[o] <- paste(tid, "_uORF", o, sep = "")
		}
		temp.uORF$uORFoverlap <- NULL
	}
	bundledDF <- rbind(bundledDF, temp.uORF, temp.CDS)
}

setwd(HOMEDIR)
setwd("RPKM_2")
setwd("summary")
write.xlsx(bundledDF, file = "bundledDF_20200312.xlsx")
setwd("../RData")
save(bundledDF, file = "bundledDF_20200312.RData")



