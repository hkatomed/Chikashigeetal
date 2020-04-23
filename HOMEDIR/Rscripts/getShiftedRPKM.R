SEQID <- commandArgs(trailingOnly=TRUE)[1]
cat("SEQID:", SEQID, "\n")
wd <- getwd()

library("openxlsx")
setwd("pombase/RData")
load(file = "all_ORFs.RData")
load(file = "CDS_ATGSTOP.RData")
setwd("../../")

for(i in 1:nrow(all.ORFs)){
	tid <- all.ORFs$transcript_id[i]
	name <- all.ORFs$name[i]
	gene_id <- substr(tid, start = 1, stop = nchar(tid)-2)
	if(name == gene_id) all.ORFs$name[i] <- paste(gene_id, "_CDS", sep = "")
	if(name != gene_id) all.ORFs$name[i] <- paste(gene_id, "_", name, sep = "")
}

CDS.ATGSTOP$mRNA.length <- 0
for(i in 1:nrow(CDS.ATGSTOP)){
	start <- 1
	end <- CDS.ATGSTOP$STOP[i]
	if(CDS.ATGSTOP$UTR5[i] == TRUE) start <- CDS.ATGSTOP$UTR5.start[i]
	if(CDS.ATGSTOP$UTR3[i] == TRUE) end <- CDS.ATGSTOP$UTR3.end[i]
	CDS.ATGSTOP$mRNA.length[i] <- end - start + 1
}

setwd("totalReads_2")
filename <- paste(SEQID, "_totalReads.txt", sep = "")
totalReads <- as.integer(readLines(con = filename))
setwd("../mappedReads_2")
filename <- paste(SEQID, "_mappedReads.txt", sep = "")
mappedReads <- as.integer(readLines(con = filename))
setwd("../bedGraph_5end_2")
filename <- paste(SEQID, "_Sense.bedGraph", sep = "")
bedGraph <- read.table(file = filename, stringsAsFactors = FALSE)
colnames(bedGraph) <- c("transcript_id", "pos", "score.5end")
setwd("../")

all.ORFs$uORFseq <- NULL
all.ORFs$allPhase.score <- as.integer(NA)
all.ORFs$phase1.score <- as.integer(NA)
all.ORFs$phase2.score <- as.integer(NA)
all.ORFs$phase3.score <- as.integer(NA)
all.ORFs$allPhase.RPM <- as.integer(NA)
all.ORFs$phase1.RPM <- as.integer(NA)
all.ORFs$phase2.RPM <- as.integer(NA)
all.ORFs$phase3.RPM <- as.integer(NA)
all.ORFs$allPhase.RPKM <- as.integer(NA)
all.ORFs$phase1.RPKM <- as.integer(NA)
all.ORFs$phase2.RPKM <- as.integer(NA)
all.ORFs$phase3.RPKM <- as.integer(NA)
cat("nrow(all.ORFs):", nrow(all.ORFs), "\n", sep = "")

for(i in 1:nrow(all.ORFs)){
	tid <- all.ORFs$transcript_id[i]
	start <- all.ORFs$start[i]
	end <- all.ORFs$end[i]
	mRNA.length <- CDS.ATGSTOP$mRNA.length[CDS.ATGSTOP$transcript_id == tid]
	tid.bedGraph <- subset(bedGraph, transcript_id == tid)
	allPhase.start <- start - 12
	allPhase.end <- end - 12
	phase1.index <- seq(start - 12, allPhase.end, 3)
	phase2.index <- seq(start - 12+1, allPhase.end, 3)
	phase3.index <- seq(start - 12+2, allPhase.end, 3)
	if(allPhase.start < 1)	allPhase.start <- 1
	phase1.index <- phase1.index[phase1.index>0]
	phase2.index <- phase2.index[phase2.index>0]
	phase3.index <- phase3.index[phase3.index>0]
	
	all.ORFs$allPhase.score[i] <- sum(subset(tid.bedGraph, 
		pos >= allPhase.start & pos <= allPhase.end)$score)
	all.ORFs$phase1.score[i] <- sum(tid.bedGraph$score[phase1.index])
	all.ORFs$phase2.score[i] <- sum(tid.bedGraph$score[phase2.index])
	all.ORFs$phase3.score[i] <- sum(tid.bedGraph$score[phase3.index])
	all.ORFs$allPhase.RPM[i] <- all.ORFs$allPhase.score[i] * 1000000/mappedReads
	all.ORFs$phase1.RPM[i] <- all.ORFs$phase1.score[i] * 1000000/mappedReads
	all.ORFs$phase2.RPM[i] <- all.ORFs$phase2.score[i] * 1000000/mappedReads
	all.ORFs$phase3.RPM[i] <- all.ORFs$phase3.score[i] * 1000000/mappedReads
	all.ORFs$allPhase.RPKM[i] <- all.ORFs$allPhase.RPM[i] /(end - start + 1) *1000
	all.ORFs$phase1.RPKM[i] <- all.ORFs$phase1.RPM[i] /(end - start + 1) *1000
	all.ORFs$phase2.RPKM[i] <- all.ORFs$phase2.RPM[i] /(end - start + 1) *1000
	all.ORFs$phase3.RPKM[i] <- all.ORFs$phase3.RPM[i] /(end - start + 1) *1000
	if(i%%1000 == 0)	cat(",", i, sep = "")
}

setwd("RPKM_2")
objname <- paste(SEQID, ".alll.ORFs", sep = "")
assign(objname, all.ORFs)
if(length(grep(pattern = "RData", list.files())) == 0) dir.create("RData")
setwd("RData")
filename <- paste(SEQID, "_all_ORFs.RData", sep = "")
save(list = objname, file = filename)
setwd("../")
if(length(grep(pattern = "TXT", list.files())) == 0) dir.create("TXT")
setwd("TXT")
filename <- paste(SEQID, "_all_ORFs.txt", sep = "")
write.table(all.ORFs, file = filename, quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
setwd("../")
if(length(grep(pattern = "TXT", list.files())) == 0) dir.create("XLSX")
setwd("XLSX")
filename <- paste(SEQID, "_all_ORFs.xlsx", sep = "")
write.xlsx(x = all.ORFs, file = filename)
setwd(wd)
