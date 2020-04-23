# Specify HOMEDIR before run.
HOMEDIR <- "XXX"

library(Biostrings)

setwd(HOMEDIR)
setwd("pombase")
gff3 <- read.table(file = "Schizosaccharomyces_pombe_all_chromosomes.gff3.gz", 
			sep = "\t", skip = 1, stringsAsFactors = FALSE)
names(gff3) <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
nrow(gff3)

nrow(subset(gff3, type == "mRNA"))
nrow(subset(gff3, type == "five_prime_UTR"))
nrow(subset(gff3, type == "three_prime_UTR"))
nrow(subset(gff3, type == "CDS"))

get.ID.attrb1 <- function(attrb1){
	out <- strsplit(x = strsplit(x = attrb1, split = ";")[[1]][1], split = "ID=")[[1]][2]
	return(out)
}

get.ID.attrb2 <- function(attrb2){
	temp <- strsplit(x = strsplit(x = strsplit(x = attrb2, split = ";")[[1]][1], 
				split = ":")[[1]][1], split = "ID=")[[1]][2]
	if(substr(temp, start = nchar(temp)-1, stop = nchar(temp)) == ".1"){
		out <- substr(temp, start = 1, stop = nchar(temp)-2)
	}
	if(substr(temp, start = nchar(temp)-1, stop = nchar(temp)) != ".1"){
		out <- "NA"
	}
	return(out)
}

gff3.CDS <- subset(gff3, type == "CDS")
gff3.CDS$gene_id <- "NA"
for(i in 1:nrow(gff3.CDS)){
	gff3.CDS$gene_id[i] <- get.ID.attrb2(gff3.CDS$attributes[i])
}

nrow(gff3.CDS)
length(unique(gff3.CDS$gene_id))

CDS.gene_id <- unique(gff3.CDS$gene_id)
CDS.transcript_id <- character(length = length(CDS.gene_id))
for(i in 1:length(CDS.gene_id)){
	CDS.transcript_id[i] <- paste(CDS.gene_id[i], ".1", sep = "")
}

gff3.gene <- subset(gff3, type == "gene")
gff3.gene$gene_id <- "NA"
for(i in 1:nrow(gff3.gene)){
	gff3.gene$gene_id[i] <- get.ID.attrb1(gff3.gene$attributes[i])
}

nrow(gff3.gene)
length(unique(gff3.gene$gene_id))

gff3.gene$transcript_id <- "NA"
for(i in 1:nrow(gff3.gene)){
	gff3.gene$transcript_id[i] <- paste(gff3.gene$gene_id[i], ".1", sep = "")
}

nrow(subset(gff3, type == "gene"))
nrow(subset(gff3, type == "mRNA"))
nrow(subset(gff3, type == "CDS"))
nrow(subset(gff3, type == "intron"))
nrow(subset(gff3, type == "five_prime_UTR"))
nrow(subset(gff3, type == "three_prime_UTR"))
nrow(subset(gff3, type == "pseudogenic_transcript"))
nrow(subset(gff3, type == "tRNA"))
nrow(subset(gff3, type == "ncRNA"))
nrow(subset(gff3, type == "snRNA"))
nrow(subset(gff3, type == "snoRNA"))
nrow(subset(gff3, type == "rRNA"))

gff3.mRNA <- subset(gff3, type == "mRNA")
gff3.mRNA$gene_id <- "NA"
for(i in 1:nrow(gff3.mRNA)){
	gff3.mRNA$gene_id[i] <- get.ID.attrb2(gff3.mRNA$attributes[i])
}

nrow(gff3.mRNA)
length(unique(gff3.mRNA$gene_id))

mRNA.gene_id <- unique(gff3.mRNA$gene_id)
mRNA.transcript_id <- character(length = length(mRNA.gene_id))
for(i in 1:length(mRNA.gene_id)){
	mRNA.transcript_id[i] <- paste(mRNA.gene_id[i], ".1", sep = "")
}

gff3.pc <- gff3[0,]
for(i in 1:length(mRNA.transcript_id)){
	id <- mRNA.transcript_id[i]
	index <- grep(pattern = id, x = gff3$attributes)
	temp <- gff3[index, ]
	gff3.pc <- rbind(gff3.pc, temp)
	if(i%%100 == 0) cat("i: ", i, "\n", sep = "")
}

gff3.pc$gene_id <- "NA"
for(i in 1:nrow(gff3.pc)){
	gff3.pc$gene_id[i] <- get.ID.attrb2(gff3.pc$attributes[i])
}

gff3.pc$transcript_id <- "NA"
for(i in 1:nrow(gff3.pc)){
	gff3.pc$transcript_id[i] <- paste(gff3.pc$gene_id[i], ".1", sep = "")
}

setwd(HOMEDIR)
setwd("pombase")
sp.genome <- readDNAStringSet(filepath = "Schizosaccharomyces_pombe_all_chromosomes.fa.gz")

chrNames <- names(sp.genome)
for(i in 1:6){
	chrNames[i] <- strsplit(chrNames[i], split = " ")[[1]][1]
}

names(sp.genome) <- chrNames

geneSeqs <- character(length = nrow(gff3.gene))
for(i in 1:nrow(gff3.gene)){
	seqid <- gff3.gene$seqid[i]
	start <- gff3.gene$start[i]
	end <- gff3.gene$end[i]
	strand <- gff3.gene$strand[i]
	Seq <- sp.genome[[seqid]][start:end]
	if(strand == "-") Seq <- reverseComplement(Seq)
	geneSeqs[i] <- as.character(Seq)
}

geneSeqs <- DNAStringSet(geneSeqs)
names(geneSeqs) <- gff3.gene$transcript_id

setwd(HOMEDIR)
setwd("reference_transcriptome")
filename <- "sp_transcripts.fasta"
writeXStringSet(geneSeqs, filepath = filename)

pcSeqs <- character(length = nrow(gff3.pc))
for(i in 1:nrow(gff3.pc)){
	seqid <- gff3.pc$seqid[i]
	start <- gff3.pc$start[i]
	end <- gff3.pc$end[i]
	strand <- gff3.pc$strand[i]
	Seq <- sp.genome[[seqid]][start:end]
	if(strand == "-") Seq <- reverseComplement(Seq)
	pcSeqs[i] <- as.character(Seq)
}

pcSeqs <- DNAStringSet(pcSeqs)
names(pcSeqs) <- gff3.pc$transcript_id

tids <- unique(gff3.pc$transcript_id)
splicedSeqs <- DNAStringSet()
splicedATGs <- integer(length = nrow(gff3.mRNA))
splicedSTOPs <- integer(length = nrow(gff3.mRNA))
for(i in 1:length(tids)){
	temp <- subset(gff3.pc, transcript_id == tids[i])
	temp$start <- as.integer(temp$start)
	temp$end <- as.integer(temp$end)
	temp.mRNA <- subset(temp, type == "mRNA")
	if(nrow(temp.mRNA) > 1) stop("nrow(temp.mRNA) > 1")
	gene.start <- temp.mRNA$start[1]
	gene.end <- temp.mRNA$end[1]
	gene.strand <- temp.mRNA$strand[1]
	gene.length <- gene.end - gene.start + 1
	 
		temp$start <- temp$start - gene.start + 1
		temp$end <- temp$end - gene.start + 1
	
	if(gene.strand == "-"){ 
		start <- temp$start
		temp$start <- gene.length - temp$end + 1
		temp$end <- gene.length - start + 1
	}
	temp.intron <- subset(temp, type == "intron")
	temp.intron <- temp.intron[order(temp.intron$start, decreasing = FALSE),]
	temp.CDS <- subset(temp, type == "CDS")
	temp.CDS <- temp.CDS[order(temp.CDS$start, decreasing = FALSE),]
	
	if(nrow(temp.intron) != 0){
		seq.index <- data.frame(start = integer(length = nrow(temp.intron)+1), 
				end = integer(length = nrow(temp.intron)+1))
		for(j in 1:(nrow(temp.intron)+1)){
			if(j == 1) seq.index$start[j] <- 1
			if(j == 1) seq.index$end[j] <- temp.intron$start[j] - 1
			if(j != 1) seq.index$start[j] <- temp.intron$end[j-1] +1
			if(j != 1) seq.index$end[j] <- temp.intron$start[j] - 1
			if(j == (nrow(temp.intron)+1)) seq.index$start[j] <- temp.intron$end[j-1] +1
			if(j == (nrow(temp.intron)+1)) seq.index$end[j] <- gene.length
		}
	
		splicedSeq <- ""
		for(j in 1:nrow(seq.index)){
			seq <- as.character(pcSeqs[[tids[i]]][seq.index$start[j]:seq.index$end[j]])
			splicedSeq <- paste(splicedSeq, seq, collapse = "", sep = "")
		}
		splicedSeqs[[i]] <- splicedSeq
		names(splicedSeqs)[i] <- tids[i]
		
		unsplicedATG <- temp.CDS$start[1]
		gap <- 0
		for(j in 1:nrow(temp.intron)){
			if(temp.intron$start[j] < unsplicedATG){
				gap <- gap + temp.intron$end[j] - temp.intron$start[j] + 1
			}
		}
		splicedATG <- unsplicedATG - gap
		splicedATGs[i] <- splicedATG
		
		unsplicedSTOP <- temp.CDS$end[nrow(temp.CDS)]
		gap <- 0
		for(j in 1:nrow(temp.intron)){
			if(temp.intron$end[j] < unsplicedSTOP){
				gap <- gap + temp.intron$end[j] - temp.intron$start[j] + 1
			}
		}
		splicedSTOP <- unsplicedSTOP - gap
		splicedSTOPs[i] <- splicedSTOP
	}

	if(nrow(temp.intron) == 0){
		splicedSeqs[[i]] <- as.character(pcSeqs[[tids[i]]])
		names(splicedSeqs)[i] <- tids[i]
		splicedATGs[i] <- temp.CDS$start[1]
		splicedSTOPs[i] <- temp.CDS$end[nrow(temp.CDS)]
	}
	if(i%%100 == 0) cat("i: ", i, "\n", sep = "")
}

CDS.ATGSTOP <- data.frame(transcript_id = tids, 
	ATG = splicedATGs, STOP = splicedSTOPs, length = width(splicedSeqs))

CDS.ATGSTOP$UTR5 <- FALSE
CDS.ATGSTOP$UTR3 <- FALSE
CDS.ATGSTOP$UTR5.start <- as.integer(NA)
CDS.ATGSTOP$UTR5.end <- as.integer(NA)
CDS.ATGSTOP$UTR3.start <- as.integer(NA)
CDS.ATGSTOP$UTR3.end <- as.integer(NA)
for(i in 1:nrow(CDS.ATGSTOP)){
	if(CDS.ATGSTOP$ATG[i] > 1){
		CDS.ATGSTOP$UTR5[i] <- TRUE
		CDS.ATGSTOP$UTR5.start[i] <- 1
		CDS.ATGSTOP$UTR5.end[i] <- CDS.ATGSTOP$ATG[i] - 1
	}
	if(CDS.ATGSTOP$STOP[i] < CDS.ATGSTOP$length[i]){
		CDS.ATGSTOP$UTR3[i] <- TRUE
		CDS.ATGSTOP$UTR3.start[i] <- CDS.ATGSTOP$STOP[i] + 1
		CDS.ATGSTOP$UTR3.end[i] <- CDS.ATGSTOP$length[i]
	}
}

setwd(HOMEDIR)
setwd("pombase")
dir.create("fasta")
dir.create("RData")
dir.create("txt")
save(pcSeqs, file = "RData/pcSeqs.RData")
save(splicedSeqs, file = "RData/splicedSeqs.RData")
save(CDS.ATGSTOP, file = "RData/CDS_ATGSTOP.RData")
save(gff3.pc, file = "RData/gff3_pc.RData")
save(gff3.mRNA, file = "RData/gff3_mRNA.RData")
writeXStringSet(pcSeqs, filepath = "fasta/pcSeqs.fasta")
writeXStringSet(splicedSeqs, filepath = "fasta/splicedSeqs.fasta")

write.table(CDS.ATGSTOP, file = "txt/CDS_ATGSTOP.txt", 
		col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
write.table(gff3.pc, file = "txt/gff3_pc.txt", 
		col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
write.table(gff3.mRNA, file = "txt/gff3_mRNA.txt", 
		col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

setwd(HOMEDIR)
setwd("pombase")
setwd("RData")
load(file = "splicedSeqs.RData")
load(file = "CDS_ATGSTOP.RData")

all.ORFs <- data.frame(transcript_id = as.character(NA), name = as.character(NA), 
			start = as.integer(NA), end = as.integer(NA), 
			length = as.integer(NA), phase = as.integer(NA), 
			uORF = FALSE, CDS = FALSE, CDSoverlap = FALSE, 
			stringsAsFactors = FALSE)
all.ORFs <- all.ORFs[0,]
for(i in 1:nrow(CDS.ATGSTOP)){
	transcript_id <- as.character(CDS.ATGSTOP$transcript_id[i])
	if(is.na(CDS.ATGSTOP$UTR5.start[i]) == FALSE){
		seq <- splicedSeqs[[transcript_id]]
		CDSATG <- CDS.ATGSTOP$ATG[i]
		CDSlength <- CDS.ATGSTOP$length[i]
		
		ORFs <- data.frame(transcript_id = as.character(NA), name = as.character(NA), 
					start = as.integer(NA), end = as.integer(NA), 
					length = as.integer(NA), phase = as.integer(NA), 
					stringsAsFactors = FALSE)
		ORFs <- ORFs[0,]
		for(j in 1:3){
			AA <- translate(subseq(seq, start = j), no.init.codon = TRUE)	# ATG のみ拾う
			ATGpos <- matchPattern(pattern = "M", subject = AA)@ranges@start
			STOPpos <- matchPattern(pattern = "*", subject = AA)@ranges@start
			ATGpos <- (ATGpos * 3) - 2 + j-1
			STOPpos <- (STOPpos * 3) + j-1
			ATGpos <- ATGpos[ATGpos <= CDSATG]
			if(length(ATGpos) > 0){
				ORFsph <- data.frame(transcript_id = as.character(NA), name = as.character(NA), 
						start = ATGpos, end = as.integer(NA), 
						length = as.integer(NA), phase = j, 
						stringsAsFactors = FALSE)
				for(k in 1:nrow(ORFsph)){
					ORFsph$end[k] <- STOPpos[STOPpos > ATGpos[k]][1]
					ORFsph$length[k] <- ORFsph$end[k] - ORFsph$start[k] + 1
				}
				ORFs <- rbind(ORFs, ORFsph)
			}
		}
		ORFs$uORF <- FALSE
		ORFs$CDS <- FALSE
		ORFs$CDSoverlap <- FALSE
		for(j in 1:nrow(ORFs)){
			if(ORFs$start[j] == CDSATG) ORFs$CDS[j] <- TRUE
			if(ORFs$start[j] != CDSATG){
				ORFs$uORF[j] <- TRUE
				if(ORFs$end[j] > CDSATG) ORFs$CDSoverlap[j] <- TRUE
			}
		}
		ORFs <- ORFs[order(ORFs$start, decreasing = FALSE),]
		ORFs$transcript_id <- transcript_id
		for(j in 1:nrow(ORFs)){
			if(ORFs$CDS[j] == FALSE){
				ORFs$name[j] <- paste("uORF", j, sep = "")
			}
			if(ORFs$CDS[j] == TRUE){
				ORFs$name[j] <- substr(transcript_id, start = 1, stop = nchar(transcript_id) - 2)
				ORFs$CDSoverlap[j] <- as.logical(NA)
			}
		}
		all.ORFs <- rbind(all.ORFs, ORFs)
	}

	if(is.na(CDS.ATGSTOP$UTR5.start[i]) == TRUE){
		ORFs <- data.frame(transcript_id = transcript_id, 
					name = substr(transcript_id, start = 1, stop = nchar(transcript_id) - 2), 
					start = CDS.ATGSTOP$ATG[i], 
					end = CDS.ATGSTOP$STOP[i], 
					length = CDS.ATGSTOP$length[i], 
					phase = 1, 
					stringsAsFactors = FALSE)
		ORFs$uORF <- FALSE
		ORFs$CDS <- TRUE
		ORFs$CDSoverlap <- as.logical(NA)
		all.ORFs <- rbind(all.ORFs, ORFs)
	}
	if(i%%100 == 0) cat("i:", i, "\n", sep = "")
}

all.ORFs$uORFseq <- as.character(NA)
for(i in 1:nrow(all.ORFs)){
	if(all.ORFs$uORF[i] == TRUE){
		transcript_id <- all.ORFs$transcript_id[i]
		start <- all.ORFs$start[i]
		end <- all.ORFs$end[i]
		all.ORFs$uORFseq[i] <- as.character(splicedSeqs[[transcript_id]][start:end])
	}
	if(i%%100 == 0) cat("i:", i, "\n", sep = "")
}

setwd(HOMEDIR)
setwd("pombase")
save(all.ORFs, file = "RData/all_ORFs.RData")
write.table(all.ORFs, file = "txt/all_ORFs.txt", 
		col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

setwd(HOMEDIR)
setwd("pombase")
load(file = "RData/all_ORFs.RData")

for(i in 1:nrow(all.ORFs)){
	tid <- all.ORFs$transcript_id[i]
	name <- all.ORFs$name[i]
	gene_id <- substr(tid, start = 1, stop = nchar(tid)-2)
	if(name == gene_id) all.ORFs$name[i] <- paste(gene_id, "_CDS", sep = "")
	if(name != gene_id) all.ORFs$name[i] <- paste(gene_id, "_", name, sep = "")
}

dir.create("xlsx")
library("openxlsx")
write.xlsx(x = all.ORFs, file = "xlsx/all_ORFs.xlsx")
