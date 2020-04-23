# Specify HOMEDIR before run.
HOMEDIR <- "XXX"


library(Biostrings)

setwd(HOMEDIR)
setwd("pombase")
gff3 <- read.table(file = "Schizosaccharomyces_pombe_all_chromosomes.gff3.gz", 
			sep = "\t", skip = 1, stringsAsFactors = FALSE)
names(gff3) <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")

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

gff3.gene <- subset(gff3, type == "gene")
gff3.gene$gene_id <- "NA"
for(i in 1:nrow(gff3.gene)){
	gff3.gene$gene_id[i] <- get.ID.attrb1(gff3.gene$attributes[i])
}
gff3.gene$transcript_id <- "NA"
for(i in 1:nrow(gff3.gene)){
	gff3.gene$transcript_id[i] <- paste(gff3.gene$gene_id[i], ".1", sep = "")
}

nrow(gff3.gene)
length(unique(gff3.gene$gene_id))
head(gff3.gene)

gff3.intron <- subset(gff3, type == "intron")
gff3.intron$gene_id <- "NA"
for(i in 1:nrow(gff3.intron)){
	gff3.intron$gene_id[i] <- get.ID.attrb2(gff3.intron$attributes[i])
}
gff3.intron$transcript_id <- "NA"
for(i in 1:nrow(gff3.intron)){
	gff3.intron$transcript_id[i] <- paste(gff3.intron$gene_id[i], ".1", sep = "")
}

nrow(gff3.intron)
length(unique(gff3.intron$gene_id))
head(gff3.intron)

setwd(HOMEDIR)
setwd("pombase")
sp.genome <- readDNAStringSet(filepath = "Schizosaccharomyces_pombe_all_chromosomes.fa.gz")
chrNames <- names(sp.genome)
for(i in 1:6){
	chrNames[i] <- strsplit(chrNames[i], split = " ")[[1]][1]
}
names(sp.genome) <- chrNames

splicedSeqs2 <- DNAStringSet()
for(i in 1:nrow(gff3.gene)){
	tid <- gff3.gene$transcript_id[i]
	temp <- subset(gff3.gene, transcript_id == tid)
	chromosome <- temp$seqid[1]
	gene.start <- temp$start[1]
	gene.end <- temp$end[1]
	gene.strand <- temp$strand[1]
	gene.length <- gene.end - gene.start + 1

	temp.intron <- 	subset(gff3.intron, transcript_id == tid)
	temp.intron <- temp.intron[order(temp.intron$start, decreasing = FALSE),]
	
	if(nrow(temp.intron) != 0){
		seq.index <- data.frame(start = integer(length = nrow(temp.intron)+1), 
				end = integer(length = nrow(temp.intron)+1))
		for(j in 1:(nrow(temp.intron)+1)){
			if(j == 1) seq.index$start[j] <- gene.start
			if(j == 1) seq.index$end[j] <- temp.intron$start[j] - 1
			if(j != 1) seq.index$start[j] <- temp.intron$end[j-1] +1
			if(j != 1) seq.index$end[j] <- temp.intron$start[j] - 1
			if(j == (nrow(temp.intron)+1)) seq.index$start[j] <- temp.intron$end[j-1] +1
			if(j == (nrow(temp.intron)+1)) seq.index$end[j] <- gene.end
		}
	
		splicedSeq <- ""
		for(j in 1:nrow(seq.index)){
			seq <- as.character(sp.genome[[chromosome]][seq.index$start[j]:seq.index$end[j]])
			splicedSeq <- paste(splicedSeq, seq, collapse = "", sep = "")
		}
		if(gene.strand == "+")	splicedSeqs2[[i]] <- splicedSeq
		if(gene.strand == "-")	splicedSeqs2[[i]] <- as.character(reverseComplement(DNAString(splicedSeq)))
		
		names(splicedSeqs2)[i] <- tid
	}

	if(nrow(temp.intron) == 0){
		unsplicedSeq <- sp.genome[[chromosome]][gene.start:gene.end]
		if(gene.strand == "-")  unsplicedSeq <- reverseComplement(unsplicedSeq)
			
		splicedSeqs2[[i]] <- as.character(unsplicedSeq)
		names(splicedSeqs2)[i] <- tid
	}
	if(i%%100 == 0) cat("i: ", i, "\n", sep = "")
}

setwd(HOMEDIR)
setwd("pombase")
save(gff3.gene, file = "RData/gff3_gene.RData")
save(gff3.intron, file = "RData/gff3_intron.RData")
save(sp.genome, file = "RData/sp_genome.RData")
save(splicedSeqs2, file = "RData/splicedSeqs2.RData")
writeXStringSet(splicedSeqs2, filepath = "fasta/splicedSeqs2.fasta")

setwd(HOMEDIR)
dir.create("reference_transcriptome")
setwd("reference_transcriptome")
filename <- "sp_spliced_transcripts.fasta"
writeXStringSet(splicedSeqs2, filepath = filename)

out <- data.frame(name = names(splicedSeqs2), width = width(splicedSeqs2), 
			stringsAsFactors = FALSE)
setwd(HOMEDIR)
setwd("genome_info")
write.table(out, file = "sp_spliced_transcriptome_info.txt", quote = FALSE, 
	col.names = FALSE, row.names = FALSE, sep = "\t")


setwd(HOMEDIR)
setwd("pombase/RData")
load(file = "splicedSeqs2.RData")

targetWidth <- 20
allLines <- 0
for(i in 1:length(splicedSeqs2)){
	seqWidth <- width(splicedSeqs2[i])
	lines <- seqWidth - targetWidth + 1
	allLines <- allLines + lines
}

seqDB <- character(length = allLines)

allLines <- 0
for(i in 1:length(splicedSeqs2)){
	seqFull <- splicedSeqs2[[i]]
	seqWidth <- width(splicedSeqs2[i])
	lines <- seqWidth - targetWidth + 1
	for(j in 1:lines){
		allLines <- allLines + 1
		SEQ <- subseq(seqFull, start = j, end = (j + targetWidth - 1))
		seqDB[allLines] <- as.character(SEQ)
	}
	if(i%%20 == 0) cat(i, ",", sep = "")
	if(i%%200 == 0) cat("/ total ", length(splicedSeqs2), "\n", sep = "")
}

duplicatedSeqs <- unique(seqDB[duplicated(seqDB)])

dupTRUE <- duplicated(seqDB)

dupTRUEfromLast <- duplicated(seqDB, fromLast = TRUE)

dupTRUEboth <- dupTRUE + dupTRUEfromLast
dupTRUEboth[dupTRUEboth > 1] <- 1

dupTRUEboth <- as.logical(dupTRUEboth)

splicedSeqsM <- splicedSeqs2

allLines <- 0
for(i in 1:length(splicedSeqsM)){
	seqFull <- splicedSeqsM[[i]]
	seqWidth <- width(splicedSeqsM[i])
	start <- allLines + 1
	lines <- seqWidth - targetWidth + 1
	end <- allLines + lines

	dupTRUEseq <- dupTRUEboth[start:end]
	
	# 20200311
	Nstart <- which(dupTRUEseq == TRUE)
	if(length(Nstart) > 0){
		for(n in 1:length(Nstart)){
			N20 <- Nstart[n]:(Nstart[n]+19)
			splicedSeqsM[[i]][N20] <- "N"
		}
	}
	
	allLines <- end
	if(i%%20 == 0) cat(i, ",", sep = "")
	if(i%%200 == 0) cat("/ total ", length(splicedSeqsM), "\n", sep = "")
}

setwd(HOMEDIR)
setwd("pombase")
save(splicedSeqsM, file = "RData/splicedSeqsM2.RData")		
writeXStringSet(splicedSeqsM, filepath = "fasta/splicedSeqsM2.fasta")

