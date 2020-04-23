## The output files are used in the indicated figures.
# output: ratioWith3AT/ratioWith3AT_RPKMthresh6_bw.pdf
#         used in S7A,C1-4
# output: ratioWith3AT/Length_RPKM6_WIDTH27_v11/ratioWith3AT_RPKMthresh6_len_bar.pdf
#         used in S7A, C5
# output: ratioWo3AT/ratioWo3AT_RPKMthresh6_bw.pdf
#         used in S7B,D1-4
# output: ratioWo3AT/Length_RPKM6_WIDTH27_v11/ratioWo3AT_RPKMthresh6_len_bar.pdf
#         used in S7B, D5


# Specify HOMEDIR before run.
HOMEDIR <- "XXX"


setwd(HOMEDIR)
setwd("pombase")
load(file = "RData/splicedSeqs2.RData")
load(file = "RData/splicedSeqsM2.RData")

setwd(HOMEDIR)
setwd("RPKM_2")
setwd("RData")
load(file = "bundledDF_20200312.RData")

analyzeRatiouORFv11 <- function(inDF, type = "ratioWT", RPKM.threshold = 6, WIDTH = 27){

	wd3 <- getwd()
	
	OUTDIR1 <- paste("RPKM", RPKM.threshold, "_WIDTH", WIDTH, sep = "")
	dir.create(OUTDIR1)
	setwd(OUTDIR1)
	
	OUTDIR2 <- type
	dir.create(OUTDIR2)
	setwd(OUTDIR2)
	wd <- getwd()

	SAMPLE1 <- "WTwo3AT"
	SAMPLE2 <- "WTwith3AT"
	SAMPLE3 <- "gcn2wo3AT"
	SAMPLE4 <- "gcn2with3AT"

	LOG1 <- character(length = 12)
	LOG1[1] <- paste("type: ", type, sep = "")
	temp <- paste("SAMPLE1: ", SAMPLE1, ", SAMPLE2: ", SAMPLE2, sep = "")
	LOG1[2] <- paste(temp, ", SAMPLE3: ", SAMPLE3, ", SAMPLE4: ", SAMPLE4, sep = "")
	LOG1[3] <- ""

	LOG1[4] <- paste("input: ", nrow(inDF), sep = "")
	LOG1[5] <- paste(" CDS: ", nrow(subset(inDF, CDS == TRUE)), sep = "")

	LOG1[6] <- paste(" CDS masked: ", nrow(subset(inDF, CDS == TRUE & masked4 == TRUE)), sep = "")
	LOG1[7] <- paste(" CDS notMasked: ", nrow(subset(inDF, CDS == TRUE & masked4 == FALSE)), sep = "")
	LOG1[8] <- paste(" uORF: ", nrow(subset(inDF, uORF == TRUE)), sep = "")
	LOG1[9] <- paste(" uORF masked: ", nrow(subset(inDF, uORF == TRUE & masked == TRUE)), sep = "")
	LOG1[10] <- paste(" uORF notMasked: ", nrow(subset(inDF, uORF == TRUE & masked == FALSE)), sep = "")

	LOG1[11] <- ""
	LOG1[12] <- paste(" uORF, CDSoverlap==TRUE: ", 
			nrow(subset(inDF, uORF == TRUE & CDSoverlap == TRUE)), sep = "")
	LOG1[13] <- paste(" uORF, CDSoverlap==FALSE: ", 
			nrow(subset(inDF, uORF == TRUE & CDSoverlap == FALSE)), sep = "")

	LOG1[14] <- paste(" uORF and CDS notMasked: ", 
			nrow(subset(inDF, uORF == TRUE & masked == FALSE & masked4 == FALSE)), sep = "")
	LOG1[15] <- paste(" uORF and CDS notMasked, CDSoverlap==TRUE: ", 
			nrow(subset(inDF, uORF == TRUE & masked == FALSE & masked4 == FALSE & CDSoverlap == TRUE)), 
			sep = "")
	LOG1[16] <- paste(" uORF and CDS notMasked, CDSoverlap==FALSE: ", 
			nrow(subset(inDF, uORF == TRUE & masked == FALSE & masked4 == FALSE & CDSoverlap == FALSE)), 
			sep = "")

	LOG1[17] <- ""
	LOG1[18] <- paste("Stepwise filtering (RPKM.threshold=", RPKM.threshold, ")", sep = "")
	
	checkuORFs <- function(inputDF, char = ""){
		uORFnames <- c("SPCC1393.08.1_uORF1", "SPCC1393.08.1_uORF4", 
			"SPAC1952.05.1_uORF3", "SPAC20G4.03c.1_uORF2", "SPAC222.07c.1_uORF3")
		out <- inputDF[0,]
		for(i in 1:length(uORFnames)){
			index <- grep(pattern = uORFnames[i], x = inputDF$name)
				out <- rbind(out, inputDF[index,])
		}
		filename <- paste("checkuORFs_", char, "_RPKMth", RPKM.threshold, ".xlsx", sep = "")
		write.xlsx(out, filename)
	}
	checkuORFs(inDF, char = "1_inDF")

	
	selectedDF <- inDF
	LOG1[19] <- paste(" Before filtering: ", nrow(selectedDF), sep = "")
	selectedDF <- subset(selectedDF, uORF == TRUE)
	LOG1[20] <- paste(" uORF==TRUE: ", nrow(selectedDF), sep = "")
	selectedDF <- subset(selectedDF, CDSoverlap == FALSE)
	LOG1[21] <- paste(" CDSoverlap==FALSE: ", nrow(selectedDF), sep = "")

	selectedDF <- subset(selectedDF, masked == FALSE & masked4 == FALSE)
	checkuORFs(selectedDF, char = "4_maskedRemoved")
	LOG1[22] <- paste(" unMasked: ", nrow(selectedDF), sep = "")

	selectedDF <- selectedDF[selectedDF[[SAMPLE1]] > RPKM.threshold,]
	checkuORFs(selectedDF, char = "5_SAMPLE1_threshold")
	LOG1[23] <- paste(" ", SAMPLE1, " more than RPKM.threshold (", RPKM.threshold, "):", nrow(selectedDF), sep = "")
	selectedDF <- selectedDF[selectedDF[[SAMPLE2]] > RPKM.threshold,]
	checkuORFs(selectedDF, char = "6_SAMPLE2_threshold")
	LOG1[24] <- paste(" ", SAMPLE2, " more than RPKM.threshold (", RPKM.threshold, "):", nrow(selectedDF), sep = "")
	selectedDF <- selectedDF[selectedDF[[SAMPLE3]] > RPKM.threshold,]
	checkuORFs(selectedDF, char = "7_SAMPLE3_threshold")
	LOG1[25] <- paste(" ", SAMPLE3, " more than RPKM.threshold (", RPKM.threshold, "):", nrow(selectedDF), sep = "")
	selectedDF <- selectedDF[selectedDF[[SAMPLE4]] > RPKM.threshold,]
	checkuORFs(selectedDF, char = "8_SAMPLE4_threshold")
	LOG1[26] <- paste(" ", SAMPLE4, " more than RPKM.threshold (", RPKM.threshold, "):", nrow(selectedDF), sep = "")

	LOG1[27] <- ""
	LOG1[28] <- "Frame filtering"

	FOLD <- 1.5
	colNamePhase2_1 <- paste("phase2", SAMPLE1, sep = "")
	selectedDF <- selectedDF[selectedDF[[SAMPLE1]] > selectedDF[[colNamePhase2_1]]*FOLD,]
	colNamePhase2_2 <- paste("phase2", SAMPLE2, sep = "")
	selectedDF <- selectedDF[selectedDF[[SAMPLE2]] > selectedDF[[colNamePhase2_2]]*FOLD,]
	colNamePhase2_3 <- paste("phase2", SAMPLE3, sep = "")
	selectedDF <- selectedDF[selectedDF[[SAMPLE3]] > selectedDF[[colNamePhase2_3]]*FOLD,]
	colNamePhase2_4 <- paste("phase2", SAMPLE4, sep = "")
	selectedDF <- selectedDF[selectedDF[[SAMPLE4]] > selectedDF[[colNamePhase2_4]]*FOLD,]

	colNamePhase3_1 <- paste("phase3", SAMPLE1, sep = "")
	selectedDF <- selectedDF[selectedDF[[SAMPLE1]] > selectedDF[[colNamePhase3_1]]*FOLD,]
	colNamePhase3_2 <- paste("phase3", SAMPLE2, sep = "")
	selectedDF <- selectedDF[selectedDF[[SAMPLE2]] > selectedDF[[colNamePhase3_2]]*FOLD,]
	colNamePhase3_3 <- paste("phase3", SAMPLE3, sep = "")
	selectedDF <- selectedDF[selectedDF[[SAMPLE3]] > selectedDF[[colNamePhase3_3]]*FOLD,]
	colNamePhase3_4 <- paste("phase3", SAMPLE4, sep = "")
	selectedDF <- selectedDF[selectedDF[[SAMPLE4]] > selectedDF[[colNamePhase3_4]]*FOLD,]
	checkuORFs(selectedDF, char = "9_phase23")
	# LOG1[29] <- paste(" phase1 RPKM is twice higher than those of phases 2 and 3: ", nrow(selectedDF), sep = "")
	LOG1[29] <- paste(" phase1 RPKM is ", FOLD, " fold-higher than those of phases 2 and 3: ", 
			nrow(selectedDF), sep = "")
	
	temp1 <- list(WT = log2(selectedDF$ratioWT), gcn2 = log2(selectedDF$ratiogcn2))
	ylab1 <- "(uORF/CDS with 3AT) / (uORF/CDS wo 3AT)"
	main1 <- paste("ratioWT vs ratiogcn2", ", RPKM=", RPKM.threshold, sep = "")
	pvalue1 <- t.test(temp1[[1]], temp1[[2]])$p.value	# Welch Two Sample t-test
	sub1 <- paste("n=", nrow(selectedDF), ", p=", formatC(pvalue1, format = "e", digits = 2), sep = "")

	temp2 <- list(wo3AT = log2(selectedDF$ratioWo3AT), with3AT = log2(selectedDF$ratioWith3AT))
	ylab2 <- "(uORF/CDS in gcn2) / (uORF/CDS in WT)"
	main2 <- paste("ratioWo3AT vs ratioWith3AT", ", RPKM=", RPKM.threshold, sep = "")
	pvalue2 <- t.test(temp2[[1]], temp2[[2]])$p.value	# Welch Two Sample t-test
	sub2 <- paste("n=", nrow(selectedDF), ", p=", formatC(pvalue2, format = "e", digits = 2), sep = "")

	ylim1 = c(min(unlist(temp1), na.rm = TRUE), max(unlist(temp1), na.rm = TRUE))
	ylim2 = c(min(unlist(temp2), na.rm = TRUE), max(unlist(temp2), na.rm = TRUE))

	filename <- paste(type, "_RPKMthresh", RPKM.threshold, "_bw.pdf", sep = "")
	pdf(file = filename, width = 8*5, height = 8)
	par(mfrow = c(1, 2*5))
		boxplot(temp1, outline = FALSE, main = main1, sub = sub1, cex.axis = 1, ylab = ylab1, ylim = ylim1)
		abline(h = 0, col = "gray")
		beeswarm(temp1, corral = "wrap", col = rgb(1, 0, 0, alpha = 0.5), pch = 20, add = TRUE)

		boxplot(temp2, outline = FALSE, main = main2, sub = sub2, cex.axis = 1, ylab = ylab2, ylim = ylim2)
		abline(h = 0, col = "gray")
		beeswarm(temp2, corral = "wrap", col = rgb(1, 0, 0, alpha = 0.5), pch = 20, add = TRUE)

	LOG1[30] <- ""
	LOG1[31] <- main1
	LOG1[32] <- sub1
	LOG1[33] <- main2
	LOG1[34] <- sub2

	colName <- "unmasked_start"
	selectedDF[[colName]] <- as.character(NA)
	colName <- "masked_start"
	selectedDF[[colName]] <- as.character(NA)
	
	colName <- "unmasked_end"
	selectedDF[[colName]] <- as.character(NA)
	colName <- "masked_end"
	selectedDF[[colName]] <- as.character(NA)
	
	colName <- "is_masked_start"
	selectedDF[[colName]] <- FALSE
	
	colName <- "is_masked_end"
	selectedDF[[colName]] <- FALSE
	
	for(i in 1:nrow(selectedDF)){
		tid <- selectedDF$transcript_id[i]
		start <- selectedDF$start[i]
		end <- selectedDF$end[i]
		if(start > WIDTH & (WIDTH+end) < (length(splicedSeqs2[[tid]]) -1) ){
			SeqStart <- as.character(splicedSeqs2[[tid]][(start-WIDTH):(start+WIDTH+2)])
			SeqStartM <- as.character(splicedSeqsM[[tid]][(start-WIDTH):(start+WIDTH+2)])
			SeqEnd <- as.character(splicedSeqs2[[tid]][(end-2-WIDTH):(end+WIDTH)])
			SeqEndM <- as.character(splicedSeqsM[[tid]][(end-2-WIDTH):(end+WIDTH)])
		}else{
			SeqStart <- as.character(NA)
			SeqStartM <- as.character(NA)
			SeqEnd <- as.character(NA)
			SeqEndM <- as.character(NA)
		}
		colName1 <- "unmasked_start"
		selectedDF[[colName1]][i] <- SeqStart
		colName2 <- "masked_start"
		selectedDF[[colName2]][i] <- SeqStartM
		colName3 <- "unmasked_end"
		selectedDF[[colName3]][i] <- SeqEnd
		colName4 <- "masked_end"
		selectedDF[[colName4]][i] <- SeqEndM
		colName5 <- "is_masked_start"
		colName6 <- "is_masked_end"
		if(is.na(SeqStart) == FALSE & is.na(SeqStartM) == FALSE & 
			is.na(SeqEnd) == FALSE & is.na(SeqEndM) == FALSE){
			if(SeqStart != SeqStartM){
				selectedDF[[colName5]][i] <- TRUE	
			}
			if(SeqEnd != SeqEndM){
				selectedDF[[colName6]][i] <- TRUE	
			}
		}
		if(is.na(SeqStart) == TRUE | is.na(SeqStartM) == TRUE){
			selectedDF[[colName5]][i] <- as.logical(NA)
		}
		if(is.na(SeqEnd) == TRUE | is.na(SeqEndM) == TRUE){
			selectedDF[[colName6]][i] <- as.logical(NA)
		}
	}

	LOG1[35] <- ""
	LOG1[36] <- "Sequence"
	LOG1[37] <- paste("LENGTH: ", "not applicable", ", ", "WIDTH: ", WIDTH, sep = "")
	LOG1[38] <- paste("nrow(selectedDF): ", nrow(selectedDF), sep = "")
	LOG1[39] <- paste("is_masked_start == FALS: ", nrow(subset(selectedDF, is_masked_start == FALSE)), sep = "")
	LOG1[40] <- paste("is_masked_start == TRUE: ", nrow(subset(selectedDF, is_masked_start == TRUE)), sep = "")
	LOG1[41] <- paste("is.na(is_masked_start) == TRUE): ", 
			nrow(subset(selectedDF, is.na(is_masked_start) == TRUE)), sep = "")
	LOG1[42] <- paste("is_masked_end == FALSE): ", nrow(subset(selectedDF, is_masked_end == FALSE)), sep = "")
	LOG1[43] <- paste("is_masked_end == TRUE: ", nrow(subset(selectedDF, is_masked_end == TRUE)), sep = "")
	LOG1[44] <- paste("is.na(is_masked_end) == TRUE: ", 
			nrow(subset(selectedDF, is.na(is_masked_end) == TRUE)), sep = "")

	qvalues <- quantile(selectedDF[[type]], na.rm = TRUE)
	TYPE.Q1 <- subset(selectedDF, selectedDF[[type]] > qvalues[4])
	TYPE.Q2 <- subset(selectedDF, selectedDF[[type]] > qvalues[3] & selectedDF[[type]] <= qvalues[4])
	TYPE.Q3 <- subset(selectedDF, selectedDF[[type]] > qvalues[2] & selectedDF[[type]] <= qvalues[3])
	TYPE.Q4 <- subset(selectedDF, selectedDF[[type]] >= qvalues[1] & selectedDF[[type]] <= qvalues[2])

	temp3 <- list(Q1 = log2(TYPE.Q1$ratioWT), Q2 = log2(TYPE.Q2$ratioWT), 
			Q3 = log2(TYPE.Q3$ratioWT), Q4 = log2(TYPE.Q4$ratioWT))
	temp4 <- list(Q1 = log2(TYPE.Q1$ratiogcn2), Q2 = log2(TYPE.Q2$ratiogcn2), 
			Q3 = log2(TYPE.Q3$ratiogcn2), Q4 = log2(TYPE.Q4$ratiogcn2))
	ylim34 <- c(min(c(unlist(temp3), unlist(temp4)), na.rm = TRUE), 
			max(c(unlist(temp3), unlist(temp4)), na.rm = TRUE))
	main3 <- "ratioWT in groups Q1-Q4"
	main4 <- "ratiogcn2 in groups Q1-Q4"
	ylab34 <- "(uORF/CDS with 3AT) / (uORF/CDS wo 3AT)"

	pQ1 <- formatC(t.test(temp3[[1]], temp4[[1]])$p.value, format = "e", digits = 2)
	pQ2 <- formatC(t.test(temp3[[2]], temp4[[2]])$p.value, format = "e", digits = 2)
	pQ3 <- formatC(t.test(temp3[[3]], temp4[[3]])$p.value, format = "e", digits = 2)
	pQ4 <- formatC(t.test(temp3[[4]], temp4[[4]])$p.value, format = "e", digits = 2)
	nQ1 <- length(temp3[[1]])
	nQ2 <- length(temp3[[2]])
	nQ3 <- length(temp3[[3]])
	nQ4 <- length(temp3[[4]])
	sub34 <- paste("n=", nQ1, ",", nQ2, ",", nQ3, ",", nQ4, sep = "")
	pvalue34 <- as.numeric(c(pQ1, pQ2, pQ3, pQ4))

	boxplot(temp3, outline = FALSE, main = main3, sub = sub34, cex.axis = 1, ylab = ylab34, ylim = ylim34)
	abline(h = 0, col = "gray")
	beeswarm(temp3, corral = "wrap", col = rgb(1, 0, 0, alpha = 0.5), pch = 20, add = TRUE)

	boxplot(temp4, outline = FALSE, main = main4, sub = sub34, cex.axis = 1, ylab = ylab34, ylim = ylim34)
	y <- c(max(temp4[[1]]), max(temp4[[2]]), max(temp4[[3]]), max(temp4[[4]])) + 0.1
	abline(h = 0, col = "gray")
	beeswarm(temp4, corral = "wrap", col = rgb(1, 0, 0, alpha = 0.5), pch = 20, add = TRUE)

	LOG1[45] <- ""
	LOG1[46] <- main3
	LOG1[47] <- main4
	LOG1[48] <- sub34
	LOG1[49] <- paste("p=", pQ1, pQ2, pQ3, pQ4, sep = ",")

	temp5 <- list(Q1 = log2(TYPE.Q1$WTwo3ATvsCDS), Q2 = log2(TYPE.Q2$WTwo3ATvsCDS), 
			Q3 = log2(TYPE.Q3$WTwo3ATvsCDS), Q4 = log2(TYPE.Q4$WTwo3ATvsCDS))
	temp6 <- list(Q1 = log2(TYPE.Q1$WTwith3ATvsCDS), Q2 = log2(TYPE.Q2$WTwith3ATvsCDS), 
			Q3 = log2(TYPE.Q3$WTwith3ATvsCDS), Q4 = log2(TYPE.Q4$WTwith3ATvsCDS))
	ylim56 <- c(min(c(unlist(temp5), unlist(temp6)), na.rm = TRUE), 
			max(c(unlist(temp5), unlist(temp6)), na.rm = TRUE))
	main5 <- "WTwo3ATvsCDS in groups Q1-Q4"
	main6 <- "WTwith3ATvsCDS in groups Q1-Q4"
	ylab5 <- "uORF/CDS wo 3AT"
	ylab6 <- "uORF/CDS with 3AT"

	pQ1 <- formatC(t.test(temp5[[1]], temp6[[1]])$p.value, format = "e", digits = 2)
	pQ2 <- formatC(t.test(temp5[[2]], temp6[[2]])$p.value, format = "e", digits = 2)
	pQ3 <- formatC(t.test(temp5[[3]], temp6[[3]])$p.value, format = "e", digits = 2)
	pQ4 <- formatC(t.test(temp5[[4]], temp6[[4]])$p.value, format = "e", digits = 2)
	nQ1 <- length(temp5[[1]])
	nQ2 <- length(temp5[[2]])
	nQ3 <- length(temp5[[3]])
	nQ4 <- length(temp5[[4]])
	sub56 <- paste("n=", nQ1, ",", nQ2, ",", nQ3, ",", nQ4, sep = "")
	pvalue56 <- as.numeric(c(pQ1, pQ2, pQ3, pQ4))

	boxplot(temp5, outline = FALSE, main = main5, sub = sub56, cex.axis = 1, ylab = ylab5, ylim = ylim56)
	abline(h = 0, col = "gray")
	beeswarm(temp5, corral = "wrap", col = rgb(1, 0, 0, alpha = 0.5), pch = 20, add = TRUE)

	boxplot(temp6, outline = FALSE, main = main6, sub = sub56, cex.axis = 1, ylab = ylab6, ylim = ylim56)
	y <- c(max(temp6[[1]]), max(temp6[[2]]), max(temp6[[3]]), max(temp6[[4]])) + 0.1
	abline(h = 0, col = "gray")
	beeswarm(temp6, corral = "wrap", col = rgb(1, 0, 0, alpha = 0.5), pch = 20, add = TRUE)

	LOG1[50] <- ""
	LOG1[51] <- main5
	LOG1[52] <- main6
	LOG1[53] <- sub56
	LOG1[54] <- paste("p=", pQ1, pQ2, pQ3, pQ4, sep = ",")

	temp7 <- list(Q1 = log2(TYPE.Q1$gcn2wo3ATvsCDS), Q2 = log2(TYPE.Q2$gcn2wo3ATvsCDS), 
			Q3 = log2(TYPE.Q3$gcn2wo3ATvsCDS), Q4 = log2(TYPE.Q4$gcn2wo3ATvsCDS))
	temp8 <- list(Q1 = log2(TYPE.Q1$gcn2with3ATvsCDS), Q2 = log2(TYPE.Q2$gcn2with3ATvsCDS), 
			Q3 = log2(TYPE.Q3$gcn2with3ATvsCDS), Q4 = log2(TYPE.Q4$gcn2with3ATvsCDS))
	ylim78 <- c(min(c(unlist(temp7), unlist(temp8)), na.rm = TRUE), 
			max(c(unlist(temp7), unlist(temp8)), na.rm = TRUE))
	main5 <- "gcn2wo3ATvsCDS in groups Q1-Q4"
	main6 <- "gcn2with3ATvsCDS in groups Q1-Q4"
	ylab7 <- "uORF/CDS wo 3AT"
	ylab8 <- "uORF/CDS with 3AT"

	pQ1 <- formatC(t.test(temp7[[1]], temp8[[1]])$p.value, format = "e", digits = 2)
	pQ2 <- formatC(t.test(temp7[[2]], temp8[[2]])$p.value, format = "e", digits = 2)
	pQ3 <- formatC(t.test(temp7[[3]], temp8[[3]])$p.value, format = "e", digits = 2)
	pQ4 <- formatC(t.test(temp7[[4]], temp8[[4]])$p.value, format = "e", digits = 2)
	nQ1 <- length(temp7[[1]])
	nQ2 <- length(temp7[[2]])
	nQ3 <- length(temp7[[3]])
	nQ4 <- length(temp7[[4]])
	sub78 <- paste("n=", nQ1, ",", nQ2, ",", nQ3, ",", nQ4, sep = "")
	pvalue78 <- as.numeric(c(pQ1, pQ2, pQ3, pQ4))

	boxplot(temp7, outline = FALSE, main = main5, sub = sub78, cex.axis = 1, ylab = ylab7, ylim = ylim78)
	abline(h = 0, col = "gray")
	beeswarm(temp7, corral = "wrap", col = rgb(1, 0, 0, alpha = 0.5), pch = 20, add = TRUE)

	boxplot(temp8, outline = FALSE, main = main6, sub = sub78, cex.axis = 1, ylab = ylab8, ylim = ylim78)
	y <- c(max(temp8[[1]]), max(temp8[[2]]), max(temp8[[3]]), max(temp8[[4]])) + 0.1
	abline(h = 0, col = "gray")
	beeswarm(temp8, corral = "wrap", col = rgb(1, 0, 0, alpha = 0.5), pch = 20, add = TRUE)

	LOG1[55] <- ""
	LOG1[56] <- main5
	LOG1[57] <- main6
	LOG1[58] <- sub78
	LOG1[59] <- paste("p=", pQ1, pQ2, pQ3, pQ4, sep = ",")

	temp9 <- list(Q1 = log2(TYPE.Q1$WTwo3ATvsCDS), Q2 = log2(TYPE.Q2$WTwo3ATvsCDS), 
			Q3 = log2(TYPE.Q3$WTwo3ATvsCDS), Q4 = log2(TYPE.Q4$WTwo3ATvsCDS))
	temp10 <- list(Q1 = log2(TYPE.Q1$gcn2wo3ATvsCDS), Q2 = log2(TYPE.Q2$gcn2wo3ATvsCDS), 
			Q3 = log2(TYPE.Q3$gcn2wo3ATvsCDS), Q4 = log2(TYPE.Q4$gcn2wo3ATvsCDS))
	ylim910 <- c(min(c(unlist(temp9), unlist(temp10)), na.rm = TRUE), 
			max(c(unlist(temp9), unlist(temp10)), na.rm = TRUE))
	main9 <- "WTwo3ATvsCDS in groups Q1-Q4"
	main10 <- "gcn2wo3ATvsCDS in groups Q1-Q4"
	ylab9 <- "uORF/CDS wo 3AT"
	ylab10 <- "uORF/CDS wo 3AT"

	pQ1 <- formatC(t.test(temp9[[1]], temp10[[1]])$p.value, format = "e", digits = 2)
	pQ2 <- formatC(t.test(temp9[[2]], temp10[[2]])$p.value, format = "e", digits = 2)
	pQ3 <- formatC(t.test(temp9[[3]], temp10[[3]])$p.value, format = "e", digits = 2)
	pQ4 <- formatC(t.test(temp9[[4]], temp10[[4]])$p.value, format = "e", digits = 2)
	nQ1 <- length(temp9[[1]])
	nQ2 <- length(temp9[[2]])
	nQ3 <- length(temp9[[3]])
	nQ4 <- length(temp9[[4]])
	sub910 <- paste("n=", nQ1, ",", nQ2, ",", nQ3, ",", nQ4, sep = "")
	pvalue910 <- as.numeric(c(pQ1, pQ2, pQ3, pQ4))

	LOG1[60] <- ""
	LOG1[61] <- main9
	LOG1[62] <- main10
	LOG1[63] <- sub910
	LOG1[64] <- paste("p=", pQ1, pQ2, pQ3, pQ4, sep = ",")

	temp11 <- list(Q1 = log2(TYPE.Q1$WTwith3ATvsCDS), Q2 = log2(TYPE.Q2$WTwith3ATvsCDS), 
			Q3 = log2(TYPE.Q3$WTwith3ATvsCDS), Q4 = log2(TYPE.Q4$WTwith3ATvsCDS))
	temp12 <- list(Q1 = log2(TYPE.Q1$gcn2with3ATvsCDS), Q2 = log2(TYPE.Q2$gcn2with3ATvsCDS), 
			Q3 = log2(TYPE.Q3$gcn2with3ATvsCDS), Q4 = log2(TYPE.Q4$gcn2with3ATvsCDS))
	ylim1112 <- c(min(c(unlist(temp11), unlist(temp12)), na.rm = TRUE), 
			max(c(unlist(temp11), unlist(temp12)), na.rm = TRUE))
	main9 <- "WTwith3ATvsCDS in groups Q1-Q4"
	main10 <- "gcn2with3ATvsCDS in groups Q1-Q4"
	ylab11 <- "uORF/CDS with 3AT"
	ylab12 <- "uORF/CDS with 3AT"

	pQ1 <- formatC(t.test(temp11[[1]], temp12[[1]])$p.value, format = "e", digits = 2)
	pQ2 <- formatC(t.test(temp11[[2]], temp12[[2]])$p.value, format = "e", digits = 2)
	pQ3 <- formatC(t.test(temp11[[3]], temp12[[3]])$p.value, format = "e", digits = 2)
	pQ4 <- formatC(t.test(temp11[[4]], temp12[[4]])$p.value, format = "e", digits = 2)
	nQ1 <- length(temp11[[1]])
	nQ2 <- length(temp11[[2]])
	nQ3 <- length(temp11[[3]])
	nQ4 <- length(temp11[[4]])
	sub1112 <- paste("n=", nQ1, ",", nQ2, ",", nQ3, ",", nQ4, sep = "")
	pvalue1112 <- as.numeric(c(pQ1, pQ2, pQ3, pQ4))

	LOG1[65] <- ""
	LOG1[66] <- main9
	LOG1[67] <- main10
	LOG1[68] <- sub1112
	LOG1[69] <- paste("p=", pQ1, pQ2, pQ3, pQ4, sep = ",")

	temp11 <- list(Q1 = log2(TYPE.Q1$ratioWith3AT), Q2 = log2(TYPE.Q2$ratioWith3AT), 
			Q3 = log2(TYPE.Q3$ratioWith3AT), Q4 = log2(TYPE.Q4$ratioWith3AT))
	temp12 <- list(Q1 = log2(TYPE.Q1$ratioWo3AT), Q2 = log2(TYPE.Q2$ratioWo3AT), 
			Q3 = log2(TYPE.Q3$ratioWo3AT), Q4 = log2(TYPE.Q4$ratioWo3AT))
	ylim1112 <- c(min(c(unlist(temp11), unlist(temp12)), na.rm = TRUE), 
			max(c(unlist(temp11), unlist(temp12)), na.rm = TRUE))
	main9 <- "ratioWith3AT in groups Q1-Q4"
	main10 <- "ratioWo3AT in groups Q1-Q4"
	ylab11 <- "(uORF/CDS in gcn2) / (uORF/CDS in WT)"
	ylab12 <- "(uORF/CDS in gcn2) / (uORF/CDS in WT)"

	pQ1 <- formatC(t.test(temp11[[1]], temp12[[1]])$p.value, format = "e", digits = 2)
	pQ2 <- formatC(t.test(temp11[[2]], temp12[[2]])$p.value, format = "e", digits = 2)
	pQ3 <- formatC(t.test(temp11[[3]], temp12[[3]])$p.value, format = "e", digits = 2)
	pQ4 <- formatC(t.test(temp11[[4]], temp12[[4]])$p.value, format = "e", digits = 2)
	nQ1 <- length(temp11[[1]])
	nQ2 <- length(temp11[[2]])
	nQ3 <- length(temp11[[3]])
	nQ4 <- length(temp11[[4]])
	sub1112 <- paste("n=", nQ1, ",", nQ2, ",", nQ3, ",", nQ4, sep = "")
	pvalue1112 <- as.numeric(c(pQ1, pQ2, pQ3, pQ4))

	boxplot(temp11, outline = FALSE, main = main9, sub = sub1112, cex.axis = 1, ylab = ylab11, ylim = ylim1112)
	abline(h = 0, col = "gray")
	beeswarm(temp11, corral = "wrap", col = rgb(1, 0, 0, alpha = 0.5), pch = 20, add = TRUE)

	boxplot(temp12, outline = FALSE, main = main10, sub = sub1112, cex.axis = 1, ylab = ylab12, ylim = ylim1112)
	y <- c(max(temp12[[1]]), max(temp12[[2]]), max(temp12[[3]]), max(temp12[[4]])) + 0.1
	abline(h = 0, col = "gray")
	beeswarm(temp12, corral = "wrap", col = rgb(1, 0, 0, alpha = 0.5), pch = 20, add = TRUE)
	dev.off() 

	LOG1[70] <- ""
	LOG1[71] <- main9
	LOG1[72] <- main10
	LOG1[73] <- sub1112
	LOG1[74] <- paste("p=", pQ1, pQ2, pQ3, pQ4, sep = ",")

	LOG1[50+25] <- ""
	LOG1[51+25] <- paste("Sequence for Q1-Q4 (", type, ")", sep = "")
	LOG1[52+25] <- paste("nrow(TYPE.Q1): ", nrow(TYPE.Q1), sep = "")
	LOG1[53+25] <- paste(" is_masked_start == FALS: ", nrow(subset(TYPE.Q1, is_masked_start == FALSE)), sep = "")
	LOG1[54+25] <- paste(" is_masked_start == TRUE: ", nrow(subset(TYPE.Q1, is_masked_start == TRUE)), sep = "")
	LOG1[55+25] <- paste(" is.na(is_masked_start) == TRUE): ", 
			nrow(subset(TYPE.Q1, is.na(is_masked_start) == TRUE)), sep = "")
	LOG1[56+25] <- paste(" is_masked_end == FALSE): ", nrow(subset(TYPE.Q1, is_masked_end == FALSE)), sep = "")
	LOG1[57+25] <- paste(" is_masked_end == TRUE: ", nrow(subset(TYPE.Q1, is_masked_end == TRUE)), sep = "")
	LOG1[58+25] <- paste(" is.na(is_masked_end) == TRUE: ", 
			nrow(subset(TYPE.Q1, is.na(is_masked_end) == TRUE)), sep = "")

	LOG1[59+25] <- paste("nrow(TYPE.Q2): ", nrow(TYPE.Q2), sep = "")
	LOG1[60+25] <- paste(" is_masked_start == FALS: ", nrow(subset(TYPE.Q2, is_masked_start == FALSE)), sep = "")
	LOG1[61+25] <- paste(" is_masked_start == TRUE: ", nrow(subset(TYPE.Q2, is_masked_start == TRUE)), sep = "")
	LOG1[62+25] <- paste(" is.na(is_masked_start) == TRUE): ", 
			nrow(subset(TYPE.Q2, is.na(is_masked_start) == TRUE)), sep = "")
	LOG1[63+25] <- paste(" is_masked_end == FALSE): ", nrow(subset(TYPE.Q2, is_masked_end == FALSE)), sep = "")
	LOG1[64+25] <- paste(" is_masked_end == TRUE: ", nrow(subset(TYPE.Q2, is_masked_end == TRUE)), sep = "")
	LOG1[65+25] <- paste(" is.na(is_masked_end) == TRUE: ", 
			nrow(subset(TYPE.Q2, is.na(is_masked_end) == TRUE)), sep = "")

	LOG1[66+25] <- paste("nrow(TYPE.Q3): ", nrow(TYPE.Q3), sep = "")
	LOG1[67+25] <- paste(" is_masked_start == FALS: ", nrow(subset(TYPE.Q3, is_masked_start == FALSE)), sep = "")
	LOG1[68+25] <- paste(" is_masked_start == TRUE: ", nrow(subset(TYPE.Q3, is_masked_start == TRUE)), sep = "")
	LOG1[69+25] <- paste(" is.na(is_masked_start) == TRUE): ", 
			nrow(subset(TYPE.Q3, is.na(is_masked_start) == TRUE)), sep = "")
	LOG1[70+25] <- paste(" is_masked_end == FALSE): ", nrow(subset(TYPE.Q3, is_masked_end == FALSE)), sep = "")
	LOG1[71+25] <- paste(" is_masked_end == TRUE: ", nrow(subset(TYPE.Q3, is_masked_end == TRUE)), sep = "")
	LOG1[72+25] <- paste(" is.na(is_masked_end) == TRUE: ", 
			nrow(subset(TYPE.Q3, is.na(is_masked_end) == TRUE)), sep = "")

	LOG1[73+25] <- paste("nrow(TYPE.Q4): ", nrow(TYPE.Q4), sep = "")
	LOG1[74+25] <- paste(" is_masked_start == FALS: ", nrow(subset(TYPE.Q4, is_masked_start == FALSE)), sep = "")
	LOG1[75+25] <- paste(" is_masked_start == TRUE: ", nrow(subset(TYPE.Q4, is_masked_start == TRUE)), sep = "")
	LOG1[76+25] <- paste(" is.na(is_masked_start) == TRUE): ", 
			nrow(subset(TYPE.Q4, is.na(is_masked_start) == TRUE)), sep = "")
	LOG1[77+25] <- paste(" is_masked_end == FALSE): ", nrow(subset(TYPE.Q4, is_masked_end == FALSE)), sep = "")
	LOG1[78+25] <- paste(" is_masked_end == TRUE: ", nrow(subset(TYPE.Q4, is_masked_end == TRUE)), sep = "")
	LOG1[79+25] <- paste(" is.na(is_masked_end) == TRUE: ", 
			nrow(subset(TYPE.Q4, is.na(is_masked_end) == TRUE)), sep = "")

	TYPE.total.start <- DNAStringSet(subset(selectedDF, is_masked_start == FALSE)$unmasked_start)
	TYPE.total.end <- DNAStringSet(subset(selectedDF, is_masked_end == FALSE)$unmasked_end)
	names(TYPE.total.start) <- subset(selectedDF, is_masked_start == FALSE)$name
	names(TYPE.total.end) <- subset(selectedDF, is_masked_end == FALSE)$name

	TYPE.Q1.start <- DNAStringSet(subset(TYPE.Q1, is_masked_start == FALSE)$unmasked_start)
	TYPE.Q1.end <- DNAStringSet(subset(TYPE.Q1, is_masked_end == FALSE)$unmasked_end)
	names(TYPE.Q1.start) <- subset(TYPE.Q1, is_masked_start == FALSE)$name
	names(TYPE.Q1.end) <- subset(TYPE.Q1, is_masked_end == FALSE)$name

	TYPE.Q2.start <- DNAStringSet(subset(TYPE.Q2, is_masked_start == FALSE)$unmasked_start)
	TYPE.Q2.end <- DNAStringSet(subset(TYPE.Q2, is_masked_end == FALSE)$unmasked_end)
	names(TYPE.Q2.start) <- subset(TYPE.Q2, is_masked_start == FALSE)$name
	names(TYPE.Q2.end) <- subset(TYPE.Q2, is_masked_end == FALSE)$name

	TYPE.Q3.start <- DNAStringSet(subset(TYPE.Q3, is_masked_start == FALSE)$unmasked_start)
	TYPE.Q3.end <- DNAStringSet(subset(TYPE.Q3, is_masked_end == FALSE)$unmasked_end)
	names(TYPE.Q3.start) <- subset(TYPE.Q3, is_masked_start == FALSE)$name
	names(TYPE.Q3.end) <- subset(TYPE.Q3, is_masked_end == FALSE)$name

	TYPE.Q4.start <- DNAStringSet(subset(TYPE.Q4, is_masked_start == FALSE)$unmasked_start)
	TYPE.Q4.end <- DNAStringSet(subset(TYPE.Q4, is_masked_end == FALSE)$unmasked_end)
	names(TYPE.Q4.start) <- subset(TYPE.Q4, is_masked_start == FALSE)$name
	names(TYPE.Q4.end) <- subset(TYPE.Q4, is_masked_end == FALSE)$name

	LOG1[80+25] <- ""
	LOG1[81+25] <- "DNAStringSet"
	LOG1[82+25] <- paste("TYPE.total.start: ", length(TYPE.total.start), sep = "")
	LOG1[83+25] <- paste("TYPE.total.end: ", length(TYPE.total.end), sep = "")
	LOG1[84+25] <- paste("TYPE.Q1.start: ", length(TYPE.Q1.start), sep = "")
	LOG1[85+25] <- paste("TYPE.Q1.end: ", length(TYPE.Q1.end), sep = "")
	LOG1[86+25] <- paste("TYPE.Q2.start: ", length(TYPE.Q2.start), sep = "")
	LOG1[87+25] <- paste("TYPE.Q2.end: ", length(TYPE.Q2.end), sep = "")
	LOG1[88+25] <- paste("TYPE.Q3.start: ", length(TYPE.Q3.start), sep = "")
	LOG1[89+25] <- paste("TYPE.Q3.end: ", length(TYPE.Q3.end), sep = "")
	LOG1[90+25] <- paste("TYPE.Q4.start: ", length(TYPE.Q4.start), sep = "")
	LOG1[91+25] <- paste("TYPE.Q4.end: ", length(TYPE.Q4.end), sep = "")

	setwd(wd)
	FASTADIR <- paste("FASTA_RPKM", RPKM.threshold, "_WIDTH", WIDTH, "_v11", sep = "")
	dir.create(FASTADIR)
	setwd(FASTADIR)

	writeXStringSet(TYPE.total.start, 
		filepath = paste(type, "_total_start_th", RPKM.threshold, "_W", WIDTH, ".fasta", sep = ""))
	writeXStringSet(TYPE.Q1.start, 
		filepath = paste(type, "_Q1_start_th", RPKM.threshold, "_W", WIDTH, ".fasta", sep = ""))
	writeXStringSet(TYPE.Q2.start, 
		filepath = paste(type, "_Q2_start_th", RPKM.threshold, "_W", WIDTH, ".fasta", sep = ""))
	writeXStringSet(TYPE.Q3.start, 
		filepath = paste(type, "_Q3_start_th", RPKM.threshold, "_W", WIDTH, ".fasta", sep = ""))
	writeXStringSet(TYPE.Q4.start, 
		filepath = paste(type, "_Q4_start_th", RPKM.threshold, "_W", WIDTH, ".fasta", sep = ""))

	writeXStringSet(TYPE.total.end, 
		filepath = paste(type, "_total_end_th", RPKM.threshold, "_W", WIDTH, ".fasta", sep = ""))
	writeXStringSet(TYPE.Q1.end, 
		filepath = paste(type, "_Q1_end_th", RPKM.threshold, "_W", WIDTH, ".fasta", sep = ""))
	writeXStringSet(TYPE.Q2.end, 
		filepath = paste(type, "_Q2_end_th", RPKM.threshold, "_W", WIDTH, ".fasta", sep = ""))
	writeXStringSet(TYPE.Q3.end, 
		filepath = paste(type, "_Q3_end_th", RPKM.threshold, "_W", WIDTH, ".fasta", sep = ""))
	writeXStringSet(TYPE.Q4.end, 
		filepath = paste(type, "_Q4_end_th", RPKM.threshold, "_W", WIDTH, ".fasta", sep = ""))

	setwd(wd)
	LogoDIR <- paste("Logo_RPKM", RPKM.threshold, "_WIDTH", WIDTH, "_v11", sep = "")
	dir.create(LogoDIR)
	setwd(LogoDIR)
	
	seqData <- list(start_total = as.character(TYPE.total.start), 
			start_Q1 = as.character(TYPE.Q1.start), 
			start_Q2 = as.character(TYPE.Q2.start), 
			start_Q3 = as.character(TYPE.Q3.start), 
			start_Q4 = as.character(TYPE.Q4.start), 
			end_total = as.character(TYPE.total.end), 
			end_Q1 = as.character(TYPE.Q1.end), 
			end_Q2 = as.character(TYPE.Q2.end), 
			end_Q3 = as.character(TYPE.Q3.end), 
			end_Q4 = as.character(TYPE.Q4.end))

	x <- ggseqlogo(seqData, ncol = 1, method = "bits")
	filename <- paste(type, "_th", RPKM.threshold, "_W", WIDTH, "_Logo_bits.pdf", sep = "")
	pdf(file = filename, width = 12, height = 12)
	print(x)
	dev.off()

	x <- ggseqlogo(seqData, ncol = 1, method = "probability")
	filename <- paste(type, "_th", RPKM.threshold, "_W", WIDTH, "_Logo_prob.pdf", sep = "")
	pdf(file = filename, width = 12, height = 12)
	print(x)
	dev.off()

	setwd(wd)
	LenDIR <- paste("Length_RPKM", RPKM.threshold, "_WIDTH", WIDTH, "_v11", sep = "")
	dir.create(LenDIR)
	setwd(LenDIR)

	tempList <- list()

	for(q in 1:4){
		objName <- paste("TYPE.Q", q, sep = "")
		tempList[[q]] <- list()
		names(tempList)[q] <- paste("Q", q, sep = "")
		tempList[[q]][[1]] <- log2(selectedDF[[type]])
		for(i in 1:9){
			tempList[[q]][[i]] <- log2(subset(get(objName), length == (i+1)*3)[[type]])
			names(tempList[[q]])[i] <- paste("aa", i, sep = "")
			if(i == 9){
				tempList[[q]][[i]] <- log2(subset(get(objName), length > 27)[[type]])
				names(tempList[[q]])[i] <- ">8"
			}
		} # i-loop
	} # q-loop
	
	ylim <- c(min(log2(selectedDF[[type]]), na.rm = TRUE), max(log2(selectedDF[[type]]), na.rm = TRUE))
	filename <- paste(type, "_RPKMthresh", RPKM.threshold, "_len.pdf", sep = "")
	pdf(file = filename, width = 8, height = 12)
	par(mfrow = c(4, 1))
	for(q in 1:4){
		if(type == "WTwo3ATvsCDS")	ylab <- "uORF/CDS wo 3AT"
		if(type == "WTwith3ATvsCDS")	ylab <- "uORF/CDS with 3AT"
		if(type == "gcn2wo3ATvsCDS")	ylab <- "uORF/CDS wo 3AT"
		if(type == "gcn2with3ATvsCDS")	ylab <- "uORF/CDS with 3AT"
		if(type == "ratioWT")	ylab <- "(uORF/CDS with 3AT) / (uORF/CDS wo 3AT)"
		if(type == "ratiogcn2")	ylab <- "(uORF/CDS with 3AT) / (uORF/CDS wo 3AT)"
		if(type == "ratioWo3AT")	ylab <- "(uORF/CDS in gcn2) / (uORF/CDS in WT)"
		if(type == "ratioWith3AT")	ylab <- "(uORF/CDS in gcn2) / (uORF/CDS in WT)"
		
			boxplot(tempList[[q]], outline = FALSE, main = paste(type, " Q", q, sep = ""), 
				cex.axis = 1, ylab = ylab, ylim = ylim)
			abline(h = 0, col = "gray")
			beeswarm(tempList[[q]], corral = "wrap", col = rgb(1, 0, 0, alpha = 0.5), 
				pch = 20, add = TRUE)
	} # q-loop
	dev.off()

	filename <- paste(type, "_RPKMthresh", RPKM.threshold, "_len_bar.pdf", sep = "")
	pdf(file = filename, width = 8, height = 12)
	par(mfrow = c(4, 1))
	scoreMatrix <- data.frame(matrix(nrow = 4, ncol = 9, data = numeric()), stringsAsFactors = FALSE)
		for(i in 1:9){
			colnames(scoreMatrix)[i] <- paste("aa", i, sep = "")
			if(i == 9)	colnames(scoreMatrix)[i] <- paste(">8", sep = "")
		}
	for(q in 1:4){
		scores <- numeric(length = 9)
		for(i in 1:9){
			scores[i] <- length(tempList[[q]][[i]])
		}
		scoreMatrix[q,] <- scores
	} # q-loop
	pvalues <- numeric(length = 9)
	for(i in 1:9){
		pvalues[i] <- formatC(chisq.test(scoreMatrix[,i])$p.value, format = "e", digits = 1)
	}
	for(q in 1:4){
		scores <- as.numeric(scoreMatrix[q,])
		names(scores) <- colnames(scoreMatrix)
		if(q != 4){
			y_coords = barplot(scores, ylab = "Occurrence", 
					main = paste(type, " Q", q, sep = ""), cex.names = 1)
			text(y_coords, scores, labels = scores, pos = 1)
		}
		if(q == 4){
			y_coords = barplot(scores, ylab = "Occurrence", 
					sub = paste(c("chisq:", pvalues), collapse = ", "), 
					main = paste(type, " Q", q, sep = ""), cex.names = 1)
			text(y_coords, scores, labels = scores, pos = 1)
		}
	} # q-loop
	dev.off()


	for(i in 1:9){
		TYPE.Q1.start <- DNAStringSet(subset(TYPE.Q1, length == i*3+3 & is_masked_start == FALSE)$unmasked_start)
		TYPE.Q1.end <- DNAStringSet(subset(TYPE.Q1, length == i*3+3 & is_masked_end == FALSE)$unmasked_end)
		TYPE.Q2.start <- DNAStringSet(subset(TYPE.Q2, length == i*3+3 & is_masked_start == FALSE)$unmasked_start)
		TYPE.Q2.end <- DNAStringSet(subset(TYPE.Q2, length == i*3+3 & is_masked_end == FALSE)$unmasked_end)
		TYPE.Q3.start <- DNAStringSet(subset(TYPE.Q3, length == i*3+3 & is_masked_start == FALSE)$unmasked_start)
		TYPE.Q3.end <- DNAStringSet(subset(TYPE.Q3, length == i*3+3 & is_masked_end == FALSE)$unmasked_end)
		TYPE.Q4.start <- DNAStringSet(subset(TYPE.Q4, length == i*3+3 & is_masked_start == FALSE)$unmasked_start)
		TYPE.Q4.end <- DNAStringSet(subset(TYPE.Q4, length == i*3+3 & is_masked_end == FALSE)$unmasked_end)
		
		names(TYPE.Q1.start) <- subset(TYPE.Q1, length == i*3+3 & is_masked_start == FALSE)$name
		names(TYPE.Q1.end) <- subset(TYPE.Q1, length == i*3+3 & is_masked_end == FALSE)$name
		names(TYPE.Q2.start) <- subset(TYPE.Q2, length == i*3+3 & is_masked_start == FALSE)$name
		names(TYPE.Q2.end) <- subset(TYPE.Q2, length == i*3+3 & is_masked_end == FALSE)$name
		names(TYPE.Q3.start) <- subset(TYPE.Q3, length == i*3+3 & is_masked_start == FALSE)$name
		names(TYPE.Q3.end) <- subset(TYPE.Q3, length == i*3+3 & is_masked_end == FALSE)$name
		names(TYPE.Q4.start) <- subset(TYPE.Q4, length == i*3+3 & is_masked_start == FALSE)$name
		names(TYPE.Q4.end) <- subset(TYPE.Q4, length == i*3+3 & is_masked_end == FALSE)$name

		LOG1[116+(i-1)*8+1] <- paste("AA", i, ", start_Q1: ", length(TYPE.Q1.start), sep = "")
		LOG1[116+(i-1)*8+2] <- paste("AA", i, ",   end_Q1: ", length(TYPE.Q1.end), sep = "")
		LOG1[116+(i-1)*8+3] <- paste("AA", i, ", start_Q2: ", length(TYPE.Q2.start), sep = "")
		LOG1[116+(i-1)*8+4] <- paste("AA", i, ",   end_Q2: ", length(TYPE.Q2.end), sep = "")
		LOG1[116+(i-1)*8+5] <- paste("AA", i, ", start_Q3: ", length(TYPE.Q3.start), sep = "")
		LOG1[116+(i-1)*8+6] <- paste("AA", i, ",   end_Q3: ", length(TYPE.Q3.end), sep = "")
		LOG1[116+(i-1)*8+7] <- paste("AA", i, ", start_Q4: ", length(TYPE.Q4.start), sep = "")
		LOG1[116+(i-1)*8+8] <- paste("AA", i, ",   end_Q4: ", length(TYPE.Q4.end), sep = "")

		setwd(wd)
		FASTADIR2 <- paste("FASTA2_RPKM", RPKM.threshold, "_WIDTH", WIDTH, "_v11", sep = "")
		dir.create(FASTADIR2)
		setwd(FASTADIR2)

		writeXStringSet(TYPE.Q1.start, filepath = paste(type, "_Q1_start_th", 
			RPKM.threshold, "_W", WIDTH, "_AA", i, ".fasta", sep = ""))
		writeXStringSet(TYPE.Q2.start, filepath = paste(type, "_Q2_start_th", 
			RPKM.threshold, "_W", WIDTH, "_AA", i, ".fasta", sep = ""))
		writeXStringSet(TYPE.Q3.start, filepath = paste(type, "_Q3_start_th", 
			RPKM.threshold, "_W", WIDTH, "_AA", i, ".fasta", sep = ""))
		writeXStringSet(TYPE.Q4.start, filepath = paste(type, "_Q4_start_th", 
			RPKM.threshold, "_W", WIDTH, "_AA", i, ".fasta", sep = ""))

		writeXStringSet(TYPE.Q1.end, filepath = paste(type, "_Q1_end_th", 
			RPKM.threshold, "_W", WIDTH, "_AA", i, ".fasta", sep = ""))
		writeXStringSet(TYPE.Q2.end, filepath = paste(type, "_Q2_end_th", 
			RPKM.threshold, "_W", WIDTH, "_AA", i, ".fasta", sep = ""))
		writeXStringSet(TYPE.Q3.end, filepath = paste(type, "_Q3_end_th", 
			RPKM.threshold, "_W", WIDTH, "_AA", i, ".fasta", sep = ""))
		writeXStringSet(TYPE.Q4.end, filepath = paste(type, "_Q4_end_th", 
			RPKM.threshold, "_W", WIDTH, "_AA", i, ".fasta", sep = ""))

		setwd(wd)
		Logo2DIR <- paste("Logo2_RPKM", RPKM.threshold, "_WIDTH", WIDTH, "_v11", sep = "")
		if(i == 1)	dir.create(Logo2DIR)
		setwd(Logo2DIR)
		logo2 <- getwd()

		seqData <- list(start_Q1 = as.character(TYPE.Q1.start), 
			start_Q2 = as.character(TYPE.Q2.start), 
			start_Q3 = as.character(TYPE.Q3.start), 
			start_Q4 = as.character(TYPE.Q4.start), 
			end_Q1 = as.character(TYPE.Q1.end), 
			end_Q2 = as.character(TYPE.Q2.end), 
			end_Q3 = as.character(TYPE.Q3.end), 
			end_Q4 = as.character(TYPE.Q4.end))

		x <- NA
		x <- ggseqlogo(seqData, ncol = 1, method = "bits")
		filename <- paste(type, "_th", RPKM.threshold, 
				"_W", WIDTH, "_AA", i, "_Logo_bits.pdf", sep = "")
		pdf(file = filename, width = 12, height = 10)
		print(x)
		dev.off()

		x <- NA
		x <- ggseqlogo(seqData, ncol = 1, method = "probability")
		filename <- paste(type, "_th", RPKM.threshold, 
				"_W", WIDTH, "_AA", i, "_Logo_prob.pdf", sep = "")
		pdf(file = filename, width = 12, height = 10)
		print(x)
		dev.off()

	} # i-loop


	setwd(wd)
	AADIR <- paste("AA_RPKM", RPKM.threshold, "_WIDTH", WIDTH, "_v11", sep = "")
	dir.create(AADIR)
	setwd(AADIR)

	for(i in 1:((WIDTH/3)-1)){
	TYPE.Q1.start <- DNAStringSet(subset(TYPE.Q1, length == i*3+3 & is_masked_start == FALSE)$unmasked_start)
	TYPE.Q2.start <- DNAStringSet(subset(TYPE.Q2, length == i*3+3 & is_masked_start == FALSE)$unmasked_start)
	TYPE.Q3.start <- DNAStringSet(subset(TYPE.Q3, length == i*3+3 & is_masked_start == FALSE)$unmasked_start)
	TYPE.Q4.start <- DNAStringSet(subset(TYPE.Q4, length == i*3+3 & is_masked_start == FALSE)$unmasked_start)
		
	temp <- oligonucleotideFrequency(DNAStringSet("AAA"), width = 3, step = 3, simplify.as = "collapsed")
	codonLetters <- names(temp)

	AAs <- 1:(i+1)
	for(a in AAs){
	codonDF <- data.frame(matrix(nrow = 4, ncol = 64, data = as.integer(0)), stringsAsFactors = FALSE)
	colnames(codonDF) <- codonLetters
	rownames(codonDF) <- c("Q1", "Q2", "Q3", "Q4")
	rownames(codonDF)[1] <- paste("len", i*3+3, "_AA", a, "_Q1", sep = "")
	rownames(codonDF)[2] <- paste("len", i*3+3, "_AA", a, "_Q2", sep = "")
	rownames(codonDF)[3] <- paste("len", i*3+3, "_AA", a, "_Q3", sep = "")
	rownames(codonDF)[4] <- paste("len", i*3+3, "_AA", a, "_Q4", sep = "")
	for(q in 1:4){
		objName <- paste("TYPE.Q", q, ".start", sep = "")
		codons <- as.character(subseq(get(objName), start = WIDTH+(a-1)*3+1, end = WIDTH+(a-1)*3+3))
		for(c in 1:length(codons)){
			codon <- codons[c]
			codonDF[[codon]][q] <- codonDF[[codon]][q] + 1
		} 
	} # q-loop
	if(a == 1) codonDF2 <- codonDF
	if(a != 1) codonDF2 <- rbind(codonDF2, codonDF)
	filename <- paste(type, "_RPKMthresh", RPKM.threshold, "_len", i*3+3, "_AA", a, ".xlsx", sep = "")
	write.xlsx(codonDF, file = filename, rowNames = TRUE)
	} # a-loop
	filename <- paste(type, "_RPKMthresh", RPKM.threshold, "_len", i*3+3, "_bundled.xlsx", sep = "")
	write.xlsx(codonDF2, file = filename, rowNames = TRUE)

	selectedQ1 <- codonDF2[grep(pattern = "Q1", rownames(codonDF2)),]
	selectedQ2 <- codonDF2[grep(pattern = "Q2", rownames(codonDF2)),]
	selectedQ3 <- codonDF2[grep(pattern = "Q3", rownames(codonDF2)),]
	selectedQ4 <- codonDF2[grep(pattern = "Q4", rownames(codonDF2)),]

	selectedQ1noSTOP <- selectedQ1[1:(nrow(selectedQ1)-1),]
	selectedQ2noSTOP <- selectedQ2[1:(nrow(selectedQ2)-1),]
	selectedQ3noSTOP <- selectedQ3[1:(nrow(selectedQ3)-1),]
	selectedQ4noSTOP <- selectedQ4[1:(nrow(selectedQ4)-1),]

	if(i != 1){
		selectedQ1noATG <- selectedQ1[2:(nrow(selectedQ1)-1),]
		selectedQ2noATG <- selectedQ2[2:(nrow(selectedQ2)-1),]
		selectedQ3noATG <- selectedQ3[2:(nrow(selectedQ3)-1),]
		selectedQ4noATG <- selectedQ4[2:(nrow(selectedQ4)-1),]
	}
	if(i == 1){
		selectedQ1noATG <- selectedQ1[0,]
		selectedQ2noATG <- selectedQ2[0,]
		selectedQ3noATG <- selectedQ3[0,]
		selectedQ4noATG <- selectedQ4[0,]
	}

	filename <- paste(type, "_RPKMthresh", RPKM.threshold, "_len", i*3+3, "_codons.pdf", sep = "")
	pdf(file = filename, width = 20, height = 10)
	par(mfrow = c(4, 1))
	for(q in 1:4){
		main <- paste(type, "_RPKMthresh", RPKM.threshold, "_len", i*3+3, ", sum (AA1-", i,") Q", q, sep = "")
		objName <- paste("selectedQ", q, "noSTOP", sep = "")
		scores <- colSums(get(objName))
		y_coords = barplot(scores, ylab = "Occurrence", main = main, cex.names = 0.8)
		text(y_coords, scores, labels = scores, pos = 1)
	}
	dev.off()

	if(i != 1){
	filename <- paste(type, "_RPKMthresh", RPKM.threshold, "_len", i*3+3, "_codonsWoSTART.pdf", sep = "")
	pdf(file = filename, width = 20, height = 10)
	par(mfrow = c(4, 1))
	for(q in 1:4){
		main <- paste(type, "_RPKMthresh", RPKM.threshold, "_len", i*3+3, ", sum (AA2-", i,") Q", q, sep = "")
		objName <- paste("selectedQ", q, "noATG", sep = "")
		scores <- colSums(get(objName))
		y_coords = barplot(scores, ylab = "Occurrence", main = main, cex.names = 0.8)
		text(y_coords, scores, labels = scores, pos = 1)
	}
	dev.off()
	}

	if(i == 1){
		selectedQ1noSTOPsum <- colSums(selectedQ1noSTOP)
		selectedQ2noSTOPsum <- colSums(selectedQ2noSTOP)
		selectedQ3noSTOPsum <- colSums(selectedQ3noSTOP)
		selectedQ4noSTOPsum <- colSums(selectedQ4noSTOP)
		selectedQ1noATGsum <- colSums(selectedQ1noATG)
		selectedQ2noATGsum <- colSums(selectedQ2noATG)
		selectedQ3noATGsum <- colSums(selectedQ3noATG)
		selectedQ4noATGsum <- colSums(selectedQ4noATG)
	}
	if(i != 1){
		selectedQ1noSTOPsum <- selectedQ1noSTOPsum + colSums(selectedQ1noSTOP)
		selectedQ2noSTOPsum <- selectedQ1noSTOPsum + colSums(selectedQ2noSTOP)
		selectedQ3noSTOPsum <- selectedQ1noSTOPsum + colSums(selectedQ3noSTOP)
		selectedQ4noSTOPsum <- selectedQ1noSTOPsum + colSums(selectedQ4noSTOP)
		selectedQ1noATGsum <- selectedQ1noATGsum + colSums(selectedQ1noATG)
		selectedQ2noATGsum <- selectedQ1noATGsum + colSums(selectedQ2noATG)
		selectedQ3noATGsum <- selectedQ1noATGsum + colSums(selectedQ3noATG)
		selectedQ4noATGsum <- selectedQ1noATGsum + colSums(selectedQ4noATG)
	}

	} # i-loop

	filename <- paste(type, "_RPKMthresh", RPKM.threshold, "_sum_codons.pdf", sep = "")
	pdf(file = filename, width = 20, height = 10)
	par(mfrow = c(4, 1))
	for(q in 1:4){
		main <- paste(type, "_RPKMthresh", RPKM.threshold, ", len6-", WIDTH-3, ", sum (AA1-", (WIDTH/3)-1, "), Q", q, sep = "")
		objName <- paste("selectedQ", q, "noSTOPsum", sep = "")
		scores <- get(objName)
		ymax <- max(scores, na.rm = TRUE) + max(scores, na.rm = TRUE) * 0.2
		y_coords = barplot(scores, ylab = "Occurrence", main = main, cex.names = 0.8, ylim = c(0, ymax))
		text(y_coords, scores, labels = scores, pos = 3)
	}
	dev.off()

	filename <- paste(type, "_RPKMthresh", RPKM.threshold, "_sum_codonsWoSTART.pdf", sep = "")
	pdf(file = filename, width = 20, height = 10)
	par(mfrow = c(4, 1))
	for(q in 1:4){
		main <- paste(type, "_RPKMthresh", RPKM.threshold, ", len6-", WIDTH-3, ", sum (AA2-", (WIDTH/3)-1, "), Q", q, sep = "")
		objName <- paste("selectedQ", q, "noATGsum", sep = "")
		scores <- get(objName)
		ymax <- max(scores, na.rm = TRUE) + max(scores, na.rm = TRUE) * 0.2
		y_coords = barplot(scores, ylab = "Occurrence", main = main, cex.names = 0.8, ylim = c(0, ymax))
		text(y_coords, scores, labels = scores, pos = 3)
	}
	dev.off()


	setwd(wd)
	LOGDIR <- paste("LOG_RPKM", RPKM.threshold, "_WIDTH", WIDTH, "_v11", sep = "")
	dir.create(LOGDIR)
	setwd(LOGDIR)
	filename <- paste(type, "_log_RPKMthresh", RPKM.threshold, "_WIDTH", WIDTH, ".txt", sep = "")
	write(LOG1, file = filename)

	setwd(wd3)
}

setwd(HOMEDIR)
dir.create("20200311_Logo_v11")
setwd("20200311_Logo_v11")

for(type in c("ratioWo3AT", "ratioWith3AT")){
	analyzeRatiouORFv11(inDF = bundledDF, type = type, WIDTH = 27, RPKM.threshold = 6)
}


