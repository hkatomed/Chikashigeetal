## The output files are used in the indicated figures.
# output: RPKMthresh6_ratioWT-ratiogcn2.pdf
#         used in Fig.4D
# output: Bases_RPKMthresh6_UPST1-4_selected_.pdf
#         used in S8A
# output: Bases_RPKMthresh6_UPST1_LEN18-24_selected_.pdf
#         used in S8B
# output: 20200422_Logo_v10_4/RPKM6/TRUEuORFs_20200311.xlsx
#         used in Table S5 (28 genes)


# Specify HOMEDIR before run.
HOMEDIR <- "XXX"


setwd(HOMEDIR)
setwd("pombase")
load(file = "RData/splicedSeqs2.RData")
load(file = "RData/splicedSeqsM2.RData")

setwd("/Volumes/Promise RAID/20190710_Ribosome_Asano")
setwd("RPKM_2")
setwd("RData")
load(file = "bundledDF_20200312.RData")


analyzeRatiouORFv10_4 <- function(inDF, RPKM.threshold = 6){

	wd3 <- getwd()
	
	setwd(wd3)
	OUTDIR1 <- paste("RPKM", RPKM.threshold, sep = "")
	dir.create(OUTDIR1)
	setwd(OUTDIR1)
	wd <- getwd()
	

	SAMPLE1 <- "WTwo3AT"
	SAMPLE2 <- "WTwith3AT"
	SAMPLE3 <- "gcn2wo3AT"
	SAMPLE4 <- "gcn2with3AT"

	LOG1 <- character(length = 12)
	# LOG1[1] <- paste("type: ", type, sep = "")
	LOG1[1] <- paste("type: ", "Not applicable", sep = "")
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
	checkuORFs(selectedDF, char = "2_uORF_selected")
	LOG1[20] <- paste(" uORF==TRUE: ", nrow(selectedDF), sep = "")
	selectedDF <- subset(selectedDF, CDSoverlap == FALSE)
	checkuORFs(selectedDF, char = "3_CDSoverlapRemoved")
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
	# LOG1[29] <- paste(" phase1 RPKM is more than twice higher than those of phases 2 and 3: ", 
			nrow(selectedDF), sep = "")
	LOG1[29] <- paste(" phase1 RPKM is ", FOLD, " fold-higher than those of phases 2 and 3: ", 
			nrow(selectedDF), sep = "")

	log2plotLen <- function(type1, type2){
		cat("log2plotLen, type1:", type1, ", type2:", type2, "\n", sep = "")
		start <- 1
		# end <- 9
		end <- 18
		filename <- paste("RPKMthresh", RPKM.threshold, "_", type1, "-", type2, ".pdf", sep = "")
		pdf(file = filename, width = 4*3, height = 4*3)
		par(mfrow = c(3, 3))
		lim <- c(min(c(selectedDF[[type1]], selectedDF[[type2]]), na.rm = TRUE), 
					max(c(selectedDF[[type1]], selectedDF[[type2]]), na.rm = TRUE))	
		xlab <- paste("log2(", type1, ")", sep = "")
		ylab <- paste("log2(", type2, ")", sep = "")
			
		for(i in start:(end+1)){
			if(i <= end){
				uORFlength <- i*3+3
				temp <- subset(selectedDF, length == uORFlength)
				main <- paste(type1, " vs ", type2, ", len=", uORFlength, ", n=", 
						nrow(temp), sep = "")
				plot(x = log2(selectedDF[[type1]]), y = log2(selectedDF[[type2]]), 
					xlim = log2(lim), ylim = log2(lim), 
					main = main, xlab = xlab, ylab = ylab, type = "n")
				abline(a = 0, b = 1, lwd = 2, col = "gray")
				# if(type1 == "WTwo3ATvsCDS" & type2 == "WTwith3ATvsCDS"){
					abline(h = log2(0.5), lwd = 2, col = "gray")
					abline(h = -log2(0.5), lwd = 2, col = "gray")
					abline(h = log2(1), lwd = 2, col = "gray")
					abline(v = log2(0.5), lwd = 2, col = "gray")
					abline(v = -log2(0.5), lwd = 2, col = "gray")
					abline(v = log2(1), lwd = 2, col = "gray")
				# }
				# if(type1 == "gcn2wo3ATvsCDS" & type2 == "gcn2with3ATvsCDS"){
					abline(h = log2(0.5), lwd = 2, col = "gray")
					abline(h = -log2(0.5), lwd = 2, col = "gray")
					abline(h = log2(1), lwd = 2, col = "gray")
					abline(v = log2(0.5), lwd = 2, col = "gray")
					abline(v = -log2(0.5), lwd = 2, col = "gray")
					abline(v = log2(1), lwd = 2, col = "gray")
				# }
				abline(a = log2(2), b = 1, lwd = 2, col = "orange")
				abline(a = -log2(2), b = 1, lwd = 2, col = "orange")
				abline(a = log2(1.5), b = 1, lwd = 2, col = "lightblue")
				abline(a = -log2(1.5), b = 1, lwd = 2, col = "lightblue")
				points(x = log2(selectedDF[[type1]]), y = log2(selectedDF[[type2]]), 
					pch = 20, col = gray(0, alpha = 0.2))
			}
			if(i == end+1){
				main <- paste(type1, " vs ", type2, ", total, n=", nrow(selectedDF), sep = "")
				plot(x = log2(selectedDF[[type1]]), y = log2(selectedDF[[type2]]), 
					xlim = log2(lim), ylim = log2(lim), 
					main = main, xlab = xlab, ylab = ylab, type = "n")
				abline(a = 0, b = 1, lwd = 2, col = "gray")
				if(type1 == "WTwo3ATvsCDS" & type2 == "WTwith3ATvsCDS"){
					abline(h = log2(0.5), lwd = 2, col = "gray")
					abline(h = -log2(0.5), lwd = 2, col = "gray")
					abline(h = log2(1), lwd = 2, col = "gray")
					abline(v = log2(0.5), lwd = 2, col = "gray")
					abline(v = -log2(0.5), lwd = 2, col = "gray")
					abline(v = log2(1), lwd = 2, col = "gray")
				}
				if(type1 == "gcn2wo3ATvsCDS" & type2 == "gcn2with3ATvsCDS"){
					abline(h = log2(0.5), lwd = 2, col = "gray")
					abline(h = -log2(0.5), lwd = 2, col = "gray")
					abline(h = log2(1), lwd = 2, col = "gray")
					abline(v = log2(0.5), lwd = 2, col = "gray")
					abline(v = -log2(0.5), lwd = 2, col = "gray")
					abline(v = log2(1), lwd = 2, col = "gray")
				}
					abline(a = log2(2), b = 1, lwd = 2, col = "orange")
					abline(a = -log2(2), b = 1, lwd = 2, col = "orange")
					abline(a = log2(1.5), b = 1, lwd = 2, col = "lightblue")
					abline(a = -log2(1.5), b = 1, lwd = 2, col = "lightblue")
					points(x = log2(selectedDF[[type1]]), y = log2(selectedDF[[type2]]), 
						pch = 20, col = gray(0, alpha = 0.2))
				
				if(type1 == "ratioWT" & type2 == "ratiogcn2"){
					plot(x = log2(selectedDF[[type1]]), y = log2(selectedDF[[type2]]), 
						xlim = log2(lim), ylim = log2(lim), 
						main = main, xlab = xlab, ylab = ylab, type = "n")
					trueuORFs <- subset(selectedDF, 
						selectedDF$ratioWT <= 0.5 & selectedDF$ratiogcn2 >= 1)
					falseuORFs <- selectedDF[(selectedDF$ratioWT <= 0.5 & selectedDF$ratiogcn2 >= 1) == FALSE,]
					cat("selectedDF: ", nrow(selectedDF), "\n", sep = "")
					cat("ratioWT <= 0.5 & ratiogcn2 >= 1: ", nrow(trueuORFs), "\n", sep = "")
					cat("ratioWT > 0.5 & ratiogcn2 < 1: ", nrow(falseuORFs), "\n", sep = "")

					abline(h = log2(1), lwd = 2, col = "gray")
					abline(v = log2(0.5), lwd = 2, col = "gray")
					abline(a = log2(2), b = 1, lwd = 2, col = "orange")
					abline(a = -log2(2), b = 1, lwd = 2, col = "orange")
					points(x = log2(trueuORFs[[type1]]), y = log2(trueuORFs[[type2]]), 
						pch = 20, col = rgb(1, 0, 0, alpha = 0.5))
					points(x = log2(falseuORFs[[type1]]), y = log2(falseuORFs[[type2]]), 
						pch = 20, col = gray(0, alpha = 0.2))

					plot(x = log2(selectedDF[[type1]]), y = log2(selectedDF[[type2]]), 
						xlim = log2(lim), ylim = log2(lim), 
						main = main, xlab = xlab, ylab = ylab, type = "n")
					abline(a = 0, b = 1, lwd = 2, col = "gray")
					abline(h = log2(1), lwd = 2, col = "gray")
					abline(v = log2(0.5), lwd = 2, col = "gray")
					points(x = log2(trueuORFs[[type1]]), y = log2(trueuORFs[[type2]]), 
						pch = 20, col = rgb(1, 0, 0, alpha = 0.5))
					points(x = log2(falseuORFs[[type1]]), y = log2(falseuORFs[[type2]]), 
						pch = 20, col = gray(0, alpha = 0.2))

					plot(x = log2(selectedDF[[type1]]), y = log2(selectedDF[[type2]]), 
						xlim = log2(lim), ylim = log2(lim), 
						main = main, xlab = xlab, ylab = ylab, type = "n")
					abline(a = 0, b = 1, lwd = 2, col = "gray")
					abline(h = log2(1), lwd = 2, col = "orange")
					abline(v = log2(0.5), lwd = 2, col = "orange")
					points(x = log2(trueuORFs[[type1]]), y = log2(trueuORFs[[type2]]), 
						pch = 20, col = rgb(1, 0, 0, alpha = 0.5))
					points(x = log2(falseuORFs[[type1]]), y = log2(falseuORFs[[type2]]), 
						pch = 20, col = gray(0, alpha = 0.2))
				}

			}

			# fil1-uORF1, pch=2, triangle
			# fil1-uORF4, pch=5, diamond
			# fil1-uORF5, pch=1, circle	
			# gcn5-uORF3, pch=6, inverted triangle
			# hri1-uORF2, pch=0, square

			index <- grep(pattern = "SPCC1393.08.1_uORF1", x = selectedDF$name)
			if(length(index) > 0){
				temp.target <- selectedDF[index,]
				points(x = log2(temp.target[[type1]]), y = log2(temp.target[[type2]]), 
					pch = 2, col = rgb(0, 1, 0, alpha = 1))
			}
			if(length(index) == 0){
				index2 <- grep(pattern = "SPCC1393.08.1_uORF1", x = inDF$name)
				temp.target <- inDF[index2,]
				points(x = log2(temp.target[[type1]]), y = log2(temp.target[[type2]]), 
					pch = 2, col = rgb(0, 0, 1, alpha = 1))
			}
			
			index <- grep(pattern = "SPCC1393.08.1_uORF4", x = selectedDF$name)
			if(length(index) > 0){
				temp.target <- selectedDF[index,]
				points(x = log2(temp.target[[type1]]), y = log2(temp.target[[type2]]), 
					pch = 5, col = rgb(0, 1, 0, alpha = 1))
			}
			if(length(index) == 0){
				index2 <- grep(pattern = "SPCC1393.08.1_uORF4", x = inDF$name)
				temp.target <- inDF[index2,]
				points(x = log2(temp.target[[type1]]), y = log2(temp.target[[type2]]), 
					pch = 5, col = rgb(0, 0, 1, alpha = 1))
			}
				
			index <- grep(pattern = "SPAC1952.05.1_uORF3", x = selectedDF$name)
			if(length(index) > 0){
				temp.target <- selectedDF[index,]
				points(x = log2(temp.target[[type1]]), y = log2(temp.target[[type2]]), 
					pch = 6, col = rgb(0, 1, 0, alpha = 1))
			}
			if(length(index) == 0){
				index2 <- grep(pattern = "SPAC1952.05.1_uORF3", x = inDF$name)
				temp.target <- inDF[index2,]
				points(x = log2(temp.target[[type1]]), y = log2(temp.target[[type2]]), 
					pch = 6, col = rgb(0, 0, 1, alpha = 1))
			}
			
			# index <- grep(pattern = "SPAC20G4.03c.1_uORF2", x = selectedDF$name)
			# if(length(index) > 0){
			# 	temp.target <- selectedDF[index,]
			# 	points(x = log2(temp.target[[type1]]), y = log2(temp.target[[type2]]), 
			# 		pch = 0, col = rgb(0, 1, 0, alpha = 1))
			# }
			# if(length(index) == 0){
			# 	index2 <- grep(pattern = "SPAC20G4.03c.1_uORF2", x = inDF$name)
			# 	temp.target <- inDF[index2,]
			# 	points(x = log2(temp.target[[type1]]), y = log2(temp.target[[type2]]), 
			# 		pch = 0, col = rgb(0, 0, 1, alpha = 1))
			# }

			if(i <= end){
				points(x = log2(temp[[type1]]), y = log2(temp[[type2]]), 
						pch = 20, col = rgb(1, 0, 0, alpha = 0.5))
	
				# ratio <- temp[[type2]]/temp[[type1]]
				ratio <- log2(temp[[type2]]/temp[[type1]])
				labels1 <- paste(length(ratio[which(ratio > log2(2))]), 
						length(ratio[which(ratio > log2(1.5))]), sep = ",")
				labels2 <- paste(length(ratio[which(ratio < -log2(2))]), 
						length(ratio[which(ratio < -log2(1.5))]), sep = ",")
				text(x = log2(lim[1])+1, y = log2(lim[2]), 
					labels = labels1, col = "red")
				text(x = log2(lim[2])-1, y = log2(lim[1]), 
					labels = labels2, col = "red")

				if(length(which(ratio > log2(1.5)))>0){
					out1 <- temp[which(ratio > log2(1.5)), c(1:9, 25:38, 47:50)]
					out1$type1 <- type1
					out1$type2 <- type2
					out1$upFOLD <- type2
					out1$dnFOLD <- "-"
					# out1$ratiolog2 <- log2(ratio[which(ratio > log2(1.5))])
					out1$ratiolog2 <- ratio[which(ratio > log2(1.5))]
					# out1 <- out1[c((ncol(out1)-4):ncol(out1), 1:27)]
					if(i == 1) out1o <- out1
					if(i != 1) out1o <- rbind(out1o, out1)
				}
				if(length(which(ratio < -log2(1.5)))>0){
					out2 <- temp[which(ratio < -log2(1.5)), c(1:9, 25:38, 47:50)]
					out2$type1 <- type1
					out2$type2 <- type2
					out2$upFOLD <- "-"
					out2$dnFOLD <- type1
					# out2$ratiolog2 <- log2(ratio[which(ratio < -log2(1.5))])
					out2$ratiolog2 <- ratio[which(ratio < -log2(1.5))]
					# out2 <- out2[c((ncol(out2)-4):ncol(out2), 1:27)]
					if(i == 1) out2o <- out2
					if(i != 1) out2o <- rbind(out2o, out2)
				}
			}
		}
		dev.off()

		for(i in (end+1):max(selectedDF$length)){
			uORFlength <- i*3+3
			temp <- subset(selectedDF, length == uORFlength)
			# ratio <- temp[[type2]]/temp[[type1]]
			ratio <- log2(temp[[type2]]/temp[[type1]])
			if(length(which(ratio > log2(1.5)))>0){
				out1 <- temp[which(ratio > log2(1.5)), c(1:9, 25:38, 47:50)]
				out1$type1 <- type1
				out1$type2 <- type2
				out1$upFOLD <- type2
				out1$dnFOLD <- "-"
				# out1$ratiolog2 <- log2(ratio[which(ratio > log2(1.5))])
				out1$ratiolog2 <- ratio[which(ratio > log2(1.5))]
				# out1 <- out1[c((ncol(out1)-4):ncol(out1)2, 1:27)]
				if(i == 1) out1o <- out1
				if(i != 1) out1o <- rbind(out1o, out1)
			}
			if(length(which(ratio < -log2(1.5)))>0){
				out2 <- temp[which(ratio < -log2(1.5)), c(1:9, 25:38, 47:50)]
				out2$type1 <- type1
				out2$type2 <- type2
				out2$upFOLD <- "-"
				out2$dnFOLD <- type1
				# out2$ratiolog2 <- log2(ratio[which(ratio < -log2(1.5))])
				out2$ratiolog2 <- ratio[which(ratio < -log2(1.5))]
				# out2 <- out2[c((ncol(out2)-4):ncol(out2), 1:27)]
				if(i == 1) out2o <- out2
				if(i != 1) out2o <- rbind(out2o, out2)
			}
		}

		out1o <- subset(out1o, is.na(transcript_id) == FALSE)
		out2o <- subset(out2o, is.na(transcript_id) == FALSE)
		if(nrow(out1o) > 0 & nrow(out2o) > 0){
			out <- rbind(out1o, out2o)
			if(type1 == "ratioWT" & type2 == "ratiogcn2"){
				outList <- list(out = out, trueuORFs = trueuORFs, falseuORFs = falseuORFs)
			}
			if(type1 != "ratioWT" | type2 != "ratiogcn2"){
				outList <- list(out = out, trueuORFs = list(), falseuORFs = list())
			}
			return(outList)
		}
		if(nrow(out1o) > 0){
			if(type1 == "ratioWT" & type2 == "ratiogcn2"){
				outList <- list(out = out1o, trueuORFs = trueuORFs, falseuORFs = falseuORFs)
			}
			if(type1 != "ratioWT" | type2 != "ratiogcn2"){
				outList <- list(out = out1o, trueuORFs = list(), falseuORFs = list())
			}
			return(outList)
		}
		if(nrow(out2o) > 0){
			if(type1 == "ratioWT" & type2 == "ratiogcn2"){
				outList <- list(out = out2o, trueuORFs = trueuORFs, falseuORFs = falseuORFs)
			}
			if(type1 != "ratioWT" | type2 != "ratiogcn2"){
				outList <- list(out = out2o, trueuORFs = list(), falseuORFs = list())
			}
			return(outList)
		}
	}


	outList <- log2plotLen(type1 = "ratioWT", type2 = "ratiogcn2")
	outList$out <- outList$out[c((ncol(outList$out)-4):ncol(outList$out), 1:(ncol(outList$out)-5))]

	filename <- "DiffFOLD_20200311.xlsx"
	write.xlsx(x = outList$out, file = filename, row.names = FALSE)

	filename <- "TRUEuORFs_20200311.xlsx"
	write.xlsx(x = outList$trueuORFs, file = filename, row.names = FALSE)

	filename <- "FALSEuORFs_20200311.xlsx"
	write.xlsx(x = outList$falseuORFs, file = filename, row.names = FALSE)

	getUPSTcodonSeq <- function(tid, start, end, length, UPSTcodon = 1, noSTART = FALSE){
		if(noSTART == FALSE & UPSTcodon*3+3 > length) return(as.character(NA))
		if(noSTART == TRUE & UPSTcodon*3+3+3 > length) return(as.character(NA))
		codonStart <- end - 2 - UPSTcodon*3
		codonEnd <- end - UPSTcodon*3
		return(as.character(splicedSeqs2[[tid]][codonStart:codonEnd]))
	}

	getDNSTcodonSeq <- function(tid, start, end, length, DNSTcodon = 1, noSTOP = FALSE){
		if(noSTOP == FALSE & DNSTcodon*3+3 > length) return(as.character(NA))
		if(noSTOP == TRUE & DNSTcodon*3+3+3 > length) return(as.character(NA))
		codonStart <- start + DNSTcodon*3
		codonEnd <- start +2 + DNSTcodon*3
		return(as.character(splicedSeqs2[[tid]][codonStart:codonEnd]))
	}

	baseFreq <- function(x){
		if(is.na(x) == TRUE) return(as.character(NA))
		bases <- character(length = 3)
		bases[1] <- substr(x, start = 1, stop = 1)
		bases[2] <- substr(x, start = 2, stop = 2)
		bases[3] <- substr(x, start = 3, stop = 3)
		numA <- length(grep(pattern = "A", bases))
		numC <- length(grep(pattern = "C", bases))
		numG <- length(grep(pattern = "G", bases))
		numT <- length(grep(pattern = "T", bases))
		return(paste(numA, numC, numG, numT, sep = ","))
	}

	checkMorePur <- function(x){
		if(is.na(x) == TRUE) return(as.character(NA))
		baseNum <- as.integer(strsplit(x, ",")[[1]])
		if(baseNum[1] + baseNum[3] >= 2) return(TRUE)	# A, G
		if(baseNum[2] + baseNum[4] >= 2) return(FALSE)	# C, T
	}

	checkMorePyr <- function(x){
		if(is.na(x) == TRUE) return(as.character(NA))
		baseNum <- as.integer(strsplit(x, ",")[[1]])
		if(baseNum[2] + baseNum[4] >= 2) return(TRUE)	# C, T
		if(baseNum[1] + baseNum[3] >= 2) return(FALSE)	# A, G
	}

	checkPurOnly <- function(x){
		if(is.na(x) == TRUE) return(as.character(NA))
		baseNum <- as.integer(strsplit(x, ",")[[1]])
		if(baseNum[1] + baseNum[3] == 3) return(TRUE)	# A, G
		if(baseNum[1] + baseNum[3] != 3) return(FALSE)
	}

	checkPyrOnly <- function(x){
		if(is.na(x) == TRUE) return(as.character(NA))
		baseNum <- as.integer(strsplit(x, ",")[[1]])
		if(baseNum[2] + baseNum[4] == 3) return(TRUE)	# C, T
		if(baseNum[2] + baseNum[4] != 3) return(FALSE)	
	}

	checkCGOnly <- function(x){
		if(is.na(x) == TRUE) return(as.character(NA))
		baseNum <- as.integer(strsplit(x, ",")[[1]])
		if(baseNum[2] + baseNum[3] == 3) return(TRUE)	# C, G
		if(baseNum[2] + baseNum[3] != 3) return(FALSE)
	}

	checkATOnly <- function(x){
		if(is.na(x) == TRUE) return(as.character(NA))
		baseNum <- as.integer(strsplit(x, ",")[[1]])
		if(baseNum[1] + baseNum[4] == 3) return(TRUE)	# A, T
		if(baseNum[1] + baseNum[4] != 3) return(FALSE)	
	}

	checkMoreA <- function(x){
		if(is.na(x) == TRUE) return(as.character(NA))
		baseNum <- as.integer(strsplit(x, ",")[[1]])
		if(baseNum[2] + baseNum[3] >= 1) return(NA)	# C, G
		if(baseNum[1] >= 2) return(TRUE)	# A
		if(baseNum[1] <= 1) return(FALSE)	
	}

	checkMoreT <- function(x){
		if(is.na(x) == TRUE) return(as.character(NA))
		baseNum <- as.integer(strsplit(x, ",")[[1]])
		if(baseNum[2] + baseNum[3] >= 1) return(NA)	# C, G
		if(baseNum[4] >= 2) return(TRUE)	# T
		if(baseNum[4] <= 1) return(FALSE)	
	}

	checkMoreC <- function(x){
		if(is.na(x) == TRUE) return(as.character(NA))
		baseNum <- as.integer(strsplit(x, ",")[[1]])
		if(baseNum[1] + baseNum[4] >= 1) return(NA)	# A, T
		if(baseNum[2] >= 2) return(TRUE)	# C
		if(baseNum[2] <= 1) return(FALSE)	
	}

	checkMoreG <- function(x){
		if(is.na(x) == TRUE) return(as.character(NA))
		baseNum <- as.integer(strsplit(x, ",")[[1]])
		if(baseNum[1] + baseNum[4] >= 1) return(NA)	# A, T
		if(baseNum[3] >= 2) return(TRUE)	# G
		if(baseNum[3] <= 1) return(FALSE)	
	}

	checkMoreCG <- function(x){
		if(is.na(x) == TRUE) return(as.character(NA))
		baseNum <- as.integer(strsplit(x, ",")[[1]])
		if(baseNum[2] + baseNum[3] >= 2) return(TRUE)	# C, G
		if(baseNum[2] + baseNum[3] <= 1) return(FALSE)	
	}

	checkMoreAT <- function(x){
		if(is.na(x) == TRUE) return(as.character(NA))
		baseNum <- as.integer(strsplit(x, ",")[[1]])
		if(baseNum[1] + baseNum[4] >= 2) return(TRUE)	# A, T
		if(baseNum[1] + baseNum[4] <= 1) return(FALSE)	
	}

	typeLists <- list() # 20200422
		
	tempDF <- selectedDF
	for(i in 1:9){
		colName <- paste("UPST", i, sep = "")
		tempDF[[colName]] <- unlist(Map(f = getUPSTcodonSeq, 
					tid = tempDF$transcript_id,
					start = tempDF$start, 
					end <- tempDF$end, 
					length <- tempDF$length, 
					UPSTcodon = i, noSTART = FALSE), use.names = FALSE)
		
		colNameBF <- paste("baseFreq_UPST", i, sep = "")
		tempDF[[colNameBF]] <- unlist(Map(f = baseFreq, 
					x = tempDF[[colName]]), use.names = FALSE)
		
		colNameMPur <- paste("MPur_UPST", i, sep = "")
		tempDF[[colNameMPur]] <- unlist(Map(f = checkMorePur, 
					x = tempDF[[colNameBF]]), use.names = FALSE)
		
		colNameMPyr <- paste("MPyr_UPST", i, sep = "")
		tempDF[[colNameMPyr]] <- unlist(Map(f = checkMorePyr, 
					x = tempDF[[colNameBF]]), use.names = FALSE)

		colNameOPur <- paste("OPur_UPST", i, sep = "")
		tempDF[[colNameOPur]] <- unlist(Map(f = checkPurOnly, 
					x = tempDF[[colNameBF]]), use.names = FALSE)
		colNameOPyr <- paste("OPyr_UPST", i, sep = "")
		tempDF[[colNameOPyr]] <- unlist(Map(f = checkPyrOnly, 
					x = tempDF[[colNameBF]]), use.names = FALSE)

		colNameOCG <- paste("OCG_UPST", i, sep = "")
		tempDF[[colNameOCG]] <- unlist(Map(f = checkCGOnly, 
					x = tempDF[[colNameBF]]), use.names = FALSE)

		colNameOAT <- paste("OAT_UPST", i, sep = "")
		tempDF[[colNameOAT]] <- unlist(Map(f = checkATOnly, 
					x = tempDF[[colNameBF]]), use.names = FALSE)

		colNameMA <- paste("MA_UPST", i, sep = "")
		tempDF[[colNameMA]] <- unlist(Map(f = checkMoreA, 
					x = tempDF[[colNameBF]]), use.names = FALSE)

		colNameMT <- paste("MT_UPST", i, sep = "")
		tempDF[[colNameMT]] <- unlist(Map(f = checkMoreT, 
					x = tempDF[[colNameBF]]), use.names = FALSE)

		colNameMC <- paste("MC_UPST", i, sep = "")
		tempDF[[colNameMC]] <- unlist(Map(f = checkMoreC, 
					x = tempDF[[colNameBF]]), use.names = FALSE)

		colNameMG <- paste("MG_UPST", i, sep = "")
		tempDF[[colNameMG]] <- unlist(Map(f = checkMoreG, 
					x = tempDF[[colNameBF]]), use.names = FALSE)

		colNameMCG <- paste("MCG_UPST", i, sep = "")
		tempDF[[colNameMCG]] <- unlist(Map(f = checkMoreCG, 
					x = tempDF[[colNameBF]]), use.names = FALSE)

		colNameMAT <- paste("MAT_UPST", i, sep = "")
		tempDF[[colNameMAT]] <- unlist(Map(f = checkMoreAT, 
					x = tempDF[[colNameBF]]), use.names = FALSE)

		MPur <- tempDF[tempDF[[colNameMPur]] == TRUE & is.na(tempDF[[colNameMPur]]) == FALSE,]
		MPyr <- tempDF[tempDF[[colNameMPyr]] == TRUE & is.na(tempDF[[colNameMPyr]]) == FALSE,]
		OPur <- tempDF[tempDF[[colNameOPur]] == TRUE & is.na(tempDF[[colNameOPur]]) == FALSE,]
		OPyr <- tempDF[tempDF[[colNameOPyr]] == TRUE & is.na(tempDF[[colNameOPyr]]) == FALSE,]
		OCG <- tempDF[tempDF[[colNameOCG]] == TRUE & is.na(tempDF[[colNameOCG]]) == FALSE,]
		OAT <- tempDF[tempDF[[colNameOAT]] == TRUE & is.na(tempDF[[colNameOAT]]) == FALSE,]
		MA <- tempDF[tempDF[[colNameMA]] == TRUE & is.na(tempDF[[colNameMA]]) == FALSE,]
		MT <- tempDF[tempDF[[colNameMT]] == TRUE & is.na(tempDF[[colNameMT]]) == FALSE,]
		MC <- tempDF[tempDF[[colNameMC]] == TRUE & is.na(tempDF[[colNameMC]]) == FALSE,]
		MG <- tempDF[tempDF[[colNameMG]] == TRUE & is.na(tempDF[[colNameMG]]) == FALSE,]
		MCG <- tempDF[tempDF[[colNameMCG]] == TRUE & is.na(tempDF[[colNameMCG]]) == FALSE,]
		MAT <- tempDF[tempDF[[colNameMAT]] == TRUE & is.na(tempDF[[colNameMAT]]) == FALSE,]


		filename <- paste("Bases_RPKMthresh", RPKM.threshold, "_UPST", i, ".pdf", sep = "")
		pdf(file = filename, width = 6*2, height = 4*2)
		par(mfrow = c(2, 2))

		# types <- c("ratioWT", "ratiogcn2", "ratioWo3AT", "ratioWith3AT", 
		# 		"WTwo3AT", "WTwith3AT", "gcn2wo3AT", "gcn2with3AT", 
		# 		"WTwo3ATvsCDS", "WTwith3ATvsCDS", "gcn2wo3ATvsCDS", "gcn2with3ATvsCDS")

		types <- c("WTwo3ATvsCDS")
		typeList <- list()
		for(t in 1:length(types)){
			type <- types[t]
			typeList$MPur <- log2(MPur[[type]])
			typeList$MPyr <- log2(MPyr[[type]])
			typeList$OPur <- log2(OPur[[type]])
			typeList$OPyr <- log2(OPyr[[type]])
			typeList$OCG <- log2(OCG[[type]])
			typeList$OAT <- log2(OAT[[type]])
			typeList$MA <- log2(MA[[type]])
			typeList$MT <- log2(MT[[type]])
			typeList$MC <- log2(MC[[type]])
			typeList$MG <- log2(MG[[type]])
			typeList$MCG <- log2(MCG[[type]])
			typeList$MAT <- log2(MAT[[type]])

			codonNumbers <- c(nrow(MPur), nrow(MPyr), nrow(OPur), nrow(OPyr), 
					nrow(OCG), nrow(OAT), nrow(MA), nrow(MT), nrow(MC), nrow(MG), nrow(MCG), nrow(MAT))
			sub <- paste("n=", paste(codonNumbers, collapse = ","), sep = "")
			main <- paste(type, " RPKMthresh", RPKM.threshold, " UPST", i, sep = "")
			ylim <- c(min(unlist(typeList), na.rm = TRUE), max(unlist(typeList), na.rm = TRUE))
			ylab <- paste(type, " (log2)", sep = "")
			boxplot(typeList, outline = FALSE, main = main, sub = sub, 
				ylab = ylab, ylim = ylim, cex.axis = 0.8)
			beeswarm(typeList, corral = "wrap", col = rgb(1, 0, 0, alpha = 0.5), pch = 20, add = TRUE)

			codonNumbers2 <- codonNumbers[c(1,2,5,6)]
			sub <- paste("n=", paste(codonNumbers2, collapse = ","), sep = "")

			boxplot(typeList[c(1,2,5,6)], outline = FALSE, main = main, sub = sub, 
				ylab = ylab, ylim = ylim, cex.axis = 0.8)
			beeswarm(typeList[c(1,2,5,6)], corral = "wrap", col = rgb(1, 0, 0, alpha = 0.5), pch = 20, add = TRUE)

			pvalues <- numeric(length = 6)
			pvalues[1:6] <- NA
			names(pvalues) <- c("MPur vs MPyr", "OPur vs OPyr", "OCG vs OAT", "MA vs MT", "MC vs MG", "MCG vs MAT")
			if(codonNumbers[1] > 5 & codonNumbers[2] > 5){
				pvalues[1] <- t.test(typeList[[1]], typeList[[2]])$p.value
			}
			if(codonNumbers[3] > 5 & codonNumbers[4] > 5){
				pvalues[2] <- t.test(typeList[[3]], typeList[[4]])$p.value
			}
			if(codonNumbers[5] > 5 & codonNumbers[6] > 5){
				pvalues[3] <- t.test(typeList[[5]], typeList[[6]])$p.value
			}
			if(codonNumbers[7] > 5 & codonNumbers[8] > 5){
				pvalues[4] <- t.test(typeList[[7]], typeList[[8]])$p.value
			}
			if(codonNumbers[9] > 5 & codonNumbers[10] > 5){
				pvalues[5] <- t.test(typeList[[9]], typeList[[10]])$p.value
			}
			if(codonNumbers[11] > 5 & codonNumbers[12] > 5){
				pvalues[6] <- t.test(typeList[[11]], typeList[[12]])$p.value
			}

			ymax <- -log10(min(pvalues, na.rm = TRUE)) + -log10(min(pvalues, na.rm = TRUE))*0.1
			if(ymax < 3) ymax <- 3.2
			ylab <- "-log10(p-value)"
			y_coords <- barplot(-log10(pvalues), main = main, ylab = ylab, ylim = c(0, ymax), cex.names = 0.6)
			text(y_coords, -log10(pvalues), labels = formatC(pvalues, format = "e", digits = 1), pos = 3)
			abline(h = -log10(0.05), col = "orange")
			abline(h = -log10(0.01), col = "red")
			abline(h = -log10(0.001), col = "purple")

			pvalues2 <- pvalues[c(1,3)]
			y_coords <- barplot(-log10(pvalues2), main = main, ylab = ylab, ylim = c(0, ymax), cex.names = 0.6)
			text(y_coords, -log10(pvalues2), labels = formatC(pvalues2, format = "e", digits = 1), pos = 3)
			abline(h = -log10(0.05), col = "orange")
			abline(h = -log10(0.01), col = "red")
			abline(h = -log10(0.001), col = "purple")

			if(i == 1)	typeLists[[type]] <- list()
			if(i == 1)	typeLists[[type]][["MPur_vs_MPyr"]] <- list()
			listName <- paste("UPST", i, sep = "")
			typeLists[[type]][["MPur_vs_MPyr"]][[listName]] <- list()
			typeLists[[type]][["MPur_vs_MPyr"]][[listName]][["typeList"]] <- typeList[c(1,2)]
			typeLists[[type]][["MPur_vs_MPyr"]][[listName]][["codonNumbers"]] <- codonNumbers[c(1,2)]
			typeLists[[type]][["MPur_vs_MPyr"]][[listName]][["pvalue"]] <- pvalues[c(1)]
			typeLists[[type]][["OCG_vs_OAT"]][[listName]] <- list()
			typeLists[[type]][["OCG_vs_OAT"]][[listName]][["typeList"]] <- typeList[c(5,6)]
			typeLists[[type]][["OCG_vs_OAT"]][[listName]][["codonNumbers"]] <- codonNumbers[c(5,6)]
			typeLists[[type]][["OCG_vs_OAT"]][[listName]][["pvalue"]] <- pvalues[c(3)]
		}
		dev.off()
	} # i-loop
	
	filename <- paste("Bases_RPKMthresh", RPKM.threshold, "_UPST1-4_selected_.pdf", sep = "")
	pdf(file = filename, width = 6*2, height = 4*2)
	par(mfrow = c(2, 2))
	for(PAIR in c("MPur_vs_MPyr", "OCG_vs_OAT")){	

		type <- "WTwo3ATvsCDS"
		codonNumbers <- integer(length = 8)
		pvalues <- numeric(length = 4)
		DATA <- list()
		for(i in 1:4){
			listName <- paste("UPST", i, sep = "")
			codonNumbers[c((i-1)*2+1, (i-1)*2+2)] <- typeLists[[type]][[PAIR]][[listName]][["codonNumbers"]]
			pvalues[i] <- typeLists[[type]][[PAIR]][[listName]][["pvalue"]]
			if(PAIR == "MPur_vs_MPyr"){
				listName2 <- paste("MPur", i, sep = "")
				DATA[[listName2]] <- typeLists[[type]][[PAIR]][[listName]][["typeList"]]$MPur
				listName3 <- paste("MPyr", i, sep = "")
				DATA[[listName3]] <- typeLists[[type]][[PAIR]][[listName]][["typeList"]]$MPyr
			}
			if(PAIR == "OCG_vs_OAT"){
				listName2 <- paste("OCG", i, sep = "")
				DATA[[listName2]] <- typeLists[[type]][[PAIR]][[listName]][["typeList"]]$OCG
				listName3 <- paste("OAT", i, sep = "")
				DATA[[listName3]] <- typeLists[[type]][[PAIR]][[listName]][["typeList"]]$OAT
			}
		}

		sub <- paste("n=", paste(codonNumbers, collapse = ","), sep = "")
		main <- paste(type, " RPKMthresh", RPKM.threshold, " UPST1-4 ", PAIR, sep = "")
		ylim <- c(min(unlist(DATA), na.rm = TRUE), max(unlist(DATA), na.rm = TRUE))
		ylab <- paste(type, " (log2)", sep = "")
		boxplot(DATA, outline = FALSE, main = main, sub = sub, 
			ylab = ylab, ylim = ylim, cex.axis = 0.8)
		beeswarm(DATA, corral = "wrap", col = rgb(1, 0, 0, alpha = 0.5), pch = 20, add = TRUE)

		ymax <- -log10(min(pvalues, na.rm = TRUE)) + -log10(min(pvalues, na.rm = TRUE))*0.1
		if(ymax < 3) ymax <- 3.2
		ylab <- "-log10(p-value)"
		y_coords <- barplot(-log10(pvalues), main = main, ylab = ylab, ylim = c(0, ymax), cex.names = 0.6)
		text(y_coords, -log10(pvalues), labels = formatC(pvalues, format = "e", digits = 1), pos = 3)
		abline(h = -log10(0.05), col = "orange")
		abline(h = -log10(0.01), col = "red")
		abline(h = -log10(0.001), col = "purple")
	} # PAIR-loop
	dev.off()

	# noSTART
	for(i in 1:9){
		colName <- paste("UPST", i, "noSTART", sep = "")
		tempDF[[colName]] <- unlist(mcMap(f = getUPSTcodonSeq, 
					tid = tempDF$transcript_id,
					start = tempDF$start, 
					end <- tempDF$end, 
					length <- tempDF$length, 
					UPSTcodon = i, noSTART = TRUE, mc.cores = 4), use.names = FALSE)
		
		colNameBF <- paste("baseFreq_UPST", i, "noSTART", sep = "")
		tempDF[[colNameBF]] <- unlist(Map(f = baseFreq, 
					x = tempDF[[colName]]), use.names = FALSE)
		
		colNameMPur <- paste("MPur_UPST", i, "noSTART", sep = "")
		tempDF[[colNameMPur]] <- unlist(Map(f = checkMorePur, 
					x = tempDF[[colNameBF]]), use.names = FALSE)
		
		colNameMPyr <- paste("MPyr_UPST", i, "noSTART", sep = "")
		tempDF[[colNameMPyr]] <- unlist(Map(f = checkMorePyr, 
					x = tempDF[[colNameBF]]), use.names = FALSE)

		colNameOPur <- paste("OPur_UPST", i, "noSTART", sep = "")
		tempDF[[colNameOPur]] <- unlist(Map(f = checkPurOnly, 
					x = tempDF[[colNameBF]]), use.names = FALSE)
		colNameOPyr <- paste("OPyr_UPST", i, "noSTART", sep = "")
		tempDF[[colNameOPyr]] <- unlist(Map(f = checkPyrOnly, 
					x = tempDF[[colNameBF]]), use.names = FALSE)

		colNameOCG <- paste("OCG_UPST", i, "noSTART", sep = "")
		tempDF[[colNameOCG]] <- unlist(Map(f = checkCGOnly, 
					x = tempDF[[colNameBF]]), use.names = FALSE)

		colNameOAT <- paste("OAT_UPST", i, "noSTART", sep = "")
		tempDF[[colNameOAT]] <- unlist(Map(f = checkATOnly, 
					x = tempDF[[colNameBF]]), use.names = FALSE)

		colNameMA <- paste("MA_UPST", i, "noSTART", sep = "")
		tempDF[[colNameMA]] <- unlist(Map(f = checkMoreA, 
					x = tempDF[[colNameBF]]), use.names = FALSE)

		colNameMT <- paste("MT_UPST", i, "noSTART", sep = "")
		tempDF[[colNameMT]] <- unlist(Map(f = checkMoreT, 
					x = tempDF[[colNameBF]]), use.names = FALSE)

		colNameMC <- paste("MC_UPST", i, "noSTART", sep = "")
		tempDF[[colNameMC]] <- unlist(Map(f = checkMoreC, 
					x = tempDF[[colNameBF]]), use.names = FALSE)

		colNameMG <- paste("MG_UPST", i, "noSTART", sep = "")
		tempDF[[colNameMG]] <- unlist(Map(f = checkMoreG, 
					x = tempDF[[colNameBF]]), use.names = FALSE)

		colNameMCG <- paste("MCG_UPST", i, "noSTART", sep = "")
		tempDF[[colNameMCG]] <- unlist(Map(f = checkMoreCG, 
					x = tempDF[[colNameBF]]), use.names = FALSE)

		colNameMAT <- paste("MAT_UPST", i, "noSTART", sep = "")
		tempDF[[colNameMAT]] <- unlist(Map(f = checkMoreAT, 
					x = tempDF[[colNameBF]]), use.names = FALSE)

		MPur <- tempDF[tempDF[[colNameMPur]] == TRUE & is.na(tempDF[[colNameMPur]]) == FALSE,]
		MPyr <- tempDF[tempDF[[colNameMPyr]] == TRUE & is.na(tempDF[[colNameMPyr]]) == FALSE,]
		OPur <- tempDF[tempDF[[colNameOPur]] == TRUE & is.na(tempDF[[colNameOPur]]) == FALSE,]
		OPyr <- tempDF[tempDF[[colNameOPyr]] == TRUE & is.na(tempDF[[colNameOPyr]]) == FALSE,]
		OCG <- tempDF[tempDF[[colNameOCG]] == TRUE & is.na(tempDF[[colNameOCG]]) == FALSE,]
		OAT <- tempDF[tempDF[[colNameOAT]] == TRUE & is.na(tempDF[[colNameOAT]]) == FALSE,]
		MA <- tempDF[tempDF[[colNameMA]] == TRUE & is.na(tempDF[[colNameMA]]) == FALSE,]
		MT <- tempDF[tempDF[[colNameMT]] == TRUE & is.na(tempDF[[colNameMT]]) == FALSE,]
		MC <- tempDF[tempDF[[colNameMC]] == TRUE & is.na(tempDF[[colNameMC]]) == FALSE,]
		MG <- tempDF[tempDF[[colNameMG]] == TRUE & is.na(tempDF[[colNameMG]]) == FALSE,]
		MCG <- tempDF[tempDF[[colNameMCG]] == TRUE & is.na(tempDF[[colNameMCG]]) == FALSE,]
		MAT <- tempDF[tempDF[[colNameMAT]] == TRUE & is.na(tempDF[[colNameMAT]]) == FALSE,]


		filename <- paste("Bases_RPKMthresh", RPKM.threshold, "_UPST", i, "_noSTART.pdf", sep = "")
		pdf(file = filename, width = 6*2, height = 4*2)
		par(mfrow = c(2, 2))

		# types <- c("ratioWT", "ratiogcn2", "ratioWo3AT", "ratioWith3AT", 
		# 		"WTwo3AT", "WTwith3AT", "gcn2wo3AT", "gcn2with3AT", 
		# 		"WTwo3ATvsCDS", "WTwith3ATvsCDS", "gcn2wo3ATvsCDS", "gcn2with3ATvsCDS")

		types <- c("WTwo3ATvsCDS")
		typeList <- list()
		for(t in 1:length(types)){
			type <- types[t]
			typeList$MPur <- log2(MPur[[type]])
			typeList$MPyr <- log2(MPyr[[type]])
			typeList$OPur <- log2(OPur[[type]])
			typeList$OPyr <- log2(OPyr[[type]])
			typeList$OCG <- log2(OCG[[type]])
			typeList$OAT <- log2(OAT[[type]])
			typeList$MA <- log2(MA[[type]])
			typeList$MT <- log2(MT[[type]])
			typeList$MC <- log2(MC[[type]])
			typeList$MG <- log2(MG[[type]])
			typeList$MCG <- log2(MCG[[type]])
			typeList$MAT <- log2(MAT[[type]])

			codonNumbers <- c(nrow(MPur), nrow(MPyr), nrow(OPur), nrow(OPyr), 
					nrow(OCG), nrow(OAT), nrow(MA), nrow(MT), 
					nrow(MC), nrow(MG), nrow(MCG), nrow(MAT))
			sub <- paste("n=", paste(codonNumbers, collapse = ","), sep = "")
			main <- paste(type, " RPKMthresh", RPKM.threshold, " UPST", i, " noSTART", sep = "")
			ylim <- c(min(unlist(typeList), na.rm = TRUE), max(unlist(typeList), na.rm = TRUE))
			ylab <- paste(type, " (log2)", sep = "")
			boxplot(typeList, outline = FALSE, main = main, sub = sub, 
				ylab = ylab, ylim = ylim, cex.axis = 0.8)
			beeswarm(typeList, corral = "wrap", col = rgb(1, 0, 0, alpha = 0.5), pch = 20, add = TRUE)

			pvalues <- numeric(length = 6)
			pvalues[1:6] <- NA
			names(pvalues) <- c("MPur vs MPyr", "OPur vs OPyr", "OCG vs OAT", 
						"MA vs MT", "MC vs MG", "MCG vs MAT")
			if(codonNumbers[1] > 5 & codonNumbers[2] > 5){
				pvalues[1] <- t.test(typeList[[1]], typeList[[2]])$p.value
			}
			if(codonNumbers[3] > 5 & codonNumbers[4] > 5){
				pvalues[2] <- t.test(typeList[[3]], typeList[[4]])$p.value
			}
			if(codonNumbers[5] > 5 & codonNumbers[6] > 5){
				pvalues[3] <- t.test(typeList[[5]], typeList[[6]])$p.value
			}
			if(codonNumbers[7] > 5 & codonNumbers[8] > 5){
				pvalues[4] <- t.test(typeList[[7]], typeList[[8]])$p.value
			}
			if(codonNumbers[9] > 5 & codonNumbers[10] > 5){
				pvalues[5] <- t.test(typeList[[9]], typeList[[10]])$p.value
			}
			if(codonNumbers[11] > 5 & codonNumbers[12] > 5){
				pvalues[6] <- t.test(typeList[[11]], typeList[[12]])$p.value
			}

			ymax <- -log10(min(pvalues, na.rm = TRUE)) + -log10(min(pvalues, na.rm = TRUE))*0.1
			if(ymax < 3) ymax <- 3.2
			ylab <- "-log10(p-value)"
			y_coords <- barplot(-log10(pvalues), main = main, ylab = ylab, 
					ylim = c(0, ymax), cex.names = 0.6)
			text(y_coords, -log10(pvalues), labels = formatC(pvalues, format = "e", digits = 1), pos = 3)
			abline(h = -log10(0.05), col = "orange")
			abline(h = -log10(0.01), col = "red")
			abline(h = -log10(0.001), col = "purple")
		}
		dev.off()
	} # i-loop

	log2plotLen_tempDF <- function(type1, type2){
		cat("log2plotLen_tempDF, type1:", type1, ", type2:", type2, "\n", sep = "")
		
			out1o <- tempDF[0,]
			out2o <- tempDF[0,]
		for(i in 1:max(tempDF$length)){
			# cat(i, ",", sep = "")
			uORFlength <- i*3+3
			temp <- subset(tempDF, length == uORFlength)
			# ratio <- temp[[type2]]/temp[[type1]]
			ratio <- log2(temp[[type2]]/temp[[type1]])
			if(length(which(ratio > log2(1.5)))>0){
				# out1 <- temp[which(ratio > log2(1.5)), c(1:9, 25:38, 47:50)]
				out1 <- temp[which(ratio > log2(1.5)), ]
				# if(nrow(out1) == 0) cat(nrow(out1))
				out1$type1 <- type1
				out1$type2 <- type2
				out1$upFOLD <- type2
				out1$dnFOLD <- "-"
				# out1$ratiolog2 <- log2(ratio[which(ratio > log2(1.5))])
				out1$ratiolog2 <- ratio[which(ratio > log2(1.5))]
				# out1 <- out1[c((ncol(out1)-4):ncol(out1), 1:(ncol(out1)-5))]
				if(i == 1) out1o <- out1
				if(i != 1) out1o <- rbind(out1o, out1)
			}
			if(length(which(ratio < -log2(1.5)))>0){
				# out2 <- temp[which(ratio < -log2(1.5)), c(1:9, 25:38, 47:50)]
				out2 <- temp[which(ratio < -log2(1.5)), ]
				# if(nrow(out2) == 0) cat(nrow(out2))
				out2$type1 <- type1
				out2$type2 <- type2
				out2$upFOLD <- "-"
				out2$dnFOLD <- type1
				# out2$ratiolog2 <- log2(ratio[which(ratio < -log2(1.5))])
				out2$ratiolog2 <- ratio[which(ratio < -log2(1.5))]
				# out2 <- out2[c((ncol(out2)-4):ncol(out2), 1:(ncol(out2)-5))]
				if(i == 1) out2o <- out2
				if(i != 1) out2o <- rbind(out2o, out2)
			}
		}

		out1o <- subset(out1o, is.na(transcript_id) == FALSE)
		out2o <- subset(out2o, is.na(transcript_id) == FALSE)
		if(nrow(out1o) > 0 & nrow(out2o) > 0){
			out <- rbind(out1o, out2o)
			return(out)
		}
		if(nrow(out1o) > 0){
			return(out1o)
		}
		if(nrow(out2o) > 0){
			return(out2o)
		}
		if(nrow(out1o) == 0 & nrow(out2o) == 0){
			return(out2o)
		}
	}

	out <- log2plotLen_tempDF(type1 = "ratioWT", type2 = "ratiogcn2")

	filename <- "Diff2FOLDBase_20200311.xlsx"
	write.xlsx(x = out, file = filename, row.names = FALSE)

	filename <- "Base_20200311.xlsx"
	write.xlsx(x = tempDF, file = filename, row.names = FALSE)
	
	typeLists <- list() # 20200422

	tempDF2 <- selectedDF
	for(i in 1:1){			
		colName <- paste("UPST", i, sep = "")
		tempDF2[[colName]] <- unlist(Map(f = getUPSTcodonSeq, 
					tid = tempDF2$transcript_id,
					start = tempDF2$start, 
					end <- tempDF2$end, 
					length <- tempDF2$length, 
					UPSTcodon = i, noSTART = FALSE), use.names = FALSE)
		
		colNameBF <- paste("baseFreq_UPST", i, sep = "")
		tempDF2[[colNameBF]] <- unlist(Map(f = baseFreq, 
					x = tempDF2[[colName]]), use.names = FALSE)
		
		colNameMPur <- paste("MPur_UPST", i, sep = "")
		tempDF2[[colNameMPur]] <- unlist(Map(f = checkMorePur, 
					x = tempDF2[[colNameBF]]), use.names = FALSE)
		
		colNameMPyr <- paste("MPyr_UPST", i, sep = "")
		tempDF2[[colNameMPyr]] <- unlist(Map(f = checkMorePyr, 
					x = tempDF2[[colNameBF]]), use.names = FALSE)

		colNameOPur <- paste("OPur_UPST", i, sep = "")
		tempDF2[[colNameOPur]] <- unlist(Map(f = checkPurOnly, 
					x = tempDF2[[colNameBF]]), use.names = FALSE)
		colNameOPyr <- paste("OPyr_UPST", i, sep = "")
		tempDF2[[colNameOPyr]] <- unlist(Map(f = checkPyrOnly, 
					x = tempDF2[[colNameBF]]), use.names = FALSE)

		colNameOCG <- paste("OCG_UPST", i, sep = "")
		tempDF2[[colNameOCG]] <- unlist(Map(f = checkCGOnly, 
					x = tempDF2[[colNameBF]]), use.names = FALSE)

		colNameOAT <- paste("OAT_UPST", i, sep = "")
		tempDF2[[colNameOAT]] <- unlist(Map(f = checkATOnly, 
					x = tempDF2[[colNameBF]]), use.names = FALSE)

		colNameMA <- paste("MA_UPST", i, sep = "")
		tempDF2[[colNameMA]] <- unlist(Map(f = checkMoreA, 
					x = tempDF2[[colNameBF]]), use.names = FALSE)

		colNameMT <- paste("MT_UPST", i, sep = "")
		tempDF2[[colNameMT]] <- unlist(Map(f = checkMoreT, 
					x = tempDF2[[colNameBF]]), use.names = FALSE)

		colNameMC <- paste("MC_UPST", i, sep = "")
		tempDF2[[colNameMC]] <- unlist(Map(f = checkMoreC, 
					x = tempDF2[[colNameBF]]), use.names = FALSE)

		colNameMG <- paste("MG_UPST", i, sep = "")
		tempDF2[[colNameMG]] <- unlist(Map(f = checkMoreG, 
					x = tempDF2[[colNameBF]]), use.names = FALSE)

		colNameMCG <- paste("MCG_UPST", i, sep = "")
		tempDF2[[colNameMCG]] <- unlist(Map(f = checkMoreCG, 
					x = tempDF2[[colNameBF]]), use.names = FALSE)

		colNameMAT <- paste("MAT_UPST", i, sep = "")
		tempDF2[[colNameMAT]] <- unlist(Map(f = checkMoreAT, 
					x = tempDF2[[colNameBF]]), use.names = FALSE)

		MPur <- tempDF2[tempDF2[[colNameMPur]] == TRUE & is.na(tempDF2[[colNameMPur]]) == FALSE,]
		MPyr <- tempDF2[tempDF2[[colNameMPyr]] == TRUE & is.na(tempDF2[[colNameMPyr]]) == FALSE,]
		OPur <- tempDF2[tempDF2[[colNameOPur]] == TRUE & is.na(tempDF2[[colNameOPur]]) == FALSE,]
		OPyr <- tempDF2[tempDF2[[colNameOPyr]] == TRUE & is.na(tempDF2[[colNameOPyr]]) == FALSE,]
		OCG <- tempDF2[tempDF2[[colNameOCG]] == TRUE & is.na(tempDF2[[colNameOCG]]) == FALSE,]
		OAT <- tempDF2[tempDF2[[colNameOAT]] == TRUE & is.na(tempDF2[[colNameOAT]]) == FALSE,]
		MA <- tempDF2[tempDF2[[colNameMA]] == TRUE & is.na(tempDF2[[colNameMA]]) == FALSE,]
		MT <- tempDF2[tempDF2[[colNameMT]] == TRUE & is.na(tempDF2[[colNameMT]]) == FALSE,]
		MC <- tempDF2[tempDF2[[colNameMC]] == TRUE & is.na(tempDF2[[colNameMC]]) == FALSE,]
		MG <- tempDF2[tempDF2[[colNameMG]] == TRUE & is.na(tempDF2[[colNameMG]]) == FALSE,]
		MCG <- tempDF2[tempDF2[[colNameMCG]] == TRUE & is.na(tempDF2[[colNameMCG]]) == FALSE,]
		MAT <- tempDF2[tempDF2[[colNameMAT]] == TRUE & is.na(tempDF2[[colNameMAT]]) == FALSE,]

		filename <- paste("Bases_RPKMthresh", RPKM.threshold, "_length_UPST", i, ".pdf", sep = "")
		pdf(file = filename, width = 6*2, height = 4*2)
		par(mfrow = c(2, 2))

		# types <- c("ratioWT", "ratiogcn2", "ratioWo3AT", "ratioWith3AT", 
		# 		"WTwo3AT", "WTwith3AT", "gcn2wo3AT", "gcn2with3AT", 
		# 		"WTwo3ATvsCDS", "WTwith3ATvsCDS", "gcn2wo3ATvsCDS", "gcn2with3ATvsCDS")

		types <- c("WTwo3ATvsCDS")
		typeList <- list()
		for(t in 1:length(types)){
		# for(l in c(9, 12, 15, 18, 21, 24, 27)){	
		for(l in c(9, 12, 15, 18, 21, 24, 27, 30, 33, 36, 39, 42, 45, 48, 51, 54, 57, 60)){	
			type <- types[t]
			typeList$MPur <- log2(MPur[[type]][MPur$length == l])	
			typeList$MPyr <- log2(MPyr[[type]][MPyr$length == l])	
			typeList$OPur <- log2(OPur[[type]][OPur$length == l])	
			typeList$OPyr <- log2(OPyr[[type]][OPyr$length == l])	
			typeList$OCG <- log2(OCG[[type]][OCG$length == l])	
			typeList$OAT <- log2(OAT[[type]][OAT$length == l])	
			typeList$MA <- log2(MA[[type]][MA$length == l])		
			typeList$MT <- log2(MT[[type]][MT$length == l])		
			typeList$MC <- log2(MC[[type]][MC$length == l])		
			typeList$MG <- log2(MG[[type]][MG$length == l])		
			typeList$MCG <- log2(MCG[[type]][MCG$length == l])	
			typeList$MAT <- log2(MAT[[type]][MAT$length == l])	

			codonNumbers <- c(length(which(is.na(typeList$MPur) == FALSE)), 
					length(which(is.na(typeList$MPyr) == FALSE)), 
					length(which(is.na(typeList$OPur) == FALSE)), 
					length(which(is.na(typeList$OPyr) == FALSE)), 
					length(which(is.na(typeList$OCG) == FALSE)), 
					length(which(is.na(typeList$OAT) == FALSE)), 
					length(which(is.na(typeList$MA) == FALSE)), 
					length(which(is.na(typeList$MT) == FALSE)), 
					length(which(is.na(typeList$MC) == FALSE)), 
					length(which(is.na(typeList$MG) == FALSE)), 
					length(which(is.na(typeList$MCG) == FALSE)), 
					length(which(is.na(typeList$MAT) == FALSE)))
			sub <- paste("n=", paste(codonNumbers, collapse = ","), sep = "")
			main <- paste(type, " RPKMthresh", RPKM.threshold, " length=", l, " UPST", i, sep = "")
			ylim <- c(min(unlist(typeList), na.rm = TRUE), max(unlist(typeList), na.rm = TRUE))
			ylab <- paste(type, " (log2)", sep = "")
			boxplot(typeList, outline = FALSE, main = main, sub = sub, 
				ylab = ylab, ylim = ylim, cex.axis = 0.8)
			beeswarm(typeList, corral = "wrap", col = rgb(1, 0, 0, alpha = 0.5), pch = 20, add = TRUE)

			pvalues <- numeric(length = 6)
			pvalues[1:6] <- NA
			names(pvalues) <- c("MPur vs MPyr", "OPur vs OPyr", "OCG vs OAT", 
						"MA vs MT", "MC vs MG", "MCG vs MAT")
			if(codonNumbers[1] > 5 & codonNumbers[2] > 5){
				pvalues[1] <- t.test(typeList[[1]], typeList[[2]])$p.value
			}
			if(codonNumbers[3] > 5 & codonNumbers[4] > 5){
				pvalues[2] <- t.test(typeList[[3]], typeList[[4]])$p.value
			}
			if(codonNumbers[5] > 5 & codonNumbers[6] > 5){
				pvalues[3] <- t.test(typeList[[5]], typeList[[6]])$p.value
			}
			if(codonNumbers[7] > 5 & codonNumbers[8] > 5){
				pvalues[4] <- t.test(typeList[[7]], typeList[[8]])$p.value
			}
			if(codonNumbers[9] > 5 & codonNumbers[10] > 5){
				pvalues[5] <- t.test(typeList[[9]], typeList[[10]])$p.value
			}
			if(codonNumbers[11] > 5 & codonNumbers[12] > 5){
				pvalues[6] <- t.test(typeList[[11]], typeList[[12]])$p.value
			}

			if(sum(pvalues, na.rm = TRUE)>0){
				ymax <- -log10(min(pvalues, na.rm = TRUE)) + -log10(min(pvalues, na.rm = TRUE))*0.1
			}
			if(ymax < 3) ymax <- 3.2
			ylab <- "-log10(p-value)"
			y_coords <- barplot(-log10(pvalues), main = main, 
					ylab = ylab, ylim = c(0, ymax), cex.names = 0.6)
			text(y_coords, -log10(pvalues), labels = formatC(pvalues, format = "e", digits = 1), pos = 3)
			abline(h = -log10(0.05), col = "orange")
			abline(h = -log10(0.01), col = "red")
			abline(h = -log10(0.001), col = "purple")
		} # l-loop
		} # t-loop
		dev.off()

		filename <- paste("Bases_RPKMthresh", RPKM.threshold, "_lengthCumulative_UPST", i, ".pdf", sep = "")
		pdf(file = filename, width = 6*2, height = 4*2)
		par(mfrow = c(2, 2))

		# types <- c("ratioWT", "ratiogcn2", "ratioWo3AT", "ratioWith3AT", 
		# 		"WTwo3AT", "WTwith3AT", "gcn2wo3AT", "gcn2with3AT", 
		# 		"WTwo3ATvsCDS", "WTwith3ATvsCDS", "gcn2wo3ATvsCDS", "gcn2with3ATvsCDS")

		types <- c("WTwo3ATvsCDS")
		typeList <- list()
		for(t in 1:length(types)){
			type <- types[t]
		for(l in c(9, 12, 15, 18, 21, 24, 27, 30, 33, 36, 39, 42, 45, 48, 51, 54, 57, 60)){
			if(l == 9){
				typeList$MPur <- log2(MPur[[type]][MPur$length == l])
				typeList$MPyr <- log2(MPyr[[type]][MPyr$length == l])
				typeList$OPur <- log2(OPur[[type]][OPur$length == l])
				typeList$OPyr <- log2(OPyr[[type]][OPyr$length == l])
				typeList$OCG <- log2(OCG[[type]][OCG$length == l])
				typeList$OAT <- log2(OAT[[type]][OAT$length == l])
				typeList$MA <- log2(MA[[type]][MA$length == l])	
				typeList$MT <- log2(MT[[type]][MT$length == l])	
				typeList$MC <- log2(MC[[type]][MC$length == l])	
				typeList$MG <- log2(MG[[type]][MG$length == l])	
				typeList$MCG <- log2(MCG[[type]][MCG$length == l])
				typeList$MAT <- log2(MAT[[type]][MAT$length == l])
			}
			if(l > 9){
				typeList$MPur <- c(typeList$MPur, log2(MPur[[type]][MPur$length == l]))
				typeList$MPyr <- c(typeList$MPyr, log2(MPyr[[type]][MPyr$length == l]))
				typeList$OPur <- c(typeList$OPur, log2(OPur[[type]][OPur$length == l]))
				typeList$OPyr <- c(typeList$OPyr, log2(OPyr[[type]][OPyr$length == l]))
				typeList$OCG <- c(typeList$OCG, log2(OCG[[type]][OCG$length == l]))
				typeList$OAT <- c(typeList$OAT, log2(OAT[[type]][OAT$length == l]))
				typeList$MA <- c(typeList$MA, log2(MA[[type]][MA$length == l]))
				typeList$MT <- c(typeList$MT, log2(MT[[type]][MT$length == l]))
				typeList$MC <- c(typeList$MC, log2(MC[[type]][MC$length == l]))
				typeList$MG <- c(typeList$MG, log2(MG[[type]][MG$length == l]))
				typeList$MCG <- c(typeList$MCG, log2(MCG[[type]][MCG$length == l]))
				typeList$MAT <- c(typeList$MAT, log2(MAT[[type]][MAT$length == l]))
			}

			codonNumbers <- c(length(which(is.na(typeList$MPur) == FALSE)), 
					length(which(is.na(typeList$MPyr) == FALSE)), 
					length(which(is.na(typeList$OPur) == FALSE)), 
					length(which(is.na(typeList$OPyr) == FALSE)), 
					length(which(is.na(typeList$OCG) == FALSE)), 
					length(which(is.na(typeList$OAT) == FALSE)), 
					length(which(is.na(typeList$MA) == FALSE)), 
					length(which(is.na(typeList$MT) == FALSE)), 
					length(which(is.na(typeList$MC) == FALSE)), 
					length(which(is.na(typeList$MG) == FALSE)), 
					length(which(is.na(typeList$MCG) == FALSE)), 
					length(which(is.na(typeList$MAT) == FALSE)))
			sub <- paste("n=", paste(codonNumbers, collapse = ","), sep = "")
			main <- paste(type, " RPKMthresh", RPKM.threshold, " length=", l, " UPST", i, sep = "")
			ylim <- c(min(unlist(typeList), na.rm = TRUE), max(unlist(typeList), na.rm = TRUE))
			ylab <- paste(type, " (log2)", sep = "")
			boxplot(typeList, outline = FALSE, main = main, sub = sub, 
				ylab = ylab, ylim = ylim, cex.axis = 0.8)
			beeswarm(typeList, corral = "wrap", col = rgb(1, 0, 0, alpha = 0.3), pch = 20, add = TRUE)

			pvalues <- numeric(length = 6)
			pvalues[1:6] <- NA
			names(pvalues) <- c("MPur vs MPyr", "OPur vs OPyr", "OCG vs OAT", 
						"MA vs MT", "MC vs MG", "MCG vs MAT")
			if(codonNumbers[1] > 5 & codonNumbers[2] > 5){
				pvalues[1] <- t.test(typeList[[1]], typeList[[2]])$p.value
			}
			if(codonNumbers[3] > 5 & codonNumbers[4] > 5){
				pvalues[2] <- t.test(typeList[[3]], typeList[[4]])$p.value
			}
			if(codonNumbers[5] > 5 & codonNumbers[6] > 5){
				pvalues[3] <- t.test(typeList[[5]], typeList[[6]])$p.value
			}
			if(codonNumbers[7] > 5 & codonNumbers[8] > 5){
				pvalues[4] <- t.test(typeList[[7]], typeList[[8]])$p.value
			}
			if(codonNumbers[9] > 5 & codonNumbers[10] > 5){
				pvalues[5] <- t.test(typeList[[9]], typeList[[10]])$p.value
			}
			if(codonNumbers[11] > 5 & codonNumbers[12] > 5){
				pvalues[6] <- t.test(typeList[[11]], typeList[[12]])$p.value
			}

			if(sum(pvalues, na.rm = TRUE)>0){
				ymax <- -log10(min(pvalues, na.rm = TRUE)) + -log10(min(pvalues, na.rm = TRUE))*0.1
			}
			if(ymax < 3) ymax <- 3.2
			ylab <- "-log10(p-value)"
			y_coords <- barplot(-log10(pvalues), main = main, ylab = ylab, 
					ylim = c(0, ymax), cex.names = 0.6)
			text(y_coords, -log10(pvalues), labels = formatC(pvalues, format = "e", digits = 1), pos = 3)
			abline(h = -log10(0.05), col = "orange")
			abline(h = -log10(0.01), col = "red")
			abline(h = -log10(0.001), col = "purple")

			# 20200422
			if(l == 18 | l == 21 | l == 24){
				if(i == 1 & l == 18)	typeLists[[type]] <- list()
				if(i == 1 & l == 18)	typeLists[[type]][["MPur_vs_MPyr"]] <- list()
				listName <- paste("UPST", i, "LEN", l, sep = "")
				typeLists[[type]][["MPur_vs_MPyr"]][[listName]] <- list()
				typeLists[[type]][["MPur_vs_MPyr"]][[listName]][["typeList"]] <- typeList[c(1,2)]
				typeLists[[type]][["MPur_vs_MPyr"]][[listName]][["codonNumbers"]] <- codonNumbers[c(1,2)]
				typeLists[[type]][["MPur_vs_MPyr"]][[listName]][["pvalue"]] <- pvalues[c(1)]
				typeLists[[type]][["OCG_vs_OAT"]][[listName]] <- list()
				typeLists[[type]][["OCG_vs_OAT"]][[listName]][["typeList"]] <- typeList[c(5,6)]
				typeLists[[type]][["OCG_vs_OAT"]][[listName]][["codonNumbers"]] <- codonNumbers[c(5,6)]
				typeLists[[type]][["OCG_vs_OAT"]][[listName]][["pvalue"]] <- pvalues[c(3)]
			}

		} # l-loop
		} # t-loop
		dev.off()
	}
	
	filename <- paste("Bases_RPKMthresh", RPKM.threshold, "_UPST1_LEN18-24_selected_.pdf", sep = "")
	pdf(file = filename, width = 6*2, height = 4*2)
	par(mfrow = c(2, 2))
	i <- 1
	
		type <- "WTwo3ATvsCDS"
		codonNumbers <- integer(length = 12)
		pvalues <- numeric(length = 6)
		DATA <- list()
		for(l in c(18, 21, 24)){
			AA <- (l-3)/3
			listName <- paste("UPST", i, "LEN", l, sep = "")
			codonNumbers[c((AA-4-1)*4+1, (AA-4-1)*4+2)] <- typeLists[[type]][["MPur_vs_MPyr"]][[listName]][["codonNumbers"]]
			codonNumbers[c((AA-4-1)*4+3, (AA-4-1)*4+4)] <- typeLists[[type]][["OCG_vs_OAT"]][[listName]][["codonNumbers"]]
			pvalues[(AA-4-1)*2+1] <- typeLists[[type]][["MPur_vs_MPyr"]][[listName]][["pvalue"]]
			pvalues[(AA-4-1)*2+2] <- typeLists[[type]][["OCG_vs_OAT"]][[listName]][["pvalue"]]

			listName2 <- paste("MPur", AA, sep = "")
			DATA[[listName2]] <- typeLists[[type]][["MPur_vs_MPyr"]][[listName]][["typeList"]]$MPur
			listName3 <- paste("MPyr", AA, sep = "")
			DATA[[listName3]] <- typeLists[[type]][["MPur_vs_MPyr"]][[listName]][["typeList"]]$MPyr
			listName4 <- paste("OCG", AA, sep = "")
			DATA[[listName4]] <- typeLists[[type]][["OCG_vs_OAT"]][[listName]][["typeList"]]$OCG
			listName5 <- paste("OAT", AA, sep = "")
			DATA[[listName5]] <- typeLists[[type]][["OCG_vs_OAT"]][[listName]][["typeList"]]$OAT
		}

		sub <- paste("n=", paste(codonNumbers, collapse = ","), sep = "")
		main <- paste(type, " RPKMthresh", RPKM.threshold, " UPST1-4 ", PAIR, sep = "")
		ylim <- c(min(unlist(DATA), na.rm = TRUE), max(unlist(DATA), na.rm = TRUE))
		ylab <- paste(type, " (log2)", sep = "")
		boxplot(DATA, outline = FALSE, main = main, sub = sub, 
			ylab = ylab, ylim = ylim, cex.axis = 0.8)
		beeswarm(DATA, corral = "wrap", col = rgb(1, 0, 0, alpha = 0.5), pch = 20, add = TRUE)

		ymax <- -log10(min(pvalues, na.rm = TRUE)) + -log10(min(pvalues, na.rm = TRUE))*0.1
		if(ymax < 3) ymax <- 3.2
		ylab <- "-log10(p-value)"
		y_coords <- barplot(-log10(pvalues), main = main, ylab = ylab, ylim = c(0, ymax), cex.names = 0.6)
		text(y_coords, -log10(pvalues), labels = formatC(pvalues, format = "e", digits = 1), pos = 3)
		abline(h = -log10(0.05), col = "orange")
		abline(h = -log10(0.01), col = "red")
		abline(h = -log10(0.001), col = "purple")
	
	dev.off()

	setwd(wd)
	LOGDIR <- paste("LOG_RPKM", RPKM.threshold, "_v10_4", sep = "")
	dir.create(LOGDIR)
	setwd(LOGDIR)
	filename <- paste("log_RPKMthresh", RPKM.threshold, ".txt", sep = "")
	write(LOG1, file = filename)

	setwd(wd3)
}

setwd(HOMEDIR)
dir.create("20200422_Logo_v10_4")
setwd("20200422_Logo_v10_4")
analyzeRatiouORFv10_4(inDF = bundledDF, RPKM.threshold = 6)

