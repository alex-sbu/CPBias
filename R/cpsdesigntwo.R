##' @description
##' Shuffle codons to generate sequences with different average codon pair scores relative to two poorly correlated reference codon pair biases.
##'
##' @details
##' An input sequence can be differentially recoded between two poorly correlated CPB references. A correlation test on two CPB reference tables can be computed directly with \code{\link{CPBcorr}}. CPSdesign.dual will attempt to create a recoded sequence characterized by two ideal codon pair scores relative to two different CPBs. See \code{\link{listCPB}} for a list of available CPB reference tables.
##'
##' To get best results with dual reference recoding it is not advised to use 'max' and 'min' inputs for ideal score. It is better to first determine the range of possible scores relative to each reference alone, and then use specific scores when recoding for two references. Use the \code{bind} argument if the 'max' or 'min' possible score is desired. The \code{bind} argument controls the preference and direction of recoding for either reference. It takes a four element numeric vector input, the first position of the vector controls the permissivity of scores \emph{less than} the ideal score for the first reference, the number in the second position controls scores \emph{greater than} ideal for the first reference, and the third and fourth numbers control scores less than and greater than ideal for the \emph{second reference}, respectively. Preference between directions and references is decided relative to each other, therefore any four identical numbers result in no biased preference. Increasing one value relative to the others will allow recoding to explore more sequences in the direction and reference designated by that position in the \code{bind} vector. Output scores can be further optimized by increasing the \code{buffer} value, if the ideal scores between two references are dramatically different or there is very little correlation between the reference CPBs. Increasing the \code{buffer} will increase the probability of alternating between references during a single permutation, whereas normally each round of codon shuffling is performed relative to one reference at a time.
##'
##' Recoded sequences can be saved as a fasta file with the \code{save} argument. Additional information about the recoded sequence is returned as an invisible list.
##'
##' @note
##' Both single and dual reference recoding offer the option of removing certain types of sequence elements by restricting which codons can be paired together. For example, restriction enzyme recognition sequences can be removed or prevented from appearing in the recoded sequence. These restricted sequences can be entered to the \code{restrictSeqs} argument as a character string with sequences seperated by commas. A comprehensive list of restriction enzyme recognition sequences expressed as R regular expressions is provided in this package, see \code{\link{REseqs}}. Restricted sequences are searched only on the single input strand, the \code{complementary} arguement allows for searching on the reverse complement strand. Depending on the type and number of sequences that are restricted this functionality can dramatically slow down recoding and restrict the CPS.
##'
##' @title Gene design with two codon pair bias references
##'
##' @param sequence Sequences can be input directly as a character string or as the file path to a fasta file. All sequences must be in the correct reading frame, stop codons or codons not defined in the translation table are not allowed and will generate an error.
##' @param reference CPB reference table (see \code{\link{CPBtable}}).
##' @param referenceTwo A second reference table for recoding relative to two CPB's.
##' @param score Ideal score relative to the first reference. Input can be numeric, `min', or `max'.
##' @param scoreTwo Ideal score relative to the second reference. Format same as \code{score}.
##' @param start Nucleotide position in the sequence.
##' @param end Nucleotide position in the sequence. If NULL, the last in frame nucleotide position is used.
##' @param cycles Optional input designating the number of recoding cycles. If empty, a minimal number of cycles is determined. Increasing the number of cycles may result in scores closer to the ideal value.
##' @param buffer Controls the probability of alternating recoding preference toward the other reference during a single round of codon shuffling. Accepts numeric input.
##' @param bind Controls the permissivity of scores greater than or less than the ideal score for each reference. Input is a four element numeric vector interpreted relative to each other. Position in the vector designates which direction and which reference. The first position controls scores less than the ideal first score, the second position controls scores greater than the ideal first score, and the third and fourth numbers control scores less than, and greater than the desired second score. Default values represent no bias toward either reference and direction.
##' @param transTable Alternative translation tables can be used.
##' @param restrictSeqs A string of comma separated sequences to remove or avoid in the input sequence while recoding. Search is performed 5' to 3' on the given strand. R regular expressions are allowed.
##' @param complementary If TRUE search for restricted sequences is also performed on the complementary strand 5' to 3'.
##' @param windowSize CPS line plots are smoothed by locally weighted polynomial regression where \code{windowSize} designates the number of nucleotides over which individual codon pair scores are smoothed. If NULL, smoothing spans 7.5\% of the sequence length.
##' @param save Save recoded sequence in fasta file.
##' @param name Output sequence name and fasta file name.
##' @param location Save location of output sequence.
##' @param draw If TRUE a line plot showing the local CPS along the length of the sequence is output to the graphics device during recoding.
##' @param silent If TRUE output to the console is suppressed.
##' @return
##' \item{firstoldCPS}{Average codon pair score of the input sequence relative to reference 1.}
##' \item{firstnewCPS}{Average codon pair score of recoded sequence relative to reference 1.}
##' \item{secondoldCPS}{Average codon pair score of the input sequence relative to reference 2.}
##' \item{secondnewCPS}{Average codon pair score of recoded sequence relative to reference 2.}
##' \item{codonchanges}{If codons were changed function will return an error.}
##' \item{mutations}{Number of point mutations generated by recoding.}
##' \item{oldCPSarray}{Vector of individual codon pair scores for the input sequence.}
##' \item{newCPSarray}{Vector of individual codon pair scores for the recoded sequence.}
##' \item{returnSeq}{The recoded sequence.}
##'
##' If restricted sequences are given an additional item is returned:
##' \item{restrSeqs}{The number of matches in the recoded sequence to the restricted sequences, value will include the complementary sequence is \code{complementary} is TRUE.}
##'
##' @examples
##' fastaLocation <- system.file('tbevns5.fasta', package = 'CPBias')
##' tbev <- importFasta(fastaLocation)[[1]]
##'
##' # A dual CPS recoding of TBE virus relative to the CPB of two of its natural hosts,
##' # Homo.sapiens and Aedes aegypti. Design strategy is to make Homo.sapiens CPS as high
##' # as possible and less than WT CPS in Aedes.aegypti.
##'
##' # If correlation between codon pair biases is too high, dual differential recoding may not
##' # be possible.
##' CPBcorr(Homo.sapiens, Aedes.aegypti)
##'
##' # First estimate the possible range of scores relative to both hosts
##' HumMax <- CPSdesign.single(tbev, Homo.sapiens, 'max', silent=TRUE, draw=FALSE)[[2]]
##' AedMin <- CPSdesign.single(tbev, Aedes.aegypti, 'min', silent=TRUE, draw=FALSE)[[2]]
##'
##' # Bind is set to prefer greater CPS in Homo.sapienss while greater CPS in Aedes
##' # is strongly prohibited. There is room to play with these settings.
##' tbevDual <- CPSdesign.dual(tbev, Homo.sapiens, Aedes.aegypti, .25, -.02,
##' bind=c(1,1000,1,.001))
##' tbevDual
CPSdesign.dual <- function(sequence, reference, referenceTwo, score, scoreTwo, start = 1, end = NULL, cycles = NULL, buffer = 0.5, bind = c(1, 1, 1, 1), transTable = standardTranslation, 
    restrictSeqs = NULL, complementary = FALSE, windowSize = NULL, save = FALSE, name = NULL, location = NULL, draw = TRUE, silent = FALSE) {
    
    if (any(grepl(".fasta", sequence))) {
        readSeq <- importFasta(sequence)[[1]]
    } else if (!is.null(sequence)) {
        readSeq <- tolower(as.character(sequence))
    }
    
    if (is.null(restrictSeqs) & complementary) {
        stop("Restricted sequence input is NULL.")
    }
    
    
    if (is.null(end)) {
        end <- nchar(readSeq)
    }
    windowSize <- windowSize
    
    seqORF <- substr(readSeq, start, end)
    seqVesselORF <- splitbyCodon(seqORF)
    seqORF <- paste0(seqVesselORF, collapse = "")
    if (is.null(cycles)) {
        cycles <- round((nchar(seqORF)/5), digits = 0)
        if (cycles < 500) {
            cycles <- 500
        }
    }
    
    seqREC <- substr(readSeq, start, end)
    DrawStart <- CPScalcTwo(seqREC, reference, referenceTwo, blank = draw, draw = draw, windowSize = windowSize)
    startingCPS <- DrawStart[1]
    startingCPSSec <- DrawStart[2]
    
    determineThresholdDual <- FALSE
    determineThreshold <- FALSE
    bind <- 1/bind
    if (is.numeric(score) & is.numeric(scoreTwo)) {
        firstThreshold <- round(score, 5)
        secondThreshold <- round(scoreTwo, 5)
        switch(which.min(c(score, scoreTwo)), {
            probFirstHost <- 100
        }, {
            probFirstHost <- 0.001
        })
        if (score < startingCPS) {
            directMin <- TRUE
        } else if (score > startingCPS) {
            directMin <- FALSE
        }
        if (scoreTwo < startingCPSSec) {
            directMinSec <- TRUE
        } else if (scoreTwo > startingCPSSec) {
            directMinSec <- FALSE
        }
        
    } else if (!is.numeric(score) & !is.numeric(scoreTwo)) {
        determineThresholdDual <- TRUE
        probFirstHost <- 100
        if (tolower(score) == "min") {
            directMin <- TRUE
        } else if (tolower(score) == "max") {
            directMin <- FALSE
        } else {
            stop("Improper input")
        }
        if (tolower(scoreTwo) == "min") {
            directMinSec <- TRUE
        } else if (tolower(scoreTwo) == "max") {
            directMinSec <- FALSE
        } else {
            stop("Improper input")
        }
        
    } else if (is.numeric(score) & !is.numeric(scoreTwo)) {
        probFirstHost <- 0.001
        determineThreshold <- TRUE
        determineFirst <- FALSE
        firstThreshold <- round(score, 5)
        if (tolower(scoreTwo) == "min") {
            directMinSec <- TRUE
        } else if (tolower(scoreTwo) == "max") {
            directMinSec <- FALSE
        } else {
            stop("Improper input")
        }
        if (score < startingCPS) {
            directMin <- TRUE
        } else if (score > startingCPS) {
            directMin <- FALSE
        }
        
    } else if (!is.numeric(score) & is.numeric(scoreTwo)) {
        probFirstHost <- 100
        determineThreshold <- TRUE
        determineFirst <- TRUE
        secondThreshold <- round(scoreTwo, 5)
        if (tolower(score) == "min") {
            directMin <- TRUE
        } else if (tolower(score) == "max") {
            directMin <- FALSE
        } else {
            stop("Improper input")
        }
        if (scoreTwo < startingCPSSec) {
            directMinSec <- TRUE
        } else if (scoreTwo > startingCPSSec) {
            directMinSec <- FALSE
        }
    }
    
    toSep <- ifelse(directMin == directMinSec, FALSE, TRUE)
    lockedDirMin <- directMin
    lockedDirMinSec <- directMinSec
    codonpairs <- as.character(reference[, 1])
    templateCPStable <- makeCPBtemplate(transTable)
    stopCodons <- tolower(as.character(transTable[which(tolower(transTable[, 2]) == "stop"), 1]))
    if (any(stopCodons %in% seqVesselORF)) {
        stop("Sequence contains stop codon.")
    }
    indivCPS <- as.vector(reference[, 2])
    codonpairsSec <- as.character(referenceTwo[, 1])
    indivCPSSec <- as.vector(referenceTwo[, 2])
    AA <- translateSeq(seqORF, templateCPStable)
    AAlist <- as.character(unique(AA))
    totalcycles <- cycles
    if (totalcycles < 200) {
        stop("Number of cycles is too small.")
    }
    loopsize <- 3
    probSep <- 0
    probBuffer <- 0
    repeati <- 0
    repeatScore <- numeric()
    repeatSeq <- list()
    repeatiSec <- 0
    repeatScoreSec <- numeric()
    repeatSeqSec <- list()
    trackREsum <- numeric()
    toRestrictSeqs <- FALSE
    ## Restrict allowed sequences
    if (!is.null(restrictSeqs)) {
        notAllowedRE <- strsplit(restrictSeqs, split = ",")[[1]]
    }
    if (!silent) {
        pb <- txtProgressBar(min = 0, max = totalcycles, initial = 0, char = ".", width = (getOption("width")/2), style = 3)
    }
    repeat {
        returnSep <- sample(c(TRUE, FALSE), 1, prob = c(probSep, 1))
        for (ll in 1:loopsize) {
            if (!is.null(restrictSeqs)) {
                toRestrictSeqs <- sample(c(TRUE, FALSE), 1, prob = c(5, 1))
            }
            FirstHost <- sample(c(TRUE, FALSE), 1, prob = c(probFirstHost, 1))
            master <- sample(1:length(AAlist), 1)
            seqVessel <- splitbyCodon(seqREC)
            
            positionA <- which(AA == AAlist[master])
            if (length(positionA) == 1) {
                next
            }
            if (positionA[1] == 1) {
                positionB <- c(2, positionA[2:length(positionA)] - 1)
                positionC <- positionA + 1
            } else if (positionA[length(positionA)] == length(seqVessel)) {
                positionB <- positionA - 1
                positionC <- c(positionA[1:(length(positionA) - 1)] + 1, positionA[length(positionA)] - 1)
            } else if (positionA[length(positionA)] == length(seqVessel) & positionA[1] == 1) {
                positionB <- c(2, positionA[2:length(positionA)] - 1)
                positionC <- c(positionA[1:(length(positionA) - 1)] + 1, positionA[length(positionA)] - 1)
            } else {
                positionB <- positionA - 1
                positionC <- positionA + 1
            }
            
            codonsA <- seqVessel[positionA]
            codonsB <- seqVessel[positionB]
            codonsC <- seqVessel[positionC]
            codonArrayBA <- outer(codonsB, codonsA, "paste0")
            codonArrayAC <- outer(codonsA, codonsC, "paste0")
            cpsBA <- sapply(match(codonArrayBA, codonpairs), function(x) indivCPS[x])
            cpsArrayBA <- matrix(cpsBA, ncol = ncol(codonArrayBA))
            cpsAC <- sapply(match(codonArrayAC, codonpairs), function(x) indivCPS[x])
            cpsArrayAC <- matrix(cpsAC, ncol = ncol(codonArrayAC))
            cpsArrayABC <- (cpsArrayBA + cpsArrayAC)/2
            cpsfinal <- (cpsArrayABC + t(cpsArrayABC))/2
            codonArrayBASec <- outer(codonsB, codonsA, "paste0")
            codonArrayACSec <- outer(codonsA, codonsC, "paste0")
            cpsBASec <- sapply(match(codonArrayBASec, codonpairsSec), function(x) indivCPSSec[x])
            cpsArrayBASec <- matrix(cpsBASec, ncol = ncol(codonArrayBASec))
            cpsACSec <- sapply(match(codonArrayACSec, codonpairsSec), function(x) indivCPSSec[x])
            cpsArrayACSec <- matrix(cpsACSec, ncol = ncol(codonArrayACSec))
            cpsArrayABCSec <- (cpsArrayBASec + cpsArrayACSec)/2
            cpsfinalSec <- (cpsArrayABCSec + t(cpsArrayABCSec))/2
            
            HOLDcpsfinal <- cpsfinal
            HOLDdirectMin <- directMin
            BUFFERcpsfinal <- cpsfinalSec
            BUFFERdirectMin <- directMinSec
            
            if (!FirstHost) {
                HOLDdirectMin <- directMinSec
                HOLDcpsfinal <- cpsfinalSec
                BUFFERcpsfinal <- cpsfinal
                BUFFERdirectMin <- directMin
            }
            
            
            if (HOLDdirectMin) {
                colRanked <- cbind(seq(1, ncol(HOLDcpsfinal), by = 1), apply(HOLDcpsfinal, 2, min))
                colOrdered <- as.numeric(colRanked[order(colRanked[, 2]), 1])
            } else {
                colRanked <- cbind(seq(1, ncol(HOLDcpsfinal), by = 1), apply(HOLDcpsfinal, 2, max))
                colOrdered <- as.numeric(colRanked[order(colRanked[, 2], decreasing = TRUE), 1])
            }
            if (BUFFERdirectMin) {
                BUFFERcolRanked <- cbind(seq(1, ncol(BUFFERcpsfinal), by = 1), apply(BUFFERcpsfinal, 2, min))
                BUFFERcolOrdered <- as.numeric(BUFFERcolRanked[order(BUFFERcolRanked[, 2]), 1])
            } else {
                BUFFERcolRanked <- cbind(seq(1, ncol(BUFFERcpsfinal), by = 1), apply(BUFFERcpsfinal, 2, max))
                BUFFERcolOrdered <- as.numeric(BUFFERcolRanked[order(BUFFERcolRanked[, 2], decreasing = TRUE), 1])
            }
            
            rowIndex <- seq(1, length(colOrdered), by = 1)
            
            
            repeat {
                if (sample(c(TRUE, FALSE), 1, prob = c(probBuffer, 1))) {
                  
                  codonCol <- as.numeric(BUFFERcolOrdered[1])
                  
                  if (toRestrictSeqs) {
                    
                    innerTrackRE <- 0
                    rowIndexRE <- rowIndex
                    FoundBadRE <- numeric()
                    if (complementary) {
                      FoundBadRErevomp <- numeric()
                    }
                    innerTrackRow <- numeric()
                    repeat {
                      innerTrackRE <- innerTrackRE + 1
                      
                      if (length(rowIndexRE) == 0) {
                        codonRow <- ifelse(BUFFERdirectMin, {
                          as.numeric(rowIndex[which.min(BUFFERcpsfinal[rowIndex, BUFFERcolOrdered[1]])])
                        }, {
                          as.numeric(rowIndex[which.max(BUFFERcpsfinal[rowIndex, BUFFERcolOrdered[1]])])
                        })
                        (break)()
                      }
                      
                      codonRow <- ifelse(BUFFERdirectMin, {
                        as.numeric(rowIndexRE[which.min(BUFFERcpsfinal[rowIndexRE, BUFFERcolOrdered[1]])])
                      }, {
                        as.numeric(rowIndexRE[which.max(BUFFERcpsfinal[rowIndexRE, BUFFERcolOrdered[1]])])
                      })
                      
                      seqVesselHolder <- seqVessel
                      seqVesselHolder[c(positionA[codonRow], positionA[codonCol])] <- seqVesselHolder[c(positionA[codonCol], positionA[codonRow])]
                      seqRECre <- paste0(seqVesselHolder, collapse = "")
                      FoundBadRE[innerTrackRE] <- sum(sapply(1:length(notAllowedRE), function(x) {
                        any(grepl(notAllowedRE[x], seqRECre))
                      }), na.rm = TRUE)
                      
                      if (complementary) {
                        seqrevcomp <- revcomp(seqRECre)
                        FoundBadRErevomp[innerTrackRE] <- sum(sapply(1:length(notAllowedRE), function(x) {
                          any(grepl(notAllowedRE[x], seqrevcomp))
                        }), na.rm = TRUE)
                        FoundBadRE <- FoundBadRE + FoundBadRErevomp
                      }
                      
                      innerTrackRow[innerTrackRE] <- codonRow
                      if (FoundBadRE[innerTrackRE] == 0) {
                        (break)()
                      }
                      rowIndexRE <- rowIndexRE[-which(rowIndexRE == codonRow)]
                      (next)()
                    }
                    codonRow <- innerTrackRow[which(FoundBadRE == min(FoundBadRE, na.rm = TRUE))[1]]
                    seqVessel[c(positionA[codonRow], positionA[codonCol])] <- seqVessel[c(positionA[codonCol], positionA[codonRow])]
                    
                  } else {
                    
                    codonRow <- ifelse(BUFFERdirectMin, {
                      as.numeric(rowIndex[which.min(BUFFERcpsfinal[rowIndex, BUFFERcolOrdered[1]])])
                    }, {
                      as.numeric(rowIndex[which.max(BUFFERcpsfinal[rowIndex, BUFFERcolOrdered[1]])])
                    })
                    seqVessel[c(positionA[codonRow], positionA[codonCol])] <- seqVessel[c(positionA[codonCol], positionA[codonRow])]
                  }
                  
                  BUFFERcolOrdered <- BUFFERcolOrdered[-which(BUFFERcolOrdered == codonRow)]
                  colOrdered <- colOrdered[-which(colOrdered == codonRow)]
                  rowIndex <- rowIndex[-which(rowIndex == codonCol)]
                  
                  if (codonRow != codonCol) {
                    BUFFERcolOrdered <- BUFFERcolOrdered[-which(BUFFERcolOrdered == codonCol)]
                    colOrdered <- colOrdered[-which(colOrdered == codonCol)]
                    rowIndex <- rowIndex[-which(rowIndex == codonRow)]
                  }
                  
                  if (length(BUFFERcolOrdered) == 0) {
                    (break)()
                  }
                } else {
                  codonCol <- as.numeric(colOrdered[1])
                  if (toRestrictSeqs) {
                    
                    innerTrackRE <- 0
                    rowIndexRE <- rowIndex
                    FoundBadRE <- numeric()
                    if (complementary) {
                      FoundBadRErevomp <- numeric()
                    }
                    innerTrackRow <- numeric()
                    repeat {
                      
                      innerTrackRE <- innerTrackRE + 1
                      
                      if (length(rowIndexRE) == 0) {
                        codonRow <- ifelse(HOLDdirectMin, {
                          as.numeric(rowIndex[which.min(HOLDcpsfinal[rowIndex, colOrdered[1]])])
                        }, {
                          as.numeric(rowIndex[which.max(HOLDcpsfinal[rowIndex, colOrdered[1]])])
                        })
                        (break)()
                      }
                      
                      
                      codonRow <- ifelse(HOLDdirectMin, {
                        as.numeric(rowIndexRE[which.min(HOLDcpsfinal[rowIndexRE, colOrdered[1]])])
                      }, {
                        as.numeric(rowIndexRE[which.max(HOLDcpsfinal[rowIndexRE, colOrdered[1]])])
                      })
                      
                      
                      
                      seqVesselHolder <- seqVessel
                      seqVesselHolder[c(positionA[codonRow], positionA[codonCol])] <- seqVesselHolder[c(positionA[codonCol], positionA[codonRow])]
                      seqRECre <- paste0(seqVesselHolder, collapse = "")
                      FoundBadRE[innerTrackRE] <- sum(sapply(1:length(notAllowedRE), function(x) {
                        any(grepl(notAllowedRE[x], seqRECre))
                      }), na.rm = TRUE)
                      
                      if (complementary) {
                        seqrevcomp <- revcomp(seqRECre)
                        FoundBadRErevomp[innerTrackRE] <- sum(sapply(1:length(notAllowedRE), function(x) {
                          any(grepl(notAllowedRE[x], seqrevcomp))
                        }), na.rm = TRUE)
                        FoundBadRE <- FoundBadRE + FoundBadRErevomp
                      }
                      
                      
                      innerTrackRow[innerTrackRE] <- codonRow
                      if (FoundBadRE[innerTrackRE] == 0) {
                        (break)()
                      }
                      rowIndexRE <- rowIndexRE[-which(rowIndexRE == codonRow)]
                      (next)()
                    }
                    codonRow <- innerTrackRow[which(FoundBadRE == min(FoundBadRE, na.rm = TRUE))[1]]
                    seqVessel[c(positionA[codonRow], positionA[codonCol])] <- seqVessel[c(positionA[codonCol], positionA[codonRow])]
                    
                  } else {
                    
                    codonRow <- ifelse(HOLDdirectMin, {
                      as.numeric(rowIndex[which.min(HOLDcpsfinal[rowIndex, colOrdered[1]])])
                    }, {
                      as.numeric(rowIndex[which.max(HOLDcpsfinal[rowIndex, colOrdered[1]])])
                    })
                    seqVessel[c(positionA[codonRow], positionA[codonCol])] <- seqVessel[c(positionA[codonCol], positionA[codonRow])]
                  }
                  
                  BUFFERcolOrdered <- BUFFERcolOrdered[-which(BUFFERcolOrdered == codonRow)]
                  colOrdered <- colOrdered[-which(colOrdered == codonRow)]
                  rowIndex <- rowIndex[-which(rowIndex == codonCol)]
                  
                  if (codonRow != codonCol) {
                    BUFFERcolOrdered <- BUFFERcolOrdered[-which(BUFFERcolOrdered == codonCol)]
                    colOrdered <- colOrdered[-which(colOrdered == codonCol)]
                    rowIndex <- rowIndex[-which(rowIndex == codonRow)]
                  }
                  
                  if (length(colOrdered) == 0) {
                    (break)()
                  }
                  
                }
            }
            seqREC <- paste0(seqVessel, collapse = "")
        }
        if (draw) {
            Sys.sleep(0.1)
        }
        drawTracking <- CPScalcTwo(seqREC, reference, referenceTwo, blank = FALSE, draw = draw, windowSize = windowSize)
        tracking <- drawTracking[1]
        repeati <- repeati + 1
        repeatScore[repeati] <- as.numeric(tracking)
        repeatSeq[[repeati]] <- seqREC
        trackingSec <- drawTracking[2]
        repeatScoreSec[repeati] <- as.numeric(trackingSec)
        if (!is.null(restrictSeqs)) {
            trackREsum[repeati] <- sum(sapply(1:length(notAllowedRE), function(x) {
                any(grepl(notAllowedRE[x], seqREC))
            }), na.rm = TRUE)
        }
        if (!silent) {
            setTxtProgressBar(pb, repeati)
        }
        stepsize <- 20
        initialStep <- ifelse(determineThresholdDual, 200, 100)
        stepladder <- seq(initialStep, totalcycles, by = stepsize)
        
        roundscores <- round(repeatScore, 5)
        roundscoresSec <- round(repeatScoreSec, 5)
        
        if (repeati > (initialStep - 1)) {
            
            thresholdDistance <- abs(roundscores - firstThreshold)
            thresholdDistanceSec <- abs(roundscoresSec - secondThreshold)
            heavyScores <- 1/(thresholdDistance^2 + 1)
            heavyScoresSec <- 1/(thresholdDistanceSec^2 + 1)
            
            heavyScores[roundscores < firstThreshold] <- 1/(bind[1] * thresholdDistance[roundscores < firstThreshold]^2 + 1)
            heavyScores[roundscores > firstThreshold] <- 1/(bind[2] * thresholdDistance[roundscores > firstThreshold]^2 + 1)
            heavyScoresSec[roundscoresSec < secondThreshold] <- 1/(bind[3] * thresholdDistanceSec[roundscoresSec < secondThreshold]^2 + 1)
            heavyScoresSec[roundscoresSec > secondThreshold] <- 1/(bind[4] * thresholdDistanceSec[roundscoresSec > secondThreshold]^2 + 1)
            heavyScoresSec <- -heavyScoresSec
            if (returnSep) {
                if (toRestrictSeqs) {
                  trackREindex <- seq(1, length(trackREsum), by = 1)
                  trackREindex <- trackREindex[which(trackREsum == min(trackREsum, na.rm = TRUE))]
                  REheavyscores <- heavyScores[trackREindex]
                  REheavyScoresSec <- heavyScoresSec[trackREindex]
                  seqREC <- repeatSeq[[as.numeric(trackREindex[which.max(abs(REheavyscores - REheavyScoresSec))])]]
                } else {
                  seqREC <- repeatSeq[[which.max(abs(heavyScores - heavyScoresSec))]]
                }
            }
        }
        
        if (any(repeati == stepladder)) {
            scoremean <- mean(abs(heavyScores[(length(heavyScores) - stepsize):length(heavyScores)]), na.rm = TRUE)
            scoremeanSec <- mean(abs(heavyScoresSec[(length(heavyScoresSec) - stepsize):length(heavyScoresSec)]), na.rm = TRUE)
            repeatmean <- mean(repeatScore[(length(repeatScore) - stepsize):length(repeatScore)], na.rm = TRUE)
            repeatmeanSec <- mean(repeatScoreSec[(length(repeatScoreSec) - stepsize):length(repeatScoreSec)], na.rm = TRUE)
            
            switch(which.min(c(abs(scoremean), abs(scoremeanSec))), {
                probFirstHost <- 10
                directMin <- ifelse(repeatmean > firstThreshold, TRUE, FALSE)
            }, {
                probFirstHost <- 0.1
                directMinSec <- ifelse(repeatmeanSec > secondThreshold, TRUE, FALSE)
            })
        }
        
        if (determineThresholdDual & repeati == (initialStep/2)) {
            probFirstHost <- 0.001
        }
        switch((repeati - (initialStep - 1)), {
            probSep <- 3
            probBuffer <- buffer
        })
        switch((repeati - (totalcycles/2)), {
            probSep <- 10
        })
        
        if (repeati == (initialStep - 2)) {
            if (determineThreshold) {
                if (determineFirst) {
                  firstThreshold <- ifelse(lockedDirMin, min(repeatScore), max(repeatScore))
                  firstThreshold <- round(firstThreshold, 5)
                } else {
                  secondThreshold <- ifelse(lockedDirMinSec, min(repeatScoreSec), max(repeatScoreSec))
                  secondThreshold <- round(secondThreshold, 5)
                }
            }
            if (determineThresholdDual) {
                firstThreshold <- ifelse(lockedDirMin, min(repeatScore), max(repeatScore))
                firstThreshold <- round(firstThreshold, 5)
                secondThreshold <- ifelse(lockedDirMinSec, min(repeatScoreSec), max(repeatScoreSec))
                secondThreshold <- round(secondThreshold, 5)
            }
        }
        
        switch((repeati - (totalcycles - 1)), {
            (break)()
        })
    }
    
    if (!silent) {
        close(pb)
    }
    roundscores <- round(repeatScore, 5)
    roundscoresSec <- round(repeatScoreSec, 5)
    thresholdDistance <- abs(roundscores - firstThreshold)
    thresholdDistanceSec <- abs(roundscoresSec - secondThreshold)
    heavyScores <- 1/(thresholdDistance^2 + 1)
    heavyScoresSec <- 1/(thresholdDistanceSec^2 + 1)
    
    heavyScores[roundscores < firstThreshold] <- 1/(bind[1] * thresholdDistance[roundscores < firstThreshold]^2 + 1)
    heavyScores[roundscores > firstThreshold] <- 1/(bind[2] * thresholdDistance[roundscores > firstThreshold]^2 + 1)
    heavyScoresSec[roundscoresSec < secondThreshold] <- 1/(bind[3] * thresholdDistanceSec[roundscoresSec < secondThreshold]^2 + 1)
    heavyScoresSec[roundscoresSec > secondThreshold] <- 1/(bind[4] * thresholdDistanceSec[roundscoresSec > secondThreshold]^2 + 1)
    heavyScoresSec <- (-heavyScoresSec)
    
    if (toRestrictSeqs) {
        trackREindex <- seq(1, length(trackREsum), by = 1)
        trackREindex <- trackREindex[which(trackREsum == min(trackREsum, na.rm = TRUE))]
        REheavyscores <- heavyScores[trackREindex]
        REheavyScoresSec <- heavyScoresSec[trackREindex]
        seqREC <- repeatSeq[[as.numeric(trackREindex[which.max(abs(REheavyscores - REheavyScoresSec))])]]
    } else {
        seqREC <- repeatSeq[[which.max(abs(heavyScores - heavyScoresSec))]]
    }
    
    finalCPS <- CPScalcTwo(seqREC, reference, referenceTwo, blank = FALSE, draw = FALSE, windowSize = windowSize)[1]
    finalCPSSec <- CPScalcTwo(seqREC, reference, referenceTwo, blank = FALSE, draw = FALSE, windowSize = windowSize)[2]
    compare.aa <- translateSeq(seqORF, transTable) == translateSeq(seqREC, transTable)
    if (any(!compare.aa)) {
        message("Error: Amino acids were changed.")
    }
    
    if (save) {
        if (is.null(location)) {
            location <- getwd()
        }
        if (is.null(name)) {
            if (any(grepl(".fasta", sequence))) {
                name <- gsub(".fasta", " CPSdesign.fasta", sequence)
            } else {
                name <- "CPSdesign.fasta"
            }
        }
        exportFasta(sequence = seqREC, name = name, location = location)
        
    }
    if (draw) {
        CPSplot(seqOne = seqREC, refOne = reference, seqTwo = seqORF, refTwo = referenceTwo, start = 1, end = NULL, windowSize = windowSize)
    }
    compare.codons <- countCodons(seqREC, transTable = transTable) == countCodons(seqORF, transTable = transTable)
    if (any(!compare.codons)) {
        numberofCodonChanges <- as.numeric(table(compare.codons)["FALSE"])
    } else {
        numberofCodonChanges <- 0
    }
    numberofMutations <- as.numeric(table(strsplit(seqORF, "")[[1]] == strsplit(seqREC, "")[[1]])["FALSE"])
    oldCPSarray <- CPScalcOne(seqORF, reference, average = FALSE, blank = FALSE, windowSize = windowSize)
    newCPSarray <- CPScalcOne(seqREC, reference, average = FALSE, blank = FALSE, windowSize = windowSize)
    
    if (is.null(restrictSeqs)) {
        ReturnTwo <- list(startingCPS, finalCPS, startingCPSSec, finalCPSSec, numberofCodonChanges, numberofMutations, oldCPSarray, newCPSarray, seqREC)
        names(ReturnTwo) <- c("firstoldCPS", "firstnewCPS", "secondoldCPS", "secondnewCPS", "codonchanges", "mutations", "oldCPSarray", "newCPSarray", "returnSeq")
    } else {
        reSeqs <- sum(sapply(1:length(notAllowedRE), function(x) grepl(notAllowedRE[x], seqREC)), na.rm = TRUE)
        ReturnTwo <- list(startingCPS, finalCPS, startingCPSSec, finalCPSSec, numberofCodonChanges, numberofMutations, oldCPSarray, newCPSarray, seqREC, reSeqs)
        names(ReturnTwo) <- c("firstoldCPS", "firstnewCPS", "secondoldCPS", "secondnewCPS", "codonchanges", "mutations", "oldCPSarray", "newCPSarray", "returnSeq", 
            "restrSeqs")
    }
    if (!silent) {
        cat("\n Reference 1 \n", "Original mean codon pair score: ", ReturnTwo$firstoldCPS, "\n", "New mean codon pair score: ", ReturnTwo$firstnewCPS, "\n", "\n Reference 2 \n", 
            "Original mean codon pair score: ", ReturnTwo$secondoldCPS, "\n", "New mean codon pair score: ", ReturnTwo$secondnewCPS, "\n \n", "Number of codon changes: ", 
            ReturnTwo$codonchanges, "\n", "Number of mutations: ", as.numeric(ReturnTwo$mutations), "\n")
    }
    
    invisible(ReturnTwo)
} 
