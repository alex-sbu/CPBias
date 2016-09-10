##' @description
##' Shuffle codons to generate sequences with altered average codon pair scores, or random permutations, relative to one reference codon pair bias.
##'
##' @details
##' This function optimizes the shuffling of existing codons in a protein coding sequence while preserving the order of amino acids. Codon usage is not changed by recoding because codons are not added or removed in the process. Shuffling is directed toward an ideal average codon pair score relative to a reference codon pair bias (CPB). CPB references are calculations of observed to expected codon pair frequencies performed on a large number of CDS sequences. See \code{\link{listCPB}} for a list of available CPB reference tables.
##'
##' The CPSdesign algorithm generates multiple sequence permutations. Returned sequences are by default those with a score (CPS) closest to the ideal (designated by the \code{score} argument). The shuffling algorithm can also be set to favor sequences dissimilar to the original sequence at any possible CPS, by setting \code{scramble} to TRUE. Scrambling will preferentially select codon positions different from the original sequence, however extreme scores may not be possible. To fully maximize the number of codon position differences without regard to the ideal score set \code{maxmutations} to TRUE.
##'
##' Scrambling the sequence is not the same as randomly shuffling codons. To generate a true random permutation of existing codons enter ``random'' for the \code{score} argument.
##'
##' Recoded sequences can be saved as a fasta file with the \code{save} argument. Additional information about the recoded sequence is returned as an invisible list.
##'
##' @note
##' Both single and dual reference recoding offer the option of removing certain types of sequence elements by restricting which codons can be paired together. For example, restriction enzyme recognition sequences can be removed or prevented from appearing in the recoded sequence. These restricted sequences can be entered to the \code{restrictSeqs} argument as a character string with sequences seperated by commas. A comprehensive list of restriction enzyme recognition sequences expressed as R regular expressions is provided in this package, see \code{\link{REseqs}}. Restricted sequences are searched only on the single input strand, the \code{complementary} arguement allows for searching on the reverse complement strand. Depending on the type and number of sequences that are restricted this functionality can dramatically slow down recoding and restrict the CPS.
##'
##'
##'
##' @title Gene design with a single codon pair bias reference
##'
##' @param sequence Sequences can be input directly as a character string or as the file path to a fasta file. All sequences must be in the correct reading frame, stop codons or codons not defined in the translation table are not allowed and will generate an error.
##' @param reference CPB reference table (see \code{\link{CPBtable}}).
##' @param score Ideal score relative to the first reference. Input can be numeric, `min', `max', or `random'.
##' @param start Nucleotide position in the sequence.
##' @param end Nucleotide position in the sequence. If NULL, the last in frame nucleotide position is used.
##' @param cycles Optional input designating the number of recoding cycles. If empty, a minimal number of cycles is determined. Increasing the number of cycles may result in scores closer to the ideal value.
##' @param scramble Optional TRUE or FALSE input designating whether priority should be given to increasing the number of mutations. If NULL, \code{scramble} is set to FALSE.
##' @param maxmutations If \code{scramble} and \code{maxmutations} are TRUE, the sequence generated with the greatest number of mutations is returned, with little control over the CPS.
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
##' \item{oldCPS}{Average codon pair score of the input sequence.}
##' \item{newCPS}{Average codon pair score of recoded sequence.}
##' \item{codonchanges}{If codons were changed function will return an error.}
##' \item{mutations}{Number of point mutations generated by recoding.}
##' \item{oldCPSarray}{Vector of individual codon pair scores for the input sequence.}
##' \item{newCPSarray}{Vector of individual codon pair scores for the recoded sequence.}
##' \item{returnSeq}{The recoded sequence.}
##'
##' If restricted sequences are given an additional item is returned:
##' \item{restrSeqs}{The number of matches in the recoded sequence to the restricted sequences, value will include the complementary sequence is \code{complementary} is TRUE.}
##' @examples
##' fastaLocation <- system.file('tbevns5.fasta', package = 'CPBias')
##' tbev <- importFasta(fastaLocation)[[1]]
##'
##' # Create a Homo.sapiens 'min' using 300 cycles, and save the recoded sequence.
##' CPSdesign.single(tbev, Homo.sapiens, 'min', cycles = 300, name = 'demoseq CPS min.fasta',
##' save = TRUE)
##'
##' # Create a scrambled sequence with a wild-type Homo.sapiens average CPS and omit some
##' # restriction enzyme sequences.
##'
##' # Find restriction enzyme recognition sequences in REseqs
##' selEnz <- which(REseqs[,1] %in% c('PfoI','SmlI','PflFI'))
##'
##' # Create comma separated string containing regex versions of the recognition sequences
##' omitRE <- paste0(REseqs[selEnz,5], collapse=',')
##'
##' # Get WT CPS relative to Homo.sapiens CPB by running CPScalc in silent mode
##' tbevScrambled <- CPSdesign.single(tbev, Homo.sapiens, CPScalc(tbev, Homo.sapiens, silent=TRUE,
##' draw= FALSE)[[1]], scramble = TRUE, restrictSeqs = omitRE)
CPSdesign.single <- function(sequence, reference, score, start = 1, end = NULL, cycles = NULL, scramble = FALSE, maxmutations = FALSE, transTable = standardTranslation, 
    restrictSeqs = NULL, complementary = FALSE, windowSize = NULL, save = FALSE, name = NULL, location = NULL, draw = TRUE, silent = FALSE) {
    
    if (any(grepl(".fasta", sequence))) {
        readSeq <- importFasta(sequence)[[1]]
    } else if (!is.null(sequence)) {
        readSeq <- tolower(as.character(sequence))
    }
    
    if (is.null(end)) {
        end <- nchar(readSeq)
    }
    
    if (is.null(restrictSeqs) & complementary) {
        stop("Restricted sequence input is NULL.")
    }
    
    seqLength <- end - start
    templateCPStable <- makeCPBtemplate(transTable)
    stopCodons <- tolower(as.character(transTable[which(tolower(transTable[, 2]) == "stop"), 1]))
    seqORF <- substr(readSeq, start, end)
    seqVesselORF <- splitbyCodon(seqORF)
    seqORF <- paste0(seqVesselORF, collapse = "")
    seqREC <- seqORF
    
    if (any(stopCodons %in% seqVesselORF)) {
        stop("Sequence contains stop codon.")
    }
    
    if (is.numeric(cycles)) {
        totalcycles <- cycles
        if (totalcycles < 100) {
            stop("Cycle number too small.")
        }
    } else {
        totalcycles <- round(seqLength/10, digits = 0)
        if (totalcycles < 500 & is.numeric(score)) {
            totalcycles <- 500
        }
        if (totalcycles < 300 & !is.numeric(score)) {
            totalcycles <- 300
        }
    }
    
    if (draw) {
        startingCPS <- CPScalcOne(seqORF, reference, blank = TRUE, windowSize = windowSize)
    } else {
        startingCPS <- CPScalcOne(seqORF, reference, blank = FALSE, windowSize = windowSize, draw = FALSE)
    }
    
    trueRandom <- FALSE
    if (is.numeric(score)) {
        scoreGiven <- TRUE
        desiredScore <- score
        directMin <- if (score < startingCPS) {
            TRUE
        } else if (score > startingCPS) {
            FALSE
        }
        stepsize <- 20
        initialstep <- 100
        stepladder <- seq(initialstep, totalcycles, by = stepsize)
        loopsize <- 5
        
    } else if (is.character(score)) {
        scoreGiven <- FALSE
        if (tolower(score) == "min") {
            directMin <- TRUE
        } else if (tolower(score) == "max") {
            directMin <- FALSE
        } else if (tolower(score) == "random") {
            trueRandom <- TRUE
            # doesn't do anything, assigned so other stuff doesnt break ...
            directMin <- FALSE
        }
        loopsize <- 10
        
    } else if (is.null(score)) {
        stop("Error: enter score")
    }
    
    
    lockedDirMin <- directMin
    codonpairs <- as.character(reference[, 1])
    indivCPS <- as.vector(reference[, 2])
    AA <- translateSeq(seqORF, templateCPStable)
    AAlist <- as.character(unique(AA))
    Trandom <- ifelse(scoreGiven, 10, 2)
    Tbest <- 0
    repeati <- 0
    repeatScore <- numeric()
    repeatMutations <- numeric()
    repeatSeq <- list()
    trackREsum <- numeric()
    toRestrictSeqs <- FALSE
    ## Restrict allowed sequences
    if (!is.null(restrictSeqs)) {
        notAllowedRE <- strsplit(restrictSeqs, split = ",")[[1]]
        if (trueRandom) {
            warning("Generated sequence is not a true random due to sequence restrictions.")
        }
    }
    if (!silent) {
        pb <- txtProgressBar(min = 0, max = totalcycles, initial = 0, char = ".", width = (getOption("width")/2), style = 3)
    }
    repeat {
        
        returnBest <- sample(c(TRUE, FALSE), 1, prob = c(Tbest, 1))
        for (ll in 1:loopsize) {
            if (!is.null(restrictSeqs)) {
                toRestrictSeqs <- sample(c(TRUE, FALSE), 1, prob = c(5, 1))
            }
            randomWalk <- sample(c(TRUE, FALSE), 1, prob = c(Trandom, 1))
            
            if (trueRandom) {
                randomWalk <- TRUE
                returnBest <- FALSE
            }
            
            directMin <- if (randomWalk & scoreGiven) {
                sample(c(TRUE, FALSE), 1, prob = c(1, 1))
            } else {
                lockedDirMin
            }
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
            
            if (randomWalk) {
                colRanked <- cbind(seq(1, ncol(cpsfinal), by = 1), apply(cpsfinal, 2, function(x) sample(x, 1)))
                colOrdered <- as.numeric(colRanked[sample(1:length(colRanked[, 2])), 1])
            } else {
                colRanked <- if (directMin) {
                  cbind(seq(1, ncol(cpsfinal), by = 1), apply(cpsfinal, 2, function(x) min(x, na.rm = TRUE)))
                } else {
                  cbind(seq(1, ncol(cpsfinal), by = 1), apply(cpsfinal, 2, function(x) max(x, na.rm = TRUE)))
                }
                colOrdered <- if (directMin) {
                  as.numeric(colRanked[order(colRanked[, 2]), 1])
                } else {
                  as.numeric(colRanked[order(colRanked[, 2], decreasing = TRUE), 1])
                }
            }
            
            rowIndex <- seq(1, nrow(cpsfinal), by = 1)
            
            if (scramble) {
                repeat {
                  codonCol <- as.numeric(colOrdered[1])
                  
                  if (toRestrictSeqs) {
                    innerTrackRE <- 0
                    rowIndexRE <- rowIndex
                    FoundBadRE <- numeric()
                    if (complementary) {
                      FoundBadRErevomp <- numeric()
                    }
                    innerTrackRow <- numeric()
                    isSame <- logical()
                    repeat {
                      innerTrackRE <- innerTrackRE + 1
                      if (length(rowIndexRE) == 0) {
                        if (trueRandom) {
                          codonRow <- as.numeric(sample(rowIndex, 1))
                          
                        } else {
                          codonRow <- ifelse(directMin, {
                            as.numeric(rowIndex[which.min(cpsfinal[rowIndex, colOrdered[1]])])
                          }, {
                            as.numeric(rowIndex[which.max(cpsfinal[rowIndex, colOrdered[1]])])
                          })
                        }
                        (break)()
                      }
                      
                      if (trueRandom) {
                        codonRow <- as.numeric(sample(rowIndexRE, 1))
                        
                      } else {
                        codonRow <- ifelse(directMin, {
                          as.numeric(rowIndexRE[which.min(cpsfinal[rowIndexRE, colOrdered[1]])])
                        }, {
                          as.numeric(rowIndexRE[which.max(cpsfinal[rowIndexRE, colOrdered[1]])])
                        })
                      }
                      
                      isSame[innerTrackRE] <- seqVessel[positionA[codonCol]] == seqVesselORF[positionA[codonRow]] | seqVessel[positionA[codonRow]] == seqVesselORF[positionA[codonCol]]
                      
                      seqVesselHolder <- seqVessel
                      seqVesselHolder[c(positionA[codonRow], positionA[codonCol])] <- seqVesselHolder[c(positionA[codonCol], positionA[codonRow])]
                      seqRECre <- paste0(seqVesselHolder, collapse = "")
                      FoundBadRE[innerTrackRE] <- sum(sapply(1:length(notAllowedRE), function(x) {
                        any(grepl(notAllowedRE[x], seqRECre))
                      }))
                      
                      if (complementary) {
                        seqrevcomp <- revcomp(seqRECre)
                        FoundBadRErevomp[innerTrackRE] <- sum(sapply(1:length(notAllowedRE), function(x) {
                          any(grepl(notAllowedRE[x], seqrevcomp))
                        }), na.rm = TRUE)
                        FoundBadRE <- FoundBadRE + FoundBadRErevomp
                      }
                      innerTrackRow[innerTrackRE] <- codonRow
                      if (FoundBadRE[innerTrackRE] == 0 & !isSame[innerTrackRE]) {
                        (break)()
                      }
                      rowIndexRE <- rowIndexRE[-which(rowIndexRE == codonRow)]
                      (next)()
                    }
                    
                    codonRows <- innerTrackRow[which(FoundBadRE == min(FoundBadRE, na.rm = TRUE))]
                    isSameCodonRows <- isSame[which(FoundBadRE == min(FoundBadRE, na.rm = TRUE))]
                    codonRows <- codonRows[!isSameCodonRows]
                    if (length(codonRows) == 0) {
                      codonRow <- innerTrackRow[which(FoundBadRE == min(FoundBadRE, na.rm = TRUE))[1]]
                    } else {
                      codonRow <- codonRows[1]
                    }
                    seqVessel[c(positionA[codonRow], positionA[codonCol])] <- seqVessel[c(positionA[codonCol], positionA[codonRow])]
                    
                  } else {
                    
                    rowIndexScramble <- rowIndex
                    repeat {
                      
                      if (length(rowIndexScramble) == 0) {
                        if (trueRandom) {
                          codonRow <- as.numeric(sample(rowIndex, 1))
                          
                        } else {
                          codonRow <- ifelse(directMin, {
                            as.numeric(rowIndex[which.min(cpsfinal[rowIndex, colOrdered[1]])])
                          }, {
                            as.numeric(rowIndex[which.max(cpsfinal[rowIndex, colOrdered[1]])])
                          })
                        }
                        (break)()
                      }
                      
                      if (trueRandom) {
                        codonRow <- as.numeric(sample(rowIndexScramble, 1))
                        
                      } else {
                        codonRow <- ifelse(directMin, {
                          as.numeric(rowIndexScramble[which.min(cpsfinal[rowIndexScramble, colOrdered[1]])])
                        }, {
                          as.numeric(rowIndexScramble[which.max(cpsfinal[rowIndexScramble, colOrdered[1]])])
                        })
                      }
                      
                      
                      isSame <- seqVessel[positionA[codonCol]] == seqVesselORF[positionA[codonRow]] | seqVessel[positionA[codonRow]] == seqVesselORF[positionA[codonCol]]
                      
                      if (isSame) {
                        rowIndexScramble <- rowIndexScramble[-which(rowIndexScramble == codonRow)]
                        (next)()
                      } else {
                        seqVessel[c(positionA[codonRow], positionA[codonCol])] <- seqVessel[c(positionA[codonCol], positionA[codonRow])]
                      }
                      (break)()
                    }
                  }
                  
                  colOrdered <- colOrdered[-which(colOrdered == codonRow)]
                  rowIndex <- rowIndex[-which(rowIndex == codonCol)]
                  if (codonRow != codonCol) {
                    colOrdered <- colOrdered[-which(colOrdered == codonCol)]
                    rowIndex <- rowIndex[-which(rowIndex == codonRow)]
                  }
                  if (length(colOrdered) == 0) {
                    (break)()
                  }
                }
                
            } else {
                repeat {
                  
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
                        if (trueRandom) {
                          codonRow <- as.numeric(sample(rowIndex, 1))
                          
                        } else {
                          codonRow <- ifelse(directMin, {
                            as.numeric(rowIndex[which.min(cpsfinal[rowIndex, colOrdered[1]])])
                          }, {
                            as.numeric(rowIndex[which.max(cpsfinal[rowIndex, colOrdered[1]])])
                          })
                        }
                        (break)()
                      }
                      
                      if (trueRandom) {
                        codonRow <- as.numeric(sample(rowIndexRE, 1))
                        
                      } else {
                        codonRow <- ifelse(directMin, {
                          as.numeric(rowIndexRE[which.min(cpsfinal[rowIndexRE, colOrdered[1]])])
                        }, {
                          as.numeric(rowIndexRE[which.max(cpsfinal[rowIndexRE, colOrdered[1]])])
                        })
                      }
                      
                      
                      seqVesselHolder <- seqVessel
                      seqVesselHolder[c(positionA[codonRow], positionA[codonCol])] <- seqVesselHolder[c(positionA[codonCol], positionA[codonRow])]
                      seqRECre <- paste0(seqVesselHolder, collapse = "")
                      FoundBadRE[innerTrackRE] <- sum(sapply(1:length(notAllowedRE), function(x) {
                        any(grepl(notAllowedRE[x], seqRECre))
                      }))
                      
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
                    if (trueRandom) {
                      codonRow <- as.numeric(sample(rowIndex, 1))
                    } else {
                      codonRow <- ifelse(directMin, {
                        as.numeric(rowIndex[which.min(cpsfinal[rowIndex, colOrdered[1]])])
                      }, {
                        as.numeric(rowIndex[which.max(cpsfinal[rowIndex, colOrdered[1]])])
                      })
                    }
                    
                    seqVessel[c(positionA[codonRow], positionA[codonCol])] <- seqVessel[c(positionA[codonCol], positionA[codonRow])]
                    
                  }
                  
                  colOrdered <- colOrdered[-which(colOrdered == codonRow)]
                  rowIndex <- rowIndex[-which(rowIndex == codonCol)]
                  if (codonRow != codonCol) {
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
        
        if (!is.null(draw)) {
            Sys.sleep(0.1)
        }
        
        tracking <- CPScalcOne(seqREC, reference, windowSize = windowSize, draw = draw)
        repeati <- repeati + 1
        repeatScore[repeati] <- as.numeric(tracking)
        repeatMutations[repeati] <- as.numeric(table(strsplit(seqORF, "")[[1]] == strsplit(seqREC, "")[[1]])["FALSE"])
        repeatSeq[[repeati]] <- seqREC
        if (!is.null(restrictSeqs)) {
            trackREsum[repeati] <- sum(sapply(1:length(notAllowedRE), function(x) {
                any(grepl(notAllowedRE[x], seqREC))
            }))
        }
        
        if (!silent) {
            setTxtProgressBar(pb, repeati)
        }
        
        if (scoreGiven) {
            
            switch((repeati - (totalcycles/5)), {
                Trandom <- 1
                Tbest <- 0.1
            })
            switch((repeati - (totalcycles/2)), {
                Trandom <- 0.1
                Tbest <- 1
            })
            switch((repeati - (totalcycles/1.5)), {
                Trandom <- 0.01
                Tbest <- 10
            })
            switch((repeati - (totalcycles/1.1)), {
                Trandom <- 0
                Tbest <- 100
            })
            if (any(repeati == stepladder)) {
                scoremean <- mean(repeatScore[(length(repeatScore) - stepsize):length(repeatScore)], na.rm = TRUE)
                directMin <- if (scoremean > desiredScore) {
                  TRUE
                } else if (scoremean < desiredScore) {
                  FALSE
                }
            }
            roundscores <- round(repeatScore, 5)
            thresholdDistance <- abs(roundscores - desiredScore)
            heavyScores <- 1/(thresholdDistance^2 + 1)
            if (returnBest & repeati > initialstep) {
                if (toRestrictSeqs) {
                  trackREindex <- seq(1, length(trackREsum), by = 1)
                  trackREindex <- trackREindex[which(trackREsum == min(trackREsum, na.rm = TRUE))]
                  REheavyscores <- heavyScores[trackREindex]
                  seqREC <- repeatSeq[[as.numeric(trackREindex[which.max(REheavyscores)])]]
                } else {
                  seqREC <- repeatSeq[[which.max(heavyScores)]]
                }
            }
            
            switch((repeati - (totalcycles - 1)), {
                (break)()
            })
        } else {
            if (returnBest) {
                if (toRestrictSeqs) {
                  trackREindex <- seq(1, length(trackREsum), by = 1)
                  
                  trackREindex <- trackREindex[which(trackREsum == min(trackREsum, na.rm = TRUE))]
                  RErepeatScore <- repeatScore[trackREindex]
                  
                  seqREC <- ifelse(directMin, {
                    repeatSeq[[as.numeric(trackREindex[which.min(RErepeatScore)])]]
                  }, {
                    repeatSeq[[as.numeric(trackREindex[which.max(RErepeatScore)])]]
                  })
                } else {
                  seqREC <- ifelse(directMin, {
                    repeatSeq[[which.min(repeatScore)]]
                  }, {
                    repeatSeq[[which.max(repeatScore)]]
                  })
                }
            }
            switch((repeati - (totalcycles/5)), {
                Trandom <- 0.5
                Tbest <- 0.1
            })
            switch((repeati - (totalcycles/2)), {
                Trandom <- 0.1
                Tbest <- 1
            })
            switch((repeati - (totalcycles/1.5)), {
                Trandom <- 0.01
                Tbest <- 10
            })
            switch((repeati - (totalcycles/1.1)), {
                Trandom <- 0
                Tbest <- 100
            })
            switch((repeati - (totalcycles - 1)), {
                break
            })
        }
        if (seqLength < 1300) {
            Sys.sleep(time = 0.2)
        }
    }
    if (!silent) {
        close(pb)
    }
    
    
    standardScores <- if (directMin) {
        (repeatScore - max(repeatScore))/sd(repeatScore)
    } else {
        (repeatScore - min(-(repeatScore)))/sd(-repeatScore)
    }
    if (scoreGiven) {
        standardheavy <- (heavyScores - max(heavyScores))/sd(heavyScores)
    }
    
    if (!is.null(restrictSeqs)) {
        trackREindex <- seq(1, length(trackREsum), by = 1)
        trackREindex <- trackREindex[which(trackREsum == min(trackREsum))]
        if (scoreGiven) {
            REheavyscores <- heavyScores[trackREindex]
        } else {
            RErepeatScore <- repeatScore[trackREindex]
        }
    }
    
    if (maxmutations) {
        trackMutindex <- seq(1, length(repeatMutations), by = 1)
        maxMuts <- max(repeatMutations, na.rm = TRUE)
        trackMutindex <- trackMutindex[which(repeatMutations == maxMuts)]
        if (!is.null(restrictSeqs)) {
            if (any(trackREsum == 0)) {
                bestRE <- which(trackREsum == 0)
            } else {
                minRE <- min(trackREsum, na.rm = TRUE)
                bestRE <- which(trackREsum == minRE)
            }
            trackMutindexCheck <- trackMutindex[trackREsum]
            if (length(trackMutindexCheck) != 0) {
                trackMutindex <- trackMutindexCheck
            } else {
                warning("Maximizing number of mutations in scramble recoding overrides sequence restriction, some sequences have not been removed.")
            }
        }
        if (scoreGiven) {
            MUTheavyscores <- heavyScores[trackMutindex]
        } else {
            MUTrepeatScore <- repeatScore[trackMutindex]
        }
    }
    
    if (scoreGiven) {
        if (maxmutations) {
            seqREC <- repeatSeq[[as.numeric(trackMutindex[which.max(MUTheavyscores)])]]
            
        } else if (!maxmutations) {
            if (!is.null(restrictSeqs)) {
                seqREC <- repeatSeq[[as.numeric(trackREindex[which.max(REheavyscores)])]]
            } else {
                seqREC <- repeatSeq[[which.max(heavyScores)]]
            }
        }
    } else {
        if (maxmutations) {
            
            if (directMin) {
                seqREC <- repeatSeq[[as.numeric(trackMutindex[which.min(MUTrepeatScore)])]]
            } else {
                seqREC <- repeatSeq[[as.numeric(trackMutindex[which.max(MUTrepeatScore)])]]
            }
            
        } else if (!maxmutations) {
            if (!is.null(restrictSeqs)) {
                trackREindex <- seq(1, length(trackREsum), by = 1)
                trackREindex <- trackREindex[trackREsum == min(trackREsum)]
                RErepeatScore <- repeatScore[trackREindex]
                seqREC <- ifelse(directMin, {
                  repeatSeq[[as.numeric(trackREindex[which.min(RErepeatScore)])]]
                }, {
                  repeatSeq[[as.numeric(trackREindex[which.max(RErepeatScore)])]]
                })
            } else {
                seqREC <- if (directMin) {
                  repeatSeq[[which.min(repeatScore)]]
                } else {
                  repeatSeq[[which.max(repeatScore)]]
                }
            }
        }
    }
    
    # override return seq with last shuffle if true random
    if (trueRandom) {
        seqREC <- repeatSeq[[length(repeatSeq)]]
        if (scramble) {
            warning("Using scramble together with random will not generate a true random sequence.")
        }
    }
    
    finalCPS <- CPScalcOne(seqREC, reference, windowSize = windowSize, draw = draw)
    
    compare.codons <- countCodons(seqREC, transTable = transTable) == countCodons(seqORF, transTable = transTable)
    
    if (any(!compare.codons)) {
        numberofCodonChanges <- as.numeric(table(compare.codons)["FALSE"])
        message("Error: Codon usage was changed.")
    } else {
        numberofCodonChanges <- 0
    }
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
        CPSplot(seqOne = seqREC, refOne = reference, seqTwo = seqORF, refTwo = NULL, start = 1, end = NULL, windowSize = windowSize)
    }
    compare.codons <- countCodons(seqREC, transTable = transTable) == countCodons(seqORF, transTable = transTable)
    if (any(!compare.codons)) {
        numberofCodonChanges <- as.numeric(table(compare.codons)["FALSE"])
    } else {
        numberofCodonChanges <- 0
    }
    numberofMutations <- as.numeric(table(strsplit(seqORF, "")[[1]] == strsplit(seqREC, "")[[1]])["FALSE"])
    oldCPSarray <- CPScalcOne(seqORF, reference, average = FALSE, windowSize = windowSize)
    newCPSarray <- CPScalcOne(seqREC, reference, average = FALSE, windowSize = windowSize)
    if (is.null(restrictSeqs)) {
        ReturnOne <- list(startingCPS, finalCPS, numberofCodonChanges, numberofMutations, oldCPSarray, newCPSarray, seqREC)
        names(ReturnOne) <- c("oldCPS", "newCPS", "codonchanges", "mutations", "oldCPSarray", "newCPSarray", "returnSeq")
    } else {
        reSeqs <- sum(sapply(1:length(notAllowedRE), function(x) grepl(notAllowedRE[x], seqREC)), na.rm = TRUE)
        ReturnOne <- list(startingCPS, finalCPS, numberofCodonChanges, numberofMutations, oldCPSarray, newCPSarray, seqREC, reSeqs)
        names(ReturnOne) <- c("oldCPS", "newCPS", "codonchanges", "mutations", "oldCPSarray", "newCPSarray", "returnSeq", "restrSeqs")
    }
    if (!silent) {
        cat("\n Original mean codon pair score: ", ReturnOne$oldCPS, "\n \n", "New mean codon pair score: ", ReturnOne$newCPS, "\n \n", "Number of codon changes: ", 
            ReturnOne$codonchanges, "\n \n", "Number of mutations: ", as.numeric(ReturnOne$mutations), "\n")
    }
    invisible(ReturnOne)
    
} 
