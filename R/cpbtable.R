##' @description
##' Calculate the observed to expected frequency of all codon pairs for a given set of protein coding gene sequences.
##'
##' @details
##'There are 3,721 coding codon pairs, if using a standard translation table. The score (CPS) of each individual codon pair is determined by,
##'\deqn{\ln \left( \frac{{\textrm{codon pair$_{ab}$} \times (\textrm{amino acid$_{a}$} \times \textrm{amino acid$_{b}$})}}{{(\textrm{codon$_{a}$} \times \textrm{codon$_{b}$}) \times \textrm{amino acid pair$_{ab}$}}} \right)}{ln(codon pair[ab] x (amino acid[a] x amino acid[b]))/((codon[a] x codon[b]) x amino acid pair[ab])}
##'Each value is measured as the relative frequency of the total. Tandem codon positions are marked \emph{a} and \emph{b}. A codon pair consists of 6 nucleotides, and counting is every three nucleotides along the sequence.
##'
##'Sequences containing nucleotides not found in the translation table and sequences not divisible by three are excluded. Codon pairs containing codons undefined in the translation table, and codon pairs containing stop codons, will generate NA's and those codon pairs will not be included in the CPB calculation. By default, the standard translation table is used (see \code{\link{standardTranslation}}). All input sequences should be in frame, protein coding (CDS) sequences.
##'
##' @note
##' Codon pair bias, also called codon context, is typically regarded to be species specific. To this end CPB reference tables have been calculated for organisms which have whole genome CDS sequence data available. See \code{\link{listCPB}} for a list of pre-calculated CPB reference tables.
##'
##' @title Determine the codon pair bias
##'
##' @param sequences Input can be the location of a fasta file, or a character string vector.
##' @param dnfControl If TRUE, III-I dinucleotide bias is factored out.
##' @param transTable Translation table to use for identifying codons. See \code{\link{standardTranslation}}.
##' @param name Name the table if save is TRUE.
##' @param save TRUE or FALSE to save the CPB reference table as a comma delimited csv file.
##' @param location File path to save the csv file.
##' @param silent If TRUE the progress bar is suppressed.
##'
##' @return
##' A list with two elements is returned invisibly:
##' \item{CPBtable}{CPB reference table containing all coding codon pairs and their individual CPS calculated with the above formula. This format can be used as the CPB reference input to the CPSdesign functions.}
##' \item{complete.CPBtable}{Larger CPBtable containing the frequencies of all components of the codon pair score calculation.}
##'
##' @references Coleman JR, et al. 2008 Virus attenuation by genome-scale changes in codon pair bias. \emph{Science} \bold{320}(5884):1784--1787.
##'
##' @examples
##' fastaLocation <- system.file('ccds.fasta', package = 'CPBias')
##' ccds <- importFasta(fastaLocation, sepSequences = TRUE)[[1]]
##'
##' ccdsCPB <- CPBtable(ccds)
##' # First element in returned list is the  CPB reference table
##' ccds.sample <- ccdsCPB[[1]]
##'
##' # CPBtable will import sequences automatically if fasta location is given
##' ccdsCPB <- CPBtable(fastaLocation)[[1]]
##' head(ccdsCPB)
##'
##' # Factor out dinucleotide bias between codons
##' dnCPB <- CPBtable(fastaLocation, dnfControl=TRUE)[[1]]
##' plot(ccdsCPB[,2], dnCPB[,2])
CPBtable <- function(sequences, dnfControl = FALSE, transTable = standardTranslation, name = NULL, save = FALSE, location = NULL, silent = FALSE) {
    
    ## Generate output table name
    if (any(grepl(".fasta", sequences)) & length(sequences) == 1) {
        genome <- importFasta(sequences, sepSequences = TRUE)[[1]]
        if (is.null(name)) {
            name <- gsub(".fasta", "", sequences)
        }
    } else {
        genome <- tolower(sequences)
        if (is.null(name)) {
            name <- "CPBtable"
        }
    }
    
    ## Sequence quality check
    CPSblanktable <- makeCPBtemplate(transTable)
    nt <- paste0("[^", paste0(unique(strsplit(paste0(CPSblanktable$codon1, collapse = ""), split = "")[[1]]), collapse = ""), "]")
    
    if (any(grepl(nt, genome))) {
        removeNA <- grep(nt, genome)
        warning(paste0("Sequence(s) removed for containing unknown nucleotides: ", length(removeNA)))
        genome <- genome[-removeNA]
    }
    
    holdchar <- numeric()
    for (i in 1:length(genome)) {
        holdchar[i] <- nchar(genome[i])%%3
    }
    
    if (any(holdchar != 0)) {
        removeNMT <- which(holdchar != 0)
        warning(paste0("Sequence(s) removed because sequence length was not a multiple of three: ", length(removeNMT)))
        genome <- genome[-removeNMT]
    }
    
    if (length(genome) < 1) {
        stop("No valid sequences found.")
    }
    
    ## Codon and codon pair raw counts
    if (!silent) {
        CPStableCounter <- 0
        pb <- txtProgressBar(min = 0, max = length(genome) * 2, initial = 0, char = ".", width = (getOption("width")/2), style = 3)
    }
    
    AllCP <- lapply(1:length(genome), function(x) {
        
        if (!silent) {
            CPStableCounter <<- CPStableCounter + 1
            setTxtProgressBar(pb, CPStableCounter)
        }
        
        splitbyCodonPair(genome[x])
    })
    
    AllCP <- unlist(AllCP)
    AllCPSUMs <- table(factor(AllCP, levels = unique(CPSblanktable$codonpair)))
    
    AllC <- lapply(1:length(genome), function(x) {
        
        if (!silent) {
            CPStableCounter <<- CPStableCounter + 1
            setTxtProgressBar(pb, CPStableCounter)
        }
        
        splitbyCodon(genome[x])
    })
    
    AllC <- unlist(AllC)
    AllCSUMs <- table(factor(AllC, levels = unique(CPSblanktable$codon1)))
    
    if (!silent) {
        close(pb)
    }
    
    ## Calculate III-I dinucleotide observed:expected frequencies for each codon pair
    if (dnfControl) {
        AlldnIII.I <- substring(AllCP, 3, 4)
        AlldnIII <- substring(AllCP, 3, 3)
        AlldnI <- substring(AllCP, 4, 4)
        dnCPSblanktable <- substring(CPSblanktable$codonpair, 3, 4)
        dnIIICPSblanktable <- substring(CPSblanktable$codonpair, 3, 3)
        dnICPSblanktable <- substring(CPSblanktable$codonpair, 4, 4)
        AlldnIII.ISUMs <- table(factor(AlldnIII.I, levels = unique(AlldnIII.I)))
        AlldnIIISUMs <- table(factor(AlldnIII, levels = unique(AlldnIII)))
        AlldnISUMs <- table(factor(AlldnI, levels = unique(AlldnI)))
        matchdnIII <- sapply(match(dnIIICPSblanktable, names(AlldnIIISUMs)), function(x) AlldnIIISUMs[x])
        matchdnI <- sapply(match(dnICPSblanktable, names(AlldnISUMs)), function(x) AlldnISUMs[x])
        matchdnIII.I <- sapply(match(dnCPSblanktable, names(AlldnIII.ISUMs)), function(x) AlldnIII.ISUMs[x])
        totaldnIII <- sum(as.numeric(AlldnIIISUMs))
        totaldnI <- sum(as.numeric(AlldnISUMs))
        totaldnIII.I <- sum(as.numeric(AlldnIII.ISUMs))
        dnIII <- as.numeric(matchdnIII)/totaldnIII
        dnI <- as.numeric(matchdnI)/totaldnI
        dnIII.Iexp <- dnIII * dnI
        dnIII.Iobs <- as.numeric(matchdnIII.I)/totaldnIII.I
        dnIII.Ifreq <- dnIII.Iobs/dnIII.Iexp
    }
    
    ## Calculate codon pair observed:expected frequencies while correcting for amino acid bias
    matchCone <- sapply(match(CPSblanktable$codon1, names(AllCSUMs)), function(x) AllCSUMs[x])
    matchCtwo <- sapply(match(CPSblanktable$codon2, names(AllCSUMs)), function(x) AllCSUMs[x])
    matchCP <- sapply(match(CPSblanktable$codonpair, names(AllCPSUMs)), function(x) AllCPSUMs[x])
    
    AAonecount <- sapply(unique(CPSblanktable$amino1), function(x) sum(matchCone[CPSblanktable$amino1 == x]))
    matchAAone <- sapply(match(CPSblanktable$amino1, unique(CPSblanktable$amino1)), function(x) AAonecount[x])
    AAtwocount <- sapply(unique(CPSblanktable$amino2), function(x) sum(matchCtwo[CPSblanktable$amino2 == x]))
    matchAAtwo <- sapply(match(CPSblanktable$amino2, unique(CPSblanktable$amino2)), function(x) AAtwocount[x])
    AAPcount <- sapply(unique(CPSblanktable$aminopair), function(x) sum(matchCP[CPSblanktable$aminopair == x]))
    matchAAP <- sapply(match(CPSblanktable$aminopair, unique(CPSblanktable$aminopair)), function(x) AAPcount[x])
    
    totalC <- sum(as.numeric(AllCSUMs))
    uniqueAA <- matchAAone[!duplicated(CPSblanktable$amino1)]
    totalAA <- sum(as.numeric(uniqueAA))
    totalCP <- sum(as.numeric(matchCP))
    uniqueAAP <- matchAAP[!duplicated(CPSblanktable$aminopair)]
    totalAAP <- sum(as.numeric(uniqueAAP))
    
    Cone <- as.numeric(matchCone)/totalC
    Ctwo <- as.numeric(matchCtwo)/totalC
    CPexp <- Cone * Ctwo
    CPobs <- as.numeric(matchCP)/totalCP
    
    AAone <- as.numeric(matchAAone)/totalAA
    AAtwo <- as.numeric(matchAAtwo)/totalAA
    AAPexp <- AAone * AAtwo
    AAPobs <- as.numeric(matchAAP)/totalAAP
    
    CPfinal <- (CPobs/AAPobs)/(CPexp/AAPexp)
    
    if (dnfControl) {
        beforeCPfinal <- CPfinal
        CPfinal <- CPfinal/dnIII.Ifreq
    }
    
    ## Natural log of final codon pair frequencies
    CPfinal <- log(CPfinal)
    
    
    ## Make CPBtable output (small table)
    outputCPB <- data.frame(CPSblanktable$codonpair, CPfinal, stringsAsFactors = FALSE)
    
    if (any(outputCPB[, 2] == "-Inf")) {
        outputCPB[, 2][which(outputCPB[, 2] == "-Inf")] <- NA
    }
    
    outputCPB <- outputCPB[order(outputCPB[, 1]), ]
    names(outputCPB) <- c("codon pair", "cps")
    name <- tolower(name)
    
    ## Save CPBtable as csv
    if (save) {
        
        if (is.null(location)) {
            location <- getwd()
        }
        
        write.csv(outputCPB, file = file.path(location, paste0(name, ".csv")), row.names = FALSE)
        
    }
    
    ## Make large CPBtable output with all raw data used in the calculations
    
    if (dnfControl) {
        largeTable <- data.frame(CPSblanktable$codon1, CPSblanktable$amino1, CPSblanktable$codon2, CPSblanktable$amino2, CPSblanktable$aminopair, CPSblanktable$codonpair, 
            dnIIICPSblanktable, dnICPSblanktable, dnCPSblanktable, Cone, AAone, Ctwo, AAtwo, AAPexp, AAPobs, CPexp, CPobs, beforeCPfinal, dnIII, dnI, dnIII.Iexp, 
            dnIII.Iobs, dnIII.Ifreq, CPfinal, stringsAsFactors = FALSE)
        names(largeTable) <- c("codon1", "amino1", "codon2", "amino2", "aminopair", "codonpair", "nucleotide.III", "nucleotide.I", "dinucleotide.III.I", "codon1 frequency", 
            "amino1 frequency", "codon2 frequency", "amino2 frequency", "amino acid pair expected", "amino acid pair observed", "codon pair expected", "codon pair observed", 
            "codon pair bias (observed:expected)", "nucleotide III frequency", "nucleotide I frequency", "dinucleotide III.I expected", "dinucleotide III.I observed", 
            "dinucleotide III.I bias (observed:expected)", "natural log of codon pair bias:dinucleotide III.I bias")
    } else {
        largeTable <- data.frame(CPSblanktable$codon1, CPSblanktable$amino1, CPSblanktable$codon2, CPSblanktable$amino2, CPSblanktable$aminopair, CPSblanktable$codonpair, 
            Cone, AAone, Ctwo, AAtwo, AAPexp, AAPobs, CPexp, CPobs, CPfinal, stringsAsFactors = FALSE)
        names(largeTable) <- c("codon1", "amino1", "codon2", "amino2", "aminopair", "codonpair", "codon1 frequency", "amino1 frequency", "codon2 frequency", "amino2 frequency", 
            "amino acid pair expected", "amino acid pair observed", "codon pair expected", "codon pair observed", "natural log of (codon pair observed:expected)")
    }
    
    ## Return invisible list with both data frames
    CPBreturnList <- list(outputCPB, largeTable)
    names(CPBreturnList) <- c("CPBtable", "complete.CPBtable")
    return(invisible(CPBreturnList))
} 
