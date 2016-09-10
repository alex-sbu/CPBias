# split by 3 nucleotides every 3 nucleotides for in frame sequence
splitbyCodon <- function(sequence) {
    sequence <- strsplit(sequence, split = "")[[1]]
    if (length(sequence) < 3) {
        return(NA)
    }
    splitsequence <- sapply(seq(from = 1, to = length(sequence), by = 3), function(x) {
        paste0(sequence[x:(x + 2)], collapse = "")
    })
    if ((length(sequence)%%3) == 1 | (length(sequence)%%3) == 2) {
        splitsequence <- splitsequence[1:(length(splitsequence) - 1)]
    }
    splitsequence
}

# split by 6 nucleotides every 3 nucleotides for in frame sequence
splitbyCodonPair <- function(sequence) {
    sequence <- strsplit(sequence, split = "")[[1]]
    if (length(sequence) < 6) {
        return(NA)
    }
    CountOrder <- seq(from = 1, to = (length(sequence) - 2), by = 3)
    splitsequence <- sapply(CountOrder, function(x) {
        paste0(sequence[x:(x + 5)], collapse = "")
    })
    splitsequence <- splitsequence[1:(length(splitsequence) - 1)]
    splitsequence
}

# count number of occurances per codon
countCodons <- function(sequence, transTable = standardTranslation) {
    if (is.na(sequence)) {
        return(NA)
    }
    CPSblanktable <- makeCPBtemplate(transTable)
    splitsequence <- splitbyCodon(sequence)
    UniquecodonCounts <- table(factor(splitsequence, levels = unique(CPSblanktable[, 1])))
    UniquecodonCounts
}

## translate sequence
translateSeq <- function(sequence, templateCPStable = standardTranslation) {
    sequence <- splitbyCodon(sequence)
    TranslatedString <- sapply(sequence, function(x) {
        tolower(as.character((templateCPStable[match(x, templateCPStable[, 1]), 4])))
    })
    PepSeq <- as.character(na.omit(TranslatedString))
    PepSeq
}

# Generate a template of codons and codon pairs for recoding and generating CPBtable
makeCPBtemplate <- function(transTable = standardTranslation) {
    transTable <- standardTranslation
    codons <- as.character(tolower(transTable[, 2]))
    stopCodons <- which(grepl("stop", codons))
    codons <- transTable[-stopCodons, 1]
    codons <- tolower(as.character(codons))
    aminoacids <- as.character(tolower(transTable[, 2]))
    aminoacids <- aminoacids[-stopCodons]
    codonPairs <- expand.grid(codons, codons, stringsAsFactors = FALSE)
    codonPairsC <- apply(codonPairs, 1, function(x) paste0(x, collapse = ""))
    aminoOne <- sapply(1:length(codonPairs[, 1]), function(x) aminoacids[match(codonPairs[x, 1], codons)])
    aminoTwo <- sapply(1:length(codonPairs[, 2]), function(x) aminoacids[match(codonPairs[x, 2], codons)])
    aminoPairs <- data.frame(aminoOne, aminoTwo, stringsAsFactors = FALSE)
    aminoPairsC <- apply(aminoPairs, 1, function(x) paste0(x, collapse = ""))
    returnTable <- data.frame(codonPairs[, 1], codonPairs[, 2], codonPairsC, aminoOne, aminoTwo, aminoPairsC, stringsAsFactors = FALSE)
    names(returnTable) <- c("codon1", "codon2", "codonpair", "amino1", "amino2", "aminopair")
    returnTable
}

# Write and export fasta file
exportFasta <- function(sequence, name, location = NULL) {
    if (is.null(location)) {
        location <- getwd()
    }
    
    if (any(grepl(".fasta", name))) {
        name <- gsub(".fasta", "", name)
    }
    
    namedLoc <- file.path(location, paste0(name, ".fasta"))
    outputFASTA <- file(namedLoc)
    FastaName <- paste0("> ", name, "\n")
    sequence <- strsplit(sequence, split = "")[[1]]
    inputLength <- length(sequence)
    seqLength <- seq(from = 1, to = length(sequence), by = 80)
    preFinalString <- sapply(seqLength[1:length(seqLength) - 1], function(x) {
        paste0(paste0(sequence[x:(x + 79)], collapse = ""), "\n")
    })
    ConStrings <- paste0(preFinalString, collapse = "")
    toAdd <- length(sequence) - seqLength[length(seqLength)]
    FinalString <- paste0(FastaName, ConStrings, paste0(sequence[seqLength[length(seqLength)]:(seqLength[length(seqLength)] + toAdd)], collapse = ""), collapse = "")
    writeLines(FinalString, outputFASTA)
    close(outputFASTA)
}

# Read and import fasta file
##' @description
##' A fasta file containing single or multiple sequences is imported as a lowercase character string, or a vector of multiple character strings.
##'
##' @details
##' Comments marked by `;' are ignored. If importing sequences for use in \code{\link{CPBtable}}, make sure to set \code{sepSequences} to TRUE otherwise sequences not divisible by three might shift the read frame.
##'
##' @title Import nucleotide sequence from fasta file
##'
##' @param fastaLocation Location of fasta file.
##' @param sepSequences If FALSE fasta files containing multiple sequences are concatenated into one string, otherwise a vector of strings is returned.
##'
##' @return
##' \item{seqs}{Sequence(s).}
##' \item{seqNames}{Names of all sequences.}
##'
##' @examples
##' # Tick borne encephalitis virus NS5 coding region
##' fastaLocation <- system.file('tbevns5.fasta', package = 'CPBias')
##' tbev <- importFasta(fastaLocation)
##' tbev[[2]]
##'
##' # Randomly selected protein coding sequences from the human
##' # CCDS data set.
##' fastaLocation <- system.file('ccds.fasta', package = 'CPBias')
##' ccds <- importFasta(fastaLocation, sepSequences = TRUE)
##' ccds[[2]]
importFasta <- function(fastaLocation, sepSequences = FALSE) {
    rawFasta <- readLines(fastaLocation)
    seqNames <- rawFasta[grep(">.+$", rawFasta)]
    seqNames <- gsub(">", "", seqNames)
    
    if (sepSequences) {
        rawFasta <- gsub(">.+$", ",", rawFasta)
        rawFasta <- gsub(";.+$", "", rawFasta)
        rawFasta <- paste0(rawFasta, collapse = "")
        rawFasta <- gsub("^,", "", rawFasta)
        rawFasta <- strsplit(rawFasta, split = ",")[[1]]
    } else {
        rawFasta <- gsub(">.+$", "", rawFasta)
        rawFasta <- gsub(";.+$", "", rawFasta)
        rawFasta <- paste0(rawFasta, collapse = "")
    }
    
    Sequence <- tolower(rawFasta)
    seqList <- list(Sequence, seqNames)
    names(seqList) <- c("seqs", "seqNames")
    seqList
}

##'@description
##' Spearman rank correlation test for comparing two CPB reference tables
##'@details
##' Uses \code{cor.test} R base function with the spearman method. Depending on the CPB's, it might be useful to standardize the data first.
##'
##' @title Codon pair bias correlation test
##' @param refOne First CPB reference table
##' @param refTwo Second CPB reference table
##' @return Output from \code{cor.test}, and a scatter plot.
##' @examples
##' CPBcorr(Homo.sapiens, Aedes.aegypti)
CPBcorr <- function(refOne, refTwo) {
    CPBcor <- cor.test(refOne[, 2], refTwo[, 2], method = "spearman")
    plot(refOne[, 2] ~ refTwo[, 2], pch = 1, col = "#2e3436", type = "p")
    return(CPBcor)
}

# generate reverse complement sequence
revcomp <- function(seq) {
    seq <- tolower(seq)
    reco <- seq
    seq <- strsplit(seq, "")[[1]]
    reco <- strsplit(reco, "")[[1]]
    reco[which(seq == "a")] <- "t"
    reco[which(seq == "t")] <- "a"
    reco[which(seq == "c")] <- "g"
    reco[which(seq == "g")] <- "c"
    reco[which(seq == "r")] <- "y"
    reco[which(seq == "y")] <- "r"
    reco[which(seq == "m")] <- "k"
    reco[which(seq == "k")] <- "m"
    reco[which(seq == "b")] <- "v"
    reco[which(seq == "v")] <- "b"
    reco[which(seq == "d")] <- "h"
    reco[which(seq == "h")] <- "d"
    reco <- paste0(rev(reco), collapse = "")
} 
