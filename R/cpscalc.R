##' @description
##' Calculate the codon pair score of a sequence relative to a reference codon pair bias.
##'
##' @title Codon pair score calculator
##'
##' @param sequence Sequences can be input directly as a character string or as the file path to a fasta file.
##' @param reference Reference CPB table. See \code{\link{CPBtable}}.
##' @param start Nucleotide position in the sequence marking the region of the sequence to use for calculation.
##' @param end Nucleotide position in the sequence marking the region of the sequence to use for calculation.
##' @param draw If TRUE, a line plot showing the local CPS along the length of the sequence is output to the graphics device.
##' @param windowSize CPS line plots are smoothed by locally weighted polynomial regression where \code{windowSize} designates the number of nucleotides over which individual codon pair scores are smoothed. If NULL, smoothing spans 7.5\% of the sequence length.
##' @param silent If TRUE the progress bar is suppressed.
##'
##' @return
##' An invisible list is returned:
##' \item{averageCPS}{Average of the individual codon pair scores of the input sequence.}
##' \item{totalCPS}{Numeric vector containing the individual codon pair scores along the sequence length.}
##'
##' @examples
##' fastaLocation <- system.file('tbevns5.fasta', package = 'CPBias')
##' tbev <- importFasta(fastaLocation)[[1]]
##' tbevCPS <- CPScalc(tbev, Homo.sapiens)
##' tbevCPS[[1]]
##'
##' # Can plot a smoothed line plot with CPScalc output as it is done in
##' # CPSplot()
##' CPSx <- 1:length(tbevCPS[[2]])
##' smoothCPS <- loess(tbevCPS[[2]] ~ CPSx, span = (.05), degree = 2)
##' lineCPS <- predict(smoothCPS)
##' plot(lineCPS, type='l')
CPScalc <- function(sequence, reference, start = 1, end = NULL, draw = TRUE, windowSize = NULL, silent = FALSE) {
    
    if (any(grepl(".fasta", sequence))) {
        userGene <- importFasta(sequence)[[1]]
    } else {
        userGene <- tolower(paste0(sequence, collapse = ""))
    }
    
    if (is.null(end)) {
        end <- nchar(userGene)
    }
    userGene <- substr(userGene, start, end)
    userGeneCP <- splitbyCodonPair(userGene)
    userGeneCPS <- sapply(match(userGeneCP, as.character(reference[, 1])), function(x) as.vector(reference[, 2])[x])
    
    userGeneCPSmean <- mean(userGeneCPS, na.rm = TRUE)
    if (draw) {
        CPSplot(seqOne = userGene, refOne = reference, start = start, end = end, windowSize = windowSize)
    }
    
    CPSlist <- list(round(userGeneCPSmean, 5), round(userGeneCPS, 5))
    names(CPSlist) <- c("averageCPS", "totalCPS")
    
    if (!silent) {
        cat("Average codon pair score: ", userGeneCPSmean, "\n")
    }
    
    return(invisible(CPSlist))
}

CPScalcOne <- function(userGene, reference, average = TRUE, blank = FALSE, windowSize = NULL, draw = TRUE) {
    userGeneCP <- splitbyCodonPair(userGene)
    userGeneCPS <- sapply(match(userGeneCP, as.character(reference[, 1])), function(x) as.vector(reference[, 2])[x])
    if (average) {
        userGeneCPSmean <- mean(userGeneCPS, na.rm = TRUE)
        if (draw) {
            rtCPSplot(refs = 1, userGeneCPS, userGeneCPSmean, blank = blank, geneLength = nchar(userGene), windowSize = windowSize)
        }
        return(userGeneCPSmean)
    } else if (!average) {
        return(userGeneCPS)
    }
}

CPScalcTwo <- function(userGene, reference, referenceTwo, blank = FALSE, draw = FALSE, windowSize = NULL) {
    userGeneCP <- splitbyCodonPair(userGene)
    userGeneCPS <- sapply(match(userGeneCP, as.character(reference[, 1])), function(x) as.vector(reference[, 2])[x])
    userGeneCPSmean <- mean(userGeneCPS, na.rm = TRUE)
    userGeneCPStwo <- sapply(match(userGeneCP, as.character(referenceTwo[, 1])), function(x) as.vector(referenceTwo[, 2])[x])
    userGeneCPStwomean <- mean(userGeneCPStwo, na.rm = TRUE)
    if (draw) {
        rtCPSplot(refs = 2, userGeneCPS, userGeneCPSmean, userGeneCPStwo, userGeneCPStwomean, blank = blank, geneLength = nchar(userGene), windowSize = windowSize)
    }
    return(c(userGeneCPSmean, userGeneCPStwomean))
    
} 
