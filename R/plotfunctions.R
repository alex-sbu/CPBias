##' @description
##' Local codon pair scores plotted as a smooth line along the length of a nucleotide sequence.
##'
##' @details
##' The individual codon pair scores across a gene sequence are plotted as a smooth curve by locally weighted quadratic regression (\code{loess}). Span is defined as the window size divided by total sequence length in nucleotides. Two different CPS encodings of a sequence can be compared with up to two reference CPB's.
##'
##' @title Codon pair score line plot
##'
##' @param seqOne Sequences can be input directly as a character string or as the file path to a fasta file. All sequences must be in the correct reading frame, stop codons or codons not defined in the translation table are not allowed and will generate an error.
##' @param refOne CPB reference table. See \code{\link{CPBtable}}.
##' @param seqTwo An alternate codon pair recoding of seqOne. seqTwo only accepts sequences with identical length to seqOne.
##' @param refTwo A second CPB reference. The sequence(s) will be plotted relative to both references.
##' @param start Nucleotide position in the sequence.
##' @param end Nucleotide position in the sequence.
##' @param windowSize Designates the number of nucleotides over which smoothing occurs. If NULL windowSize is 7.5\% of the sequence length.
##'
##' @return Line plot in the graphics device.
##' @examples
##' fastaLocation <- system.file('tbevns5.fasta', package = 'CPBias')
##' tbev <- importFasta(fastaLocation)[[1]]
##'
##' humCPS <- CPScalc(tbev, Homo.sapiens, silent=TRUE, draw=FALSE)[[1]]
##' tbevhmax <- CPSdesign.dual(tbev, Aedes.aegypti, Homo.sapiens, humCPS, 'max',
##' draw=FALSE, silent=TRUE)
##' CPSplot(tbev, Homo.sapiens, tbevhmax[[9]], Homo.sapiens)
CPSplot <- function(seqOne, refOne, seqTwo = NULL, refTwo = NULL, start = 1, end = NULL, windowSize = NULL) {
    if (any(grepl(".fasta", seqOne))) {
        seqI <- importFasta(seqOne)[[1]]
        if (any(grepl(".fasta", seqTwo))) {
            seqII <- importFasta(seqTwo)[[1]]
        } else if (!is.null(seqTwo)) {
            seqII <- tolower(paste0(seqTwo, collapse = ""))
        }
    } else if (!is.null(seqOne)) {
        seqI <- tolower(paste0(seqOne, collapse = ""))
        if (any(grepl(".fasta", seqTwo))) {
            seqII <- importFasta(seqTwo)[[1]]
        } else if (!is.null(seqTwo)) {
            seqII <- tolower(paste0(seqTwo, collapse = ""))
        }
    }
    
    if (is.null(end)) {
        end <- nchar(seqI)
    }
    
    seqlength <- end - start
    lengthbycodon <- floor(seqlength/3)
    
    if (is.null(windowSize)) {
        windowSize <- seqlength * 0.075
    }
    
    xraw <- list()
    yraw <- list()
    ymean <- numeric()
    ysmooth <- list()
    yline <- list()
    listref <- list()
    listseq <- character()
    
    if (!is.null(seqTwo)) {
        if (!is.null(refTwo)) {
            listseq[1] <- seqI
            listseq[2] <- seqII
            listref[[1]] <- refOne
            listref[[2]] <- refTwo
        } else {
            listseq[1] <- seqI
            listseq[2] <- seqII
            listref[[1]] <- refOne
        }
    } else {
        if (!is.null(refTwo)) {
            listseq[1] <- seqI
            listref[[1]] <- refOne
            listref[[2]] <- refTwo
        } else {
            listseq[1] <- seqI
            listref[[1]] <- refOne
        }
    }
    
    linecol <- c("red", "#ef2929", "blue", "#729fcf")
    linetypes <- c(1, 2, 1, 2)
    
    counter <- 0
    for (n in 1:length(listseq)) {
        for (m in 1:length(listref)) {
            counter <- counter + 1
            yraw[[counter]] <- as.numeric(CPScalc(listseq[n], listref[[m]], start, end, draw = FALSE, silent = TRUE)[[2]])
            ymean[counter] <- as.numeric(CPScalc(listseq[n], listref[[m]], start, end, draw = FALSE, silent = TRUE)[[1]])
            xraw[[counter]] <- seq(1, length(yraw[[counter]]))
            ysmooth[[counter]] <- loess(yraw[[counter]] ~ xraw[[counter]], span = (windowSize/seqlength), degree = 2)
            yline[[counter]] <- predict(ysmooth[[counter]])
        }
    }
    
    titles <- character()
    if (!is.null(seqTwo)) {
        if (!is.null(refTwo)) {
            titles[1:4] <- c(paste0("seq1 ref1 CPS: ", round(ymean[1], 4)), paste0("seq1 ref2 CPS: ", round(ymean[2], 4)), paste0("seq2 ref1 CPS: ", round(ymean[3], 
                4)), paste0("seq2 ref2 CPS: ", round(ymean[4], 4)))
        } else {
            titles[1:2] <- c(paste0("seq1 ref1 CPS: ", round(ymean[1], 4)), paste0("seq2 ref1 CPS: ", round(ymean[2], 4)))
        }
    } else {
        if (!is.null(refTwo)) {
            titles[1:2] <- c(paste0("seq1 ref1 CPS: ", round(ymean[1], 4)), paste0("seq1 ref2 CPS: ", round(ymean[2], 4)))
        } else {
            titles[1] <- c(paste0("seq1 ref1 CPS: ", round(ymean[1], 4)))
        }
    }
    
    ylimit <- c(floor(min(as.numeric(unlist(yline)))), ceiling(max(as.numeric(unlist(yline)))))
    xlimit <- c(0, length(yline[[1]]))
    
    par(mar = c(4, 4, 4, 1))
    plot(1, 1, type = "l", col = "white", xlim = xlimit, ylim = ylimit, yaxt = "n", xaxt = "n", ann = FALSE)
    
    axis(1, at = xlimit, labels = c(start, end))
    axis(2, las = 1)
    title(ylab = "CPS")
    title(xlab = "Sequence Length (nucleotides)", line = 0.5)
    
    for (i in 1:counter) {
        lines(yline[[i]], col = linecol[i], lty = linetypes[i], lwd = 1.5)
    }
    if (counter > 1) {
        legend("top", legend = titles, box.col = "black", xpd = TRUE, lwd = 2, lty = linetypes, col = linecol, inset = c(0, -0.05), ncol = 2, x.intersp = 1.2, bg = "white")
    } else {
        legend("top", legend = titles, box.col = "black", xpd = TRUE, lwd = 2, lty = linetypes, col = linecol, inset = c(0, -0.05), x.intersp = 1.2, bg = "white")
    }
}


rtCPSplot <- function(refs, seqORF, score, Second.seqORF, Second.score, windowSize, geneLength, blank = FALSE) {
    if (refs == 1) {
        x <- seq(1, length(seqORF), by = 1)
        y <- seqORF
        
        if (is.null(windowSize)) {
            windowSize <- geneLength * 0.075
        }
        
        smooth.line.ORF <- loess(y ~ x, span = (windowSize/geneLength), degree = 2)
        ylimit <- c(floor(min(as.numeric(y))) - 0.2, ceiling(max(as.numeric(y))) + 0.2)
        xlimit <- c(0, length(x))
        
        if (blank) {
            par(mar = c(4, 4, 1, 1))
            plot(1, 1, type = "l", xlim = xlimit, col = "white", ylim = ylimit, yaxt = "n", xaxt = "n", ann = FALSE)
            axis(1, at = c(1, length(x)), labels = c(1, length(x) * 3))
            axis(2, las = 1)
            title(ylab = "CPS")
            title(xlab = "Sequence Length (nucleotides)", line = 0.5)
        } else if (!blank) {
            polygon(c(-1e+06, -1e+06, 1e+06, 1e+06, -1e+06), c(100, -100, -100, 100, 100), col = "white", xpd = FALSE)
            box()
            lines(predict(smooth.line.ORF), type = "l")
            mtext(side = 3, line = -1, paste0("CPS: ", round(score, 5)))
        }
    } else if (refs == 2) {
        if (is.null(windowSize)) {
            windowSize <- geneLength * 0.075
        }
        x <- seq(1, length(seqORF), by = 1)
        y <- seqORF
        smooth.line.ORF <- loess(y ~ x, span = (windowSize/geneLength), degree = 2)
        Second.y <- Second.seqORF
        Second.smooth.line.ORF <- loess(Second.y ~ x, span = (windowSize/geneLength), degree = 2)
        ylimit <- c(floor(min(as.numeric(c(y, Second.y)))) - 0.2, ceiling(max(as.numeric(c(y, Second.y)))) + 0.2)
        if (blank) {
            plot(1, 1, type = "l", xlim = c(0, length(x)), lwd = 1.5, col = "white", ylim = ylimit, yaxt = "n", xaxt = "n", ann = FALSE)
            axis(1, at = c(1, length(x)), labels = c(1, length(x) * 3))
            axis(2, las = 1)
            title(ylab = "CPS")
            title(xlab = "Sequence Length (nucleotides)", line = 0.5)
            
        } else if (!blank) {
            polygon(c(-1e+06, -1e+06, 1e+06, 1e+06, -1e+06), c(100, -100, -100, 100, 100), col = "white", xpd = FALSE)
            box()
            lines(predict(smooth.line.ORF), lwd = 1.5, col = "red")
            lines(predict(Second.smooth.line.ORF), lwd = 1.5, col = "blue")
            mtext(side = 3, line = -1, paste0("CPS [reference 1]: ", round(score, 5)))
            mtext(side = 3, line = -2, paste0("CPS [reference 2]: ", round(Second.score, 5)))
            
        }
    }
}

newCPSplot <- function(Sequences, References, windowSize, Colors, Limits, Legend) {
    xraw <- list()
    yraw <- list()
    ymean <- numeric()
    ysmooth <- list()
    yline <- list()
    plot(1, 1, type = "l", col = "white", xlim = Limits[[1]], ylim = Limits[[2]], yaxt = "n", xaxt = "n", ann = FALSE)
    for (i in 1:length(Sequences)) {
        start <- 1
        end <- nchar(Sequences[[i]])
        seqlength <- end - start
        yraw[[i]] <- as.numeric(CPScalc(Sequences[[i]], References[[i]], start, end, draw = FALSE, silent = TRUE)[[2]])
        ymean[i] <- as.numeric(CPScalc(Sequences[[i]], References[[i]], start, end, draw = FALSE, silent = TRUE)[[1]])
        xraw[[i]] <- seq(1, length(yraw[[i]]))
        ysmooth[[i]] <- loess(yraw[[i]] ~ xraw[[i]], span = (windowSize/seqlength), degree = 2)
        yline[[i]] <- predict(ysmooth[[i]])
        lines(yline[[i]], col = Colors[i], lwd = ifelse(i == 1, 1.5, 0.8))
    }
    axis(1, at = Limits[[1]], labels = Limits[[1]] * 3)
    axis(2, at = Limits[[2]], las = 1)
    title(ylab = "CPS")
    title(xlab = "Sequence Length (nucleotides)", line = 0.5)
    legend("topleft", legend = Legend, box.col = "black", xpd = TRUE, lwd = 2, col = Colors, inset = c(0, -0.05), ncol = 2, x.intersp = 1.2, bg = "white")
} 
