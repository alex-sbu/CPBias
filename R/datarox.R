#' @name listCPB
#' @description
#' Table of species specific CPB reference tables and information about the sequences used to create them.
#' @details
#' Sequences were collected for those organisms which have whole genome sequence data available in the EMBL database (ENA Release 121). Only protein coding sequences were used, mitochondrion and chloroplast sequences were omitted.
#'
#' @note
#' Only the Homo.sapiens and Aedes.aegypti CPB reference table data sets are included in this package, the rest of the references listed can be downloaded at \url{https://github.com/zayets} as rda files and loaded directly in R.
#'
#' @title Codon pair bias reference tables
#' @format
#' \describe{
#' \item{Species}{ Species scientific name.}
#' \item{EMBL.ID}{ EMBL WGS project ID number.}
#' \item{Taxonomy}{ Species lineage.}
#' \item{Number.of.Sequences}{ Number of sequences used in CPB calculation.}
#' \item{CPB.reference}{ Name of R object containing the CPB reference table.}
#' }
#' @source \url{http://www.ebi.ac.uk/ena}
#' @examples
#' # write to CSV
#' data(listCPB)
#' write.csv(listCPB, 'listCPB.csv', row.names=FALSE)
#'
NULL

#' @name REseqs
#' @description
#' Restriction enzyme recognition sequences in the form of R regular expressions.
#' @title Restriction enzyme table
#' @format
#' \describe{
#' \item{Enzyme}{ Restriction enzyme.}
#' \item{Prototype}{ Prototype enzyme for the recognition sequence.}
#' \item{Organism}{ Organism in which enzyme was isolated.}
#' \item{Recognition.sequence}{ 5' to 3' restriction enzyme recognition sequence, with cleavage marked where applicable.}
#' \item{Regex.sequence}{ 5' to 3' recognition sequence represented as a regular expression.}
#' }
#' @source \url{http://rebase.neb.com}
NULL


#' @name standardTranslation
#' @description
#' Standard code translation table as a two column data frame containing all unique coding and non-coding codons in column one and a single letter identifier for the encoded amino acid in column two.
#' @title Standard translation table
#' @format
#' \describe{
#' \item{codon}{ Individual codons.}
#' \item{amino}{ One letter identifiers for amino acids or `STOP' for stop codons.}
#' }
#' @source \url{http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi}
NULL 
