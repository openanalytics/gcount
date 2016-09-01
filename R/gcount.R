#' gcount: Obtain Raw Read Counts from RNASeq Data
#' 
#' \code{gcount} is an R-package that allows to obtain \emph{raw} read counts 
#' quickly and easily from RNASeq data to be used in downstream analyses.
#' 
#' The current functionality is as follows:
#' 
#' \itemize{
#' \item Provide a \code{bam/bed} file or object along with gtf/gff 
#' file or object as input.
#' \item Specify sequencing type: \code{paired} or \code{single} end data.
#' \item Specify library type: \code{unstranded}, \code{first-strand} or 
#' \code{second-strand} specific.
#' \item Filter reads based on number of mismatches prior to counting if 
#' necessary.
#' \item Consider or ignore reads that map to \emph{overlapping genes}.
#' \item Count reads that overlap feature of interest -- \code{genes}, 
#' \code{exons}, \code{introns} etc.
#' }
#'  
#' @aliases gcount gcount-package
#' @docType package
#' @name gcount
#' @importFrom utils packageVersion capture.output head
#' @import methods
#' @import gread
#' @import data.table
#' @import BiocGenerics
#' @importFrom GenomicRanges GRanges split mcols
#' @importMethodsFrom GenomicRanges disjoin reduce intersect findOverlaps 
#' @importMethodsFrom GenomicRanges countOverlaps seqnames start end strand
#' @importFrom GenomicAlignments findOverlaps summarizeOverlaps 
#' @importFrom GenomicAlignments readGAlignments readGAlignmentPairs cigar 
#' @importFrom GenomicAlignments qwidth width njunc rname
#' @importFrom S4Vectors queryHits subjectHits
#' @importMethodsFrom S4Vectors Rle elementMetadata
#' @import bitops
#' @seealso \code{\link{get_counts}}
NULL
