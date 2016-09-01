#' @title Find overlapping indices of two gtf/gff/bed/bam objects
#'
#' @description For internal use only. Function for finding overlaps between 
#' two objects of class \code{gtf/gff/bed/bam} using 
#' \code{GenomicRanges::findOverlaps}.
#' 
#' @param x,y An object of class \code{gtf}, \code{gff}, \code{bed} or 
#' \code{bam}.
#' @param ignore_redundant Should redundant overlaps be ignored?
#' @param ignore_strand Logical argument to pass to \code{GRanges} function. 
#' Indicates whether \code{strand} should be ignored when constructing 
#' \code{GRanges} object or not. Default is \code{FALSE}.
#' @param ... Additional arguments passed to 
#' \code{GenomicRanges::findOverlaps}.
#' @return A \code{data.table} containing overlapping indices.
find_overlaps <- function(x, y, ignore_redundant=FALSE, 
                    ignore_strand=FALSE, ...) {
    # this function is copied from gread
    stopifnot(inherits(x, "GRanges"), inherits(y, "GRanges"), 
                ignore_redundant %in% c(FALSE, TRUE), 
                ignore_strand %in% c(FALSE, TRUE))
    olaps = GenomicRanges::findOverlaps(x, y, 
                ignore.strand=ignore_strand,  ...)
    olaps = setDT(list(queryHits = queryHits(olaps), 
                    subjectHits = subjectHits(olaps)))
    # findOverlaps for GRanges objects doesn't seem to have ignoreRedundant
    # argument. so mimicing that functionality below.
    if (ignore_redundant) {
        olaps = olaps[, `:=`(queryHits = pmin(queryHits, subjectHits), 
                            subjectHits = pmax(queryHits, subjectHits))]
        olaps = unique(olaps, by=names(olaps))
    }
    olaps[]
}

is.gtf <- function(x) inherits(x, 'gtf')
is.gff <- function(x) inherits(x, 'gff')
is.bed <- function(x) inherits(x, 'bed')
is.bam <- function(x) inherits(x, 'bam')
