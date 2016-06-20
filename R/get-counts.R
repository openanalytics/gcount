#' Obtain read counts from \code{bam} object.
#'
#' @description \code{get_counts} allows to obtain read counts quickly and 
#' easily from RNASeq data to be used in downstream analyses.
#' 
#' @param reads A \code{bam/bed} object (see \code{gread::read_format}) or 
#' complete path to a \code{bam/bed} file.
#' @param annotation A \code{gtf/gff} object (see \code{gread::read_format}) 
#' or complete path to \code{gtf/gff} file.
#' @param transcript_id Column name in \code{x} corresponding to transcript 
#' id. Default value is \code{"transcript_id"}.
#' @param gene_id Column name in \code{x} corresponding to gene id. Default 
#' value is \code{"gene_id"}.
#' @param mismatches Default \code{-1} is to ignore mismatches. If \code{> 0}, 
#' reads with mismatches \code{<= mismatches} are alone retained.
#' @param minoverlap Argument that's passed to 
#' \code{GenomicRanges::findOverlaps}.
#' Default is \code{1L}, i.e., select all overlapping reads.  
#' @param feature How to count reads?
#' 
#' \code{"gene_exon"} counts only exonic reads within genes.
#' 
#' \code{"gene"} counts any/all reads overlapping the gene.
#' 
#' \code{"exon"} returns read counts for each exon separately.
#' 
#' \code{"intron"} returns read counts for each intron separately.
#' 
#' See \code{type} argument for more advanced operations on extracting 
#' feature coordinates.
#' @param type Same as \code{gread::extract}. See \code{?extract} in 
#' \code{gread}.
#' @param library Either \code{"unstranded"} (default), \code{"firststrand"} 
#' or \code{"secondstrand"}.
#' @param paired Default is \code{FALSE}. If the library is \code{paired end}, 
#' set it to \code{TRUE}
#' @param multiple_feature_overlaps logical. Should reads that overlap 
#' multiple features be counted. Default is \code{FALSE}, i.e., to discard. If 
#' \code{TRUE}, reads across overlapping features will be counted.
#' @param verbose logical. Default is \code{FALSE}. If \code{TRUE}, provides 
#' helpful messages to the console. 
#' @return A data.table with calculated raw counts of overlapping reads 
#' for each feature.
#' @aliases get_counts
#' @export
#' @examples
#' path = system.file("tests", package="gcount")
#' gtf_file = file.path(path, "sample.gtf")
#' bam_file = file.path(path, "sample.bam")
#' bam_counts = get_counts(bam_file, gtf_file, feature="gene_exon", 
#'              type="union", paired=FALSE, library="unstranded", 
#'              verbose=TRUE)
get_counts <- function(reads, annotation, transcript_id="transcript_id", 
                gene_id="gene_id", mismatches=-1L, minoverlap=1L, 
                feature=c("gene_exon", "gene", "exon", "intron"), 
                type=c("default", "union", "disjoin", "intersect", 
                        "longest", "shortest", "overlap"), 
                library=c("unstranded", "firststrand", "secondstrand"), 
                paired=FALSE, multiple_feature_overlaps=FALSE, 
                verbose=FALSE) {

    paired = as.logical(paired[1L])
    multiple_feature_overlaps = as.logical(multiple_feature_overlaps[1L])
    verbose = as.logical(verbose[1L])
    minoverlap = as.integer(minoverlap[1L])
    mismatches = as.integer(mismatches[1L])
    feature = match.arg(feature)
    type = match.arg(type)
    library = match.arg(library)

    # to please R CMD CHECK
    bam=NM=NULL
    if (!paired %in% c(TRUE, FALSE))
        stop("'paired' must be logical TRUE/FALSE.")
    if (!multiple_feature_overlaps %in% c(TRUE, FALSE))
        stop("'multiple_feature_overlaps' must be logical TRUE/FALSE.")
    if (!verbose %in% c(TRUE, FALSE))
        stop("'verbose' must be logical TRUE/FALSE.")
    if (is.na(minoverlap) || minoverlap < 1L) {
        warning("invalid 'minoverlap' value = ", minoverlap, 
            "; set back to default 1.")
        minoverlap = 1L
    }
    if (!inherits(annotation, "gtf") && !inherits(annotation, "gff")) {
        if (is.character(annotation)) {
            if (verbose) 
                cat("Argument 'annotation' is of type character. Assuming ", 
                    "it is a path to a gtf/gff file and attempting to load ", 
                    "using gread::read_format.\n", sep="")
            annotation = gread::read_format(annotation)
        } else {
            stop("Argument 'annotation' should be either a path to gtf/gff ", 
                "file or an object of class 'gtf' or 'gff' object; see ", 
                "gread::read_format.")
        }
    }
    if (!inherits(reads, "bed") && !inherits(reads, "bam")) {
        if (is.character(reads)) {
            if (verbose)
                cat("'reads' argument is of type character. Assuming ", 
                    "it is a path to bed/bam file and attempting to load ", 
                    "using gread::read_format.\n", sep="")
            reads = gread::read_format(reads)
        } else {
            stop("Argument 'reads' should be either a path to a bed/bam file", 
                " or an object of class 'bam' or 'bed; see ", 
                "gread::read_format.")
        }
    }
    if (is.na(mismatches) || mismatches < 0L) {
        if (verbose) cat("Reads are not be filtered based on mismatches.\n")
    } else {
        if (inherits(reads, "bam")) {
            if (verbose) {
                cat("Retaining only those reads with <= ", 
                        mismatches, " mismatches.\n", sep="")
            }
            stopifnot("NM" %chin% names(bam))
            reads = reads[NM <= mismatches]
        } else {
            cat("No reads are filtered from bed file. ", 
                "It doesn't contain mismatch info.\n", sep="")
        }
    }
    this = gread::extract(annotation, feature=feature, type=type, 
                    ignore_strand=FALSE, transcript_id=transcript_id, 
                    gene_id=gene_id)
    switch(feature,
      gene_exon = {
        # this would compute lengths incorrectly, but gene_exon takes 
        # care of this by computng length from 'this' (=exon) and 
        # replacing the length column.
        genes = gread::extract(annotation, feature=feature, type="default", 
                    ignore_strand=FALSE, transcript_id=transcript_id, 
                    gene_id=gene_id)
          rc = gene_exon(reads, this, genes, minoverlap, library, 
                        paired, multiple_feature_overlaps, verbose)
      },
      gene = {
        stop("Not yet implemented.")
      }, 
      exon = {
        rc = exon_only(reads, this, minoverlap, library, paired, verbose)
      },
      intron = {
        rc = intron_only(reads, this, minoverlap, library, paired, verbose)
      }
    )
    set(rc, j=c("first_strand", "second_strand", "unstranded"), value=NULL)
    # # just provide the chromosomes that were input'd
    # chrs = as.vector(attributes(seqnames(bam))$values)
    # if (!identical(sort(rc$seqname), sort(chrs))) {
    #     rc = rc[seqname %in% chrs]
    # }
    rc[]
}
