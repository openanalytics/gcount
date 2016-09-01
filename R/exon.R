exon_only <- function(reads, exons, minoverlap, library, paired, verbose) {

    if (verbose) message("Getting read counts for feature type 'exon'")
    # to please R CMD CHECK
    g_strand=m_strand=fs=ss=first_strand=second_strand=unstranded=NULL
    seqname=NULL
    olaps = find_overlaps(reads, exons, ignore_strand=TRUE, 
                type="any", minoverlap=minoverlap)
    olaps[, `:=`(m_strand = as.character(strand(reads))[queryHits], 
                 g_strand = as.character(strand(exons))[subjectHits])]
    if (!paired) {
      out = olaps[, list(g_strand = g_strand[1], 
                     first_strand = sum(m_strand != g_strand[1]), 
                     second_strand = sum(m_strand == g_strand[1]), 
                     unstranded = .N), 
           by = subjectHits]
    } else {
        olaps[, "reads_flag" := reads$flag[queryHits]]
        # read and mate in reverse strand
        rs.idx  = bitAnd(olaps$reads_flag, 48) == 48
        # first in pair
        fip.idx = bitAnd(olaps$reads_flag, 64) == 64
        # second in pair
        sip.idx = bitAnd(olaps$reads_flag, 128) == 128
        # first in pair and reverse strand
        idx1 = olaps$m_strand != olaps$g_strand & fip.idx
        # second in pair and forward strand
        idx2 = olaps$m_strand == olaps$g_strand & sip.idx
        # either or, but not both in reverse strand
        idx = (idx1 | idx2) & !rs.idx
        # first strand
        olaps[, "fs" := idx]
        # first in pair and first strand
        idx1 = olaps$m_strand == olaps$g_strand & fip.idx
        # second in pair and reverse strand
        idx2 = olaps$m_strand != olaps$g_strand & sip.idx
        # either or, but not both in reverse strand
        idx = (idx1 | idx2) & !rs.idx
        # second strand
        olaps[, "ss" := idx]
        out = olaps[, list(g_strand = g_strand[1],
                        first_strand = sum(fs),
                        second_strand = sum(ss)),
              by=subjectHits]
        # don't use unstranded = .N it counts 'bad' reads as well. 
        out[, "unstranded" := first_strand + second_strand]
    }
    exons = data.table::setDT(as(exons, "data.frame"))
    del_cols = which(names(exons) %in% c("overlaps", "width"))
    if (length(del_cols)) set(exons, j=del_cols, value=NULL)
    setnames(exons, "seqnames", "seqname")
    exons[, "length" := end-start+1L
        ][, "seqname" := as.character(seqname)
        ][, "strand" := as.character(strand)]
    exons[, c("first_strand", "second_strand", "unstranded") := 0L]
    exons[out$subjectHits, c("first_strand", "second_strand", "unstranded") := 
                            list(first_strand, second_strand, unstranded)]
    if (library == "firststrand") {
        exons[, "reads" := first_strand]
    } else if (library == "secondstrand") {
        exons[, "reads" := second_strand]
    } else if (library == "unstranded") {
        exons[, "reads" := unstranded]
    }
    setcolorder(exons, c("seqname", "start", "end", "length", "strand", 
        "transcript_id", "gene_id",  "first_strand", "second_strand", 
        "unstranded", "reads"))
    exons[]
}
