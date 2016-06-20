gene_exon <- function(reads, exons, genes, minoverlap, library, paired, 
                multiple_feature_overlaps, verbose) {

    if (verbose) cat("Getting read counts for feature type 'gene_exon'.\n")
    # to please R CMD CHECK
    N=g_strand=m_strand=hits=fs=ss=first_strand=second_strand=gene_id=NULL
    i.first_strand=i.second_strand=unstranded=bam=NM=i.unstranded=NULL
    olaps = find_overlaps(reads, exons, ignore_strand=TRUE, type="any", 
                    minoverlap=minoverlap)
    olaps[, `:=`(m_strand = reads$strand[queryHits], 
                 g_strand = exons$strand[subjectHits], 
                 hits = match(exons$gene_id[subjectHits], genes$gene_id))]
    olaps = olaps[!(duplicated(olaps, by=c("queryHits", "hits")))]
    if (!multiple_feature_overlaps) {
      if (verbose) cat("Filtering reads that overlap more than 1 feature.\n")
      olaps = olaps[, "N" := .N, by="queryHits"][N == 1L][, "N" := NULL][]
    }
    if (!paired) {
        out = olaps[, list(g_strand = g_strand[1], 
                        first_strand = sum(m_strand != g_strand[1]), 
                        second_strand = sum(m_strand == g_strand[1]), 
                        unstranded = .N), by = hits]
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
              by=hits]
        out[, "unstranded" := first_strand + second_strand]
        # don't use unstranded = .N it counts 'bad' reads as well. 
    }
    exons[, "length" := end - start + 1L]
    transcripts.length = exons[, list(length = sum(length)), by = gene_id]

    genes = copy(genes)
    genes[, c("first_strand", "second_strand", "unstranded") := 0L]
    out[, "gene_id" := genes$gene_id[hits]]
    genes[out, `:=`(first_strand=i.first_strand, 
                     second_strand=i.second_strand, 
                     unstranded=i.unstranded), 
                     on="gene_id"]
    if (library == "firststrand") {
        genes[, "reads" := first_strand]
    } else if (library == "secondstrand") {
        genes[, "reads" := second_strand]
    } else if (library == "unstranded") {
        genes[, "reads" := unstranded]
    }
    i.length = NULL
    genes[, "length" := NA
        ][transcripts.length, "length" := i.length, on="gene_id"]
    setcolorder(genes, c("seqname", "start", "end", "length", "strand", 
        "transcript_id", "gene_id", "overlaps",  "first_strand", 
        "second_strand", "unstranded", "reads"))
    genes[]
}

# gene_all <- function(reads, genes, exons, minoverlap, library, 
#                   paired, verbose) {
# 
#     if (verbose) cat("Getting read counts for feature type 'gene_all'.\n")
#     olaps = find_overlaps(reads, genes, ignore_strand=TRUE, 
#                   type = "any", minoverlap = minoverlap)
#     olaps[, `:=`(m_strand = as.vector(strand(reads)[queryHits]), 
#                  g_strand = genes$strand[subjectHits])]
#     if (!paired) {
#         out = olaps[, list(g_strand = g_strand[1], 
#                         first_strand = sum(m_strand != g_strand[1]), 
#                         second_strand = sum(m_strand == g_strand[1]), 
#                         unstranded = .N), 
#               by = subjectHits]
#     } else {
#         olaps[, reads_flag := elementMetadata(reads)$flag[queryHits]]
#         rs.idx  = bitAnd(olaps$reads_flag, 48) == 48
#         fip.idx = bitAnd(olaps$reads_flag, 64) == 64 
#         sip.idx = bitAnd(olaps$reads_flag, 128) == 128
#         idx1 = olaps$m_strand != olaps$g_strand & fip.idx 
#         idx2 = olaps$m_strand == olaps$g_strand & sip.idx
#         idx = (idx1 | idx2) & !rs.idx
#         olaps[, fs := idx]
#         idx1 = olaps$m_strand == olaps$g_strand & fip.idx 
#         idx2 = olaps$m_strand != olaps$g_strand & sip.idx
#         idx = (idx1 | idx2) & !rs.idx
#         olaps[, ss := idx]
#         out = olaps[, list(g_strand = g_strand[1], 
#                         first_strand = sum(fs), 
#                         second_strand = sum(ss), 
#                         unstranded = .N),
#               by=subjectHits]
#     }
#     exons[, length := end - start + 1]
#     transcripts.length = exons[, list(length = sum(length)), by = gene_id]
#     
#     genes[out[, subjectHits], first_strand := out[, first_strand]
#           ][is.na(first_strand), first_strand := 0L]
#     genes[out[, subjectHits], second_strand := out[, second_strand]
#           ][is.na(second_strand), second_strand := 0L]
#     genes[out[, subjectHits], unstranded := out[, unstranded]
#           ][is.na(unstranded), unstranded := 0L]
#     if (library == "firststrand") {
#         genes[, reads := first_strand]
#     } else if (library == "secondstrand") {
#         genes[, reads := second_strand]
#     } else if (library == "unstranded") {
#         genes[, reads := unstranded]
#     }
#     # genes[, length := end - start + 1]
#     setkeyv(genes, "gene_id")
#     setkeyv(transcripts.length, "gene_id")
#     genes = transcripts.length[genes][is.na(length), length := 0L]   
#     setcolorder(genes, c("seqname", "start", "end", "length", 
#               "strand", "gene_id", "gene.type", "first_strand", 
#               "second_strand", "unstranded", "reads"))
#     genes
# }
