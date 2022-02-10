#' Creating a GRanges of exon-exon junction breakpoints
#' @details 
#' This is an internal function which returns a GRanges of candidate
#' exon-exon junction breakpoints with some simple filtering
#' @param hits.tx hits dataframe containing breakpoints overlapping with
#' exon edges.
#' @param exons GRanges object with transcript info.
#' @param gr SV calls in a GRanges object
#' @keywords internal
#' @return a GRanges of candidate exon-exon junctions
annotate_rt_gr <- function(exons, hits.tx, gr){
    tx_name1 = exons[hits.tx$subjectHits.x]$tx_name
    tx_name2 = exons[hits.tx$subjectHits.y]$tx_name
    exon_id1 = exons[hits.tx$subjectHits.x]$exon_id
    exon_id2 = exons[hits.tx$subjectHits.y]$exon_id
    #annotate with exon ids
    rt.gr<- c(gr[hits.tx$queryHits], partner(gr)[hits.tx$queryHits])
    rt.gr$exon <- c(exon_id1, exon_id2)
    #annoatate with tx names
    txs <- mapply(intersect, tx_name1, tx_name2, SIMPLIFY = F)
    rt.gr$txs <- c(IRanges::CharacterList(txs), IRanges::CharacterList(txs))
    #remove bps without tx names
    rt.gr <- rt.gr[!vapply(rt.gr$txs, rlang::is_empty, logical(1))]
    #combine matching exons and transcripts of the same breakend
    names <- unique(names(rt.gr))
    #update tx names
    rt.txs <- lapply(names, function(x) union_items(x, rt.gr, "txs"))
    names(rt.txs) <- names
    
    rt.gr$txs <- rt.txs[names(rt.gr)]
    #update exon ids
    rt.exons <- lapply(names, function(x) union_items(x, rt.gr, "exon"))
    names(rt.exons) <- names
    rt.gr$exons <- rt.exons[names(rt.gr)]
    #remove duplicate breakend records
    #unique() and duplicated() for granges compare RANGES, not names
    rt.gr <- rt.gr[!duplicated(names(rt.gr))]
    return(rt.gr)              
}

#' Initial filtering a GRanges of exon-exon junction breakpoints
#' @details 
#' This is an internal function which returns a filtered GRanges of candidate
#' exon-exon junction breakpoints filtered by number of matching exons
#' @param rt.gr GRanges of RT breakpoints
#' @keywords internal
#' @return a GRanges of RT insertion site breakpoints
filterByExon_rt_gr <- function(rt.gr){
    #RT filter 1: breakpoint should have at least one set of matching exon
    rt.gr <- rt.gr[!mapply(identical, partner(rt.gr)$exons, rt.gr$exons) | 
                       (mapply(identical, partner(rt.gr)$exons, rt.gr$exons) & 
                            vapply(rt.gr$exons, length, numeric(1))>1)]
    return(rt.gr)
}

#' Filtering a GRanges of RT insertion site breakpoints
#' @details 
#' This is an internal function which returns a filtered GRanges of candidate
#' exon-exon junction breakpoints filtered by matching scores
#' @param insSite.gr GRanges of candidate RT insertion site breakpoints
#' @param tx.rank A dataframe of RT transcripts with matching scores
#' @param genes TxDb object of genes.
#' @param minscore The minimum proportion of intronic deletions of a 
#' transcript should be identified.
#' @keywords internal
#' @return a GRanges of RT insertion site breakpoints
filterByScore_rt_gr <- function(rt.gr, tx.rank, genes, minscore){
    #RT filter 2:minimal proportion of exon-exon detected for a transcript
    #remove rows and transcripts which are not in the tx.rank
    rt.gr <- rt.gr[stringr::str_detect(unstrsplit(rt.gr$txs), 
                                       paste(tx.rank$tx_name, collapse = "|"))]
    rt.gr$txs <- mapply('[', 
                        rt.gr$txs, 
                        mapply(stringr::str_detect, 
                               rt.gr$txs, 
                               paste(tx.rank$tx_name, collapse = "|"), SIMPLIFY = F), SIMPLIFY = F)
    return(rt.gr)
}



