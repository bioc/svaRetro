#' Creating a GRanges of RT insertion site breakpoints
#' @details 
#' This is an internal function which returns a GRanges of candidate
#' RT insertion site breakpoints with some initial filtering
#' @param hits.tx hits dataframe containing breakpoints overlapping with
#' exon edges.
#' @param exons TxDb object with transcript info.
#' @param gr SV calls in a GRanges object
#' @param same.tx A logical vector of whether the paired breakpoints belong to 
#' the same transcript
#' @param idx A dataframe of junctions with only one side overlapping with exons
#' @keywords internal
#' @return a GRanges of RT insertion site breakpoints
annotate_is_gr <- function(exons, hits, gr, same.tx, idx){
    #insSite.gr include bps with both AND one bnd mapped to exon boundry
    insSite.gr <- c(gr[hits[!same.tx,]$queryHits], 
                    partner(gr)[hits[!same.tx,]$queryHits], 
                    gr[idx$queryHits])
    #annotate with exon ids
    exon_id1 = exons[hits[!same.tx,]$subjectHits.x]$exon_id
    exon_id2 = exons[hits[!same.tx,]$subjectHits.y]$exon_id
    exon_id3 = exons[idx$subjectHits]$exon_id
    insSite.gr$exons <- c(exon_id1, exon_id2, exon_id3)
    #annotate with transcript names
    tx_name1 = exons[hits[!same.tx,]$subjectHits.x]$tx_name
    tx_name2 = exons[hits[!same.tx,]$subjectHits.y]$tx_name
    tx_name3 = exons[idx$subjectHits]$tx_name
    insSite.gr$txs <- c(tx_name1, tx_name2, tx_name3)
    #remove bps without tx names
    insSite.gr <- insSite.gr[!vapply(insSite.gr$txs, 
                                     rlang::is_empty, logical(1))]
    #return(insSite.gr)
    #combine matching exons and transcripts of the same breakend
    names <- unique(names(insSite.gr))
    #print(names)
    #print(insSite.gr[names(insSite.gr)==names[1]])
    #insSite.txs <- lapply(names, function(x) {Reduce(union, insSite.gr[names(insSite.gr)==x]$txs)})
    insSite.txs <- lapply(names, 
                          function(x) union_items(x, insSite.gr, "txs"))
    names(insSite.txs) <- names
    insSite.exons <- lapply(names, 
                            function(x) union_items(x, insSite.gr, "exons"))
    names(insSite.exons) <- names
    insSite.gr$txs <- insSite.txs[names(insSite.gr)]
    insSite.gr$exons <- insSite.exons[names(insSite.gr)]
    return(insSite.gr)
}

#' Filtering a GRanges of RT insertion site breakpoints
#' @details 
#' This is an internal function which returns a filtered GRanges of candidate
#' RT insertion site breakpoints
#' @param insSite.gr GRanges of candidate RT insertion site breakpoints
#' @param rt.gr GRanges of RT breakpoints
#' @param gr SV calls in a GRanges object
#' @param tx.rank A dataframe of RT transcripts with matching scores
#' @keywords internal
#' @return a GRanges of RT insertion site breakpoints
filterByRT_is_gr <- function(insSite.gr, rt.gr, gr, tx.rank){
    #remove duplicates and bps in rt.gr
    insSite.gr <- insSite.gr[!duplicated(names(insSite.gr))]
    insSite.gr <- insSite.gr[!names(insSite.gr) %in% names(rt.gr)]
    #add the missing partner bp to insSite.gr
    l <- !insSite.gr$partner %in% names(insSite.gr)
    insSite.gr <- c(insSite.gr, gr[insSite.gr[l]$partner])
    insSite.gr$rtFound <- mapply(stringr::str_detect, 
                                 insSite.gr$txs, 
                                 paste(tx.rank$tx_name, collapse = "|"))
    insSite.gr$rtFoundSum <- vapply(insSite.gr$rtFound, 
                                    function(x) {sum(x) > 0}, logical(1))
    return(insSite.gr)
}