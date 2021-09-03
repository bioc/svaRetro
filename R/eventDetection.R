#' Detecting retrotranscript insertion in nuclear genomes.
#'
#' @details
#' This function searches for retroposed transcripts by identifying breakpoints 
#' supporting intronic deletions and fusions between exons and remote loci.
#' Only BND notations are supported at the current stage.
#' @param gr A GRanges object
#' @param genes TxDb object of genes. hg19 and hg38 are supported in the 
#' current version.
#' @param maxgap The maxium distance allowed on the reference genome between 
#' the paired exon boundries.
#' @param minscore The minimum proportion of intronic deletions of a 
#' transcript should be identified.
#' @return A GRangesList object, named insSite and rt, reporting breakpoints 
#' supporting insert sites and 
#' retroposed transcripts respectively. 'exon' and 'txs' in the metadata 
#' columns report exon_id and transcript_name from the 'genes' object.
#' @examples
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' genes <- TxDb.Hsapiens.UCSC.hg19.knownGene
#' vcf.file <- system.file("extdata", "diploidSV.vcf",
#'                          package = "svaRetro")
#' vcf <- VariantAnnotation::readVcf(vcf.file, "hg19")
#' gr <- breakpointRanges(vcf, nominalPosition=TRUE)
#' rt <- rtDetect(gr, genes, maxgap=30, minscore=0.6)
#' @export
#' 
rtDetect <- function(gr, genes, maxgap=100, minscore=0.4){
    #message("rtDetect")
    #check args
    assertthat::assert_that(is(gr, "GRanges"), 
                            msg = "gr should be a GRanges object")
    assertthat::assert_that(!isEmpty(gr), 
                            msg = "gr can't be empty")
    assertthat::assert_that(is(genes, "TxDb"), 
                            msg = "genes should be a TxDb object")
    
    #prepare annotation exons
    GenomeInfoDb::seqlevelsStyle(genes) <- GenomeInfoDb::seqlevelsStyle(gr)[1]
    exons <- exons(genes, columns=c("exon_id", "tx_id", "tx_name","gene_id"))
    
    #find exon-SV overlaps:
    hits.start <- findOverlaps(gr, exons, maxgap = maxgap, type = "start", 
                               ignore.strand = TRUE)
    hits.end <- findOverlaps(partner(gr), exons, maxgap = maxgap, type = "end", 
                             ignore.strand = TRUE)
    
    # 1.return breakpoints overlapping with exons on both ends (>=2 exons)
    hits <- dplyr::inner_join(dplyr::as_tibble(hits.start), 
                              dplyr::as_tibble(hits.end), by="queryHits")
    same.tx <- vapply(Reduce(BiocGenerics::intersect, 
                             list(mcols(exons)[hits$subjectHits.x, 'tx_id'], 
                                  mcols(exons)[hits$subjectHits.y, 'tx_id'])), 
                      length, numeric(1))!=0
    hits.tx <- hits[same.tx,]
    
    # 2.return breakpoints of insertionSite-exon 
    hits.insSite <- hits[!same.tx,] %>%
        dplyr::bind_rows(dplyr::anti_join(dplyr::as_tibble(hits.start), 
                                          dplyr::as_tibble(hits.end), 
                                          by='queryHits')) %>%
        dplyr::bind_rows(dplyr::anti_join(dplyr::as_tibble(hits.end), 
                                          dplyr::as_tibble(hits.start), 
                                          by='queryHits'))
    
    if (nrow(hits.tx)==0 & nrow(hits.insSite)==0) {
        message("There is no retroposed gene detected.")
        return(GRanges())
    }else{
        # 3.filter exon-exon junctions by minscore(>=2 exons)
        rt.gr <- annotate_rt_gr(exons, hits.tx, gr)
        rt.gr <- filterByExon_rt_gr(rt.gr)
        tx.rank <- .scoreByTranscripts(genes, unlist(rt.gr$txs)) 
        #dataframe of valid retro transcripts
        tx.rank <- tx.rank[tx.rank$score >= minscore,]
        rt.gr <- filterByScore_rt_gr(rt.gr, tx.rank, genes, minscore)
 
        # 4.filter insertion site junctions, reduce duplications
        #junctions with only one side overlapping with exons:
        idx <- dplyr::bind_rows(dplyr::anti_join(dplyr::as_tibble(hits.start),
                                                 dplyr::as_tibble(hits.end),
                                                 by='queryHits'),
                                dplyr::anti_join(dplyr::as_tibble(hits.end),
                                                 dplyr::as_tibble(hits.start),
                                                 by='queryHits'))
        insSite.gr <- annotate_is_gr(exons, hits, gr, same.tx, idx)
        insSite.gr <- filterByRT_is_gr(insSite.gr, rt.gr, gr, tx.rank)

        # 5.create one GrangesList per gene
        #get all genes detected
        rt.gr$gene_symbol <- .txs2genesym(rt.gr$txs)
        insSite.gr$gene_symbol <- .txs2genesym(insSite.gr$txs)
        # unlisted gene_symbols are factors which are converted to numbers
        # if not converted to characters prior to unique()
        l_gene_symbol <- unique(c(as.character(unlist(rt.gr$gene_symbol)),
                                  as.character(unlist(insSite.gr$gene_symbol))))

        #RT GRangesList by gene
        rt.gr.idx <- lapply(l_gene_symbol, 
                            function(gs) vapply(rt.gr$gene_symbol, 
                                                function(x) gs %in% x, 
                                                logical(1)))
        rt.grlist <- stats::setNames(lapply(rt.gr.idx,
                                     function(i) list(rt=rt.gr[i])), 
                                     l_gene_symbol)

        #InsSite GRangesList by gene
        insSite.gr.idx <- lapply(l_gene_symbol, 
                                 function(gs) vapply(insSite.gr$gene_symbol, 
                                                     function(x) gs %in% x, 
                                                     logical(1)))
        #include partnered insSite bnds which don't have genesym labeling (NA)
        insSite.grlist <- stats::setNames(
            lapply(insSite.gr.idx,
                   function(i) list(insSite=c(insSite.gr[i], 
                                              partner(insSite.gr)[i]))),
            l_gene_symbol)

        #group inssite and rt as one GRangesList
        gr.list <- S4Vectors::pc(rt.grlist, insSite.grlist)
        gr.list <- lapply(gr.list, 
                          function(x) stats::setNames(x, 
                                                      c('junctions', 
                                                        'insSite')))

        #TODO: add L1/Alu annotation for insertion site filtering.
        #TODO: single exon tx awareness in filtering
        return(gr.list)
        
    }
}



