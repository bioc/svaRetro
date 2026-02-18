#' Adding gene symbol annotations
#' @details
#' This is an internal function which takes a list of txs and converts
#' them to gene symbols. Supports both UCSC (kgID) and Ensembl transcript IDs.
#' For UCSC IDs, uses a static lookup table. For Ensembl IDs, uses TxDb
#' gene_id mapping and org.Hs.eg.db for Entrez-to-symbol conversion.
#' @param txs A list of transcript ids.
#' @param genes Optional TxDb object, needed for Ensembl transcript ID lookup.
#' @param unique.genesyms TRUE or FALSE. If TRUE, the converted gene symbols
#' will remove duplicates.
#' @keywords internal
#' @return A list of names in gene symbols
.txs2genesym <- function(txs, genes = NULL, unique.genesyms = TRUE) {
    assertthat::assert_that(is(txs, "list") | is(txs, "vector"),
        msg = "txs should be a list or vector"
    )
    # First try: static UCSC kgID lookup table
    gene_symbol <- utils::read.delim(
        system.file("extdata", "gene_symbol.txt",
            package = "svaRetro"
        ),
        header = TRUE, comment.char = "#"
    )
    gene_symbol <- dplyr::bind_rows(
        gene_symbol,
        data.frame(kgID = NA, geneSymbol = NA)
    )

    gene_syms <- lapply(
        txs,
        function(x) {
            gene_symbol$geneSymbol[gene_symbol$kgID %in% x]
        }
    )

    # Check if UCSC lookup found anything
    has_matches <- any(vapply(gene_syms, length, numeric(1)) > 0)

    # Fallback: dynamic lookup for Ensembl transcript IDs via TxDb + org.Hs.eg.db
    if (!has_matches && !is.null(genes)) {
        # Build tx_name -> gene_id mapping from TxDb
        tx2gene <- tryCatch(
            {
                AnnotationDbi::select(genes,
                    keys = unique(unlist(txs)),
                    keytype = "TXNAME",
                    columns = "GENEID"
                )
            },
            error = function(e) NULL
        )

        if (!is.null(tx2gene) && requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
            # Map Entrez Gene IDs -> gene symbols
            entrez_ids <- unique(stats::na.omit(tx2gene$GENEID))
            if (length(entrez_ids) > 0) {
                id2sym <- tryCatch(
                    {
                        AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
                            keys = entrez_ids,
                            keytype = "ENTREZID",
                            column = "SYMBOL",
                            multiVals = "first"
                        )
                    },
                    error = function(e) NULL
                )

                if (!is.null(id2sym)) {
                    # Build tx_name -> symbol mapping
                    tx2sym <- stats::setNames(
                        id2sym[tx2gene$GENEID],
                        tx2gene$TXNAME
                    )
                    gene_syms <- lapply(txs, function(x) {
                        syms <- unique(stats::na.omit(unname(tx2sym[x])))
                        if (length(syms) == 0) character(0) else syms
                    })
                }
            }
        }
    }

    if (unique.genesyms) {
        gene_syms <- lapply(gene_syms, unique)
    }
    return(gene_syms)
}
