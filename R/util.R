#' Replaces the NA values in a with corresponding values in b
#' @param a,b objects to be tested or coerced.
#' @keywords internal
#' @return The altered object.
'%na%' <- function(a, b) {
	if (is.null(a) || length(a) == 0) return(b)
	if (is.null(b) || length(b) == 0) return(a)
	return(ifelse(is.na(a), b, a))
}

#' Uses b if a is NULL
#' @param a,b objects to be tested or coerced.
#' @keywords internal
#' @return An un-null object.
'%null%' <- function(a, b) {
	if (is.null(a)) return(b)
	return (a)
}

#' Merge elements of a column from a subset of GRanges elements
#' @param name Names of GRanges elements.
#' @param gr The GRanges set
#' @param item_col The GRanges colomn of interest
#' @keywords internal
#' @return A vector.
union_items <- function(name, gr, item_col){
    Reduce(union, mcols(gr[names(gr)==name])[[item_col]])
}


