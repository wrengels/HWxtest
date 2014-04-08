# genepop.to.genind

#' Imports a \code{.txt} file in \code{GenePop} format into an object of type \code{genind}

#' 
#' The main work is done by the function \code{adegenet::import2genind}. However, that function requires text files with an extension of \code{.gen}, whereas such files usually have extension \code{.txt}. The sole purpose of this function is to work around the \dQuote{.gen} requirement.
#' 
#' @param name the name of a file in \code{GenePop} format
#' 
#' @return an object of class \code{genind}
#' 
#' @export
genepop.to.genind <- function(name){
	require(adegenet)
	tempfile <- file(name)
	tmp <- readLines(tempfile)
	writeLines(tmp, "tempgenepop.gen")
	ind <- import2genind("tempgenepop.gen")
	unlink("tempgenepop.gen")
	ind
}