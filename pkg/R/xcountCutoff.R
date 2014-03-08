# xcountCutoff
# (c) William R. Engels, 2014
# 'Exact Tests For Hardy-Weinberg Proportions', 2009, Genetics, 183, pp1431-1441


#' Determine immediately whether number of tables is over a limit
#' 
#' Calling \code{scountCutoff} gives you a quick answer to whether the number of tables is over a given cutoff. It is useful in deciding whether to analyze a data set with \code{xtest} or \code{mtest}. This function is used by \code{hw.test} and not normally called directly by the user.
#' 
#' 
#' @param m vector containing the numbers of alleles of each type. Length must be at least 2 and all must be positive integers.
#' @param cutoff Is the number of tables above or below this value?
#' 
#' @method generic xcountCutoff
#' 
#' @return TRUE or FALSE depending on whether the table count is above or below \code{cutoff}
#' 
#' @examples
#' #
#' alleles <- c(15, 14, 11, 12, 2, 2, 1, 3)
#' if(xcountCutoff(alleles)) cat("There are too many tables")


#' @useDynLib HWxtest

#' @export xcountCutoff
xcountCutoff <- 
function(x, ...){
	UseMethod("xcountCutoff")
}

#' @method xcountCutoff integer
#' @S3method xcountCutoff integer
xcountCutoff.integer <- 
function(m, cutoff=1e7) {
	if(length(m) < 2) stop("\nThere must be at least two alleles\n");
	if(any(m < 1)) stop("\nThere must be at least one copy of each allele\n");
	if(acount(m) > 1e10) return(TRUE);
		value <- .C("xcount",
		counts=as.integer(sort(m, decreasing=T)),
		nAlleles = as.integer(length(m)),
		tableCount=as.double(cutoff),
		safeSecs=as.integer(5)
		,PACKAGE="HWxtest"
		);
		n <- value$tableCount;
		if(n < 0) return(TRUE);
		return(FALSE)
}

#' @method xcountCutoff matrix
#' @S3method xcountCutoff matrix
xcountCutoff.matrix <- 
function(c, ...) {
	m <- alleleCounts(c);
	xcountCutoff.integer(m,...)
}

#' @method xcountCutoff numeric
#' @S3method xcountCutoff numeric
xcountCutoff.numeric <-
function(c, ...) {
	m <- as.integer(c);
	xcountCutoff.integer(c,...)
}