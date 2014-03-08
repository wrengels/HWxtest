# xcount
# (c) William R. Engels, 2014
# 'Exact Tests For Hardy-Weinberg Proportions', 2009, Genetics, 183, pp1431-1441


#' Find Exact Number of Genotype Tables
#' 
#' Use \code{xcount} to determine the exact number of tables (i.e., genotype numbers) for a given set of allele counts. This method enumerates all tables, and is best when the total number is less than 10^10 or so. This function is mostly called by \code{hw.test} rather than directly by the user.
#' 
#' @param m vector containing the numbers of alleles of each type. Length must be at least 2 and all must be positive integers.
#' @param safety Stop execution if the approximate table number obtained from \code{acount()} is more than this cutoff.
#' @param safeSecs Time limit in seconds. Another safety feature to prevent getting stuck in a too-long computation
#' 
#' @return The exact number of tables
#' 
#' @examples
#' # Allele counts from human Rh locus. Guo and Thompson, 1992, Figure 1
#' #
#' alleles <- c(15, 14, 11, 12, 2, 2, 1, 3)
#' xcount(alleles)



# dyn.load("~/DropBox/HWxtest/pkg/src/HWxcount.so")


#' @useDynLib HWxtest
#' @export xcount
xcount <- 
function(x, ...) {
	UseMethod("xcount")
}


#' @method xcount integer
#' @S3method xcount integer
xcount.integer <- 
function(m, safety = 1e10, safeSecs = 10) {
	if(class(m)=="table") m <- unclass(m);
	if(class(m)=="matrix") m <- alleleCounts(m);
	if(length(m) < 2) stop("\nThere must be at least two alleles\n");
	if(any(m < 1)) stop("\nThere must be at least one copy of each allele\n");
	ac <- acount(m);
	if(ac > safety) stop("\nToo many to count. Try increasing safety paramater\n");
	xc <- -1;
	value <- .C("xcount",
		counts=as.integer(sort(m, decreasing=T)),
		nAlleles = as.integer(length(m)),
		tableCount=as.double(xc),
		safeSecs=as.integer(safeSecs)
		,PACKAGE="HWxtest"
		);
	n <- value$tableCount;
	if(n < 0) {
		warning("\nOperation timed out after ", -n, 
		"tables.\nThe true count is greater than that.\nNegative value for count reported.\nTime limit is defined by parameter safeSecs.")
	}
	value$tableCount;
}

#' @method xcount matrix
#' @S3method xcount matrix
xcount.matrix <- 
function(c, ...) {
	m <- alleleCounts(c);
	xcount.integer(m,...)
}

#' @method xcount numeric
#' @S3method xcount numeric
xcount.numeric <- 
function(c, ...) {
	m <- as.integer(c);
	xcount.integer(m,...)
}

#' @method xcount genotype
#' @S3method xcount genotype
xcount.genotype <- 
function(x, ...) {
	tab <- table(factor(allele(x, 1), levels = allele.names(x)), factor(allele(x, 2), levels = allele.names(x)));
	m <- alleleCounts(unclass(t(tab)));
	xcount.integer(m,...)
}


