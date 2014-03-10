# Functions for arranging data used in HW test
# (c) William R. Engels, 2014


#' Make a matrix symmetrical, using only the lower-left half
#' 
#' In a matrix of genotype counts, \code{a[i,j]} and \code{a[j,i]} both represent the same
#' heterozygote if \code{i!=j}. So only half the matrix need be used. If \code{i==j} then
#' it is a homozygote.
#' 
#' @param t a matrix of non-negative integers representing genotype counts. At least 2x2
#' 
#' @return a matrix of the same dimensions as \code{t} but symmetrical. Upper-right is replaced.
#' 
#' 
#' @examples
#' d <- matrix(1:25,5,5)
#' fillUpper(d)
#' 
#' @export
fillUpper <- 
function(t){
	if(!(is.matrix(t) && (dim(t)[1]==dim(t)[2]) && dim(t)[1] > 1)) stop("Must be square matrix at least 2x2")
	k <- dim(t)[1];
	for(j in 2:k) {t[1:(j-1), j] <- t[j,1:(j-1)]};
	t
}



#' Extract allele counts from matrix of genotypes
#' 
#' Create a vector of the number of each allele. The sum will be double that of genotpes.
#' 
#' @param g a matrix of non-negative integers representing genotype counts
#' 
#' @return a vector of allele counts
#' 
#' @examples
#' d <- matrix(1:25,5,5)
#' alleleCounts(d)
#' 
#' 
#' @export
alleleCounts <- 
function(g) {
	t <- fillUpper(g);
	k <- dim(t)[1];
	m <- integer(k);
	for(i in 1:k) {m[i] <- sum(t[i,]) + t[i,i]};
	m
}



#' Convert a by-rows vector to a genotype matrix
#' 
#' If there are \code{k} alleles, there are \code{k(k+1)/2} possible genotypes. Convert a row-
#' wise vector of genotype counts into a \code{k x k} symmetrical matrix of genotype counts.
#' Genotype counts should be in the order: \code{a11, a21, a22, a31, a32, ..., a}
#' 
#' @param cv vector containing \code{k(k+1)/2} genotype counts. All non-negative integers.
#' 
#' @return a \code{k x k} symmetrical matrix of genotype counts.
#' 
#' 
#' @examples
#' vec.to.matrix(c(0,3,1,5,18,1,3,7,5,2))
#' 


#' @export
vec.to.matrix <- 
function(cv){
	nGenotypes <- length(cv)
	nAlleles <- as.integer((sqrt(8*nGenotypes + 1) - 1)/2)
	if(nGenotypes != nAlleles*(nAlleles + 1)/2) stop("\nWrong number of genotype counts")
	t <- matrix(NA, nAlleles, nAlleles)
	for(i in 1:nAlleles){t[i, 1:i] <- cv[(i*(i-1)/2 + 1):(i*(i+1)/2)]}
	t	
}


#' Extract a by-rows vector from a genotype matrix
#' 
#' Given a square matrix of size \eqn{k x k}, generate a vector of genotype counts of
#' size \eqn{k(k-1)/2}
#' 
#' @param t square matrix of genotype counts. Only the lower-left half is used.
#' 
#' @return a vector of length \eqn{k(k-1)/2} containing the genotype counts
#' 
#' @examples
#' t <- vec.to.matrix(c(0,3,1,5,18,1,3,7,5,2))
#' v <- matrix.to.vec(t)
#' 


#' @export
matrix.to.vec <- 
function(t){
	if(!(is.matrix(t) && (dim(t)[1]==dim(t)[2]) && dim(t)[1] > 1)) stop("Must be square matrix at least 2x2")
	v <- c();
	k <- dim(t)[1]
	for(i in 1:k){v <- append(v,t[i,1:i])}
	names(v) <- NULL;
	v	
}
NULL


# #' Convert a list of genotypes into a genotype count matrix
# #' 
# #' Genotype lists as are used by packages `genetics` and `adegenet` are converted to an array of genotype counts.
# #' This function requires package `genetics`
# #' 
# #' @param g List of text objects indicating genotypes. Alleles are separated by \dQuote{/}
# #' 
# #' @return matrix of \eqn{k x k} genotype counts
# #' 
# #' @examples
# #' g <- c(rep("a/b",3),
# #'		"b/b",
# #'		rep("a/c", 5),
# #'		rep("b/c", 18),
# #'		"c/c",
# #'		rep("a/d",3),
# #'		rep("b/d", 7),
# #'		rep("c/d", 5),
# #'		rep("d/d", 2))
# #' genotypeList.to.matrix(g)

# #' @export
# genotypeList.to.matrix <- 
# function(g){
	# if(!require(genetics)) stop("\ngenetics package is required")
	# x <- genotype(g);
	# tab <- table(factor(allele(x, 1), levels = allele.names(x)), factor(allele(x, 2), levels = allele.names(x)));
	# t(tab)
# }

