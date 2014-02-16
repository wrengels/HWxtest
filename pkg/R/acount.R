# acount
# Bill Engels
# 'Exact Tests For Hardy-Weinberg Proportions', 2009, Genetics, 183, pp1431-1441


#' Find Approximate Number of Genotype Tables
#' 
#' Use \code{acount} to obtain the approximate number of genotype tables for a given set of allele counts. This method uses a normal approximation and is much faster than enumerating the tables.
#' 
#' @param m vector containing the numbers of alleles of each type. Length must be at least 2. All items are positive integers.
#' 
#' @return The approximate number of tables.
#' 
#' @examples
#' # Allele counts from human Rh locus. Guo and Thompson, 1992, Figure 1
#' #
#' alleles <- c(15, 14, 11, 12, 2, 2, 1, 3)
#' acount(alleles)
#' # This approximation may be compared with the exact value of 250552020
#' #
#' ld <- c(6329, 319, 47, 2773, 75, 6702, 14, 2, 333)
#' acount(ld)
#' #
#' # This is an example where the number of tables is too large for a full enumeration.




#' @export
acount <- 
function(m) {
	if(length(m) < 2) stop("\nThere must be at least two alleles\n");
	if(any(m < 1)) stop("\nThere must be at least one copy of each allele\n");
	n <- sum(m)/2;
	k <- length(m);
	summ2 <- sum(m^2);
	b <- k * (k+1)/2 -1;
	va <- n * b * (n+b+1)/((b+1)*(b+1)*(b+2));
	vm <- (k+1) * va;
	q <- ((k-1)/(vm * k)) * (summ2 - (4 * n * n/k));
	lnPm <- log(sqrt(k))  +  ((k-1)/2) * log((k-1)/(2*pi*k*vm)) - q/2;
	lns <- lchoose(n + b, b);
	exp(lns + lnPm)
}