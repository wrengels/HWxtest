# hw.test  Generic function to test for HW by either full enumeration or by Monte Carlo
# (c) William R. Engels, 2014
# Accepts matrix, vector, table, genotype, etc. 


#' Test for HW by either full enumeration or Monte Carlo.
#' 
#' The \code{hw.test()} function is the main function of the \code{HWxtest} package. This function produces a valid test for Hardy-Weinberg frequencies for virtually any set of genotype counts. It will use either a full-enumeration method in which all possible tables with the same allele numbers are examined, or a Monte Carlo test where a large number of random tables is examined. To decide which to use, it calls \code{xcountCutoff} to determine whether the number of tables to examine is greater than \code{cutoff}. If it is, then \code{mtest} is used. Otherwise \code{xtest} is used. The result is a robust test which will always provide a meaningful and accurate P value. Each table examined is compared with the observed counts according to each of four measures of fit:  \dQuote{LLR}, \dQuote{Prob}, \dQuote{U}, or \dQuote{Chisq} corresponding to the log-likelihood ratio, the null-hypothesis probability, the U-score or the Pearson X^2 value. It can also plot a histogram showing the distribution of any of these statistics.
#' 
#' 
#' @aliases hw.test.matrix hw.test.integer
#' 
#' @param x The genotype counts. You must provide the number of each genotype. So if there are \eqn{k} alleles, you need to include the number of each of the \eqn{k(k+1)/2} genotypes. The format of \code{x} is somewhat flexible: It can be a square matrix, but only the lower-left half is used. It can be a vector of the observations in the order \eqn{a_11, a_21, a_22, a_31, ..., a_kk}. For compatability with the packages \code{genetics} and \code{adegenet}, it can also be an object of class \code{genotype}.
#' @param method Can be \dQuote{auto}, \dQuote{exact} or \dQuote{monte} to indicate the method to use. If \dQuote{auto}, the \code{hw.test} will first check to see whether the total number of tables exceeds a cutoff specified by the parameter \code{cutoff}.
#' @param cutoff If \code{method} is set to \dQuote{auto}, then \code{cutoff} is used to decide whether to perform the test via the full enumeration or Monte Carlo method. If the number of tables is less than \code{cutoff}, then a full enumeration is performed. Otherwise the method will be Monte Carlo with \code{B} random trials.
#' @param B The number of trials to perform if Monte Carlow method is used
#' @param statName can be \dQuote{LLR}, \dQuote{Prob}, \dQuote{U}, or \dQuote{Chisq} depending on which one is to be ploted. Note that P values for all four are computed regardless of which one is specified with this parameter.
#' @param histobins If 0, no histogram is plotted. If 1 or \code{TRUE} a histogram with 500 bins is plotted. If \code{histobins} is set to a number greater than 1, a histogram with \code{histobins} bins is plotted.
#' @param histobounds A vector containing the left and right boundaries for the histogram's x axis. If you leave this as the default, \code{c(0,0)}, then \code{hw.test} will compute reasonable bounds to include most of the distribution.
#' @param showCurve whether to show a blue curve indicating the asymptotic (chi squared) distribution. This only works for \code{LLR} and \code{Chisq}
#' @param safeSecs After this many seconds the calculation will be aborted. This is a safety valve to prevent attempts to compute impossibly large sets of tables.
#' @param detail Determines how much detail is printed. If it is set to 0, nothing is printed (useful if you use \code{hw.test} programmatically.)
#' 

#' 
#' @return Returns a list of class \code{hwtest} which includes the following items:
#' \item{$ Pvalues}{The four computed P values corresponding to the test statistics: \code{LLR}, \code{Prob}, \code{U} and \code{Chisq} in that order.}
#' \item{$ p.value}{The P value corresponding to the specified test}
#' \item{$ observed}{The four observed statistics in the same order as above}
#' \item{$ ntrials}{The number of tables examined during the calculation if done by Monte Carlo}
#' \item{$ tableCount}{The total number of tables if done by full enumeration}
#' \item{$ genotypes}{The input matrix of genotype counts}
#' \item{$ alleles}{The allele counts \eqn{m} corresponding to the input genotype counts}
#' \item{$ statName}{Which statistic to use for the histogram and in the \code{p.value} item}
#' \item{$ method}{Which method was used, \dQuote{exact} or \dQuote{monte}}
#' \item{$ detail}{An integer indicating how much detail to print. Use 0 for no printing}

#' 
#' @examples
#' # Data from Louis and Dempster 1987 Table 2 and Guo and Thompson 1992 Figure 2:
#' c <- c(0,3,1,5,18,1,3,7,5,2)
#' hw.test(c)
#' # To see a histogram of the LLR statistic:
#' hw.test(c, histobins=TRUE)
#' # For a histogram of the U statistic and other details of the result:
#' hw.test(c, statName="U", histobins=TRUE, detail=3)
#' 

#' @export hw.test

hw.test <- 
function(x, ...) {
	UseMethod("hw.test")
}

#' @method hw.test matrix
#' @S3method hw.test matrix
hw.test.matrix <- 
function(c, method="auto", cutoff=1e7, B=100000, statName="LLR", histobins=0, histobounds=c(0,0), showCurve=T, safeSecs=100, detail=2, ...) {
	statNames <- c("LLR", "Prob", "U", "Chisq");
	statID <- which(statNames==statName);
	if(method=="auto") method <- if(xcountCutoff(c, cutoff))method <- "monte" else method <- "exact"
	if(method=="monte") 
		value <- mtest(c, ntrials=B, statName, histobins, histobounds, showCurve, safeSecs, detail)
	else
		value <- xtest(c, statName, histobins, histobounds, showCurve, safeSecs, detail);
	value$method=method;
	value$p.value=value$Pvalues[statID];
	value$statName=statName;
	value$detail=detail;
	class(value) <- "hwtest"
	return(value)
}

#' @method hw.test integer
#' @S3method hw.test integer
hw.test.integer <- 
function(c, ...) {
	c <- vec.to.matrix(c);
	hw.test(c, ...)
}
#' @method hw.test numeric
#' @S3method hw.test numeric
hw.test.numeric <- 
function(c, ...) {
	c <- vec.to.matrix(c);
	hw.test(c, ...)
}

#' @method hw.test table
#' @S3method hw.test table
hw.test.table <- 
function(c, ...) {
	hw.test(unclass(c), ...)
}

#' @method hw.test genotype
#' @S3method hw.test genotype
hw.test.genotype <- 
function(x, ...) {
	tab <- table(factor(allele(x, 1), levels = allele.names(x)), factor(allele(x, 2), levels = allele.names(x)));
	hw.test(unclass(t(tab)), ...)
}
#' @method hw.test list
#' @S3method hw.test list
hw.test.list <- 
function(s, ...){
	lapply(s, hw.test, ...)
}

#' @method hw.test logical
#' @S3method hw.test logical
hw.test.logical <- 
function(s, ...) {
	a <- list(p.value=NA);
	class(a) <- "hwtest"
	return(a)
}