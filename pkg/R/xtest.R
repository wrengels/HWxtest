
#' Performs an exact test with full enumeration for Hardy-Weinberg proportions
#' 
#' Given a set of genotype counts, \code{xtest} examines all possible outcomes with the same set of allele counts. For each table, it computes four test statistics and compares them with the observed values. It returns the total probability of all tables with test statistics as \dQuote{extreme} or more so than the observed. It can also plot a histogram of one of the statitistics if \code{histobins} is greater than zero. More about these four test statistics and other information can be found in the vignette.
#' 
#' @param c A matrix containing the genotype counts. It should be a square matrix, but only the lower-left half is used.
#' @param statName can be \dQuote{LLR}, \dQuote{Prob}, \dQuote{U}, or \dQuote{Chisq} depending on which one is to be ploted. Note that P values for all four are computed regardless of which one is specified with this parameter.
#' @param histobins If 0 no histogram is plotted. If 1 or \code{TRUE} a histogram with 500 bins is plotted. If set to a number greater than 1, a histogram with \code{histobins} is plotted.
#' @param histobounds A vector containing the left and right boundaries for the histogram's x axis. If you leave this as the default, \code{c(0,0)}, then \code{xtest} will compute reasonable bounds to include most of the distribution.
#' @param showCurve whether to show a blue curve indicating the asymptotic (chi squared) distribution. This only works for \code{LLR} and \code{Chisq}
#' @param safeSecs After this many seconds the calculation will be aborted. This is a safety valve to prevent attempts to compute impossibly large sets of tables.
#' @param detail Determines how much detail is printed. If set to 0, nothing is printed (useful if you use \code{xtest} programmatically.).
#' 
#' @return \code{xtest} returns a list components
#' \item{$ Pvalues}{The four computed P values corresponding to the test statistics: \code{LLR}, \code{Prob}, \code{U} and \code{Chisq} in that order.}
#' \item{$ observed}{The four observed statistics in the same order as above}
#' \item{$ tableCount}{The number of tables examined during the calculation}
#' \item{$ genotypes}{The input matrix of genotype counts}
#' \item{$ alleles}{The allele counts \eqn{m} corresponding to the input genotype counts}

#' 
#' @examples
#' # Data from Louis and Dempster 1987 Table 2 and Guo and Thompson 1992 Figure 2:
#' c <- vec.to.matrix(c(0,3,1,5,18,1,3,7,5,2))
#' xtest(c)
#' # To see a histogram of the LLR statistic:
#' xtest(c, histobins=TRUE)
#' # For a histogram of the U statistic and other details of the result:
#' xtest(c, statName="U", histobins=TRUE, detail=4)
#' 

#' @useDynLib HWxtest
#' @export
xtest <- 
function(c, statName="LLR", histobins=0, histobounds=c(0,0), showCurve=T, safeSecs=100, detail=2) {
	statNames <- c("LLR", "Prob", "U", "Chisq");
	statID <- which(statNames==statName);
	m <- alleleCounts(c);
	if(histobins == 1){ histobins  <- 500};   #The default is 500 bins
	ostats <- c(observedLLR(c), observedProb(c), observedU(c), observedX2(c));
	if(histobounds[1]==histobounds[2] && histobins) histobounds <- defaultHistobounds(ostats, statID, m);
	x <- .C("xtest",
			m=as.integer(sort(m, decreasing=T)),
			nAlleles=as.integer(length(m)),
			observed=as.double(ostats),
			Pvalues=as.double(double(4)),
			statID=as.integer(statID-1),
			histobins=as.integer(histobins),
			histobounds=as.double(histobounds),
			histoData=as.double(double(histobins + 1)),
			safeSecs=as.integer(safeSecs),
			tableCount=as.double(0)
			,PACKAGE="HWxtest"
			);
	comments = character(4);
	if(ostats[3]>=0) {comments[[3]] <- " (test for homozygote excess)"}
	if(ostats[3]<0) {comments[[3]] <- " (test for heterozygote excess)"};
	if(x$tableCount < 0) stop("Calculation timed out. You can change the time limit by setting parameter 'safeSecs'");
	if(detail==1) cat("\nP value (", statName,") = ", formatC(x$Pvalues[statID]), comments[statID], sep="");
	if(detail>=2) for(i in 1:4){cat("\nP value (",statNames[i],")",strtrim("     ", 6-nchar(statNames[i])),"= ", formatC(x$Pvalues[i]), comments[i], sep="")};
	if(detail>=3) cat("\n\nExamined", x$tableCount, "tables\n");
	if(detail>=4) {cat("\nObserved Test Statistics:\n");
		for(i in 1:4){cat("\n",strtrim("     ", 6-nchar(statNames[i])), statNames[i],"  :  ", formatC(ostats[i]))}}
	if(histobins) plotHistogram(ostats, statID, m, histobins, histobounds, x$histoData, showCurve);
	return = list(Pvalues=x$Pvalues, observed=ostats, tableCount= x$tableCount, genotypes=c, alleles=m)

}