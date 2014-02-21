
#' Performs an exact test with full enumeration for Hardy-Weinberg proportions
#' 
#' Given a set of genotype counts, \code{xtest} examines all possible outcomes with the same set of allele counts. For each table, it computes four test statistics and compares them with the observed values. It returns the total probability of all tables with test statistics as \dQuote{extreme} or more so than the observed. It can also plot a histogram of one of the statitistics if \code{histobins} is greater than zero.
#' 
#' @param c A matrix containing the genotype counts. It should be a square matrix, but only the lower-left half is used.
#' @param statname can be \dQuote{LLR}, \dQuote{Prob}, \dQuote{U}, or \dQuote{Chisq} depending on which one is to be ploted. Note that P values for all four are computed regardless of which one is specified with this parameter.
#' @param histobins If 0 no histogram is plotted. If 1 or \code{TRUE} a histogram with 500 bins is plotted. If set to a number greater than 1, a histogram with \code{histobins} is plotted.
#' @param histobounds A vector containing the left and right boundaries for the histogram's x axis. If you leave this as the default, \code{c(0,0)}, then \code{xtest} will compute reasonable bounds to include most of the distribution.
#' @param safeSecs After this many seconds the calculation will be aborted. This is a safety valve to prevent attempts to compute impossibly large sets of tables.
#' @param detail Determines how much detail is printed. If set to 0, nothing is printed (useful if you use \code{xtest} programmatically.).
#' 


#' @export
xtest <- 
function(c, statName="LLR", histobins=0, histobounds=c(0,0), safeSecs=100, detail=2) {
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
			tableCount=as.double(0));
	if(x$tableCount < 0) stop("Calculation timed out. You can change the time limit by setting parameter 'safeSecs'");
	if(detail==1) cat("\nP value (", statName,") = ", formatC(x$Pvalues[statID]), sep="");
	if(detail>=2) for(i in 1:4){cat("\nP value (",statNames[i],")",strtrim("     ", 6-nchar(statNames[i])),"= ", formatC(x$Pvalues[i]), sep="")};
	if(detail>=3) cat("\n\nExamined", x$tableCount, "tables\n");
	if(detail>=4) {cat("\nObserved Test Statistics:\n");
		for(i in 1:4){cat("\n",strtrim("     ", 6-nchar(statNames[i])), statNames[i],"  :  ", formatC(ostats[i]))}}
	x
}