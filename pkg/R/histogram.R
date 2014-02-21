
#' Functions to plot a histogram of test statistic
#' 
#' Running  \code{xtest()} or \code{mtest()} can create data for a frequency distribution plot of one of the four test statistics. Then these functions can be used to make the plot. The user will not normally call these functions. Instead, they will be called by \code{xtest} or \code{mtest} provided \code{histobins} is positive.
#' 
#' @param obstats Observed statistics for the 4 test measures, \code{LLR}, \code{Prob}, \code{U} and \code{Chisq}.
#' @param statID Value 1-4 indicating which statistic to use for the plot.
#' @m A vector of the allele counts

#' 
#' @return \code{defaultHistobounds} returns a vector containing the left and right boundaries for the x axis.
#' 
#' @return \code{plotHistogram} does not return a value, but plots the histogram.


#' @rdname histogramFunctions
#' @export
defaultHistobounds <- 
function(ostats, statID, m) {
	k <- length(m);
	b <- double(2);
	df <- k *(k+1)/2 - k - 1;
	if(statID==1 || statID==4) b[2] <- qchisq(.999, df);
	if(statID==2) {
		# find maximum probability
		n <- sum(m)/2;
		ae <- matrix(0,k,k);
		for(i in 1:k) for(j in 1:k) ae[i,j] <- m[i]*m[j]/(2*n);
		for(i in 1:k)ae[i,i] <- ae[i,i]/2;
		b[1] <- (-2)*log(observedProb(ae));
		b[2] <- b[[1]] + qchisq(.999, df);
	}
	if(statID==3) {
		b[[1]] <- -abs(ostats[3]) -5;
		b[[2]] <- abs(ostats[3]) + 5;
	}
	return(b);	
}



# plot(seq(0, 16.3, length.out = 500),c$histoData, type = "h", xlab="The Test Statistic", col=c(rep("gray20", 200), rep("lightcoral", 400)))
