# Plot histogram of HW test statistics
# (c) William R. Engels, 2014


#' Functions to plot a histogram of test statistic
#' 
#' Running the function \code{\link{hwx.test}}, \code{\link{xtest}} or \code{\link{mtest}} can create data for a frequency distribution plot of one of the four test statistics. Then these functions can be used to make the plot. The user will not normally call these functions. Instead, they will be called by \code{\link{hwx.test}} provided \code{histobins} is positive. 
#' 
#' @param ostats Observed statistics for the 4 test measures, \code{LLR}, \code{Prob}, \code{U} and \code{Chisq}.
#' @param statID Value 1-4 indicating which statistic to use for the plot.
#' @param m A vector of the allele counts
#' @param histobins the number of bins for the histogram
#' @param histobounds The left and right boundary of the histogram x-axis
#' @param histoData A vector of probabilities (or counts) returned from the \code{xtest} or \code{mtest} function
#' @param showCurve Whether to draw a blue curve indicating the asymptotic distribution (if known)
#' @param color1 The color for outcomes fitting the null distribution better than the observed
#' @param color2 The color for outcomes deviating from the null at leasst as much as observed. Area of \code{color2} is the P value.
#' @param ntrials If greater than 0, this is the number of trials used in a Monte Carlo test.

#' 
#' @return \code{defaultHistobounds} returns a vector containing the left and right boundaries for the x axis.
#' 
#' @return \code{plotHistogram} does not return a value, but plots the histogram.
#' 
#' @seealso \code{\link{hwx.test}}


#' @rdname histogramFunctions
#' @export
defaultHistobounds <- 
function(ostats, statID, m) {
	k <- length(m);
	b <- double(2);
	n <- sum(m)/2;
	df <- k *(k-1)/2;
	if(statID==1 || statID==4) b[2] <- qchisq(.999, df);
	if(statID==2) {
		# find maximum probability
		ae <- matrix(0,k,k);
		for(i in 1:k) for(j in 1:k) ae[i,j] <- m[i]*m[j]/(2*n);
		for(i in 1:k)ae[i,i] <- ae[i,i]/2;
		b[1] <- (-2)*log(observedProb(ae));
		b[2] <- b[[1]] + qchisq(.999, df);
	}
	if(statID==3) {
		seu <- sqrt(n * (k-1));
		b[[1]] <- -4 * seu;
		b[[2]] <- 4 * seu;
	}
	return(b);	
}


#' @rdname histogramFunctions
#' @export
plotHistogram <- 
function(ostats, statID, m, histobins, histobounds, histoData, showCurve=TRUE, color1="gray40", color2="lightcoral", ntrials=0) {
	k <- length(m);
	labels <- c("-2 ln(LR)", "-2 ln(Probability)", "U Score (test for homozygote excess)", "Pearson's Chi Squared");
	os <- ostats[statID];
	if(statID==1) os <- -2*os;
	if(statID==2) os <- -2*log(os);
	nleft <- as.integer(histobins*(os-histobounds[1])/(histobounds[2]-histobounds[1]));
	if(nleft < 0) nleft <- 0;
	nright <- histobins - nleft + 5;
	if(nright < 0) nright <- 0;
	if(statID==3 && os < 0) {  #swap colors
		temp <- color1;
		color1 <- color2;
		color2 <- temp;
		labels[3] <- "U Score (test for heterozygote excess)"
	}
	colrs <- c(rep(color1, nleft), rep(color2, nright));
	lwd <- 500/histobins
	plot(seq(histobounds[1], histobounds[2], length.out=histobins+1), 
		histoData, 
		type="h",
		xlab=labels[statID],
		ylab="Frequency",
		col=colrs,
		lend=1,
		lwd=lwd);
	if(showCurve && (statID==1 || statID==4)){
		if(ntrials==0) ntrials <- 1;
		dx <- seq(from=histobounds[[1]], to=histobounds[[2]], length.out=histobins);
		dy <- dchisq(dx,k *(k-1)/2);
		lines(dx,dy * (ntrials) * (histobounds[[2]] - histobounds[[1]])/histobins, col="blue", lwd=2)
	}
}
