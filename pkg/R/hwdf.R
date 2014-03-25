# converts results of hw.test into a list of hwtest objects
#' Convert results of \code{\link{hw.test}} to a single list of \code{hwtest} objects.
#' 
#' There are two main uses of \code{listify}. You can simplify a complex result from \code{\link{hw.test}} containing multiple populations and multiple loci into a simple list of \code{hwtest} objects. At the same time, you have a chance to change the parameters \code{detail} and \code{statName}. Useful to get output from a test.
#' 
#' @param hwlist the results of a call to \code{\link{hw.test}}. It can be an \code{hwtest} object, a list of them or a list of lists of them.
#' @param detail Used only if you wish to reset the \code{detail} of each object.
#' @param statName Used only if you want to rest the \code{statName} of each object
#' 
#' 
#' @return a list of \code{hwtest} objects, possibly with their \code{detail} and \code{statName} parameters reset
#' 
#' @examples
#' data(HWcases)
#' outcome <- hw.test(HWcases, detail=4, statName="LLR")
#' listify(outcome, detail=1, statName="U")

#' @export
listify <- function(hwlist, detail = NA, statName = NA) {
	# If it's a plain old hwtest, make it a list
	if (class(hwlist) == "hwtest") 
		hwlist <- list(result = hwlist)
	# It should now be a list
	if (class(hwlist) != "list") 
		stop("Format must be list of hwtest")
	# If it is a list of lists, flatten it by 1
	if (all(lapply(hwlist, class) == "list", na.rm = T)) 
		hwlist <- unlist(hwlist, recursive = F)
	# All elements should now be hwtest
	if (any(lapply(hwlist, class) != "hwtest", na.rm = T)) 
		stop("Format must be a list of hwtest")
	if (!is.na(detail)) 
		hwlist <- lapply(hwlist, function(z) {
			z$detail <- detail
			z
		})
	if (!is.na(statName)) {
		statNames <- c("LLR", "Prob", "U", "Chisq");
		statID <- which(statNames==statName);
		hwlist <- lapply(hwlist, function(z) {
			z$statName <- statName
			z$p.value <- z$Pvalues[[statID]]
			z
		})
	}
	isnum <- sapply(hwlist, function(x) is.numeric(x$p.value))
	hwlist <- hwlist[isnum]
	hwlist
}


#' Construct a data frame from \code{\link{hw.test}} output
#' 
#' If the \code{\link{hw.test}} output has multiple populations and/or multiple loci, use this function to make a data frame to display the results in tabular form. 

#' 
#' @param hwlist The output from a call to \code{\link{hw.test}}
#' @param statName gives you the option of changing which statistic's P value is reported
#' @param showN whether to show a column of sample size (number of diploids in the sample)
#' @param showk whether to show the number of alleles
#' @param showMethod whether to show whether the exact or Monte Carlo method was used
#' @param showSE whether to include the standard error for those tests which used the Monte Carlo method
#' @param showTables whether to show the total number of tables examined when full enumeration (exact) method is used
#' @param showTrials whether to show the number of random trials when Monte Carlo method is used
#' @param showStat whether to show the observed statistic

#' @export
hwdf <- 
function(hwlist, statName=NA, showN=TRUE, showk=TRUE, showMethod=TRUE, showSE=TRUE, showTables=TRUE, showTrials=TRUE, showStat=TRUE) {
	hwlist <- listify(hwlist, statName=statName)
	statNames <- c("LLR", "Prob", "U", "Chisq")
	statID <- which(statNames==hwlist[[1]]$statName)
	P.Value <- unlist(sapply(hwlist, function(x) x$Pvalues[[statID]]))
	method <- sapply(hwlist, function(x) x$method)
	exact <- method=="exact"
	monte <- method=="monte"
	ntrials <- sapply(hwlist, function(x) x$ntrials)
	ntrials <- unlist(ifelse(monte, ntrials, NA))
	observedStat <- unlist(sapply(hwlist, function(x) x$observed[[statID]]))
	tableCount <- sapply(hwlist, function(x) x$tableCount)
	N <- sapply(hwlist, function(x) sum(x$alleles)/2)
	k <- sapply(hwlist, function(x) length(x$alleles))
	denom <- ifelse(monte, ntrials, 1)
	denom <- unlist(denom)
	se.all <- sqrt(abs(P.Value * (1-P.Value)/denom))
	SE <- ifelse(monte, se.all, NA)
	pm <- ifelse(monte, "\u00b1","")
	columns <- list(P.Value)
	names <- paste("P value (", statNames[[statID]],")", sep="")
	if(showSE && any(monte)){
		names <- c(names," ", "\u00b1 SE")
		columns$pm <- pm
		columns$SE <- SE
	}
	if(showStat) {
		names <- c(names, paste("observed ", statNames[[statID]], sep=""))
		columns$observedStat <- observedStat
	}
	if(showMethod){
		names <- c(names, "Method")
		columns$method <- unlist(method)
	}
	if(showTables && any(method=="exact")){
		names <- c(names, "All Tables")
		columns$Tables <- unlist(ifelse(exact, tableCount, NA))
	}
	if(showTrials && any(method=="monte")){
		names <- c(names, "Random Trials")
		columns$Trials <- ntrials
	}
	if(showN){
		names <- c(names, "N Diploids")
		columns$N <- N
	}
	if(showk) {
		names <- c(names, "Alleles")
		columns$k <- k
	}
	df <- as.data.frame(columns)
	names(df) <- names
	names(columns) <- names
	df
	
	
}