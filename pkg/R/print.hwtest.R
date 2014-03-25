# For printing the output of hw.test
# (c) William R. Engels, 2014


#' S3 Method for printing \code{hwtest} objects

#' 
#' Prints test results (\code{hwtest}) objects depending on how much detail is provided.
#' 
#' @param x the results from a call to \code{\link{hw.test}}
#' @param ... other parameters passed to \code{print}. You can specify \code{detail} or \code{statName}


#' @method print hwtest
#' @S3method print hwtest
print.hwtest <- 
function(x,...) show.hwtest(x,...)


#' Display results from \code{\link{hw.test}}
#' 
#' \code{show.hwtest} is called by \code{print}, but it can also be called directly by the user if the object is a list of \code{hwtest} objects, or even a list of such lists.
#' 
#' 
#' @param h the results of a call to \code{\link{hw.test}}
#' @param detail may be used to change the detail of the output
#' @param statName may be used to specify which statistic to use if \code{detail} is 1.
#' 
#' @return none

#' @export
show.hwtest <- 
function(h, detail=NA, statName=NA) {
	UseMethod("show.hwtest")
}

#' @export
show.hwtest.hwtest <- 
function(h, detail=NA, statName=NA) {
	if(is.na(h$method)) return()
	if(!is.na(detail)) h$detail=detail
	if(!is.na(statName))h$statName=statName
	if(h$method=="exact" && h$tableCount < 0) stop("Calculation timed out. You can change the time limit by setting parameter 'safeSecs'");
	p <- h$Pvalues;
	ob <- h$observed;
	statNames <- names(h$SE)
	comments <- character(4);
	if(ob[3]>=0) {comments[[3]] <- " (test for homozygote excess)"}
	if(ob[3]<0) {comments[[3]] <- " (test for heterozygote excess)"};
	names(comments) <- statNames;
	detail <- h$detail
	if(detail==1) {
		cat("P value (", h$statName,") = ", formatC(p[h$statName]), sep="");
		if(h$method=="monte") cat( " \u00b1 ",formatC(h$SE[h$statName], digits=5), sep="");
		cat(comments[h$statName],"\n");
	}
	if(detail >= 2) {
		cat("\n*****    Sample of ", sum(h$alleles)/2," diploids with ", length(h$alleles), " alleles", sep="")
		if(h$method=="monte") cat("\nMonte Carlo test for HW with ", h$ntrials," trials.\n", sep="");
		if(h$method=="exact") cat("\nFull enumeration of ",h$tableCount, " tables to test for HW\n", sep="");
		for(i in 1:4){
			cat("\nP value (",statNames[i],")",strtrim("     ", 6-nchar(statNames[i])),"= ", formatC(p[i], digits=6, format="f"), sep="");
			if(h$method=="monte") cat( " \u00b1 ",formatC(h$SE[i], digits=5), sep="");
			cat(comments[i]);
		}
	}
	if(detail>=3) {
		cat("\n\nObserved Test Statistics:\n");
		for(i in 1:4){cat("\n",strtrim("     ", 6-nchar(statNames[i])), statNames[i],"  :  ", formatC(ob[i]))}
	}
	if(detail>=4){
		cat("\n\nObserved Allele Counts: ", h$alleles);
		cat("\n\nObserved Genotype Counts\n");
		print(clearUpper(h$genotypes), na.print="");
	}
	if(detail>=5){
		cat("\n\nExpected Genotype Counts\n");
		ecounts <- observedX2(h$genotypes, returnExpected=T);
		rownames(ecounts) <- rownames(h$genotypes);
		colnames(ecounts) <- rownames(h$genotypes);
		print(clearUpper(ecounts), digits=3, na.print="");
	}
	cat("\n", sep="");
}

#' @export
show.hwtest.list <- 
function(h, detail=NA, statName=NA) {
	listify(h, detail=detail, statName=statName) 
}