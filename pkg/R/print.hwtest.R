# For printing the output of hw.test
# (c) William R. Engels, 2014


#' S3 Method for printing \code{hwtest} objects

#' 
#' Prints test results (\code{hwtest}) objects depending on how much detail is provided.
#' 
#' @param h the results from a call to \code{\link{hw.test}}
#' @param detail If unspecified, the level of detail used in the \code{hw.test} call will be applied.
#' @param ... other parameters passed to \code{print}


#' @method print hwtest
#' @S3method print hwtest
print.hwtest <- 
function(h, detail=h$detail, ...) {
	if(length(h) < 8) return();
	if(h$method=="exact" && h$tableCount < 0) stop("Calculation timed out. You can change the time limit by setting parameter 'safeSecs'");
	statNames <- c("LLR", "Prob", "U", "Chisq");
	statID <- which(statNames==h$statName);
	p <- h$Pvalues;
	ob <- h$observed;
	if(h$method=="monte") se <- sqrt(p * (1-p)/h$ntrials);
	comments <- character(4);
	if(ob[3]>=0) {comments[[3]] <- " (test for homozygote excess)"}
	if(ob[3]<0) {comments[[3]] <- " (test for heterozygote excess)"};
	if(detail==1) {
		cat("\nP value (", h$statName,") = ", formatC(p[statID]), sep="");
		if(h$method=="monte") cat( " \u00b1 ",formatC(se[statID], digits=5), sep="");
		cat(comments[statID],"\n");
	}
	if(detail >= 2) {
		if(h$method=="monte") cat("\nMonte Carlo test for HW with ", h$ntrials," trials.\n", sep="");
		if(h$method=="exact") cat("\nFull enumeration of ",h$tableCount, " tables to test for HW\n", sep="");
		for(i in 1:4){
			cat("\nP value (",statNames[i],")",strtrim("     ", 6-nchar(statNames[i])),"= ", formatC(p[i], digits=6, format="f"), sep="");
			if(h$method=="monte") cat( " \u00b1 ",formatC(se[i], digits=5), sep="");
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
		print(h$genotypes, na.print="");
	}
	if(detail>=5){
		cat("\n\nExpected Genotype Counts\n");
		print(observedX2(h$genotypes, returnExpected=T), digits=3, na.print="");
	}
	cat("\n", sep="");
}