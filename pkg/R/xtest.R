xtesttest <- 
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