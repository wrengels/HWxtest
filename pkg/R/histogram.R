# plot(seq(0, 16.3, length.out = 500),c$histoData, type = "h", xlab="The Test Statistic", col=c(rep("gray20", 200), rep("lightcoral", 400)))

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