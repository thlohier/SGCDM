mvnrnd <- function(inMu, inSigma, inPhi, inCases)
{	
	
	#######################################################################################################################################
	#																																	  #
	#				Returns an n-by-d matrix "r" of random vectors chosen from the multivariate normal distribution with				  #
	#				mean "inMu", and covariance "inSigma". 	 																			  #
	#																																	  #
	#######################################################################################################################################
	
	
	d <- dim(inMu)[1];
	
	T <- chol(inSigma);
	M <- matrix(inMu,ncol=d,nrow=inCases,byrow=TRUE);
	
	tmp <- matrix(0,(inCases+200),d);
	V   <- (1 - inPhi^2)*1^2;

	for(sp in 1:d)
	{
		eps <- rnorm((inCases+200),0,sqrt(V));
		
		tmp[1,sp] <- 0;
		for(i in 2:length(eps))
			tmp[i,sp] <- inPhi*tmp[i-1,sp] + eps[i];
	}
	
	G <- tmp[201:(inCases+200),];
	
	r <- G%*%T + M;
	
	return(r);
}
