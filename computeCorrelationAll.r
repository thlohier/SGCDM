computeCorrelationAll <- function(inTrait, inParam)
{

	#######################################################################################################################################
	#																																	  #
	#		Compute the correlations coefficients between the trait "inTrait" and each parameter included in "inParam" as well as		  #
	#		the significance of these correlations.																				 		  #
	#																																	  #
	#######################################################################################################################################
	
	
	R      <- numeric();
	p_val  <- numeric();
	
	to_del <- which(is.na(inTrait));
	
	if(length(to_del) > 0)
	{
		trait <- inTrait[-to_del];
		param <- list(inParam[1,-to_del],inParam[2,-to_del],inParam[3,-to_del],inParam[4,-to_del],inParam[5,-to_del]);
	}
	else
	{
		trait <- inTrait;
		param <- list(inParam[1,],inParam[2,],inParam[3,],inParam[4,],inParam[5,]);
	}
		
	for(i in 1:length(param))
	{
		R     <- c(R,round(cor(param[[i]],trait),2));
		p_val <- c(p_val, round(cor.test(param[[i]],trait)[[3]],3));
	}
	
	return(cbind(R,p_val));
}
