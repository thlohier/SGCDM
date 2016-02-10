extractTraits <- function(inTrait, inNames)
{

	#######################################################################################################################################
	#																																	  #
	#			Extract the values of the traits "inTrait" corresponding to the species whose name is included in "inNames". 	 		  #
	#																																	  #
	#######################################################################################################################################
	
	
	sp_names  <- numeric();
	trait     <- numeric();
	
	for(sp in 1:length(inTrait$species))
	{
		tmp           <- c(inTrait$genus[sp],inTrait$species[sp]);
		sp_name       <- paste(substr(tmp[1],1,3),".",substr(tmp[2],1,3),sep="");
		sp_names      <- c(sp_names,sp_name);
	}
	
	for(sp in 1:length(inNames))
	{
		ind <- which(sp_names_jena[sp]==sp_names);
		
		if(length(ind)!=0)
		{
			trait_tmp    <- numeric();

			for(i in ind)
				trait_tmp <- c(trait_tmp,inTrait$value[i]);
			
			trait <- c(trait,mean(trait_tmp,na.rm=TRUE));
		}
		else
			trait <- c(trait,NA);
	}
	
	return(trait);
}
