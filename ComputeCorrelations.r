rm(list=ls(all=TRUE))

source("ExtractTraits.r");
source("ComputeCorrelationAll.r");


#######################################################################################################################################
#																																	  #
#		Compute the correlations between the different species model parameters and between species model parameters and  	 		  #
#		functional traits when all species are pooled together and when they are sorted according to their functional groups		  #
#																																	  #
#######################################################################################################################################


## Read species parameters
param_tmp     <- read.table("Full/parameters.txt",h=T);
param         <- matrix(unlist(param_tmp),ncol=length(param_tmp));
sp_names      <- names(param_tmp);


## Exclude species which have been rejected in the model checking step
sp_names  <- names(param_tmp);
names_del <- c("Car.pra", "Ger.pra", "Luz.cam", "Tri.fra");
to_del    <- numeric();

for(i in 1:length(names_del))
	to_del <- c(to_del, which(sp_names == names_del[i]));

N_species     <- (1:length(param))[-to_del];
sp_names_jena <- sp_names[-to_del];
param_clean   <- param[,-to_del];



######################################### Relationship between model parameters and leaf traits #########################################

sla    <- read.table("Data/leda/SLA.txt",h=T,stringsAsFactors=F);
height <- read.table("Data/leda/canopy_height.txt",h=T,stringsAsFactors=F);
ldmc   <- read.table("Data/leda/LDMC.txt",h=T,stringsAsFactors=F);


H    <- extractTraits(height, sp_names_jena);
SLA  <- extractTraits(sla, sp_names_jena);
LDMC <- extractTraits(ldmc, sp_names_jena);


cor_traits_all <- cbind(computeCorrelationAll(H,param_clean),computeCorrelationAll(SLA,param_clean),computeCorrelationAll(LDMC,param_clean));



####################### Relationship between model parameters and leaf traits according to plant functional types #######################

pft         <- read.table("Data/pft.csv",sep=";");
plant_types <- unique(as.vector(pft[[3]]));

param_pft   <- list();
SLA_pft     <- list();
LDMC_pft    <- list();
H_pft       <- list();


for(i in 1:length(plant_types))
{
	type_sp   <- which(as.vector(pft[[3]])==plant_types[i]);
	ind_param <- numeric();
	
	for(j in 1:length(type_sp))
		ind_param <- c(ind_param,which(sp_names_jena==pft[[1]][type_sp[j]]));
	
	
	param_pft <- c(param_pft,list(rbind(param_clean[1,ind_param],param_clean[2,ind_param],param_clean[3,ind_param],param_clean[4,ind_param],param_clean[5,ind_param])));
	SLA_pft   <- c(SLA_pft,list(SLA[ind_param]));
	LDMC_pft  <- c(LDMC_pft,list(LDMC[ind_param]));
	H_pft     <- c(H_pft,list(H[ind_param]));
}

cor_traits_pft   <- matrix(0,0,6);

for(i in 1:length(param_pft))
	cor_traits_pft <- rbind(cor_traits_pft,cbind(computeCorrelationAll(H_pft[[i]],param_pft[[i]]),computeCorrelationAll(SLA_pft[[i]],param_pft[[i]]),computeCorrelationAll(LDMC_pft[[i]],param_pft[[i]])));



########################################## Relationship between the different model parameters ##########################################

cor_param_all <- matrix(NA,dim(param)[1],dim(param)[1]);
p_param_all   <- matrix(NA,dim(param)[1],dim(param)[1]);


for(i in 1:dim(param)[1])
{
	for(j in 1:dim(param)[1])
	{
		cor_param_all[i,j] <- round(cor(param_clean[i,],param_clean[j,]),2);
		p_param_all[i,j]   <- round(cor.test(param_clean[i,],param_clean[j,])[[3]],3);
	}
}



################################## Relationship between the different model parameters according to PFT #################################

cor_param_pft <- list();
p_param_pft   <- list();

for(pft in 1:length(param_pft))
{
	cor_param_tmp <- matrix(NA,dim(param)[1],dim(param)[1]);
	p_param_tmp   <- matrix(NA,dim(param)[1],dim(param)[1]);
		
	for(i in 1:dim(param)[1])
	{
			for(j in 1:dim(param)[1])
			{
				cor_param_tmp[i,j] <- round(cor(param_pft[[pft]][i,],param_pft[[pft]][j,]),2);
				p_param_tmp[i,j]   <- round(cor.test(param_pft[[pft]][i,],param_pft[[pft]][j,])[[3]],3);
			}
	}
	
	cor_param_pft <- c(cor_param_pft, list(cor_param_tmp));
	p_param_pft   <- c(p_param_pft,list(p_param_tmp));
}



###################################### Write the correlation coefficient and p-values in .csv file ######################################

to_write <- 1;

if(to_write)
{
	## Correlation between model parameters and plant functional traits
	tab_2     <- matrix(NA,0,6);
	order_tab <- c(0,5,10,15);
	
	for(i in 1:dim(cor_traits_all)[1])
		tab_2 <- rbind(tab_2,cor_traits_all[i,],cor_traits_pft[(order_tab+i),]);
		
	write.csv2(tab_2,"Tables/Table_2.csv");
	
	
	## Correlation between the different model parameters
	tab_S3    <- matrix(NA,0,6);
	order_row <- c(1,3,4);
	order_col <- c(3,4,5);
	
	for(i in 1:length(order_row))
	{
		line_tmp <- numeric();
		
		for(j in 1:length(order_col))
			line_tmp <- c(line_tmp,cor_param_all[order_row[i],order_col[j]],p_param_all[order_row[i],order_col[j]]);
		
		tab_S3 <- rbind(tab_S3,line_tmp);
			
		
		for(pft in 1:length(cor_param_pft))
		{
			line_pft_tmp <- numeric();
			
			for(j in 1:length(order_col))
				line_pft_tmp <- c(line_pft_tmp,cor_param_pft[[pft]][order_row[i],order_col[j]],p_param_pft[[pft]][order_row[i],order_col[j]]);
			
			tab_S3 <- rbind(tab_S3,line_pft_tmp);
		}
	}
	
	write.csv2(tab_S3,"Tables/Table_S3.csv", row.names=FALSE);
}
