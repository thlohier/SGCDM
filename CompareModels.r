rm(list=ls(all=TRUE))


log_like <- rep(0,4);



############################################################## Full model ###############################################################

param  <- read.table("Full/parameters.txt",h=TRUE);


## Exclude species which have been rejected in the model checking step
sp_names  <- names(param);
names_del <- c("Car.pra", "Ger.pra", "Luz.cam", "Tri.fra");
to_del    <- numeric();

for(i in 1:length(names_del))
	to_del <- c(to_del, which(sp_names == names_del[i]));

N_species <- (1:length(param))[-to_del];


like_full <- array(NA,c(17,2,length(param)));

for(sp in N_species)
{
	fic_name  <- paste("Full/Biomass/bio_obs_", sp_names[sp], ".txt", sep="");
	B_tmp     <- read.table(fic_name);
	B_obs_mat <- matrix(unlist(B_tmp), ncol=length(B_tmp));
	B_obs     <- list(B_obs_mat[,1:2],B_obs_mat[,3:4]);
	
	fic_name  <- paste("Full/Biomass/bio_real_", sp_names[sp], ".txt", sep="");
	B_tmp     <- read.table(fic_name);
	B_real    <- matrix(unlist(B_tmp), ncol=length(B_tmp));
	
	fic_name  <- paste("Full/Growth/r_", sp_names[sp], ".txt", sep="");
	r_tmp     <- read.table(fic_name);
	r         <- matrix(unlist(r_tmp), ncol=length(r_tmp));
	
	N_season  <- dim(B_real)[1];
	N_plot    <- dim(B_real)[2];
	
	for(season in 1:N_season)
	{	
		to_keep <- which(!is.na(r[season,]));
		N_p     <- length(to_keep);
		
		if(N_p > 0)
		{
			A     <- (N_p*param[[sp]][3]^2 + param[[sp]][4]^2)/((param[[sp]][3]*param[[sp]][4])^2);
			B     <- (param[[sp]][3]^2)*sum(r[season,to_keep] - param[[sp]][1])/(N_p*param[[sp]][3]^2 + param[[sp]][4]^2);
			C     <- (1/(param[[sp]][4]^2))*sum((r[season,to_keep] - param[[sp]][1])^2);
		
			log_like[1] <- log_like[1] - log((2*pi)^(N_p/2)) - log(sqrt(A)) - log(param[[sp]][3]*(param[[sp]][4]^N_p)) - (1/2)*(C - A*B^2);
		}
		
		for(pl in 1:N_plot)
		{
			to_keep <- which(!is.na(B_obs[[pl]][season,]) & B_obs[[pl]][season,] != 0);
			
			if(length(to_keep) > 0)
			{
				if(B_real[season,pl] != 0 && !is.na(B_real[season,pl]))
					log_like[1]  <- log_like[1] + sum(log(dnorm((log(B_obs[[pl]][season,to_keep]) - log(B_real[season,pl])), 0, param[[sp]][5])),na.rm=TRUE);
			}
		}
	}	
}



######################################################## Model with sigma_d = 0 #########################################################

B_real_tmp <- read.table("Models/M2/bio_real.txt", h=TRUE);
B_real     <- matrix(unlist(B_real_tmp), ncol=length(B_real_tmp));
r          <- read.table("Models/M2/r.txt", h=TRUE);
param      <- read.table("Models/M2/parameters.txt",h=TRUE);


for(sp in N_species)
{
	fic_name  <- paste("Full/Biomass/bio_obs_", sp_names[sp], ".txt", sep="");
	B_tmp     <- read.table(fic_name);
	B_obs     <- matrix(unlist(B_tmp), ncol=length(B_tmp));
	
	N_season  <- dim(B_obs)[1];
	
	for(season in 1:N_season)
	{
		if(!is.na(r[[sp]][season]))
			log_like[2]  <- log_like[2] + log(dnorm((r[[sp]][season] - param[[sp]][1]), 0, param[[sp]][3]));
		
		
		to_keep <- which(!is.na(B_obs[season,]) & B_obs[season,] != 0);
		
		if(length(to_keep) > 0)
		{
			if(B_real[season,sp] != 0 && !is.na(B_real[season,sp]))
				log_like[2]  <- log_like[2] + sum(log(dnorm((log(B_obs[season,to_keep]) - log(B_real[season,sp])), 0, param[[sp]][4])),na.rm=TRUE);
		}
	}
}



################################################## Model with sigma_d = 0 & sigma_e = 0 #################################################

param  <- read.table("Models/M1/parameters.txt",h=TRUE);


for(sp in N_species)
{
	fic_name  <- paste("Full/Biomass/bio_obs_", sp_names[sp], ".txt", sep="");
	B_tmp     <- read.table(fic_name);
	B_obs     <- matrix(unlist(B_tmp), ncol=length(B_tmp));
	
	N_season  <- dim(B_obs)[1];
	
	for(season in 1:N_season)
	{
		to_keep <- which(!is.na(B_obs[season,]) & B_obs[season,] != 0);
		
		if(length(to_keep) > 0)
		{
			if(param[[sp]][1] != 0 && !is.na(param[[sp]][1]))
				log_like[3]  <- log_like[3] + sum(log(dnorm((log(B_obs[season,to_keep]) - log(param[[sp]][1])), 0, param[[sp]][4])),na.rm=TRUE);
		}
	}
}



############################### Model with sigma_d = 0 & sigma_e = 0 & all species are assumed identical ################################

param  <- read.table("Models/M0/parameters.txt",h=TRUE);


for(sp in N_species)
{
	fic_name  <- paste("Full/Biomass/bio_obs_", sp_names[sp], ".txt", sep="");
	B_tmp     <- read.table(fic_name);
	B_obs     <- matrix(unlist(B_tmp), ncol=length(B_tmp));
	
	N_season  <- dim(B_obs)[1];
	
	for(season in 1:N_season)
	{
		to_keep <- which(!is.na(B_obs[season,]) & B_obs[season,] != 0);
		
		if(length(to_keep) > 0)
		{
			if(param[[sp]][1] != 0 && !is.na(param[[sp]][1]))
				log_like[4]  <- log_like[4] + sum(log(dnorm((log(B_obs[season,to_keep]) - log(param[[sp]][1])), 0, param[[sp]][4])),na.rm=TRUE);
		}
	}
}



################################################# Compute AIC and DAIC from likelihood ##################################################

## Number of parameters to estimate for each model
N_sp     <- c(rep(length(N_species),3),1);
nb_param <- c(5,4,3,3);
k_m      <- N_sp*nb_param;

## Aikake information criterion (AIC)
AIC  <- 2*k_m - 2*log_like;



#################################### Write the results of model comparison in a .csv file (Tab. S2) #####################################

to_write <- 1;

if(to_write)
{
	## Create the directory if it does not exixt yet
	main_dir        <- getwd();
	sub_dir         <- "Tables";

	if(!file.exists(paste(main_dir, sub_dir, sep = "/", collapse = "/")))
		dir.create(file.path(main_dir, sub_dir), showWarnings = FALSE);

	
	## Write Tab. S2 in a .csv file
	tab_S2           <- rbind(log_like, k_m, AIC);
	colnames(tab_S2) <- c("Full","M2","M1","M0");
	rownames(tab_S2) <- c("L","k","AIC");

	if(file.exists("Tables/Table_S2.csv"))
	{
		x <- readline("The file Table_S2.csv already exists. Do you want to overwritte it ? (y/n):")
		
		if(x == "y")
			write.csv2(round(tab_S2,1), "Tables/Table_S2.csv", row.names=TRUE);
	}
	else
		write.csv2(round(tab_S2,1), "Tables/Table_S2.csv", row.names=TRUE);
}
