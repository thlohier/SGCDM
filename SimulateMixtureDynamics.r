rm(list=ls(all=TRUE))


## Generate random vectors chosen from the multivariate normal distribution
source("Tools/mvnrnd.r");


#######################################################################################################################################
#																																	  #
#				Simulate the community dynamics of mixtures varying in their diversity level (2,4,8 and 16 species)					  #
#				The environmental variables u_e(T) are those estimated with the inference method (EstimateParameters.r)				  #
#				The demographic and observation variable u_d(T,p) and u_o(T,p,s) are drawn from multivariate normal distributions	  #
#				For each mixtures, 100 replicated simualtions are run.																  #
#																																	  #
#######################################################################################################################################


## Simulation details
N_sp        <- c(2,4,8,16);
N_season    <- 14;
N_subsample <- 4;
N_days      <- 60;
N_replica   <- 100;


## Read observed biomass, parameters and environmental variables in .txt files
biomass <- read.table("Data/biomass.csv",sep="\t",h=T,stringsAsFactors=F);
param   <- read.table("Full/parameters.txt",sep="\t",h=T,stringsAsFactors=F);
env     <- read.table("Full/env.txt",h=TRUE);


## Mixture codes
plotcode   <- as.vector(unique(biomass$plotcode));


for(n in N_sp)
{
	mix_codes  <- numeric();
	mix_names  <- list();
	ig_names   <- list();
	ignored    <- numeric();

	## Observed statistics
	M_obs     <- numeric();
	SD_obs    <- numeric();
	l_obs     <- numeric();
	
	## Averaged statistics across replicated simulations
	M_sim     <- numeric();
	SD_sim    <- numeric();
	l_sim     <- numeric();
	
	## Statistics for each replicated simulation
	M_rep       <- matrix(NA,N_replica,0);
	SD_rep      <- matrix(NA,N_replica,0);
	l_rep       <- matrix(NA,N_replica,0);
	
	
	for(p in 1:length(plotcode))
	{	
		######################################## Extract biomass data for species grown in mixture ########################################
		
		t1 <- biomass[biomass$plotcode==plotcode[p],];
		
		if(t1$sowndiv[1]==n)
		{
			
			## Names of species present in the mixture
			ind          <- which(t1[1,22:81]==1);
			N_species    <- length(ind);
			sp_names     <- numeric(0);
			
			for(i in 1:N_species)
			{
				name_t   <- colnames(t1)[ind[i]+21];
				sp_names <- c(sp_names,substr(name_t,2,nchar(name_t)));
			}
			
			if("Car.pra" %in% sp_names || "Ger.pra" %in% sp_names || "Luz.cam" %in% sp_names || "Tri.fra" %in% sp_names)
			{
				ignored  <- c(ignored,p);
				ig_names <- c(ig_names,list(sp_names));
			}
			else
			{
				mix_codes <- c(mix_codes,plotcode[p]);
				mix_names <- c(mix_names,list(sp_names));
				
				## Biomass of species present in the mixture
				years        <- unique(t1$year);
				annual_data  <- array(NA,c(N_species,N_season,N_subsample));
				harvest_date <- matrix(NA,50,4);
				
				ind  <- 1;
				ind1 <- 1;
				for(i in 1:length(years))
				{
					while(ind<length(t1$year) && t1$year[ind]==years[i])
					{
						ind2 <- t1$subsample[ind];
						
						if(t1$month[ind]=="May" || t1$month[ind]=="Jun")
						{
							harvest_date[ind,ind2] <- paste(t1$month[ind],t1$year[ind]);
							for(sp in 1:N_species)
								annual_data[sp,ind1,ind2] <- t1[[sp_names[sp]]][ind];
						}
						else if(t1$month[ind]=="Aug" || t1$month[ind]=="Sept")
						{
							harvest_date[ind,ind2] <- paste(t1$month[ind],t1$year[ind]);
							for(sp in 1:N_species)
								annual_data[sp,(ind1+1),ind2] <- t1[[sp_names[sp]]][ind];
						}

						ind  <- ind+1;
					}
					
					ind1 <- ind1+2;
				}
				
				annual_data[annual_data<0] <- NA;
				
				B_tot <- matrix(NA,N_season,N_subsample);
				
				for(subsample in 1:N_subsample)
					B_tot[,subsample] <- apply(annual_data[,,subsample],2,sum);
					
				ind_na <- which(is.na(B_tot),arr.ind=TRUE);
				
				
				
				############################### Compute the community statistics from Jena polyculture biomass data ###############################
				
				M_mix      <- numeric(0);
				SD_mix     <- numeric(0);
				lembda_mix <- numeric(0);
				prop_mix   <- matrix(NA,N_species,N_subsample);
				
				## Compute communtiy statistics for each subsample
				for(subsample in 1:N_subsample)
				{
					M_mix  <- c(M_mix,mean(B_tot[,subsample],na.rm=TRUE));
					SD_mix <- c(SD_mix,sd(B_tot[,subsample],na.rm=TRUE));
					
					for(sp in 1:dim(annual_data)[1])
						prop_mix[sp,subsample] <- mean(annual_data[sp,,subsample],na.rm=TRUE)/M_mix[subsample];
					
					lembda_mix <- c(lembda_mix,sum(prop_mix[,subsample]^2));
				}
				
				## Compute the average community statistics
				M_obs    <- c(M_obs,mean(M_mix,na.rm=TRUE));
				SD_obs   <- c(SD_obs,mean(SD_mix,na.rm=TRUE));
				l_obs    <- c(l_obs,mean(lembda_mix,na.rm=TRUE));
				
				
				
				######################################### Simulate community dynamics with the full model #########################################
				
				M_rep_mix    <- numeric();
				SD_rep_mix   <- numeric();
				l_rep_mix    <- numeric();
				
				## Initial biomass
				B_ini  <- rep(2/N_species,N_species);
				
				## Species parameters (r_m, K, s_e, s_d, s_o), environmental variables (u_e) and competition coefficient (alpha)
				r_m    <- numeric();
				K      <- numeric();
				s_e    <- numeric();
				s_d    <- numeric();
				s_o    <- numeric();
				u_e    <- matrix(NA,N_season,N_species);
				alpha  <- matrix(1,N_species,N_species);
				
				for(sp in 1:N_species)
				{
					r_m <- c(r_m,param[[sp_names[sp]]][1]);
					K   <- c(K,param[[sp_names[sp]]][2]);
					s_e <- c(s_e,param[[sp_names[sp]]][3]);
					s_d <- c(s_d,sqrt(2)*param[[sp_names[sp]]][4]);
					s_o <- c(s_o,param[[sp_names[sp]]][5]);
					
					u_e[,sp] <- env[[sp_names[sp]]][1:N_season];
				}
				
				
				for(replica in 1:N_replica)
				{	
					B       <- array(NA,c(N_species,N_season,N_subsample));
					
					## Generation of random variable with zero mean and unit variance (demographic fluctuations)
					set.seed(replica);
					u_d     <- matrix(NA,N_season,N_species);
					mu_d    <- matrix(0,N_species,1);
					sig_d   <- diag(N_species);
					u_d     <- mvnrnd(mu_d,sig_d,0.0,N_season);
					
					
					## Generation of random variable with zero mean and unit variance (observation error)
					set.seed(N_replica + replica);
					u_o_tmp <- matrix(NA,N_season*N_subsample,N_species);
					mu_o    <- matrix(0,N_species,1);
					sig_o   <- diag(N_species);
					u_o_tmp <- mvnrnd(mu_o,sig_o,0.0,N_season*N_subsample);
					u_o     <- array(u_o_tmp,c(N_season,N_subsample,N_species))
					
					
					for(samp in 1:N_subsample)
					{
						bio  <- array(0,c(N_species,N_days,N_season));
						
						for(season in 1:N_season)
						{
							B_prec <- B_ini
						
							for(day in 1:N_days)
							{		
								for(sp in 1:N_species)
								{	
									## Growth rate of the species
									r    <- max(0,r_m[sp] + s_e[sp]*u_e[season,sp] + (s_d[sp]/sqrt(B_ini[sp]))*u_d[season,sp]);

									## Competition term
									comp <- max(0,(1-((alpha%*%matrix(B_prec,N_species,1))[sp]/K[sp])));
								
									## Biomass in the next day
									bio[sp,day,season] <- B_prec[sp] + B_prec[sp]*r*comp;
									B_prec[sp]         <- bio[sp,day,season];
								}
							}
						
							## Observation error
							non_null                <- which(bio[,N_days,season] > 0)
							B[non_null,season,samp] <- exp(log(bio[non_null,N_days,season]) + s_o[sp]*u_o[season,samp,non_null]);
						}
					}
				
					B_tot <- matrix(NA,N_season,N_subsample);
				
					for(subsample in 1:N_subsample)
						B_tot[,subsample] <- apply(B[,,subsample],2,sum);
					
					B_tot[ind_na] <- NA;
					for(ii in 1:dim(ind_na)[1])
						B[,ind_na[ii,1],ind_na[ii,2]] <- NA;
					
					
					
					##################################### Compute the CV from simulated polyculture biomass data ######################################
				
					M_mix      <- numeric();
					SD_mix     <- numeric();
					lembda_mix <- numeric();
					prop_mix   <- matrix(NA,N_species,N_subsample);
				
					
					## Compute communtiy statistics for each subsample
					for(subsample in 1:N_subsample)
					{
						M_mix   <- c(M_mix,mean(B_tot[,subsample],na.rm=TRUE));
						SD_mix  <- c(SD_mix,sd(B_tot[,subsample],na.rm=TRUE));
					
						for(sp in 1:dim(B)[1])
							prop_mix[sp,subsample] <- mean(B[sp,,subsample],na.rm=T)/M_mix[subsample];
					
						lembda_mix <- c(lembda_mix,sum(prop_mix[,subsample]^2));
					}
					
					
					M_rep_mix  <- c(M_rep_mix,mean(M_mix));
					SD_rep_mix <- c(SD_rep_mix,mean(SD_mix));
					l_rep_mix  <- c(l_rep_mix,mean(lembda_mix));
				}
				
				
				## Stock community statistics obtained for each replica
				M_rep     <- cbind(M_rep,M_rep_mix);
				SD_rep    <- cbind(SD_rep,SD_rep_mix);
				l_rep     <- cbind(l_rep,l_rep_mix);
				
				
				## Compute the average community statistics
				M_sim     <- c(M_sim,mean(M_rep_mix));
				SD_sim    <- c(SD_sim,mean(SD_rep_mix));
				l_sim     <- c(l_sim,mean(l_rep_mix));
			}
		}
	}
	
	
	
	################################# Write the observed and simulated community statistics in .txt files #################################
	
	to_write <- 1;
	
	if(to_write)
	{
		main_dir        <- getwd();
		sub_dir         <- "Full";
		new_dir         <- "Mixtures";
		
		if(!file.exists(paste(main_dir, sub_dir, sep = "/", collapse = "/")))
		{
			dir.create(file.path(main_dir, sub_dir), showWarnings = FALSE);
			setwd(file.path(main_dir, sub_dir));
			
			if(!file.exists(paste(main_dir, sub_dir, new_dir, sep = "/", collapse = "/")))
				dir.create(file.path(main_dir, sub_dir, new_dir), showWarnings = FALSE);
		}
		
		
		## Save community statistics obtained for each replica in .txt files
		
		main_dir        <- getwd();
		sub_dir         <- "Full";
		new_dir         <- "Mixtures";
		dir.create(file.path(main_dir, sub_dir, new_dir), showWarnings = FALSE);
		
		
		## Save community statistics obtained for each replica
		colnames(M_rep) <- mix_codes;
		write.table(round(M_rep,4), paste("Full/Mixtures/mean_rep_",n,".txt",sep=""), sep="\t", row.names=FALSE, col.names=TRUE);

		colnames(SD_rep) <- mix_codes;
		write.table(round(SD_rep,4), paste("Full/Mixtures/deviation_rep_",n,".txt",sep=""), sep="\t", row.names=FALSE, col.names=TRUE);

		colnames(l_rep) <- mix_codes;
		write.table(round(l_rep,4), paste("Full/Mixtures/simpson_rep_",n,".txt",sep=""), sep="\t", row.names=FALSE, col.names=TRUE);
		
		
		## Save averaged community satistics
		stat_obs           <- rbind(M_obs,SD_obs,l_obs);
		rownames(stat_obs) <- mix_codes;
		write.table(round(stat_obs,4), paste("Full/Mixtures/stat_obs_",n,".txt",sep=""), sep="\t", row.names=FALSE, col.names=TRUE);

		stat_sim           <- rbind(M_sim,SD_sim,l_sim);
		rownames(stat_sim) <- mix_codes;
		write.table(round(stat_sim,4), paste("Full/Mixtures/stat_sim_",n,".txt",sep=""), sep="\t", row.names=FALSE, col.names=TRUE);
	}
}
