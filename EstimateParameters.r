rm(list=ls(all=TRUE))


#######################################################################################################################################
#																																	  #
#		Compute the full model parameters (r_m, K, sigma_e, sigma_d, sigma_o) as well as the environmental, demographic and 		  #
#       observation variables u_e(T), u_d(T,p) and u_o(T,p,s).																		  #
#																																	  #
#######################################################################################################################################


## Biomass in the beginning of the growing season
B_i     <- 2;


######################################## Extract biomass data for species grown in monoculture ########################################

info       <- read.table("Data/info_small.csv",sep=";",h=T,stringsAsFactors=F);
bio        <- read.table("Data/biomass_mono.csv",sep=";",h=T,stringsAsFactors=F);

headers    <- names(info);
plotcode   <- as.vector(unique(bio[[1]]));

mono       <- matrix(0,length(plotcode),2);

for(p in 1:length(plotcode))
{
	t2_info  <- as.vector(info[info$plot==plotcode[p],]);
	
	if(length(which(t2_info==1))==0)
		mono[p,] <- c("",plotcode[p])
	else
		mono[p,] <- c(headers[which(t2_info==1)],plotcode[p]);
}

mono     <- mono[order(mono[,1]),];
sp_names <- unique(mono[,1])[-1];
biomass  <- list();

for(sp in 1:length(sp_names))
{
	bio_sp   <- list();
	plot_cur <- mono[mono[,1]==sp_names[sp],2];
	for(p in 1:length(plot_cur))
	{
		t2_bio   <- as.vector(bio[bio$plot==plot_cur[p],]);
		year     <- t2_bio$year;
		harvest  <- t2_bio$harvest;
		subsamp  <- t2_bio$subsample;
		bio_cur  <- matrix(NA,max(harvest),max(subsamp));
		
		for(i in 1:length(harvest))
			bio_cur[harvest[i],subsamp[i]] <- t2_bio$mass.gm2[i];
		
		bio_sp <- c(bio_sp,list(bio_cur));
	}
	
	biomass <- c(biomass,list(bio_sp));
}

names(biomass) <- sp_names;



########################## Computation of observation variance (sigma_o) and observation variables u_o(T,p,s) #########################

N_species <- length(biomass);

bio_real  <- list();
s_o       <- numeric();
u_o       <- list();


for(sp in 1:N_species)
{
	N_plot       <- length(biomass[[sp]]);
	B_real       <- matrix(NA,17,N_plot);
	e_obs_tmp    <- numeric()
	e_obs        <- list();
	u_o_tmp      <- array(NA,c(17,N_plot,2));
	
	for(p in 1:N_plot)
	{
		N_season     <- dim(biomass[[sp]][[p]])[1];
		N_subsample  <- dim(biomass[[sp]][[p]])[2];
		e_obs_p      <- matrix(NA,N_season,N_subsample);
		s_o_tmp      <- matrix(NA,N_season,N_plot);
		
		for(season in 1:N_season)
		{
			bio_tmp <- biomass[[sp]][[p]][season,];
			to_keep <- which(bio_tmp!=0 & !is.na(bio_tmp));
			
			if(length(to_keep)!=0)
			{
				B_real[season,p]        <- exp(mean(log(bio_tmp[to_keep])));
				e_obs_p[season,to_keep] <- log(bio_tmp[to_keep]) - log(B_real[season,p]);
			}
			else
			{
				B_real[season,p] <- NA;
				e_obs_p[season,] <- rep(NA,N_subsample);
			}
		}
		
		e_obs_tmp <- c(e_obs_tmp,e_obs_p);
		e_obs     <- c(e_obs,list(e_obs_p));
	}
	
	s_o <- c(s_o,sd(e_obs_tmp,na.rm=TRUE));
	
	for(p in 1:N_plot)
		u_o_tmp[1:dim(e_obs[[p]])[1],p,] <- e_obs[[p]]/s_o[sp];
	
	u_o <- c(u_o,list(u_o_tmp));
	
	bio_real <- c(bio_real,list(B_real));
}



############################################### Computation of the carrying capacity K  ###############################################

K         <- numeric();

for(sp in 1:N_species)
	K <- c(K,2*max(bio_real[[sp]],na.rm=T));



################################################ Computation of the growth rate r(T,p) ################################################

N_days  <- 60;
r       <- list();

## The function to minimize is the absolute value of the difference between simulated and observed biomass 
computeR <- function(inR)
{	
	B_prec <- B_i;
	
	for(day in 1:N_days)
	{
		B_cour <- B_prec + inR*B_prec*(1-(B_prec/K[sp]));
		B_prec <- B_cour;
	}
	
	R <- abs(B_emp - B_cour);
}


for(sp in 1:N_species)
{
	N_season <- dim(bio_real[[sp]])[1];
	N_plot   <- dim(bio_real[[sp]])[2];
	
	r_spec <- matrix(NA,N_season,N_plot);
	
	for(season in 1:N_season)
	{
		for(p in 1:N_plot)
		{
			B_emp <- bio_real[[sp]][season,p];
				
			if(!is.na(B_emp))
			{
				if(B_emp > B_i)
				{
					fit_r <- optimize(f=computeR,interval=c(-0.1,0.5),tol=1e-20)
					
					r_spec[season,p] <- fit_r$minimum;
				}
				else
					r_spec[season,p] <- 0;
			}
				
		}
	}
	
	r <- c(r,list(r_spec));
}



######################################## Verification of r values founded ########################################

bio_test <- list();
error    <- numeric();

for(sp in 1:N_species)
{	
	B_test <- matrix(-1,N_season,N_plot);
	
	for(p in 1:N_plot)
	{
		for(season in 1:N_season)
		{
			B_prec <- B_i;
			
			for(day in 1:N_days)
			{
				B_cour <- B_prec + r[[sp]][season,p]*B_prec*(1-(B_prec/K[sp]));
				B_prec <- B_cour;
			}
			
			B_test[season,p] <- B_cour;
		}
	}
	
	bio_test <- c(bio_test,list(B_test));
}

for(sp in 1:N_species)
	error <- c(error,max(abs(bio_real[[sp]] - bio_test[[sp]]),na.rm=T));



######################################## Estimation of model parameters (r_m, sigma_e, sigma_d) #######################################

# Esimated model parameters
r_m  <- numeric(0);
s_e  <- numeric(0);
s_d  <- numeric(0); 
u_e  <- matrix(0,N_season,N_species);
ac_e <- numeric(0);

## Log likelihood function
F <- function(inP)
{
	R     <- numeric(0);

	for(season in 1:N_season)
	{
		r_use <- r[[sp]][season,][!is.na(r[[sp]][season,])];
		N_p   <- length(r_use);
		
		if(N_p>0)
		{
			t1 <- (N_p*inP[2]^2+inP[3]^2)/((inP[2]*inP[3])^2);
			t2 <- (inP[2]^2)*sum(r_use - inP[1])/(N_p*inP[2]^2 + inP[3]^2);
			t3 <- (1/(inP[3]^2))*sum((r_use - inP[1])^2);
			
			R   <- c(R,(-log((2*pi)^(N_p/2)) - log(sqrt(t1)) - log(inP[2]*(inP[3]^N_p)) - (1/2)*(t3 - t1*t2^2)));
		}
	}
	
	return(-sum(R))
}


## Estimation of model parameters
for(sp in 1:N_species)
{
	## Maximum likelihood estimation with quasi-Newton optimization algorithm
	fit   <- optim(c(0.05,0.025,0.025),F);
	
	## Model parameters
	r_m <- c(r_m,fit$par[1]);
	s_e <- c(s_e,fit$par[2]);
	s_d <- c(s_d,fit$par[3]);
}



############################################ Estimation of environmental variables u_e(T) #############################################

u_e     <- matrix(NA,17,N_species);
u_d     <- list();
error_r <- numeric();

G <- function(inU_e)
{
	R <- (1/2)*inU_e^2;
	for(p in 1:N_p)
		R <- R + (1/2)*((r_use[p] - r_m[sp] - s_e[sp]*inU_e)/s_d[sp])^2;
	
	return(R);
}

for(sp in 1:N_species)
{
	for(season in 1:N_season)
	{
		r_use <- r[[sp]][season,][!is.na(r[[sp]][season,])];
		N_p   <- length(r_use);
		
		if(N_p>0)
		{
			fit_env <- optimize(f=G,interval=c(-10,10),tol=1e-20)
					
			u_e[season,sp] <- fit_env$minimum;
		}
	}
	
	u_d <- c(u_d,list((r[[sp]] - r_m[sp] - s_e[[sp]]*u_e[,sp])/s_d[sp]));
	
	r_test <- matrix(NA,dim(u_d[[sp]])[1],dim(u_d[[sp]])[2]);
	for(p in 1:dim(u_d[[sp]])[2])
		r_test[,p]  <- r_m[sp] + s_e[sp]*u_e[,sp] + s_d[sp]*u_d[[sp]][,p];
	
	error_r <- c(error_r,max(abs(r[[sp]] - r_test),na.rm=T));
}



###################################### Save model parameters and simulation details in .txt files #####################################

to_write <- 0;

if(to_write==TRUE)
{
	main_dir        <- getwd();
	sub_dir         <- "Full";
	dir.create(file.path(main_dir, sub_dir), showWarnings = FALSE);
	setwd(file.path(main_dir, sub_dir));
	
	sec_dir         <- getwd();
	sub_dir_1       <- "Growth";
	dir.create(file.path(sec_dir, sub_dir_1), showWarnings = FALSE);
	sub_dir_2       <- "Demo";
	dir.create(file.path(sec_dir, sub_dir_2), showWarnings = FALSE);
	sub_dir_3       <- "Obs";
	dir.create(file.path(sec_dir, sub_dir_3), showWarnings = FALSE);
	sub_dir_4       <- "Biomass";
	dir.create(file.path(sec_dir, sub_dir_4), showWarnings = FALSE);
	
	setwd(main_dir);
	
	
	param           <- cbind(r_m,K,s_e,s_d,s_o);
	rownames(param) <- sp_names;
	colnames(u_e)   <- sp_names;
	
	
	for(sp in 1:N_species)
		write.table(round(r[[sp]],8),paste("Full/Growth/r_",sp_names[sp],".txt",sep=""), sep="\t", row.names=FALSE, col.names=FALSE);
	
	
	write.table(round(t(param),8), "Full/parameters.txt", sep="\t", row.names=FALSE, col.names=TRUE);
	write.table(round(u_e,8),"Full/env.txt", sep="\t", row.names=FALSE, col.names=TRUE);
	
	
	for(sp in 1:N_species)
	{
		write.table(round(u_d[[sp]],8),paste("Full/Demo/demo_",sp_names[sp],".txt",sep=""), sep="\t", row.names=FALSE, col.names=FALSE);
		
		u_o_mat <- cbind(u_o[[sp]][,1,],u_o[[sp]][,2,]);
		write.table(round(u_o_mat,8),paste("Full/Obs/obs_",sp_names[sp],".txt",sep=""), sep="\t", row.names=FALSE, col.names=FALSE);		
	}
	
	
	for(sp in 1:N_species)
	{
		bio_tmp       <- matrix(NA,17,4);
		bio_tmp[1:dim(biomass[[sp]][[1]])[1],1:2] <- biomass[[sp]][[1]];
		bio_tmp[1:dim(biomass[[sp]][[2]])[1],3:4] <- biomass[[sp]][[2]];
		
		write.table(round(bio_tmp,8),paste("Full/Biomass/bio_obs_",sp_names[sp],".txt",sep=""), sep="\t", row.names=FALSE, col.names=FALSE);
		
		write.table(round(bio_real[[sp]],8),paste("Full/Biomass/bio_real_",sp_names[sp],".txt",sep=""), sep="\t", row.names=FALSE, col.names=FALSE);
	}
}
