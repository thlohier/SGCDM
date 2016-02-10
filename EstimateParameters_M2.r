rm(list=ls(all=TRUE))


#######################################################################################################################################
#																																	  #
#								Compute the model M2 parameters (r_mi, K_i, sigma_ei, sigma_oi) 							 		  #
#																																	  #
#######################################################################################################################################


## Biomass at the beginning of the season [g.m-2]
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



######################################## Computation of observation variance ########################################

N_season    <- 17;
N_subsample <- 4;
N_species   <- length(biomass);
bio_real    <- matrix(NA,N_season,N_species);
s_o         <- numeric();


for(sp in 1:N_species)
{	
	bio_tmp       <- matrix(NA,N_season,N_subsample);
	bio_tmp[1:dim(biomass[[sp]][[1]])[1],1:2] <- biomass[[sp]][[1]];
	bio_tmp[1:dim(biomass[[sp]][[2]])[1],3:4] <- biomass[[sp]][[2]];

	e_obs         <- matrix(NA,N_season,N_subsample);
		
	for(season in 1:N_season)
	{
		to_keep <- which(bio_tmp[season,]!=0 & !is.na(bio_tmp[season,]));
		
		if(length(to_keep)!=0)
		{
			bio_real[season,sp]    <- exp(mean(log(bio_tmp[season,to_keep])));
			e_obs[season,to_keep]  <- log(bio_tmp[season,to_keep]) - log(bio_real[season,sp]);
		}
	}
	
	s_o      <- c(s_o,sd(e_obs,na.rm=TRUE));
}



######################################## Computation of K as maximum biomass through plot and year ########################################

K         <- numeric();

for(sp in 1:N_species)
	K <- c(K,2*max(bio_real[,sp],na.rm=TRUE));



######################################## Computation of r(T,p) with a gradient descent ########################################

N_days  <- 60;
r       <- matrix(NA,N_season,N_species);

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
	for(season in 1:N_season)
	{
		B_emp <- bio_real[season,sp];
		
		if(!is.na(B_emp))
		{
			if(B_emp > B_i)
			{
				fit_r        <- optimize(f=computeR,interval=c(-0.1,0.5),tol=1e-20)
				r[season,sp] <- fit_r$minimum;
			}
			else
				r[season,sp] <- 0;
		}
	}
}



######################################## Estimation of model parameters with maximum likelihood ########################################

r_m  <- numeric();
s_e  <- numeric();

for(sp in 1:N_species)
{
	r_m  <- c(r_m, mean(r[,sp],na.rm=TRUE));
	s_e  <- c(s_e, sd((r[,sp] - r_m[sp]),na.rm=TRUE));
}



######################################## Save model parameters in a .txt file ########################################

to_write <- 1;

if(to_write)
{
	## Create the directories if it does not exixt yet
	main_dir <- getwd();
	sub_dir  <- "Models";
	
	if(!file.exists(paste(main_dir, sub_dir, sep = "/", collapse = "/")))
		dir.create(file.path(main_dir, sub_dir), showWarnings = FALSE);
	
	setwd(file.path(main_dir, sub_dir));
	
	sec_dir <- getwd();
	sub_dir <- "M2";
	
	if(!file.exists(paste(sec_dir, sub_dir, sep = "/", collapse = "/")))
		dir.create(file.path(sec_dir, sub_dir), showWarnings = FALSE);
	
	setwd(main_dir);
	
	param              <- rbind(r_m,K,s_e,s_o);
	colnames(param)    <- sp_names;
	colnames(r)        <- sp_names;
	colnames(bio_real) <- sp_names;
	
	file_names <- c("Models/M2/parameters.txt","Models/M2/r.txt","Models/M2/bio_real.txt");
	files      <- list(round(param,8),round(r,8),round(bio_real,8));
	
	for(i in 1:length(file_names))
	{
		if(file.exists(file_names[i]))
		{
			x <- readline(paste("The file", file_names[i], "already exists. Do you want to overwritte it ? (y/n):"))
		
			if(x == "y")
				write.table(files[[i]], file_names[i], sep="\t", row.names=FALSE, col.names=TRUE);
		}
		else
			write.table(files[[i]], file_names[i], sep="\t", row.names=FALSE, col.names=TRUE);
	}
}
