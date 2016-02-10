rm(list=ls(all=TRUE))


#######################################################################################################################################
#																																	  #
#								Compute the minimal model (M0) parameters (r_m, K, sigma_o) 										  #
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

N_species <- length(biomass);

## Compute the average biomass across seasons, plots, samples and species (B_real)
bio_tmp <- numeric();

for(sp in 1:N_species)
{
	for(p in 1:length(biomass[[sp]]))	
		bio_tmp <- c(bio_tmp, c(biomass[[sp]][[p]]));
}

to_keep <- which(bio_tmp!=0 & !is.na(bio_tmp));
B_real  <- exp(mean(log(bio_tmp[to_keep])));


## Compute the observation error
e_obs <- log(bio_tmp[to_keep]) - log(B_real);
s_o   <- sd(e_obs, na.rm=TRUE);



######################################## Computation of K as maximum biomass through plot and year ########################################

K         <- 2*B_real;



######################################## Computation of r(T,p) with a gradient descent ########################################

N_days  <- 60;
r       <- list();
r_m     <- numeric();

computeR <- function(inR)
{	
	B_prec <- B_i;
	
	for(day in 1:N_days)
	{
		B_cour <- B_prec + inR*B_prec*(1-(B_prec/K));
		B_prec <- B_cour;
	}
	
	R <- abs(B_emp - B_cour);
}


B_emp <- B_real;		
fit_r <- optimize(f=computeR,interval=c(-0.1,0.5),tol=1e-20);
r_m   <- fit_r$minimum;



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
	sub_dir <- "M0";
	
	if(!file.exists(paste(sec_dir, sub_dir, sep = "/", collapse = "/")))
		dir.create(file.path(sec_dir, sub_dir), showWarnings = FALSE);
	
	setwd(main_dir);
	
	
	## Write model parameters in a .txt file
	param           <- rbind(rep(B_real),rep(r_m,N_species),rep(K,N_species),rep(s_o,N_species));
	colnames(param) <- sp_names;
	
	if(file.exists("Models/M0/parameters.txt"))
	{
		x <- readline("The file parameters.txt already exists. Do you want to overwritte it ? (y/n):")
		
		if(x == "y")
			write.table(round(param,8), "Models/M0/parameters.txt", sep="\t", row.names=FALSE, col.names=TRUE);
	}
	else
		write.table(round(param,8), "Models/M0/parameters.txt", sep="\t", row.names=FALSE, col.names=TRUE);
}
