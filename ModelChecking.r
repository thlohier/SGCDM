rm(list=ls(all=TRUE))


#######################################################################################################################################
#																																	  #
#		Model checking including:																									  #
#				1. Normality test (Kolmogorov-Smirnov) for the environmental, demographic and observation variables					  #
#				2. Test of the model ability to reproduce monospecific community dynamics											  #
#																																	  #
#######################################################################################################################################


param_tmp   <- read.table("Full/parameters.txt",h=TRUE);
param       <- matrix(unlist(param_tmp),ncol=length(param_tmp));

u_e_tmp     <- read.table("Full/env.txt",h=TRUE);
u_e         <- matrix(unlist(u_e_tmp),ncol=length(u_e_tmp));

N_species <- length(param_tmp);
sp_names  <- names(param_tmp);


########################################## Normality test for environmental variables u_e(T) ##########################################

kolmo_env  <- numeric();

for(sp in 1:N_species)
{	
	u_e_sp    <- u_e[which(!is.na(u_e[,sp])),sp];
	
	if(length(u_e_sp) > 3)
	{
		kolmo_tmp <- try(ks.test(u_e_sp,"pnorm",0,1),silent=TRUE);
		
		if(class(kolmo_tmp)!="try-error")
		kolmo_env <- c(kolmo_env,kolmo_tmp$p.value)
		else
			kolmo_env <- c(kolmo_env,NA);
	}
	else
		kolmo_env <- c(kolmo_env,NA);
}

kolmo_env_plot <- kolmo_env[which(!is.na(kolmo_env))];



########################################## Normality test for demographic variables u_d(T,p) ##########################################

kolmo_demo <- numeric();

for(sp in 1:N_species)
{	
	fic_name <- paste("Full/Demo/demo_",sp_names[sp],".txt",sep="");
	u_d_tmp  <- read.table(fic_name);
	u_d      <- c(matrix(unlist(u_d_tmp),ncol=length(u_d_tmp)));
	
	if(length(u_d[which(!is.na(u_d))]) > 3)
	{
		kolmo_tmp <- try(ks.test(u_d[which(!is.na(u_d))],"pnorm",0,1),silent=TRUE);
		
		if(class(kolmo_tmp)!="try-error")
			kolmo_demo <- c(kolmo_demo,kolmo_tmp$p.value)
		else
			kolmo_demo <- c(kolmo_demo,NA);
	}
	else
		kolmo_demo <- c(kolmo_demo,NA);
}

kolmo_demo_plot <- kolmo_demo[which(!is.na(kolmo_demo))];



######################################### Normality test for observation variables u_o(T,p,s) #########################################

kolmo_obs <- numeric();

for(sp in 1:N_species)
{	
	fic_name <- paste("Full/Obs/obs_",sp_names[sp],".txt",sep="");
	u_o_tmp  <- read.table(fic_name);
	u_o      <- c(matrix(unlist(u_o_tmp),ncol=length(u_o_tmp)));
				
	if(length(u_o[which(!is.na(u_o))]) > 3)
	{
		kolmo_tmp <- try(ks.test(u_o[which(!is.na(u_o))],"pnorm",0,1),silent=TRUE);
		
		if(class(kolmo_tmp)!="try-error")
			kolmo_obs <- c(kolmo_obs,kolmo_tmp$p.value)
		else
			kolmo_obs <- c(kolmo_obs,NA);
	}
	else
		kolmo_obs <- c(kolmo_obs,NA);
}

kolmo_obs_plot <- kolmo_obs[which(!is.na(kolmo_obs))];
	


############################# Ability of the model to reproduce observed statistics with random variables #############################

N_plot     <- 2;
N_sample   <- 2;
N_season   <- 17;
N_days     <- 60;
N_replica  <- 100;

B_obs_m  <- numeric();
B_obs_sd <- numeric();

p_val_m  <- numeric();
p_val_sd <- numeric();

for(sp in 1:N_species)
{
	fic_name  <- paste("Full/Biomass/bio_obs_",sp_names[sp],".txt",sep="");
	B_obs_tmp <- read.table(fic_name);
	B_obs     <- array(unlist(B_obs_tmp),c(N_season,N_sample,N_plot));
	
	to_na     <- which(is.na(B_obs),arr.ind=T);
	
	B_obs_m   <- c(B_obs_m, mean(B_obs,na.rm=TRUE));
	B_obs_sd  <- c(B_obs_sd, sd(B_obs,na.rm=TRUE));
	
	B_sim_m   <- numeric();
	B_sim_sd  <- numeric();

	for(r in 1:N_replica)
	{
		## Environmental variables u_e(T)
		set.seed(r);
		u_e  <- rnorm(N_season);
		
		## Demographic variables u_d(T,p)
		set.seed(N_replica + r);
		u_d  <- matrix(rnorm(N_season*N_plot,0,1),N_season,N_plot);
		
		## Observation variables u_o(T,p,s)
		set.seed(2*N_replica + r);
		u_o  <- array(rnorm(N_season*N_plot*N_sample,0,1),c(N_season,N_sample,N_plot));
		
		## Species biomass
		B_sp <- array(NA,c(N_season,N_sample,N_plot));
		
		for(p in 1:N_plot)
		{
			for(samp in 1:N_sample)
			{
				for(s in 1:N_season)
				{
					B_prec <- 2;
					
					for(day in 1:N_days)
					{
						B_cour <- B_prec + (param[1,sp] + param[3,sp]*u_e[s] + param[4,sp]*u_d[s,p])*B_prec*(1-(B_prec/param[2,sp]));
						B_prec <- B_cour;
					}
					
					B_sp[s,samp,p] <- exp(log(B_cour) + param[5,sp]*u_o[s,samp,p]);
				}
			}
		}
		
		B_sp[to_na] <- NA;
			
		B_sim_m  <- c(B_sim_m,mean(B_sp,na.rm=TRUE));
		B_sim_sd <- c(B_sim_sd,sd(B_sp,na.rm=TRUE));
	}
	
	p_val_m  <- c(p_val_m,length(which(B_sim_m > B_obs_m[sp]))/length(B_sim_m));
	p_val_sd <- c(p_val_sd,length(which(B_sim_sd > B_obs_sd[sp]))/length(B_sim_sd));
}



######################################### Write the results of model checking in a .csv fils ##########################################

to_write <- 1;

if(to_write)
{
	## Create the directory if it does not exist yet
	main_dir        <- getwd();
	sub_dir         <- "Tables";
	
	if(!file.exists(paste(main_dir, sub_dir, sep = "/", collapse = "/")))
		dir.create(file.path(main_dir, sub_dir), showWarnings = FALSE);
	
	
	## Write Tab. S4 in a .csv file
	tab_S4           <- cbind(round(cbind(kolmo_env,kolmo_demo,kolmo_obs),3),round(cbind(p_val_m,p_val_sd),2));
	rownames(tab_S4) <- sp_names;
	colnames(tab_S4) <- c("u_e","u_d","u_o","B_av","B_sd");
	
	if(file.exists("Tables/Table_S4.csv"))
	{
		x <- readline("The file Table_S4.csv already exists. Do you want to overwritte it ? (y/n):")
		
		if(x == "y")
			write.csv2(tab_S4, "Tables/Table_S4.csv", row.names=TRUE);
	}
	else
		write.csv2(tab_S4, "Tables/Table_S4.csv", row.names=TRUE);
}



###################################################### Plot the results (Fig. 2) ######################################################

to_plot <- 1;

if(to_plot)
{
	## Create the directory if it does not exist yet
	main_dir        <- getwd();
	sub_dir         <- "Figures";
	
	if(!file.exists(paste(main_dir, sub_dir, sep = "/", collapse = "/")))
		dir.create(file.path(main_dir, sub_dir), showWarnings = FALSE);
	
	
	## Function for the drawing of Fig. 2
	plotFig2 <- function()
	{
		png("Figures/Figure_2.png",width=900,height=600);

		par(mfrow=c(2,3))
		par(cex.axis=1.75)
		par(cex.lab=2)
		par(mar=c(5,6,4,3))

		br <- seq(0,round(max(kolmo_env_plot)+0.05,1),0.05);
		hist(kolmo_env_plot,breaks=br,xlab="p-values",main="");
		lines(c(0.05,0.05),c(-5,20),col="blue",lwd=2);
		mtext("(A)",side=3,line=1.5,adj=0.0,cex=1.5)

		br <- seq(0,round(max(kolmo_demo_plot)+0.05,1),0.05);
		hist(kolmo_demo_plot,breaks=br,xlab="p-values",main="");
		lines(c(0.05,0.05),c(-5,20),col="blue",lwd=2);
		mtext("(B)",side=3,line=1.5,adj=0.0,cex=1.5)

		br <- seq(0,round(max(kolmo_obs_plot)+0.05,1),0.05);
		hist(kolmo_obs_plot,breaks=br,xlab="p-values",main="");
		lines(c(0.05,0.05),c(-5,20),col="blue",lwd=2);
		mtext("(C)",side=3,line=1.5,adj=0.0,cex=1.5)

		br <- seq(0,round(max(p_val_m),1),0.025);
		hist(p_val_m,breaks=br,xlab="p-values",main="");
		lines(c(0.025,0.025),c(-5,20),col="blue",lwd=2);
		lines(c(0.975,0.975),c(-5,20),col="blue",lwd=2);
		mtext("(D)",side=3,line=1.5,adj=0.0,cex=1.5)

		br <- seq(0,round(max(p_val_sd),1),0.025);
		hist(p_val_sd,breaks=br,xlab="p-values",main="");
		lines(c(0.025,0.025),c(-5,20),col="blue",lwd=2);
		lines(c(0.975,0.975),c(-5,20),col="blue",lwd=2);
		mtext("(E)",side=3,line=1.5,adj=0.0,cex=1.5)

		dev.off();
	}
	
	
	## Check if Fig. 3 has already been drawn
	if(file.exists("Figures/Figure_2.png"))
	{
		x <- readline("The file Figure_2.png already exists. Do you want to overwritte it ? (y/n):")
		
		if(x == "y")
			plotFig2();
	}
	else
		plotFig2();	
}
