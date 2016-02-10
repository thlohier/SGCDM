rm(list=ls(all=TRUE))


#######################################################################################################################################
#																																	  #
#		Compute three different indicators enabling to assess the model prediction accuracy regarding biomass average and standard	  #
#		deviation as well as Simpson concentration index for mixtures varying in their diversity level.								  #
#																																	  #
#######################################################################################################################################


N_species   <- c(2,4,8,16);
N_replica   <- 100;

plotcode    <- numeric();
N_mixtures  <- c(0);

## Correlation coefficients for each diversity level
R     <- matrix(0,3,(length(N_species)+1));

## Proportion of mixtures significantly deviating from model predictions for diversity level (N_dev)
N_dev <- matrix(0,3,(length(N_species)+1));

## Overall deviation significance from predictions for each diversity level(p_dev)
p_dev <- matrix(0,3,(length(N_species)+1));


rho_mean <- matrix(NA,N_replica,length(N_species));
rho_sd   <- matrix(NA,N_replica,length(N_species));
rho_l    <- matrix(NA,N_replica,length(N_species));

M_sim_vect  <- numeric();
SD_sim_vect <- numeric();
L_sim_vect  <- numeric();

M_obs_vect  <- numeric();
SD_obs_vect <- numeric();
L_obs_vect  <- numeric();

Mean_rnd    <- matrix(NA,N_replica,0);
Sd_rnd      <- matrix(NA,N_replica,0);
L_rnd       <- matrix(NA,N_replica,0);


for(i in 1:length(N_species))
{
	## Observed community statistics
	fic_names    <- paste("Full/Mixtures/stat_obs_",N_species[i],".txt",sep="");
	stat_obs_tmp <- read.table(fic_names,h=TRUE);
	stat_obs     <- matrix(unlist(stat_obs_tmp),ncol=length(stat_obs_tmp));

	## Simulated community statistics
	fic_names    <- paste("Full/Mixtures/stat_sim_",N_species[i],".txt",sep="");
	stat_sim_tmp <- read.table(fic_names,h=TRUE);
	stat_sim     <- matrix(unlist(stat_sim_tmp),ncol=length(stat_sim_tmp));

	## Simulated average community biomass
	fic_names <- paste("Full/Mixtures/mean_rep_",N_species[i],".txt",sep="");
	Mean_tmp  <- read.table(fic_names,h=TRUE);
	Mean_rep  <- matrix(unlist(Mean_tmp),ncol=length(Mean_tmp));
	mix_codes <- names(Mean_tmp);
	
	## Simulated standard deviation of community biomass
	fic_names <- paste("Full/Mixtures/deviation_rep_",N_species[i],".txt",sep="");
	Sd_tmp    <- read.table(fic_names,h=TRUE);
	Sd_rep    <- matrix(unlist(Sd_tmp),ncol=length(Sd_tmp));
	
	## Simulated  Simpson's concentration index  of community biomass
	fic_names <- paste("Full/Mixtures/simpson_rep_",N_species[i],".txt",sep="");
	L_tmp     <- read.table(fic_names,h=TRUE);
	L_rep     <- matrix(unlist(L_tmp),ncol=length(L_tmp));
	
	
	## Correlation coefficient between observed and predicted values (R^2)
	R[1,(i+1)] <- round(cor(stat_obs[1,],stat_sim[1,])^2,2);
	R[2,(i+1)] <- round(cor(stat_obs[2,],stat_sim[2,])^2,2);
	R[3,(i+1)] <- round(cor(stat_obs[3,],stat_sim[3,])^2,2);
	
	
	## Proportion of mixtures significantly deviating from model predictions (N_dev)
	plotcode <- c(plotcode,mix_codes);
	
	p_val_m  <- numeric();
	p_val_sd <- numeric();
	p_val_l  <- numeric();

	for(m in 1:length(mix_codes))
	{
		p_val_m  <- c(p_val_m,(length(which(Mean_rep[,m] > stat_obs[1,m]))/length(Mean_rep[,m])));
		p_val_sd <- c(p_val_sd,(length(which(Sd_rep[,m] > stat_obs[2,m]))/length(Sd_rep[,m])));
		p_val_l  <- c(p_val_l,(length(which(L_rep[,m] > stat_obs[3,m]))/length(L_rep[,m])));
	}
	
	N_dev[1,(i+1)] <- round(length(which(p_val_m < 0.025 | p_val_m > 0.975)),2);
	N_dev[2,(i+1)] <- round(length(which(p_val_sd < 0.025 | p_val_sd > 0.975)),2);
	N_dev[3,(i+1)] <- round(length(which(p_val_l < 0.025 | p_val_l > 0.975)),2);
	
	
	## Overall deviation significance from predictions (p_dev)
	Mean_rnd_tmp <- matrix(NA,N_replica,length(mix_codes));
	Sd_rnd_tmp   <- matrix(NA,N_replica,length(mix_codes));
	L_rnd_tmp    <- matrix(NA,N_replica,length(mix_codes));
	
	rho_mean_tmp <- numeric();
	rho_sd_tmp   <- numeric();
	rho_l_tmp    <- numeric();
	
	for(r in 1:N_replica)
	{
		set.seed(r);
		to_draw <- sample(1:dim(Mean_rep)[1],length(mix_codes));
		
		for(k in 1:length(to_draw))
		{
			Mean_rnd_tmp[r,k] <- Mean_rep[to_draw[k],k];
			Sd_rnd_tmp[r,k]   <- Sd_rep[to_draw[k],k];
			L_rnd_tmp[r,k]    <- L_rep[to_draw[k],k];
		}
		
		rho_mean_tmp <- c(rho_mean_tmp,cor(Mean_rnd_tmp[r,],stat_sim[1,])^2);
		rho_sd_tmp   <- c(rho_sd_tmp,cor(Sd_rnd_tmp[r,],stat_sim[2,])^2);
		rho_l_tmp    <- c(rho_l_tmp,cor(L_rnd_tmp[r,],stat_sim[3,])^2);
	}

	Mean_rnd <- cbind(Mean_rnd,Mean_rnd_tmp);
	Sd_rnd   <- cbind(Sd_rnd,Sd_rnd_tmp);
	L_rnd    <- cbind(L_rnd,L_rnd_tmp);

	rho_mean[,i] <- rho_mean_tmp;
	rho_sd[,i]   <- rho_sd_tmp;
	rho_l[,i]    <- rho_l_tmp;

	p_dev[1,(i+1)] <- min(length(which(rho_mean[,i] < R[1,(i+1)])), length(which(rho_mean[,i] > R[1,(i+1)])))/100;
	p_dev[2,(i+1)] <- min(length(which(rho_sd[,i] < R[2,(i+1)])), length(which(rho_sd[,i] > R[2,(i+1)])))/100;
	p_dev[3,(i+1)] <- min(length(which(rho_l[,i] < R[3,(i+1)])), length(which(rho_l[,i] > R[3,(i+1)])))/100;
	
	
	## Number of mixtures for each diversity level
	N_mixtures <- c(N_mixtures, length(mix_codes));
		
	
	## Data for computing previous indicators across mixtures
	M_obs_vect  <- c(M_obs_vect,stat_obs[1,]);
	SD_obs_vect <- c(SD_obs_vect,stat_obs[2,]);
	L_obs_vect  <- c(L_obs_vect,stat_obs[3,]);
	
	M_sim_vect  <- c(M_sim_vect,stat_sim[1,]);
	SD_sim_vect <- c(SD_sim_vect,stat_sim[2,]);
	L_sim_vect  <- c(L_sim_vect,stat_sim[3,]);
}

## Correlation coefficient between observed and predicted values (R^2)
R[1,1]  <- round(cor(M_obs_vect,M_sim_vect)^2,2);
R[2,1]  <- round(cor(SD_obs_vect,SD_sim_vect)^2,2);
R[3,1]  <- round(cor(L_obs_vect,L_sim_vect)^2,2);


## Proportion of mixtures significantly deviating from model predictions (N_dev)
N_dev[1,1]  <- sum(N_dev[1,])/sum(N_mixtures);
N_dev[2,1]  <- sum(N_dev[2,])/sum(N_mixtures);
N_dev[3,1]  <- sum(N_dev[3,])/sum(N_mixtures);


## Overall deviation significance from predictions (p_dev)
rho_mean_tot <- numeric();
rho_sd_tot   <- numeric();
rho_l_tot    <- numeric();

for(r in 1:N_replica)
{
	rho_mean_tot <- c(rho_mean_tot,cor(Mean_rnd[r,],M_sim_vect)^2);
	rho_sd_tot   <- c(rho_sd_tot,cor(Sd_rnd[r,],SD_sim_vect)^2);
	rho_l_tot    <- c(rho_l_tot,cor(L_rnd[r,],L_sim_vect)^2);
}

p_dev[1,1]  <- min(length(which(rho_mean_tot < R[1,1])), length(which(rho_mean_tot > R[1,1])))/100;
p_dev[2,1]  <- min(length(which(rho_sd_tot < R[2,1])), length(which(rho_sd_tot > R[2,1])))/100;
p_dev[3,1]  <- min(length(which(rho_l_tot < R[3,1])), length(which(rho_l_tot > R[3,1])))/100;



#################################### Write the indicators used to assess model accuracy in .csv file ####################################

to_write <- 1;

if(to_write)
{
	tab_3 <- matrix(NA,0,(length(N_species)+1));
	
	for(i in 1:dim(R)[1])
		tab_3 <- rbind(tab_3,rbind(R[i,],N_dev[i,],p_dev[i,]));
		
	write.csv2(round(tab_3,2),"Tables/Table_3.csv");
}


###################################################### Plot the results (Fig. 4) ######################################################

to_plot <- 1;

if(to_plot)
{
	stat_obs_vect <- rbind(M_obs_vect, SD_obs_vect, L_obs_vect);
	stat_sim_vect <- rbind(M_sim_vect, SD_sim_vect, L_sim_vect);
	
	fig_nb <- c("(A)","(B)","(C)");
	x_lab  <- c(expression("log(Observed "*bar(B)*")"),expression("log(Observed "*sigma[B]*")"),expression("Observed "*lambda));
	y_lab  <- c(expression("log(Predicted "*bar(B)*")"),expression("log(Predicted "*sigma[B]*")"),expression("Predicted "*lambda));
	my_col <- c("black","red","blue","green");

	
	png("Figures/predictions.png",width=900,height=300)
	par(mfrow=c(1,3))
	par(cex.axis=1.7)
	par(cex.lab=1.9)
	par(mar=c(5,6,4,3))
	
	for(i in 1:dim(stat_obs_vect)[1])
	{
		y_min <- log(min(c(stat_obs_vect[i,],stat_sim_vect[i,])));
		y_max <- log(max(c(stat_obs_vect[i,],stat_sim_vect[i,])));
		
		plot(log(stat_obs_vect[i,]),log(stat_obs_vect[i,]),ylim=c(y_min,y_max),xlab=x_lab[i],ylab=y_lab[i],type="l");

		for(k in 1:length(N_species))
			lines(log(stat_obs_vect[i,(sum(N_mixtures[1:k])+1):(sum(N_mixtures[1:(k+1)]))]),log(stat_sim_vect[i,(sum(N_mixtures[1:k])+1):(sum(N_mixtures[1:(k+1)]))]),type="p",pch=20+k,col=my_col[k],cex=2,lwd=2);
		
		my_legend <- bquote(R^2 == .(R[i,1]));
		legend("topleft",as.expression(my_legend),bty="n",cex=1.75);
		mtext(fig_nb[i],side=3,line=1.5,adj=0.0,cex=1.5);
	}
	
	dev.off();
}
