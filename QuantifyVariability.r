rm(list=ls(all=TRUE))


#######################################################################################################################################
#																																	  #
#		Compute the effects of environmental and demographic variability as well as observation error on the obseved variance of	  #
#		monospecific community biomass.																								  #
#																																	  #
#######################################################################################################################################


## Read species parameters
param     <- read.table("Full/parameters.txt",sep="\t",h=TRUE,stringsAsFactors=FALSE);
param_mat <- matrix(unlist(param), ncol=length(param));

## Simulation details
B_ini     <- 2;
N_seasons <- 1000;
N_days    <- 60;

## Exclude species which have been rejected in the model checking step
sp_names  <- names(param);
names_del <- c("Car.pra", "Ger.pra", "Luz.cam", "Tri.fra");
to_del    <- numeric();

for(i in 1:length(names_del))
	to_del <- c(to_del, which(sp_names == names_del[i]));

N_species <- (1:length(param))[-to_del];



#################################### Estimate environmental, demographic and observation effects ######################################

## Draw random variables in normal distribution with zero mean and unit variance
set.seed(0)
u_e <- rnorm(N_seasons);
set.seed(1)
u_d <- rnorm(N_seasons);
set.seed(2)
u_o <- rnorm(N_seasons);


B_ref  <- matrix(NA,N_seasons,length(param));
B_env  <- matrix(NA,N_seasons,length(param));
B_demo <- matrix(NA,N_seasons,length(param));

V_ref   <- numeric();
V_env   <- numeric();
V_demo  <- numeric();


## Simulate the community dynamics with the complete model and with the shuting down of observation error and/or demographic stochasticity
for(sp in N_species)
{
	for(s in 1:N_seasons)
	{
		B_prec     <- B_ini;
		B_prec_env <- B_ini;
		
		for(d in 1:N_days)
		{
			B_cour_env <- B_prec_env + (param[[sp]][1] + param[[sp]][3]*u_e[s])*B_prec_env*(1 - (B_prec_env/param[[sp]][2]));
			B_cour     <- B_prec + (param[[sp]][1] + param[[sp]][3]*u_e[s] + param[[sp]][4]*u_d[s])*B_prec*(1 - (B_prec/param[[sp]][2]));
			
			B_prec     <- B_cour;
			B_prec_env <- B_cour_env;
		}
		
		B_ref[s,sp]  <- exp(log(B_cour) + param[[sp]][5]*u_o[s]);
		B_env[s,sp]  <- B_cour_env;
		B_demo[s,sp] <- B_cour;
	}
	
	V_ref  <- c(V_ref,var(B_ref[,sp]));
	V_env  <- c(V_env,var(B_env[,sp]));
	V_demo <- c(V_demo,var(B_demo[,sp]));
}

## Compute the effects of environmental and demographic stochasticity as well as observation error
effect_env  <- V_env/V_ref;
effect_demo <- (V_demo - V_env)/V_ref;
effect_obs  <- (V_ref - V_demo)/V_ref;

M_effects   <- c(mean(effect_env),mean(effect_demo),mean(effect_obs));


## Number of species for which environmental, demographic and observation variability respectively dominates
env_dom  <- numeric();
demo_dom <- numeric();
obs_dom  <- numeric();

for(sp in 1:length(effect_env))
{
	if(effect_env[sp] > effect_demo[sp] && effect_env[sp] > effect_obs[sp])
		env_dom <- c(env_dom,sp)
	else if(effect_demo[sp] > effect_env[sp] && effect_demo[sp] > effect_obs[sp])
		demo_dom <- c(demo_dom,sp)
	else
		obs_dom <- c(obs_dom,sp);
}



########################################### Relative importance of effects according to PFT ###########################################

pft           <- read.table("Data/pft.csv",sep=";");
plant_types   <- unique(as.vector(pft[[3]]));
sp_names_jena <- sp_names[-to_del];
param_pft     <- matrix(0,0,5);
domina        <- list();

for(i in 1:length(plant_types))
{
	type_sp   <- which(as.vector(pft[[3]])==plant_types[i]);
	ind_param <- numeric();
	for(j in 1:length(type_sp))
		ind_param <- c(ind_param,which(sp_names_jena==pft[[1]][type_sp[j]]));
	
	names_pft  <- c(names_pft, sp_names[ind_param]) 
	param_pft  <- rbind(param_pft,t(param_mat[,ind_param]))
	domina_tmp <- cbind(effect_env,effect_demo,effect_obs)[ind_param,];
	domina     <- c(domina,list(domina_tmp));
}



######################################## Write the relative importance of effects in .csv file ########################################

to_write <- 1;

if(to_write)
{
	tab_1 <- cbind(round(param_pft,3),round(rbind(domina[[1]],domina[[2]],domina[[3]],domina[[4]]),2));
	write.csv2(tab_1,"Full/Table_1.csv")
}



###################################################### Plot the results (Fig. 3) ######################################################

to_plot <- 0;


if(to_plot)
{
	png("Figures/Variability.png",width=800,height=800);

	par(mfrow=c(2,2))
	par(cex.axis=1.75)
	par(cex.lab=2)
	par(mar=c(5,6,4,3))
	par(mgp=c(3,2,0))

	my_letters <- c("(A)","(B)","(C)","(D)");

	for(ft in 1:length(domina))
	{
		boxplot(domina[[ft]],ylim=c(0,1),names=c(expression(Sigma[e]),expression(Sigma[d]),expression(Sigma[o])));
		mtext(my_letters[ft],side=3,line=1.5,adj=0.0,cex=1.5)
	}

	dev.off();
}
