rm(list=ls(all=TRUE))


#######################################################################################################################################
#																																	  #
#									Perform the analysis describe in Lohier et al. (2015).											  #
#																																	  #
#######################################################################################################################################


## 2.3. Inference based on monoculture data
source("EstimateParameters.r");


## 2.4. Model checking
source("ModelChecking.r");


## 2.5. Assessing the need for the current level of model complexity
source("EstimateParameters_M0.r");
source("EstimateParameters_M1.r");
source("EstimateParameters_M2.r");

source("CompareModels.r");


## 2.6. Quantifying the sources of biomass variability
source("QuantifyVariability.r");


## 2.7. Analyzing species demographical characteristics
source("ComputeCorrelations.r");


## 2.8. Predicting multi-specific community dynamics based on the assumption of competitive symmetry
source("SimulateMixtureDynamics.r");

source("ModelPredictions.r");
