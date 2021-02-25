'''
The purpose of this script is to see if regression inferences change based on regression on different subsets of the scale-free dataset.  I have found in its parent script that repression_rate can sometimes have heterogenous effects, when low enough and sample is small enough.

Forked from: /Users/Zack/Documents/UCLA/Research/ProtestTheory/Scripts/Simulations_Regression_v11.R.
Fork date: 07.29.2019.

/Users/Zack/Documents/UCLA/Research/ProtestTheory/Scripts/Simulations_Regression_v11.R is cluttered with other stuff, so I need to start with a fresh slate.  Changes made before starting:
	- I keep everything up to the data section.
	- Delete everything except for work having to do with exp7, which was the scale-free network name in earlier versions of the dataset.

First version.
	1. Update runRegressions to be fixed sample size, not percentage
	2. Get data in form for plotting

v2. 1. Update functions
	2. Results with scaling=2.3
v3. 1. Add facet_wrap by sample size
v4. 1. Lots of subsets
	2. Show main results, xlim(0,.05)
v5. 1. Add code to make figure comparing SW, SF, and HK graphs. - DONE
	2. Delete old functions: makePlot, makeHist, pipeline - DONE
	3. Delete old work - DONE
	4. Debug just below "Keep any start that is equal or below median size and there is change of size" - DONE, ADDED AN IF STATEMENT
	5. Make version of combined graph shortened to 0.05-.15 on old data while wait for UvAmsterdam new data to finish
	6. Update the 3way graph types with scale_color_discrete(limits=c('Watts-Strogatz', 'Holme-Kim', 'Barabasi-Albert')), which makes order of legend match order of results in chart.
v6. Finally have W-S, H-K with repression [0-.15]
	1. Read in new data
	2. Subset scale-free, repeat results
	3. Remove testing_ from file names now that have better small world, holme-kim models
v7. 1. x-axis, "Maximum REpression Rate" -> "Repression Rate"
	2. Update model m, scaling_parameter -> networkStructureControl
v8. 1. Make additional figures that look at other regression results
		1. Update makeResults
		2. Update processResults
v9. 1. Change 'Scale-Free' -> 'Watts-Strogatz'
		NB: 1 a.m. on 11.08.2019, need to rerun some code to get B-A as line label, requires some loading and processing.  Need to do for varying network structure control, others are done.
		NB: 1 a.m. on 11.08.2019, need to do some processResults.  Need to see results when restrict to size change, comparing all 3.  Still have to do at 1:20 a.m., when boarding.  See line 625.
v10. 1. Create chart where initial_size less than median, for all 3 at once.
v11. 1. R1 from PNAS rejection wants robustness check on parameters.
		- Only repression_rate
v12. Version for Network Science R&R.
v13. 1. Update facet sample size
	 2. Add vertical gray
v14. 1. Updating figures for SM with new aesthetics, fewer facets
v15. 1. Update makeResults so that it writes its results to file, is useful to have the raw samples and not just aggregated of it.
	 2. Update processResults w/ new name ot keep separate from other file\'s processResults

'''
##########################
##
##	GLOBALS
##
##########################
set.seed(480193)

library(dplyr)  # For sample_n
library(reshape)
library(ggplot2)

setwd('/Users/Zack/Documents/UCLA/Research/ProtestTheory')


# Runs regression, returns t-statistic of each coefficient on final size
runRegressions <- function(data, model, sampleSize, family, scaling=FALSE, repression=TRUE){
	#keep <- round(nrow(data)*percentage)
	sub <- sample_n(data, sampleSize, replace=TRUE)

	meanRepression <- mean(sub$repression_rate)
	if(family != 'negbinom'){
		reg <- summary(lm(model, sub))$coefficients
		temp <- reg[,3]  # t-statistics
		temp2 <- reg[,2] # standard error
		temp3 <- reg[,1] # coefficients
	}

	# The below if statement is because for one particular iteration, there is sometimes not a value for initial_neighborhood_clustering.  I will therefore add just to keep from erroring out.
	if(length(grep('initial_neighborhood_clustering', names(temp))) == 0){ # If there is no neighborhood clustering variable in temp, do the below.  There will also not be one in temp2 or temp3
		temp[['initial_neighborhood_clustering']] <- NA
		temp2[['initial_neighborhood_clustering']] <- NA
		temp3[['initial_neighborhood_clustering']] <- NA

		}

	# The below if statement is because for one particular iteration, there is sometimes not a value for initial_median_degree.  I will therefore add just to keep from erroring out.
	if(length(grep('initial_median_degree', names(temp))) == 0){ # If there is no neighborhood clustering variable in temp, do the below.  There will also not be one in temp2 or temp3
		temp[['initial_median_degree']] <- NA
		temp2[['initial_median_degree']] <- NA
		temp3[['initial_median_degree']] <- NA

		}


	# The below if statement is because for one particular iteration, there is sometimes not a value for initial_median_degree.  I will therefore add just to keep from erroring out.
	if(length(grep('initial_size', names(temp))) == 0){ # If there is no neighborhood clustering variable in temp, do the below.  There will also not be one in temp2 or temp3
		temp[['initial_size']] <- NA
		temp2[['initial_size']] <- NA
		temp3[['initial_size']] <- NA

		}

	if(length(grep('initial_density', names(temp))) == 0){ # If there is no neighborhood clustering variable in temp, do the below.  There will also not be one in temp2 or temp3
		temp[['initial_density']] <- NA
		temp2[['initial_density']] <- NA
		temp3[['initial_density']] <- NA

		}

	if(length(grep('initial_nodes_clustering', names(temp))) == 0){ # If there is no neighborhood clustering variable in temp, do the below.  There will also not be one in temp2 or temp3
		temp[['initial_nodes_clustering']] <- NA
		temp2[['initial_nodes_clustering']] <- NA
		temp3[['initial_nodes_clustering']] <- NA

		}


	# Commented out for now because fails to converge, not worth worrying about at this early stage.
	# if(family == 'negbinom'){
	# 	temp <- summary(glm.nb(model, data=sub, start=rep(0, )))$coefficients[,3]  # t-statistics

	# }

	if(scaling==FALSE & repression==FALSE){
		df <- data.frame(initial_size=temp[['initial_size']], initial_density=temp[['initial_density']], initial_global_clustering=temp[['initial_global_clustering']], initial_mean_degree=temp[['initial_mean_degree']], initial_median_degree=temp[['initial_median_degree']], initial_neighborhood_clustering=temp[['initial_neighborhood_clustering']], initial_nodes_clustering=temp[['initial_nodes_clustering']], initial_size=temp[['initial_size']], initial_size_se=temp2[['initial_size']], initial_density_se=temp2[['initial_density']], initial_global_clustering_se=temp2[['initial_global_clustering']], initial_mean_degree_se=temp2[['initial_mean_degree']], initial_median_degree_se=temp2[['initial_median_degree']], initial_neighborhood_clustering_se=temp2[['initial_neighborhood_clustering']], initial_nodes_clustering_se=temp2[['initial_nodes_clustering']], initial_size_se=temp2[['initial_size']],initial_size_c=temp3[['initial_size']], initial_density_c=temp3[['initial_density']], initial_global_clustering_c=temp3[['initial_global_clustering']], initial_mean_degree=temp3[['initial_mean_degree']], initial_median_degree_c=temp3[['initial_median_degree']], initial_neighborhood_clustering_c=temp3[['initial_neighborhood_clustering']], initial_nodes_clustering_c=temp3[['initial_nodes_clustering']], initial_size_c=temp3[['initial_size']], meanRepressionRate=meanRepression, sampleSize=sampleSize)
	}

	if(scaling==TRUE & repression==FALSE){
	df <- data.frame(scaling_parameter=temp[['scaling_parameter']], initial_size=temp[['initial_size']], initial_density=temp[['initial_density']], initial_global_clustering=temp[['initial_global_clustering']], initial_mean_degree=temp[['initial_mean_degree']], initial_median_degree=temp[['initial_median_degree']], initial_neighborhood_clustering=temp[['initial_neighborhood_clustering']], initial_nodes_clustering=temp[['initial_nodes_clustering']], initial_size=temp[['initial_size']], scaling_parameter_se=temp2[['scaling_parameter']], initial_size_se=temp2[['initial_size']], initial_density_se=temp2[['initial_density']], initial_global_clustering_se=temp2[['initial_global_clustering']], initial_mean_degree_se=temp2[['initial_mean_degree']], initial_median_degree_se=temp2[['initial_median_degree']], initial_neighborhood_clustering_se=temp2[['initial_neighborhood_clustering']], initial_nodes_clustering_se=temp2[['initial_nodes_clustering']], initial_size_se=temp2[['initial_size']], scaling_parameter_c=temp3[['scaling_parameter']], initial_size_c=temp3[['initial_size']], initial_density_c=temp3[['initial_density']], initial_global_clustering_c=temp3[['initial_global_clustering']], initial_mean_degree=temp3[['initial_mean_degree']], initial_median_degree_c=temp3[['initial_median_degree']], initial_neighborhood_clustering_c=temp3[['initial_neighborhood_clustering']], initial_nodes_clustering_c=temp3[['initial_nodes_clustering']], initial_size_c=temp3[['initial_size']], meanRepressionRate=meanRepression, sampleSize=sampleSize)
	}

	if(scaling==TRUE & repression==TRUE){
	df <- data.frame(repression_rate=temp[['repression_rate']], scaling_parameter=temp[['scaling_parameter']], initial_size=temp[['initial_size']], initial_density=temp[['initial_density']], initial_global_clustering=temp[['initial_global_clustering']], initial_mean_degree=temp[['initial_mean_degree']], initial_median_degree=temp[['initial_median_degree']], initial_neighborhood_clustering=temp[['initial_neighborhood_clustering']], initial_nodes_clustering=temp[['initial_nodes_clustering']], initial_size=temp[['initial_size']], repression_rate_se=temp2[['repression_rate']], scaling_parameter_se=temp2[['scaling_parameter']], initial_size_se=temp2[['initial_size']], initial_density_se=temp2[['initial_density']], initial_global_clustering_se=temp2[['initial_global_clustering']], initial_mean_degree_se=temp2[['initial_mean_degree']], initial_median_degree_se=temp2[['initial_median_degree']], initial_neighborhood_clustering_se=temp2[['initial_neighborhood_clustering']], initial_nodes_clustering_se=temp2[['initial_nodes_clustering']], initial_size_se=temp2[['initial_size']], repression_rate_c=temp3[['repression_rate']], scaling_parameter_c=temp3[['scaling_parameter']], initial_size_c=temp3[['initial_size']], initial_density_c=temp3[['initial_density']], initial_global_clustering_c=temp3[['initial_global_clustering']], initial_mean_degree=temp3[['initial_mean_degree']], initial_median_degree_c=temp3[['initial_median_degree']], initial_neighborhood_clustering_c=temp3[['initial_neighborhood_clustering']], initial_nodes_clustering_c=temp3[['initial_nodes_clustering']], initial_size_c=temp3[['initial_size']], meanRepressionRate=meanRepression, sampleSize=sampleSize)
	}
	if(scaling==FALSE & repression==TRUE){
	df <- data.frame(repression_rate=temp[['repression_rate']], initial_size=temp[['initial_size']], initial_density=temp[['initial_density']], initial_global_clustering=temp[['initial_global_clustering']], initial_mean_degree=temp[['initial_mean_degree']], initial_median_degree=temp[['initial_median_degree']], initial_neighborhood_clustering=temp[['initial_neighborhood_clustering']], initial_nodes_clustering=temp[['initial_nodes_clustering']], initial_size=temp[['initial_size']], repression_rate_se=temp2[['repression_rate']], initial_size_se=temp2[['initial_size']], initial_density_se=temp2[['initial_density']], initial_global_clustering_se=temp2[['initial_global_clustering']], initial_mean_degree_se=temp2[['initial_mean_degree']], initial_median_degree_se=temp2[['initial_median_degree']], initial_neighborhood_clustering_se=temp2[['initial_neighborhood_clustering']], initial_nodes_clustering_se=temp2[['initial_nodes_clustering']], initial_size_se=temp2[['initial_size']], repression_rate_c=temp3[['repression_rate']], initial_size_c=temp3[['initial_size']], initial_density_c=temp3[['initial_density']], initial_global_clustering_c=temp3[['initial_global_clustering']], initial_mean_degree=temp3[['initial_mean_degree']], initial_median_degree_c=temp3[['initial_median_degree']], initial_neighborhood_clustering_c=temp3[['initial_neighborhood_clustering']], initial_nodes_clustering_c=temp3[['initial_nodes_clustering']], initial_size_c=temp3[['initial_size']], meanRepressionRate=meanRepression, sampleSize=sampleSize)
	}

	return(df)
}




# 01.07.2020.  R1 from PNAS has asked for robustness check of different variable combinations, so I need to figure out how to do this dynamically.
runRegressions2 <- function(data, model, sampleSize, family, scaling=FALSE, repression=TRUE){
	#keep <- round(nrow(data)*percentage)
	sub <- sample_n(data, sampleSize, replace=TRUE)

	meanRepression <- mean(sub$repression_rate)
	if(family != 'negbinom'){
		reg <- summary(lm(model, sub))$coefficients
		temp <- reg[,3]  # t-statistics
		temp2 <- reg[,2] # standard error
		temp3 <- reg[,1] # coefficients

		# Do not want the intercept
		temp <- temp[2:length(temp)]
		temp2 <- temp2[2:length(temp2)]
		temp3 <- temp3[2:length(temp3)]

		# rename columns, temp stays as is
		names(temp2) <- paste(names(temp2), 'se', sep='_')
		names(temp3) <- paste(names(temp3), 'c', sep='_')

		# Make one dataframe. transpose because data.frame makes 1 column, want 1 row
		temp <- t(data.frame(temp))
		temp2 <- t(data.frame(temp2))
		temp3 <- t(data.frame(temp3))

		df <- data.frame(cbind(temp, temp2, temp3), meanRepressionRate=meanRepression, sampleSize=sampleSize)
	}



	# # The below if statement is because for one particular iteration, there is sometimes not a value for initial_neighborhood_clustering.  I will therefore add just to keep from erroring out.
	# if(length(grep('initial_neighborhood_clustering', names(temp))) == 0){ # If there is no neighborhood clustering variable in temp, do the below.  There will also not be one in temp2 or temp3
	# 	temp[['initial_neighborhood_clustering']] <- NA
	# 	temp2[['initial_neighborhood_clustering']] <- NA
	# 	temp3[['initial_neighborhood_clustering']] <- NA

	# 	}

	# # The below if statement is because for one particular iteration, there is sometimes not a value for initial_median_degree.  I will therefore add just to keep from erroring out.
	# if(length(grep('initial_median_degree', names(temp))) == 0){ # If there is no neighborhood clustering variable in temp, do the below.  There will also not be one in temp2 or temp3
	# 	temp[['initial_median_degree']] <- NA
	# 	temp2[['initial_median_degree']] <- NA
	# 	temp3[['initial_median_degree']] <- NA

	# 	}


	# # The below if statement is because for one particular iteration, there is sometimes not a value for initial_median_degree.  I will therefore add just to keep from erroring out.
	# if(length(grep('initial_size', names(temp))) == 0){ # If there is no neighborhood clustering variable in temp, do the below.  There will also not be one in temp2 or temp3
	# 	temp[['initial_size']] <- NA
	# 	temp2[['initial_size']] <- NA
	# 	temp3[['initial_size']] <- NA

	# 	}

	# if(length(grep('initial_density', names(temp))) == 0){ # If there is no neighborhood clustering variable in temp, do the below.  There will also not be one in temp2 or temp3
	# 	temp[['initial_density']] <- NA
	# 	temp2[['initial_density']] <- NA
	# 	temp3[['initial_density']] <- NA

	# 	}

	# if(length(grep('initial_nodes_clustering', names(temp))) == 0){ # If there is no neighborhood clustering variable in temp, do the below.  There will also not be one in temp2 or temp3
	# 	temp[['initial_nodes_clustering']] <- NA
	# 	temp2[['initial_nodes_clustering']] <- NA
	# 	temp3[['initial_nodes_clustering']] <- NA

	# 	}


	# # Commented out for now because fails to converge, not worth worrying about at this early stage.
	# # if(family == 'negbinom'){
	# # 	temp <- summary(glm.nb(model, data=sub, start=rep(0, )))$coefficients[,3]  # t-statistics

	# # }

	# if(scaling==FALSE & repression==FALSE){
	# 	df <- data.frame(initial_size=temp[['initial_size']], initial_density=temp[['initial_density']], initial_global_clustering=temp[['initial_global_clustering']], initial_mean_degree=temp[['initial_mean_degree']], initial_median_degree=temp[['initial_median_degree']], initial_neighborhood_clustering=temp[['initial_neighborhood_clustering']], initial_nodes_clustering=temp[['initial_nodes_clustering']], initial_size=temp[['initial_size']], initial_size_se=temp2[['initial_size']], initial_density_se=temp2[['initial_density']], initial_global_clustering_se=temp2[['initial_global_clustering']], initial_mean_degree_se=temp2[['initial_mean_degree']], initial_median_degree_se=temp2[['initial_median_degree']], initial_neighborhood_clustering_se=temp2[['initial_neighborhood_clustering']], initial_nodes_clustering_se=temp2[['initial_nodes_clustering']], initial_size_se=temp2[['initial_size']],initial_size_c=temp3[['initial_size']], initial_density_c=temp3[['initial_density']], initial_global_clustering_c=temp3[['initial_global_clustering']], initial_mean_degree=temp3[['initial_mean_degree']], initial_median_degree_c=temp3[['initial_median_degree']], initial_neighborhood_clustering_c=temp3[['initial_neighborhood_clustering']], initial_nodes_clustering_c=temp3[['initial_nodes_clustering']], initial_size_c=temp3[['initial_size']], meanRepressionRate=meanRepression, sampleSize=sampleSize)
	# }

	# if(scaling==TRUE & repression==FALSE){
	# df <- data.frame(scaling_parameter=temp[['scaling_parameter']], initial_size=temp[['initial_size']], initial_density=temp[['initial_density']], initial_global_clustering=temp[['initial_global_clustering']], initial_mean_degree=temp[['initial_mean_degree']], initial_median_degree=temp[['initial_median_degree']], initial_neighborhood_clustering=temp[['initial_neighborhood_clustering']], initial_nodes_clustering=temp[['initial_nodes_clustering']], initial_size=temp[['initial_size']], scaling_parameter_se=temp2[['scaling_parameter']], initial_size_se=temp2[['initial_size']], initial_density_se=temp2[['initial_density']], initial_global_clustering_se=temp2[['initial_global_clustering']], initial_mean_degree_se=temp2[['initial_mean_degree']], initial_median_degree_se=temp2[['initial_median_degree']], initial_neighborhood_clustering_se=temp2[['initial_neighborhood_clustering']], initial_nodes_clustering_se=temp2[['initial_nodes_clustering']], initial_size_se=temp2[['initial_size']], scaling_parameter_c=temp3[['scaling_parameter']], initial_size_c=temp3[['initial_size']], initial_density_c=temp3[['initial_density']], initial_global_clustering_c=temp3[['initial_global_clustering']], initial_mean_degree=temp3[['initial_mean_degree']], initial_median_degree_c=temp3[['initial_median_degree']], initial_neighborhood_clustering_c=temp3[['initial_neighborhood_clustering']], initial_nodes_clustering_c=temp3[['initial_nodes_clustering']], initial_size_c=temp3[['initial_size']], meanRepressionRate=meanRepression, sampleSize=sampleSize)
	# }

	# if(scaling==TRUE & repression==TRUE){
	# df <- data.frame(repression_rate=temp[['repression_rate']], scaling_parameter=temp[['scaling_parameter']], initial_size=temp[['initial_size']], initial_density=temp[['initial_density']], initial_global_clustering=temp[['initial_global_clustering']], initial_mean_degree=temp[['initial_mean_degree']], initial_median_degree=temp[['initial_median_degree']], initial_neighborhood_clustering=temp[['initial_neighborhood_clustering']], initial_nodes_clustering=temp[['initial_nodes_clustering']], initial_size=temp[['initial_size']], repression_rate_se=temp2[['repression_rate']], scaling_parameter_se=temp2[['scaling_parameter']], initial_size_se=temp2[['initial_size']], initial_density_se=temp2[['initial_density']], initial_global_clustering_se=temp2[['initial_global_clustering']], initial_mean_degree_se=temp2[['initial_mean_degree']], initial_median_degree_se=temp2[['initial_median_degree']], initial_neighborhood_clustering_se=temp2[['initial_neighborhood_clustering']], initial_nodes_clustering_se=temp2[['initial_nodes_clustering']], initial_size_se=temp2[['initial_size']], repression_rate_c=temp3[['repression_rate']], scaling_parameter_c=temp3[['scaling_parameter']], initial_size_c=temp3[['initial_size']], initial_density_c=temp3[['initial_density']], initial_global_clustering_c=temp3[['initial_global_clustering']], initial_mean_degree=temp3[['initial_mean_degree']], initial_median_degree_c=temp3[['initial_median_degree']], initial_neighborhood_clustering_c=temp3[['initial_neighborhood_clustering']], initial_nodes_clustering_c=temp3[['initial_nodes_clustering']], initial_size_c=temp3[['initial_size']], meanRepressionRate=meanRepression, sampleSize=sampleSize)
	# }
	# if(scaling==FALSE & repression==TRUE){
	# df <- data.frame(repression_rate=temp[['repression_rate']], initial_size=temp[['initial_size']], initial_density=temp[['initial_density']], initial_global_clustering=temp[['initial_global_clustering']], initial_mean_degree=temp[['initial_mean_degree']], initial_median_degree=temp[['initial_median_degree']], initial_neighborhood_clustering=temp[['initial_neighborhood_clustering']], initial_nodes_clustering=temp[['initial_nodes_clustering']], initial_size=temp[['initial_size']], repression_rate_se=temp2[['repression_rate']], initial_size_se=temp2[['initial_size']], initial_density_se=temp2[['initial_density']], initial_global_clustering_se=temp2[['initial_global_clustering']], initial_mean_degree_se=temp2[['initial_mean_degree']], initial_median_degree_se=temp2[['initial_median_degree']], initial_neighborhood_clustering_se=temp2[['initial_neighborhood_clustering']], initial_nodes_clustering_se=temp2[['initial_nodes_clustering']], initial_size_se=temp2[['initial_size']], repression_rate_c=temp3[['repression_rate']], initial_size_c=temp3[['initial_size']], initial_density_c=temp3[['initial_density']], initial_global_clustering_c=temp3[['initial_global_clustering']], initial_mean_degree=temp3[['initial_mean_degree']], initial_median_degree_c=temp3[['initial_median_degree']], initial_neighborhood_clustering_c=temp3[['initial_neighborhood_clustering']], initial_nodes_clustering_c=temp3[['initial_nodes_clustering']], initial_size_c=temp3[['initial_size']], meanRepressionRate=meanRepression, sampleSize=sampleSize)
	# }

	return(df)
}





##########################
##
##
##	DATA
##
##
##########################
# Below is added on 07.29.2019.  Is from /Users/Zack/Documents/UCLA/Research/ProtestTheory/Scripts/PolNet_01_ProcessSimulationData_v2.R
all <- read.csv('Data/NetworkSimulation/Summer2019/processedData/01_experimentsCombined_repressionNarrow_eigenCore.csv', stringsAsFactors=FALSE)
all$normalized <- round(1000*all$active_nodes/all$total_nodes)

# Make names to match model m
all$scaling_parameter <- all$networkStructureControl

# Separate out the networks
sf <- all[all$graph_type=='scale_free_graph',]
sworld <- all[all$graph_type=='watts_strogatz_graph',]
hk <- all[all$graph_type=='powerlaw_cluster_graph',]



##########################
##
##
##	MODELS
##
##  #NB: If logging, make sure to log initial size but keep variable name the same
##########################
m <- log(normalized+1, 10) ~ initial_density + initial_global_clustering + initial_mean_degree + initial_median_degree + initial_neighborhood_clustering + initial_nodes_clustering + initial_size + scaling_parameter + repression_rate

# Where scaling parameter not vary
ma <- log(normalized+1, 10) ~ initial_density + initial_global_clustering + initial_mean_degree + initial_median_degree + initial_neighborhood_clustering + initial_nodes_clustering + initial_size  + repression_rate

# Only repression_rate
mb <- log(normalized+1, 10) ~ repression_rate


# Only repression_rate, initial_size
mc <- log(normalized+1, 10) ~ repression_rate + initial_size


# Only one of clustering, degree
md <- log(normalized+1, 10) ~ initial_density + initial_mean_degree + initial_neighborhood_clustering + initial_size  + repression_rate




################# WHAT IF I DON'T SELECT A MAX REPRESSION RATE BUT INSTEAD SAMPLE ALL REPRESSION?  THEN LOOK AT SIG. SIGN?


##########################
##
##
##	WORK, ONE GRAPH FOR NETWORK TYPE, 3 LINES PER GRAPH
##
##
##########################

# data = data to use
# repressionRate = vector of repressionRates to examine
# sampleSizes = vector of sample sizes to examine
# trials = how many regressions to run
# regressionMOdel = formula for model
# scalingBool = TRUE, FALSE to control results dataframe
# repressionBool = TRUE, FALSE to control results dataframe
makeResults_PolNet_05 <- function(data, repressionRates, sampleSizes, trials=100, regressionModel, scalingBool=TRUE, repressionBool=TRUE, networkType){
	results <- NULL


	for(i in 2:length(repressionRates)){  # Use 2 so that will always have variation of repression rate
		print(paste0('Repression rate is ', repressionRates[i]))
		temp <- data[data$repression_rate <= repressionRates[i],]

		for(j in 1:length(sampleSizes)){
			for(k in 1:trials){
				r <- runRegressions(data=temp, model=regressionModel, sampleSize=sampleSizes[j], family='regular', scaling=scalingBool, repression=repressionBool)
				r$repression_rate_max <- repressionRates[i]
				#r$percSampled <- percUse[j]
				results <- rbind(results, r)
			}
		print(paste0('Just finished sample size of ', sampleSizes[j]))
		}
	}

	fileout <- paste0('Data/NetworkSimulation/Summer2019/processedData/', networkType, '_filteredSimulationData.csv')
	write.csv(results, fileout, row.names=FALSE)
	return(results)
}


makeResults3 <- function(data, repressionRates=NULL, networkSCRates=NULL, sampleSizes, trials=100, regressionModel, scalingBool=TRUE, repressionBool=TRUE, networkType){
	results <- NULL

	if(!is.null(repressionRates)){
		for(i in 2:length(repressionRates)){  # Use 2 so that will always have variation of repression rate
			print(paste0('Repression rate is ', repressionRates[i]))
			temp <- data[data$repression_rate <= repressionRates[i],]

			for(j in 1:length(sampleSizes)){
				for(k in 1:trials){
					r <- runRegressions(data=temp, model=regressionModel, sampleSize=sampleSizes[j], family='regular', scaling=scalingBool, repression=repressionBool)
					r$repression_rate_max <- repressionRates[i]
					#r$percSampled <- percUse[j]
					results <- rbind(results, r)
				}
			print(paste0('Just finished sample size of ', sampleSizes[j]))
			}
		}
	}
	if(!is.null(networkSCRates)){
		for(i in 2:length(networkSCRates)){  # Use 2 so that will always have variation of repression rate
			print(paste0('Network control rate is ', networkSCRates[i]))
			temp <- data[data$networkStructureControl<=networkSCRates[i],]

			for(j in 1:length(sampleSizes)){
				for(k in 1:trials){
					r <- runRegressions(data=temp, model=regressionModel, sampleSize=sampleSizes[j], family='regular', scaling=scalingBool, repression=repressionBool)
					r$scaling_parameter_max <- networkSCRates[i]
					#r$percSampled <- percUse[j]
					results <- rbind(results, r)
				}
			print(paste0('Just finished sample size of ', sampleSizes[j]))
			}
		}
	}
	fileout <- paste0('Data/NetworkSimulation/Summer2019/processedData/', networkType, '_filteredSimulationData.csv')
	write.csv(results, fileout, row.names=FALSE)
	return(results)
}


# NB: on 01.07.2020, change to runRegressions2
makeResults2 <- function(data, repressionRates, sampleSizes, trials=100, regressionModel, scalingBool=TRUE, repressionBool=TRUE, networkType){
	results <- NULL


	for(i in 2:length(repressionRates)){  # Use 2 so that will always have variation of repression rate
		print(paste0('Repression rate is ', repressionRates[i]))
		#temp <- data[data$repression_rate <= repressionRates[i],]
		temp <- subset(data, repression_rate <= repressionRates[i])

		for(j in 1:length(sampleSizes)){
			for(k in 1:trials){
				r <- runRegressions2(data=temp, model=regressionModel, sampleSize=sampleSizes[j], family='regular', scaling=scalingBool, repression=repressionBool)
				r$repression_rate_max <- repressionRates[i]
				#r$percSampled <- percUse[j]
				results <- rbind(results, r)
			}
		print(paste0('Just finished sample size of ', sampleSizes[j]))
		}
	}
		fileout <- paste0('Data/NetworkSimulation/Summer2019/processedData/', networkType, '_filteredSimulationData.csv')
	write.csv(results, fileout, row.names=FALSE)
	return(results)
}


# data = the results from a big loop, all the regressions
# tscore = level of significance wanted, make sure is positive
# networkType = string describing network being analyzed
processResults_PolNet_05 <- function(data, tscore, networkType){
	tscore <- tscore

	data$repressionSigSign <- cut(data$repression_rate, breaks=c(-Inf, -1*tscore, tscore, Inf), labels=c(-1,0,1), right=FALSE, include.lowest=TRUE)

	fileout <- paste0('Data/NetworkSimulation/Summer2019/processedData/', networkType, '_regressionData.csv')
	write.csv(data, fileout, row.names=FALSE)

	temp <- data.frame(data %>% group_by(repression_rate_max, sampleSize) %>% summarize(Neg = sum(repressionSigSign==-1)/n(), Zero = sum(repressionSigSign==0)/n(), Pos = sum(repressionSigSign==1)/n()))

	temp <- melt(temp, id.vars=c('repression_rate_max', 'sampleSize'))
	temp$Significance <- temp$variable

	fileout_agg <- paste0('Data/NetworkSimulation/Summer2019/processedData/', networkType, '_regressionData_toplot.csv')
	write.csv(temp, fileout_agg, row.names=FALSE)

	return(temp)
}


# data = the results from a big loop, all the regressions
# tscore = level of significance wanted, make sure is positive
# networkType = string describing network being analyzed
processResults2 <- function(data, tscore, networkType){
	tscore <- tscore

	data$SigSign <- cut(data$scaling_parameter, breaks=c(-Inf, -1*tscore, tscore, Inf), labels=c(-1,0,1), right=FALSE, include.lowest=TRUE)

	fileout <- paste0('Data/NetworkSimulation/Summer2019/processedData/', networkType, '_scaling_parameter_regressionData.csv')
	write.csv(data, fileout, row.names=FALSE)

	#data$grouping <- data[[variable]]  #ggplot will not read string from variable

	temp <- data.frame(data %>% group_by(scaling_parameter_max, sampleSize) %>% summarize(Neg = sum(SigSign==-1)/n(), Zero = sum(SigSign==0)/n(), Pos = sum(SigSign==1)/n()))

	temp <- melt(temp, id.vars=c('scaling_parameter_max', 'sampleSize'))  # variable is a user supplied string, so will not be same as variable in following line
	temp$Significance <- temp$variable

	fileout_agg <- paste0('Data/NetworkSimulation/Summer2019/processedData/', networkType, '_networkStructureControl_regressionData_toplot.csv')
	write.csv(temp, fileout_agg, row.names=FALSE)

	return(temp)
}



###################################
#
#	DO WORK, ONE NETWORK AT A TIME
#
###################################

###################################
### Regardless of network scaling
###################################
rates <- unique(sf$repression_rate)  # First will be 0
n <- c(100, 250, 500, 1000, 1500, 2000)

results <- makeResults(data=sf, repressionRates=rates, sampleSizes=n, trials=100, regressionModel=m)

toPlot <- processResults(data=results, tscore=1.96, networkType='scaleFree')

## LIMIT XAXIS
toPlot <- read.csv('Data/NetworkSimulation/Summer2019/processedData/scaleFree_regressionData_toplot.csv', stringsAsFactors=FALSE)
#toPlot$facets <- paste0('Sample Size: ', toPlot$sampleSize, '\n', 'Scaling Parameter = 2.3')
ggplot(toPlot, aes(x=repression_rate_max, y=value, linetype=Significance)) + geom_line() + facet_wrap(~sampleSize, ncol=3, dir='h') + theme_classic()  + xlab('Repression Rate') + ylab('Percent of Regressions') + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlim(0,.06)
ggsave('Figures/June2019/signficance_facet_alphaAll_xlim.jpg', plot=last_plot())



# Below is to read if have not run toPlot line above
toPlot <- read.csv('Data/NetworkSimulation/Summer2019/processedData/scaleFree_regressionData_toplot.csv', stringsAsFactors=FALSE)
#toPlot$facets <- paste0('Sample Size: ', toPlot$sampleSize, '\n', 'Scaling Parameter = 2.3')
ggplot(toPlot, aes(x=repression_rate_max, y=value, linetype=Significance)) + geom_line() + facet_wrap(~sampleSize, ncol=3, dir='h') + theme_classic()  + xlab('Repression Rate') + ylab('Percent of Regressions') + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave('Figures/June2019/signficance_facet_alphaAll.jpg', plot=last_plot())



###################################
### Network Scaling of 2.3
###################################

sf2 <- sf[sf$scaling_parameter==2.3,]

rates <- unique(sf2$repression_rate)
n <- c(100, 250, 500, 1000, 1500, 2000)
#n <- c(100)

results <- makeResults(data=sf2, repressionRates=rates, sampleSizes=n, trials=100, regressionModel=m, scalingBool=FALSE, repressionBool=TRUE)

toPlot <- processResults(data=results, tscore=1.96, networkType='scaleFree_alpha23')

# Below is to read if have not run toPlot line above
toPlot <- read.csv('Data/NetworkSimulation/Summer2019/processedData/scaleFree_alpha23_regressionData_toplot.csv', stringsAsFactors=FALSE)

#toPlot$facets <- paste0('Sample Size: ', toPlot$sampleSize, '\n', 'Scaling Parameter = 2.3')
ggplot(toPlot, aes(x=repression_rate_max, y=value, linetype=Significance)) + geom_line() + facet_wrap(~sampleSize, ncol=3, dir='h') + theme_classic()  + xlab('Repression Rate') + ylab('Percent of Regressions') + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave('Figures/June2019/signficance_facet_alpha23.jpg', plot=last_plot())




###################################
### Keep any start that is equal or below median size
###################################
sf2 <- sf[sf$normalized <= median(sf$normalized),]
rates <- unique(sf2$repression_rate)
n <- c(100, 250, 500, 1000, 1500, 2000)
#n <- c(100)

results <- makeResults(data=sf2, repressionRates=rates, sampleSizes=n, trials=100, regressionModel=m, scalingBool=TRUE, repressionBool=TRUE)

toPlot <- processResults(data=results, tscore=1.96, networkType='scaleFree_initialSize_lteMedian')

# Below is to read if have not run toPlot line above

#toPlot$facets <- paste0('Sample Size: ', toPlot$sampleSize, '\n', 'Scaling Parameter = 2.3')
ggplot(toPlot, aes(x=repression_rate_max, y=value, linetype=Significance)) + geom_line() + facet_wrap(~sampleSize, ncol=3, dir='h') + theme_classic()  + xlab('Repression Rate') + ylab('Percent of Regressions') + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave('Figures/June2019/signficance_facet_initialSize_lteMedian.jpg', plot=last_plot())


###################################
### Keep any start that is equal or below median size and there is change of size
###################################
sf2 <- sf[sf$normalized <= median(sf$normalized),]
sf2 <- sf2[sf2$initial_size != sf2$active_nodes,]

rates <- unique(sf2$repression_rate)
n <- c(100, 250, 500, 1000, 1500, 2000)
#n <- c(100)

### NB: WEIRD, WHEN I RUN MANUALLY IT WORKS, I.E. GO UP AND USE LOOP, NOT FUNCTION.  SO ON 08.08.2019, I WAS ABLE TO USE processResults() on the correct results dataframe
results <- makeResults(data=sf2, repressionRates=rates, sampleSizes=n, trials=100, regressionModel=m, scalingBool=TRUE, repressionBool=TRUE)

toPlot <- processResults(data=results, tscore=1.96, networkType='scaleFree_initialSize_lteMedian_sizeChanged')

ggplot(toPlot, aes(x=repression_rate_max, y=value, linetype=Significance)) + geom_line() + facet_wrap(~sampleSize, ncol=3, dir='h') + theme_classic()  + xlab('Repression Rate') + ylab('Percent of Regressions') + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave('Figures/June2019/signficance_facet_initialSize_lteMedian_sizeChanged.jpg', plot=last_plot())





###################################
### Keep any start that is equal or below median size and realistic network scaling
###################################
sf2 <- sf[sf$normalized <= median(sf$normalized),]
sf2 <- sf2[sf2$scaling_parameter == 2.3,]

rates <- unique(sf2$repression_rate)
n <- c(100, 250, 500, 1000, 1500, 2000)
#n <- c(100)

results <- makeResults(data=sf2, repressionRates=rates, sampleSizes=n, trials=100, regressionModel=m, scalingBool=FALSE, repressionBool=TRUE)

toPlot <- processResults(data=results, tscore=1.96, networkType='scaleFree_initialSize_lteMedian_alpha23')

ggplot(toPlot, aes(x=repression_rate_max, y=value, linetype=Significance)) + geom_line() + facet_wrap(~sampleSize, ncol=3, dir='h') + theme_classic()  + xlab('Repression Rate') + ylab('Percent of Regressions') + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave('Figures/June2019/signficance_facet_initialSize_lteMedian_alpha23.jpg', plot=last_plot())


###################################
### Keep any start that is equal or below median size and there is change and network scaling is 2.3
###################################

sf2 <- sf[sf$normalized <= median(sf$normalized),]
sf2 <- sf2[sf2$initial_size != sf2$active_nodes,]
sf2 <- sf2[sf2$scaling_parameter == 2.3,]

rates <- unique(sf2$repression_rate)
n <- c(100, 250, 500, 1000, 1500, 2000)
#n <- c(100)

results <- makeResults(data=sf2, repressionRates=rates, sampleSizes=n, trials=100, regressionModel=m, scalingBool=FALSE, repressionBool=TRUE)

toPlot <- processResults(data=results, tscore=1.96, networkType='scaleFree_initialSize_lteMedian_sizeChanged_alpha23')

# Below is to read if have not run toPlot line above

#toPlot$facets <- paste0('Sample Size: ', toPlot$sampleSize, '\n', 'Scaling Parameter = 2.3')
ggplot(toPlot, aes(x=repression_rate_max, y=value, linetype=Significance)) + geom_line() + facet_wrap(~sampleSize, ncol=3, dir='h') + theme_classic()  + xlab('Repression Rate') + ylab('Percent of Regressions') + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave('Figures/June2019/signficance_facet_initialSize_lteMedian_sizeChanged_alpha23.jpg', plot=last_plot())


###################################
### Keep any start that is above median size
###################################

sf2 <- sf[sf$normalized > median(sf$normalized),]
rates <- unique(sf2$repression_rate)
n <- c(100, 250, 500, 1000, 1500, 2000)
#n <- c(100)

results <- makeResults(data=sf2, repressionRates=rates, sampleSizes=n, trials=100, regressionModel=m, scalingBool=TRUE, repressionBool=TRUE)

toPlot <- processResults(data=results, tscore=1.96, networkType='scaleFree_initialSize_gtMedian')

# Below is to read if have not run toPlot line above

#toPlot$facets <- paste0('Sample Size: ', toPlot$sampleSize, '\n', 'Scaling Parameter = 2.3')
ggplot(toPlot, aes(x=repression_rate_max, y=value, linetype=Significance)) + geom_line() + facet_wrap(~sampleSize, ncol=3, dir='h') + theme_classic()  + xlab('Repression Rate') + ylab('Percent of Regressions') + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave('Figures/June2019/signficance_facet_initialSize_gtMedian.jpg', plot=last_plot())

###################################
### Keep those where size actually changed
###################################

sf2 <- sf[sf$initial_size != d$active_nodes,]
rates <- unique(sf2$repression_rate)
n <- c(100, 250, 500, 1000, 1500, 2000)
#n <- c(100)

results <- makeResults(data=sf2, repressionRates=rates, sampleSizes=n, trials=100, regressionModel=m, scalingBool=TRUE, repressionBool=TRUE)

toPlot <- processResults(data=results, tscore=1.96, networkType='scaleFree_sizeChanged')

# Below is to read if have not run toPlot line above

#toPlot$facets <- paste0('Sample Size: ', toPlot$sampleSize, '\n', 'Scaling Parameter = 2.3')
ggplot(toPlot, aes(x=repression_rate_max, y=value, linetype=Significance)) + geom_line() + facet_wrap(~sampleSize, ncol=3, dir='h') + theme_classic()  + xlab('Repression Rate') + ylab('Percent of Regressions') + theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text=element_text(size=16), axis.title=element_text(size=16), legend.text=element_text(size=16), legend.title = element_text(size=16), strip.text=element_text(size=16), legend.position='bottom', legend.box='vertical')
ggsave('Figures/June2019/signficance_facet3_sizeChanged.jpg', plot=last_plot())




###################################
### Keep those where initial size > 2, those are dyads so very unlikely to lead to anything
###################################
sf2 <- sf[sf$initial_size > 2,]
rates <- unique(sf2$repression_rate)
n <- c(100, 250, 500, 1000, 1500, 2000)
#n <- c(100)

results <- makeResults(data=sf2, repressionRates=rates, sampleSizes=n, trials=100, regressionModel=m, scalingBool=TRUE, repressionBool=TRUE)

toPlot <- processResults(data=results, tscore=1.96, networkType='scaleFree_initialSize_gt2')

# Below is to read if have not run toPlot line above

#toPlot$facets <- paste0('Sample Size: ', toPlot$sampleSize, '\n', 'Scaling Parameter = 2.3')
ggplot(toPlot, aes(x=repression_rate_max, y=value, linetype=Significance)) + geom_line() + facet_wrap(~sampleSize, ncol=3, dir='h') + theme_classic()  + xlab('Repression Rate') + ylab('Percent of Regressions') + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave('Figures/June2019/signficance_facet_initialSize_gt2.jpg', plot=last_plot())


###################################
### Initial size > 2, is change in protest size.
###################################
sf2 <- sf[sf$initial_size > 2,]
sf2 <- sf2[sf2$initial_size != sf2$active_nodes,]
rates <- unique(sf2$repression_rate)
n <- c(100, 250, 500, 1000, 1500, 2000)
#n <- c(100)

results <- makeResults(data=sf2, repressionRates=rates, sampleSizes=n, trials=100, regressionModel=m, scalingBool=TRUE, repressionBool=TRUE)

toPlot <- processResults(data=results, tscore=1.96, networkType='scaleFree_initialSize_gt2_sizeChanged')

# Below is to read if have not run toPlot line above

#toPlot$facets <- paste0('Sample Size: ', toPlot$sampleSize, '\n', 'Scaling Parameter = 2.3')
ggplot(toPlot, aes(x=repression_rate_max, y=value, linetype=Significance)) + geom_line() + facet_wrap(~sampleSize, ncol=3, dir='h') + theme_classic()  + xlab('Repression Rate') + ylab('Percent of Regressions') + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave('Figures/June2019/signficance_facet_initialSize_gt2_sizeChanged.jpg', plot=last_plot())




##########################
##
##
##	WORK, ONE GRAPH WITH 1 LINE PER NETWORK TYPE
##
##
##########################
### Standardize some information
sf$initial_size <- log(sf$initial_size+1, 10)  # Use log since skewed
sf$scaling_parameter <- sf$networkStructureControl  # Change name to match model variable name

sworld$initial_size <- log(sworld$initial_size+1, 10)
sworld$scaling_parameter <- sworld$networkStructureControl

hk$initial_size <- log10(hk$initial_size+1)
hk$scaling_parameter <- hk$networkStructureControl


###################################
### Regardless of network scaling
###################################

# rates <- sort(unique(all$repression_rate))
# rates <- rates[rates <= .15]  # Keep .15 to match what new results will be.  Only necessary on 08.08.2019 while wait for Amsterdam cluster to finish
rates <- unique(all$repression_rate)
  # On 08.08.2019, hardcode because sf has different rates
#n <- c(1000, 1500, 2000)
n <- c(100, 250, 500, 1000, 1500, 2000)
#n <- c(100)

### COMMENT OUT BELOW IF NECESSARY.
### HAVE RUN AS OF 08.08.2019
# make results
results_sf <- makeResults(data=sf, repressionRates=rates, sampleSizes=n, trials=100, regressionModel=m)
results_sworld <- makeResults(data=sworld, repressionRates=rates, sampleSizes=n, trials=100, regressionModel=m)
results_hk <- makeResults(data=hk, repressionRates=rates, sampleSizes=n, trials=100, regressionModel=m)

# get data to plot
toPlot_sf <- processResults(data=results_sf, tscore=1.96, networkType='3wayCompare_scaleFree')
toPlot_sworld <- processResults(data=results_sworld, tscore=1.96, networkType='3wayCompare_smallWorld')
toPlot_hk <- processResults(data=results_hk, tscore=1.96, networkType='3wayCompare_holmeKim')

toPlot_sf$graph_type <- 'Barabasi-Albert'
toPlot_sworld$graph_type <- 'Watts-Strogatz'
toPlot_hk$graph_type <- 'Holme-Kim'

# combine
toPlot <- rbind(toPlot_sf, toPlot_sworld, toPlot_hk)


###
### RUN BELOW IF NEED TO LOAD DATA LATER, LIKE WHEN MAKING NARROW XLIM
###
toPlot_sf <- read.csv('Data/NetworkSimulation/Summer2019/processedData/3wayCompare_scaleFree_regressionData_toplot.csv', stringsAsFactors=FALSE)
toPlot_sworld <- read.csv('Data/NetworkSimulation/Summer2019/processedData/3wayCompare_smallWorld_regressionData_toplot.csv', stringsAsFactors=FALSE)
toPlot_hk <- read.csv('Data/NetworkSimulation/Summer2019/processedData/3wayCompare_holmeKim_regressionData_toplot.csv', stringsAsFactors=FALSE)


# toPlot_sf <- read.csv('Data/NetworkSimulation/Summer2019/processedData/3wayCompare_scaleFree_regressionData.csv', stringsAsFactors=FALSE)
# toPlot_sworld <- read.csv('Data/NetworkSimulation/Summer2019/processedData/3wayCompare_smallWorld_regressionData.csv', stringsAsFactors=FALSE)
# toPlot_hk <- read.csv('Data/NetworkSimulation/Summer2019/processedData/3wayCompare_holmeKim_regressionData.csv', stringsAsFactors=FALSE)


toPlot_sf$graph_type <- 'Barabasi-Albert'
toPlot_sworld$graph_type <- 'Watts-Strogatz'
toPlot_hk$graph_type <- 'Holme-Kim'
# ###
# ### RUN ABOVE IF NEED TO LOAD DATA LATER, LIKE WHEN MAKING NARROW XLIM
# ###

# combine
toPlot <- rbind(toPlot_sf, toPlot_sworld, toPlot_hk)

# smaller size, comment out as necessary
toPlot <- subset(toPlot, sampleSize >= 1000)

### Plot, all
ggplot(toPlot, aes(x=repression_rate_max, y=value, linetype=Significance, color=graph_type)) + geom_rect(aes(xmin=.054, xmax=.1005, ymin=0, ymax=1), fill='grey95', alpha=.05, inherit.aes=FALSE) + geom_line() + facet_wrap(~sampleSize, ncol=3, dir='h') + theme_classic() +  xlab('Repression Rate') + ylab('Percent of Regressions') + theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text=element_text(size=16), axis.title=element_text(size=16), legend.text=element_text(size=16), legend.title = element_text(size=16), strip.text=element_text(size=16), legend.position='bottom', legend.box='veritical') + labs(color='Graph Type') + scale_color_discrete(limits=c('Watts-Strogatz', 'Holme-Kim', 'Barabasi-Albert'))
ggsave('Figures/August2020/compareNetworks_significance_facet3_color.jpg', plot=last_plot())

ggplot(toPlot, aes(x=repression_rate_max, y=value, linetype=Significance, color=graph_type))+ geom_rect(aes(xmin=.054, xmax=.1005, ymin=0, ymax=1), fill='grey95', alpha=.05, inherit.aes=FALSE) + geom_line() + facet_wrap(~sampleSize, ncol=3, dir='h') + theme_classic() +  xlab('Repression Rate') + ylab('Percent of Regressions') + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_color_grey(limits=c('Watts-Strogatz', 'Holme-Kim', 'Barabasi-Albert')) + labs(color='Graph Type')
ggsave('Figures/August2020/compareNetworks_significance_facet3_grey.jpg', plot=last_plot())




###################################
### ROBUSTNESS, 01.07.2020
###	Regardless of network scaling, different parameter combinations. model is mb.  response to pnas R1
###################################

rates <- unique(all$repression_rate)
  # On 08.08.2019, hardcode because sf has different rates
n <- c(100, 250, 500, 1000, 1500, 2000)

# make results
results_sf <- makeResults2(data=sf, repressionRates=rates, sampleSizes=n, trials=100, regressionModel=mb, scalingBool=FALSE)
results_sworld <- makeResults2(data=sworld, repressionRates=rates, sampleSizes=n, trials=100, regressionModel=mb, scalingBool=FALSE)
results_hk <- makeResults2(data=hk, repressionRates=rates, sampleSizes=n, trials=100, regressionModel=mb, scalingBool=FALSE)

# get data to plot
toPlot_sf <- processResults(data=results_sf, tscore=1.96, networkType='3wayCompare_repressionOnly_scaleFree')
toPlot_sworld <- processResults(data=results_sworld, tscore=1.96, networkType='3wayCompare_repressionOnly_smallWorld')
toPlot_hk <- processResults(data=results_hk, tscore=1.96, networkType='3wayCompare_repressionOnly_holmeKim')

toPlot_sf$graph_type <- 'Barabasi-Albert'
toPlot_sworld$graph_type <- 'Watts-Strogatz'
toPlot_hk$graph_type <- 'Holme-Kim'

# combine
toPlot <- rbind(toPlot_sf, toPlot_sworld, toPlot_hk)


toPlot_sf$graph_type <- 'Barabasi-Albert'
toPlot_sworld$graph_type <- 'Watts-Strogatz'
toPlot_hk$graph_type <- 'Holme-Kim'
# ###
# ### RUN ABOVE IF NEED TO LOAD DATA LATER, LIKE WHEN MAKING NARROW XLIM
# ###

# combine
toPlot <- rbind(toPlot_sf, toPlot_sworld, toPlot_hk)

### Plot, all
ggplot(toPlot, aes(x=repression_rate_max, y=value, linetype=Significance, color=graph_type)) + geom_line() + facet_wrap(~sampleSize, ncol=3, dir='h') + theme_classic() +  xlab('Repression Rate') + ylab('Percent of Regressions') + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(color='Graph Type') + scale_color_discrete(limits=c('Watts-Strogatz', 'Holme-Kim', 'Barabasi-Albert'))
ggsave('Figures/June2019/compareNetworks_significance_facet_color_repressionOnly.jpg', plot=last_plot())

ggplot(toPlot, aes(x=repression_rate_max, y=value, linetype=Significance, color=graph_type)) + geom_line() + facet_wrap(~sampleSize, ncol=3, dir='h') + theme_classic() +  xlab('Repression Rate') + ylab('Percent of Regressions') + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_color_grey(limits=c('Watts-Strogatz', 'Holme-Kim', 'Barabasi-Albert')) + labs(color='Graph Type')
ggsave('Figures/June2019/compareNetworks_significance_facet_grey_repressionOnly.jpg', plot=last_plot())



###################################
### ROBUSTNESS, 01.07.2020
###	Regardless of network scaling, different parameter combinations. model is mc.  response to pnas R1
###################################

rates <- unique(all$repression_rate)
  # On 08.08.2019, hardcode because sf has different rates
n <- c(100, 250, 500, 1000, 1500, 2000)

# make results
results_sf <- makeResults2(data=sf, repressionRates=rates, sampleSizes=n, trials=100, regressionModel=mc, scalingBool=FALSE)
results_sworld <- makeResults2(data=sworld, repressionRates=rates, sampleSizes=n, trials=100, regressionModel=mc, scalingBool=FALSE)
results_hk <- makeResults2(data=hk, repressionRates=rates, sampleSizes=n, trials=100, regressionModel=mc, scalingBool=FALSE)

# get data to plot
toPlot_sf <- processResults(data=results_sf, tscore=1.96, networkType='3wayCompare_repressionAndSize_scaleFree')
toPlot_sworld <- processResults(data=results_sworld, tscore=1.96, networkType='3wayCompare_repressionAndSize_smallWorld')
toPlot_hk <- processResults(data=results_hk, tscore=1.96, networkType='3wayCompare_repressionAndSize_holmeKim')

toPlot_sf$graph_type <- 'Barabasi-Albert'
toPlot_sworld$graph_type <- 'Watts-Strogatz'
toPlot_hk$graph_type <- 'Holme-Kim'

# combine
toPlot <- rbind(toPlot_sf, toPlot_sworld, toPlot_hk)


toPlot_sf$graph_type <- 'Barabasi-Albert'
toPlot_sworld$graph_type <- 'Watts-Strogatz'
toPlot_hk$graph_type <- 'Holme-Kim'
# ###
# ### RUN ABOVE IF NEED TO LOAD DATA LATER, LIKE WHEN MAKING NARROW XLIM
# ###

# combine
toPlot <- rbind(toPlot_sf, toPlot_sworld, toPlot_hk)

### Plot, all
ggplot(toPlot, aes(x=repression_rate_max, y=value, linetype=Significance, color=graph_type)) + geom_line() + facet_wrap(~sampleSize, ncol=3, dir='h') + theme_classic() +  xlab('Repression Rate') + ylab('Percent of Regressions') + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(color='Graph Type') + scale_color_discrete(limits=c('Watts-Strogatz', 'Holme-Kim', 'Barabasi-Albert'))
ggsave('Figures/June2019/compareNetworks_significance_facet_color_repressionAndSize.jpg', plot=last_plot())

ggplot(toPlot, aes(x=repression_rate_max, y=value, linetype=Significance, color=graph_type)) + geom_line() + facet_wrap(~sampleSize, ncol=3, dir='h') + theme_classic() +  xlab('Repression Rate') + ylab('Percent of Regressions') + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_color_grey(limits=c('Watts-Strogatz', 'Holme-Kim', 'Barabasi-Albert')) + labs(color='Graph Type')
ggsave('Figures/June2019/compareNetworks_significance_facet_grey_repressionAndSize.jpg', plot=last_plot())




###################################
### ROBUSTNESS, 01.07.2020
###	Regardless of network scaling, different parameter combinations. model is md.  response to pnas R1
###################################

rates <- unique(all$repression_rate)
  # On 08.08.2019, hardcode because sf has different rates
n <- c(100, 250, 500, 1000, 1500, 2000)

# make results
results_sf <- makeResults2(data=sf, repressionRates=rates, sampleSizes=n, trials=100, regressionModel=md, scalingBool=FALSE)
results_sworld <- makeResults2(data=sworld, repressionRates=rates, sampleSizes=n, trials=100, regressionModel=md, scalingBool=FALSE)
results_hk <- makeResults2(data=hk, repressionRates=rates, sampleSizes=n, trials=100, regressionModel=md, scalingBool=FALSE)

# get data to plot
toPlot_sf <- processResults(data=results_sf, tscore=1.96, networkType='3wayCompare_networkReduced_scaleFree')
toPlot_sworld <- processResults(data=results_sworld, tscore=1.96, networkType='3wayCompare_networkReduced_smallWorld')
toPlot_hk <- processResults(data=results_hk, tscore=1.96, networkType='3wayCompare_networkReduced_holmeKim')

toPlot_sf$graph_type <- 'Barabasi-Albert'
toPlot_sworld$graph_type <- 'Watts-Strogatz'
toPlot_hk$graph_type <- 'Holme-Kim'

# combine
toPlot <- rbind(toPlot_sf, toPlot_sworld, toPlot_hk)


toPlot_sf$graph_type <- 'Barabasi-Albert'
toPlot_sworld$graph_type <- 'Watts-Strogatz'
toPlot_hk$graph_type <- 'Holme-Kim'
# ###
# ### RUN ABOVE IF NEED TO LOAD DATA LATER, LIKE WHEN MAKING NARROW XLIM
# ###

# combine
toPlot <- rbind(toPlot_sf, toPlot_sworld, toPlot_hk)

### Plot, all
ggplot(toPlot, aes(x=repression_rate_max, y=value, linetype=Significance, color=graph_type)) + geom_line() + facet_wrap(~sampleSize, ncol=3, dir='h') + theme_classic() +  xlab('Repression Rate') + ylab('Percent of Regressions') + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(color='Graph Type') + scale_color_discrete(limits=c('Watts-Strogatz', 'Holme-Kim', 'Barabasi-Albert'))
ggsave('Figures/June2019/compareNetworks_significance_facet_color_networkReduced.jpg', plot=last_plot())

ggplot(toPlot, aes(x=repression_rate_max, y=value, linetype=Significance, color=graph_type)) + geom_line() + facet_wrap(~sampleSize, ncol=3, dir='h') + theme_classic() +  xlab('Repression Rate') + ylab('Percent of Regressions') + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_color_grey(limits=c('Watts-Strogatz', 'Holme-Kim', 'Barabasi-Albert')) + labs(color='Graph Type')
ggsave('Figures/June2019/compareNetworks_significance_facet_grey_networkReduced.jpg', plot=last_plot())





###################################
### Regardless of network scaling, size changes
### 11.08.2019
###################################
sf2 <- sf[sf$initial_size != sf$active_nodes,]
sworld2 <- sworld[sworld$initial_size != sworld$active_nodes,]
hk2 <- hk[hk$initial_size != hk$active_nodes,]

# rates <- sort(unique(all$repression_rate))
# rates <- rates[rates <= .15]  # Keep .15 to match what new results will be.  Only necessary on 08.08.2019 while wait for Amsterdam cluster to finish
rates <- unique(all$repression_rate)
  # On 08.08.2019, hardcode because sf has different rates
n <- c(100, 250, 500, 1000, 1500, 2000)
#n <- c(100)

# make results
results_sf <- makeResults(data=sf2, repressionRates=rates, sampleSizes=n, trials=100, regressionModel=m)
results_sworld <- makeResults(data=sworld2, repressionRates=rates, sampleSizes=n, trials=100, regressionModel=m)
results_hk <- makeResults(data=hk2, repressionRates=rates, sampleSizes=n, trials=100, regressionModel=m)

# get data to plot
toPlot_sf <- processResults(data=results_sf, tscore=1.96, networkType='3wayCompare_scaleFree_sizeChange')
toPlot_sworld <- processResults(data=results_sworld, tscore=1.96, networkType='3wayCompare_smallWorld_sizeChange')
toPlot_hk <- processResults(data=results_hk, tscore=1.96, networkType='3wayCompare_holmeKim_sizeChange')

###
### RUN BELOW IF NEED TO LOAD DATA LATER, LIKE WHEN MAKING NARROW XLIM
###
toPlot_sf <- read.csv('Data/NetworkSimulation/Summer2019/processedData/3wayCompare_scaleFree_sizeChange_regressionData_toplot.csv', stringsAsFactors=FALSE)
toPlot_sworld <- read.csv('Data/NetworkSimulation/Summer2019/processedData/3wayCompare_smallWorld_sizeChange_regressionData_toplot.csv', stringsAsFactors=FALSE)
toPlot_hk <- read.csv('Data/NetworkSimulation/Summer2019/processedData/3wayCompare_holmeKim_sizeChange_regressionData_toplot.csv', stringsAsFactors=FALSE)


toPlot_sf$graph_type <- 'Barabasi-Albert'
toPlot_sworld$graph_type <- 'Watts-Strogatz'
toPlot_hk$graph_type <- 'Holme-Kim'
# ###
# ### RUN ABOVE IF NEED TO LOAD DATA LATER, LIKE WHEN MAKING NARROW XLIM
# ###

# combine
toPlot <- rbind(toPlot_sf, toPlot_sworld, toPlot_hk)
toPlot <- subset(toPlot, sampleSize >= 1000)  # for 3 facets only

### Plot, all
ggplot(toPlot, aes(x=repression_rate_max, y=value, linetype=Significance, color=graph_type)) + geom_rect(aes(xmin=.054, xmax=.1005, ymin=0, ymax=1), fill='grey95', alpha=.05, inherit.aes=FALSE) + geom_line() + facet_wrap(~sampleSize, ncol=3, dir='h') + theme_classic() +  xlab('Repression Rate') + ylab('Percent of Regressions') + theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text=element_text(size=16), axis.title=element_text(size=16), legend.text=element_text(size=14), legend.title = element_text(size=14), strip.text=element_text(size=16), legend.position='bottom', legend.box='vertical') + labs(color='Graph Type') + scale_color_discrete(limits=c('Watts-Strogatz', 'Holme-Kim', 'Barabasi-Albert'))
ggsave('Figures/August2020/compareNetworks_significance_facet3_color_sizeChanged.jpg', plot=last_plot())




ggplot(toPlot, aes(x=repression_rate_max, y=value, linetype=Significance, color=graph_type)) + geom_rect(aes(xmin=.054, xmax=.1005, ymin=0, ymax=1), fill='grey95', alpha=.05, inherit.aes=FALSE) + geom_line() + facet_wrap(~sampleSize, ncol=3, dir='h') + theme_classic() +  xlab('Repression Rate') + ylab('Percent of Regressions') + theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text=element_text(size=16), axis.title=element_text(size=16), legend.text=element_text(size=16), legend.title = element_text(size=14), strip.text=element_text(size=14), legend.position='bottom', legend.box='vertical') + scale_color_grey(limits=c('Watts-Strogatz', 'Holme-Kim', 'Barabasi-Albert')) + labs(color='Graph Type')
ggsave('Figures/August2020/compareNetworks_significance_facet3_grey_sizeChanged.jpg', plot=last_plot())





###################################
### AS NETWORK STRUCTURE CONTROL CHANGES
###
###################################
  # On 08.08.2019, hardcode because sf has different rates
n <- c(100, 250, 500, 1000, 1500, 2000)


rates <- sort(unique(sf$networkStructureControl))
results_sf <- makeResults3(data=sf, repressionRates=NULL, networkSCRates=rates, sampleSizes=n, trials=100, regressionModel=m)
toPlot_sf <- processResults2(data=results_sf, tscore=1.96, networkType='3wayCompare_scaleFree_networkStructureControl')  ## GOT ERROR HERE, 11.08.2019, 10:15 PM.  I LET THE NEXT TWO SEGMENTS RUN, SHOWER, AND SLEEP, SO WILL NEED TO CHECK IF WORKED.
toPlot_sf$graph_type <- 'Barabasi-Albert'

rates <- sort(unique(sworld$networkStructureControl))
results_sworld <- makeResults3(data=sworld, repressionRates=NULL, networkSCRates=rates, sampleSizes=n, trials=100, regressionModel=m)
toPlot_sworld <- processResults2(data=results_sworld, tscore=1.96, networkType='3wayCompare_smallWorld_networkStructureControl')
toPlot_sworld$graph_type <- 'Watts-Strogatz'


rates <- sort(unique(hk$networkStructureControl))
results_hk <- makeResults3(data=hk, repressionRates=NULL, networkSCRates=rates, sampleSizes=n, trials=100, regressionModel=m)
toPlot_hk <- processResults2(data=results_hk, tscore=1.96, networkType='3wayCompare_holmeKim_networkStructureControl')
toPlot_hk$graph_type <- 'Holme-Kim'

toPlot_sf <- read.csv('Data/NetworkSimulation/Summer2019/processedData/3wayCompare_scaleFree_networkStructureControl_networkStructureControl_regressionData_toplot.csv', stringsAsFactors=FALSE)
toPlot_sworld <- read.csv('Data/NetworkSimulation/Summer2019/processedData/3wayCompare_smallWorld_networkStructureControl_networkStructureControl_regressionData_toplot.csv', stringsAsFactors=FALSE)
toPlot_hk <- read.csv('Data/NetworkSimulation/Summer2019/processedData/3wayCompare_holmeKim_networkStructureControl_scalingParameter_regressionData_toplot.csv', stringsAsFactors=FALSE)


## Standardize rates so can be on own plot
toPlot_sf$scaling_parameter_max <- toPlot_sf$scaling_parameter_max - 2  # -2 because range for the others is [0-1], this is [2-3]


toPlot_sf$graph_type <- 'Barabasi-Albert'
toPlot_sworld$graph_type <- 'Watts-Strogatz'
toPlot_hk$graph_type <- 'Holme-Kim'

toPlot <- rbind(toPlot_sf, toPlot_sworld, toPlot_hk)

toPlot <- subset(toPlot, sampleSize>=1000)

### Plot, all, scaling parameter
ggplot(toPlot, aes(x=scaling_parameter_max, y=value, linetype=Significance, color=graph_type)) + geom_line() + facet_wrap(~sampleSize, ncol=3, dir='h') + theme_classic() +  xlab('Network Structure Control') + ylab('Percent of Regressions') + theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text=element_text(size=16), axis.title=element_text(size=16), legend.text=element_text(size=14), legend.title = element_text(size=14), strip.text=element_text(size=16), legend.position='bottom', legend.box='vertical') + labs(color='Graph Type') + scale_color_discrete(limits=c('Watts-Strogatz', 'Holme-Kim', 'Barabasi-Albert'))
ggsave('Figures/August2020/compareNetworks_significance_facet3_color_networkStructureControl.jpg', plot=last_plot())




###################################
### Realistic values
###################################

# rates <- sort(unique(all$repression_rate))
# rates <- rates[rates <= .15]  # Keep .15 to match what new results will be.  Only necessary on 08.08.2019 while wait for Amsterdam cluster to finish
rates <- unique(all$repression_rate)
  # On 08.08.2019, hardcode because sf has different rates
n <- c(100, 250, 500, 1000, 1500, 2000)
#n <- c(100)

# Subset
sf2 <- sf[sf$scaling_parameter==2.3,]
sworld$scaling_parameter <- round(sworld$scaling_parameter, 2)
sworlsf2 <- sworld[sworld$scaling_parameter == .34,]  # sw, rounded global clustering and p, looked for p where global clustering = .16. Looked at Watts-Strogatz original paper, eyeballed .4.  .34 is closest value
hk2 <- hk[hk$scaling_parameter == .315,]

### COMMENT OUT BELOW IF NECESSARY.
### HAVE RUN AS OF 08.08.2019
# make results
# results_sf <- makeResults(data=sf2, repressionRates=rates, sampleSizes=n, trials=100, regressionModel=ma, scalingBool=FALSE, repressionBool=TRUE)
# results_sworld <- makeResults_PolNet_05(data=sworlsf2, repressionRates=rates, sampleSizes=n, trials=100, regressionModel=ma, scalingBool=FALSE, repressionBool=TRUE, networkType='3wayCompare_smallWorld_realisticValue')  # Could not find old data for sworld, so had to run again in August 2020.  Updated function name, so taht's why this function name is different.
# results_hk <- makeResults(data=hk2, repressionRates=rates, sampleSizes=n, trials=100, regressionModel=ma, scalingBool=FALSE, repressionBool=TRUE)

# # get data to plot
# toPlot_sf <- processResults(data=results_sf, tscore=1.96, networkType='3wayCompare_scaleFree_realisticValue')
# toPlot_sworld <- processResults_PolNet_05(data=results_sworld, tscore=1.96, networkType='3wayCompare_smallWorld_realisticValue') # Could not find old data for sworld, so had to run again in August 2020.  Updated function name, so taht's why this function name is different.
# toPlot_hk <- processResults(data=results_hk, tscore=1.96, networkType='3wayCompare_holmeKim_realisticValue')

# toPlot_sf$graph_type <- 'Barabasi-Albert'
# toPlot_sworld$graph_type <- 'Watts-Strogatz'
# toPlot_hk$graph_type <- 'Holme-Kim'

# # combine
# toPlot <- rbind(toPlot_sf, toPlot_sworld, toPlot_hk)

# ### RUN BELOW IF NEED TO LOAD DATA LATER, LIKE WHEN MAKING NARROW XLIM
toPlot_sf <- read.csv('Data/NetworkSimulation/Summer2019/processedData/3wayCompare_scaleFree_realisticValue_regressionData_toplot.csv', stringsAsFactors=FALSE)
toPlot_sworld <- read.csv('Data/NetworkSimulation/Summer2019/processedData/3wayCompare_smallWorld_realisticValue_regressionData_toplot.csv', stringsAsFactors=FALSE)
toPlot_hk <- read.csv('Data/NetworkSimulation/Summer2019/processedData/3wayCompare_holmeKim_realisticValue_regressionData_toplot.csv', stringsAsFactors=FALSE)

toPlot_sf$graph_type <- 'Barabasi-Albert'
toPlot_sworld$graph_type <- 'Watts-Strogatz'
toPlot_hk$graph_type <- 'Holme-Kim'

# # combine
toPlot <- rbind(toPlot_sf, toPlot_sworld, toPlot_hk)
toPlot <- subset(toPlot, sampleSize >= 1000)

### Plot, all
ggplot(toPlot, aes(x=repression_rate_max, y=value, linetype=Significance, color=graph_type)) + geom_rect(aes(xmin=.054, xmax=.1005, ymin=0, ymax=1), fill='grey95', alpha=.05, inherit.aes=FALSE) + geom_line() + facet_wrap(~sampleSize, ncol=3, dir='h') + theme_classic() +  theme_classic() + xlab('Repression Rate') + ylab('Percent of Regressions') + theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text=element_text(size=16), axis.title=element_text(size=16), legend.text=element_text(size=16), legend.title = element_text(size=16), strip.text=element_text(size=16), legend.position='bottom', legend.box='vertical') + labs(color='Graph Type') + scale_color_discrete(limits=c('Watts-Strogatz', 'Holme-Kim', 'Barabasi-Albert'))
ggsave('Figures/August2020/compareNetworks_significance_realisticValue_facet3_color.jpg', plot=last_plot())

ggplot(toPlot, aes(x=repression_rate_max, y=value, linetype=Significance, color=graph_type)) + geom_rect(aes(xmin=.054, xmax=.1005, ymin=0, ymax=1), fill='grey95', alpha=.05, inherit.aes=FALSE) + geom_line() + facet_wrap(~sampleSize, ncol=3, dir='h') + theme_classic() +  xlab('Repression Rate') + ylab('Percent of Regressions') + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_color_grey(limits=c('Watts-Strogatz', 'Holme-Kim', 'Barabasi-Albert')) + labs(color='Graph Type')
ggsave('Figures/August2020/compareNetworks_significance_realisticValue_facet3_grey.jpg', plot=last_plot())


###################################
### Random network from Watts-Strogatz, random from Holme-Kim, scale-free at scaling=3 (closest to H-K)
###################################

# rates <- sort(unique(all$repression_rate))
# rates <- rates[rates <= .15]  # Keep .15 to match what new results will be.  Only necessary on 08.08.2019 while wait for Amsterdam cluster to finish
rates <- unique(all$repression_rate)
  # On 08.08.2019, hardcode because sf has different rates
n <- c(100, 250, 500, 1000, 1500, 2000)

# Subset, have loaded at beginning of script
sf2 <- sf[sf$scaling_parameter==3,]
sworld$scaling_parameter <- round(sworld$scaling_parameter, 2)
sworlsf2 <- sworld[sworld$scaling_parameter == 1,]  # total rewiring
hk2 <- hk[hk$scaling_parameter == 0,]  # no triad clustering, so as close to random as possible


# make results
results_sf <- makeResults(data=sf2, repressionRates=rates, sampleSizes=n, trials=100, regressionModel=ma, scalingBool=FALSE, repressionBool=TRUE)
results_sworld <- makeResults(data=sworlsf2, repressionRates=rates, sampleSizes=n, trials=100, regressionModel=ma, scalingBool=FALSE, repressionBool=TRUE)
results_hk <- makeResults(data=hk2, repressionRates=rates, sampleSizes=n, trials=100, regressionModel=ma, scalingBool=FALSE, repressionBool=TRUE)

# get data to plot
toPlot_sf <- processResults(data=results_sf, tscore=1.96, networkType='3wayCompare_scaleFree_mostSimilar')
toPlot_sworld <- processResults(data=results_sworld, tscore=1.96, networkType='3wayCompare_smallWorld_mostSimilar')
toPlot_hk <- processResults(data=results_hk, tscore=1.96, networkType='3wayCompare_holmeKim_mostSimilar')

toPlot_sf$graph_type <- 'Barabasi-Albert'
toPlot_sworld$graph_type <- 'Watts-Strogatz'
toPlot_hk$graph_type <- 'Holme-Kim'

# combine
toPlot <- rbind(toPlot_sf, toPlot_sworld, toPlot_hk)

### Plot, all + restrict x [0.05,.15]
ggplot(toPlot, aes(x=repression_rate_max, y=value, linetype=Significance, color=graph_type)) + geom_line() + facet_wrap(~sampleSize, ncol=3, dir='h') + theme_classic() +  xlab('Repression Rate') + ylab('Percent of Regressions') + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(color='Graph Type') + scale_color_discrete(limits=c('Watts-Strogatz', 'Holme-Kim', 'Barabasi-Albert'))
ggsave('Figures/June2019/compareNetworks_significance_facet_mostSimilar_color.jpg', plot=last_plot())

ggplot(toPlot, aes(x=repression_rate_max, y=value, linetype=Significance, color=graph_type)) + geom_line() + facet_wrap(~sampleSize, ncol=3, dir='h') + theme_classic() +  xlab('Repression Rate') + ylab('Percent of Regressions') + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_color_grey(limits=c('Watts-Strogatz', 'Holme-Kim', 'Barabasi-Albert')) + labs(color='Graph Type')
ggsave('Figures/June2019/compareNetworks_significance_facet_mostSimilar_grey.jpg', plot=last_plot())




###################################
### Random network from Watts-Strogatz, random from Holme-Kim, scale-free at scaling=2 (most skew)
###################################

# rates <- sort(unique(all$repression_rate))
# rates <- rates[rates <= .15]  # Keep .15 to match what new results will be.  Only necessary on 08.08.2019 while wait for Amsterdam cluster to finish
rates <- unique(all$repression_rate)
  # On 08.08.2019, hardcode because sf has different rates
n <- c(100, 250, 500, 1000, 1500, 2000)

# Subset, have loaded at beginning of script
sf2 <- sf[sf$scaling_parameter==2,]
sworld$scaling_parameter <- round(sworld$scaling_parameter, 2)
sworlsf2 <- sworld[sworld$scaling_parameter == 1,]  # total rewiring
hk2 <- hk[hk$scaling_parameter == 0,]  # no triad clustering, so as close to random as possible


# make results
results_sf <- makeResults(data=sf2, repressionRates=rates, sampleSizes=n, trials=100, regressionModel=ma, scalingBool=FALSE, repressionBool=TRUE)
results_sworld <- makeResults(data=sworlsf2, repressionRates=rates, sampleSizes=n, trials=100, regressionModel=ma, scalingBool=FALSE, repressionBool=TRUE)
results_hk <- makeResults(data=hk2, repressionRates=rates, sampleSizes=n, trials=100, regressionModel=ma, scalingBool=FALSE, repressionBool=TRUE)

# get data to plot
toPlot_sf <- processResults(data=results_sf, tscore=1.96, networkType='3wayCompare_scaleFree_mostSimilarV2')
toPlot_sworld <- processResults(data=results_sworld, tscore=1.96, networkType='3wayCompare_smallWorld_mostSimilarV2')
toPlot_hk <- processResults(data=results_hk, tscore=1.96, networkType='3wayCompare_holmeKim_mostSimilarV2')

toPlot_sf$graph_type <- 'Barabasi-Albert'
toPlot_sworld$graph_type <- 'Watts-Strogatz'
toPlot_hk$graph_type <- 'Holme-Kim'

# combine
toPlot <- rbind(toPlot_sf, toPlot_sworld, toPlot_hk)

### Plot, all + restrict x [0.05,.15]
ggplot(toPlot, aes(x=repression_rate_max, y=value, linetype=Significance, color=graph_type)) + geom_line() + facet_wrap(~sampleSize, ncol=3, dir='h') + theme_classic() +  xlab('Repression Rate') + ylab('Percent of Regressions') + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(color='Graph Type') + scale_color_discrete(limits=c('Watts-Strogatz', 'Holme-Kim', 'Barabasi-Albert'))
ggsave('Figures/June2019/compareNetworks_significance_facet_mostSimilarV2_color.jpg', plot=last_plot())

ggplot(toPlot, aes(x=repression_rate_max, y=value, linetype=Significance, color=graph_type)) + geom_line() + facet_wrap(~sampleSize, ncol=3, dir='h') + theme_classic() +  xlab('Repression Rate') + ylab('Percent of Regressions') + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_color_grey(limits=c('Watts-Strogatz', 'Holme-Kim', 'Barabasi-Albert')) + labs(color='Graph Type')
ggsave('Figures/June2019/compareNetworks_significance_facet_mostSimilarV2_grey.jpg', plot=last_plot())




###################################
### Random network from Watts-Strogatz, normal Holme-Kim, scale-free at scaling=2 (most skew)
###################################

# rates <- sort(unique(all$repression_rate))
# rates <- rates[rates <= .15]  # Keep .15 to match what new results will be.  Only necessary on 08.08.2019 while wait for Amsterdam cluster to finish
rates <- unique(all$repression_rate)
  # On 08.08.2019, hardcode because sf has different rates
n <- c(100, 250, 500, 1000, 1500, 2000)

# Subset, have loaded at beginning of script
sf2 <- sf[sf$scaling_parameter==2,]
sworld$scaling_parameter <- round(sworld$scaling_parameter, 2)
sworld2 <- sworld[sworld$scaling_parameter == 1,]  # total rewiring
hk2 <- hk # keep all


# make results
# results_sf <- makeResults(data=sf2, repressionRates=rates, sampleSizes=n, trials=100, regressionModel=ma, scalingBool=FALSE, repressionBool=TRUE)
# results_sworld <- makeResults(data=sworlsf2, repressionRates=rates, sampleSizes=n, trials=100, regressionModel=ma, scalingBool=FALSE, repressionBool=TRUE)
# results_hk <- makeResults(data=hk2, repressionRates=rates, sampleSizes=n, trials=100, regressionModel=ma, scalingBool=FALSE, repressionBool=TRUE)

# # get data to plot
# toPlot_sf <- processResults(data=results_sf, tscore=1.96, networkType='3wayCompare_scaleFree_mostSimilarV3')
# toPlot_sworld <- processResults(data=results_sworld, tscore=1.96, networkType='3wayCompare_smallWorld_mostSimilarV3')
# toPlot_hk <- processResults(data=results_hk, tscore=1.96, networkType='3wayCompare_holmeKim_mostSimilarV3')

# toPlot_sf$graph_type <- 'Barabasi-Albert'
# toPlot_sworld$graph_type <- 'Watts-Strogatz'
# toPlot_hk$graph_type <- 'Holme-Kim'

### RUN BELOW IF NEED TO LOAD DATA LATER, LIKE WHEN MAKING NARROW XLIM
toPlot_sf <- read.csv('Data/NetworkSimulation/Summer2019/processedData/3wayCompare_scaleFree_mostSimilarV3_regressionData_toplot.csv', stringsAsFactors=FALSE)
toPlot_sworld <- read.csv('Data/NetworkSimulation/Summer2019/processedData/3wayCompare_smallWorld_mostSimilarV3_regressionData_toplot.csv', stringsAsFactors=FALSE)
toPlot_hk <- read.csv('Data/NetworkSimulation/Summer2019/processedData/3wayCompare_holmeKim_mostSimilarV3_regressionData_toplot.csv', stringsAsFactors=FALSE)

toPlot_sf$graph_type <- 'Barabasi-Albert'
toPlot_sworld$graph_type <- 'Watts-Strogatz'
toPlot_hk$graph_type <- 'Holme-Kim'

# combine
toPlot <- rbind(toPlot_sf, toPlot_sworld, toPlot_hk)

### Plot, all + restrict x [0.05,.15]
ggplot(toPlot, aes(x=repression_rate_max, y=value, linetype=Significance, color=graph_type)) + geom_line() + facet_wrap(~sampleSize, ncol=3, dir='h') + theme_classic() +  xlab('Repression Rate') + ylab('Percent of Regressions') + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(color='Graph Type') + scale_color_discrete(limits=c('Watts-Strogatz', 'Holme-Kim', 'Barabasi-Albert'))
ggsave('Figures/June2019/compareNetworks_significance_facet_mostSimilarV3_color.jpg', plot=last_plot())

ggplot(toPlot, aes(x=repression_rate_max, y=value, linetype=Significance, color=graph_type)) + geom_line() + facet_wrap(~sampleSize, ncol=3, dir='h') + theme_classic() +  xlab('Repression Rate') + ylab('Percent of Regressions') + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_color_grey(limits=c('Watts-Strogatz', 'Holme-Kim', 'Barabasi-Albert')) + labs(color='Graph Type')
ggsave('Figures/June2019/compareNetworks_significance_facet_mostSimilarV3_grey.jpg', plot=last_plot())





###################################
### Random network from Watts-Strogatz, random from Holme-Kim, scale-free at scaling=3 (closest to H-K)
###################################
####
# NB: ON 08.14.2020 AT 1;30 P.M., I RERUN THIS BUT WITH AN UPDATED makeResults FROM Scripts/PolNet_08_RegressionResultsOverTime_v8.R THAT SAVES OUT DATA.  THAT WRITE.CSV CAPABILITY IS ADDED SO DO NOT HAVE CHANGE FUNCTION IN THAT FILE BUT CAN USE THAT FILE TO READ IN DATA MADE HERE.  NB AS WELL THAT I WRITE OUT TO OLD DIRECTORY, SUMMER2019, IS FINE SINCE ALL THAT CHANGVES IS WRITING TO FILE.
# don't forget model ma
# rates <- sort(unique(all$repression_rate))
# rates <- rates[rates <= .15]  # Keep .15 to match what new results will be.  Only necessary on 08.08.2019 while wait for Amsterdam cluster to finish
rates <- unique(all$repression_rate)
  # On 08.08.2019, hardcode because sf has different rates
n <- c(100, 250, 500, 1000, 1500, 2000)

# Subset, have loaded at beginning of script
sf2 <- sf[sf$scaling_parameter==3,]
# sworld$scaling_parameter <- round(sworld$scaling_parameter, 2)
# sworld2 <- sworld[sworld$scaling_parameter == 1,]  # total rewiring
# hk2 <- hk[hk$scaling_parameter == 0,]  # no triad clustering, so as close to random as possible
sworld2 <- sworld
hk2 <- hk


# make results
results_sf <- makeResults_PolNet_05(data=sf2, repressionRates=rates, sampleSizes=n, trials=100, regressionModel=ma, scalingBool=FALSE, repressionBool=TRUE, networkType='3wayCompare_scaleFree_mostSimilar')
results_sworld <- makeResults_PolNet_05(data=sworld2, repressionRates=rates, sampleSizes=n, trials=100, regressionModel=ma, scalingBool=FALSE, repressionBool=TRUE, networkType='3wayCompare_smallWorld_mostSimilar')
results_hk <- makeResults_PolNet_05(data=hk2, repressionRates=rates, sampleSizes=n, trials=100, regressionModel=ma, scalingBool=FALSE, repressionBool=TRUE, networkType='3wayCompare_holmeKim_mostSimilar')

# get data to plot
toPlot_sf <- processResults_PolNet_05(data=results_sf, tscore=1.96, networkType='3wayCompare_scaleFree_mostSimilar')
toPlot_sworld <- processResults_PolNet_05(data=results_sworld, tscore=1.96, networkType='3wayCompare_smallWorld_mostSimilar')
toPlot_hk <- processResults_PolNet_05(data=results_hk, tscore=1.96, networkType='3wayCompare_holmeKim_mostSimilar')

### RUN BELOW IF NEED TO LOAD DATA LATER, LIKE WHEN MAKING NARROW XLIM
toPlot_sf <- read.csv('Data/NetworkSimulation/Summer2019/processedData/3wayCompare_scaleFree_mostSimilar_regressionData_toplot.csv', stringsAsFactors=FALSE)
toPlot_sworld <- read.csv('Data/NetworkSimulation/Summer2019/processedData/3wayCompare_smallWorld_mostSimilar_regressionData_toplot.csv', stringsAsFactors=FALSE)
toPlot_hk <- read.csv('Data/NetworkSimulation/Summer2019/processedData/3wayCompare_holmeKim_mostSimilar_regressionData_toplot.csv', stringsAsFactors=FALSE)



toPlot_sf$graph_type <- 'Barabasi-Albert'
toPlot_sworld$graph_type <- 'Watts-Strogatz'
toPlot_hk$graph_type <- 'Holme-Kim'

# combine
toPlot <- rbind(toPlot_sf, toPlot_sworld, toPlot_hk)
toPlot <- subset(toPlot, sampleSize >= 1000)  # will ensure facets are only 3

### Plot, all + restrict x [0.05,.15]
ggplot(toPlot, aes(x=repression_rate_max, y=value, linetype=Significance, color=graph_type)) + geom_rect(aes(xmin=.054, xmax=.1005, ymin=0, ymax=1), fill='grey95', alpha=.05, inherit.aes=FALSE) + geom_line() + facet_wrap(~sampleSize, ncol=3, dir='h') + theme_classic() +  xlab('Repression Rate') + ylab('Percent of Regressions') + theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text=element_text(size=16), axis.title=element_text(size=16), legend.text=element_text(size=14), legend.title = element_text(size=14), strip.text=element_text(size=16), legend.position='bottom', legend.box='vertical') + labs(color='Graph Type') + scale_color_discrete(limits=c('Watts-Strogatz', 'Holme-Kim', 'Barabasi-Albert'))
ggsave('Figures/August2020/compareNetworks_significance_facet3_mostSimilarV4_color.jpg', plot=last_plot())

ggplot(toPlot, aes(x=repression_rate_max, y=value, linetype=Significance, color=graph_type)) + geom_line() + facet_wrap(~sampleSize, ncol=3, dir='h') + theme_classic() +  xlab('Repression Rate') + ylab('Percent of Regressions') + theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text=element_text(size=16), axis.title=element_text(size=16), legend.text=element_text(size=14), legend.title = element_text(size=14), strip.text=element_text(size=16), legend.position='bottom', legend.box='vertical') + scale_color_grey(limits=c('Watts-Strogatz', 'Holme-Kim', 'Barabasi-Albert')) + labs(color='Graph Type')
ggsave('Figures/August2020/compareNetworks_significance_facet3_mostSimilarV4_grey.jpg', plot=last_plot())



###################################
### Keep any start that is equal or below median size
###################################
sf2 <- sf[sf$normalized <= median(sf$normalized),]
hk2 <- hk[hk$normalized <= median(hk$normalized),]
sworld2 <- sworld[sworld$normalized <= median(sworld$normalized),]


rates <- unique(sf2$repression_rate)
n <- c(100, 250, 500, 1000, 1500, 2000)


# make results
results_sf <- makeResults(data=sf2, repressionRates=rates, sampleSizes=n, trials=100, regressionModel=m, scalingBool=FALSE, repressionBool=TRUE)
results_sworld <- makeResults(data=sworld2, repressionRates=rates, sampleSizes=n, trials=100, regressionModel=m, scalingBool=FALSE, repressionBool=TRUE)
results_hk <- makeResults(data=hk2, repressionRates=rates, sampleSizes=n, trials=100, regressionModel=m, scalingBool=FALSE, repressionBool=TRUE)

# # get data to plot
toPlot_sf <- processResults(data=results_sf, tscore=1.96, networkType='3wayCompare_scaleFree_lteMedian')
toPlot_sworld <- processResults(data=results_sworld, tscore=1.96, networkType='3wayCompare_smallWorld_lteMedian')
toPlot_hk <- processResults(data=results_hk, tscore=1.96, networkType='3wayCompare_holmeKim_lteMedian')


# Below is to read if have not run toPlot line above
# toPlot_sf <- read.csv('Data/NetworkSimulation/Summer2019/processedData/3wayCompare_3wayCompare_scaleFree_lteMedian_regressionData_toplot.csv', stringsAsFactors=FALSE)
# toPlot_sworld <- read.csv('Data/NetworkSimulation/Summer2019/processedData/3wayCompare_smallWorld_lteMedian_mostSimilar_regressionData_toplot.csv', stringsAsFactors=FALSE)
# toPlot_hk <- read.csv('Data/NetworkSimulation/Summer2019/processedData/3wayCompare_holmeKim_lteMedian_regressionData_toplot.csv', stringsAsFactors=FALSE)



toPlot_sf$graph_type <- 'Barabasi-Albert'
toPlot_sworld$graph_type <- 'Watts-Strogatz'
toPlot_hk$graph_type <- 'Holme-Kim'

# combine
toPlot <- rbind(toPlot_sf, toPlot_sworld, toPlot_hk)


#toPlot$facets <- paste0('Sample Size: ', toPlot$sampleSize, '\n', 'Scaling Parameter = 2.3')
ggplot(toPlot, aes(x=repression_rate_max, y=value, linetype=Significance, color=graph_type)) + geom_line() + facet_wrap(~sampleSize, ncol=3, dir='h') + theme_classic()  + xlab('Repression Rate') + ylab('Percent of Regressions') + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave('Figures/June2019/compareNetworks_signficance_facet_initialSize_lteMedian_color.jpg', plot=last_plot())

ggplot(toPlot, aes(x=repression_rate_max, y=value, linetype=Significance, color=graph_type)) + geom_line() + facet_wrap(~sampleSize, ncol=3, dir='h') + theme_classic()  + xlab('Repression Rate') + ylab('Percent of Regressions') + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_color_grey(limits=c('Watts-Strogatz', 'Holme-Kim', 'Barabasi-Albert'))
ggsave('Figures/June2019/compareNetworks_signficance_facet_initialSize_lteMedian_grey.jpg', plot=last_plot())





###################################
### Keep any start that is equal or below median size and size changes
### See phase shift for scale-free, watts-strogatz, no repression effect for holme-kim
###################################
sf2 <- sf[sf$normalized <= median(sf$normalized),]
hk2 <- hk[hk$normalized <= median(hk$normalized),]
sworld2 <- sworld[sworld$normalized <= median(sworld$normalized),]

sf2 <- sf2[sf2$initial_size != sf2$active_nodes,]
hk2 <- hk2[hk2$initial_size != hk2$active_nodes,]
sf2 <- sworld2[sworld2$initial_size != sworld2$active_nodes,]



rates <- unique(sf2$repression_rate)
n <- c(100, 250, 500, 1000, 1500, 2000)


# make results
results_sf <- makeResults(data=sf2, repressionRates=rates, sampleSizes=n, trials=100, regressionModel=m, scalingBool=FALSE, repressionBool=TRUE)
results_sworld <- makeResults(data=sworld2, repressionRates=rates, sampleSizes=n, trials=100, regressionModel=m, scalingBool=FALSE, repressionBool=TRUE)
results_hk <- makeResults(data=hk2, repressionRates=rates, sampleSizes=n, trials=100, regressionModel=m, scalingBool=FALSE, repressionBool=TRUE)

# # get data to plot
toPlot_sf <- processResults(data=results_sf, tscore=1.96, networkType='3wayCompare_scaleFree_lteMedian_sizeChanged')
toPlot_sworld <- processResults(data=results_sworld, tscore=1.96, networkType='3wayCompare_smallWorld_lteMedian_sizeChanged')
toPlot_hk <- processResults(data=results_hk, tscore=1.96, networkType='3wayCompare_holmeKim_lteMedian_sizeChanged')




ggplot(toPlot, aes(x=repression_rate_max, y=value, linetype=Significance, color=graph_type)) + geom_line() + facet_wrap(~sampleSize, ncol=3, dir='h') + theme_classic()  + xlab('Repression Rate') + ylab('Percent of Regressions') + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_color_grey(limits=c('Watts-Strogatz', 'Holme-Kim', 'Barabasi-Albert'))
ggsave('Figures/June2019/compareNetworks_signficance_facet_initialSize_lteMedian_sizeChanged_grey.jpg', plot=last_plot())


ggplot(toPlot, aes(x=repression_rate_max, y=value, linetype=Significance, color=graph_type)) + geom_rect(aes(xmin=.054, xmax=.1005, ymin=0, ymax=1), fill='grey95', alpha=.05, inherit.aes=FALSE) +



###################################
#
# 	USING NEW TYPES OF REPRESSION
#
###################################
#### NEW WEIGHTING
all <- read.csv('Data/NetworkSimulation/Summer2020/processedData/01_experimentsCombined_eigenCore_newRepression.csv', stringsAsFactors=FALSE)
all$normalized <- round(1000*all$active_nodes/all$total_nodes)

# Make names to match model m
all$scaling_parameter <- all$networkStructureControl

rates <- unique(all$repression_rate)
  # On 08.08.2019, hardcode because sf has different rates
n <- c(100, 250, 500, 1000, 1500, 2000)


# Separate out the networks
sf <- all[all$graph_type=='scale_free_graph',]
sworld <- all[all$graph_type=='watts_strogatz_graph',]
hk <- all[all$graph_type=='powerlaw_cluster_graph',]

# make results
results_sf <- makeResults_PolNet_05(data=sf, repressionRates=rates, sampleSizes=n, trials=100, regressionModel=ma, scalingBool=FALSE, repressionBool=TRUE, networkType='3wayCompare_forOldPresentation_scaleFree_newRepression')
results_sworld <- makeResults_PolNet_05(data=sworld, repressionRates=rates, sampleSizes=n, trials=100, regressionModel=ma, scalingBool=FALSE, repressionBool=TRUE, networkType='3wayCompare_forOldPresentation_smallWorld_newRepression')
results_hk <- makeResults_PolNet_05(data=hk, repressionRates=rates, sampleSizes=n, trials=100, regressionModel=ma, scalingBool=FALSE, repressionBool=TRUE, networkType='3wayCompare_forOldPresentation_holmeKim_newRepression')

# get data to plot
toPlot_sf <- processResults_PolNet_05(data=results_sf, tscore=1.96, networkType='3wayCompare_forOldPresentation_scaleFree_newRepression')
toPlot_sworld <- processResults_PolNet_05(data=results_sworld, tscore=1.96, networkType='3wayCompare_forOldPresentation_smallWorld_newRepression')
toPlot_hk <- processResults_PolNet_05(data=results_hk, tscore=1.96, networkType='3wayCompare_forOldPresentation_holmeKim_newRepression')

### RUN BELOW IF NEED TO LOAD DATA LATER, LIKE WHEN MAKING NARROW XLIM
# toPlot_sf <- read.csv('Data/NetworkSimulation/Summer2019/processedData/3wayCompare_forOldPresentation_scaleFree_newRepression_regressionData_toplot.csv', stringsAsFactors=FALSE)
# toPlot_sworld <- read.csv('Data/NetworkSimulation/Summer2019/processedData/3wayCompare_forOldPresentation_smallWorld_newRepression_regressionData_toplot.csv', stringsAsFactors=FALSE)
# toPlot_hk <- read.csv('Data/NetworkSimulation/Summer2019/processedData/3wayCompare_forOldPresentation_holmeKim_newRepression_regressionData_toplot.csv', stringsAsFactors=FALSE)



toPlot_sf$graph_type <- 'Barabasi-Albert'
toPlot_sworld$graph_type <- 'Watts-Strogatz'
toPlot_hk$graph_type <- 'Holme-Kim'

# combine
toPlot <- rbind(toPlot_sf, toPlot_sworld, toPlot_hk)
toPlot <- subset(toPlot, sampleSize >= 1000)  # will ensure facets are only 3

### Plot, all + restrict x [0.05,.15]
ggplot(toPlot, aes(x=repression_rate_max, y=value, linetype=Significance, color=graph_type)) + geom_rect(aes(xmin=.054, xmax=.1005, ymin=0, ymax=1), fill='grey95', alpha=.05, inherit.aes=FALSE) + geom_line() + facet_wrap(~sampleSize, ncol=3, dir='h') + theme_classic() +  xlab('Repression Rate') + ylab('Percent of Regressions') + theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text=element_text(size=16), axis.title=element_text(size=16), legend.text=element_text(size=14), legend.title = element_text(size=14), strip.text=element_text(size=16), legend.position='bottom', legend.box='vertical') + labs(color='Graph Type') + scale_color_discrete(limits=c('Watts-Strogatz', 'Holme-Kim', 'Barabasi-Albert'))
ggsave('Figures/August2020/compareNetworks_significance_facet3_newRepression_color.jpg', plot=last_plot())

ggplot(toPlot, aes(x=repression_rate_max, y=value, linetype=Significance, color=graph_type)) + geom_line() + facet_wrap(~sampleSize, ncol=3, dir='h') + theme_classic() +  xlab('Repression Rate') + ylab('Percent of Regressions') + theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text=element_text(size=16), axis.title=element_text(size=16), legend.text=element_text(size=14), legend.title = element_text(size=14), strip.text=element_text(size=16), legend.position='bottom', legend.box='vertical') + scale_color_grey(limits=c('Watts-Strogatz', 'Holme-Kim', 'Barabasi-Albert')) + labs(color='Graph Type')
ggsave('Figures/August2020/compareNetworks_significance_facet3_newRepression_grey.jpg', plot=last_plot())




#### EDGE REPRESSION
all <- read.csv('Data/NetworkSimulation/Summer2018/processedData/01_experimentsCombined_edgeRepression.csv', stringsAsFactors=FALSE)
all$normalized <- round(1000*all$final_size/all$total_nodes)

# Make names to match model m
all$scaling_parameter <- all$networkStructureControl

rates <- unique(all$repression_rate)
  # On 08.08.2019, hardcode because sf has different rates
n <- c(100, 250, 500, 1000, 1500, 2000)


# Separate out the networks
sf <- all[all$graph_type=='scale_free_graph',]
sworld <- all[all$graph_type=='watts_strogatz_graph',]
hk <- all[all$graph_type=='powerlaw_cluster_graph',]



# make results
results_sf <- makeResults_PolNet_05(data=sf, repressionRates=rates, sampleSizes=n, trials=100, regressionModel=ma, scalingBool=FALSE, repressionBool=TRUE, networkType='3wayCompare_forOldPresentation_scaleFree_edgeRepression')
results_sworld <- makeResults_PolNet_05(data=sworld, repressionRates=rates, sampleSizes=n, trials=100, regressionModel=ma, scalingBool=FALSE, repressionBool=TRUE, networkType='3wayCompare_forOldPresentation_smallWorld_edgeRepression')
results_hk <- makeResults_PolNet_05(data=hk, repressionRates=rates, sampleSizes=n, trials=100, regressionModel=ma, scalingBool=FALSE, repressionBool=TRUE, networkType='3wayCompare_forOldPresentation_holmeKim_edgeRepression')

# get data to plot
toPlot_sf <- processResults_PolNet_05(data=results_sf, tscore=1.96, networkType='3wayCompare_forOldPresentation_scaleFree_edgeRepression')
toPlot_sworld <- processResults_PolNet_05(data=results_sworld, tscore=1.96, networkType='3wayCompare_forOldPresentation_smallWorld_edgeRepression')
toPlot_hk <- processResults_PolNet_05(data=results_hk, tscore=1.96, networkType='3wayCompare_forOldPresentation_holmeKim_edgeRepression')

### RUN BELOW IF NEED TO LOAD DATA LATER, LIKE WHEN MAKING NARROW XLIM
# toPlot_sf <- read.csv('Data/NetworkSimulation/Summer2019/processedData/3wayCompare_forOldPresentation_scaleFree_edgeRepression_regressionData_toplot.csv', stringsAsFactors=FALSE)
# toPlot_sworld <- read.csv('Data/NetworkSimulation/Summer2019/processedData/3wayCompare_forOldPresentation_smallWorld_edgeRepression_regressionData_toplot.csv', stringsAsFactors=FALSE)
# toPlot_hk <- read.csv('Data/NetworkSimulation/Summer2019/processedData/3wayCompare_forOldPresentation_holmeKim_edgeRepression_regressionData_toplot.csv', stringsAsFactors=FALSE)



toPlot_sf$graph_type <- 'Barabasi-Albert'
toPlot_sworld$graph_type <- 'Watts-Strogatz'
toPlot_hk$graph_type <- 'Holme-Kim'

# combine
toPlot <- rbind(toPlot_sf, toPlot_sworld, toPlot_hk)
toPlot <- subset(toPlot, sampleSize >= 1000)  # will ensure facets are only 3

### Plot, all + restrict x [0.05,.15]
ggplot(toPlot, aes(x=repression_rate_max, y=value, linetype=Significance, color=graph_type)) + geom_rect(aes(xmin=.054, xmax=.1005, ymin=0, ymax=1), fill='grey95', alpha=.05, inherit.aes=FALSE) + geom_line() + facet_wrap(~sampleSize, ncol=3, dir='h') + theme_classic() +  xlab('Repression Rate') + ylab('Percent of Regressions') + theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text=element_text(size=16), axis.title=element_text(size=16), legend.text=element_text(size=14), legend.title = element_text(size=14), strip.text=element_text(size=16), legend.position='bottom', legend.box='vertical') + labs(color='Graph Type') + scale_color_discrete(limits=c('Watts-Strogatz', 'Holme-Kim', 'Barabasi-Albert')) + xlim(c(0,.15))
ggsave('Figures/August2020/compareNetworks_significance_facet3_edgeRepression_color.jpg', plot=last_plot())

ggplot(toPlot, aes(x=repression_rate_max, y=value, linetype=Significance, color=graph_type)) + geom_line() + facet_wrap(~sampleSize, ncol=3, dir='h') + theme_classic() +  xlab('Repression Rate') + ylab('Percent of Regressions') + theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text=element_text(size=16), axis.title=element_text(size=16), legend.text=element_text(size=14), legend.title = element_text(size=14), strip.text=element_text(size=16), legend.position='bottom', legend.box='vertical') + scale_color_grey(limits=c('Watts-Strogatz', 'Holme-Kim', 'Barabasi-Albert')) + labs(color='Graph Type')
ggsave('Figures/August2020/compareNetworks_significance_facet3_edgeRepression_grey.jpg', plot=last_plot())

