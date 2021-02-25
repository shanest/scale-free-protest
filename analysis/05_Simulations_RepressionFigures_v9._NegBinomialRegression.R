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
v9_NegBinomialRegression. Is a fork of v9 to replicate main results but not using OLS, as robustness.
	1. Change outfile names. - DONE
		- CHANGED OUTPATH OF write.csv in processResults, processResults2
		- CHANGED GGSAVE OUTPATH TO June2019/negBinom_
	2. Update regression functions
	3. Keep only parts of code that make figures in paper as of 11.11.2019. - DONE

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
#library(MASS)  # For glm.nb.  Fails to converge a lot.
library(arm)  # bayes glm


setwd('/Users/Zack/Documents/UCLA/Research/ProtestTheory')


# Runs regression, returns t-statistic of each coefficient on final size
runRegressions <- function(data, model, sampleSize, family, scaling=FALSE, repression=TRUE){
	#keep <- round(nrow(data)*percentage)
	sub <- sample_n(data, sampleSize, replace=TRUE)

	meanRepression <- mean(sub$repression_rate)
	if(family != 'negbinom'){
		#reg <- summary(lm(model, sub))$coefficients
		reg <- summary(glm(model, data=sub))$coefficients  # bayes glm, from arm package
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
m <- normalized ~ initial_density + initial_global_clustering + initial_mean_degree + initial_median_degree + initial_neighborhood_clustering + initial_nodes_clustering + initial_size + scaling_parameter + repression_rate

# Where scaling parameter not vary
ma <- normalized ~ initial_density + initial_global_clustering + initial_mean_degree + initial_median_degree + initial_neighborhood_clustering + initial_nodes_clustering + initial_size  + repression_rate



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
makeResults <- function(data, repressionRates, sampleSizes, trials=100, regressionModel, scalingBool=TRUE, repressionBool=TRUE){
	results <- NULL


	for(i in 2:length(repressionRates)){  # Use 2 so that will always have variation of repression rate
		print(paste0('Repression rate is ', repressionRates[i]))
		temp <- data[data$repression_rate<=repressionRates[i],]

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
	return(results)
}


makeResults3 <- function(data, repressionRates=NULL, networkSCRates=NULL, sampleSizes, trials=100, regressionModel, scalingBool=TRUE, repressionBool=TRUE){
	results <- NULL

	if(!is.null(repressionRates)){
		for(i in 2:length(repressionRates)){  # Use 2 so that will always have variation of repression rate
			print(paste0('Repression rate is ', repressionRates[i]))
			temp <- data[data$repression_rate<=repressionRates[i],]

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

	return(results)
}



makeResults2 <- function(data, repressionRates=NULL, networkSCRates=NULL sampleSizes, trials=100, regressionModel, scalingBool=TRUE, repressionBool=TRUE){
	results <- NULL

	#if(!is.null(repressionRates)){
		for(i in 2:length(repressionRates)){  # Use 2 so that will always have variation of repression rate
			print(paste0('Repression rate is ', repressionRates[i]))
			temp <- data[data$repression_rate<=repressionRates[i],]

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
			return(results)
	#}
}

	#if(!is.null(networkSCRates)){
		# for(i in 2:length(networkSCRates)){  # Use 2 so that will always have variation of repression rate
		# 	print(paste0('Network structure control rate is ', networkSCRates[i]))
		# 	temp <- data[data$networkStructureControl<=networkSCRates[i],]

		# 	for(j in 1:length(sampleSizes)){
		# 		for(k in 1:trials){
		# 			r <- runRegressions(data=temp, model=regressionModel, sampleSize=sampleSizes[j], family='regular', scaling=scalingBool, repression=repressionBool)
		# 			r$networkSCRate_max <- networkSCRates[i]
		# 			#r$percSampled <- percUse[j]
		# 			results <- rbind(results, r)
		# 		}
		# 	print(paste0('Just finished sample size of ', sampleSizes[j]))
		# 	}
		# }
	#}
}


# data = the results from a big loop, all the regressions
# tscore = level of significance wanted, make sure is positive
# networkType = string describing network being analyzed
processResults <- function(data, tscore, networkType){
	tscore <- tscore

	data$repressionSigSign <- cut(data$repression_rate, breaks=c(-Inf, -1*tscore, tscore, Inf), labels=c(-1,0,1), right=FALSE, include.lowest=TRUE)

	fileout <- paste0('Data/NetworkSimulation/Summer2019/processedData/negBinom_', networkType, '_regressionData.csv')
	write.csv(data, fileout, row.names=FALSE)

	temp <- data.frame(data %>% group_by(repression_rate_max, sampleSize) %>% summarize(Neg = sum(repressionSigSign==-1)/n(), Zero = sum(repressionSigSign==0)/n(), Pos = sum(repressionSigSign==1)/n()))

	temp <- melt(temp, id.vars=c('repression_rate_max', 'sampleSize'))
	temp$Significance <- temp$variable

	fileout_agg <- paste0('Data/NetworkSimulation/Summer2019/processedData/negBinom_', networkType, '_regressionData_toplot.csv')
	write.csv(temp, fileout_agg, row.names=FALSE)

	return(temp)
}


# data = the results from a big loop, all the regressions
# tscore = level of significance wanted, make sure is positive
# networkType = string describing network being analyzed
# 11.08.2019, seems to be the same as processResults but for scaling paraemter, not for repression rate
processResults2 <- function(data, tscore, networkType){
	tscore <- tscore

	data$SigSign <- cut(data$scaling_parameter, breaks=c(-Inf, -1*tscore, tscore, Inf), labels=c(-1,0,1), right=FALSE, include.lowest=TRUE)

	fileout <- paste0('Data/NetworkSimulation/Summer2019/processedData/negBinom_', networkType, '_scalingParameter_regressionData.csv')
	write.csv(data, fileout, row.names=FALSE)

	#data$grouping <- data[[variable]]  #ggplot will not read string from variable

	temp <- data.frame(data %>% group_by(scaling_parameter_max, sampleSize) %>% summarize(Neg = sum(SigSign==-1)/n(), Zero = sum(SigSign==0)/n(), Pos = sum(SigSign==1)/n()))

	temp <- melt(temp, id.vars=c('scaling_parameter_max', 'sampleSize'))  # variable is a user supplied string, so will not be same as variable in following line
	temp$Significance <- temp$variable

	fileout_agg <- paste0('Data/NetworkSimulation/Summer2019/processedData/negBinom_', networkType, '_scalingParameter_regressionData_toplot.csv')
	write.csv(temp, fileout_agg, row.names=FALSE)

	return(temp)
}





##########################
##
##
##	WORK, ONE GRAPH WITH 1 LINE PER NETWORK TYPE
##
##
##########################
### Standardize some information
#sf$initial_size <- log(sf$initial_size+1, 10)  # Use log since skewed
sf$scaling_parameter <- sf$networkStructureControl  # Change name to match model variable name

#sworld$initial_size <- log(sworld$initial_size+1, 10)
sworld$scaling_parameter <- sworld$networkStructureControl

#hk$initial_size <- log10(hk$initial_size+1)
hk$scaling_parameter <- hk$networkStructureControl


###################################
### Regardless of network scaling
###################################

# rates <- sort(unique(all$repression_rate))
# rates <- rates[rates <= .15]  # Keep .15 to match what new results will be.  Only necessary on 08.08.2019 while wait for Amsterdam cluster to finish
rates <- unique(all$repression_rate)
  # On 08.08.2019, hardcode because sf has different rates
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
# toPlot_sf <- read.csv('Data/NetworkSimulation/Summer2019/processedData/negBinom_3wayCompare_scaleFree_regressionData_toplot.csv', stringsAsFactors=FALSE)
# toPlot_sworld <- read.csv('Data/NetworkSimulation/Summer2019/processedData/negBinom_3wayCompare_smallWorld_regressionData_toplot.csv', stringsAsFactors=FALSE)
# toPlot_hk <- read.csv('Data/NetworkSimulation/Summer2019/processedData/negBinom_3wayCompare_holmeKim_regressionData_toplot.csv', stringsAsFactors=FALSE)

toPlot_sf$graph_type <- 'Barabasi-Albert'
toPlot_sworld$graph_type <- 'Watts-Strogatz'
toPlot_hk$graph_type <- 'Holme-Kim'
# ###
# ### RUN ABOVE IF NEED TO LOAD DATA LATER, LIKE WHEN MAKING NARROW XLIM
# ###

# combine
toPlot <- rbind(toPlot_sf, toPlot_sworld, toPlot_hk)

### Plot, all
ggplot(toPlot, aes(x=repression_rate_max, y=value, linetype=Significance, color=graph_type)) + geom_line() + facet_wrap(~sampleSize, ncol=3, dir='h') + theme_bw() +  xlab('Repression Rate') + ylab('Percent of Regressions') + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(color='Graph Type') + scale_color_discrete(limits=c('Watts-Strogatz', 'Holme-Kim', 'Barabasi-Albert'))
ggsave('Figures/June2019/negBinom_compareNetworks_significance_facet_color.jpg', plot=last_plot())

ggplot(toPlot, aes(x=repression_rate_max, y=value, linetype=Significance, color=graph_type)) + geom_line() + facet_wrap(~sampleSize, ncol=3, dir='h') + theme_bw() +  xlab('Repression Rate') + ylab('Percent of Regressions') + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_color_grey(limits=c('Watts-Strogatz', 'Holme-Kim', 'Barabasi-Albert')) + labs(color='Graph Type')
ggsave('Figures/June2019/negBinom_compareNetworks_significance_facet_grey.jpg', plot=last_plot())



###################################
### Regardless of network scaling, size changes
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
# toPlot_sf <- read.csv('Data/NetworkSimulation/Summer2019/processedData/3wayCompare_scaleFree_sizeChange_regressionData_toplot.csv', stringsAsFactors=FALSE)
# toPlot_sworld <- read.csv('Data/NetworkSimulation/Summer2019/processedData/3wayCompare_smallWorld_sizeChange_regressionData_toplot.csv', stringsAsFactors=FALSE)
# toPlot_hk <- read.csv('Data/NetworkSimulation/Summer2019/processedData/3wayCompare_holmeKim_sizeChange_regressionData_toplot.csv', stringsAsFactors=FALSE)


toPlot_sf$graph_type <- 'Barabasi-Albert'
toPlot_sworld$graph_type <- 'Watts-Strogatz'
toPlot_hk$graph_type <- 'Holme-Kim'
# ###
# ### RUN ABOVE IF NEED TO LOAD DATA LATER, LIKE WHEN MAKING NARROW XLIM
# ###

# combine
toPlot <- rbind(toPlot_sf, toPlot_sworld, toPlot_hk)

### Plot, all
ggplot(toPlot, aes(x=repression_rate_max, y=value, linetype=Significance, color=graph_type)) + geom_line() + facet_wrap(~sampleSize, ncol=3, dir='h') + theme_bw() +  xlab('Repression Rate') + ylab('Percent of Regressions') + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(color='Graph Type') + scale_color_discrete(limits=c('Watts-Strogatz', 'Holme-Kim', 'Barabasi-Albert'))
ggsave('Figures/June2019/negBinom_compareNetworks_significance_facet_color_sizeChanged.jpg', plot=last_plot())

ggplot(toPlot, aes(x=repression_rate_max, y=value, linetype=Significance, color=graph_type)) + geom_line() + facet_wrap(~sampleSize, ncol=3, dir='h') + theme_bw() +  xlab('Repression Rate') + ylab('Percent of Regressions') + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_color_grey(limits=c('Watts-Strogatz', 'Holme-Kim', 'Barabasi-Albert')) + labs(color='Graph Type')
ggsave('Figures/June2019/negBinom_compareNetworks_significance_facet_grey_sizeChanged.jpg', plot=last_plot())





###################################
### AS NETWORK STRUCTURE CONTROL CHANGES
###################################
rates <- sort(unique(sf$networkStructureControl))
n <- c(100, 250, 500, 1000, 1500, 2000)

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

# toPlot_sf <- read.csv('Data/NetworkSimulation/Summer2019/processedData/3wayCompare_scaleFree_regressionData_networkStructureControl_toplot.csv', stringsAsFactors=FALSE)
# toPlot_sworld <- read.csv('Data/NetworkSimulation/Summer2019/processedData/3wayCompare_smallWorld_regressionData_networkStructureControl_toplot.csv', stringsAsFactors=FALSE)
# toPlot_hk <- read.csv('Data/NetworkSimulation/Summer2019/processedData/3wayCompare_holmeKim_regressionData_networkStructureControl_toplot.csv', stringsAsFactors=FALSE)


## Standardize rates so can be on own plot
toPlot_sf$scaling_parameter_max <- toPlot_sf$scaling_parameter_max - 2  # -2 because range for the others is [0-1], this is [2-3]


toPlot <- rbind(toPlot_sf, toPlot_sworld, toPlot_hk)

### Plot, all, scaling parameter
ggplot(toPlot, aes(x=scaling_parameter_max, y=value, linetype=Significance, color=graph_type)) + geom_line() + facet_wrap(~sampleSize, ncol=3, dir='h') + theme_bw() +  xlab('Network Structure Control') + ylab('Percent of Regressions') + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(color='Graph Type') + scale_color_discrete(limits=c('Watts-Strogatz', 'Holme-Kim', 'Barabasi-Albert'))
ggsave('Figures/June2019/negBinom_compareNetworks_significance_facet_color_networkStructureControl.jpg', plot=last_plot())


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
results_sf <- makeResults(data=sf2, repressionRates=rates, sampleSizes=n, trials=100, regressionModel=ma, scalingBool=FALSE, repressionBool=TRUE)
results_sworld <- makeResults(data=sworlsf2, repressionRates=rates, sampleSizes=n, trials=100, regressionModel=ma, scalingBool=FALSE, repressionBool=TRUE)
results_hk <- makeResults(data=hk2, repressionRates=rates, sampleSizes=n, trials=100, regressionModel=ma, scalingBool=FALSE, repressionBool=TRUE)

# get data to plot
toPlot_sf <- processResults(data=results_sf, tscore=1.96, networkType='3wayCompare_scaleFree_realisticValue')
toPlot_sworld <- processResults(data=results_sworld, tscore=1.96, networkType='3wayCompare_smallWorld_realisticValue')
toPlot_hk <- processResults(data=results_hk, tscore=1.96, networkType='3wayCompare_holmeKim_realisticValue')

toPlot_sf$graph_type <- 'Barabasi-Albert'
toPlot_sworld$graph_type <- 'Watts-Strogatz'
toPlot_hk$graph_type <- 'Holme-Kim'

# combine
toPlot <- rbind(toPlot_sf, toPlot_sworld, toPlot_hk)

### RUN BELOW IF NEED TO LOAD DATA LATER, LIKE WHEN MAKING NARROW XLIM
# toPlot_sf <- read.csv('Data/NetworkSimulation/Summer2019/processedData/3wayCompare_scaleFree_realisticValue_regressionData_toplot.csv', stringsAsFactors=FALSE)
# toPlot_sworld <- read.csv('Data/NetworkSimulation/Summer2019/processedData/3wayCompare_smallWorld_realisticValue_regressionData_toplot.csv', stringsAsFactors=FALSE)
# toPlot_hk <- read.csv('Data/NetworkSimulation/Summer2019/processedData/3wayCompare_holmeKim_realisticValue_regressionData_toplot.csv', stringsAsFactors=FALSE)

toPlot_sf$graph_type <- 'Barabasi-Albert'
toPlot_sworld$graph_type <- 'Watts-Strogatz'
toPlot_hk$graph_type <- 'Holme-Kim'

# combine
toPlot <- rbind(toPlot_sf, toPlot_sworld, toPlot_hk)

### Plot, all
ggplot(toPlot, aes(x=repression_rate_max, y=value, linetype=Significance, color=graph_type)) + geom_line() + facet_wrap(~sampleSize, ncol=3, dir='h') + theme_bw() +  xlab('Repression Rate') + ylab('Percent of Regressions') + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(color='Graph Type') + scale_color_discrete(limits=c('Watts-Strogatz', 'Holme-Kim', 'Barabasi-Albert'))
ggsave('Figures/June2019/negBinom_compareNetworks_significance_realisticValue_facet_color.jpg', plot=last_plot())

ggplot(toPlot, aes(x=repression_rate_max, y=value, linetype=Significance, color=graph_type)) + geom_line() + facet_wrap(~sampleSize, ncol=3, dir='h') + theme_bw() +  xlab('Repression Rate') + ylab('Percent of Regressions') + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_color_grey(limits=c('Watts-Strogatz', 'Holme-Kim', 'Barabasi-Albert')) + labs(color='Graph Type')
ggsave('Figures/June2019/negBinom_compareNetworks_significance_realisticValue_facet_grey.jpg', plot=last_plot())




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
ggplot(toPlot, aes(x=repression_rate_max, y=value, linetype=Significance, color=graph_type)) + geom_line() + facet_wrap(~sampleSize, ncol=3, dir='h') + theme_bw() +  xlab('Repression Rate') + ylab('Percent of Regressions') + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(color='Graph Type') + scale_color_discrete(limits=c('Watts-Strogatz', 'Holme-Kim', 'Barabasi-Albert'))
ggsave('Figures/June2019/negBinom_compareNetworks_significance_facet_mostSimilar_color.jpg', plot=last_plot())

ggplot(toPlot, aes(x=repression_rate_max, y=value, linetype=Significance, color=graph_type)) + geom_line() + facet_wrap(~sampleSize, ncol=3, dir='h') + theme_bw() +  xlab('Repression Rate') + ylab('Percent of Regressions') + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_color_grey(limits=c('Watts-Strogatz', 'Holme-Kim', 'Barabasi-Albert')) + labs(color='Graph Type')
ggsave('Figures/June2019/negBinom_compareNetworks_significance_facet_mostSimilar_grey.jpg', plot=last_plot())



