'''
This script processes the edge removal data and makes a figure for the main paper.
Done in response to R1 from Network Science.
It combines 01_ProcessSimulationData with 08_RegressionResultsOvertime
'''
set.seed(480193)

library(dplyr)

setwd('/Users/Zack/Documents/UCLA/Research/ProtestTheory')


processRawData <- function(file_path){
	data <- read.csv(file_path, stringsAsFactors=FALSE)
	# data$param_idx <- sprintf("%03d", data$param_idx)  # Doing this makes the numbers nicer.  None greater than 999
	# data$trial_idx <- sprintf("%03d", data$trial_idx)  # Doing this makes the numbers nicer.  None greater than 999
	# data$key <- paste(data$graph_type, data$param_idx, data$trial_idx, sep='_')

	data$active_percent <- data$initial_size/data$final_size

	data$growth_percent <- (data$final_size-data$initial_size)/data$initial_size

	return(data)
}


### The edge repression data are so old that they do not have step by step, so no need to aggregated.
# aggRawData <- function(data){
# 	endTrials <- data.frame(data %>% group_by(key) %>% slice(n()))
# 	startTrials <- data.frame(data %>% group_by(key) %>% slice(1))

# 	endTrials <- merge(endTrials, startTrials[,c('key', 'communities_with_protesters', 'mean_community_protest_percent')], by='key')

# 	names(endTrials)[names(endTrials)=='communities_with_protesters.x'] <- 'ending_communities_with_protesters'
# 	names(endTrials)[names(endTrials)=='communities_with_protesters.y'] <- 'starting_communities_with_protesters'
# 	names(endTrials)[names(endTrials)=='mean_community_protest_percent.y'] <- 'starting_mean_community_protest_percent'
# 	names(endTrials)[names(endTrials)=='mean_community_protest_percent.x'] <- 'ending_mean_community_protest_percent'

# 	return(endTrials)

# }


assignCoreAndSuccess <- function(data){
	data$successDouble <- ifelse(data$growth_percent >= 2, 1, 0)
	data$successGreater35 <- ifelse(data$final_size >= 35, 1, 0)
	#core_thresh <- quantile(data$initial_mean_degree, .95)
	#data$coreStart <- ifelse(data$initial_mean_degree >= core_thresh, 1, 0)

	core_degree_thresh <- quantile(data$initial_median_degree, .95)
	data$coreDegree <- ifelse(data$initial_median_degree >= core_degree_thresh, 1, 0)

	eigen_threshold <- quantile(data$initial_mean_eigen_centrality, .95)
#	data$coreEigen <- ifelse(data$initial_mean_eigen_centrality >= eigen_threshold, 1, 0)  # Did not look at this that early, so comment out

	# data$coreStart <- data$coreEigen
	# data$coreStart2 <- data$coreEigen*data$coreDegree

	data$repressionWorked <- ifelse(data$final_size < data$initial_size & data$final_size < 35, 1, 0)

	return(data)
}

########
##
## PROCESS EXPERIMENTS' RESULTS, MERGE, SAVE OUT
##	# Subset by repression rate because original models went all the way to 1, not .15 like later.
########

#### HOLME-KIM
hk <- processRawData(file_path='Data/NetworkSimulation/Summer2018/exp10-threshold_type=uniform-repression_type=edge_removal-m=3-graph_type=powerlaw_cluster_graph-num_nodes=3000.csv')
hk2 <- assignCoreAndSuccess(data=hk)
hk2 <- subset(hk2, repression_rate <= .1501)


#### WATTS-STROGATZ
ws <- processRawData(file_path='Data/NetworkSimulation/Summer2018/exp14-threshold_type=uniform-num_nodes=1000-repression_type=edge_removal-graph_type=watts_strogatz_graph-k=50.csv')
ws2 <- assignCoreAndSuccess(data=ws)
ws2 <- subset(ws2, repression_rate <= .1501)

#### SCALE-FREE
sf <- processRawData(file_path='Data/NetworkSimulation/Summer2018/exp16-threshold_type=uniform-num_nodes=1000-repression_type=edge_removal-graph_type=scale_free_graph.csv')
sf2 <- assignCoreAndSuccess(data=sf)
sf2 <- subset(sf2, repression_rate <= .1501)

hk2$k <- NA
ws2$m <- NA

hk2$scaling_parameter <- NA
ws2$scaling_parameter <- NA

sf2$p <- NA
sf2$k <- NA
sf2$m <- NA

# Since p in ws2 and hk2 and scaling_parameter in sf are what control network structure, let's make a new variable with an easy to remember name.
hk2$networkStructureControl <- hk2$p
ws2$networkStructureControl <- ws2$p
sf2$networkStructureControl <- sf2$scaling_parameter

all <- rbind(hk2, ws2, sf2)
all <- all[,!(names(all) %in% c('X', 'p', 'k', 'm', 'scaling_parameter', 'repression_type'))]
write.csv(all, 'Data/NetworkSimulation/Summer2018/processedData/01_experimentsCombined_edgeRepression.csv', row.names=FALSE)

