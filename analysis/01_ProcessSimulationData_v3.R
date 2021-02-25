'''
This script is to process the network data.

Designed to work on R 3.5.2

v2.  1. Trying new measures of core.
	 2. Change file out name to reflect new core
	 	A. based on mean eigen centrality
	 	B. based on core according to mean eigen and median threshold
	 3. Change repressionWorked to account for final_size
v3.  1. Use data with new repression type for R&R

'''
set.seed(480193)

library(dplyr)

setwd('/Users/Zack/Documents/UCLA/Research/ProtestTheory')


processRawData <- function(file_path){
	data <- read.csv(file_path, stringsAsFactors=FALSE)
	data$param_idx <- sprintf("%03d", data$param_idx)  # Doing this makes the numbers nicer.  None greater than 999
	data$trial_idx <- sprintf("%03d", data$trial_idx)  # Doing this makes the numbers nicer.  None greater than 999
	data$key <- paste(data$graph_type, data$param_idx, data$trial_idx, sep='_')

	data$active_percent <- data$active_nodes/data$total_nodes

	data$growth_percent <- (data$active_nodes-data$initial_size)/data$initial_size

	return(data)
}

aggRawData <- function(data){
	endTrials <- data.frame(data %>% group_by(key) %>% slice(n()))
	startTrials <- data.frame(data %>% group_by(key) %>% slice(1))

	endTrials <- merge(endTrials, startTrials[,c('key', 'communities_with_protesters', 'mean_community_protest_percent')], by='key')

	names(endTrials)[names(endTrials)=='communities_with_protesters.x'] <- 'ending_communities_with_protesters'
	names(endTrials)[names(endTrials)=='communities_with_protesters.y'] <- 'starting_communities_with_protesters'
	names(endTrials)[names(endTrials)=='mean_community_protest_percent.y'] <- 'starting_mean_community_protest_percent'
	names(endTrials)[names(endTrials)=='mean_community_protest_percent.x'] <- 'ending_mean_community_protest_percent'

	return(endTrials)

}


assignCoreAndSuccess <- function(data){
	data$successDouble <- ifelse(data$growth_percent >= 2, 1, 0)
	data$successGreater35 <- ifelse(data$active_nodes >= 35, 1, 0)
	#core_thresh <- quantile(data$initial_mean_degree, .95)
	#data$coreStart <- ifelse(data$initial_mean_degree >= core_thresh, 1, 0)

	core_degree_thresh <- quantile(data$initial_median_degree, .95)
	data$coreDegree <- ifelse(data$initial_median_degree >= core_degree_thresh, 1, 0)

	eigen_threshold <- quantile(data$initial_mean_eigen_centrality, .95)
	data$coreEigen <- ifelse(data$initial_mean_eigen_centrality >= eigen_threshold, 1, 0)

	data$coreStart <- data$coreEigen
	data$coreStart2 <- data$coreEigen*data$coreDegree

	data$repressionWorked <- ifelse(data$active_nodes < data$initial_size & data$active_nodes < 35, 1, 0)

	return(data)
}


# selectTrialParameters <- function(data, column, value){
# 	temp <- data[data[[column]] == value,]

# 	return(temp)
# }

########
##
## PROCESS EXPERIMENTS' RESULTS, MERGE, SAVE OUT
##
########

#### HOLME-KIM
hk1 <- processRawData(file_path='Data/NetworkSimulation/Summer2020/exp3-1-graph_type=powerlaw_cluster_graph-repression_type=node_removal-threshold_type=uniform-num_nodes=1000-m=3.csv')
hk12 <- aggRawData(data=hk1)
hk12 <- assignCoreAndSuccess(data=hk12)

hk2 <- processRawData(file_path='Data/NetworkSimulation/Summer2020/exp3-2-graph_type=powerlaw_cluster_graph-repression_type=node_removal-threshold_type=uniform-num_nodes=1000-m=3.csv')
hk22 <- aggRawData(data=hk2)
hk22 <- assignCoreAndSuccess(data=hk22)

hk <- rbind(hk12, hk22)

#### WATTS-STROGATZ
ws1 <- processRawData(file_path='Data/NetworkSimulation/Summer2020/exp2-1-graph_type=watts_strogatz_graph-repression_type=node_removal-threshold_type=uniform-num_nodes=1000-k=15.csv')
ws12 <- aggRawData(data=ws1)
ws12 <- assignCoreAndSuccess(data=ws12)

ws2 <- processRawData(file_path='Data/NetworkSimulation/Summer2020/exp2-2-graph_type=watts_strogatz_graph-repression_type=node_removal-threshold_type=uniform-num_nodes=1000-k=15.csv')
ws22 <- aggRawData(data=ws2)
ws22 <- assignCoreAndSuccess(data=ws22)

ws <- rbind(ws12, ws22)

#### SCALE-FREE
## HAD TO SPLIT INTO FOUR RUNS, SO NEED TO LOAD THOSE AND THEN RBIND.
sf1 <- processRawData(file_path='Data/NetworkSimulation/Summer2020/exp1-1-graph_type=scale_free_graph-repression_type=node_removal-threshold_type=uniform-num_nodes=1000.csv')
sf12 <- aggRawData(data=sf1)
sf12 <- assignCoreAndSuccess(data=sf12)

sf2 <- processRawData(file_path='Data/NetworkSimulation/Summer2020/exp1-2-graph_type=scale_free_graph-repression_type=node_removal-threshold_type=uniform-num_nodes=1000.csv')
sf22 <- aggRawData(data=sf2)
sf22 <- assignCoreAndSuccess(data=sf22)

sf3 <- processRawData(file_path='Data/NetworkSimulation/Summer2020/exp1-3-graph_type=scale_free_graph-repression_type=node_removal-threshold_type=uniform-num_nodes=1000.csv')
sf32 <- aggRawData(data=sf3)
sf32 <- assignCoreAndSuccess(data=sf32)


sf4 <- processRawData(file_path='Data/NetworkSimulation/Summer2020/exp1-4-graph_type=scale_free_graph-repression_type=node_removal-threshold_type=uniform-num_nodes=1000.csv')
sf42 <- aggRawData(data=sf4)
sf42 <- assignCoreAndSuccess(data=sf42)


sf <- rbind(sf12, sf22, sf32, sf42)




########### COMBINE
### Make column names the name
names(ws)[(names(ws) %in% names(hk)==FALSE)]  # k.  all are 15, is degree of all starting nodes.
names(hk)[(names(hk) %in% names(ws)==FALSE)]  # m.  all are 3, is number edges to add per new node

hk$k <- NA
ws$m <- NA

names(ws)[(names(ws) %in% names(sf)==FALSE)]  # p, k, m.  Will be same for hk2 since has same columns as ws2
names(sf)[(names(sf) %in% names(ws)==FALSE)]  # scaling_parameter

hk$scaling_parameter <- NA
ws$scaling_parameter <- NA

sf$p <- NA
sf$k <- NA
sf$m <- NA

# Since p in ws2 and hk2 and scaling_parameter in sf are what control network structure, let's make a new variable with an easy to remember name.
hk$networkStructureControl <- hk$p
ws$networkStructureControl <- ws$p
sf$networkStructureControl <- sf$scaling_parameter

all <- rbind(hk, ws, sf)

# Result will be 797mb when saved, so let's trim.
# p and scaling_parameter vary, but are now captured in networkStructureControl variable
all <- all[,!(names(all) %in% c('X', 'p', 'k', 'm', 'scaling_parameter', 'repression_type'))]


########### SAVE
write.csv(all, 'Data/NetworkSimulation/Summer2020/processedData/01_experimentsCombined_eigenCore_newRepression.csv', row.names=FALSE)



