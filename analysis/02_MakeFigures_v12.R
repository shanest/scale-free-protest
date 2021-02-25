'''
This script is to make the large faceted figures that may be the main result to show.

Designed to work on R 3.5.2

v1. Get the basic structure working.  Get up to faceting by outcome and distinguishing by core start.  So save and more to v2.
v2. 1. Add logic to work on all network types at once.
v3. 1. Got the logic working for code to repeat, so will now put in function.  The key was learning how to pass strings as object names to ggplot.  Here: https://stackoverflow.com/questions/46372691/ggplot-and-dplyr-and-column-name-as-string.
	2. Add second funtion for naming outfile
v4. 1. Experimenting with changing y-axis range by facet, gets complicated so starting new version.  See: https://stackoverflow.com/questions/51735481/ggplot2-change-axis-limits-for-each-individual-facet-panel.
v5. Got the code working, so now clean up the scratch work.
v6. Basically all graphs are just black lines; only ones that are not are what I tested code on.  Need to find out if that is bug in code or actual results.
	1. The bug was that I did not round values, so things like Eigenvector centrality or clusteirng would have 10, 15 significant figures and therefore be their own groups.  Am now rounding everything to 2 decimal places.
v7. Making more charts with different subsets.  Want to update makePlots to not aggregate to core success based on argument to function.
	1. Add + .00001 to value in results_melted so that do not get error messages about infinite values when log.
v8. 1. Use new core measure based on Eigenvector centrality
v9.	1. Change outfile name so group by variable looking at, makes easier to compare across data subsets.
	2. Round lots of measures so that grouping is clearer
	3. Comnent out line in function that does that
v10. CHANGE FILE NAME LATE AT NIGHT.  THIS WILL BE TAKE UP WORK FROM SIEGEL 2018, SEE MOLESKINE.  MAIN THING IS TO MAKE DV ACTIVE_PERCENT, THEN SLOWLY BUILD MODELS UP BY EXPLORING DIFFERENT PARTS OF PARAMETER SPACE.  SO, LOTS OF DIFFERENT GROUPINGS AND SUBSETS, MAY REQUIRE TOTALLY NEW FUNCTION OR AT LEAST HEAVY REWORK.

	ONE REASON WILL NEED HEAVY REWORK IS BECAUSE WANT TO SEPARATE BY GRAPH TYPE BEFORE STARTING ANALYSIS.  THAT IS, GRAPH TYPE WILL NOT BE A VARIABLE, WILL BE DIFFERENT SUBSETS OF DATA ON WHICH TO WORK.
	0. Get rid of "Repression Works" plotting - DONE
	1. Make active_percent the outcome variable - DONE
	2. Separate out the 3 graph types - DONE
	3. Show change in active_percent as initial_mean_eigen_centrality increases - DONE
	4. Then, introduce change in networkStructureControl (clustering for hk, sw; skew for scale free)
	5. Then, introduce repression_rate

	On morning of 05.28.2019, full model will be active_percent = mean_eigen + network structure + (repression_rate OR another network measure).  Need to make initial subsets based on that final model.
v11. v10 got very messy.  Have settled on m0 is seed_eigen_centrality, then m1 adds initial_median_degree.  Clean up the rest from v10 that is not that.
v12. IS MAJOR UPDATE, FIRST UPDATE AFTER POLNET.  SEE PAGE 121 OF MOLESKINE NOTEBOOK.
	M0 - mean eigenvector centrality
			- scale free parameter is 2.3
			- holme kim parameter is .315
			- watts strogatz parameter is .4

	M1 - network parameter control
	M2 - repression rate

	All done separately by network type.  Because need to identify regions of common behavior per model, cannot really do a function, will have to manually do each section fo each network type.




NB: If m0 is seed_eigen_centrality and m1 is initial_mean_eigen_centrality, then the story is starter needs to be very central but not surrounded by very central people.  That is not necessarily a bad story either...
'''
set.seed(480193)

library(dplyr)
library(ggplot2)
library(reshape)
library(scales)

# Below is for custom y-axis to facet grid.
# See https://stackoverflow.com/questions/51735481/ggplot2-change-axis-limits-for-each-individual-facet-panel
#library(devtools)
#devtools::install_github("zeehio/facetscales")
library(facetscales)



setwd('/Users/Zack/Documents/UCLA/Research/ProtestTheory')



########
##
## DATA
##
########
all <- read.csv('Data/NetworkSimulation/Spring2019/processedData/01_experimentsCombined_eigenCore.csv', stringsAsFactors=FALSE)


#temp <- all

# Round to 2 decimals
toRound <- c('ending_mean_commmunity_protest_percent', 'growth_percent', 'initial_density', 'initial_global_clustering', 'initial_median_degree', 'initial_mean_eigen_centrality', 'initial_neighborhood_clustering', 'initial_nodes_clustering', 'seed_eigen_centrality', 'starting_mean_commmunity_protest_percent', 'active_percent')
for(i in 1:length(toRound)){
	all[toRound[i]] <- round(all[toRound[i]], 2)
}

# Should be whole number
toRound <- c('initial_median_degree', 'initial_mean_degree')
for(i in 1:length(toRound)){
	all[toRound[i]] <- round(all[toRound[i]])
}

########
##
## FUNCTIONS
##
########
# data = what loaded initially.
# outfile = the variable name to plot as the x-axis, will form part of filename too.
# outfileMore = descriptor to add
# xlabPretty = Make x-axis label human readable
# groupByCore = Whether or not to distinguish between core and periphery, matters because some x-axis variables are essentially the core anyway
makePlots <- function(data, outfile, outfileMore=NULL, xlabPretty=NULL, groupByCore=TRUE){
	data[[outfile]] <- round(data[[outfile]], 2)

	if(groupByCore==TRUE){
		#results <- data.frame(data %>% group_by(graph_type, UQ(as.name(outfile)), coreStart) %>% summarize(successDouble=mean(successDouble), successGreater35=mean(successGreater35), repressionWorked=mean(repressionWorked), active_percent=mean(active_percent, na.rm=TRUE)))
		results <- data.frame(data %>% group_by(graph_type, UQ(as.name(outfile)), coreStart) %>% summarize(successDouble=mean(successDouble), successGreater35=mean(successGreater35), active_percent=mean(active_percent, na.rm=TRUE)))
		results_melted <- melt(results, id.vars= c('graph_type', outfile, 'coreStart'))
		results_melted$coreStart <- ifelse(results_melted$coreStart==1, "Yes", "No")

	}

	if(groupByCore==FALSE){
		#results <- data.frame(data %>% group_by(graph_type, UQ(as.name(outfile))) %>% summarize(successDouble=mean(successDouble), successGreater35=mean(successGreater35), repressionWorked=mean(repressionWorked), active_percent=mean(active_percent, na.rm=TRUE)))
		results <- data.frame(data %>% group_by(graph_type, UQ(as.name(outfile))) %>% summarize(successDouble=mean(successDouble), successGreater35=mean(successGreater35), active_percent=mean(active_percent, na.rm=TRUE)))
		results_melted <- melt(results, id.vars= c('graph_type', outfile))

		outfileMore <- paste0('noCoreDistinction_', outfileMore)

	}

	# Make names pretty
	results_melted$graph_type[results_melted$graph_type=='powerlaw_cluster_graph'] <- 'Holme-Kim'
	results_melted$graph_type[results_melted$graph_type=='watts_strogatz_graph'] <- 'Small World'
	results_melted$graph_type[results_melted$graph_type=='scale_free_graph'] <- 'Scale Free'

	results_melted$variable <- as.character(results_melted$variable)
	results_melted$variable[results_melted$variable=='successDouble'] <- 'Protest Doubles'
	results_melted$variable[results_melted$variable=='successGreater35'] <- 'Final Protest >= 35'
	#results_melted$variable[results_melted$variable=='repressionWorked'] <- 'Repression Works'
	results_melted$variable[results_melted$variable=='active_percent'] <- 'Percent Mobilized'



	# Make custom y scale
	minPD <- min(results_melted$value[results_melted$variable=='Protest Doubles'])
	minPD <- ifelse(minPD==0, .0000001, minPD)
	maxPD <- max(results_melted$value[results_melted$variable=='Protest Doubles'])+.1

	minFP <- min(results_melted$value[results_melted$variable=='Final Protest >= 35'])
	minFP <- ifelse(minFP==0, .0000001, minFP)
	maxFP <- max(results_melted$value[results_melted$variable=='Final Protest >= 35'])+.1

	# minRW <- min(results_melted$value[results_melted$variable=='Repression Works'])
	# minRW <- ifelse(minRW==0, .0000001, minRW)
	# maxRW <- max(results_melted$value[results_melted$variable=='Repression Works'])  # Don't need +.1 because repression rate is always high

	minPM <- min(results_melted$value[results_melted$variable=='Percent Mobilized'])
	minPM <- ifelse(minPM==0, .0000001, minPM)
	maxPM <- max(results_melted$value[results_melted$variable=='Percent Mobilized'])

	results_melted$value <- results_melted$value + .0000001

	scales_y <- list(
			'Final Protest >= 35'=scale_y_log10(limits=c(minFP,maxFP)),
			'Protest Doubles' = scale_y_log10(limits=c(minPD, maxPD)),
			#'Repression Works' = scale_y_log10(limits=c(minRW, maxRW)),
			'Percent Mobilized' = scale_y_log10(limits=c(minPM, maxPM))
			#'Repression Works' = scale_y_log10(limits=c(0.45,1), breaks=seq(.45, 1, .1))
			)

	if(groupByCore==TRUE){
		thePlot <- ggplot(results_melted, aes(x=UQ(as.name(outfile)), y=value, linetype=factor(coreStart))) + geom_line() + facet_grid_sc(variable ~ graph_type, scales=list(y=scales_y, x='free')) + theme_bw() + xlab(xlabPretty) + ylab("") + scale_y_continuous(trans=log10_trans(), breaks=trans_breaks('log10', function(x) 10^x), labels=trans_format('log10', math_format(10^.x))) + scale_linetype_discrete(name='Core Start')
	}

	if(groupByCore==FALSE){
		thePlot <- ggplot(results_melted, aes(x=UQ(as.name(outfile)), y=value)) + geom_line() + facet_grid_sc(variable ~ graph_type, scales=list(y=scales_y, x='free')) + theme_bw() + xlab(xlabPretty) + ylab("") + scale_y_continuous(trans=log10_trans(), breaks=trans_breaks('log10', function(x) 10^x), labels=trans_format('log10', math_format(10^.x)))
	}

	# Below is trying to facet with more variables
	jpeg(paste0('Figures/June2019/', outfile, '_eigenCoreMeasure_', outfileMore, '.jpg'), res=300, width=8, height=8, units='in')
	print(thePlot)
	dev.off()
}


############################################################
############################################################
#
#
#		FOLLOW METHODS IN SEIGEL 2018, ANALYZING COMPUTATIONAL MODELS
#			NB: I AM NOT SURE HOW YOU COULD AUTOMATE THE BELOW SINCE IT REQUIRES VISUAL INSPECTION TO FIND DIFFERENT AREAS OF BEHAVIOR.  YOU COULD WITH SOME PROGRAMMED RULES, BUT NOT WORTH IT FOR THIS PROJECT.
#			NB: MEANS WILL HAVE TO MANUALLY DO THE SAME FOR HOLME-KIM, SMALL-WORLD NETWORKS.
#			NB: COMPLETED SCALE-FREE AT 1:30 P.M. ON 05.28.2019, WILL FOCUS ON POSTER MAKING NOW.
#
############################################################
############################################################

## FIRST, WILL KEEP THOSE WHERE INITIAL SIZE LEAST 2
## DECIDE NOT TO DO BECAUSE ADDS A RANDOM SPIKE AT HIGH LEVELS OF EIGENVECTOR CENTRALITY THAT I THINK IS JUST NOISE
# all2 <- all[all$initial_size >= 3,]

# ## SECOND, WILL SPLIT BY GRAPH TYPE
# sw <- all2[all2$graph_type=='watts_strogatz_graph',]
# hk <- all2[all2$graph_type=='powerlaw_cluster_graph',]
# sf <- all2[all2$graph_type=='scale_free_graph',]


## SECOND, WILL SPLIT BY GRAPH TYPE
sw <- all[all$graph_type=='watts_strogatz_graph',]
sw$networkStructureControl <- round(sw$networkStructureControl, 2)

hk <- all[all$graph_type=='powerlaw_cluster_graph',]
sf <- all[all$graph_type=='scale_free_graph',]



###########
###########
# SCALE_FREE
###########
###########
sf_m0 <- subset(sf, repression_rate==0 & networkStructureControl==2.3)

#### M0: eigen centrality
#m0 <- sf[sf$networkStructureControl==2.3 & sf$repression_rate==0,]
m0 <- data.frame(sf_m0 %>% group_by(initial_mean_eigen_centrality) %>% summarize(active_percent=mean(active_percent)))
jpeg('Figures/June2019/Siegel2018_scale_free_a_m0_meanEigenCentrality.jpeg', width=5, height=5, res=300, units='in')
ggplot(m0, aes(x=initial_mean_eigen_centrality, y=active_percent)) + geom_line() + xlab('Mean Eigenvector Centrality of First Protesters') + ylab('Avg. Percent Mobilization of Network') + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


#### M1: M0 + network control parameter
# Visual inspect of M0 shows an increasing and then decreasing zone
sf_m1 <- subset(sf, repression_rate==0)
zone1 <- m0$initial_mean_eigen_centrality[which.max(m0$active_percent)]
m1 <- data.frame(sf_m1 %>% group_by(initial_mean_eigen_centrality, networkStructureControl) %>% summarize(active_percent=mean(active_percent)))

## The two regions of m0
to_keep <- seq(from=0, to=zone1, length.out=4)  # Not from 0 because get weird behavior there
to_keep <- round(to_keep, 2)

to_keep2 <- seq(from=zone1+.01, to=max(m1$initial_mean_eigen_centrality), length.out=4)
to_keep2 <- round(to_keep2, 2)

to_keep <- c(to_keep, to_keep2)

# Will look at just a few points from m0
m1 <- m1[which(m1$initial_mean_eigen_centrality %in% to_keep),]

# Prepare for plotting
m1$facets <- paste0('Mean Eigen Centrality: ', m1$initial_mean_eigen_centrality)  # Do this line to get pretty labels for facets

jpeg('Figures/June2019/Siegel2018_scale_free_a_m1_networkStructureControl.jpeg', width=6, height=6, res=300, units='in')
ggplot(m1, aes(x=networkStructureControl, y=active_percent)) + geom_line() + facet_wrap(~facets, ncol=3, dir='v') + xlab('Scaling Parameter') + ylab('Avg. Percent Mobilization of Network') + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


##### Visual inspection of above shows 3 zones: increasing (0, .2), decreasing, (.05, .07, .08), no change (.18, .28, .3).

#### M2: M1 + repression_rate
increasing <- subset(sf, initial_mean_eigen_centrality==.02 & networkStructureControl==2.3)
decreasing <- subset(sf, initial_mean_eigen_centrality==.07 & networkStructureControl==2.3)
constant <- subset(sf, initial_mean_eigen_centrality==.18 & networkStructureControl==2.3)

increasing <- data.frame(increasing %>% group_by(repression_rate) %>% summarize(active_percent=mean(active_percent), initial_mean_eigen_centrality=mean(initial_mean_eigen_centrality)))
decreasing <- data.frame(decreasing %>% group_by(repression_rate) %>% summarize(active_percent=mean(active_percent), initial_mean_eigen_centrality=mean(initial_mean_eigen_centrality)))
constant <- data.frame(constant %>% group_by(repression_rate) %>% summarize(active_percent=mean(active_percent), initial_mean_eigen_centrality=mean(initial_mean_eigen_centrality)))

increasing$facets <- paste0('Mean Eigen Centrality: ', increasing$initial_mean_eigen_centrality, '\n', 'Scaling Parameter = 2.3')
decreasing$facets <- paste0('Mean Eigen Centrality: ', decreasing$initial_mean_eigen_centrality, '\n', 'Scaling Parameter = 2.3')
constant$facets <- paste0('Mean Eigen Centrality: ', constant$initial_mean_eigen_centrality, '\n', 'Scaling Parameter = 2.3')

# Plot first zone, increasing
ggplot(increasing, aes(x=repression_rate, y=active_percent)) + geom_line() + xlab('Repression Rate') + ylab('Avg. Percent Mobilization of Network') + theme_bw() + facet_wrap(~facets) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave('Figures/June2019/Siegel2018_scale_free_a_m2_repressionRate_zoneIncreasing.jpeg', plot=last_plot())

# Plot second zone, decreasing
ggplot(decreasing, aes(x=repression_rate, y=active_percent)) + geom_line() + xlab('Repression Rate') + ylab('Avg. Percent Mobilization of Network') + theme_bw() + facet_wrap(~facets) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave('Figures/June2019/Siegel2018_scale_free_a_m2_repressionRate_zoneDecreasing.jpeg', plot=last_plot())

# Plot first zone, flat
ggplot(constant, aes(x=repression_rate, y=active_percent)) + geom_line() + xlab('Repression Rate') + ylab('Avg. Percent Mobilization of Network') + theme_bw() + facet_wrap(~facets) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave('Figures/June2019/Siegel2018_scale_free_a_m2_repressionRate_zoneConstant.jpeg', plot=last_plot())




###########
###########
# HOLME-HIM
###########
###########
hk_m0 <- subset(hk, repression_rate==0 & networkStructureControl==.315)

#### M0: eigen centrality
#m0 <- sf[sf$networkStructureControl==2.3 & sf$repression_rate==0,]
m0 <- data.frame(hk_m0 %>% group_by(initial_mean_eigen_centrality) %>% summarize(active_percent=mean(active_percent)))
ggplot(m0, aes(x=initial_mean_eigen_centrality, y=active_percent)) + geom_line() + xlab('Mean Eigenvector Centrality of First Protesters') + ylab('Avg. Percent Mobilization of Network') + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave('Figures/June2019/Siegel2018_holme_kim_a_m0_meanEigenCentrality.jpeg', plot=last_plot())



#### M1: M0 + network control parameter
# Noisy, overall is decreasing
hk_m1 <- subset(hk, repression_rate==0)

zone1 <- m0$initial_mean_eigen_centrality[which.max(m0$initial_mean_eigen_centrality)]
m1 <- data.frame(hk_m1 %>% group_by(initial_mean_eigen_centrality, networkStructureControl) %>% summarize(active_percent=mean(active_percent)))

## The two regions of m0
to_keep <- round(seq(from=0, to=zone1, length.out=9), 2)  # Not from 0 because get weird behavior there

# Will look at just a few points from m0
m1 <- m1[which(m1$initial_mean_eigen_centrality %in% to_keep),]

# Prepare for plotting
m1$facets <- paste0('Mean Eigen Centrality: ', m1$initial_mean_eigen_centrality)  # Do this line to get pretty labels for facets


ggplot(m1, aes(x=networkStructureControl, y=active_percent)) + geom_line() + facet_wrap(~facets, ncol=3, dir='v') + xlab('Triad Formation Rate') + ylab('Avg. Percent Mobilization of Network')   + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave('Figures/June2019/Siegel2018_holme_kim_a_m1_networkStructureControl.jpeg', plot=last_plot())


##### Visual inspection of above shows 2 zones: flat but above 0 (mean eigen centraliy 0-.13), flat and no diffusion (mean eigen centrality .16 and above)

#### M2: M1 + repression_rate
stablePositive <- subset(hk, initial_mean_eigen_centrality==.02 & networkStructureControl==.315)
stableZero <- subset(hk, initial_mean_eigen_centrality==.07 & networkStructureControl==.315)

stablePositive <- data.frame(stablePositive %>% group_by(repression_rate) %>% summarize(active_percent=mean(active_percent), initial_mean_eigen_centrality=mean(initial_mean_eigen_centrality)))
stableZero <- data.frame(stableZero %>% group_by(repression_rate) %>% summarize(active_percent=mean(active_percent), initial_mean_eigen_centrality=mean(initial_mean_eigen_centrality)))

stablePositive$facets <- paste0('Mean Eigen Centrality: ', stablePositive$initial_mean_eigen_centrality, '\n', 'Triad Formation Rate = .315')
stableZero$facets <- paste0('Mean Eigen Centrality: ', stableZero$initial_mean_eigen_centrality, '\n', 'Triad Formation Rate = .315')

# Plot first zone, stablePositive
ggplot(stablePositive, aes(x=repression_rate, y=active_percent)) + geom_line() + xlab('Repression Rate') + ylab('Avg. Percent Mobilization of Network') + theme_bw() + facet_wrap(~facets) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave('Figures/June2019/Siegel2018_holme_kim_a_m2_repressionRate_zoneStablePositive.jpeg', plot=last_plot())

# Plot second zone, stableZero
ggplot(stableZero, aes(x=repression_rate, y=active_percent)) + geom_line() + xlab('Repression Rate') + ylab('Avg. Percent Mobilization of Network') + theme_bw() + facet_wrap(~facets) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave('Figures/June2019/Siegel2018_holme_kim_a_m2_repressionRate_zoneStableZero.jpeg', plot=last_plot())




###########
###########
# SMALL_WORLD, USING MEAN EIGENVECTOR
#  - Very little range in this value, results ugly
###########
###########
sw_m0 <- subset(sw, repression_rate==0 & (networkStructureControl==.48 | networkStructureControl == .34))  # Global clustering of .16 is not well estimated from our networkStructureControl, so this straddles it.

#### M0: eigen centrality
#m0 <- sw[sw$networkStructureControl==.4 & sw$repression_rate==0,]
m0 <- data.frame(sw_m0 %>% group_by(initial_mean_eigen_centrality) %>% summarize(active_percent=mean(active_percent)))
ggplot(m0, aes(x=initial_mean_eigen_centrality, y=active_percent)) + geom_line() + xlab('Mean Eigenvector Centrality of First Protesters') + ylab('Avg. Percent Mobilization of Network') + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave('Figures/June2019/Siegel2018_small_world_a_m0_meanEigenCentrality.jpeg', plot=last_plot())



#### M1: M0 + network control parameter
# Visual inspect of M0 shows only increasing.  Was only 3 points
sw_m1 <- subset(sw, repression_rate==0)
m1 <- data.frame(sw_m1 %>% group_by(initial_mean_eigen_centrality, networkStructureControl) %>% summarize(active_percent=mean(active_percent)))

## The two regions of m0
to_keep <- seq(from=0, to=max(sw$initial_mean_eigen_centrality), by=.01)  # Not from 0 because get weird behavior there

# Will look at just a few points from m0
m1 <- m1[which(m1$initial_mean_eigen_centrality %in% to_keep),]

# Prepare for plotting
m1$facets <- paste0('Mean Eigen Centrality: ', m1$initial_mean_eigen_centrality)  # Do this line to get pretty labels for facets

ggplot(m1, aes(x=networkStructureControl, y=active_percent)) + geom_line() + facet_wrap(~facets, ncol=3, dir='v') + xlab('Rewiring Percent') + ylab('Avg. Percent Mobilization of Network')   + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave('Figures/June2019/Siegel2018_small_world_a_m1_networkStructureControl.jpeg', plot=last_plot())


##### Visual inspection makes it look like no clear pattern, only .02, .03, .04 have enough data

#### M2: M1 + repression_rate
one <- subset(sw, initial_mean_eigen_centrality==.02 & (networkStructureControl==.48 | networkStructureControl == .34))
two <- subset(sw, initial_mean_eigen_centrality==.03 & (networkStructureControl==.48 | networkStructureControl == .34))
three <- subset(sw, initial_mean_eigen_centrality==.04 & (networkStructureControl==.48 | networkStructureControl == .34))

one <- data.frame(one %>% group_by(repression_rate) %>% summarize(active_percent=mean(active_percent), initial_mean_eigen_centrality=mean(initial_mean_eigen_centrality)))
two <- data.frame(two %>% group_by(repression_rate) %>% summarize(active_percent=mean(active_percent), initial_mean_eigen_centrality=mean(initial_mean_eigen_centrality)))
three <- data.frame(three %>% group_by(repression_rate) %>% summarize(active_percent=mean(active_percent), initial_mean_eigen_centrality=mean(initial_mean_eigen_centrality)))

one$facets <- paste0('Mean Eigen Centrality: ', one$initial_mean_eigen_centrality, '\n', 'Scaling Parameter = .4')
two$facets <- paste0('Mean Eigen Centrality: ', two$initial_mean_eigen_centrality, '\n', 'Scaling Parameter = .4')
three$facets <- paste0('Mean Eigen Centrality: ', three$initial_mean_eigen_centrality, '\n', 'Scaling Parameter = .4')



all_sw <- rbind(one, two, three)

# Plot first zone, increasing
ggplot(all_sw, aes(x=repression_rate, y=active_percent)) + geom_line() + xlab('Repression Rate') + ylab('Avg. Percent Mobilization of Network') + theme_bw() + facet_wrap(~facets) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave('Figures/June2019/Siegel2018_small_world_a_m2_repressionRate.jpeg', plot=last_plot())



###########
###########
# SMALL_WORLD, USING MEDIAN DEGREE CENTRALITY
#  - Very little range in this value, results ugly
###########
###########
sw_m0 <- subset(sw, repression_rate==0 & (networkStructureControl==.48 | networkStructureControl == .34))  # Global clustering of .16 is not well estimated from our networkStructureControl, so this straddles it.

#### M0: eigen centrality
#m0 <- sw[sw$networkStructureControl==.4 & sw$repression_rate==0,]
m0 <- data.frame(sw_m0 %>% group_by(initial_median_degree) %>% summarize(active_percent=mean(active_percent)))
ggplot(m0, aes(x=initial_median_degree, y=active_percent)) + geom_line() + xlab('Median Degree of First Protesters') + ylab('Avg. Percent Mobilization of Network') + theme_bw()
ggsave('Figures/June2019/Siegel2018_small_world_b_m0_initialMedianDegree.jpeg', plot=last_plot())


#### M1: M0 + network control parameter
# Visual inspect of M0 shows increasing up to 16, then decreasing
sw_m1 <- subset(sw, repression_rate==0)
zone1 <- m0$initial_median_degree[which.max(m0$active_percent)]
m1 <- data.frame(sw_m1 %>% group_by(initial_median_degree, networkStructureControl) %>% summarize(active_percent=mean(active_percent)))

## The two regions of m0
to_keep <- c(1, 3, 5, 7, 9, 11, 13, 15, 17)

# Will look at just a few points from m0
m1 <- m1[which(m1$initial_median_degree %in% to_keep),]

# Prepare for plotting
m1$facets <- paste0('Initial Median Degree: ', m1$initial_median_degree)  # Do this line to get pretty labels for facets

ggplot(m1, aes(x=networkStructureControl, y=active_percent)) + geom_line() + facet_wrap(~facets, ncol=3, dir='v') + xlab('Rewiring Percent') + ylab('Avg. Percent Mobilization of Network')   + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave('Figures/June2019/Siegel2018_small_world_b_m1_networkStructureControl.jpeg', plot=last_plot())



##### Visual inspection of above shows 2 zones: increasing when mean initial degree is 13, increasing then flat for mean degree of 15

#### M2: M1 + repression_rate
increasing <- subset(sw, initial_median_degree==13)
el_increasing <- subset(sw, initial_median_degree==15 & networkStructureControl==.08)
el_flat <- subset(sw, initial_median_degree==15 & networkStructureControl==.34)

increasing <- data.frame(increasing %>% group_by(repression_rate) %>% summarize(active_percent=mean(active_percent), initial_median_degree=mean(initial_median_degree)))
el_increasing <- data.frame(el_increasing %>% group_by(repression_rate) %>% summarize(active_percent=mean(active_percent), initial_median_degree=mean(initial_median_degree)))
el_flat <- data.frame(el_flat %>% group_by(repression_rate) %>% summarize(active_percent=mean(active_percent), initial_median_degree=mean(initial_median_degree)))

increasing$facets <- paste0('Median Degree of First Protesters: ', increasing$initial_median_degree, '\n', 'Rewiring Percent = .02')
el_increasing$facets <- paste0('Median Degree of First Protesters: ', el_increasing$initial_median_degree, '\n', 'Rewiring Percent = .02')
el_flat$facets <- paste0('Median Degree of First Protesters: ', el_flat$initial_median_degree, '\n', 'Rewiring Percent = .34')

all_sw <- rbind(increasing, el_increasing, el_flat)


# Plot first zone, increasing
ggplot(all_sw, aes(x=repression_rate, y=active_percent)) + geom_line() + xlab('Repression Rate') + ylab('Avg. Percent Mobilization of Network') + theme_bw() + facet_wrap(~facets) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ylim(0, .03)
ggsave('Figures/June2019/Siegel2018_small_world_b_m2_repressionRate.jpeg', plot=last_plot())









########
##
## WORK
##
########
makePlots(data=all, outfile='repression_rate', xlabPretty='Repression Rate')
makePlots(data=all, outfile='networkStructureControl', xlabPretty='Network Control Parameter')
makePlots(data=all, outfile='initial_size', xlabPretty='Initial Protest Size')
makePlots(data=all, outfile='initial_median_degree', xlabPretty='Median Degree of Initial Protesters')
makePlots(data=all, outfile='initial_mean_degree', xlabPretty='Mean Degree of Initial Protesters')
makePlots(data=all, outfile='num_communities', xlabPretty='Number of Starting Communities')
makePlots(data=all, outfile='initial_density', xlabPretty='Initial Protester Density')
makePlots(data=all, outfile='initial_neighborhood_clustering', xlabPretty='Within Protester Clustering')
makePlots(data=all, outfile='initial_nodes_clustering', xlabPretty='Overall Protester Clustering')
makePlots(data=all, outfile='initial_mean_eigen_centrality', xlabPretty='Average Eigenvector Centrality of Initial Protesters')
makePlots(data=all, outfile='starting_mean_commmunity_protest_percent', xlabPretty='Community Penetration of Initial Protesters')


########
##
## WORK, DO NOT DISTINGUISH BY CORE, INITIAL SIZE GT 4
##
########
all_short <- all[all$initial_size >= 5,]

makePlots(data=all_short, outfile='repression_rate', xlabPretty='Repression Rate', outfileMore='initialSizegt4', groupByCore=TRUE)
makePlots(data=all_short, outfile='networkStructureControl', xlabPretty='Network Control Parameter', outfileMore='initialSizegt4', groupByCore=TRUE)
makePlots(data=all_short, outfile='initial_size', xlabPretty='Initial Protest Size', outfileMore='initialSizegt4', groupByCore=TRUE)
makePlots(data=all_short, outfile='initial_median_degree', xlabPretty='Median Degree of Initial Protesters', outfileMore='initialSizegt4', groupByCore=TRUE)
makePlots(data=all_short, outfile='initial_mean_degree', xlabPretty='Mean Degree of Initial Protesters', outfileMore='initialSizegt4', groupByCore=TRUE)
makePlots(data=all_short, outfile='num_communities', xlabPretty='Number of Starting Communities', outfileMore='initialSizegt4', groupByCore=TRUE)
makePlots(data=all_short, outfile='initial_density', xlabPretty='Initial Protester Density', outfileMore='initialSizegt4', groupByCore=TRUE)
makePlots(data=all_short, outfile='initial_neighborhood_clustering', xlabPretty='Within Protester Clustering', outfileMore='initialSizegt4', groupByCore=TRUE)
makePlots(data=all_short, outfile='initial_nodes_clustering', xlabPretty='Protester Clustering', outfileMore='initialSizegt4', groupByCore=TRUE)
makePlots(data=all_short, outfile='starting_mean_commmunity_protest_percent', xlabPretty='Community Penetration of Initial Protesters', outfileMore='initialSizegt4', groupByCore=TRUE)
makePlots(data=all_short, outfile='initial_mean_eigen_centrality', xlabPretty='Average Eigenvector Centrality of Initial Protesters', outfileMore='initialSizegt4', groupByCore=TRUE)


########
##
## WORK, SUBSET FOR INITIAL SIZE BETWEEN 5, 15
##
########
all_short <- all[all$initial_size >= 5 & all$initial_size <= 15,]

makePlots(data=all_short, outfile='repression_rate', xlabPretty='Repression Rate', outfileMore='initialSize5to15')
makePlots(data=all_short, outfile='networkStructureControl', xlabPretty='Network Control Parameter', outfileMore='initialSize5to15')
makePlots(data=all_short, outfile='initial_size', xlabPretty='Initial Protest Size', outfileMore='initialSize5to15')
makePlots(data=all_short, outfile='initial_median_degree', xlabPretty='Median Degree of Initial Protesters', outfileMore='initialSize5to15')
makePlots(data=all_short, outfile='initial_mean_degree', xlabPretty='Mean Degree of Initial Protesters', outfileMore='initialSize5to15')
makePlots(data=all_short, outfile='num_communities', xlabPretty='Number of Starting Communities', outfileMore='initialSize5to15')
makePlots(data=all_short, outfile='initial_density', xlabPretty='Initial Protester Density', outfileMore='initialSize5to15')
makePlots(data=all_short, outfile='initial_neighborhood_clustering', xlabPretty='Within Protester Clustering', outfileMore='initialSize5to15')
makePlots(data=all_short, outfile='initial_nodes_clustering', xlabPretty='Overall Protester Clustering', outfileMore='initialSize5to15')
makePlots(data=all_short, outfile='initial_mean_eigen_centrality', xlabPretty='Average Eigenvector Centrality of Initial Protesters', outfileMore='initialSize5to15')
makePlots(data=all_short, outfile='starting_mean_commmunity_protest_percent', xlabPretty='Community Penetration of Initial Protesters', outfileMore='initialSize5to15')


########
##
## WORK, SUBSET FOR INITIAL SIZE LESS THAN 35
##
########
all_short <- all[all$initial_size <= 35,]

makePlots(data=all_short, outfile='repression_rate', xlabPretty='Repression Rate', outfileMore='initialSizelt35')
makePlots(data=all_short, outfile='networkStructureControl', xlabPretty='Network Control Parameter', outfileMore='initialSizelt35')
makePlots(data=all_short, outfile='initial_size', xlabPretty='Initial Protest Size', outfileMore='initialSizelt35')
makePlots(data=all_short, outfile='initial_median_degree', xlabPretty='Median Degree of Initial Protesters', outfileMore='initialSizelt35')
makePlots(data=all_short, outfile='initial_mean_degree', xlabPretty='Mean Degree of Initial Protesters', outfileMore='initialSizelt35')
makePlots(data=all_short, outfile='num_communities', xlabPretty='Number of Starting Communities', outfileMore='initialSizelt35')
makePlots(data=all_short, outfile='initial_density', xlabPretty='Initial Protester Density', outfileMore='initialSizelt35')
makePlots(data=all_short, outfile='initial_neighborhood_clustering', xlabPretty='Within Protester Clustering', outfileMore='initialSizelt35')
makePlots(data=all_short, outfile='initial_nodes_clustering', xlabPretty='Protester Clustering', outfileMore='initialSizelt35')
makePlots(data=all_short, outfile='initial_mean_eigen_centrality', xlabPretty='Average Eigenvector Centrality of Initial Protesters', outfileMore='initialSizelt35')
makePlots(data=all_short, outfile='starting_mean_commmunity_protest_percent', xlabPretty='Community Penetration of Initial Protesters', outfileMore='initialSizelt35')


########
##
## WORK, SUBSET BASED ON REALISTIC VALUES
##
##  scale-free, see citations from my paper with shane
##	hk, see citations from my paper with shane
##	sw, rounded global clustering and p, looked for p where global clustering = .16.  .3 comes from value used to find hk metric as well
########
sf <- all[all$graph_type=='scale_free_graph',]
sf <- sf[sf$networkStructureControl==2.3,]

hk <- all[all$graph_type=='powerlaw_cluster_graph',]
hk <- hk[hk$networkStructureControl==.315,]

sw <- all[all$graph_type=='watts_strogatz_graph',]
sw <- 	data[[outfile]] <- round(data[[outfile]], 2)

sw <- sw[sw$networkStructureControl==.4,]

all_short <- rbind(sf, hk, sf)

makePlots(data=all_short, outfile='repression_rate', xlabPretty='Repression Rate', outfileMore='realisticNetworkControlValues')
makePlots(data=all_short, outfile='networkStructureControl', xlabPretty='Network Control Parameter', outfileMore='realisticNetworkControlValues')
makePlots(data=all_short, outfile='initial_size', xlabPretty='Initial Protest Size', outfileMore='realisticNetworkControlValues')
makePlots(data=all_short, outfile='initial_median_degree', xlabPretty='Median Degree of Initial Protesters', outfileMore='realisticNetworkControlValues')
makePlots(data=all_short, outfile='initial_mean_degree', xlabPretty='Mean Degree of Initial Protesters', outfileMore='realisticNetworkControlValues')
makePlots(data=all_short, outfile='num_communities', xlabPretty='Number of Starting Communities', outfileMore='realisticNetworkControlValues')
makePlots(data=all_short, outfile='initial_density', xlabPretty='Initial Protester Density', outfileMore='realisticNetworkControlValues')
makePlots(data=all_short, outfile='initial_neighborhood_clustering', xlabPretty='Within Protester Clustering', outfileMore='realisticNetworkControlValues')
makePlots(data=all_short, outfile='initial_nodes_clustering', xlabPretty='Protester Clustering', outfileMore='realisticNetworkControlValues')
makePlots(data=all_short, outfile='initial_mean_eigen_centrality', xlabPretty='Average Eigenvector Centrality of Initial Protesters', outfileMore='realisticNetworkControlValues')
makePlots(data=all_short, outfile='starting_mean_commmunity_protest_percent', xlabPretty='Community Penetration of Initial Protesters', outfileMore='realisticNetworkControlValues')



############################################################
############################################################
#
#		DO NOT DISTINGUISH BY CORE
#
############################################################
############################################################



########
##
## WORK, DO NOT DISTINGUISH BY CORE
##
########
makePlots(data=all, outfile='repression_rate', xlabPretty='Repression Rate', groupByCore=FALSE)
makePlots(data=all, outfile='networkStructureControl', xlabPretty='Network Control Parameter', groupByCore=FALSE)
makePlots(data=all, outfile='initial_size', xlabPretty='Initial Protest Size',groupByCore=FALSE)
makePlots(data=all, outfile='initial_median_degree', xlabPretty='Median Degree of Initial Protesters', groupByCore=FALSE)
makePlots(data=all, outfile='initial_mean_degree', xlabPretty='Mean Degree of Initial Protesters', groupByCore=FALSE)
makePlots(data=all, outfile='num_communities', xlabPretty='Number of Starting Communities', groupByCore=FALSE)
makePlots(data=all, outfile='initial_density', xlabPretty='Initial Protester Density', groupByCore=FALSE)
makePlots(data=all, outfile='initial_neighborhood_clustering', xlabPretty='Within Protester Clustering', groupByCore=FALSE)
makePlots(data=all, outfile='initial_nodes_clustering', xlabPretty='Protester Clustering', groupByCore=FALSE)
makePlots(data=all, outfile='initial_mean_eigen_centrality', xlabPretty='Average Eigenvector Centrality of Initial Protesters', groupByCore=FALSE)
makePlots(data=all, outfile='starting_mean_commmunity_protest_percent', xlabPretty='Community Penetration of Initial Protesters', groupByCore=FALSE)


########
##
## WORK, DO NOT DISTINGUISH BY CORE, INITIAL SIZE GT 4
##
########
all_short <- all[all$initial_size >= 5,]

makePlots(data=all_short, outfile='repression_rate', xlabPretty='Repression Rate', outfileMore='initialSizegt4', groupByCore=FALSE)
makePlots(data=all_short, outfile='networkStructureControl', xlabPretty='Network Control Parameter', outfileMore='initialSizegt4', groupByCore=FALSE)
makePlots(data=all_short, outfile='initial_size', xlabPretty='Initial Protest Size', outfileMore='initialSizegt4', groupByCore=FALSE)
makePlots(data=all_short, outfile='initial_median_degree', xlabPretty='Median Degree of Initial Protesters', outfileMore='initialSizegt4', groupByCore=FALSE)
makePlots(data=all_short, outfile='initial_mean_degree', xlabPretty='Mean Degree of Initial Protesters', outfileMore='initialSizegt4', groupByCore=FALSE)
makePlots(data=all_short, outfile='num_communities', xlabPretty='Number of Starting Communities', outfileMore='initialSizegt4', groupByCore=FALSE)
makePlots(data=all_short, outfile='initial_density', xlabPretty='Initial Protester Density', outfileMore='initialSizegt4', groupByCore=FALSE)
makePlots(data=all_short, outfile='initial_neighborhood_clustering', xlabPretty='Within Protester Clustering', outfileMore='initialSizegt4', groupByCore=FALSE)
makePlots(data=all_short, outfile='initial_nodes_clustering', xlabPretty='Protester Clustering', outfileMore='initialSizegt4', groupByCore=FALSE)
makePlots(data=all_short, outfile='starting_mean_commmunity_protest_percent', xlabPretty='Community Penetration of Initial Protesters', outfileMore='initialSizegt4', groupByCore=FALSE)
makePlots(data=all_short, outfile='initial_mean_eigen_centrality', xlabPretty='Average Eigenvector Centrality of Initial Protesters', outfileMore='initialSizegt4', groupByCore=FALSE)



########
##
## WORK, SUBSET FOR INITIAL SIZE BETWEEN 5, 15
##
########
all_short <- all[all$initial_size >= 5 & all$initial_size <= 15,]

makePlots(data=all_short, outfile='repression_rate', xlabPretty='Repression Rate', outfileMore='initialSize5to15', groupByCore=FALSE)
makePlots(data=all_short, outfile='networkStructureControl', xlabPretty='Network Control Parameter', outfileMore='initialSize5to15', groupByCore=FALSE)
makePlots(data=all_short, outfile='initial_size', xlabPretty='Initial Protest Size', outfileMore='initialSize5to15', groupByCore=FALSE)
makePlots(data=all_short, outfile='initial_median_degree', xlabPretty='Median Degree of Initial Protesters', outfileMore='initialSize5to15', groupByCore=FALSE)
makePlots(data=all_short, outfile='initial_mean_degree', xlabPretty='Mean Degree of Initial Protesters', outfileMore='initialSize5to15', groupByCore=FALSE)
makePlots(data=all_short, outfile='num_communities', xlabPretty='Number of Starting Communities', outfileMore='initialSize5to15', groupByCore=FALSE)
makePlots(data=all_short, outfile='initial_density', xlabPretty='Initial Protester Density', outfileMore='initialSize5to15', groupByCore=FALSE)
makePlots(data=all_short, outfile='initial_neighborhood_clustering', xlabPretty='Within Protester Clustering', outfileMore='initialSize5to15', groupByCore=FALSE)
makePlots(data=all_short, outfile='initial_nodes_clustering', xlabPretty='Overall Protester Clustering', outfileMore='initialSize5to15', groupByCore=FALSE)
makePlots(data=all_short, outfile='initial_mean_eigen_centrality', xlabPretty='Average Eigenvector Centrality of Initial Protesters', outfileMore='initialSize5to15', groupByCore=FALSE)
makePlots(data=all_short, outfile='starting_mean_commmunity_protest_percent', xlabPretty='Community Penetration of Initial Protesters', outfileMore='initialSize5to15', groupByCore=FALSE)


########
##
## WORK, DO NOT DISTINGUISH BY CORE, 4 < INITIAL SIZE < 35
##
########
all_short <- all[all$initial_size >= 5 & all$initial_size <= 34,]

makePlots(data=all_short, outfile='repression_rate', xlabPretty='Repression Rate', outfileMore='initialSizegt4lt35', groupByCore=FALSE)
makePlots(data=all_short, outfile='networkStructureControl', xlabPretty='Network Control Parameter', outfileMore='initialSizegt4lt35', groupByCore=FALSE)
makePlots(data=all_short, outfile='initial_size', xlabPretty='Initial Protest Size', outfileMore='initialSizegt4lt35', groupByCore=FALSE)
makePlots(data=all_short, outfile='initial_median_degree', xlabPretty='Median Degree of Initial Protesters', outfileMore='initialSizegt4lt35', groupByCore=FALSE)
makePlots(data=all_short, outfile='initial_mean_degree', xlabPretty='Mean Degree of Initial Protesters', outfileMore='initialSizegt4lt35', groupByCore=FALSE)
makePlots(data=all_short, outfile='num_communities', xlabPretty='Number of Starting Communities', outfileMore='initialSizegt4lt35', groupByCore=FALSE)
makePlots(data=all_short, outfile='initial_density', xlabPretty='Initial Protester Density', outfileMore='initialSizegt4lt35', groupByCore=FALSE)
makePlots(data=all_short, outfile='initial_neighborhood_clustering', xlabPretty='Within Protester Clustering', outfileMore='initialSizegt4lt35', groupByCore=FALSE)
makePlots(data=all_short, outfile='initial_nodes_clustering', xlabPretty='Protester Clustering', outfileMore='initialSizegt4lt35', groupByCore=FALSE)
makePlots(data=all_short, outfile='initial_mean_eigen_centrality', xlabPretty='Average Eigenvector Centrality of Initial Protesters', outfileMore='initialSizegt4lt35', groupByCore=FALSE)
makePlots(data=all_short, outfile='starting_mean_commmunity_protest_percent', xlabPretty='Community Penetration of Initial Protesters', outfileMore='initialSizegt4lt35', groupByCore=FALSE)
makePlots(data=all_short, outfile='initial_mean_degree', xlabPretty='Mean Degree of Initial Protesters', outfileMore='initialSizegt4lt35', groupByCore=FALSE)




########
##
## WORK, SUBSET BASED ON REALISTIC VALUES BUT DO NOT DISTINGUISH BY CORE
##
##  scale-free, see citations from my paper with shane
##	hk, see citations from my paper with shane
##	sw, rounded global clustering and p, looked for p where global clustering = .16.  .3 comes from value used to find hk metric as well
########
sf <- all[all$graph_type=='scale_free_graph',]
sf <- sf[sf$networkStructureControl==2.3,]

hk <- all[all$graph_type=='powerlaw_cluster_graph',]
hk <- hk[hk$networkStructureControl==.315,]

sw <- all[all$graph_type=='watts_strogatz_graph',]
sw <- sw[sw$networkStructureControl==.4,]

all_short <- rbind(sf, hk, sf)

makePlots(data=all_short, outfile='repression_rate', xlabPretty='Repression Rate', outfileMore='realisticNetworkControlValues_noCoreDistinction', groupByCore=FALSE)
makePlots(data=all_short, outfile='networkStructureControl', xlabPretty='Network Control Parameter', outfileMore='realisticNetworkControlValues_noCoreDistinction', groupByCore=FALSE)
makePlots(data=all_short, outfile='initial_size', xlabPretty='Initial Protest Size', outfileMore='realisticNetworkControlValues_noCoreDistinction', groupByCore=FALSE)
makePlots(data=all_short, outfile='initial_mean_degree', xlabPretty='Mean Degree of Initial Protesters', outfileMore='realisticNetworkControlValues_noCoreDistinction', groupByCore=FALSE)
makePlots(data=all_short, outfile='initial_median_degree', xlabPretty='Median Degree of Initial Protesters', outfileMore='realisticNetworkControlValues_noCoreDistinction', groupByCore=FALSE)
makePlots(data=all_short, outfile='num_communities', xlabPretty='Number of Starting Communities', outfileMore='realisticNetworkControlValues_noCoreDistinction', groupByCore=FALSE)
makePlots(data=all_short, outfile='initial_density', xlabPretty='Initial Protester Density', outfileMore='realisticNetworkControlValues_noCoreDistinction', groupByCore=FALSE)
makePlots(data=all_short, outfile='initial_neighborhood_clustering', xlabPretty='Within Protester Clustering', outfileMore='realisticNetworkControlValues_noCoreDistinction', groupByCore=FALSE)
makePlots(data=all_short, outfile='initial_nodes_clustering', xlabPretty='Protester Clustering', outfileMore='realisticNetworkControlValues_noCoreDistinction', groupByCore=FALSE)
makePlots(data=all_short, outfile='initial_mean_eigen_centrality', xlabPretty='Average Eigenvector Centrality of Initial Protesters', outfileMore='realisticNetworkControlValues_noCoreDistinction', groupByCore=FALSE)
makePlots(data=all_short, outfile='starting_mean_commmunity_protest_percent', xlabPretty='Community Penetration of Initial Protesters', outfileMore='realisticNetworkControlValues_noCoreDistinction', groupByCore=FALSE)































#### as of 05.21.2019 at 8:30 p.m., above work is quite plenty.  Need to spend time prettifying and that will probably get tricky, and certainly be annoying to test, but it a distinct enough task that I can pause here.
###### networkStructureControl, for changing p (rewiring, clustering) or alpha (scaling)
outfile <- 'initial_size'
all[[outfile]] <- round(all[[outfile]], 2)
results <- data.frame(all %>% group_by(graph_type, UQ(as.name(outfile)), coreStart) %>% summarize(successDouble=mean(successDouble), successGreater35=mean(successGreater35), repressionWorked=mean(repressionWorked)))
results_melted <- melt(results, id.vars= c('graph_type', outfile, 'coreStart'))

results_melted$graph_type[results_melted$graph_type=='powerlaw_cluster_graph'] <- 'Holme-Kim'
results_melted$graph_type[results_melted$graph_type=='watts_strogatz_graph'] <- 'Small World'
results_melted$graph_type[results_melted$graph_type=='scale_free_graph'] <- 'Scale Free'

results_melted$variable <- as.character(results_melted$variable)

results_melted$variable[results_melted$variable=='successDouble'] <- 'Protest Doubles'
results_melted$variable[results_melted$variable=='successGreater35'] <- 'Final Protest >= 35'
results_melted$variable[results_melted$variable=='repressionWorked'] <- 'Repression Works'

results_melted$coreStart <- ifelse(results_melted$coreStart==1, "Yes", "No")




jpeg(paste0('Figures/June2019/coreVperiphSuccess', outfile, '.jpg'), res=300, width=8, height=8, units='in')
ggplot(results_melted, aes(x=UQ(as.name(outfile)), y=value, linetype=factor(coreStart))) + geom_line() + facet_grid(variable ~ graph_type) + theme_bw() + xlab("Network Control Parameter") + ylab("Percent of Protests Successful") + scale_y_continuous(trans=log10_trans(), breaks=trans_breaks('log10', function(x) 10^x), labels=trans_format('log10', math_format(10^.x))) + scale_linetype_discrete(name='Core Start')
dev.off()


###### repression_rate, for changing p (rewiring, clustering) or alpha (scaling)
outfile <- 'repression_rate'
results <- data.frame(all %>% group_by(graph_type, UQ(as.name(outfile)), coreStart) %>% summarize(successDouble=mean(successDouble), successGreater35=mean(successGreater35), repressionWorked=mean(repressionWorked)))
results_melted <- melt(results, id.vars= c('graph_type', outfile, 'coreStart'))

# Below is trying to facet with more variables
jpeg(paste0('Figures/June2019/coreVperiphSuccess', outfile, '.jpg'), res=300, width=8, height=8, units='in')
ggplot(results_melted, aes(x=UQ(as.name(outfile)), y=value, linetype=factor(coreStart))) + geom_line() + facet_grid(variable ~ graph_type)
dev.off()




