'''
The purpose of this script is to see if regression inferences change based on regression on different subsets of data.
I do not know at this point (01.22.2018) if this will be for the network simulation paper or the methods paper, but it will be somewhere to show the effect of scaling.

Results after analyzing exp2 as regression:
1. Scaling matters by affecting size of initial neighborhood, not in distribution of connections itself.  So, dictator cares most about initial size, pushes repression so that harder to get initial sizes.

v2. Started June 5, 2018, with new experimental results.
	1. Delete plots of SEs, coefficients, and histograms
v3. 1. Added pipeline for functions, deleted old code that did manually
	2. Delete regressions where add more data because can just do by adding new loops with percentages
v4. 1. Updated exp7 code to get the scaling parameter, will do others as needed.  But other calls to pipeline() will not work until updated because lacking two arguments.
v5. 1. Make the main plots histograms instead of lines.  x-axis is t-statistic - DONE
v6. 1. Made v5, realized that did not save regression results by percentage.  So need to rerun entire pipeline.
v8. 1. Make font size larger in plots - DONE
	2. Remove extra underscore from file names - DONE
v9. 1. Log transformation of protest size for regressions
		- See that more consistent results for scale-free than small world
		- See that more consistent results as sample size increases
v10. 1. Add exp13 where p=1. Results: Most are 90% or above not significant, but then cool that sig. can be + or -.
	 2. Examine exp13 where global clustering between .1 and .2. Results: same randomness as random.
'''
##########################
##
##	GLOBALS
##
##########################
set.seed(480193)

library(dplyr)  # For sample_n
library(MASS)  # For glm.nb

setwd('/Users/Zack/Documents/UCLA/Research/ProtestTheory')


# Runs regression, returns t-statistic of each coefficient on final size
runRegressions <- function(data, model, percentage, family, scaling=FALSE, repression=FALSE){
	keep <- round(nrow(data)*percentage)
	sub <- sample_n(data, keep, replace=TRUE)

	if(family != 'negbinom'){
		reg <- summary(lm(model, sub))$coefficients
		temp <- reg[,3]  # t-statistics
		temp2 <- reg[,2] # standard error
		temp3 <- reg[,1] # coefficients
	}

	# Commented out for now because fails to converge, not worth worrying about at this early stage.
	# if(family == 'negbinom'){
	# 	temp <- summary(glm.nb(model, data=sub, start=rep(0, )))$coefficients[,3]  # t-statistics

	# }
	if(scaling==FALSE & repression==FALSE){
		df <- data.frame(initial_size=temp[['initial_size']], initial_density=temp[['initial_density']], initial_global_clustering=temp[['initial_global_clustering']], initial_mean_degree=temp[['initial_mean_degree']], initial_median_degree=temp[['initial_median_degree']], initial_neighborhood_clustering=temp[['initial_neighborhood_clustering']], initial_nodes_clustering=temp[['initial_nodes_clustering']], initial_size=temp[['initial_size']], initial_size_se=temp2[['initial_size']], initial_density_se=temp2[['initial_density']], initial_global_clustering_se=temp2[['initial_global_clustering']], initial_mean_degree_se=temp2[['initial_mean_degree']], initial_median_degree_se=temp2[['initial_median_degree']], initial_neighborhood_clustering_se=temp2[['initial_neighborhood_clustering']], initial_nodes_clustering_se=temp2[['initial_nodes_clustering']], initial_size_se=temp2[['initial_size']],initial_size_c=temp3[['initial_size']], initial_density_c=temp3[['initial_density']], initial_global_clustering_c=temp3[['initial_global_clustering']], initial_mean_degree=temp3[['initial_mean_degree']], initial_median_degree_c=temp3[['initial_median_degree']], initial_neighborhood_clustering_c=temp3[['initial_neighborhood_clustering']], initial_nodes_clustering_c=temp3[['initial_nodes_clustering']], initial_size_c=temp3[['initial_size']])
	}

	if(scaling==TRUE & repression==FALSE){
	df <- data.frame(scaling_parameter=temp[['scaling_parameter']], initial_size=temp[['initial_size']], initial_density=temp[['initial_density']], initial_global_clustering=temp[['initial_global_clustering']], initial_mean_degree=temp[['initial_mean_degree']], initial_median_degree=temp[['initial_median_degree']], initial_neighborhood_clustering=temp[['initial_neighborhood_clustering']], initial_nodes_clustering=temp[['initial_nodes_clustering']], initial_size=temp[['initial_size']], scaling_parameter_se=temp2[['scaling_parameter']], initial_size_se=temp2[['initial_size']], initial_density_se=temp2[['initial_density']], initial_global_clustering_se=temp2[['initial_global_clustering']], initial_mean_degree_se=temp2[['initial_mean_degree']], initial_median_degree_se=temp2[['initial_median_degree']], initial_neighborhood_clustering_se=temp2[['initial_neighborhood_clustering']], initial_nodes_clustering_se=temp2[['initial_nodes_clustering']], initial_size_se=temp2[['initial_size']], scaling_parameter_c=temp3[['scaling_parameter']], initial_size_c=temp3[['initial_size']], initial_density_c=temp3[['initial_density']], initial_global_clustering_c=temp3[['initial_global_clustering']], initial_mean_degree=temp3[['initial_mean_degree']], initial_median_degree_c=temp3[['initial_median_degree']], initial_neighborhood_clustering_c=temp3[['initial_neighborhood_clustering']], initial_nodes_clustering_c=temp3[['initial_nodes_clustering']], initial_size_c=temp3[['initial_size']])
	}

	if(scaling==TRUE & repression==TRUE){
	df <- data.frame(repression_rate=temp[['repression_rate']], scaling_parameter=temp[['scaling_parameter']], initial_size=temp[['initial_size']], initial_density=temp[['initial_density']], initial_global_clustering=temp[['initial_global_clustering']], initial_mean_degree=temp[['initial_mean_degree']], initial_median_degree=temp[['initial_median_degree']], initial_neighborhood_clustering=temp[['initial_neighborhood_clustering']], initial_nodes_clustering=temp[['initial_nodes_clustering']], initial_size=temp[['initial_size']], repression_rate_se=temp2[['repression_rate']], scaling_parameter_se=temp2[['scaling_parameter']], initial_size_se=temp2[['initial_size']], initial_density_se=temp2[['initial_density']], initial_global_clustering_se=temp2[['initial_global_clustering']], initial_mean_degree_se=temp2[['initial_mean_degree']], initial_median_degree_se=temp2[['initial_median_degree']], initial_neighborhood_clustering_se=temp2[['initial_neighborhood_clustering']], initial_nodes_clustering_se=temp2[['initial_nodes_clustering']], initial_size_se=temp2[['initial_size']], repression_rate_c=temp3[['repression_rate']], scaling_parameter_c=temp3[['scaling_parameter']], initial_size_c=temp3[['initial_size']], initial_density_c=temp3[['initial_density']], initial_global_clustering_c=temp3[['initial_global_clustering']], initial_mean_degree=temp3[['initial_mean_degree']], initial_median_degree_c=temp3[['initial_median_degree']], initial_neighborhood_clustering_c=temp3[['initial_neighborhood_clustering']], initial_nodes_clustering_c=temp3[['initial_nodes_clustering']], initial_size_c=temp3[['initial_size']])
	}
	return(df)
}

makePlot <- function(data, filecustomizer, xlab='Regression Number', ylab='T-statistic', ylimits, significance, lines=TRUE, legend=TRUE){
	ymin <- round(min(data))-1
	ymax <- round(max(data))*1.1
	pdf(filecustomizer)
	# First two if statements control use of y range
	if(abs(max(data))-abs(min(data)) > 1){
		plot(data, type='l', xlab=xlab, ylab=ylab, ylim=c(ymin,ymax), cex.main=1.5, cex.lab=1.5, cex.axis=1.5)
	}
	if(abs(max(data))-abs(min(data)) <= 1){
		plot(data, type='l', xlab=xlab, ylab=ylab, cex.main=1.5, cex.lab=1.5, cex.axis=1.5)
	}
	if(lines==TRUE){
		abline(a=significance, b=0, lty=3)
		abline(a=-1*significance, b=0, lty=3)
	}
	abline(a=0, b=0, lty=3)
	if(legend==TRUE){
		legend('topleft', paste0(round(sum(abs(data) < significance)/length(data)*100,2), '% are Insignificant'), cex=1.5)
	}
	dev.off()
}

makeHist <- function(data, filecustomizer, xlabel='', significance, legend=TRUE, ymax=NULL){

	xmin <- min(data)
	xmax <- max(data)


	# if the maximum cofficient is below significance to the left
	if(xmax <= significance*-1){
		xmin <- xmin*1.1
		xmax <- xmax*.9
	}

	# if the minimum coefficient is above significance to the right
	if(xmin >= significance){
		xmin <- xmin*.9
		xmax <- xmax*1.1
	}

	# if min and max span significance to the left
	if(xmin <= significance*-1 & xmax >= significance*1){
		xmin <- xmin*1.1
		xmax <- 4
	}

	# if min and max span significance to the right
	if(xmax >= significance & xmin <= significance){
		xmin <- -4
		xmax <- xmax*1.1
	}



	temp <- density(data)
	height <- max(temp$y)  # Vestigial, use if decide want to dynamically change text, arrow height

	pdf(filecustomizer)
	plot(temp, xlab=xlabel, main='', xlim=c(xmin, xmax), ylim=c(0,.5), cex.main=1.5, cex.lab=1.5, cex.axis=1.5)
	lines(x=rep(significance, 50), y=seq(0,.4, length.out=50), lty=3)
	lines(x=rep(significance*-1, 50), y=seq(0,.4, length.out=50), lty=3)

	#abline(v=c(significance, significance*-1), lty=3)
	text(x=significance, y=.38, 'Significant', pos=4, cex=1.5)
	arrows(x0=significance, y0=.35, x1=significance+.52, y1=.35, length=.1)
	text(x=significance*-1, y=.38, 'Significant', pos=2, cex=1.5)
	arrows(x0=significance*-1, y0=.35, x1=significance*-1-.52, y1=.35, length=.1)
	if(legend==TRUE){
		legend('topleft', paste0(round(sum(abs(data) < significance)/length(data)*100,2), '% are Insignificant'), cex=1.5)
	}
	dev.off()
}



pipeline <- function(data, model, expName, percUse, theseData = NULL, tscore, scaling, repression){
	if(is.null(theseData)){
		results <- data.frame(initial_size=NULL, initial_density=NULL, initial_global_clustering=NULL, initial_mean_degree=NULL, initial_median_degree=NULL, initial_neighborhood_clustering=NULL, initial_nodes_clustering=NULL, initial_size=NULL, initial_size_se=NULL, initial_density_se=NULL, initial_global_clustering_se=NULL, initial_mean_degree_se=NULL, initial_median_degree_se=NULL, initial_neighborhood_clustering_se=NULL, initial_nodes_clustering_se=NULL, initial_size_se=NULL,initial_size_c=NULL, initial_density_c=NULL, initial_global_clustering_c=NULL, initial_mean_degree=NULL, initial_median_degree_c=NULL, initial_neighborhood_clustering_c=NULL, initial_nodes_clustering_c=NULL, initial_size_c=NULL)



		#results_nb <- results

		for(i in 1:1000){
			print(paste('On loop ', i))
			r <- runRegressions(data=data, model=model, percentage=percUse, family='regular', scaling=scaling, repression=repression)
			results <- rbind(results, r)

			# Often fails to converge, so comment out for now
			# r <- runRegressions(data=exp1, model=m1, percentage=.5, family='negbinom')
			# results_nb <- rbind(results_nb, r)
		}

		write.csv(results, paste0('Data/NetworkSimulation/Summer2018/logged_', expName, '_regressions', percUse*100, 'percentage.csv'))

	}

	if(is.null(theseData)==FALSE){
		results <- theseData
	}

	sample <- round(nrow(data)*percUse)

	# t-statistics
	makePlot(results$initial_size, filecustomizer=paste0('Figures/Summer2018/', expName, 'Regressions_initial_size_tstat_n', sample, '_', percUse*100, 'percentage.pdf'), significance=tscore)
	makePlot(results$initial_density, filecustomizer=paste0('Figures/Summer2018/', expName, 'Regressions_initial_density_tstat_n', sample,  '_', percUse*100, 'percentage.pdf'), significance=tscore)
	makePlot(results$initial_global_clustering, filecustomizer=paste0('Figures/Summer2018/', expName, 'Regressions_initial_global_clustering_tstat_n', sample,  '_', percUse*100, 'percentage.pdf'), significance=tscore)
	makePlot(results$initial_mean_degree, filecustomizer=paste0('Figures/Summer2018/', expName, 'Regressions_initial_mean_degree_tstat_n', sample,  '_', percUse*100, 'percentage.pdf'), significance=tscore)
	makePlot(results$initial_median_degree, filecustomizer=paste0('Figures/Summer2018/', expName, 'Regressions_initial_median_degree_tstat_n', sample,  '_', percUse*100, 'percentage.pdf'), significance=tscore)
	makePlot(results$initial_neighborhood_clustering, filecustomizer=paste0('Figures/Summer2018/', expName, 'Regressions_initial_neighborhood_clustering_tstat_n', sample,  '_', percUse*100, 'percentage.pdf'), significance=tscore)
	makePlot(results$initial_nodes_clustering, filecustomizer=paste0('Figures/Summer2018/', expName, 'Regressions_initial_nodes_clustering_tstat_n', sample,  '_', percUse*100, 'percentage.pdf'), significance=tscore)
	if(scaling==TRUE){
		makePlot(results$scaling_parameter, filecustomizer=paste0('Figures/Summer2018/', expName, 'Regressions_scaling_parameter_tstat_n', sample,  '_', percUse*100, 'percentage.pdf'), significance=tscore)
	}
	if(repression==TRUE){
		makePlot(results$repression_rate, filecustomizer=paste0('Figures/Summer2018/', expName, 'Regressions_repression_rate_tstat_n', sample,  '_', percUse*100, 'percentage.pdf'), significance=tscore)
	}
	# Histogram
	# t-statistics
	makeHist(results$initial_size, filecustomizer=paste0('Figures/Summer2018/', expName, 'Histogram_Regressions_initial_size_tstat_n', sample, '_', percUse*100, 'percentage.pdf'), significance=tscore)
	makeHist(results$initial_density, filecustomizer=paste0('Figures/Summer2018/', expName, 'Histogram_Regressions_initial_density_tstat_n', sample,  '_', percUse*100, 'percentage.pdf'), significance=tscore)
	makeHist(results$initial_global_clustering, filecustomizer=paste0('Figures/Summer2018/', expName, 'Histogram_Regressions_initial_global_clustering_tstat_n', sample,  '_', percUse*100, 'percentage.pdf'), significance=tscore)
	makeHist(results$initial_mean_degree, filecustomizer=paste0('Figures/Summer2018/', expName, 'Histogram_Regressions_initial_mean_degree_tstat_n', sample,  '_', percUse*100, 'percentage.pdf'), significance=tscore)
	makeHist(results$initial_median_degree, filecustomizer=paste0('Figures/Summer2018/', expName, 'Histogram_Regressions_initial_median_degree_tstat_n', sample,  '_', percUse*100, 'percentage.pdf'), significance=tscore)
	makeHist(results$initial_neighborhood_clustering, filecustomizer=paste0('Figures/Summer2018/', expName, 'Histogram_Regressions_initial_neighborhood_clustering_tstat_n', sample,  '_', percUse*100, 'percentage.pdf'), significance=tscore)
	makeHist(results$initial_nodes_clustering, filecustomizer=paste0('Figures/Summer2018/', expName, 'Histogram_Regressions_initial_nodes_clustering_tstat_n', sample,  '_', percUse*100, 'percentage.pdf'), significance=tscore)
	if(scaling==TRUE){
		makeHist(results$scaling_parameter, filecustomizer=paste0('Figures/Summer2018/', expName, 'Histogram_Regressions_scaling_parameter_tstat_n', sample,  '_', percUse*100, 'percentage.pdf'), significance=tscore)
	}
	if(repression==TRUE){
		makeHist(results$repression_rate, filecustomizer=paste0('Figures/Summer2018/', expName, 'Histogram_Regressions_repression_rate_tstat_n', sample,  '_', percUse*100, 'percentage.pdf'), significance=tscore)
	}


}

##########################
##
##
##	DATA
##		NB: From v5 on, do not need anything below because the regressions have been run.  Only uncomment if have new model results.
##		NB2: Uncommented for v6 because the pipeilne code was not saving new files for each percentage of regression kept.
##
##########################

# # Experiment 3
# exp3_path <- 'Data/NetworkSimulation/Summer2018/exp3-threshold_type=normal-num_nodes=1000-repression_type=edge_removal-graph_type=scale_free_graph.csv'
# exp3 <- read.csv(exp3_path, stringsAsFactors=FALSE)
# exp3$normalized <- round(1000*exp3$final_size/exp3$total_nodes)

# # Experiment 1
# exp1 <- exp3[exp3$repression_rate==0 & exp3$scaling==2.3,]

# # Experiment 2
# exp2 <- exp3[exp3$repression_rate==0,]

# # # Experiment 4: No longer doing.  That is where varied network size.

# # # Experiment 5: PowerLaw Cluster, node removal, normal threshold.
# # #  Not run because later decided to run with varying p.

# Experiment 6
exp6_path <- 'Data/NetworkSimulation/Summer2018/exp6-threshold_type=normal-num_nodes=1000-repression_type=node_removal-graph_type=watts_strogatz_graph-k=50.csv'
exp6 <- read.csv(exp6_path, stringsAsFactors=FALSE)
exp6$normalized <- round(1000*exp6$final_size/exp6$total_nodes)

# # Experiment 7
exp7_path <- 'Data/NetworkSimulation/Summer2018/exp7-graph_type=scale_free_graph-threshold_type=uniform-repression_type=node_removal-num_nodes=1000.csv'
exp7 <- read.csv(exp7_path, stringsAsFactors=FALSE)
exp7$normalized <- round(1000*exp7$final_size/exp7$total_nodes)

# Experiment 8
exp8_path <- 'Data/NetworkSimulation/Summer2018/exp8-graph_type=scale_free_graph-threshold_type=normal-repression_type=node_removal-num_nodes=1000.csv'
exp8 <- read.csv(exp8_path, stringsAsFactors=FALSE)
exp8$normalized <- round(1000*exp8$final_size/exp8$total_nodes)

# # Experiment 9
# exp9_path <- 'Data/NetworkSimulation/Summer2018/exp9-threshold_type=uniform-m=3-num_nodes=3000-repression_type=node_removal-graph_type=powerlaw_cluster_graph.csv'
# exp9 <- read.csv(exp9_path, stringsAsFactors=FALSE)
# exp9$normalized <- round(1000*exp9$final_size/exp9$total_nodes)

# Experiment 10
# exp10_path <- 'Data/NetworkSimulation/Summer2018/exp10-threshold_type=uniform-repression_type=edge_removal-m=3-graph_type=powerlaw_cluster_graph-num_nodes=3000.csv'
# exp10 <- read.csv(exp10_path, stringsAsFactors=FALSE)
# exp10$normalized <- round(1000*exp10$final_size/exp10$total_nodes)

# Experiment 11
# exp11_path <- 'Data/NetworkSimulation/Summer2018/exp11-threshold_type=normal-m=3-num_nodes=3000-repression_type=node_removal-graph_type=powerlaw_cluster_graph.csv'
# exp11 <- read.csv(exp11_path, stringsAsFactors=FALSE)
# exp11$normalized <- round(1000*exp11$final_size/exp11$total_nodes)

# # Experiment 12
# exp12_path <- 'Data/NetworkSimulation/Summer2018/exp12-m=3-graph_type=powerlaw_cluster_graph-threshold_type=normal-num_nodes=3000-repression_type=edge_removal.csv'
# exp12 <- read.csv(exp12_path, stringsAsFactors=FALSE)
# exp12$normalized <- round(1000*exp12$final_size/exp12$total_nodes)

# Experiment 13
exp13_path <- 'Data/NetworkSimulation/Summer2018/exp13-threshold_type=uniform-num_nodes=1000-repression_type=node_removal-graph_type=watts_strogatz_graph-k=50.csv'
exp13 <- read.csv(exp13_path, stringsAsFactors=FALSE)
exp13$normalized <- round(1000*exp13$final_size/exp13$total_nodes)

exp13_lattice <- exp13[exp13$p==0,]
exp13_random <- exp13[exp13$p==1,]

exp13_clustered <- exp13[exp13$initial_global_clustering >=.1 & exp13$initial_global_clustering <= .2,]



# Experiment 14
exp14_path <- 'Data/NetworkSimulation/Summer2018/exp14-threshold_type=uniform-num_nodes=1000-repression_type=edge_removal-graph_type=watts_strogatz_graph-k=50.csv'
exp14 <- read.csv(exp14_path, stringsAsFactors=FALSE)
exp14$normalized <- round(1000*exp14$final_size/exp14$total_nodes)

# Experiment 15
exp15_path <- 'Data/NetworkSimulation/Summer2018/exp15-threshold_type=normal-num_nodes=1000-repression_type=edge_removal-graph_type=watts_strogatz_graph-k=50.csv'
exp15 <- read.csv(exp15_path, stringsAsFactors=FALSE)
exp15$normalized <- round(1000*exp15$final_size/exp15$total_nodes)

# Experiment 16
exp16_path <- 'Data/NetworkSimulation/Summer2018/exp16-threshold_type=uniform-num_nodes=1000-repression_type=edge_removal-graph_type=scale_free_graph.csv'
exp16 <- read.csv(exp16_path, stringsAsFactors=FALSE)
exp16$normalized <- round(1000*exp16$final_size/exp16$total_nodes)

##########################
##
##
##	WORK
##
##  #NB: If logging, make sure to log initial size but keep variable name the same
##########################
m <- log(normalized+1, 10) ~ initial_density + initial_global_clustering + initial_mean_degree + initial_median_degree + initial_neighborhood_clustering + initial_nodes_clustering + initial_size+ scaling_parameter
mb <- log(normalized+1, 10) ~ initial_density + initial_global_clustering + initial_mean_degree + initial_median_degree + initial_neighborhood_clustering + initial_nodes_clustering + initial_size  # No scaling_parameter, same as
m2 <- log(normalized+1, 10) ~ initial_density + initial_global_clustering + initial_mean_degree + initial_median_degree + initial_neighborhood_clustering + initial_nodes_clustering + initial_size + scaling_parameter + repression_rate



###### BELOW IS FROM V4, WHEN MADE INITIAL REGRESSION DATA.
### EXPERIMENT 1
# pipeline(data=exp1, model=m, expName='logged_exp1', percUse=.5, theseData=NULL, tscore=1.96, scaling=FALSE, repression=FALSE)
# pipeline(data=exp1, model=m, expName='logged_exp1', percUse=.25, theseData=NULL, tscore=1.96, scaling=FALSE, repression=FALSE)
# pipeline(data=exp1, model=m, expName='logged_exp1', percUse=.1, theseData=NULL, tscore=1.96, scaling=FALSE, repression=FALSE)

# ### EXPERIMENT 2
# pipeline(data=exp2, model=m, expName='logged_exp2', percUse=.5, theseData=NULL, tscore=1.96, scaling=FALSE, repression=FALSE)
# pipeline(data=exp2, model=m, expName='logged_exp2', percUse=.25, theseData=NULL, tscore=1.96, scaling=FALSE, repression=FALSE)
# pipeline(data=exp2, model=m, expName='logged_exp2', percUse=.1, theseData=NULL, tscore=1.96, scaling=FALSE, repression=FALSE)
# pipeline(data=exp2, model=m, expName='logged_exp2', percUse=.01, theseData=NULL, tscore=1.96, scaling=FALSE, repression=FALSE)

# ## EXPERIMENT 3
# ## Less percentage because more data
# pipeline(data=exp3, model=m2, expName='logged_exp3', percUse=.1, theseData=NULL, tscore=1.96, scaling=TRUE, repression=TRUE)
# pipeline(data=exp3, model=m2, expName='logged_exp3', percUse=.01, theseData=NULL, tscore=1.96, scaling=TRUE, repression=TRUE)


### EXPERIMENT 6
exp6$initial_size <- log(exp6$initial_size+1, 10)
exp6$scaling_parameter <- exp6$p  # No scaling parameter, but the variation I care about is p
pipeline(data=exp6, model=m, expName='logged_exp6', percUse=.5, theseData=NULL, tscore=1.96, scaling=TRUE, repression=FALSE)
pipeline(data=exp6, model=m, expName='logged_exp6', percUse=.25, theseData=NULL, tscore=1.96, scaling=TRUE, repression=FALSE)
# pipeline(data=exp6, model=m, expName='logged_exp6', percUse=.1, theseData=NULL, tscore=1.96, scaling=TRUE, repression=FALSE)
pipeline(data=exp6, model=m, expName='logged_exp6', percUse=.01, theseData=NULL, tscore=1.96, scaling=TRUE, repression=FALSE)

# ### EXPERIMENT 7
exp7$initial_size <- log(exp7$initial_size+1, 10)
pipeline(data=exp7, model=m, expName='logged_exp7', percUse=.5, theseData=NULL, tscore=1.96, scaling=TRUE, repression=FALSE)
pipeline(data=exp7, model=m, expName='logged_exp7', percUse=.25, theseData=NULL, tscore=1.96, scaling=TRUE, repression=FALSE)
pipeline(data=exp7, model=m, expName='logged_exp7', percUse=.1, theseData=NULL, tscore=1.96, scaling=TRUE, repression=FALSE)
pipeline(data=exp7, model=m, expName='logged_exp7', percUse=.01, theseData=NULL, tscore=1.96, scaling=TRUE, repression=FALSE)

# ### EXPERIMENT 7a
exp7a <- exp7[exp7$scaling_parameter==2.3,]
exp7a$initial_size <- log(exp7a$initial_size+1, 10)
pipeline(data=exp7a, model=mb, expName='logged_exp7a', percUse=.5, theseData=NULL, tscore=1.96, scaling=FALSE, repression=FALSE)
pipeline(data=exp7a, model=mb, expName='logged_exp7a', percUse=.25, theseData=NULL, tscore=1.96, scaling=FALSE, repression=FALSE)
pipeline(data=exp7a, model=mb, expName='logged_exp7a', percUse=.1, theseData=NULL, tscore=1.96, scaling=FALSE, repression=FALSE)


# ### EXPERIMENT 7, DROP LARGE INITIAL SIZE
initialSize <- 50
exp7_init <- exp7[exp7$initial_size < initialSize,]
exp7_init$initial_size <- log(exp7_init$initial_size+1, 10)
pipeline(data=exp7, model=m, expName='logged_exp7_SmallInitialSize', percUse=.5, theseData=NULL, tscore=1.96, scaling=TRUE, repression=FALSE)
pipeline(data=exp7, model=m, expName='logged_exp7_SmallInitialSize', percUse=.25, theseData=NULL, tscore=1.96, scaling=TRUE, repression=FALSE)
pipeline(data=exp7, model=m, expName='logged_exp7_SmallInitialSize', percUse=.1, theseData=NULL, tscore=1.96, scaling=TRUE, repression=FALSE)
pipeline(data=exp7, model=m, expName='logged_exp7_SmallInitialSize', percUse=.01, theseData=NULL, tscore=1.96, scaling=TRUE, repression=FALSE)

# ### EXPERIMENT 7a, DROP LARGE INITIAL SIZE
initialSize <- 50
exp7a <- exp7[exp7$scaling_parameter==2.3,]
exp7a_init <- exp7a[exp7a$initial_size < initialSize,]
exp7a_init$initial_size <- log(exp7a_init$initial_size+1, 10)
pipeline(data=exp7a, model=m, expName='logged_exp7a_SmallInitialSize', percUse=.5, theseData=NULL, tscore=1.96, scaling=FALSE, repression=FALSE)
pipeline(data=exp7a, model=m, expName='logged_exp7a_SmallInitialSize', percUse=.25, theseData=NULL, tscore=1.96, scaling=FALSE, repression=FALSE)
pipeline(data=exp7a, model=m, expName='logged_exp7a_SmallInitialSize', percUse=.1, theseData=NULL, tscore=1.96, scaling=FALSE, repression=FALSE)




### EXPERIMENT 8
# pipeline(data=exp8, model=m, expName='logged_exp8', percUse=.5, theseData=NULL, tscore=1.96, scaling=FALSE, repression=FALSE)
# pipeline(data=exp8, model=m, expName='logged_exp8', percUse=.25, theseData=NULL, tscore=1.96, scaling=FALSE, repression=FALSE)
# pipeline(data=exp8, model=m, expName='logged_exp8', percUse=.1, theseData=NULL, tscore=1.96, scaling=FALSE, repression=FALSE)
# pipeline(data=exp8, model=m, expName='logged_exp8', percUse=.01, theseData=NULL, tscore=1.96, scaling=FALSE, repression=FALSE)

### EXPERIMENT 9
# exp9$scaling_parameter <- exp9$p # Scale-free cluster, thing to vary is p, probability of triangle
# pipeline(data=exp9, model=m, expName='logged_exp9', percUse=.5, theseData=NULL, tscore=1.96, scaling=FALSE, repression=FALSE)
# pipeline(data=exp9, model=m, expName='logged_exp9', percUse=.25, theseData=NULL, tscore=1.96, scaling=FALSE, repression=FALSE)
# pipeline(data=exp9, model=m, expName='logged_exp9', percUse=.1, theseData=NULL, tscore=1.96, scaling=FALSE, repression=FALSE)
# pipeline(data=exp9, model=m, expName='logged_exp9', percUse=.01, theseData=NULL, tscore=1.96, scaling=FALSE, repression=FALSE)

# ## EXPERIMENT 10
# exp10$scaling_parameter <- exp10$p # Scale-free cluster, thing to vary is p, probability of triangle
# pipeline(data=exp10, model=m2, expName='logged_exp10', percUse=.1, theseData=NULL, tscore=1.96, scaling=TRUE, repression=TRUE)
# pipeline(data=exp10, model=m2, expName='logged_exp10', percUse=.01, theseData=NULL, tscore=1.96, scaling=TRUE, repression=TRUE)

### EXPERIMENT 11
# exp11$scaling_parameter <- exp11$p # Scale-free cluster, thing to vary is p, probability of triangle
# pipeline(data=exp11, model=m, expName='logged_exp11', percUse=.5, theseData=NULL, tscore=1.96, scaling=TRUE, repression=FALSE)
# pipeline(data=exp11, model=m, expName='logged_exp11', percUse=.25, theseData=NULL, tscore=1.96, scaling=TRUE, repression=FALSE)
# pipeline(data=exp11, model=m, expName='logged_exp11', percUse=.1, theseData=NULL, tscore=1.96, scaling=TRUE, repression=FALSE)
# pipeline(data=exp11, model=m, expName='logged_exp11', percUse=.01, theseData=NULL, tscore=1.96, scaling=TRUE, repression=FALSE)

# # ### EXPERIMENT 12
# exp12$scaling_parameter <- exp12$p # Scale-free cluster, thing to vary is p, probability of triangle
# pipeline(data=exp12, model=m2, expName='logged_exp12', percUse=.1, theseData=NULL, tscore=1.96, scaling=TRUE, repression=TRUE)
# pipeline(data=exp12, model=m2, expName='logged_exp12', percUse=.01, theseData=NULL, tscore=1.96, scaling=TRUE, repression=TRUE)

# ### EXPERIMENT 13
exp13$scaling_parameter <- exp13$p # Watts-Strogatz, so variation is p, probability of rewiring
exp13$initial_size <- log(exp13$initial_size+1, 10)
pipeline(data=exp13, model=m, expName='logged_exp13', percUse=.5, theseData=NULL, tscore=1.96, scaling=FALSE, repression=FALSE)
pipeline(data=exp13, model=m, expName='logged_exp13', percUse=.25, theseData=NULL, tscore=1.96, scaling=FALSE, repression=FALSE)
pipeline(data=exp13, model=m, expName='logged_exp13', percUse=.1, theseData=NULL, tscore=1.96, scaling=FALSE, repression=FALSE)
pipeline(data=exp13, model=m, expName='logged_exp13', percUse=.01, theseData=NULL, tscore=1.96, scaling=FALSE, repression=FALSE)



# ### EXPERIMENT 13, random
exp13_random$scaling_parameter <- exp13_random$p # Watts-Strogatz, so variation is p, probability of rewiring
exp13_random$initial_size <- log(exp13_random$initial_size+1, 10)
pipeline(data=exp13_random, model=m, expName='logged_exp13_random', percUse=.5, theseData=NULL, tscore=1.96, scaling=FALSE, repression=FALSE)
pipeline(data=exp13_random, model=m, expName='logged_exp13_random', percUse=.25, theseData=NULL, tscore=1.96, scaling=FALSE, repression=FALSE)
pipeline(data=exp13_random, model=m, expName='logged_exp13_random', percUse=.1, theseData=NULL, tscore=1.96, scaling=FALSE, repression=FALSE)
pipeline(data=exp13_random, model=m, expName='logged_exp13_random', percUse=.01, theseData=NULL, tscore=1.96, scaling=FALSE, repression=FALSE)


# ### EXPERIMENT 13, clustered
exp13_clustered$scaling_parameter <- exp13_clustered$p # Watts-Strogatz, so variation is p, probability of rewiring
exp13_clustered$initial_size <- log(exp13_clustered$initial_size+1, 10)
pipeline(data=exp13_clustered, model=m, expName='logged_exp13_clustered', percUse=.5, theseData=NULL, tscore=1.96, scaling=FALSE, repression=FALSE)
pipeline(data=exp13_clustered, model=m, expName='logged_exp13_clustered', percUse=.25, theseData=NULL, tscore=1.96, scaling=FALSE, repression=FALSE)
pipeline(data=exp13_clustered, model=m, expName='logged_exp13_clustered', percUse=.1, theseData=NULL, tscore=1.96, scaling=FALSE, repression=FALSE)
pipeline(data=exp13_clustered, model=m, expName='logged_exp13_clustered', percUse=.01, theseData=NULL, tscore=1.96, scaling=FALSE, repression=FALSE)



### EXPERIMENT 14
exp14$scaling_parameter <- exp14$p # Watts-Strogatz, so variation is p, probability of rewiring
pipeline(data=exp14, model=m2, expName='logged_exp14', percUse=.1, theseData=NULL, tscore=1.96, scaling=FALSE, repression=TRUE)
pipeline(data=exp14, model=m2, expName='logged_exp14', percUse=.01, theseData=NULL, tscore=1.96, scaling=FALSE, repression=TRUE)

### EXPERIMENT 15
exp15$scaling_parameter <- exp15$p # Watts-Strogatz, so variation is p, probability of rewiring
pipeline(data=exp15, model=m2, expName='logged_exp15', percUse=.1, theseData=NULL, tscore=1.96, scaling=FALSE, repression=TRUE)
pipeline(data=exp15, model=m2, expName='logged_exp15', percUse=.01, theseData=NULL, tscore=1.96, scaling=FALSE, repression=TRUE)

### EXPERIMENT 16, is scale-free graph
pipeline(data=exp16, model=m2, expName='logged_exp16', percUse=.1, theseData=NULL, tscore=1.96, scaling=TRUE, repression=TRUE)
pipeline(data=exp16, model=m2, expName='logged_exp16', percUse=.01, theseData=NULL, tscore=1.96, scaling=TRUE, repression=TRUE)

