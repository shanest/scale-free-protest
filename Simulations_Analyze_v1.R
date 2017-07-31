'''
The purpose of this script is to analyze the data created from the simulations.

Designed to work on R 3.2.2
'''

##########################
##
##	GLOBALS
##
##########################
library(poweRlaw)

setwd('path/to/working/directory')

# Copied from TwitterProtestSize/Scripts/TwitterProtestSize_TestingDistributions_v4.R
compareDistributions <- function(data, distribution){
	m1 <- displ$new(data)
	est <- estimate_xmin(m1)
	m1$setXmin(est$xmin)
	m1$setPars(estimate_pars(m1))
	gof <- est$gof

	if(tolower(distribution)=='poisson'){
		m2 <- dispois$new(data)
		m2$setXmin(m1$getXmin())  # Use xmin from power law, so compare over range power law holds.
		m2$setPars(estimate_pars(m2))

		comp <- compare_distributions(m1, m2)
		comp2 <- compare_distributions(m2, m1)

	}

	if(tolower(distribution)=='exponential'){
		m2 <- disexp$new(data)
		m2$setXmin(m1$getXmin())  # Use xmin from power law, so compare over range power law holds.
		m2$setPars(estimate_pars(m2))

		comp <- compare_distributions(m1, m2)
		comp2 <- compare_distributions(m2, m1)
	}

	if(tolower(distribution)=='lognormal'){
		m2 <- dislnorm$new(data)
		m2$setXmin(m1$getXmin())  # Use xmin from power law, so compare over range power law holds.
		m2$setPars(estimate_pars(m2))

		comp <- compare_distributions(m1, m2)
		comp2 <- compare_distributions(m2, m1)
	}

	return(list(data.frame(xmin=est$xmin, alpha=m1$pars, pl_gof=gof, test_stat = comp$test_statistic, two_sided=comp$p_two_sided, one_sided_vsother=comp$p_one_sided, onesided_othervs=comp2$p_one_sided), m1))
	# comp$p_two_sided  # Very low p, suggests difference between power law and poisson.  Null is they are the same, so we reject.
	# comp$p_one_sided  # p still low, confirms power law better
}


# Data is the model result from compareDistributions
makeFigure_pl <- function(data, filename){
	pdf(filename)
	par(mar=c(5,4,.5,2))
	plot(data, xlab='Diffusion Size', ylab='P(X>x)', pch=20, col='grey80')
	lines(data, col='black', lwd=3)
	text(x=max(exp1$final_size)*.1, y=.75, paste('Size Threshold: ', data$xmin), pos=2)
	text(x=max(exp1$final_size)*.1, y=.6, paste('Scaling Parameter:', round(data$pars, 2)), pos=2)
	dev.off()
}
##########################
##
##	DATA
##
##########################
# Experiment 1
exp1 <- read.csv('Data/NetworkSimulation/exp1.csv')

# Standard deviation by initial density
exp1$initial_density <- round(exp1$initial_density, 2)
sd_density <- data.frame(exp1 %>% group_by(initial_density) %>% summarize(diffusion_sd = sd(final_size)))
sd_clustering <- data.frame(exp1 %>% group_by(initial_clustering) %>% summarize(diffusion_sd = sd(final_size)))
sd_initial_size <- data.frame(exp1 %>% group_by(initial_size) %>% summarize(diffusion_sd = sd(final_size)))



##########################
##
##	ANALYZE
##
##########################
result_exp1 <- compareDistributions(exp1$final_size, 'lognormal')  # xmin is 1, alpha 1.214, gof .264, LR is -3.2, oneside is .999, .0006



##########################
##
##	FIGURES
##
##########################
makeFigure_pl(data=result_exp1[[2]], filename='Figures/NetworkSimulation_exp1.pdf')  # Result: Huge average diffusion, like Centola finds

####### Final diffusion size
pdf('Figures/NetworkSimulation_FinalSize_InitialDensity.pdf')
par(mar=c(5,4,.5,2))
plot(exp1$initial_density, exp1$final_size, xlab='Initial Density of Neighborhood', ylab='Final Diffusion Size', pch=20)
dev.off()

pdf('Figures/NetworkSimulation_FinalSize_InitialClustering.pdf')
par(mar=c(5,4,.5,2))
plot(exp1$initial_clustering, exp1$final_size, xlab='Initial Clustering of Neighborhood', ylab='Final Diffusion Size', pch=20)
dev.off()

pdf('Figures/NetworkSimulation_FinalSize_InitialSize.pdf')
par(mar=c(5,4,.5,2))
plot(exp1$initial_size, exp1$final_size, xlab='Initial Size of Neighborhood', ylab='Final Diffusion Size', pch=20)
dev.off()


####### Standard deviation by initial conditions
pdf('Figures/NetworkSimulation_FinalSizeSD_InitialDensity.pdf')
par(mar=c(5,4,.5,2))
plot(sd_density$initial_density, sd_density$diffusion_sd, xlab='Initial Density of Neighborhood', ylab='Standard Deviation, Final Diffusion Size', pch=20)
dev.off()

pdf('Figures/NetworkSimulation_FinalSizeSD_InitialClustering.pdf')
par(mar=c(5,4,.5,2))
plot(sd_clustering$initial_clustering, sd_clustering$diffusion_sd, xlab='Initial Clustering of Neighborhood', ylab='Standard Deviation, Final Diffusion Size', pch=20)
dev.off()

pdf('Figures/NetworkSimulation_FinalSizeSD_InitialNeighborhood.pdf')
par(mar=c(5,4,.5,2))
plot(sd_initial_size$initial_size, sd_initial_size$diffusion_sd, xlab='Initial Size of Neighborhood', ylab='Standard Deviation, Final Diffusion Size', pch=20)
dev.off()



