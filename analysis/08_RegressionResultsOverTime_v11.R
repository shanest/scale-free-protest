'''
The purpose of this script is to see if regression results change as we get new data about an event.

The idea is that it replicates how our knowledge changes as new data comes in with history.

It is a fork of /Users/Zack/Documents/UCLA/Research/ProtestTheory/Scripts/Regressions_VarianceEffects_v5.R.  That one used multiple sets of data, but this script will focus on Mass Mobilization and Mass Mobilization in Autocracies.

v1. 1. Update runRegressions_grow to also be by year
	2. File output, August2020 -> August2020
	NB: As of 08.19.2019, results for Mass Mobilization suggest that get consistent as time increases, so I stop exploring there.  The takeaway is time = data.

v2. Update as result of Network Science R&R.  Want to make histograms for the network models.
	1. June2019 -> August2020
v3. Saving after scratch work on v2
v4. 1. Use results with new repression type
	2. Make plots prettier
v5. 1. Sample trials using repression weights from NAVCO
v6. 1. New just to make sure have old version if mess up samling
v7. Continuing to make new figures for results for revision, now on different subsets.
	1. Make a plot function finally, makeNewPlot
v8. 1. I made plots for realistic values of network and repression wrong.  I read in regression results but had not made the raw data, need to make that raw data first.
v9. Make results for when do edge repression
v10. Update pipeline function to work with threshold variables
v11. Remove some identity scratchwork
	 Update makeHist_gg to take threshold variables
'''


##########################
##
##	GLOBALS
##
##########################
set.seed(01282017)

library(dplyr)  # For sample_n, bind_rows
library(ggplot2)
library(reshape)
library(zoo)  # to replace NAs, na.locf

setwd('/Users/Zack/Documents/UCLA/Research/ProtestTheory')


runRegressions_grow <- function(data, model, startRow=FALSE){
	tstats <- NULL
	coefs <- NULL
	ses <- NULL

	# Parse by year if startRow not given
	if(startRow==FALSE){

		data$year <- as.numeric(data$year)
		data <- subset(data, is.na(year) == FALSE)

		years <- unique(data$year)
		for(i in 2:length(years)){
			print(paste('On year', i, 'of', length(years)))

			temp <- data[data$year <= years[i],]
			m_temp <- summary(lm(model, data=temp))

			t <- as.data.frame(t(as.matrix(m_temp$coefficients[,3])))
			t$n <- i

			c <- as.data.frame(t(as.matrix(m_temp$coefficients[,1])))
			c$n <- i

			s <- as.data.frame(t(as.matrix(m_temp$coefficients[,2])))
			s$n <- i

			tstats <- bind_rows(tstats, t)  # Handles rows with different # of columns.  From dplyr
			coefs <-bind_rows(coefs, c)
			ses <- bind_rows(ses, s)

		}
	}

	# If a row starting number is given
	if(startRow != FALSE){
		for(i in startRow:nrow(data)){
			print(paste('On row', i, 'of', nrow(data)))

			temp <- data[1:i,]
			m_temp <- summary(lm(model, data=temp))

			t <- as.data.frame(t(as.matrix(m_temp$coefficients[,3])))
			t$n <- i

			c <- as.data.frame(t(as.matrix(m_temp$coefficients[,1])))
			c$n <- i

			s <- as.data.frame(t(as.matrix(m_temp$coefficients[,2])))
			s$n <- i

			tstats <- bind_rows(tstats, t)  # Handles rows with different # of columns.  From dplyr
			coefs <-bind_rows(coefs, c)
			ses <- bind_rows(ses, s)
			}
		}
	return (list(tstats=tstats, coefs=coefs, ses=ses))
}

# data is the tstatistics, data2 is the coefficients, data3 is standard errors
makePlots <- function(data, xlabel='Number of Observations', filemodifier, tscore){
	these <- names(data)[grep('factor', names(data), invert=TRUE)]  # For each variable that is not a factor variable
	for(i in 1:length(these)){
		if(is.null(data)==FALSE){
			ymax <- round(max(data[,these[i]], na.rm=TRUE)*1.05,2)
			if(ymax<tscore){
				ymax <- round(tscore*1.05, 2)
			}
			ymin <- round(min(data[,these[i]], na.rm=TRUE)*1.05, 2)
			if(ymin > -1*tscore){
				ymin <- round(-1.05*tscore, 2)
			}
			jpeg(paste0('Figures/August2020/Regression_', filemodifier, '_InferenceChange_', these[i], '.jpeg'), width=5, height=5, units='in', res=300)
			plot(x=data$n, y=data[,these[i]], xlab=xlabel, ylab='T-statisic', type='l', ylim=c(ymin,ymax))
			abline(h=tscore, lty=3)
			abline(h=0, lty=3)
			abline(h=-1*tscore, lty=3)
			legend('topleft', paste0(round(sum(abs(data[,these[i]]) < tscore)/nrow(data)*100,2), '% are Insignificant'))
			dev.off()
		}
	}
}


makeHist <- function(data, filecustomizer, xlabel='', significance, y_max=.5, legend=TRUE, ymax=NULL){
	these <- names(data)[grep('factor', names(data), invert=TRUE)]  # For each variable that is not a factor variable
	for(i in 1:length(these)){

		temp <- density(data[,these[i]], na.rm=TRUE)
		height <- max(temp$y)  # Vestigial, use if decide want to dynamically change text, arrow height


		xmin <- min(temp$x)
		xmax <- max(temp$x)


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



		pdf(paste0('Figures/August2020/Regression_Histogram_', filecustomizer, '_InferenceChange_', these[i], '.pdf'))
		plot(temp, xlab=xlabel, main='', xlim=c(xmin, xmax), ylim=c(0,y_max), cex.main=1.5, cex.lab=1.5, cex.axis=1.5)
		lines(x=rep(significance, 50), y=seq(0,.4, length.out=50), lty=3)
		lines(x=rep(significance*-1, 50), y=seq(0,.4, length.out=50), lty=3)

		#abline(v=c(significance, significance*-1), lty=3)
		text(x=significance, y=.38, 'Significant', pos=4, cex=1.5)
		arrows(x0=significance, y0=.35, x1=significance+.52, y1=.35, length=.1)
		text(x=significance*-1, y=.38, 'Significant', pos=2, cex=1.5)
		arrows(x0=significance*-1, y0=.35, x1=significance*-1-.52, y1=.35, length=.1)
		if(legend==TRUE){
			legend('topleft', paste0(round(sum(abs(temp$x) < significance)/length(temp$x)*100,2), '% are Insignificant'), cex=1.5)
		}
		dev.off()
	}
}

# Take sample of data, make regressions
runRegressions_sample <- function(data, model, percentage, n, weighting, weight_col){
	keep <- round(nrow(data)*percentage)
	if(weighting==FALSE){
		sub <- sample_n(data, keep, replace=TRUE)
	}
	if(weighting==TRUE){
		data$weight_col <- data[[weight_col]]
		sub <- sample_n(data, keep, replace=TRUE, weight=weight_col)
	}
	reg <- summary(lm(model, sub))$coefficients

	these <- rownames(reg)[grep('factor', rownames(reg), invert=TRUE)]
	tresults <- reg[rownames(reg) %in% these,][,3]
	names(tresults)[names(tresults) == '(Intercept)'] <- 'Intercept'
	tresults <- t(data.frame(tresults))
	tresults <- data.frame(tresults)
	tresults$n <- n

	return(tresults)  # Return t-statistics
}


# Repeat the regression sampling a bunch of times
pipeline <- function(data, model, percentage, trials, weighting=FALSE, weight_col=NULL, threshold_vars=FALSE){
	t_stats <- data.frame(Intercept=NULL, protesterviolence=NULL, response_accomodate=NULL, response_arrest=NULL, response_beat=NULL, response_kill=NULL, resopnse_shoot=NULL, n=NULL)
	for(i in 1:trials){
		print(paste0('On ', i))
		temp <- runRegressions_sample(data=data, model=model, percentage=percentage, n=i, weighting=weighting, weight_col=weight_col)

		t_stats <- rbind(t_stats, temp)
	}

	return(t_stats)

}

########
# BELOW THREE ARE MADE AS PART OF REVISION PROCESS
########
# work on each network type (makeResults)
# add sign stuff (processResults)
# manually plot
makeResults <- function(data, sampleSizes, trials=1000, regressionModel=m, weighting, weight_col, networkType){

	results <- NULL
	for(size in sampleSizes){
		temp <- pipeline(data=data, model=regressionModel, percentage=size/nrow(data), trials=1000, weighting=weighting, weight_col=weight_col)
		temp$sampleSize <- size

		results <- rbind(results, temp)
	}
	# sf_sampling <- pipeline(data=sf, model=m, percentage = size/nrow(sf), trials=1000)
	# sworld_sampling <- pipeline(data=sworld, model=m, percentage = size/nrow(sworld), trials=1000)
	# hk_sampling <- pipeline(data=hk, model=m, percentage=size/nrow(hk), trials=1000)

	# sf_sampling$graph_type <- 'Barabasi-Albert'
	# sworld_sampling$graph_type <- 'Watts-Strogatz'
	# hk_sampling$graph_type <- 'Holme-Kim'

	# # combine
	# toPlot <- rbind(sf_sampling, sworld_sampling, hk_sampling)

	fileout <- paste0('Data/NetworkSimulation/August2020/processedData/', networkType, '_filteredSimulationData.csv')
	write.csv(results, fileout, row.names=FALSE)

	return(results)
}


# data = the results from a big loop, all the regressions
# tscore = level of significance wanted, make sure is positive
# networkType = string describing network being analyzed
processResults <- function(data, tscore, networkType){
	tscore <- tscore

	data$repressionSigSign <- cut(data$repression_rate, breaks=c(-Inf, -1*tscore, tscore, Inf), labels=c(-1,0,1), right=FALSE, include.lowest=TRUE)

	fileout <- paste0('Data/NetworkSimulation/August2020/processedData/', networkType, '_regressionData.csv')
	write.csv(data, fileout, row.names=FALSE)

	temp <- data.frame(data %>% group_by(sampleSize) %>% summarize(Neg = sum(repressionSigSign==-1)/n(), Zero = sum(repressionSigSign==0)/n(), Pos = sum(repressionSigSign==1)/n()))

	temp <- melt(temp, id.vars=c('sampleSize'))
	temp$Significance <- temp$variable

	fileout_agg <- paste0('Data/NetworkSimulation/August2020/processedData/', networkType, '_regressionData_toplot.csv')
	write.csv(temp, fileout_agg, row.names=FALSE)

	return(temp)
}

## Make overlapping histograms
makeHist_gg <- function(data, filecustomizer, significance=1.96, y_max=200, legend=TRUE, ymax=NULL, thresholds){
	these <- c('initial_density', 'initial_global_clustering', 'initial_mean_degree', 'initial_median_degree', 'initial_neighborhood_clustering', 'initial_nodes_clustering', 'initial_size', 'scaling_parameter', 'repression_rate')

	if(thresholds==TRUE){
		these <- c(these, 'initial_mean_threshold', 'initial_neighbors_mean_threshold')
	}

	for(i in 1:length(these)){

		fileout <- paste0('Figures/August2020/Regression_Histogram_', filecustomizer, '_InferenceChange_', these[i], '.jpeg')
		data$this <- data[[these[i]]]
		ggplot(data, aes(x=this, fill=graph_type)) +
		geom_histogram(alpha=.5, bins=100, position='identity') +
		geom_vline(xintercept=c(-1*significance, significance), linetype= 'dotted', size=1) +
		 theme_bw() +  xlab('Test Statistic') +
		 ylab('Number of Regressions')  +
		 labs(color='Graph Type') +
		 scale_fill_manual(limits = c('Watts-Strogatz', 'Holme-Kim', 'Barabasi-Albert'), values=c('red', 'green', 'blue')) + theme(axis.text = element_text(size=16), axis.title=element_text(size=16), legend.text=element_text(size=16), legend.title=element_blank(), legend.position=c(.5,.85), panel.grid.minor = element_blank(),panel.grid.major=element_blank(), panel.border=element_blank(), axis.line=element_line(color='black')) + labs(fill='Graph Type') + scale_y_continuous(expand=c(0,.05), limit=c(0,ymax))
		ggsave(fileout, plot=last_plot())
	}
}


### Designed to show results as sample size changes.
makeNewPlot <- function(data, fileCustomizer){
	fileout <- paste0('Figures/August2020/compareNetworks_significance_color_repressionAndSize_', fileCustomizer, '.jpeg')
	ggplot(data, aes(x=sampleSize, y=value, linetype=Significance, color=graph_type)) + geom_line() + theme_bw() +  xlab('Sample Size') + ylab('Percent of Regressions') + theme(axis.text=element_text(size=16), axis.title=element_text(size=16), panel.grid.minor = element_blank(),panel.grid.major=element_blank(), panel.border=element_blank(), axis.line=element_line(color='black'), legend.text=element_text(size=16), legend.position='bottom', legend.box='horizontal', legend.direction='vertical') + labs(color='Graph Type') + scale_color_discrete(limits=c('Watts-Strogatz', 'Holme-Kim', 'Barabasi-Albert')) + scale_y_continuous(expand=c(0,.05))
	ggsave(fileout, plot=last_plot())
}

##########################
##
##	DATA
##
##########################
mm <- read.csv('Data/MassMobilization_Binghampton/mm_public_120615.csv', stringsAsFactors=FALSE)
autocracies <- read.csv('Data/MassMobilizationsAutocracyDataset/EventLevel/events2.csv', stringsAsFactors=FALSE)

# Below is for weights, based on NAVCO
fatalRates <- read.csv('Data/NAVCO/navco3_fatalRates_weight.csv', stringsAsFactors=FALSE)
casualtyRates <- read.csv('Data/NAVCO/navco3_casualtyRates_weight.csv', stringsAsFactors=FALSE)

##########################
##
##	BINGHAMPTON MASS MOBLIZATION
##
##########################

# Add leading 0
mm$startday <- sprintf('%02d', mm$startday)
mm$startmonth <- sprintf('%02d', mm$startmonth)
mm$endday <- sprintf('%02d', mm$endday)
mm$endmonth <- sprintf('%02d', mm$endmonth)

# Make date
#mm$startdate <- paste0(mm$startyear, mm$startmonth, mm$startday)
mm$startdate <- as.Date(paste(mm$startyear, mm$startmonth, mm$startday, sep='-'))
mm$enddate <- as.Date(paste(mm$endyear, mm$endmonth, mm$endday, sep='-'))

# Protest duration
mm$duration <- as.numeric(mm$enddate - mm$startdate)
mm$duration <- mm$duration + 1  # A protest starting and ending on same day is 1 days, not 0.

mm$duration <- ifelse(mm$duration < 1, 1, mm$duration)  # I noticed there were 3 observation where duration < 0.  Inspection showed typo on coders' part and duration should == 0, so update duration to be 0.


# Remove missing participant numbers
mm$participants <- gsub('s|+', '', mm$participants)
mm$participants <- as.numeric(mm$participants)  # NA will occur for strings taht can't be numbers, which is fine

# Order by start date of protest
mm <- mm[order(mm$startdate),]
mm <- mm[is.na(mm$participants) == FALSE,]
mm <- mm[2:nrow(mm),]  # The first row is a bad one, drop

m1 <- log(participants, 10) ~  protesterviolence + response_accomodate + response_arrest + response_beat + response_kill + response_shoot + factor(year) + factor(country)

results <- runRegressions_grow(data=mm, model=m1, startRow=FALSE)
mm_t <- results[['tstats']]
mm_c <- results[['coefs']]
mm_s <- results[['ses']]

makeHist(data=mm_t, filecustomizer='MassMobilizationBinghampton', significance=1.96, y_max=1)

makePlots(data=mm_t, filemodifier='MassMobilizationBinghampton', tscore=1.96)




##########################
##
##	BINGHAMPTON MASS MOBLIZATION, RESAMPLE
##
##########################
m1 <- log(participants, 10) ~  protesterviolence + response_accomodate + response_arrest + response_beat + response_kill + response_shoot + factor(year) + factor(country)

mm_sampling <- pipeline(data=mm, model=m1, percentage=.5, trials=1000)

makeHist(data=mm_sampling, filecustomizer='MassMobilizationBinghampton_50percent_sampling', significance=1.96)



###### ADD VARIABLES FOR TOPIC AND POPULATION.
## Clean, not all values are given as 0, 1.  These three are the ones that need cleaning.
# Protester identity, inspect.  Probably too many values to make factor.
identity <- data.frame(table(mm$protesteridentity))
table(identity$Freq)
identity <- identity[order(identity$Freq, decreasing=TRUE),]
identity[identity$Freq >= 3,]


# NB: The identity variables are highly colinear with the demand variables.  For example, if the demand is labor, the identity is going to involve a string including "worker", "laborer", etc.  Will therefore just make a model with demand variables.
m2 <- log(participants, 10) ~  protesterviolence + response_accomodate + response_arrest + response_beat + response_kill + response_shoot + factor(year) + factor(country) + demand_labor + demand_land + demand_policebrutality + demand_political + demand_price + demand_removal + demand_social

mm_sampling <- pipeline(data=mm, model=m2, percentage=.5, trials=1000, weighting=FALSE, weight_col=NULL)

makeHist(data=mm_sampling, filecustomizer='MassMobilizationBinghampton_50percent_sampling_demands', significance=1.96)

#makePlots(data=mm_sampling, filemodifier='MassMobilization_Binghampton_50percent_sampling', tscore=1.96, xlabel='Regression')


###### MAKE PROTEST DURATION THE OUTCOME
m3 <- log(duration, 10) ~  protesterviolence + response_accomodate + response_arrest + response_beat + response_kill + response_shoot + factor(year) + factor(country)

mm_sampling <- pipeline(data=mm, model=m3, percentage=.5, trials=1000, weighting=FALSE, weight_col=NULL)

makeHist(data=mm_sampling, filecustomizer='MassMobilizationBinghampton_50percent_sampling_DurationDV', significance=1.96)


###### MAKE PROTEST DURATION THE OUTCOME WITH DEMANDS
m4 <- log(duration, 10) ~  protesterviolence + response_accomodate + response_arrest + response_beat + response_kill + response_shoot + factor(year) + factor(country) + demand_labor + demand_land + demand_policebrutality + demand_political + demand_price + demand_removal + demand_social

mm_sampling <- pipeline(data=mm, model=m4, percentage=.5, trials=1000, weighting=FALSE, weight_col=NULL)

makeHist(data=mm_sampling, filecustomizer='MassMobilizationBinghampton_50percent_sampling_demandsDurationDV', significance=1.96)

##########################
##
##	MASS MOBILIZATION AUTOCRACIES, GROW
##
##########################
autocracies$year <- as.numeric(substr(autocracies$event_date, 1,4))


# Remove missing participant numbers
autocracies <- autocracies[is.na(autocracies$avgnumparticipants) == FALSE,]

# Order by start date of protest
autocracies <- autocracies[order(autocracies$event_date),]

m1 <- log(avgnumparticipants, 10) ~  maxscope + maxpartviolence + maxsecengagement + factor(year) + factor(cowcode)
results <- runRegressions_grow(data=autocracies, model=m1, startRow=1000)
autocracies_t <- results[['tstats']]


makeHist(data=autocracies_t, filecustomizer='MassMobilizationsAutocracies', significance=1.96, y_max=1)

#makePlots(data=autocracies_t, filemodifier='MassMobilizationsAutocracies', tscore=1.96, xlabel='Sample Size')


##########################
##
##	MASS MOBILIZATION AUTOCRACIES, RESAMPLE
##
##########################
m1 <- log(avgnumparticipants, 10) ~  maxscope + maxpartviolence + maxsecengagement + as.numeric(year) + factor(cowcode)

autocracies_sampling <- pipeline(data=autocracies, model=m1, percentage=.5, trials=1000)
makeHist(data=autocracies_sampling, filecustomizer='MassMobilizationAutocracies_50percent_sampling', significance=1.96)

#makePlots(data=autocracies_sampling, filemodifier='MassMobilizationAutocracies_50percent_sampling', tscore=1.96, xlabel='Regression')


##########################
##
##	NETWORK EXPERIMENTS, RESAMPLE.  ORIGINAL REPRESSION
##   I copy the code to load data from Scripts/PolNet_05_Simulations_RepressionFigures_v12.R on 08.10.2020
##
##########################
# Below is added on 07.29.2019.  Is from /Users/Zack/Documents/UCLA/Research/ProtestTheory/Scripts/PolNet_01_ProcessSimulationData_v2.R
all <- read.csv('Data/NetworkSimulation/Summer2019/processedData/01_experimentsCombined_repressionNarrow_eigenCore.csv', stringsAsFactors=FALSE)
all$normalized <- round(1000*all$active_nodes/all$total_nodes)

# Make names to match model m
all$scaling_parameter <- all$networkStructureControl

# Add fatality, casualty weights based on NAVCO
rates <- data.frame(RepressionRate=unique(all$repression_rate))
# interpolate missing values
fatalRates <- merge(fatalRates, rates, by.x='RepressionRate', by.y='RepressionRate', all.y=TRUE)
fatalRates$fatality_weight <- na.locf(fatalRates$Percentage)  # replaces NAs with previous observation

casualtyRates <- merge(casualtyRates, rates, by.y='RepressionRate', all.y=TRUE)
casualtyRates$casualty_weight <- na.locf(casualtyRates$Percentage)


all <- merge(all, fatalRates[,c('RepressionRate', 'fatality_weight')], by.x='repression_rate', by.y='RepressionRate', all.y=FALSE)
all <- merge(all, casualtyRates[,c('RepressionRate', 'casualty_weight')], by.x='repression_rate', by.y='RepressionRate', all.y=FALSE)


# Separate out the networks
sf <- all[all$graph_type=='scale_free_graph',]
sworld <- all[all$graph_type=='watts_strogatz_graph',]
hk <- all[all$graph_type=='powerlaw_cluster_graph',]

m <- log(normalized+1, 10) ~ initial_density + initial_global_clustering + initial_mean_degree + initial_median_degree + initial_neighborhood_clustering + initial_nodes_clustering + initial_size + scaling_parameter + repression_rate




# Using full range of repression
n <- c(100, 250, 500, 1000, 1500, 2000)
for(size in n){
	sf_sampling <- pipeline(data=sf, model=m, percentage = size/nrow(sf), trials=1000)
	sworld_sampling <- pipeline(data=sworld, model=m, percentage = size/nrow(sworld), trials=1000)
	hk_sampling <- pipeline(data=hk, model=m, percentage=size/nrow(hk), trials=1000)

	sf_sampling$graph_type <- 'Barabasi-Albert'
	sworld_sampling$graph_type <- 'Watts-Strogatz'
	hk_sampling$graph_type <- 'Holme-Kim'

	# combine
	toPlot <- rbind(sf_sampling, sworld_sampling, hk_sampling)
	makeHist_gg(data=toPlot, filecustomizer=paste0('experiments_sampleSize', size), ymax=150)
}


# Using realistic range of repression, 5.75% fatality rate from NAVCO
for(size in n){
	sf_short <- subset(sf, repression_rate <= .0575)
	sworld_short <- subset(sworld, repression_rate <= .0575)
	hk_short <- subset(hk, repression_rate <= .0575)

	sf_sampling <- pipeline(data=sf_short, model=m, percentage = size/nrow(sf_short), trials=1000)
	sworld_sampling <- pipeline(data=sworld_short, model=m, percentage = size/nrow(sworld_short), trials=1000)
	hk_sampling <- pipeline(data=hk_short, model=m, percentage=size/nrow(hk_short), trials=1000)

	sf_sampling$graph_type <- 'Barabasi-Albert'
	sworld_sampling$graph_type <- 'Watts-Strogatz'
	hk_sampling$graph_type <- 'Holme-Kim'

	# combine
	toPlot <- rbind(sf_sampling, sworld_sampling, hk_sampling)
	makeHist_gg(data=toPlot, filecustomizer=paste0('experiments_sampleSize', size, '_repression575Perc'), ymax=150)
}


# Using realistic range of repression, 10.05% casualty rate from NAVCO
for(size in n){
	sf_short <- subset(sf, repression_rate <= .105)
	sworld_short <- subset(sworld, repression_rate <= .105)
	hk_short <- subset(hk, repression_rate <= .105)

	sf_sampling <- pipeline(data=sf_short, model=m, percentage = size/nrow(sf_short), trials=1000)
	sworld_sampling <- pipeline(data=sworld_short, model=m, percentage = size/nrow(sworld_short), trials=1000)
	hk_sampling <- pipeline(data=hk_short, model=m, percentage=size/nrow(hk_short), trials=1000)

	sf_sampling$graph_type <- 'Barabasi-Albert'
	sworld_sampling$graph_type <- 'Watts-Strogatz'
	hk_sampling$graph_type <- 'Holme-Kim'

	# combine
	toPlot <- rbind(sf_sampling, sworld_sampling, hk_sampling)
	makeHist_gg(data=toPlot, filecustomizer=paste0('experiments_sampleSize', size, '_casualty105Perc'), ymax=150)
}


# Sampling using weights based on NAVCO fatality data
n <- c(100, 250, 500, 1000, 1500, 2000)
for(size in n){
	sf_sampling <- pipeline(data=sf, model=m, percentage = size/nrow(sf), trials=1000, weighting=TRUE, weight_col='fatality_weight')
	sworld_sampling <- pipeline(data=sworld, model=m, percentage = size/nrow(sworld), trials=1000, weighting=TRUE, weight_col='fatality_weight')
	hk_sampling <- pipeline(data=hk, model=m, percentage=size/nrow(hk), trials=1000, weighting=TRUE, weight_col='fatality_weight')

	sf_sampling$graph_type <- 'Barabasi-Albert'
	sworld_sampling$graph_type <- 'Watts-Strogatz'
	hk_sampling$graph_type <- 'Holme-Kim'

	# combine
	toPlot <- rbind(sf_sampling, sworld_sampling, hk_sampling)
	makeHist_gg(data=toPlot, filecustomizer=paste0('experiments_sampleSize', size, '_fatalityWeighted'), ymax=150)
}


# Sampling using weights based on NAVCO casualty data
# As of 08.13.2020, have not run.
n <- c(100, 250, 500, 1000, 1500, 2000)
for(size in n){
	sf_sampling <- pipeline(data=sf, model=m, percentage = size/nrow(sf), trials=1000, weighting=TRUE, weight_col='casualty_weight')
	sworld_sampling <- pipeline(data=sworld, model=m, percentage = size/nrow(sworld), trials=1000, weighting=TRUE, weight_col='casualty_weight')
	hk_sampling <- pipeline(data=hk, model=m, percentage=size/nrow(hk), trials=1000, weighting=TRUE, weight_col='casualty_weight')

	sf_sampling$graph_type <- 'Barabasi-Albert'
	sworld_sampling$graph_type <- 'Watts-Strogatz'
	hk_sampling$graph_type <- 'Holme-Kim'

	# combine
	toPlot <- rbind(sf_sampling, sworld_sampling, hk_sampling)
	makeHist_gg(data=toPlot, filecustomizer=paste0('experiments_sampleSize', size, '_casualtyWeighted'), ymax=150)
}


####################################################
####################################################
####################################################
######
##	NETWORK EXPERIMENTS, RESAMPLE, SHOW RESULTS AS SAMPLE SIZE CHANGES
##   I copy the structure from Scripts/PolNet_05_Simulations_RepressionFigures_v12.R on 08.10.2020
######
####################################################
####################################################
####################################################
### Sequence of sample sizes
n <- seq(from=50, to=2000, by=50)

#### For all repression rates
results_sf <- makeResults(data=sf, sampleSizes=n, trials=1000, regressionModel=m)
results_sworld <- makeResults(data=sworld, sampleSizes=n, trials=1000, regressionModel=m)
results_hk <- makeResults(data=hk, sampleSizes=n, trials=1000, regressionModel=m)

toPlot_sf <- processResults(data=results_sf, tscore=1.96, networkType='3wayCompare_hist_scaleFree')
toPlot_sworld <- processResults(data=results_sworld, tscore=1.96, networkType='3wayCompare_hist_smallWorld')
toPlot_hk <- processResults(data=results_hk, tscore=1.96, networkType='3wayCompare_hist_holmeKim')

toPlot_sf$graph_type <- 'Barabasi-Albert'
toPlot_sworld$graph_type <- 'Watts-Strogatz'
toPlot_hk$graph_type <- 'Holme-Kim'

# Read below to save time
toPlot_sf <- read.csv('Data/NetworkSimulation/August2020/processedData/3wayCompare_hist_scaleFree_regressionData_toplot.csv')
toPlot_sworld <- read.csv('Data/NetworkSimulation/August2020/processedData/3wayCompare_hist_smallWorld_regressionData_toplot.csv')
toPlot_hk <- read.csv('Data/NetworkSimulation/August2020/processedData/3wayCompare_hist_holmeKim_regressionData_toplot.csv')

toPlot_sf$graph_type <- 'Barabasi-Albert'
toPlot_sworld$graph_type <- 'Watts-Strogatz'
toPlot_hk$graph_type <- 'Holme-Kim'


# Combine
toPlot <- rbind(toPlot_sf, toPlot_sworld, toPlot_hk)


# Plot
ggplot(toPlot, aes(x=sampleSize, y=value, linetype=Significance, color=graph_type)) + geom_line() + theme_bw() +  xlab('Sample Size') + ylab('Percent of Regressions') + labs(color='Graph Type') + scale_color_discrete(limits=c('Watts-Strogatz', 'Holme-Kim', 'Barabasi-Albert')) + theme(axis.text=element_text(size=16), axis.title=element_text(size=16), panel.grid.minor = element_blank(),panel.grid.major=element_blank(), panel.border=element_blank(), axis.line=element_line(color='black'), legend.text=element_text(size=16), legend.title=element_text(size=16), legend.position='bottom', legend.box='horizontal', legend.direction='vertical') + scale_y_continuous(expand=c(0,.05))
ggsave('Figures/August2020/compareNetworks_significance_color_repressionAndSize.jpg', plot=last_plot())


#### For realistic casualty rates
sf_short <- subset(sf, repression_rate <= .105)
sworld_short <- subset(sworld, repression_rate <= .105)
hk_short <- subset(hk, repression_rate <= .105)

results_sf <- makeResults(data=sf_short, sampleSizes=n, trials=1000, regressionModel=m)
results_sworld <- makeResults(data=sworld_short, sampleSizes=n, trials=1000, regressionModel=m)
results_hk <- makeResults(data=hk_short, sampleSizes=n, trials=1000, regressionModel=m)

toPlot_sf <- processResults(data=results_sf, tscore=1.96, networkType='3wayCompare_hist_scaleFree_repression105')
toPlot_sworld <- processResults(data=results_sworld, tscore=1.96, networkType='3wayCompare_hist_smallWorld_repression105')
toPlot_hk <- processResults(data=results_hk, tscore=1.96, networkType='3wayCompare_hist_holmeKim_repression105')

# Read below to save time
toPlot_sf <- read.csv('Data/NetworkSimulation/August2020/processedData/3wayCompare_hist_scaleFree_repression105_regressionData_toplot.csv')
toPlot_sworld <- read.csv('Data/NetworkSimulation/August2020/processedData/3wayCompare_hist_smallWorld_repression105_regressionData_toplot.csv')
toPlot_hk <- read.csv('Data/NetworkSimulation/August2020/processedData/3wayCompare_hist_holmeKim_repression105_regressionData_toplot.csv')

toPlot_sf$graph_type <- 'Barabasi-Albert'
toPlot_sworld$graph_type <- 'Watts-Strogatz'
toPlot_hk$graph_type <- 'Holme-Kim'

# Combine
toPlot <- rbind(toPlot_sf, toPlot_sworld, toPlot_hk)

# Plot
ggplot(toPlot, aes(x=sampleSize, y=value, linetype=Significance, color=graph_type)) + geom_line() + theme_bw() +  xlab('Sample Size') + ylab('Percent of Regressions') + labs(color='Graph Type') + scale_color_discrete(limits=c('Watts-Strogatz', 'Holme-Kim', 'Barabasi-Albert')) + theme(axis.text=element_text(size=16), axis.title=element_text(size=16), panel.grid.minor = element_blank(),panel.grid.major=element_blank(), panel.border=element_blank(), axis.line=element_line(color='black'), legend.text=element_text(size=16), legend.title=element_text(size=16), legend.position='bottom', legend.box='horizontal', legend.direction='vertical') + scale_y_continuous(expand=c(0,.05))
ggsave('Figures/August2020/compareNetworks_significance_color_repressionAndSize_repression105.jpg', plot=last_plot())


### For casualty rates, weighted by NAVCO
results_sf <- makeResults(data=sf, sampleSizes=n, trials=1000, regressionModel=m, weighting=TRUE, weight_col='casualty_weight')
results_sworld <- makeResults(data=sworld, sampleSizes=n, trials=1000, regressionModel=m, weighting=TRUE, weight_col='casualty_weight')
results_hk <- makeResults(data=hk, sampleSizes=n, trials=1000, regressionModel=m, weighting=TRUE, weight_col='casualty_weight')

toPlot_sf <- processResults(data=results_sf, tscore=1.96, networkType='3wayCompare_hist_scaleFree_repressionCasualtyWeight')
toPlot_sworld <- processResults(data=results_sworld, tscore=1.96, networkType='3wayCompare_hist_smallWorld_repressionCasualtyWeight')
toPlot_hk <- processResults(data=results_hk, tscore=1.96, networkType='3wayCompare_hist_holmeKim_repressionCasualtyWeight')

# Read below to save time
toPlot_sf <- read.csv('Data/NetworkSimulation/August2020/processedData/3wayCompare_hist_scaleFree_repressionCasualtyWeight_regressionData_toplot.csv')
toPlot_sworld <- read.csv('Data/NetworkSimulation/August2020/processedData/3wayCompare_hist_smallWorld_repressionCasualtyWeight_regressionData_toplot.csv')
toPlot_hk <- read.csv('Data/NetworkSimulation/August2020/processedData/3wayCompare_hist_holmeKim_repressionCasualtyWeight_regressionData_toplot.csv')

toPlot_sf$graph_type <- 'Barabasi-Albert'
toPlot_sworld$graph_type <- 'Watts-Strogatz'
toPlot_hk$graph_type <- 'Holme-Kim'

# Combine
toPlot <- rbind(toPlot_sf, toPlot_sworld, toPlot_hk)

# Plot
ggplot(toPlot, aes(x=sampleSize, y=value, linetype=Significance, color=graph_type)) + geom_line() + theme_bw() +  xlab('Sample Size') + ylab('Percent of Regressions') + labs(color='Graph Type') + scale_color_discrete(limits=c('Watts-Strogatz', 'Holme-Kim', 'Barabasi-Albert')) + theme(axis.text=element_text(size=16), axis.title=element_text(size=16), panel.grid.minor = element_blank(),panel.grid.major=element_blank(), panel.border=element_blank(), axis.line=element_line(color='black'), legend.text=element_text(size=16), legend.title=element_text(size=16), legend.position='bottom', legend.box='horizontal', legend.direction='vertical') + scale_y_continuous(expand=c(0,.05))
ggsave('Figures/August2020/compareNetworks_significance_color_repressionAndSize_repressionCasualtyWeight.jpg', plot=last_plot())


#### For realistic fatality rates
sf_short <- subset(sf, repression_rate <= .0575)
sworld_short <- subset(sworld, repression_rate <= .0575)
hk_short <- subset(hk, repression_rate <= .0575)

results_sf <- makeResults(data=sf_short, sampleSizes=n, trials=1000, regressionModel=m)
results_sworld <- makeResults(data=sworld_short, sampleSizes=n, trials=1000, regressionModel=m)
results_hk <- makeResults(data=hk_short, sampleSizes=n, trials=1000, regressionModel=m)

toPlot_sf <- processResults(data=results_sf, tscore=1.96, networkType='3wayCompare_hist_scaleFree_repression0575')
toPlot_sworld <- processResults(data=results_sworld, tscore=1.96, networkType='3wayCompare_hist_smallWorld_repression0575')
toPlot_hk <- processResults(data=results_hk, tscore=1.96, networkType='3wayCompare_hist_holmeKim_repression0575')

# Read below to save time
toPlot_sf <- read.csv('Data/NetworkSimulation/August2020/processedData/3wayCompare_hist_scaleFree_repression0575_regressionData_toplot.csv')
toPlot_sworld <- read.csv('Data/NetworkSimulation/August2020/processedData/3wayCompare_hist_smallWorld_repression0575_regressionData_toplot.csv')
toPlot_hk <- read.csv('Data/NetworkSimulation/August2020/processedData/3wayCompare_hist_holmeKim_repression0575_regressionData_toplot.csv')

toPlot_sf$graph_type <- 'Barabasi-Albert'
toPlot_sworld$graph_type <- 'Watts-Strogatz'
toPlot_hk$graph_type <- 'Holme-Kim'

# Combine
toPlot <- rbind(toPlot_sf, toPlot_sworld, toPlot_hk)

# Plot
ggplot(toPlot, aes(x=sampleSize, y=value, linetype=Significance, color=graph_type)) + geom_line() + theme_bw() +  xlab('Sample Size') + ylab('Percent of Regressions') + theme(axis.text=element_text(size=16), axis.title=element_text(size=16), panel.grid.minor = element_blank(),panel.grid.major=element_blank(), panel.border=element_blank(), axis.line=element_line(color='black'), legend.text=element_text(size=16), legend.position='bottom', legend.box='horizontal', legend.direction='vertical') + labs(color='Graph Type') + scale_color_discrete(limits=c('Watts-Strogatz', 'Holme-Kim', 'Barabasi-Albert')) + scale_y_continuous(expand=c(0,.05))
ggsave('Figures/August2020/compareNetworks_significance_color_repressionAndSize_represion0575.jpg', plot=last_plot())


###################################
###################################
###################################
###################################
###################################
### Realistic values
###################################
###################################
###################################
###################################
###################################
### BELOW IS MADE ON MORNING OF 08.14 AS WORKING THROUGH NEW PLOTS FOR PAPER, MAINLY SM

#### MAKE WITH REALISTIC NETWORK STRUCTURE AND NAVCO WEIGHTING
#### USING ORIGINAL REPRESSION STYLE
# SUBSET RAW DATA ON NETWORK VALUES
sf2 <- sf[sf$scaling_parameter==2.3,]
sworld$scaling_parameter <- round(sworld$scaling_parameter, 2)
sworld2 <- sworld[sworld$scaling_parameter == .34,]  # sw, rounded global clustering and p, looked for p where global clustering = .16. Looked at Watts-Strogatz original paper, eyeballed .4.  .34 is closest value
hk2 <- hk[hk$scaling_parameter == .315,]

results_sf <- makeResults(data=sf2, sampleSizes=n, trials=1000, regressionModel=m, weighting=TRUE, weight_col='casualty_weight', networkType='3wayCompare_hist_scaleFree_realisticNetworkValue_NAVCOCasualtyWeight')
results_sworld <- makeResults(data=sworld2, sampleSizes=n, trials=1000, regressionModel=m, weighting=TRUE, weight_col='casualty_weight', networkType='3wayCompare_hist_smallWorld_realisticNetworkValue_NAVCOCasualtyWeight')
results_hk <- makeResults(data=hk2, sampleSizes=n, trials=1000, regressionModel=m, weighting=TRUE, weight_col='casualty_weight', networkType='3wayCompare_hist_holmeKim_realisticNetworkValue_NAVCOCasualtyWeight')

toPlot_sf <- processResults(data=results_sf, tscore=1.96, networkType='3wayCompare_hist_scaleFree_realisticNetworkValue_NAVCOCasualtyWeight')
toPlot_sworld <- processResults(data=results_sworld, tscore=1.96, networkType='3wayCompare_hist_smallWorld_realisticNetworkValue_NAVCOCasualtyWeight')
toPlot_hk <- processResults(data=results_hk, tscore=1.96, networkType='3wayCompare_hist_holmeKim_realisticNetworkValue_NAVCOCasualtyWeight')

toPlot_sf$graph_type <- 'Barabasi-Albert'
toPlot_sworld$graph_type <- 'Watts-Strogatz'
toPlot_hk$graph_type <- 'Holme-Kim'


# Combine
toPlot <- rbind(toPlot_sf, toPlot_sworld, toPlot_hk)

makeNewPlot(data=toPlot, fileCustomizer='realisticNetworkValue_NAVCOCasualtyWeight')


##### MAKE WITH REALISTIC NETWORK WEIGHTS ONLY, SIMULATION REPRESSION AS IS
## Use subset data from above
## Key is weighting=FALSE in function
results_sf <- makeResults(data=sf2, sampleSizes=n, trials=1000, regressionModel=m, weighting=FALSE, weight_col=NULL, networkType='3wayCompare_hist_scaleFree_realisticNetworkValue')
results_sworld <- makeResults(data=sworld2, sampleSizes=n, trials=1000, regressionModel=m, weighting=FALSE, weight_col=NULL, networkType='3wayCompare_hist_smallWorld_realisticNetworkValue')
results_hk <- makeResults(data=hk2, sampleSizes=n, trials=1000, regressionModel=m, weighting=FALSE, weight_col=NULL, networkType='3wayCompare_hist_holmeKim_realisticNetworkValue')

toPlot_sf <- processResults(data=results_sf, tscore=1.96, networkType='3wayCompare_hist_scaleFree_realisticNetworkValue')
toPlot_sworld <- processResults(data=results_sworld, tscore=1.96, networkType='3wayCompare_hist_smallWorld_realisticNetworkValue')
toPlot_hk <- processResults(data=results_hk, tscore=1.96, networkType='3wayCompare_hist_holmeKim_realisticNetworkValue')

toPlot_sf$graph_type <- 'Barabasi-Albert'
toPlot_sworld$graph_type <- 'Watts-Strogatz'
toPlot_hk$graph_type <- 'Holme-Kim'

# Combine
toPlot <- rbind(toPlot_sf, toPlot_sworld, toPlot_hk)

makeNewPlot(data=toPlot, fileCustomizer='realisticNetworkValue')



###################################
###################################
###################################
###################################
### VARYING NETWORK VALUES
###################################
###################################
###################################
###################################


#### MAKE FOR ONLY WHEN SIZE CHANGES
#### USING ORIGINAL REPRESSION STYLE
# SUBSET RAW DATA ON NETWORK VALUES
# UPDATE 08.14.2020: HAVE NOT RUN YET BECAUSE WILL PROBABLY REMOVE BECAUSE WE DO NOT REALLY DISCUSS AND R1 DID NOT MENTION
sf2 <- sf[sf$initial_size != sf$active_nodes,]
sworld2 <- sworld[sworld$initial_size != sworld$active_nodes,]
hk2 <- hk[hk$initial_size != hk$active_nodes,]


## Use subset data from above
## Key is weighting=FALSE in function
results_sf <- makeResults(data=sf2, sampleSizes=n, trials=1000, regressionModel=m, weighting=FALSE, weight_col=NULL, networkType='3wayCompare_hist_scaleFree_sizeChanged')
results_sworld <- makeResults(data=sworld2, sampleSizes=n, trials=1000, regressionModel=m, weighting=FALSE, weight_col=NULL, networkType='3wayCompare_hist_smallWorld_sizeChanged')
results_hk <- makeResults(data=hk2, sampleSizes=n, trials=1000, regressionModel=m, weighting=FALSE, weight_col=NULL, networkType='3wayCompare_hist_holmeKim_sizeChanged')

toPlot_sf <- processResults(data=results_sf, tscore=1.96, networkType='3wayCompare_hist_scaleFree_sizeChanged')
toPlot_sworld <- processResults(data=results_sworld, tscore=1.96, networkType='3wayCompare_hist_smallWorld_sizeChanged')
toPlot_hk <- processResults(data=results_hk, tscore=1.96, networkType='3wayCompare_hist_holmeKim_sizeChanged')

toPlot_sf$graph_type <- 'Barabasi-Albert'
toPlot_sworld$graph_type <- 'Watts-Strogatz'
toPlot_hk$graph_type <- 'Holme-Kim'

# Combine
toPlot <- rbind(toPlot_sf, toPlot_sworld, toPlot_hk)

makeNewPlot(data=toPlot, fileCustomizer='sizeChanged')



#### MAKE FOR WHEN SCALE FREE MOST LIKE HOLME-KIM
#### USING ORIGINAL REPRESSION STYLE
#### MADE IN /Users/Zack/Documents/UCLA/Research/ProtestTheory/Scripts/PolNet_05_Simulations_RepressionFigures_v15.R, AROUND LINE 1180-1220
## NB: Have looked at results on 08.14.2020 at 3:45.  They seem too different from old way of looking at stuff, but I can't figure out if that's wrong right now.  Need to move on, will think as write paper.
results_sf <- read.csv('Data/NetworkSimulation/Summer2019/processedData/3wayCompare_scaleFree_mostSimilar_filteredSimulationData.csv')
results_sworld <- read.csv('Data/NetworkSimulation/Summer2019/processedData/3wayCompare_smallWorld_mostSimilar_filteredSimulationData.csv')
results_hk <- read.csv('Data/NetworkSimulation/Summer2019/processedData/3wayCompare_holmeKim_mostSimilar_filteredSimulationData.csv')

toPlot_sf <- processResults(data=results_sf, tscore=1.96, networkType='3wayCompare_hist_scaleFree_mostSimilar')
toPlot_sworld <- processResults(data=results_sworld, tscore=1.96, networkType='3wayCompare_hist_smallWorld_mostSimilar')
toPlot_hk <- processResults(data=results_hk, tscore=1.96, networkType='3wayCompare_hist_holmeKim_mostSimilar')

toPlot_sf$graph_type <- 'Barabasi-Albert'
toPlot_sworld$graph_type <- 'Watts-Strogatz'
toPlot_hk$graph_type <- 'Holme-Kim'

# Combine
toPlot <- rbind(toPlot_sf, toPlot_sworld, toPlot_hk)

makeNewPlot(data=toPlot, fileCustomizer='mostSimilarV4')


#### MAKE FOR WHEN NETWORK STRUCTURE CONTROL CHANGES
#### ORIGINAL REPRESSION STYLE
#### MADE IN /Users/Zack/Documents/UCLA/Research/ProtestTheory/Scripts/PolNet_05_Simulations_RepressionFigures_v15.R, AROUND LINE 950-980
## As of 08.14.2020 at 4 p.m., the results are being made in R in Terminal.  Below code is just prepping
toPlot_sf <- read.csv('Data/NetworkSimulation/Summer2019/processedData/3wayCompare_scaleFree_networkStructureControl_networkStructureControl_regressionData_toplot.csv', stringsAsFactors=FALSE)
toPlot_sworld <- read.csv('Data/NetworkSimulation/Summer2019/processedData/3wayCompare_smallWorld_networkStructureControl_networkStructureControl_regressionData_toplot.csv', stringsAsFactors=FALSE)
toPlot_hk <- read.csv('Data/NetworkSimulation/Summer2019/processedData/3wayCompare_holmeKim_networkStructureControl_scalingParameter_regressionData_toplot.csv', stringsAsFactors=FALSE)

toPlot_sf$scaling_parameter_max <- toPlot_sf$scaling_parameter_max - 2  # -2 because range for the others is [0-1], this is [2-3]


# toPlot_sf <- processResults(data=toPlot_sf, tscore=1.96, networkType='3wayCompare_hist_scaleFree_mostSimilar')
# toPlot_sworld <- processResults(data=toPlot_sworld, tscore=1.96, networkType='3wayCompare_hist_smallWorld_mostSimilar')
# toPlot_hk <- processResults(data=toPlot_hk, tscore=1.96, networkType='3wayCompare_hist_holmeKim_mostSimilar')

toPlot_sf$graph_type <- 'Barabasi-Albert'
toPlot_sworld$graph_type <- 'Watts-Strogatz'
toPlot_hk$graph_type <- 'Holme-Kim'

# Combine
toPlot <- rbind(toPlot_sf, toPlot_sworld, toPlot_hk)
## Standardize rates so can be on own plot


makeNewPlot(data=toPlot, fileCustomizer='networkStructureControl')

ggplot(toPlot, aes(x=scaling_parameter_max, y=value, linetype=Significance, color=graph_type)) + geom_line() + facet_wrap(~sampleSize, ncol=3, dir='h') + theme_classic() +  xlab('Network Structure Control') + ylab('Percent of Regressions') + theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text=element_text(size=16), axis.title=element_text(size=16), legend.text=element_text(size=14), legend.title = element_text(size=14), strip.text=element_text(size=16), legend.position='bottom', legend.box='vertical') + labs(color='Graph Type') + scale_color_discrete(limits=c('Watts-Strogatz', 'Holme-Kim', 'Barabasi-Albert'))
ggsave('Figures/August2020/compareNetworks_significance_facet3_color_networkStructureControl.jpg', plot=last_plot())





####################################################
####################################################
####################################################
##
##	NETWORK EXPERIMENTS, RESAMPLE.  REVISION REPRESSION
##
####################################################
####################################################
####################################################
all <- read.csv('Data/NetworkSimulation/Summer2020/processedData/01_experimentsCombined_eigenCore_newRepression.csv', stringsAsFactors=FALSE)
all$normalized <- round(1000*all$active_nodes/all$total_nodes)

# Make names to match model m
all$scaling_parameter <- all$networkStructureControl

# Separate out the networks
sf <- all[all$graph_type=='scale_free_graph',]
sworld <- all[all$graph_type=='watts_strogatz_graph',]
hk <- all[all$graph_type=='powerlaw_cluster_graph',]



# Using full range of repression
n <- c(100, 250, 500, 1000, 1500, 2000)
for(size in n){
	sf_sampling <- pipeline(data=sf, model=m, percentage = size/nrow(sf), trials=1000)
	sworld_sampling <- pipeline(data=sworld, model=m, percentage = size/nrow(sworld), trials=1000)
	hk_sampling <- pipeline(data=hk, model=m, percentage=size/nrow(hk), trials=1000)

	sf_sampling$graph_type <- 'Barabasi-Albert'
	sworld_sampling$graph_type <- 'Watts-Strogatz'
	hk_sampling$graph_type <- 'Holme-Kim'

	# combine
	toPlot <- rbind(sf_sampling, sworld_sampling, hk_sampling)
	makeHist_gg(data=toPlot, filecustomizer=paste0('newRepression_experiments_sampleSize', size), ymax=150)
}


# Using realistic range of repression, 5.75% fatality rate from NAVCO
for(size in n){
	sf_short <- subset(sf, repression_rate <= .0575)
	sworld_short <- subset(sworld, repression_rate <= .0575)
	hk_short <- subset(hk, repression_rate <= .0575)

	sf_sampling <- pipeline(data=sf_short, model=m, percentage = size/nrow(sf_short), trials=1000)
	sworld_sampling <- pipeline(data=sworld_short, model=m, percentage = size/nrow(sworld_short), trials=1000)
	hk_sampling <- pipeline(data=hk_short, model=m, percentage=size/nrow(hk_short), trials=1000)

	sf_sampling$graph_type <- 'Barabasi-Albert'
	sworld_sampling$graph_type <- 'Watts-Strogatz'
	hk_sampling$graph_type <- 'Holme-Kim'

	# combine
	toPlot <- rbind(sf_sampling, sworld_sampling, hk_sampling)
	makeHist_gg(data=toPlot, filecustomizer=paste0('newRepression_experiments_sampleSize', size, '_repression0575Perc'), ymax=150)
}


# Using realistic range of repression, 10.05% casualty rate from NAVCO
for(size in n){
	sf_short <- subset(sf, repression_rate <= .105)
	sworld_short <- subset(sworld, repression_rate <= .105)
	hk_short <- subset(hk, repression_rate <= .105)

	sf_sampling <- pipeline(data=sf_short, model=m, percentage = size/nrow(sf_short), trials=1000)
	sworld_sampling <- pipeline(data=sworld_short, model=m, percentage = size/nrow(sworld_short), trials=1000)
	hk_sampling <- pipeline(data=hk_short, model=m, percentage=size/nrow(hk_short), trials=1000)

	sf_sampling$graph_type <- 'Barabasi-Albert'
	sworld_sampling$graph_type <- 'Watts-Strogatz'
	hk_sampling$graph_type <- 'Holme-Kim'

	# combine
	toPlot <- rbind(sf_sampling, sworld_sampling, hk_sampling)
	makeHist_gg(data=toPlot, filecustomizer=paste0('newRepression_experiments_sampleSize', size, '_casualty1050Perc'), ymax=150)
}


#### Sampling using weights based on NAVCO fatality data


#### Sampling using weights based on NAVCO casualty data


### Sequence of sample sizes
n <- seq(from=50, to=2000, by=50)

#### For all repression rates
# results_sf <- makeResults(data=sf, sampleSizes=n, trials=1000, regressionModel=m, weighting=FALSE, weight_col=NULL, networkType='3wayCompare_hist_scaleFree_newRepression')
# results_sworld <- makeResults(data=sworld, sampleSizes=n, trials=1000, regressionModel=m, weighting=FALSE, weight_col=NULL, networkType='3wayCompare_hist_smallWorld_newRepression')
# results_hk <- makeResults(data=hk, sampleSizes=n, trials=1000, regressionModel=m, weighting=FALSE, weight_col=NULL, networkType='3wayCompare_hist_holmeKim_newRepression')

# toPlot_sf <- processResults(data=results_sf, tscore=1.96, networkType='3wayCompare_hist_scaleFree_newRepression')
# toPlot_sworld <- processResults(data=results_sworld, tscore=1.96, networkType='3wayCompare_hist_smallWorld_newRepression')
# toPlot_hk <- processResults(data=results_hk, tscore=1.96, networkType='3wayCompare_hist_holmeKim_newRepression')

toPlot_sf <- read.csv('Data/NetworkSimulation/August2020/processedData/3wayCompare_hist_scaleFree_newRepression_regressionData_toplot.csv')
toPlot_sworld <- read.csv('Data/NetworkSimulation/August2020/processedData/3wayCompare_hist_smallWorld_newRepression_regressionData_toplot.csv')
toPlot_hk <- read.csv('Data/NetworkSimulation/August2020/processedData/3wayCompare_hist_holmeKim_newRepression_regressionData_toplot.csv')

toPlot_sf$graph_type <- 'Barabasi-Albert'
toPlot_sworld$graph_type <- 'Watts-Strogatz'
toPlot_hk$graph_type <- 'Holme-Kim'

# Combine
toPlot <- rbind(toPlot_sf, toPlot_sworld, toPlot_hk)

# Plot
makeNewPlot(data=toPlot, fileCustomizer='newRepression')



#### For realistic casualty rates
## As of 08.13.2020, have not run
sf_short <- subset(sf, repression_rate <= .105)
sworld_short <- subset(sworld, repression_rate <= .105)
hk_short <- subset(hk, repression_rate <= .105)

results_sf <- makeResults(data=sf_short, sampleSizes=n, trials=1000, regressionModel=m)
results_sworld <- makeResults(data=sworld_short, sampleSizes=n, trials=1000, regressionModel=m)
results_hk <- makeResults(data=hk_short, sampleSizes=n, trials=1000, regressionModel=m)

toPlot_sf <- processResults(data=results_sf, tscore=1.96, networkType='3wayCompare_hist_scaleFree_newRepression_casualtyRate')
toPlot_sworld <- processResults(data=results_sworld, tscore=1.96, networkType='3wayCompare_hist_smallWorld_newRepression_casualtyRate')
toPlot_hk <- processResults(data=results_hk, tscore=1.96, networkType='3wayCompare_hist_holmeKim_newRepression_casualtyRate')

toPlot_sf$graph_type <- 'Barabasi-Albert'
toPlot_sworld$graph_type <- 'Watts-Strogatz'
toPlot_hk$graph_type <- 'Holme-Kim'

# Combine
toPlot <- rbind(toPlot_sf, toPlot_sworld, toPlot_hk)

# Plot
ggplot(toPlot, aes(x=sampleSize, y=value, linetype=Significance, color=graph_type)) + geom_line() + theme_bw() +  xlab('Sample Size') + ylab('Percent of Regressions') + labs(color='Graph Type') + scale_color_discrete(limits=c('Watts-Strogatz', 'Holme-Kim', 'Barabasi-Albert')) + theme(axis.text=element_text(size=16), axis.title=element_text(size=16), panel.grid.minor = element_blank(),panel.grid.major=element_blank(), panel.border=element_blank(), axis.line=element_line(color='black'), legend.text=element_text(size=16), legend.position='bottom', legend.box='horizontal', legend.direction='vertical') + scale_y_continuous(expand=c(0,.05))
ggsave('Figures/August2020/compareNetworks_significance_color_repressionAndSize_repressionCasualtyWeight_newRepression.jpg', plot=last_plot())




### For casualty rates, weighted by NAVCO
results_sf <- makeResults(data=sf, sampleSizes=n, trials=1000, regressionModel=m, weighting=TRUE, weight_col='casualty_weight')
results_sworld <- makeResults(data=sworld, sampleSizes=n, trials=1000, regressionModel=m, weighting=TRUE, weight_col='casualty_weight')
results_hk <- makeResults(data=hk, sampleSizes=n, trials=1000, regressionModel=m, weighting=TRUE, weight_col='casualty_weight')

toPlot_sf <- processResults(data=results_sf, tscore=1.96, networkType='3wayCompare_hist_scaleFree_repressionCasualtyWeight')
toPlot_sworld <- processResults(data=results_sworld, tscore=1.96, networkType='3wayCompare_hist_smallWorld_repressionCasualtyWeight')
toPlot_hk <- processResults(data=results_hk, tscore=1.96, networkType='3wayCompare_hist_holmeKim_repressionCasualtyWeight')

# Read below to save time
# toPlot_sf <- read.csv('Data/NetworkSimulation/August2020/processedData/3wayCompare_hist_scaleFree_repressionCasualtyWeight_regressionData_toplot.csv')
# toPlot_sworld <- read.csv('Data/NetworkSimulation/August2020/processedData/3wayCompare_hist_smallWorld_repressionCasualtyWeight_regressionData_toplot.csv')
# toPlot_hk <- read.csv('Data/NetworkSimulation/August2020/processedData/3wayCompare_hist_holmeKim_repressionCasualtyWeight_regressionData_toplot.csv')

toPlot_sf$graph_type <- 'Barabasi-Albert'
toPlot_sworld$graph_type <- 'Watts-Strogatz'
toPlot_hk$graph_type <- 'Holme-Kim'

# Combine
toPlot <- rbind(toPlot_sf, toPlot_sworld, toPlot_hk)

# Plot
ggplot(toPlot, aes(x=sampleSize, y=value, linetype=Significance, color=graph_type)) + geom_line() + theme_bw() +  xlab('Sample Size') + ylab('Percent of Regressions') + labs(color='Graph Type') + scale_color_discrete(limits=c('Watts-Strogatz', 'Holme-Kim', 'Barabasi-Albert')) + theme(axis.text=element_text(size=16), axis.title=element_text(size=16), panel.grid.minor = element_blank(),panel.grid.major=element_blank(), panel.border=element_blank(), axis.line=element_line(color='black'), legend.text=element_text(size=16)) + scale_y_continuous(expand=c(0,.05))
ggsave('Figures/August2020/compareNetworks_significance_color_repressionAndSize_represionCasualtyWeight.jpg', plot=last_plot())




#### For realistic fatality rates
sf_short <- subset(sf, repression_rate <= .0575)
sworld_short <- subset(sworld, repression_rate <= .0575)
hk_short <- subset(hk, repression_rate <= .0575)

results_sf <- makeResults(data=sf_short, sampleSizes=n, trials=1000, regressionModel=m)
results_sworld <- makeResults(data=sworld_short, sampleSizes=n, trials=1000, regressionModel=m)
results_hk <- makeResults(data=hk_short, sampleSizes=n, trials=1000, regressionModel=m)

toPlot_sf <- processResults(data=results_sf, tscore=1.96, networkType='3wayCompare_hist_scaleFree_represion0575')
toPlot_sworld <- processResults(data=results_sworld, tscore=1.96, networkType='3wayCompare_hist_smallWorld_represion0575')
toPlot_hk <- processResults(data=results_hk, tscore=1.96, networkType='3wayCompare_hist_holmeKim_represion0575')

toPlot_sf$graph_type <- 'Barabasi-Albert'
toPlot_sworld$graph_type <- 'Watts-Strogatz'
toPlot_hk$graph_type <- 'Holme-Kim'

# Combine
toPlot <- rbind(toPlot_sf, toPlot_sworld, toPlot_hk)

# Plot
ggplot(toPlot, aes(x=sampleSize, y=value, linetype=Significance, color=graph_type)) + geom_line() + theme_bw() +  xlab('Sample Size') + ylab('Percent of Regressions') + labs(color='Graph Type') + scale_color_discrete(limits=c('Watts-Strogatz', 'Holme-Kim', 'Barabasi-Albert')) + theme(axis.text=element_text(size=16), axis.title=element_text(size=16), panel.grid.minor = element_blank(),panel.grid.major=element_blank(), panel.border=element_blank(), axis.line=element_line(color='black'), legend.text=element_text(size=16), legend.position='bottom', legend.box='horizontal', legend.direction='vertical') + scale_y_continuous(expand=c(0,.05))
ggsave('Figures/August2020/compareNetworks_significance_color_repressionAndSize_represion0575_newRepression.jpg', plot=last_plot())



########
# Include thresholds in regression
########
m_threshold <- log(normalized+1, 10) ~ initial_density + initial_global_clustering + initial_mean_degree + initial_median_degree + initial_neighborhood_clustering + initial_nodes_clustering + initial_size + scaling_parameter + repression_rate + initial_mean_threshold + initial_neighbors_mean_threshold


results_sf <- makeResults(data=sf, sampleSizes=n, trials=1000, regressionModel=m_threshold, weighting=FALSE, weight_col=NULL, networkType='3wayCompare_hist_scaleFree_newRepression_thresholds')
results_sworld <- makeResults(data=sworld, sampleSizes=n, trials=1000, regressionModel=m_threshold, weighting=FALSE, weight_col=NULL, networkType='3wayCompare_hist_smallWorld_newRepression_thresholds')
results_hk <- makeResults(data=hk, sampleSizes=n, trials=1000, regressionModel=m_threshold, weighting=FALSE, weight_col=NULL, networkType='3wayCompare_hist_holmeKim_newRepression_thresholds')

toPlot_sf <- processResults(data=results_sf, tscore=1.96, networkType='3wayCompare_hist_scaleFree_newRepression_thresholds')
toPlot_sworld <- processResults(data=results_sworld, tscore=1.96, networkType='3wayCompare_hist_smallWorld_newRepression_thresholds')
toPlot_hk <- processResults(data=results_hk, tscore=1.96, networkType='3wayCompare_hist_holmeKim_newRepression_thresholds')

# toPlot_sf <- read.csv('Data/NetworkSimulation/August2020/processedData/3wayCompare_hist_scaleFree_newRepression_thresholds_regressionData_toplot.csv')
# toPlot_sworld <- read.csv('Data/NetworkSimulation/August2020/processedData/3wayCompare_hist_smallWorld_newRepression_thresholds_regressionData_toplot.csv')
# toPlot_hk <- read.csv('Data/NetworkSimulation/August2020/processedData/3wayCompare_hist_holmeKim_newRepression_thresholds_regressionData_toplot.csv')

toPlot_sf$graph_type <- 'Barabasi-Albert'
toPlot_sworld$graph_type <- 'Watts-Strogatz'
toPlot_hk$graph_type <- 'Holme-Kim'

# Combine
toPlot <- rbind(toPlot_sf, toPlot_sworld, toPlot_hk)

# Plot
makeNewPlot(data=toPlot, fileCustomizer='newRepression_thresholdModel')


### Make histograms with threshold variables
# Using full range of repression
n <- c(100, 250, 500, 1000, 1500, 2000)
for(size in n){
	sf_sampling <- pipeline(data=sf, model=m_threshold, percentage = size/nrow(sf), trials=1000)
	sworld_sampling <- pipeline(data=sworld, model=m_threshold, percentage = size/nrow(sworld), trials=1000)
	hk_sampling <- pipeline(data=hk, model=m_threshold, percentage=size/nrow(hk), trials=1000)

	sf_sampling$graph_type <- 'Barabasi-Albert'
	sworld_sampling$graph_type <- 'Watts-Strogatz'
	hk_sampling$graph_type <- 'Holme-Kim'

	# combine
	toPlot <- rbind(sf_sampling, sworld_sampling, hk_sampling)
	makeHist_gg(data=toPlot, filecustomizer=paste0('newRepression_experiments_sampleSize', size, '_thresholds'), ymax=150, thresholds=TRUE)
}


####################################################
####################################################
####################################################
##
##	NETWORK EXPERIMENTS, RESAMPLE.  EDGE REPRESSION
##
####################################################
####################################################
####################################################
all <- read.csv('Data/NetworkSimulation/Summer2018/processedData/01_experimentsCombined_edgeRepression.csv', stringsAsFactors=FALSE)
all$normalized <- round(1000*all$final_size/all$total_nodes)

# Make names to match model m
all$scaling_parameter <- all$networkStructureControl

# Separate out the networks
sf <- all[all$graph_type=='scale_free_graph',]
sworld <- all[all$graph_type=='watts_strogatz_graph',]
hk <- all[all$graph_type=='powerlaw_cluster_graph',]

### Sequence of sample sizes
n <- seq(from=50, to=2000, by=50)

results_sf <- makeResults(data=sf, sampleSizes=n, trials=1000, regressionModel=m, weighting=FALSE, weight_col=NULL, networkType='3wayCompare_scaleFree_edgeRepression')
results_sworld <- makeResults(data=sworld, sampleSizes=n, trials=1000, regressionModel=m, weighting=FALSE, weight_col=NULL, networkType='3wayCompare_smallWorld_edgeRepression')
results_hk <- makeResults(data=hk, sampleSizes=n, trials=1000, regressionModel=m, weighting=FALSE, weight_col=NULL, networkType='3wayCompare_holmeKim_edgeRepression')


toPlot_sf <- processResults(data=results_sf, tscore=1.96, networkType='3wayCompare_scaleFree_edgeRepression')
toPlot_sworld <- processResults(data=results_sworld, tscore=1.96, networkType='3wayCompare_smallWorld_edgeRepression')
toPlot_hk <- processResults(data=results_hk, tscore=1.96, networkType='3wayCompare_holmeKim_edgeRepression')

toPlot_sf$graph_type <- 'Barabasi-Albert'
toPlot_sworld$graph_type <- 'Watts-Strogatz'
toPlot_hk$graph_type <- 'Holme-Kim'

# Combine
toPlot <- rbind(toPlot_sf, toPlot_sworld, toPlot_hk)


makeNewPlot(data=toPlot, fileCustomizer='edgeRepression')


