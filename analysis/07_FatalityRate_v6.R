'''
The purpose of this script is to generate an estimate of the fatality rate at protests.  It uses ACLED and NAVCO.

v2. Make tables.
v3. Add vertical lines for averages
v4. Looking at again during R&R for Network Science
	- Use newest ACLED - Too much of a pain to get numbers, not worth the squeeze.
	- Use newest NAVCO - Nothing new
v5. Show NAVCO in new way - x is size, y is fatality rate, show hline as average
'''
library(stringr)
library(ggplot2)
library(reshape)

set.seed(480193)
setwd('/Users/Zack/Documents/UCLA/Research/ProtestTheory')

#### ACLED
acled <- read.csv('Data/ACLED/ACLED_protests_1900-01-01-2020-08-10.csv', stringsAsFactors=FALSE)

mean(acled$fatalities)  # .016659
mean(acled$fatalities[acled$fatalities > 0])  # 3.9065

acled$size <- str_extract_all(acled$notes, "(?<=\\[).+?(?=\\])")
acled$size <- sapply(acled$size, function(x) x[1])

#### FRANCISCO MASSACRES
massacre <- read.csv('/Users/Zack/Documents/UCLA/Data/Massacres/Massacres_Table1.csv', stringsAsFactors=FALSE)
massacre$fatal_rate <- massacre$Fatality/massacre$Size
massacre$casualty_rate <- (massacre$Fatality+massacre$Injury)/massacre$Size

cumsum(prop.table(table(round(massacre$fatal_rate, 2))))
cumsum(prop.table(table(round(massacre$casualty_rate, 2))))


#### NAVCO
navco <- read.csv('Data/NAVCO/navco3.csv', stringsAsFactors=FALSE)
navco <- subset(navco, verb_10==14)
navco <- navco[is.na(navco$num_partic_event) == FALSE,]
navco$partic_split <- grepl('-', navco$num_partic_event)  # For simplicity, ignore all where estimate has a dash, is a range.


navco2 <- navco[navco$partic_split == FALSE,]  # Keep only precise ones
navco2$count <- as.numeric(navco2$num_partic_event)
navco2 <- navco2[is.na(navco2$count) == FALSE,]

navco2$fatal_casualty <- as.numeric(navco2$fatal_casu)
navco2$injuries <- as.numeric(navco2$injuries)

mean(navco2$fatal_casualty/navco2$count, na.rm=TRUE)  # 5.4%
mean((navco2$fatal_casualty+navco2$injuries)/navco2$count, na.rm=TRUE)  # 10.05%

quantile(navco2$fatal_casualty/navco2$count,probs=seq(0,1,.05),na.rm=TRUE)  # 0 through 80th percentile≤ .0033 is 90th percentile
quantile((navco2$fatal_casualty+navco2$injuries)/navco2$count,probs=seq(0,1,.05),na.rm=TRUE)  # 0 through 75th percentile, .025 at 85th percentile


#### NAVCO, COUNT NA AS 0 FOR CASUALTIES
navco2$fatal_casualty2 <- navco2$fatal_casualty
navco2$fatal_casualty2[is.na(navco2$fatal_casualty2)==TRUE] <- 0

navco2$injuries2 <- navco2$injuries
navco2$injuries2[is.na(navco2$injuries2)==TRUE] <- 0

mean(navco2$fatal_casualty2/navco2$count)  # 1.06%
mean((navco2$fatal_casualty2+navco2$injuries2)/navco2$count)  # 1.81%


#### MAKE HISTOGRAM
navco2$fatal_rate <- navco2$fatal_casualty/navco2$count
navco2$injury_rate <- navco2$injuries/navco2$count
navco2$casualty_rate <- (navco2$fatal_casualty + navco2$injuries)/navco2$count


# Tables
table(round(navco2$fatal_rate,2))
prop.table(table(round(navco2$fatal_rate,2)))
cumsum(prop.table(table(round(navco2$fatal_rate,2))))

table(round(navco2$injury_rate,2))
prop.table(table(round(navco2$injury_rate,2)))
cumsum(prop.table(table(round(navco2$injury_rate,2))))

table(round(navco2$casualty_rate,2))
prop.table(table(round(navco2$casualty_rate,2)))
cumsum(prop.table(table(round(navco2$casualty_rate,2))))

## MELT so can make two lines on one plot
use <- subset(navco2, navco2$fatal_rate <= 1)  # 1 observation had a fatality rate of 2
use <- subset(use, use$casualty_rate <= 1)  # same as above
hope <- melt(use, id.vars='count' ,measure.vars=c('fatal_rate', 'casualty_rate'))
hope$count_norm <- hope$count/max(hope$count)

ggplot(hope, aes(x=log(count, 10), y=value, color=variable)) + geom_point(alpha=.5) + ylab('Rate') + xlab('Protest Size') + theme_classic() + geom_smooth(method='loess', se = FALSE) + scale_x_continuous(limits = c(0, 6), breaks = seq(0,6,1), labels = c("1", "10", "100", "1000", "10000", "100000", "1000000")) + scale_y_continuous(expand=c(.02, .02)) + theme(axis.text.x=element_text(angle=45, hjust=1), axis.text=element_text(size=16), axis.title=element_text(size=16), legend.position=c(.8, .95), legend.title=element_blank(), legend.text=element_text(size=16)) + scale_color_manual(values=c('red', 'black'), breaks=c('fatal_rate', 'casualty_rate'), labels=c('Fatality Rate', 'Casualty Rate'))
ggsave('Figures/August2020/navco_fatalityCasualty.jpeg', plot=last_plot())


## Get weights for repression so can use to sample simulations
fatalRates <- data.frame(prop.table(table(round_any(use$fatal_rate,.005))))
names(fatalRates) <- c('RepressionRate', 'Percentage')  # Can do .9 for rate = 0, .1 for all else

casualtyRates <- data.frame(prop.table(table(round_any(use$casualty_rate,.005))))  # .8 for 0 repression rate, .02 for .005, .02 for .01; .16 for all else rates
names(casualtyRates) <- c('RepressionRate', 'Percentage	')


## FATALITIES
temp <- subset(navco2, navco2$fatal_rate <= 1)  # 1 observation had a fatality rate of 2
histData <- hist(temp$fatal_rate, plot=F, breaks=100)
# histData$counts[histData$counts==1] <- 10
# histData$counts[histData$counts==0] <- 1
histData$counts <- log10(histData$counts+1)
yAxis <- c(0,1, 10, 100, 1000, 10000)

jpeg('Figures/June2019/histogram_fatality.jpeg', width=7, height=7, units='in', res=300)
plot(histData, yaxt='n', xlab='Fatality Rate', ylim=c(0.00001,4), main="", cex.axis=1.5, cex.lab=1.5)
axis(2, at=log10(yAxis+1), labels=yAxis, cex.axis=1.5)
abline(v=c(median(temp$fatal_rate), mean(temp$fatal_rate)), lty=c('longdash', 'solid'))
legend('top', legend=c('Median', 'Mean'), lty=c('longdash', 'solid'), box.lty=0, cex=1.5)
dev.off()



## INJURIES
temp <- subset(navco2, navco2$injury_rate <= 1)  # 1 observation had a fatality rate of 2
histData <- hist(temp$injury_rate, plot=F, breaks=100)
# histData$counts[histData$counts==1] <- 10
# histData$counts[histData$counts==0] <- 1
histData$counts <- log10(histData$counts+1)
yAxis <- c(0,1, 10, 100, 1000, 10000)

jpeg('Figures/June2019/histogram_injuries.jpeg', width=7, height=7, units='in', res=300)
plot(histData, yaxt='n', xlab='Injury Rate', ylim=c(0.00001,4), main="", cex.axis=1.5, cex.lab=1.5)
axis(2, at=log10(yAxis+1), labels=yAxis, cex.axis=1.5)
abline(v=c(median(temp$injury_rate), mean(temp$injury_rate)), lty=c('longdash', 'solid'))
legend('top', legend=c('Median', 'Mean'), lty=c('longdash', 'solid'), box.lty=0, cex=1.5)
dev.off()

## CASUALTIES
temp <- subset(navco2, navco2$casualty_rate <= 1)  # 1 observation had a fatality rate of 2
histData <- hist(temp$casualty_rate, plot=F, breaks=100)
# histData$counts[histData$counts==1] <- 10
# histData$counts[histData$counts==0] <- 1
histData$counts <- log10(histData$counts+1)
yAxis <- c(0,1, 10, 100, 1000, 10000)

jpeg('Figures/June2019/histogram_casualties.jpeg', width=7, height=7, units='in', res=300)
plot(histData, yaxt='n', xlab='Casualty Rate', ylim=c(0.00001,4), main="", cex.axis=1.5, cex.lab=1.5)
axis(2, at=log10(yAxis+1), labels=yAxis, cex.axis=1.5)
abline(v=c(median(temp$casualty_rate), mean(temp$casualty_rate)), lty=c('longdash', 'solid'))
legend('top', legend=c('Median', 'Mean'), lty=c('longdash', 'solid'), box.lty=0, cex=1.5)
dev.off()


### WHAT ARE THE EVENTS WHERE EVERYONE IS INJURED OR KILLED?
mean(navco2$count[navco2$fatal_rate == 1], na.rm=TRUE) # 1.295
mean(navco2$count[navco2$injury_rate == 1], na.rm=TRUE)  # 24.04167


### WHAT IS FIRST DATE AND LAST DATE
temp <- navco2[order(navco2$date),]
