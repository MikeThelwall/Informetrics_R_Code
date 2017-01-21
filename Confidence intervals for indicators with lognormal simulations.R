#Discretised lognormal simulations for citation data to test the precision of the MNLCS indicator
#See the end for the location of the file to save the results in
library(SDMTools)	#for wt.mean; wt.sd http://cran.r-project.org/web/packages/SDMTools/SDMTools.pdf
library(poweRlaw)

iterations <- 10000	# Change this for a different number of iterations (less than 100 can cause errors)
bootstrapCount <- 1000 # Change this for a different number of bootstrap iterations (less than 100 can cause errors)
simulationMethod <- 1; # 1 is discrete with powerlaw, 0 is continuous lognormal rounded.

confidenceInterval <- function(x, pc) {
	pc <- pc/2
	sortindex <- order(x, decreasing=TRUE)
	cutoff.index <- (length(x)+1) * pc/100; cutoff.l <- floor(cutoff.index);
	cutoff.index <- (length(x)+1) * (100-pc)/100; cutoff.u <- ceiling(cutoff.index) 
	cutoff.lvalue <- x[sortindex][cutoff.l]  
	cutoff.uvalue <- x[sortindex][cutoff.u]  
	return (c(cutoff.uvalue, cutoff.lvalue))
}

MainFunction <- function() {
cat("mean0", "mean1", "meanRest", "sd0", "sd1", "sampleSizeAll", "subSampleSize1", "subSampleProp1",
	"mean(MNCSSample1)", "MNCSSample1.L", "MNLSSample1.U",  
	"mean(MNLCSSample1)", "MNLCSSample1.L", "MNLCSSample1.U", 
	"MNCSSample1BootMean.Mean", "MNCSSample1BootMean.L", "MNCSSample1BootMean.U",
	"MNLCSSample1BootMean.Mean", "MNLCSSample1BootMean.L", "MNLCSSample1BootMean.U",
	"MNLCSSample1BootMedian.Mean", "MNLCSSample1BootMedian.L", "MNLCSSample1BootMedian.U",
	"MNCSinBootCICount", "MNLCSinFormulaCICount", "MNLCSinBootCICount",
	"mean(MNLCSSample1Formula.L)","mean(MNLCSSample1Formula.U)","mean(MNLCSSample1Boot.L)","mean(MNLCSSample1Boot.U)",
	"\n", sep="\t")
bootMNCS1 <- rep(-999,bootstrapCount);bootMNCS1.L <- rep(-999,bootstrapCount); bootMNCS1.U <- rep(-999,bootstrapCount);
bootMNLCS1 <- rep(-999,bootstrapCount);bootMNLCS1.L <- rep(-999,bootstrapCount); bootMNLCS1.U <- rep(-999,bootstrapCount);
lnRand = dislnorm$new(); lnRand$setXmin(0);
for (mean0 in seq(5,5,by=0.1)) {				# Change this for a different range of means; note that 0.5 can cause problems if mean1 and mean2 are high
	for (sd0 in seq(0.75,1.5,by=0.25)) {				# Change this for a different range of standard deviations (if the 1st and 2nd number are the same it does that number) seq(0.5,2,by=0.5)
		sd1 <- sd0
		for (mean1 in seq(mean0,mean0,by=1)) {				# Change this for a different range of mean1s
			for (sampleSize in seq(10000,10000,by=5000)) {		# Change this for a different range of sample sizes
				for (subSampleProp1 in seq(0.01,0.25,by=0.01)) {			# Change this for a different range of sample sizes 1
					###Initialise results variables###
					topN <- rep(-999,iterations); 
					MNCSSample1 <- topN; MNCSSample1Boot.Mean <- topN; MNCSSample1Boot.Median <- topN; MNCSSample1Boot.L <- topN; MNCSSample1Boot.U <- topN;
					MNLCSSample1 <- topN; MNLCSSample1Formula.L <- topN; MNLCSSample1Formula.U <- topN; MNLCSSample1Boot.Mean <- topN; MNLCSSample1Boot.Median <- topN; MNLCSSample1Boot.L <- topN; MNLCSSample1Boot.U <- topN;
					subSampleSize1 <- round(sampleSize*subSampleProp1)
					###Set up data frame with one vector, inSub1, containing 1s for subSample1, then 0s for subSample2. A data frame is used for storing data tables. It is a list of vectors of equal length.###
					ldata <- data.frame(inSub1 = c(rep(1,subSampleSize1), rep(0,sampleSize-subSampleSize1))) #subSample1 1s, then all 0s
					###Add second vector, outSub, to data frame, which 0 for Sub1 and 1 for the rest###
					ldata <- transform(ldata, outSub = 1 - inSub1) 
					###Add third vector, all, to data frame, which is just 1s ###
					ldata <- transform(ldata, all = outSub + inSub1)
					bdata <- ldata # copy of ldata for bootstr
					###Calculate mean of Out data?###
					meanOut <- log( (exp(mean0) - subSampleProp1*exp(mean1)) / (1 - subSampleProp1) ) 
					###Add fourth vector, Emean, to data frame, which is expected mean for each sample###
					ldata <- transform(ldata, Emean = inSub1 * mean1 + outSub * meanOut)  #Final bit is adjustment to make overall mean constant at mean0, assuming equal sd
					###Add fifth vector, Esd, to data frame, which is just expacted SD?###
					ldata <- transform(ldata, Esd = sd0 + inSub1 * (sd1 - sd0)) 
					###Start the simulations###
					for (i in 1:iterations) {
						###Add sixth, y, to data frame, which is sampled from the lognormal distribution, with appropriate Mean and SD parameters (varying by sample and world)###
						if (simulationMethod == 0) { # Continous lognormal rounded
							ldata <- transform(ldata, y = rlnorm(sampleSize, mean = Emean, sd = Esd)) #random continuous lognormal @@@@@@@@ CHANGE TO DIRECT DISCRETISED LOGNORMAL IF POSSIBLE
							ldata <- transform(ldata, y = round(y))
						} else { # Discrete lognormal package
							lnRand$setPars(c(mean1,sd1)); randSample1 <- dist_rand(lnRand, subSampleSize1)
							lnRand$setPars(c(meanOut,sd0)); randOut <- dist_rand(lnRand, sampleSize - subSampleSize1)
							ldata <- transform(ldata, y = c(randSample1,randOut))							
						}
						ldata <- transform(ldata, yGeo = log(1+y))
						groupMeanLog <- wt.mean(ldata$yGeo,ldata$inSub1); allMeanLog <- wt.mean(ldata$yGeo,ldata$all) 
						MNCSSample1[i] <- wt.mean(ldata$y,ldata$inSub1) / wt.mean(ldata$y,ldata$all)				
						MNLCSSample1[i] <- groupMeanLog / allMeanLog
						## Calculate 95% confidence interval for MNLCS using Fieller formula.
						SEsample1Log <- wt.sd(ldata$yGeo,ldata$inSub1)/sqrt(subSampleSize1); SEallLog <- wt.sd(ldata$yGeo,ldata$all)/sqrt(sampleSize) 
						h <- (1.96*SEallLog/mean(allMeanLog))^2
						SEmnlcs <- MNLCSSample1[i]/(1-h)*sqrt((1-h)*(SEsample1Log/groupMeanLog)^2+(SEallLog/allMeanLog)^2)
						MNLCSSample1Formula.L[i] <- MNLCSSample1[i]/(1-h) - 1.96*SEmnlcs
						MNLCSSample1Formula.U[i] <- MNLCSSample1[i]/(1-h) + 1.96*SEmnlcs
						vecOut <- subset(ldata, outSub==1, select=y)		# http://www.statmethods.net/management/subset.html #newdata <- subset(mydata, sex=="m" & age > 25, select=weight:income)
						vecIn1 <- subset(ldata, inSub1==1, select=y)
						for (j in 1:bootstrapCount) {		# Only need to resample from distribution: don't generate from rlnorm again http://www.statmethods.net/advstats/bootstrapping.html
							vecOutBoot <- sample(vecOut$y,length(vecOut$y), replace = TRUE) #Sample same length vector with replacement
							vecIn1Boot <- sample(vecIn1$y,length(vecIn1$y), replace = TRUE) #Sample same length vector with replacement
							bootMNCS1[j] <- mean(vecIn1Boot) / mean(c(vecIn1Boot,vecOutBoot))
							vecIn1Boot <-log(1+vecIn1Boot); vecOutBoot <-log(1+vecOutBoot)
							bootMNLCS1[j] <- mean(vecIn1Boot) / mean(c(vecIn1Boot,vecOutBoot))
						}
						#Calculate bootstrap confidence intervals, mean and median for MNCS and MNLCS
						ci <- confidenceInterval(bootMNCS1, 5); MNCSSample1Boot.L[i] <- ci[1]; MNCSSample1Boot.U[i] <- ci[2]; MNCSSample1Boot.Mean[i] <- mean(bootMNCS1); MNCSSample1Boot.Median[i] <- median(bootMNCS1)
						ci <- confidenceInterval(bootMNLCS1, 5); MNLCSSample1Boot.L[i] <- ci[1]; MNLCSSample1Boot.U[i] <- ci[2]; MNLCSSample1Boot.Mean[i] <- mean(bootMNLCS1); MNLCSSample1Boot.Median[i] <- median(bootMNLCS1)
					}
					#Got 1000 bootstrapped and formula confidence intervals for MNLCS
					MNCSEstimate <- mean(MNCSSample1)	# Estimate from mean of all simulated samples generated from function
					MNLCSEstimate <- mean(MNLCSSample1)	# Estimate from mean of all simulated samples generated from function
					MNLCSinBootCICount <- 0; MNLCSinFormulaCICount <- 0; ; MNCSinBootCICount <- 0;
					for (i in 1:iterations) {
						if ( (MNCSSample1Boot.L[i] <= MNCSEstimate) && (MNCSEstimate <= MNCSSample1Boot.U[i]) )  MNCSinBootCICount <- MNCSinBootCICount + 1
						if ( (MNLCSSample1Boot.L[i] <= MNLCSEstimate) && (MNLCSEstimate <= MNLCSSample1Boot.U[i]) )  MNLCSinBootCICount <- MNLCSinBootCICount + 1
						if ( (MNLCSSample1Formula.L[i] <= MNLCSEstimate) && (MNLCSEstimate <= MNLCSSample1Formula.U[i]) )  MNLCSinFormulaCICount <- MNLCSinFormulaCICount + 1
					}
					ci <- confidenceInterval(MNCSSample1, 5); MNCSSample1.L <- ci[1]; MNCSSample1.U <- ci[2]
					ci <- confidenceInterval(MNLCSSample1, 5); MNLCSSample1.L <- ci[1]; MNLCSSample1.U <- ci[2]
					ci <- confidenceInterval(MNCSSample1Boot.Mean, 5); MNCSSample1BootMean.L <- ci[1]; MNCSSample1BootMean.U <- ci[2]; MNCSSample1BootMean.Mean <- mean(MNCSSample1Boot.Mean)
					ci <- confidenceInterval(MNLCSSample1Boot.Mean, 5); MNLCSSample1BootMean.L <- ci[1]; MNLCSSample1BootMean.U <- ci[2]; MNLCSSample1BootMean.Mean <- mean(MNLCSSample1Boot.Mean)
					ci <- confidenceInterval(MNLCSSample1Boot.Median, 5); MNLCSSample1BootMedian.L <- ci[1]; MNLCSSample1BootMedian.U <- ci[2]; MNLCSSample1BootMedian.Mean <- mean(MNCSSample1Boot.Median)
					cat(mean0, mean1, meanOut, sd0, sd1, sampleSize, subSampleSize1, subSampleProp1, 
							mean(MNCSSample1), MNCSSample1.L, MNCSSample1.U,
							MNLCSEstimate, MNLCSSample1.L, MNLCSSample1.U,
							MNCSSample1BootMean.Mean, MNCSSample1BootMean.L, MNCSSample1BootMean.U,
							MNLCSSample1BootMean.Mean, MNLCSSample1BootMean.L, MNLCSSample1BootMean.U,
							MNLCSSample1BootMedian.Mean, MNLCSSample1BootMedian.L, MNLCSSample1BootMedian.U,
							MNCSinBootCICount, MNLCSinFormulaCICount, MNLCSinBootCICount,
							mean(MNLCSSample1Formula.L),mean(MNLCSSample1Formula.U),mean(MNLCSSample1Boot.L),mean(MNLCSSample1Boot.U), 
						 	"\n", sep="\t")
				}
			}
		}
	}
}
}#MainFunction()

sink(paste("C:/results/ConfidenceIntervalFormulaAccuracyFixedH 6 10k var p 5 ",runif(1,0,1)," dis", simulationMethod, ".txt", sep="")) 
MainFunction()
sink()

