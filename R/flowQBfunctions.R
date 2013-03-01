
# Read file and Doublet Discriminations

ReadDD <- function(File, Chan1DD, Chan2DD, P, ChanGiven, ChanCompanion)
{
  # File: the path with the name of the file.
  # Chan1DD and Chan2DD: the two channel indices to be used for detecting the doublet events.
  # P: the doublet discrimination level (e.g. 95%).
  # ChanGiven:the index of the channel under study.
  # ChanCompanion: the index of the side scattering channel.
  
  # The function returns singlet positive fluorescent log intensities of ChanGiven and ChanCompanion channels.

  # Needed library flowCore to read a FCS file.
  library(flowCore)
  # READ THE FILE
  RAWDATA <- read.FCS(File, transformation = FALSE)
  # get expression
  ExpData <- data.frame(exprs(RAWDATA))
  # Extract the raw data
  RDChan <- data.frame(ExpData[, ChanGiven], ExpData[, ChanCompanion])
  # Only positive fluorescent intensities and singlets to be considered
  # Make P as a ratio
  P1 <- P/100
  P2 <- 2-P1
  RatioDD <- ExpData[, Chan1DD]/(1+ExpData[, Chan2DD])
  # Singlets with positive fluorescent intensities
  PosEv <- which(ExpData[,ChanGiven] > 0 & ExpData[,ChanCompanion] > 0 & P1 <= RatioDD & RatioDD < P2)
  # Extract the logarithm of the positive fluorescent intensities
  # of the founded singlets and return the data as a 2D array
  LogData <- data.frame(ChanGiven=log(ExpData[PosEv, ChanGiven]), ChanCompanion=log(ExpData[PosEv, ChanCompanion]))
  return(LogData)
}


# K-means Means and SDs

KmeansMeanSD <- function(transformed2Darray, nClusters, nstart, itermax, Vis)
{   

    # transformed2Darray: 2D data to cluster.
    # nClusters:  number of clusters. 
    # nstart: if centers is a number, how many random.
    # itermax: the maximum number of iterations allowed.
    # sets should be chosen?  
    # For more details, see kmeans function in R.
    # Only if  Vis equals one, the code plots the 2D Kmeans results.

# Output for each cluster the statistics: Mean and standard deviation (SD) following a percentile strategy  to exclude events out of 2% and 99% :

probs <- c(2,99)/100

if(nClusters >= 1) 
{
  cl <- kmeans(transformed2Darray, nClusters, iter.max = itermax, nstart = nstart)
  # Only if Vis equals one, the code plots the 2D Kmeans results
  if(Vis==1)   
  {
    plot(transformed2Darray, xlab="Given Channel", "Companion Channel", col = cl$cluster+1) 
    points(cl$centers, col = "black" , pch = 16, cex=1.75)
  }

  SDs <- rep(0, nClusters)
  Means <- rep(0, nClusters)

  for(jc in 1:nClusters)
  {
    Events <- which(cl$cluster==jc)
    vect <- transformed2Darray[Events, 1]
    percentiles <- as.double(quantile(vect, probs))
    te <- Events[which( percentiles[1] <= vect & vect <= percentiles[2])]
    if(length(te) > 0)  Means[jc] <- mean(exp(transformed2Darray[te, 1]))
    if(length(te) > 0)  SDs[jc] <- sd(exp(transformed2Darray[te, 1]))
    if(length(te) == 0) Means[jc] <- exp(cl$centers[jc, 1])
    if(length(te) == 0) SDs[jc] <- 0.0
  }
  # Order the means from lower intensity to higher intensity 
  OrderedMeans <- order(Means)
  return(data.frame(Means=Means[OrderedMeans], SDs=SDs[OrderedMeans]))
}
}


# Median and Robust standard deviation

KmeansMedianrSD <- function(transformed2Darray, nClusters, nstart, itermax, Vis)
{   

	# transformed2Darray: 2D data to cluster.
	# nClusters:  number of clusters. 
	# nstart: if centers is a number, how many random.
        # itermax: the maximum number of iterations allowed.
	# sets should be chosen?  
	# For more details, see kmeans function in R.
        # Only if  Vis equals one, the code plots the 2D Kmeans results.

	# Output for each cluster the robust statistics: Median and Robust standard deviation (rSD) 
        # following CS&T proposition:
	# rSD = (Median of{|Xi − Medianx |} × 1.4826.


if(nClusters >= 1) 
{
  cl <- kmeans(transformed2Darray, nClusters, iter.max = itermax, nstart = nstart)
  # Only if Vis equals one, the code plots the 2D Kmeans results
  if(Vis == 1)   
  {
    plot(transformed2Darray, xlab="Given Channel", "Companion Channel", col = cl$cluster+1) 
    points(cl$centers, col = "black" , pch = 16, cex=1.75)
  }

  rSDs <- rep(0, nClusters)
  Medians <- rep(0, nClusters)

  for(jc in 1:nClusters)
  {
    Events <- which(cl$cluster == jc)
    medianx=median(exp(transformed2Darray[Events, 1]))
    rSDs[jc] <- median(abs(medianx - exp(transformed2Darray[Events, 1]))) *1.4826
    Medians[jc] <- medianx
}


  # Order the means from lower intensity to higher intensity 
  OrderedMeans <- order(Medians)
  return(data.frame(Means=Medians[OrderedMeans], rSDs=rSDs[OrderedMeans]))
}
}




# Calibrartion
MFI2MESF <- function(MeansSDs, p, IllCorrCV)
{   

  # MeansSDs: 2D array with MFI Means and SDs.
  # p: the value used to convert MFI to MESF, MESF = p*MFI
  # Function IllCorrCV() to be used for Illumination-Correction. 
  # The SDs are corrected for Illumination-Correction.
  # If IllCorrCV is equal to zero no correction is processed.

  MESFMEAN <- MeansSDs[,1]*p
  # Compute  CVs
  CV <- MeansSDs[, 2]/(1 + MeansSDs[, 1])
  # SD^2 Illumination-Correction
  MESFV <- (MESFMEAN^2)*(CV^2 - IllCorrCV^2)
  return(data.frame(MESF=MESFMEAN, MESFV=MESFV))
}


# Linear Regression

lrMESF <- function(MESFmeanssd2, Peak1, Peak2)
{  

  #  MESFmeanssd2: 2D array having the means and thier variance.
  #  From Peak1: index of the first peak to be used consequtively 
  #  To Peak2:   index of the last peak to be used 

  # get means & variances

  MESFmeans <- MESFmeanssd2[Peak1:Peak2, 1]
  MESFvariances <- MESFmeanssd2[Peak1:Peak2, 2]

  # linear regression between means and variances
  QB <- lm(MESFvariances ~ MESFmeans)
  if(QB$coeff[2] != 0)
  {
    B=QB$coeff[1]/QB$coeff[2]
    Q=1/QB$coeff[2]
    # Note  QB@coeff[2] =1/Q and QB@coeff[1]= B*QB@coeff[2]
    return(c(Q=Q, B=B, Rsquared=round(summary(QB)$r.squared,2)))
  }
  if(QB$coeff[2] == 0) 
  {
    stop("No linear term to calculate Q and B")
  } 

}


# Quadratic Rgeression
qrMESF <- function(MESFmeanssd2, Peak1, Peak2)
{  

  #  MESFmeanssd2: 2D array having the means and thier variance.
  #  From Peak1: index of the first peak to be used. 
  #  To Peak2:   index of the last peak  to be used. 

  # get means & variances

  MESFmeans <- MESFmeanssd2[Peak1:Peak2, 1]
  MESFvariances <- MESFmeanssd2[Peak1:Peak2, 2]

  # quadratic regression between means and variances

  QB <- lm(MESFvariances ~ MESFmeans + I(MESFmeans^2))

  if(QB$coeff[2] != 0)
  {
    B <- QB$coeff[1]/QB$coeff[2]
    Q <- 1/QB$coeff[2]

    # Note  QB@coeff[2] =1/Q and QB@coeff[1]= B*QB@coeff[2]
    # We need to return only QB$coeff[3]

    return(c(Q=Q, B=B, Rsquared=round(summary(QB)$r.squared, 2), sigmaS2=QB$coeff[3]))
  }
  if(QB$coeff[2] == 0) 
  {
    stop("No linear term to calculate Q and B")
  }

}

# Discriminant examination

DiscriminantExamination <- function(Q, B, sigmaS2)
{  

  #  Q: the first value in the output of the function qrMESF .
  #  B: the second  value in the output of the function qrMESF. 
  #  sigmaS2: the fourth  value in the output of the function qrMESF. 



	if(Q == 0) 
	{
    	stop("Cytometry is not efficient because Q is equal to zero.")
   	} 

	if(B < 0) B <- 0

	if(Q != 0)
	{
	# The values of c0(sigmaE2) , c1(1/Q) and c2(sigmaS2) are:
 	c0 <- B/Q
 		c1 <- 1/Q
 	  		c2 <- sigmaS2
				DELTA <- c1^2 - 4*c0*c2

if(1==0)
{
if(DELTA >= 0) 
{
cat(paste("The sign of the discriminant is positive with the value",round(DELTA,2),"\n"))
cat("The larger the variation product of (sigmaE2) and (sigmaS2) \n")
cat("the lower the upper bound on the detection efficiency Q.\n")
}

if(DELTA < 0) 
{
cat(paste("The sign of the discriminant is negative with the value",round(DELTA,2),"\n"))
cat("The lower the variation product of (sigmaE2) and (sigmaS2) \n")
cat("the greater the upper bound on the detection efficiency Q. \n")
}
}

# Return the values of c0(sigmaE2) , c1(1/Q) and c2(sigmaS2) and the discriminant value

Coeffs <- c(c0, c1, c2, DELTA)
names(Coeffs) <- c("c0","c1","c2","Delta")

return(Coeffs)
}

}




