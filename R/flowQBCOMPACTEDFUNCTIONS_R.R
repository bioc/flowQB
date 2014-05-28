BEADflowQBCalculation <- function(file,R,width,fraction,numberoftotalepeaks,scatters,CHANNELTOBEPROCESSED,pi,pf,Viz,output)
{


# file: path plus the name of the FCS file if R is running from a folder different from where the file is stored.
# R: gate by an ellipse of radius R on the scatter data.
# width: a parameter used to find the lower bound (default width at half maximum) and the upper bound (default .5).
# fraction: a parameter used to find the densest fraction of cells (default .1).
# numberoftotalepeaks: Global Number of the peaks in the data.
# scatters: the indices of Forward Scatter (FS) and Side scatter (SS) in the FCS file.
# CHANNELTOBEPROCESSED : the indices of channels of interest in the FCS file.
# Peaks to be used in the regression,  ( pi < pf ) where
# pi: First DIMMEST peak to be used in the regression
# pf: Last peak to be used in the regression
# output: where the Kmeans image results should be stored, see Viz for more details.
# Viz: Parameter to control the visualization of Kmeans results only if Viz=1, a subfolder called flowQB in the output folder has subfolder "PLOTSKMEANS_BEADS_" with the visualizations.

# TO READ FCS FILE
library(flowCore)
# TO EXTRACT GAUSS WITHOUT OUTLIERS
# library(extremevalues)
# TO EXTRACT GAUSS WITH OUTLIERS
library("MASS")
# THE COLORS USED FOR THE PEAKS
colors=c("skyblue","grey","blue","green","pink","yellow","red","gold","violet")

# CREATE THE FIRST SUB-FOLDER WHERE THE RESULTS WILL BE STORED 
# NAME THIS DATA WITH ITS SPECIFIC NAME GIVEN BY THE LAB-INSTRUMENT
folderx=sub(".fcs","_fcs",file)
if( Viz ==1)
{
# CREATE flowQB FOLDER FOR THIS INPUT 
setwd(output)
print(output)
output0=paste(output,"/flowQB/",sep="")
dir.create(output0)
setwd(output0)
print(output0)
# CREATE THE SUBFOLDERS WHERE Kmeans Image RESULTS WILL BE STORED 
output4=paste(output0,"/","PLOTSKMEANS_BEADS_",numberoftotalepeaks,"X",sep="")
dir.create(output4)
print(output4)
}

# THIS LIST STORES THE RAW STATISTICS AND THE COEFFICIENTS, Q and B values with their "Pvalue" & "Std-Error" 
RESULTSLIST=list()

# NOW GO THE THE INPUT 
fcs <- read.FCS(file)


        for (scatter in scatters)
        {
            data <- exprs(fcs[,scatter])
            fcs <- fcs[data>0]
        }
        scatter.gated <- fcs[scatter_gate(fcs, scatters, R)]
 

# Express the data
Data=scatter.gated
orgdata=data.frame(exprs(Data))

# PRINT SINGLETS PROPORTION:
#print(paste("Singlet Proportion: ", 100*round(as.numeric(dim(Data)[1])/as.numeric(dim(fcs)[1]),2),"%"))

# GET THE DESCRIPTION
fcs <- read.FCS(file)
# PRINT SINGLETS PROPORTION:
DD=description(fcs)
# A PRIORI WE DONT KNOW THE NUMBER OF THE CHANNELS SET THE NUMBER AS 100.
# THIS PROCESS WILL GET THE NAME OF THE CHANNELS.
markers=1:100
k=0
for( i in 1:100)
{
PN=paste("P",i,"N",sep="")
MARK=grep(PN,names(DD))
if( length(MARK)>0 ) 
{
k=k+1
markers[k]=DD[MARK]
}
}
# GET ONLY ALL CHANNELS.
CHANNELSx=unlist(markers[1:k])

###################################################
### PREPARATION OF THE DATA WITHOUT DOUBLETS
###################################################
# Process only CHANNELS OF INTERST.
CHANNELS=CHANNELSx[CHANNELTOBEPROCESSED]
DATASIN=scatter.gated[,CHANNELTOBEPROCESSED]
nm=length(CHANNELTOBEPROCESSED)
# TURN MARKERSdata TO A DATAFRAME
# GET THE ORIGINAL MFI 
OriginalMFI=data.frame(exprs(DATASIN))
minx=apply(OriginalMFI,2,min)
# print(minx)
newda=OriginalMFI
for( jxi in 1:dim(OriginalMFI)[2])
{
newda[,jxi]=OriginalMFI[,jxi]-minx[jxi]
}
minx=apply(newda,2,min)
#print(minx)
MARKERSdata=data.frame(log(1+newda))
names(MARKERSdata)=CHANNELS

# CHECK THE names(OriginalMFI) IF ANY DOUBT 
# CLUSTERING THE TRANSFORMED DATA
# DEFINE THE NUMBER OF THE PEAKS OR CLUSTERS
# STEP KMEANS MULTIPROCESSING
x <- MARKERSdata
# PLEASE USE HELP TO CHECK ON kmeans METHOD AND THE NEEDED PARAMETRIZATION
nClusters=numberoftotalepeaks
cl <- kmeans(x, nClusters,500,500)


# STEP STATISTICS
Markers=names(MARKERSdata)
MK=0
for( j in 1:dim(MARKERSdata)[2])
{
MK=MK+1
# ROBUST STATISTICS 
MEANSr=rep(0,nClusters)
SDsr=rep(0,nClusters)
Nesr=rep(0,nClusters)
# GAUSS STATISTICS USING MASS LIBRARY
MEANSg=rep(0,nClusters)
SDsg=rep(0,nClusters)
Nesg=rep(0,nClusters)
# GAUSS STATISTICS USING MASS EXTREMES ( REMOVING OUTLIERS)
MEANSno=rep(0,nClusters)
SDsno=rep(0,nClusters)
Nesno=rep(0,nClusters)
# USING KMEANS RESULTS TO EXTRACT STATISTICS FOR EACH PEAK 
for(jc in 1:nClusters) 
{
eventsjc=which(cl$cluster==jc)
vect=OriginalMFI[eventsjc,j]

# ROBUST
MEANSr[jc]=median(vect)
SDsr[jc]=median(abs(vect-median(vect)))*1.4826
if(SDsr[jc]==0) SDsr[jc]=.1
Nesr[jc]=length(eventsjc)

# GAUSS MASS
fitgmass <- unlist(fitdistr(vect, 'normal'))
MEANSg[jc]=fitgmass[1]
SDsg[jc]=fitgmass[2]
if(SDsg[jc]==0) SDsg[jc]=.1
Nesg[jc]=length(eventsjc)

### remove outliers
p <- c(2,98)/100
y <- quantile(vect,p)
QQx=which( vect < y[2] & vect> y[1])
vect=vect[QQx]
fitgmass <- unlist(fitdistr(vect, 'normal'))
# GAUSS EXTREMES
MEANSno[jc]=fitgmass[1]
SDsno[jc]= fitgmass[2]
if(SDsno[jc]==0) SDsno[jc]=.1
Nesno[jc]=length(vect)

}

# SORT THE PEAKS USING ROBUST STATISTICS
ORD=order(MEANSr)
# STORE THE RESULTS OF KMEANS AS A PLOT FOR EACH CHANNEL OR MARKER
if( Viz ==1)
{
setwd(output4)
naf0=sub(" " ,"", Markers[j])
naf=paste(naf0,"-MultiBEAD8X.jpeg",sep="")

# TO HAVE A 2D PLOT ONE CAN USE A COMPANION CHANNEL
j1=1
# AS JPEG
jpeg(naf)
plot(MARKERSdata[,j],MARKERSdata[,j1],xlab=names(MARKERSdata)[j],ylab=names(MARKERSdata)[j1],main="K-means")
gg=0
for( r in ORD)
{
eventsjc=which(cl$cluster==r)
vect=OriginalMFI[eventsjc,j]
vectx=x[eventsjc,c(j,j1)]
gg=gg+1
points(vectx[,1],vectx[,2], col = colors[gg])
}
dev.off()
}
###################################################
### STATISTICS ROBUST-GAUSS FROM MASS AND GAUSS FROM EXTREME ARE STORED IN A FILE.
###################################################
datarg=data.frame(NE=Nesr[ORD],mfiRS=MEANSr[ORD],mfiGS=MEANSg[ORD],mfino=MEANSno[ORD],mfirSD=SDsr[ORD],mfigSD=SDsg[ORD],mfiSDno=SDsno[ORD],Nesno[ORD])
naf0x=sub(" " ,"", Markers[j])
row.names(datarg)=paste(naf0x,"-Peak",1:nClusters,sep="")

RESULTSLIST[[MK]]=datarg


#setwd(output2)
#write.csv(datarg,paste(Markers[j],"MultiCPeakMeanandSDvaluesBEAD8X.csv",sep=""))

# THE WEIGHTS IN THE WEIGHTED QUADRATIC REGRESSION DEPEND ON THE NUMBER OF EVENTS IN EACH PEAK. 
NUGAUSS=Nesg[ORD]
NUno=Nesno[ORD]
NUROBUST=Nesr[ORD]

# STEP WEIGHTED QUADRATIC REGRESSION
###################################################
### ROBUST
###################################################
rPSdataB515=datarg[,c(1,2,5)]
CVs=rPSdataB515[,3]/rPSdataB515[,2]
###################################################
### DATA AS INPUT FOR THE REGRESSION
###################################################
Q2D=data.frame(rPSdataB515[,2],rPSdataB515[,3]^2)
 # Q2D[1,1]=140
###################################################
# NUMBER OF EVENTS
NE=NUROBUST[pi:pf]
# MFI 
mfino=Q2D[pi:pf,1]
# MFI VARIATION
mfiVarno=Q2D[pi:pf,2]

# GUESSING THE WEIGHTS
WWs=0.5*(NE-1)/(mfiVarno^2)
# THE NUMBER OF ITERATION IS SET TO 15
#WWs=1/WWs
hlf1=round(length(pi:pf)/2,0)+1
hlf=round(length(pi:pf)/2,0)
#WWs[pf:hlf1]=WWs[pi:hlf]
for( iter in 1:15)
{
QQBw=lm(formula = mfiVarno ~ mfino + I(mfino^2),weights = WWs)
Pvalues=summary(QQBw)$coef[,4]
StdError=summary(QQBw)$coef[,2]
coefsx=summary(QQBw)$coef[,1]
c0=as.numeric(coefsx)[1]
c1=as.numeric(coefsx)[2]
c2=as.numeric(coefsx)[3]
COFS=c(c0,as.numeric(Pvalues)[1],as.numeric(StdError)[1],c1,as.numeric(Pvalues)[2],as.numeric(StdError)[2],c2,as.numeric(Pvalues)[3],as.numeric(StdError)[3])
if(iter==2) COEFFS=COFS
if(iter>2) COEFFS=rbind(COEFFS,COFS)
WWs=0.5*(NE-1)/((c2*mfino^2+c1*mfino+c0)^2)
#WWs=1/WWs
hlf1=round(length(pi:pf)/2,0)+1
hlf=round(length(pi:pf)/2,0)
#WWs[pf:hlf1]=WWs[pi:hlf]
}

c0=as.numeric(c0)
c1=as.numeric(c1) 
c2=as.numeric(c2)
B=c0/c1
Q=1/c1
COFSrx=data.frame(t(COFS))
names(COFSrx)=c("c0","Pvalue","Std-Error","c1","Pvalue","Std-Error","c2","Pvalue","Std-Error")
coefsr=c("ROBUSTStatsBeads",Markers[j],COFSrx,B,Q)
if(MK==1) MARKERSRESULTS=unlist(coefsr)
if(MK>1) MARKERSRESULTS=rbind(MARKERSRESULTS,unlist(coefsr))


###################################################
### GAUSS-MASS
###################################################
rPSdataB515=datarg[,c(1,3,6)]
CVs=rPSdataB515[,3]/rPSdataB515[,2]
###################################################
### DATA AS INPUT FOR THE REGRESSION
###################################################
Q2D=data.frame(rPSdataB515[,2],rPSdataB515[,3]^2)
###################################################
### GET THE NUMBER OF EVENTS
###################################################
NE=NUGAUSS[pi:pf]
###################################################
# MFI 
mfino=Q2D[pi:pf,1]
# MFI VARIATION
mfiVarno=Q2D[pi:pf,2]

# GUESSING THE WEIGHTS
WWs=0.5*(NE-1)/(mfiVarno^2)
# THE NUMBER OF ITERATION IS SET TO 15

for( iter in 1:15)
{
QQBw=lm(formula = mfiVarno ~ mfino + I(mfino^2),weights = WWs)
Pvalues=summary(QQBw)$coef[,4]
StdError=summary(QQBw)$coef[,2]
coefsx=summary(QQBw)$coef[,1]
c0=as.numeric(coefsx)[1]
c1=as.numeric(coefsx)[2]
c2=as.numeric(coefsx)[3]
COFS=c(c0,as.numeric(Pvalues)[1],as.numeric(StdError)[1],c1,as.numeric(Pvalues)[2],as.numeric(StdError)[2],c2,as.numeric(Pvalues)[3],as.numeric(StdError)[3])
if(iter==2) COEFFS=COFS
if(iter>2) COEFFS=rbind(COEFFS,COFS)
WWs=0.5*(NE-1)/((c2*mfino^2+c1*mfino+c0)^2)
}

c0=as.numeric(c0)
c1=as.numeric(c1) 
c2=as.numeric(c2)
B=c0/c1
Q=1/c1
COFSrx=data.frame(t(COFS))
names(COFSrx)=c("c0","Pvalue","Std-Error","c1","Pvalue","Std-Error","c2","Pvalue","Std-Error")
coefsg=c("GAUSSMASSStatsBeads",Markers[j],COFSrx,B,Q)
MARKERSRESULTS=rbind(MARKERSRESULTS,unlist(coefsg))


###################################################
### GAUSS-EXTREMES
###################################################
rPSdataB515=datarg[,c(8,4,7)]
CVs=rPSdataB515[,3]/rPSdataB515[,2]
###################################################
### DATA AS INPUT FOR THE REGRESSION
###################################################
Q2D=data.frame(rPSdataB515[,2],rPSdataB515[,3]^2)
###################################################
### GET THE NUMBER OF EVENTS
###################################################
NE=NUno[pi:pf]
###################################################
# MFI 
mfino=Q2D[pi:pf,1]
# MFI VARIATION
mfiVarno=Q2D[pi:pf,2]

# GUESSING THE WEIGHTS
WWs=0.5*(NE-1)/(mfiVarno^2)
# THE NUMBER OF ITERATION IS SET TO 15

for( iter in 1:15)
{
QQBw=lm(formula = mfiVarno ~ mfino + I(mfino^2),weights = WWs)
Pvalues=summary(QQBw)$coef[,4]
StdError=summary(QQBw)$coef[,2]
coefsx=summary(QQBw)$coef[,1]
c0=as.numeric(coefsx)[1]
c1=as.numeric(coefsx)[2]
c2=as.numeric(coefsx)[3]
COFS=c(c0,as.numeric(Pvalues)[1],as.numeric(StdError)[1],c1,as.numeric(Pvalues)[2],as.numeric(StdError)[2],c2,as.numeric(Pvalues)[3],as.numeric(StdError)[3])
if(iter==2) COEFFS=COFS
if(iter>2) COEFFS=rbind(COEFFS,COFS)
WWs=0.5*(NE-1)/((c2*mfino^2+c1*mfino+c0)^2)
}

c0=as.numeric(c0)
c1=as.numeric(c1) 
c2=as.numeric(c2)
B=c0/c1
Q=1/c1
COFSrx=data.frame(t(COFS))
names(COFSrx)=c("c0","Pvalue","Std-Error","c1","Pvalue","Std-Error","c2","Pvalue","Std-Error")
coefsg=c("GAUSSEXTREMESStatsBeads",Markers[j],COFSrx,B,Q)
MARKERSRESULTS=rbind(MARKERSRESULTS,unlist(coefsg))
# print(paste(Markers[j],"DONE"))

}


results=data.frame(MARKERSRESULTS)
names(results)=c("StatsProcedure","MARKER",names(COFSrx),"B","Q")
row.names(results)=1:dim(results)[1]

RESULTSLISTx=list()
RESULTSLISTx[[1]]=RESULTSLIST
RESULTSLISTx[[2]]=results
# print(paste(folderx,"DONE",sep=" "))
return(RESULTSLISTx)



}




find_peak <- function( data, width=.5, fraction=.1 )
{
    x <- sort(data)
    N <- length(data)
    M <- ceiling(N*fraction)
    M2 <- floor((M+1)/2)
    for (first in 1:(N-1)) if (x[first] < x[first+1]) break;
    for (last in N:2) if (x[last-1] < x[last]) break;
    i <- first
    dx <- x[last] - x[first]
    lo <- first
    hi <- last - M

    #find the densest fraction (default 10%) of cells
    for (j in first:(last-M+1))
    {
        d <- x[j+M-1]-x[j];d
        if (d < dx)
        {
            i <- j
            dx <- d
        }
    }
    #find the lower bound (default width at half maximum)
    for (j in i:first)
    {
        d <- x[j+M-1]-x[j]
        if (d > dx/width)
        {
            lo <- j
            break
        }
    }
    #find the upper bound (default width at half maximum)
    for (j in i:(last-M+1))
    {
        d <- x[j+M-1]-x[j]
        if (d > dx/width)
        {
            hi <- j
            break
        }
    }
    list( lo=x[lo+M2], hi=x[hi+M2] )
}

scatter_gate <- function( fcs, scatters, R=1 )
{
    #gate by an ellipse of radius R on the scatter data
    R2 <- vector( mode="double", nrow(fcs) )
    for (scatter in scatters)
    {
        data <- log(exprs(fcs[,scatter]))
        fwhm <- find_peak(data)
        center <- (fwhm$hi + fwhm$lo)/2
        radius <- (fwhm$hi - fwhm$lo)/2
        R2 = R2 + ((data - center)/radius)^2
    }
  # print(R)
    R2 <= R^2
}




LEDflowQBCalculation <- function(input,CHANNELTOBEPROCESSED,minF,maxF)
{
# input: Path with the names of the files.
# CHANNELTOBEPROCESSED : the names of channels of interest in the FCS file.
# minF: Minimum MFI to consider  of the events.
# maxF: Maximum MFI to consider  of the events.

# TO READ FCS FILE
library(flowCore)
# TO EXTRACT GAUSS WITHOUT OUTLIERS
# library(extremevalues)
# TO EXTRACT GAUSS WITH OUTLIERS
library("MASS")

folderx="LED"
# THIS LIST STORES THE RAW STATISTICS AND THE COEFFICIENTS, Q and B values with their "Pvalue" & "Std-Error" 
RESULTSLIST=list()
ALLFILES=input
# GET THE LED OF INTEREST
 
files=ALLFILES[grep(".fcs",ALLFILES)]
if( length(files) ==0) print("Sorry, no LED file in this LED folder")

if( length(files)>0)
{
nClusters=length(files)
MK=0
for( CHANNEL in CHANNELTOBEPROCESSED)
{
MK=MK+1
# CHANNEL="FITC.A"
MEANSr=rep(0,nClusters)
SDsr=rep(0,nClusters)
Nesr=rep(0,nClusters)

MEANSg=rep(0,nClusters)
SDsg=rep(0,nClusters)
Nesg=rep(0,nClusters)

MEANSno=rep(0,nClusters)
SDsno=rep(0,nClusters)
Nesno=rep(0,nClusters)

jc=0
for( file in files)
{
jc=jc+1
Data=read.FCS(file)
orgdata=data.frame(exprs(Data))
i=which(CHANNEL==names(orgdata))
if( length(i)>0)
{
vect=orgdata[,i]
# ROBUST
MEANSr[jc]=median(vect)
SDsr[jc]=median(abs(vect-median(vect)))*1.4826
if(SDsr[jc]==0) SDsr[jc]=.1
Nesr[jc]=length(vect)

# GAUSS MASS
fitgmass <- unlist(fitdistr(vect, 'normal'))
MEANSg[jc]=fitgmass[1]
SDsg[jc]=fitgmass[2]
if(SDsg[jc]==0) SDsg[jc]=.1
Nesg[jc]=length(vect)

### remove outliers
p <- c(2,98)/100
y <- quantile(vect,p)
QQx=which( vect < y[2] & vect> y[1])
vect=vect[QQx]
fitgmass <- unlist(fitdistr(vect, 'normal'))
# GAUSS EXTREMES
MEANSno[jc]=fitgmass[1]
SDsno[jc]= fitgmass[2]
if(SDsno[jc]==0) SDsno[jc]=.1
Nesno[jc]=length(vect)

}
}

# SORT THE PEAKS USING ROBUST STATISTICS
ORDx=order(MEANSr)
ORD=ORDx[which( MEANSr[ORDx] >= minF & MEANSr[ORDx] <= maxF )]

# STORE THE RESULTS OF KMEANS AS A PLOT FOR EACH CHANNEL OR MARKER

###################################################
### STATISTICS ROBUST-GAUSS FROM MASS AND GAUSS FROM EXTREME ARE STORED IN A FILE.
###################################################
datarg=data.frame(NE=Nesr[ORD],mfiRS=MEANSr[ORD],mfiGS=MEANSg[ORD],mfino=MEANSno[ORD],mfirSD=SDsr[ORD],mfigSD=SDsg[ORD],mfiSDno=SDsno[ORD],Nesno[ORD])

row.names(datarg)=paste(CHANNEL,"-Peak",1:length(ORD),sep="")
RESULTSLIST[[MK]]=datarg


# THE WEIGHTS IN THE WEIGHTED QUADRATIC REGRESSION DEPEND ON THE NUMBER OF EVENTS IN EACH PEAK. 
NUGAUSS=Nesg[ORD]
NUno=Nesno[ORD]
NUROBUST=Nesr[ORD]

# STEP WEIGHTED QUADRATIC REGRESSION

###################################################
### PEAKS TO BE USED IN THE REGRESSION
###################################################
pi=1
pf=length(ORD)
###################################################
### ROBUST
###################################################
rPSdataB515=datarg[,c(1,2,5)]
CVs=rPSdataB515[,3]/rPSdataB515[,2]
###################################################
### DATA AS INPUT FOR THE REGRESSION
###################################################
Q2D=data.frame(rPSdataB515[,2],rPSdataB515[,3]^2)
###################################################
# NUMBER OF EVENTS
NE=NUROBUST[pi:pf]
# MFI 
mfino=Q2D[pi:pf,1]
# MFI VARIATION
mfiVarno=Q2D[pi:pf,2]

# GUESSING THE WEIGHTS
WWs=0.5*(NE-1)/(mfiVarno^2)
# THE NUMBER OF ITERATION IS SET TO 15

for( iter in 1:15)
{
QQBw=lm(formula = mfiVarno ~ mfino + I(mfino^2),weights = WWs)
Pvalues=summary(QQBw)$coef[,4]
StdError=summary(QQBw)$coef[,2]
coefsx=summary(QQBw)$coef[,1]
c0=as.numeric(coefsx)[1]
c1=as.numeric(coefsx)[2]
c2=as.numeric(coefsx)[3]
COFS=c(c0,as.numeric(Pvalues)[1],as.numeric(StdError)[1],c1,as.numeric(Pvalues)[2],as.numeric(StdError)[2],c2,as.numeric(Pvalues)[3],as.numeric(StdError)[3])
if(iter==2) COEFFS=COFS
if(iter>2) COEFFS=rbind(COEFFS,COFS)
WWs=0.5*(NE-1)/((c2*mfino^2+c1*mfino+c0)^2)
}

c0=as.numeric(c0)
c1=as.numeric(c1) 
c2=as.numeric(c2)
B=c0/c1
Q=1/c1
COFSrx=data.frame(t(COFS))
names(COFSrx)=c("c0","Pvalue","Std-Error","c1","Pvalue","Std-Error","c2","Pvalue","Std-Error")
coefsr=c("ROBUSTStatsBeads",CHANNEL,COFSrx,B,Q)
if(MK==1) MARKERSRESULTS=unlist(coefsr)
if(MK>1) MARKERSRESULTS=rbind(MARKERSRESULTS,unlist(coefsr))
###################################################
### GAUSS-MASS
###################################################
rPSdataB515=datarg[,c(1,3,6)]
CVs=rPSdataB515[,3]/rPSdataB515[,2]
###################################################
### DATA AS INPUT FOR THE REGRESSION
###################################################
Q2D=data.frame(rPSdataB515[,2],rPSdataB515[,3]^2)
###################################################
### GET THE NUMBER OF EVENTS
###################################################
NE=NUGAUSS[pi:pf]
###################################################
# MFI 
mfino=Q2D[pi:pf,1]
# MFI VARIATION
mfiVarno=Q2D[pi:pf,2]

# GUESSING THE WEIGHTS
WWs=0.5*(NE-1)/(mfiVarno^2)
# THE NUMBER OF ITERATION IS SET TO 15

for( iter in 1:15)
{
QQBw=lm(formula = mfiVarno ~ mfino + I(mfino^2),weights = WWs)
Pvalues=summary(QQBw)$coef[,4]
StdError=summary(QQBw)$coef[,2]
coefsx=summary(QQBw)$coef[,1]
c0=as.numeric(coefsx)[1]
c1=as.numeric(coefsx)[2]
c2=as.numeric(coefsx)[3]
COFS=c(c0,as.numeric(Pvalues)[1],as.numeric(StdError)[1],c1,as.numeric(Pvalues)[2],as.numeric(StdError)[2],c2,as.numeric(Pvalues)[3],as.numeric(StdError)[3])
if(iter==2) COEFFS=COFS
if(iter>2) COEFFS=rbind(COEFFS,COFS)
WWs=0.5*(NE-1)/((c2*mfino^2+c1*mfino+c0)^2)
}

c0=as.numeric(c0)
c1=as.numeric(c1) 
c2=as.numeric(c2)
B=c0/c1
Q=1/c1
COFSrx=data.frame(t(COFS))
names(COFSrx)=c("c0","Pvalue","Std-Error","c1","Pvalue","Std-Error","c2","Pvalue","Std-Error")
coefsg=c("GAUSSMASSStatsBeads",CHANNEL,COFSrx,B,Q)
MARKERSRESULTS=rbind(MARKERSRESULTS,unlist(coefsg))


###################################################
### GAUSS-EXTREMES
###################################################
rPSdataB515=datarg[,c(8,4,7)]
CVs=rPSdataB515[,3]/rPSdataB515[,2]
###################################################
### DATA AS INPUT FOR THE REGRESSION
###################################################
Q2D=data.frame(rPSdataB515[,2],rPSdataB515[,3]^2)
###################################################
### GET THE NUMBER OF EVENTS
###################################################
NE=NUno[pi:pf]
###################################################
# MFI 
mfino=Q2D[pi:pf,1]
# MFI VARIATION
mfiVarno=Q2D[pi:pf,2]

# GUESSING THE WEIGHTS
WWs=0.5*(NE-1)/(mfiVarno^2)
# THE NUMBER OF ITERATION IS SET TO 15

for( iter in 1:15)
{
QQBw=lm(formula = mfiVarno ~ mfino + I(mfino^2),weights = WWs)
Pvalues=summary(QQBw)$coef[,4]
StdError=summary(QQBw)$coef[,2]
coefsx=summary(QQBw)$coef[,1]
c0=as.numeric(coefsx)[1]
c1=as.numeric(coefsx)[2]
c2=as.numeric(coefsx)[3]
COFS=c(c0,as.numeric(Pvalues)[1],as.numeric(StdError)[1],c1,as.numeric(Pvalues)[2],as.numeric(StdError)[2],c2,as.numeric(Pvalues)[3],as.numeric(StdError)[3])
if(iter==2) COEFFS=COFS
if(iter>2) COEFFS=rbind(COEFFS,COFS)
WWs=0.5*(NE-1)/((c2*mfino^2+c1*mfino+c0)^2)
}

c0=as.numeric(c0)
c1=as.numeric(c1) 
c2=as.numeric(c2)
B=c0/c1
Q=1/c1

COFSrx=data.frame(t(COFS))
names(COFSrx)=c("c0","Pvalue","Std-Error","c1","Pvalue","Std-Error","c2","Pvalue","Std-Error")
coefsg=c("GAUSSEXTREMESStatsBeads",CHANNEL,COFSrx,B,Q)
MARKERSRESULTS=rbind(MARKERSRESULTS,unlist(coefsg))

# print(paste(CHANNEL,"DONE"))

}

results=data.frame(MARKERSRESULTS)
names(results)=c("StatsProcedure","MARKER",names(COFSrx),"B","Q")
row.names(results)=1:dim(results)[1]
RESULTSLISTx=list()
RESULTSLISTx[[1]]=RESULTSLIST
RESULTSLISTx[[2]]=results
# print(paste(folderx,"DONE",sep=" "))
return(RESULTSLISTx)
}

}






