library(raadtools)
library(raadfiles)
library(lubridate)
library(terra)
# read files for PAR and convert to dates

Date1<-"2002-07-04"
Years<-2

PARdates<-as.Date(as.vector(par_files()[,"date"])[[1]])
DaysOfYear<-sort(unique(yday(PARdates)))
nPeriods<-length(DaysOfYear)

Start<-which(PARdates %in% Date1)
End<-Years*nPeriods-1+Start
subPAR<-PARdates[c(Start:End)]
Ydays<-yday(subPAR)
rPAR<-rast(sapply(DaysOfYear,function(d,Ydays,PAR){
               sPAR<-PAR[Ydays==d]
               mean(rast(read_par(sPAR,xylim=c(-180,180,-90,-30))),na.rm=TRUE)
                },Ydays,subPAR))
writeRaster(rPAR,filename="rasterbrickPAR.tiff")

