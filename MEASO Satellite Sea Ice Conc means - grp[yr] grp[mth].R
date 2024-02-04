# MEASO Satellite Sea ice means - grp[yr] grp[mth].R

# generate raster of mean sea ice concentration
# in following sequence
#  1. mean sea ice concentration across months as designated in input vector
#. 2. mean sea ice concentration across years as designated in input vector
# 

library(raadtools)
library(raadfiles)
library(terra)

# 1. Input Parameters

DaysStart<-as.Date("2011-07-01")
DaysEnd  <-as.Date("2021-06-30")   
dates<-seq(DaysStart,DaysEnd,by="1 month")

MeanForMths<-c(1:12)
refMonths<-sprintf("%02d",MeanForMths)

Arena<-c(-180,180,-90,-40) # West, East, South, North boundaries

rCICE<-rast(lapply(refMonths,function(m,dates){
  print(paste0("Start Time - Month ",m,": ",Sys.time()))
  r<-rast(readice_monthly(dates[substr(dates,6,7)==m],setNA=FALSE))
  set.crs(r, "EPSG:3031")
  values(r)[values(r)>100]<--1000 # need to remove continents
  rMn<-mean(r,na.rm=TRUE)
  names(rMn)<-m
  return(rMn)
    },dates))

targetgrid <- rast(ext(Arena), res = c(0.25, 0.25), crs = "EPSG:4326")
rCICE2<-project(rCICE, targetgrid) #, by_util = TRUE
values(rCICE2)[is.nan(values(rCICE2))]<-0  # reprojection leaves NaN in northern areas - convert to 0 sea ice
values(rCICE2)[values(rCICE2)<0]<-NA
writeRaster(rCICE2,file="rCICE.tif",overwrite=TRUE)

tmp<-rast("rCICE.tif")