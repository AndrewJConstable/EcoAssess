# HIMI Temp BRAN2020 Temp decadal mean (top-mid) - summer & winter

# Generate raster brick containing rasters of the mean difference between top and mid depth layers
# for the summers and winters in the decades 1993-2003 and 2010-2020, along with rasters of the mean in the top layer

# Important notes for processing BRAN data
# 1. Using Terra package to import netCDF files.
# 2. Rasters are complete x-y matrices, including when subsetting
# 3. When converting rasters to dataframes using as.data.frame, missing values are dropped using defaults.
#          Important to not remove cells with NA by setting na.rm=F i.e. as.data.frame(r, na.rm=F)
#          and manage NAs in code.

#########################################################################################################
# 1. Libraries and Input Parameters
##################################

library(RNetCDF)
library(terra)

# 1.1 Directories
 BRANrootRasterReadStatic<-"/vsicurl/https://dapds00.nci.org.au/thredds/fileServer/gb6/BRAN/BRAN2020/static/"
 BRANrootRasterReadMonthData<-"/vsicurl/https://dapds00.nci.org.au/thredds/fileServer/gb6/BRAN/BRAN2020/month/ocean_temp_mth_"
 BRANrootNetCDFreadMonthData<-"https://dapds00.nci.org.au/thredds/dodsC/gb6/BRAN/BRAN2020/month/ocean_temp_mth_"

# 1.2 BRAN2022 Parameters
NAtemp        <- -32768  # missing values for temperature
Date.Origin   <- "1979-01-01"
Years         <-  seq(1993,2020,1)
Months        <-  sprintf("%02d",seq(1,12,1))

DecadeEarly<-c(1993:2002) # sequence winter to summer with last summer crossing into 2003
DecadeRecent<-c(2010:2019)# sequence winter to summer with last summer crossing into 2020

WinterMths<-c(4,5,6,7,8,9)  # April-September
SummerMths<-c(10,11,12,1,2,3) # October-March

# GIS
Arena<-c(60,90,-60,-45) # West, East, South, North boundaries

# Depth

dStrata<-list( epipelagic    = c(0,100)     # Surface conditions (approximately in mixed layer) - surface 100m
               ,mesopelagic   = c(300,800)  # Mesopelagic conditions - (300-800m)
) # end list

dBins_BRAN<-c(   2.5,  7.5, 12.5, 17.51539, 22.66702
                 , 28.16938, 34.21801, 40.95498, 48.45498, 56.71801
                 , 65.66938, 75.16702, 85.01539, 95.0, 105.0
                 , 115.0, 125.00000, 135.0, 145.0, 155.0
                 , 165.0, 175.0, 185.0, 195.0, 205.1899
                 , 217.05449, 233.19432, 255.88423, 286.60898, 325.88422
                 , 373.19434, 427.05447, 485.18991, 545.51111, 610.41565
                 , 685.92676, 775.92676, 880.41565, 995.51111, 1115.31335
                 , 1238.35388, 1368.15747, 1507.73389, 1658.15747, 1818.35388
                 , 1985.31335, 2165.18018, 2431.10107, 2894.84180, 3603.10107
                 , 4509.18018)

# Output
zGrpTempPC<-c(0.2,0.5,0.8) # spatial percentiles of Temperature within zone groups
# file for saving output
outName<-"/perm_storage/home/acon/MEASO/HIMI_Temp_10-S-W_rbrick "
#########################################################################################################
# 2. Functions
##################################

#########################################################################################################
# 4. Cell attributes
##################################

# x_T: geographic longitude, units - degrees E
# y_T: geographic latitude, units - degrees N
# area_T: Area of T_cell, units: m2

dfile<-paste(BRANrootRasterReadStatic,"grid_spec.nc",sep="")
r<-rast(dfile,"area_T")
rCellArea<-crop(r, ext(Arena))
CellNumbers<-row.names(as.data.frame(rCellArea,na.rm=F))  # for assigning to matrices then used in apply fns

# depth
dfile<-paste(BRANrootRasterReadStatic,"ocean_grid.nc",sep="")
# ht - ocean depth on t-cells
r<-rast(dfile,"ht")
rCellDepth<-crop(r, ext(Arena))

# kmt - number of depth levels on t-grid
r<-rast(dfile,"kmt")
rCellDepthLevels<-crop(r, ext(Arena))

dfCellAttr<-cbind(as.data.frame(rCellArea,na.rm=F),as.data.frame(rCellDepth,na.rm=F),as.data.frame(rCellDepthLevels,na.rm=F))


#########################################################################################################
# 6. Depth Strata for averaging Temperature
##################################

NdBins<-length(dBins_BRAN)
depthMask<-lapply(dStrata,function(d){
  res<-dBins_BRAN>=d[1] & dBins_BRAN<d[2]
  if(sum(is.na(res))==NdBins) return(!is.na(res)) else return(res)
})

# determine depth intervals of each depth bin for weighting calculations of mean etc.
# i. determine intervals between elements in dBins_BRAN
dBins_diff<-c(dBins_BRAN[1],dBins_BRAN[2:NdBins]-dBins_BRAN[1:(NdBins-1)])

# ii. loop through the vector of differences to calculate the interval    
prev<-dBins_diff[1]
depthBins_int<-prev*2
for (d in seq(2,NdBins,1)){
  prev<-dBins_diff[d]-prev
  depthBins_int<-c(depthBins_int,prev*2)
}

rm(dBins_diff,prev)

#########################################################################################################
# 7. File Names for each Decade x Season period
##################################

# Files for generating rasters - decade (early, recent) x season (winter, summer)

dEarly_WinterFiles<-cbind(as.character(unlist(lapply(DecadeEarly,function(y,n){rep(y,n)},length(WinterMths)))),Months[rep(WinterMths,length(DecadeEarly))])
colnames(dEarly_WinterFiles)<-c("Year","Month")
dEarly_SummerFiles<-cbind(as.character(unlist(lapply(DecadeEarly,function(y,m){yrs<-rep(y,length(m)); yrs[m<10]<-y+1; return(yrs)},SummerMths))),Months[rep(SummerMths,length(DecadeEarly))])
colnames(dEarly_SummerFiles)<-c("Year","Month")

dRecent_WinterFiles<-cbind(as.character(unlist(lapply(DecadeRecent,function(y,n){rep(y,n)},length(WinterMths)))),Months[rep(WinterMths,length(DecadeRecent))])
colnames(dRecent_WinterFiles)<-c("Year","Month")
dRecent_SummerFiles<-cbind(as.character(unlist(lapply(DecadeRecent,function(y,m){yrs<-rep(y,length(m)); yrs[m<10]<-y+1; return(yrs)},SummerMths))),Months[rep(SummerMths,length(DecadeRecent))])
colnames(dRecent_SummerFiles)<-c("Year","Month")

perFiles<-list(
  eWinter = apply(dEarly_WinterFiles,1,function(d) paste(d[1],"_",d[2],sep=""))
 ,eSummer = apply(dEarly_SummerFiles,1,function(d) paste(d[1],"_",d[2],sep=""))
 ,rWinter = apply(dRecent_WinterFiles,1,function(d) paste(d[1],"_",d[2],sep=""))
 ,rSummer = apply(dRecent_SummerFiles,1,function(d) paste(d[1],"_",d[2],sep=""))
)

PeriodNames<-c("eWinter","eSummer","rWinter","rSummer")

#########################################################################################################
# 9. Analysis - for each decade x season period, 
#                generate rasters for top depth layer, mid depth layer and top-depth of 
#                the mean over the season then mean over years

##################################
# develop raster template
FilePointerRaster<-paste(BRANrootRasterReadMonthData,perFiles[[1]][1],".nc",sep="")
r           <- crop(rast(FilePointerRaster), ext(Arena))
rTemplate<-subset(r, 1)
names(rTemplate)<-"Temp"

# by period
resPeriod<-lapply(seq(1,length(PeriodNames)),function(p,Name,f){
  Files<-f[[p]]
  # by year
     dYearTop<-NULL; dYearMid<-NULL; dYearDiff<-NULL
    for(y in c(1:10)){  
      dSeasonTop<-NULL; dSeasonMid<-NULL; dSeasonDiff<-NULL
       for(m in c(1:6)){ 
         File_n<-(y-1)*6+m
         print(Files[File_n])
         
         # read file & make into data frame
         FilePointerRaster<-paste(BRANrootRasterReadMonthData,Files[File_n],".nc",sep="")
                   # delete? FilePointerNetCDF<-paste(BRANrootNetCDFreadMonthData,FileDate,".nc",sep="")
         
         r           <- rast(FilePointerRaster)
         rMonthTemp  <- crop(r, ext(Arena))
         
         dMTdf       <-as.data.frame(rMonthTemp,na.rm=F)  # note that missing values are dropped with default na.rm=T; come out as NA in vector
         dMT<-as.matrix(dMTdf)
         dMT[dMT==NAtemp]<-NA # make missing values = NA
         dimnames(dMT)<-NULL
         dimnames(dMT)[[1]]<-row.names(dMTdf) 
         rm(dMTdf)
         
         # mean temp in each cell in the top layer
         d<-depthMask[[1]]      # depth mask for top layer
         res_top<- apply(dMT,1,function(dT,d,DI){
                    dTsub<-dT[d]
                    DIsub<-DI[d]
                    if(sum(!is.na(dTsub))>0){
                          DIsub<-DIsub[!is.na(dTsub)]
                          dTsub<-dTsub[!is.na(dTsub)]
                         return(sum(dTsub*DIsub)/sum(DIsub))
                         }  else {return(NA)}}
                        ,d,depthBins_int)

         # mean temp in each cell in the mid layer
         d<-depthMask[[2]]      # depth mask for top layer
         res_mid<- apply(dMT,1,function(dT,d,DI){
                   dTsub<-dT[d]
                   DIsub<-DI[d]
                   if(sum(!is.na(dTsub))>0){
                           DIsub<-DIsub[!is.na(dTsub)]
                          dTsub<-dTsub[!is.na(dTsub)]
                          return(sum(dTsub*DIsub)/sum(DIsub))
                      }  else {return(NA)}}
                   ,d,depthBins_int)
         
         res_diff<-res_top-res_mid
         
         dSeasonTop<-cbind(dSeasonTop,res_top)
         dSeasonMid<-cbind(dSeasonMid,res_mid)
         dSeasonDiff<-cbind(dSeasonDiff,res_diff)
    } # end m
      
      #       after processing each file, calculate three mean values for each cell for the season to be saved in respective matrices
      #                          i) mean temperature in the top layer
      #                         ii) mean temperature in the middle layer
      #                        iii) for each month calculate the difference between top and middle layers and then take the mean difference for the season
      dYearTop<-cbind(dYearTop,apply(dSeasonTop,1,mean))
      names
      dYearMid<-cbind(dYearMid,apply(dSeasonMid,1,mean))
      dYearDiff<-cbind(dYearDiff,apply(dSeasonDiff,1,mean))
      } # end y
  # after processing all years, take thee mean of each of the three metrics in each cell and
  # add them to raster template as values
     rTopMean<-rTemplate; rMidMean<-rTemplate; rDiffMean<-rTemplate
     values(rTopMean)<-apply(dYearTop,1,mean)
     names(rTopMean)<-paste(Name[p],"_Top",sep="")
     values(rMidMean)<-apply(dYearMid,1,mean)
     names(rMidMean)<-paste(Name[p],"_Mid",sep="")
     values(rDiffMean)<-apply(dYearDiff,1,mean)
     names(rDiffMean)<-paste(Name[p],"_Diff",sep="")
     return(c(rTopMean, rMidMean, rDiffMean))
########  
    },PeriodNames,perFiles) # end lapply periods
resPeriod<-do.call("c",resPeriod)
print(resPeriod)

writeRaster(resPeriod, filename=paste(outName,Sys.time(),".tiff",sep=""))

print("File Saved and Finished")
