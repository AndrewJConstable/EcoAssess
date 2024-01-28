# BRAN 2022 HIMI Temp Zgroup x Depth Layer x Month PC20-50-80

# Summarise 0.2, 0.5, 0.8 quantiles (cells in space) of monthly mean temperature in zone groups for
#     three depth ranges - epi, meso, bentho pelagic

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

# 1. Inputs ####

#    1.1 Directories ####
 ShapeRoot<-"/perm_storage/home/acon/Shapefiles/"
 ShapeName<-"HIMI Zones 2500.shp"
 BRANrootRasterReadStatic<-"/vsicurl/https://dapds00.nci.org.au/thredds/fileServer/gb6/BRAN/BRAN2020/static/"
 BRANrootRasterReadMonthData<-"/vsicurl/https://dapds00.nci.org.au/thredds/fileServer/gb6/BRAN/BRAN2020/month/ocean_temp_mth_"
 BRANrootNetCDFreadMonthData<-"https://dapds00.nci.org.au/thredds/dodsC/gb6/BRAN/BRAN2020/month/ocean_temp_mth_"

#    1.2 BRAN2022 Parameters ####
NAtemp        <- -32768  # missing values for temperature
Date.Origin   <- "1979-01-01"
Years         <-  seq(1993,2022,1)
Months        <-  sprintf("%02d",seq(1,12,1))

#    1.3 Spatial Scope ####
#Domain<-list(Coords=c(63,83,-57,-45) ,Name="HIMI")
Domain<-list(Coords=c(60,90,-60,-45) ,Name="Kerguelen Plateau")
AssessmentExtent<-ext(Domain$Coords)

#    1.4 Depth layers  ####

dStrata<-list( epipelagic    = c(0,100)     # Surface conditions (approximately in mixed layer) - surface 120m
               ,winter_water = c(110,200)   # representative of a subsurface cooler winter water 
               ,mesopelagic   = c(300,800)  # Mesopelagic conditions - (300-800m)
               ,benthopelagic = NA          # Bottom cell   NULL indicates whole water column; NA = bottom
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

#    1.5 Temporal Scope ####

# Periods for assessment 

PeriodYears<-list(
  Early  = c(1993:2002) # sequence winter to summer with last summer crossing into 2010
  ,Recent = c(2011:2020)# sequence winter to summer with last summer crossing into 2020
) # end periods

# seasons need to be in sequence for a comparison over a year.  
#if months are earlier than the first month in the first season then data from the next year will be assembled
SeasonMonths<-list(  
   Winter = c(4:9)   
  ,Summer = c(10:12,1:3)
) # end seasons



#    1.5 Regional Raster or Polygon analysis ####
doPolygons<-FALSE

#    1.6 Percentiles if polygon analysis  ####
zGrpTempPC<-c(0.2,0.5,0.8) # spatial percentiles of Temperature within zone groups

#    1.7 Output  ####
# file for saving output
#outName<-"/perm_storage/home/acon/MEASO/HIMI_Temp_Zgrp_Depth_Monthly_PC20-50-80 "
outName<-"/perm_storage/home/acon/MEASO/Kerg_PLat_Temp_Raster "

# 2. Functions #######

# 3. Preparation ####

#    3.1 BRAN Cell attributes #####

# x_T: geographic longitude, units - degrees E
# y_T: geographic latitude, units - degrees N
# area_T: Area of T_cell, units: m2

dfile<-paste(BRANrootRasterReadStatic,"grid_spec.nc",sep="")
r<-rast(dfile,"area_T")
rCellArea<-crop(r, AssessmentExtent)
CellNumbers<-row.names(as.data.frame(rCellArea,na.rm=F))  # for assigning to matrices then used in apply fns

# depth
dfile<-paste(BRANrootRasterReadStatic,"ocean_grid.nc",sep="")
# ht - ocean depth on t-cells
r<-rast(dfile,"ht")
rCellDepth<-crop(r, AssessmentExtent)

# kmt - number of depth levels on t-grid
r<-rast(dfile,"kmt")
rCellDepthLevels<-crop(r, AssessmentExtent)

dfCellAttr<-cbind(as.data.frame(rCellArea,na.rm=F),as.data.frame(rCellDepth,na.rm=F),as.data.frame(rCellDepthLevels,na.rm=F))

rCellAttr<-c(rCellArea,rCellDepth,rCellDepthLevels)
# plot(rCellAttr, col = hcl.colors(256))
crs(rCellAttr)<-"+proj=longlat +datum=WGS84"

rCellAttr <- terra::project(rCellAttr, CRSinUse)  # reproject to CRS being used for analyses

#    3.2 Zones (shapefiles) and Masks (if needed)  ######
if(doPolygons){
HIMIzones<-vect(paste(ShapeRoot,ShapeName,sep=""))
CRSinUse<-crs(HIMIzones)


# Zone Group Masks

# Zone                         Group

Zgrps<-list(rep(NULL,9))

# <300m around Heard Island
#  5            Heard Island     1
Zgrps[[1]]<-c(5)

# 0-500m (not including Heard)

# 12       South Shelf Inner     2
# 13       South Shelf Outer     2
#  6             North Shelf     2
Zgrps[[2]]<-c(6,12,13)

#  8              Shell Bank     3
Zgrps[[3]]<-c(8)

# 500 - 1000m

#  7         NorthEast Shelf     4
#  3          Eastern Trough     4
#  4    Eastern Trough South     4
Zgrps[[4]]<-c(3,4,7)

# 10           South Canyons     5
Zgrps[[5]]<-c(10)

# 18  Western Trough & Banks     6
Zgrps[[6]]<-c(18)

# > 1000m

#  1                    East     7
#  2               East Deep     7
# 19          Williams Ridge     7
# 20     Williams Ridge Deep     7
Zgrps[[7]]<-c(1,2,19,20)

#  9                   South     8
# 11              South Deep     8
# 14               SouthEast     8
# 15          SouthEast Deep     8
Zgrps[[8]]<-c(9,11,14,15)

# 16                    West     9
# 17               West Deep     9
Zgrps[[9]]<-c(16,17)


rZgrpMasks<-lapply(Zgrps,function(zg) {
  rTmp<-rasterize(HIMIzones[zg],rCellAttr,values=1)
  rTmp[!is.finite(rTmp)]<-NA
  return(rTmp)
})

dZgrpMasks<-sapply(Zgrps,function(zg) {
  rTmp<-rasterize(HIMIzones[zg],rCellAttr,values=1)
  rTmp[!is.finite(rTmp)]<-NA
  dTmp<-values(rTmp)
  return(!is.na(values(rTmp)))
})
row.names(dZgrpMasks)<-CellNumbers

# head(dZgrpMasks); nrow(dZgrpMasks); ncol(dZgrpMasks)
} # end if doPolygons

# 4. Depth Strata for averaging Temperature  ##################################

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

# 5. Run-time Parameters ##################################

if(doPolygons) nZgrps<-length(Zgrps)

# 6. Test routines prior to analysis  ##################################
doTests<-FALSE
#    6.1 extract cell areas for zones group 1 ####
if(doPolygons & doTests) plot(mask(rCellAttr,ZgrpMasks[[2]]))

#    6.2 plotting rasters ####
if(doTests){
FilePointerRaster<-paste(BRANrootRasterReadMonthData,FileDate[1],".nc",sep="")
r           <- rast(FilePointerRaster)
rMonthTemp  <- crop(r, AssessmentExtent)
plot(subset(rMonthTemp,36))
} # end doTests

# 7. Analysis - for each period ##################################


# If analysis by zones
# Result is a data.frame with the following columns
# zgrp (integer)
# pelagic layer (character)
# date
# 20 percentile
# 50 percentile
# 80 percentile

FileDate<-cbind(as.character(as.vector(sapply(Years,function(y)rep(y,12)))),rep(Months,length(Years)))
FileDate<-apply(FileDate,1,function(d) paste(d[1],"_",d[2],sep=""))

# header variables
# xt_ocean  (longitude)
# yt_ocean  (latitude)
# zt_ocean  (depth)
# Time (the time of the data)

# 9.1 sapply(FileDate
resPC<-lapply(FileDate[c(1:length(FileDate))],function(FileDate,depthMask,dZgrpMasks,dfCellAttr){

  FilePointerRaster<-paste(BRANrootRasterReadMonthData,FileDate,".nc",sep="")
  FilePointerNetCDF<-paste(BRANrootNetCDFreadMonthData,FileDate,".nc",sep="")
  
  #    9.1.1 Reading date of the netCDF directly ####
  nc <- RNetCDF::open.nc(FilePointerNetCDF)
  Time<-var.get.nc(nc, "Time")
  fDate<-as.Date(Time,Date.Origin)
print(fDate) # running check if there are problems
  #    9.1.2 Reading data as a raster ####
  
  r           <- rast(FilePointerRaster)
  rMonthTemp  <- crop(r, AssessmentExtent)
  
  dMTdf       <-as.data.frame(rMonthTemp,na.rm=F)  # note that missing values are dropped with default na.rm=T; come out as NA in vector
  dMT<-as.matrix(dMTdf)
  dMT[dMT==NAtemp]<-NA # make missing values = NA
  dimnames(dMT)<-NULL
  dimnames(dMT)[[1]]<-row.names(dMTdf) 
  rm(dMTdf)

  #head(dMT); nrow(dMT); ncol(dMT)
  
  #    9.1.3 cycle through depth x zone groups masks to obtain data and extract percentiles ####
  res<-lapply(depthMask,function(d,z,Temp,Cell){ # for each depth stratum
    
    if (sum(d)==0) {  # take temperature from deepest stratum in each cell 
      # 9.1.3.1 Bottom layer ####

      # 9.1.3.1.1 apply zgroup ####
      res<-apply(z,2,function(z,d,Temp,Cell){ # for each zone group
        dT          <- Temp[z,] # extract cells in zone
        dTrows      <- as.numeric(row.names(dT)) # row numbers to identify cells
        res1        <- rep(NA,length(dTrows))
        msk         <- !is.na(Cell[dTrows,"kmt"])
        
        dTcols<-Cell[dTrows,"kmt"]
        res1   <- sapply(seq(1,nrow(dT),1),function(r,d,c,m){
                  if(m[r] & !is.na(c[r])) return(d[r,c[r]]) else return(NA)
                       },dT,dTcols,msk) # read off bottom layer
        
        names(res1) <- dTrows # need cell numbers for extracting cell areas
                # 9.1.3.1.1.1 determine percentiles, weighted by cell area ####
        cellArea    <-Cell[dTrows,"area_T"]

        sdf         <-cbind(res1,cellArea)
        sdf         <-sdf[!is.na(sdf[,1]),]
        if(length(sdf)==2) return(rep(sdf[1],length(zGrpTempPC))) else {
          sdf<-sdf[order(sdf[,1]),] # sort by Temperature in ascending order
          cumArea<-sapply(seq(1,nrow(sdf),1),function(i,df){sum(df[c(1:i),2])},sdf)
          propArea<-cumArea/sum(sdf[,2]) # calculate cumulative area and make proportion
          sdf<-cbind(sdf,propArea)
          #      read area percentiles for temp according to vector of percentiles
          #               (determine index of row given each the percentile )
          res2<-sapply(seq(1,length(zGrpTempPC),1),function(i,df){
            #   print(c(i,nrow(df),(sum(zGrpTempPC[i]>df[,propArea])+1)))
            # print(df[,propArea])
            return(df[(sum(zGrpTempPC[i]>df[,"propArea"])+1),1])
          },sdf)
          return(res2)} # end else
      },d,Temp,Cell) # end apply dZgrpMasks  
      return(res)
      
    }  else {  # do non-benthic depths - mean temperature weighted by the depth interval in the depth stratum in each cell
      # 9.1.3.2 Layers other than Bottom ####

      # 9.1.3.2.1 apply zgroup ####
      res<-apply(z,2,function(z,d,Temp,Cell){ # for each zone group
        dT<-Temp[z,] # extract cells in zone and depth stratum
        dTmeans<-apply(dT,1,function(dT,d,DI){
          dTsub<-dT[d]
          DIsub<-DI[d]
          if(sum(!is.na(dTsub))>0){
            DIsub<-DIsub[!is.na(dTsub)]
            dTsub<-dTsub[!is.na(dTsub)]
            return(sum(dTsub*DIsub)/sum(DIsub))
          }  else {return(NA)}}
          ,d,depthBins_int)
        
        dTrows<-as.numeric(row.names(dT)) # row numbers to identify cells
        res1<-rep(NA,length(dTrows))
        msk<-!is.na(Cell[dTrows,"kmt"])
        res1[msk]<-dTmeans[msk] # read off temps
        names(res1)<-dTrows # need cell numbers for extracting cell areas
        # 9.1.3.2.1.1 function to determine percentiles, accounting for cell area
        cellArea<-Cell[dTrows,"area_T"]
        sdf<-cbind(res1,cellArea)
        sdf<-sdf[!is.na(sdf[,1]),]
        res2<-rep(NA,length(zGrpTempPC))
        if(length(sdf)>0){ # ensure at least one cell has data
          if(length(sdf)==2) res2<-rep(sdf[1],length(zGrpTempPC)) else {  # sdf matrix has two columns
            sdf<-sdf[order(sdf[,1]),] # sort by Temperature in ascending order
            cumArea<-sapply(seq(1,nrow(sdf),1),function(i,df){sum(df[c(1:i),2])},sdf)
            propArea<-cumArea/sum(sdf[,2]) # calculate cumulative area and make proportion
            sdf<-cbind(sdf,propArea)
            #      read area percentiles for temp according to vector of percentiles
            #               (determine index of row given each the percentile )
            res2<-sapply(seq(1,length(zGrpTempPC),1),function(i,df){
              #   print(c(i,nrow(df),(sum(zGrpTempPC[i]>df[,propArea])+1)))
              # print(df[,propArea])
              return(df[(sum(zGrpTempPC[i]>df[,"propArea"])+1),1])
            },sdf)
          }} # end if and end else
        return(res2)
      },d,Temp,Cell) # end apply dZgrpMasks  
      return(res)
    } # end else
  },dZgrpMasks,dMT,dfCellAttr) # end lapply depthMask
  # 9.1.4 create data.frame  #####
  dfPC<-do.call("rbind", 
                lapply(c(1:length(depthMask)),function(i,m,d,fDate){
                  dRep<-ncol(d[[i]])
                  if(!is.null(d[[i]])) {
                    dm <- data.frame(rep(fDate,dRep),rep(substr(names(m)[i],1,1),dRep),c(1:dRep),t(d[[i]]))
                    colnames(dm)<-c("Date","Depth","ZoneGroup","PC20","PC50","PC80")
                    return(dm)
                  } else return(NULL)
                },depthMask,res,fDate)
         ) # end do.call rbind (combines list of data frames)
  
  # 9.1.5 return result for the file
  return(dfPC)
  
},depthMask,dZgrpMasks,dfCellAttr)# end sapply FileDate
resPCdf<-do.call("rbind",resPC)

save(resPCdf,file=paste(outName,Sys.time(),".Rdata",sep=""))

