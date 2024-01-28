# Routines to summarise model and satellite data
#
# by Andrew J Constable
#
# 1. Background information ####
#
# Outputs as summaries in rasters (cell summaries) or dataframes (polygons)
# Factors are 'Period' and 'Season'
#      Periods comprise a sequence of years with minimum vectors of 1 year
#      Seasons comprise a sequence of months with minimum vectors of one month
#
# Requires functions (used in lapply) for
#    Month (any processing of one or more rasters to deliver output for one month)
#    Year (any processing to combine months in a season to deliver output for one year)
#    Period (summarising the year outputs into the summary for the period)
# Each function requires
#         a numeric index for which period, year or month
#         all inputs in a list
# Output of a function is a single object suitable for use in the parent function
#
# Important notes for processing BRAN data
# 1. Using Terra package to import netCDF files.
# 2. Rasters are complete x-y matrices, including when subsetting
# 3. When converting rasters to dataframes using as.data.frame, missing values are dropped using defaults.
#          Important to not remove cells with NA by setting na.rm=F i.e. as.data.frame(r, na.rm=F)
#          and manage NAs in code.

# 1. Libraries and Functions ####

library(raadtools)
library(raadfiles)
library(blueant)
library(RNetCDF)
library(ggplot2)
library(gplots)
library(terra)

# 2. Functions ########

#    2.1 Function to generate list of PeriodxSeason dataframe (Year,Month)
fnPeriod_x_SeasonTextForFilenames<-function(PeriodYears,SeasonMonths){
  nPeriods<-length(PeriodYears)
  nSeasons<-length(SeasonMonths)
  Names_Period_By_Season<-NULL
#  Period_YrsMths<-NULL   ??????????? not found elsewhere in code - delete?
  FilesYrMth<-NULL
  for(p in seq(1,nPeriods,1)){
      YrMonth_1<-SeasonMonths[[1]][1]
    for(s in seq(1,nSeasons,1)){
      SeasonMonth_1 <- SeasonMonths[[s]][1]
      
      Names_Period_By_Season <- c(Names_Period_By_Season,paste(names(PeriodYears)[p],names(SeasonMonths)[s],sep="_") )
 #     Period_YrsMths<-c(Period_YrsMths,list(Years = length(PeriodYears[[p]]), Months = length(SeasonMonths[[s]])))
      
      fPS<-lapply(PeriodYears[[p]],function(y,Mths,YrMonth_1,SeasonMonth_1){
        yrs<-rep(y,length(Mths))
        yrs[Mths<SeasonMonth_1 | Mths<YrMonth_1]<-y+1
        f<-cbind(yrs,refMonths[Mths])
        colnames(f)<-c("Year","Month")
        return(f)}
        ,SeasonMonths[[s]],YrMonth_1,SeasonMonth_1) # end lapply for each year in period
      names(fPS)<-paste("Y",PeriodYears[[p]],sep="")
      FilesYrMth<-c(FilesYrMth,list(fPS))
    }} # end loops
  names(FilesYrMth)<-Names_Period_By_Season
#  names(Period_YrsMths)<-Names_Period_By_Season
  return(FilesYrMth)#list(FilesYrMth,Period_YrsMths))
} # end fnPeriodSeasonTextForFilenames

# 3. Initialisation Parameters  ##################################
#      3.1 Directories ####
Analysis<-"KergSST"

dRoot<-"/perm_storage/home/acon/"
setwd(dRoot)
dBRANall<-list(
           dBRAN_static = list(Raster = "/vsicurl/https://dapds00.nci.org.au/thredds/fileServer/gb6/BRAN/BRAN2020/static/")
          ,dBRAN_MLD_mth = list(Raster = "/vsicurl/https://dapds00.nci.org.au/thredds/fileServer/gb6/BRAN/BRAN2020/month/ocean_mld_mth_"
                                  ,NetCDF = "https://dapds00.nci.org.au/thredds/dodsC/gb6/BRAN/BRAN2020/month/ocean_mld_mth_"
                                  ) # end MLD list
          ,dBRAN_Temp_mth = list(Raster = "/vsicurl/https://dapds00.nci.org.au/thredds/fileServer/gb6/BRAN/BRAN2020/month/ocean_temp_mth_"
                                  ,NetCDF = "https://dapds00.nci.org.au/thredds/dodsC/gb6/BRAN/BRAN2020/month/ocean_temp_mth_"
                                  ) # end MLD list
) # end dBran list


#      3.3 Reference Months ####

     refMonths     <-  sprintf("%02d",seq(1,12,1))

#      3.4 Chlorophyll input data ####

#FilePrefix<-"Chla"
#fDateStart<-as.Date("2002-07-01")  # earliest month available
#fDateEnd  <-as.Date("2023-06-30")   

FilePrefix<-"SST"
fDateStart<-as.Date("1992-07-01")  # earliest month available
fDateEnd  <-as.Date("2023-06-30")   

#      3.5 BRAN global parameters ####

# Depth layers

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

# Missing Values

NAtemp        <- -32768  # missing values for temperature
NAmld        <- -1.00000002004088e+20  # missing values for temperature

# Date of origin
Date.Origin   <- "1979-01-01"


# 4. Inputs for Analyses (always after Initialisation Parameters####
source(paste("MEASO/Model-Satellite Data Summaries - ",Analysis," - Input Data.R",sep=""), local = environment()) # need to ensure the local environment is what is being used when 
                                                                                                                  # running in background
source(paste("MEASO/Model-Satellite Data Summaries - ",Analysis," - Fns.R",sep=""), local = environment())

# 5. Pre-processing of Inputs (note that globally unique variables do not need general assignment) ####

#     5.1 Extent of Assessment ####
AssessmentExtent<-ext(InputData$Domain)


#     5.2 BRAN Cell attributes and template (if needed) ####
if(InputData$useBRAN){

  #    3.1 Cell area  # this serves as the template as well as part of cell attributes (hence outside the wrapper to avoid calculating data)
  # x_T: geographic longitude, units - degrees E
  # y_T: geographic latitude, units - degrees N
  # area_T: Area of T_cell, units: m2
  
  dfile<-paste(dBRANall$dBRAN_static$Raster,"grid_spec.nc",sep="")
  r<-rast(dfile,"area_T")
  rCellArea<-crop(r, AssessmentExtent)
  CellNumbers<-row.names(as.data.frame(rCellArea,na.rm=F))  # for assigning to matrices then used in apply fns
  
  
  
# wrapper to save recalculating data    
   # objectName<-"dfCellAttr"
    # if(!exists(objectName)){  # wrapper line 1
      # if(!file.exists(paste(dRoot,"MEASO/",objectName,".rds",sep=""))) { # wrapper line 2
#  end first part of wrapper
  
        
#    3.2 Cell depth
dfile<-paste(dBRANall$dBRAN_static$Raster,"ocean_grid.nc",sep="")
# ht - ocean depth on t-cells
r<-rast(dfile,"ht")
rCellDepth<-crop(r, AssessmentExtent)

#    3.3 Cell depth levels
# kmt - number of depth levels on t-grid
r<-rast(dfile,"kmt")
rCellDepthLevels<-crop(r, AssessmentExtent)

dfCellAttr<-cbind(as.data.frame(rCellArea,na.rm=F),as.data.frame(rCellDepth,na.rm=F),as.data.frame(rCellDepthLevels,na.rm=F))

# end wrapper to avoid redoing calculations
  # saveRDS(dfCellAttr,paste(dRoot,"MEASO/",objectName,".rds",sep="")) # dataframe
    # } else {
  # readRDS(paste(dRoot,"MEASO/",objectName,".rds",sep=""))
    # } # end if file.exists
  # } # end if exists
  # end wrapper 

} # end input useBRAN

# 6. If summarising by polygons - generate spatVector of polygons ####

if(InputData$doPolygons){

zPolygons <- vect(InputData$PolygonFile)
pZones<-as.data.frame(values(zPolygons))[,"Zone"]
mPmasks<-lapply(as.data.frame(values(zPolygons))[,"ID"],function(p,zP) {
    rTmp<-rasterize(zP[values(zP["ID"])==p],rCellArea,values=1)
    rTmp[!is.finite(rTmp)]<-NA
    dTmp<-values(rTmp)
    return(!is.na(values(rTmp)))
  },zPolygons)
mPmasks<-do.call(cbind,mPmasks)
dimnames(mPmasks)[[1]]<-CellNumbers
dimnames(mPmasks)[[2]]<-  pZones
} # end if doPolygons


# 6. File Names for each Decade x Season period  #####

FilesYrMth<-fnPeriod_x_SeasonTextForFilenames(InputData$PeriodYears,InputData$SeasonMonths)
perFiles<-generateFileNames(FilesYrMth)
Names_Period_By_Season<-names(FilesYrMth)

# 7. Analysis - for each period x season ############

#     7.1 get spatRaster Template ####

# develop raster template
if(InputData$useBRAN){
  rTemplate<-rCellArea
  names(rTemplate)<-InputData$outVariableName
  crs(rTemplate) <- "epsg:4326"
} else { # end if useBRAN
  rTemplate<-get_rTemplate(1) # from the first raster of raad data
  names(rTemplate)<-InputData$outVariableName
}

#     7.2 analyses ####

resPeriod<-lapply(seq(1,length(Names_Period_By_Season)),fnAnalysis,perFiles,InData)

if(InputData$doPolygons){ 
  resPeriod<-do.call(rbind,resPeriod)
  save(resPeriod,file=paste(InputData$outFileRoot,Sys.time(),".Rdata",sep=""))
  } else {
    resPeriod<-do.call("c",resPeriod)
    writeRaster(resPeriod, filename=paste(InputData$outFileRoot,Sys.time(),".tiff",sep=""))
}

print(resPeriod)
print("File Saved and Finished")
