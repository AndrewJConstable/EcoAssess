# HIMI Chlorophyll a 
# surface density scaled up to cell area multiplied by MLD from BRAN 2020
# decadal means - summer, winter
# decades - 2000-2009; 2010-2019



# Important notes for processing BRAN data
# 1. Using Terra package to import netCDF files.
# 2. Rasters are complete x-y matrices, including when subsetting
# 3. When converting rasters to dataframes using as.data.frame, missing values are dropped using defaults.
#          Important to not remove cells with NA by setting na.rm=F i.e. as.data.frame(r, na.rm=F)
#          and manage NAs in code.

# 1. Libraries and Input Parameters  ##################################

library(raadtools)
library(raadfiles)
library(blueant)
library(RNetCDF)
library(ggplot2)
library(gplots)
library(terra)

#      1.1 Directories ####
 FileDir<-"/perm_storage/home/acon/"
 
 BRANrootRasterReadStatic<-"/vsicurl/https://dapds00.nci.org.au/thredds/fileServer/gb6/BRAN/BRAN2020/static/"
 BRANrootRasterReadMonthData<-"/vsicurl/https://dapds00.nci.org.au/thredds/fileServer/gb6/BRAN/BRAN2020/month/ocean_mld_mth_"
 BRANrootNetCDFreadMonthData<-"https://dapds00.nci.org.au/thredds/dodsC/gb6/BRAN/BRAN2020/month/ocean_mld_mth_"

#      1.2 Spatial Scope ####

# GIS set domain to be sampled (LongMin,LongMax,LatMin,LatMax)
 Domain<-list(Coords=c(63,83,-57,-45) ,Name="HIMI")
 AssessmentExtent<-ext(Domain$Coords)
 
 
#      1.3 Temporal Scope ####
 
    Years         <-  seq(2000,2023,1)
    refMonths     <-  sprintf("%02d",seq(1,12,1))

# Periods for assessment 

PeriodYears<-list(
   Early  = c(2003:2012) # sequence winter to summer with last summer crossing into 2010
  ,Recent = c(2013:2022)# sequence winter to summer with last summer crossing into 2020
  ) # end periods

# seasons need to be in sequence for a comparison over a year.  
#if months are earlier than the first month in the first season then data from the next year will be assembled
SeasonMonths<-list(  
   Spring = c(10,11,12)   
  ,Summer = c(1,2,3)
  ) # end seasons
 

#      1.4 BRAN2022 Parameters ####
NAmld        <- -1.00000002004088e+20  # missing values for temperature
Date.Origin   <- "1979-01-01"



#      1.5 Chlorophyll input data ####

FilePrefix<-"Chla"
fChlDaysStart<-as.Date("2002-07-01")  # earliest month available
fChlDaysEnd  <-as.Date("2023-06-30")   
fChlDates<-seq(fChlDaysStart,fChlDaysEnd,by="1 month")
# vectors for subsetting
YearsChlFiles  <- format(fChlDates,format="%Y")
MonthsChlFiles <- format(fChlDates,format="%m")

# reference years and months
refYears<-unique(Years)

# 2. Functions ########

# 3. BRAN Cell attributes and template  #####

#    3.1 Cell area ####
# x_T: geographic longitude, units - degrees E
# y_T: geographic latitude, units - degrees N
# area_T: Area of T_cell, units: m2

dfile<-paste(BRANrootRasterReadStatic,"grid_spec.nc",sep="")
r<-rast(dfile,"area_T")
rCellArea<-crop(r, AssessmentExtent)
CellNumbers<-row.names(as.data.frame(rCellArea,na.rm=F))  # for assigning to matrices then used in apply fns

#    3.2 Cell depth ####
dfile<-paste(BRANrootRasterReadStatic,"ocean_grid.nc",sep="")
# ht - ocean depth on t-cells
r<-rast(dfile,"ht")
rCellDepth<-crop(r, AssessmentExtent)

#    3.3 Cell depth levels ####
# kmt - number of depth levels on t-grid
r<-rast(dfile,"kmt")
rCellDepthLevels<-crop(r, AssessmentExtent)

dfCellAttr<-cbind(as.data.frame(rCellArea,na.rm=F),as.data.frame(rCellDepth,na.rm=F),as.data.frame(rCellDepthLevels,na.rm=F))

#    3.4 spatRaster Template ####

# develop raster template

rTemplate<-rCellDepthLevels
names(rTemplate)<-"Chla"
crs(rTemplate) <- "epsg:4326"


# 4. File Names for each Decade x Season period  #####

#    4.1 Files for generating rasters - periods x seasons ####

nPeriods<-length(PeriodYears)
nSeasons<-length(SeasonMonths)
month_1 <- unlist(SeasonMonths)[1]
Names_Period_By_Season<-NULL
Period_YrsMths<-NULL
FilesYrMth<-NULL
for(p in seq(1,nPeriods,1)){
  for(s in seq(1,nSeasons,1)){
       Names_Period_By_Season <- c(Names_Period_By_Season,paste(names(PeriodYears)[p],names(SeasonMonths)[s],sep="_") )
       Period_YrsMths<-c(Period_YrsMths,list(Years = length(PeriodYears[[p]]), Months = length(SeasonMonths[[s]])))
       
       fPS<-lapply(PeriodYears[[p]],function(y,Mths){
                  yrs<-rep(y,length(Mths))
                  yrs[Mths<month_1]<-y+1
                  f<-cbind(yrs,refMonths[Mths])
              colnames(f)<-c("Year","Month")
              return(f)}
              ,SeasonMonths[[s]]) # end lapply for each year in period
       names(fPS)<-paste("Y",PeriodYears[[p]],sep="")
    FilesYrMth<-c(FilesYrMth,list(fPS))
  }} # end loops
names(FilesYrMth)<-Names_Period_By_Season
names(Period_YrsMths)<-Names_Period_By_Season

#     4.2 MLD ####
perFilesMLD<-lapply(FilesYrMth,function(p){
                     lapply(p,function(y){
                       apply(y,1,function(d) paste(d[1],"_",d[2],sep=""))})
                     })

#     4.3 Chlorophyll  ####

returnDate<-function(d){
  which((YearsChlFiles%in%d[1] & MonthsChlFiles%in%d[2]))
} # end return date

perFilesChla<-lapply(FilesYrMth,function(p){
  lapply(p,function(y){apply(y,1,returnDate)})})



# 5. Designating outputs ####

#      5.1 general ####

outName<-"/perm_storage/home/acon/MEASO/HIMI_Chla_10-S-W_rbrick "

#      5.2 specific runs
AdjustByMLD <- TRUE
OnlyMLD <- FALSE # overrides chlorophyll calculations

saveRasterEachPeriod<-FALSE

# temporary adjustment of MLD files to save mean state for use on chlorophyll in summer of last year

   # perFilesMLD[[4]]<-perFilesMLD[[4]][c(1:(length(perFilesMLD[[4]])-1))]

# 5. Analysis - for each period x season ############


Chl_From_Chla_x_MLD<-function(fMLD,fChla,AdjustByMLD,OnlyMLD){
  if(OnlyMLD) print(paste("MLD : ",fMLD,sep="")) else if (AdjustByMLD) {
    print(paste("MLD : ",fMLD,";   Chla : ",fChla,sep=""))} else print(paste("Chla : ",fChla,sep=""))

  if(AdjustByMLD | OnlyMLD) {
       # Step 1 extract MLD
    # read file & make into data frame
    if(as.numeric(substring(fMLD,1,4))>2022){
      FilePointerRaster<-"MEASO/HIMI_MLD Recent_Summer.tiff"  # mean for the previous 9 years
    } else FilePointerRaster<-paste(BRANrootRasterReadMonthData,fMLD,".nc",sep="")
    r           <- rast(FilePointerRaster)
    crs(r)<-"epsg:4326"
    rMonthMLD  <- crop(r, AssessmentExtent)
    values(rMonthMLD)[values(rMonthMLD)==NAmld]<-NA   # make missing values = NA
    } # end if use MLD
  if(!OnlyMLD){
  # Step 2 extract Chla
  subChla<-rast(readCHL_month(fChlDates[fChla], xylim = extent(Domain$Coords)))  # raad tools to get chlorophyll  - note it uses package raster so need to convert to spatRaster
  crs(subChla)<-"epsg:4326"
  } # end if not OnlyMLD
  
  if(OnlyMLD) rMonthChla<-rMonthMLD else if (AdjustByMLD) rMonthChla<-project(subChla,rMonthMLD) else rMonthChla<-subChla
  # Step 3 MLD x Chla      
  if(OnlyMLD) rMonthC<-rMonthMLD else if(AdjustByMLD) rMonthC<-rMonthMLD*rMonthChla else rMonthC<-rMonthChla
  
# Step 4 convert to matrix and save      
dChldf <-  as.data.frame(rMonthC,na.rm=F)  # note that missing values are not dropped with default na.rm=T; come out as NA in vector
dChl   <-  as.matrix(dChldf)
rm(dChldf,rMonthChla,rMonthC)
return(dChl)

} # end Chl_From_Chla_x_MLD

meanChl<-function(d){mean(d,na.rm=TRUE)}


resPeriod<-lapply(seq(1,length(Names_Period_By_Season)),function(p,fM,fC,AdjustByMLD,OnlyMLD){
  print(paste("Period ",p,sep=""))
  if(AdjustByMLD | OnlyMLD)  FilesMLD<-fM[[p]] else FilesMLD<-NULL # list of years with each year have files for the months to be averaged
  FilesChla<-fC[[p]]

  dPeriodChl<-lapply(c(1:length(FilesMLD)),function(y,fMLD,fChla,AdjustByMLD,OnlyMLD){
               fMthsMLD<-fMLD[[y]]
               fMthsChla<-fChla[[y]]
               dYearChl<-lapply(c(1:length(fMthsMLD)),function(m,fMLD,fChla,AdjustByMLD,OnlyMLD){
                             Chl_From_Chla_x_MLD(fMLD[m],fChla[m],AdjustByMLD,OnlyMLD)  
                              },fMthsMLD,fMthsChla,AdjustByMLD,OnlyMLD)
               return(apply(do.call(cbind, dYearChl),1,meanChl))
               },FilesMLD,FilesChla,AdjustByMLD,OnlyMLD)
  rChlmean<-rTemplate
  values(rChlmean)<-apply(do.call(cbind, dPeriodChl),1,meanChl)
  
  if(OnlyMLD) rName<-paste(names(fM[p]),"_MLD",sep="") else rName<-paste(names(fM[p]),"_Chla",sep="")
  names(rChlmean)<-rName
  
  outfile<-paste(outName,rName,Sys.time(),".tiff",sep="")
  print(class(rChlmean))
  print(outfile)
  if(saveRasterEachPeriod) writeRaster(rChlmean, filename=outfile)
  return(rChlmean)
  ########  
},perFilesMLD,perFilesChla,AdjustByMLD,OnlyMLD) # end lapply periods


resPeriod<-do.call("c",resPeriod)
print(resPeriod)

writeRaster(resPeriod, filename=paste(outName,Sys.time(),".tiff",sep=""))

print("File Saved and Finished")
