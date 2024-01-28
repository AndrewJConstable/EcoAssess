# HIMI Mixed Layer Depth BRAN2020 Temp decadal mean (top-mid) - summer & winter

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
 BRANrootRasterReadMonthData<-"/vsicurl/https://dapds00.nci.org.au/thredds/fileServer/gb6/BRAN/BRAN2020/month/ocean_mld_mth_"
 BRANrootNetCDFreadMonthData<-"https://dapds00.nci.org.au/thredds/dodsC/gb6/BRAN/BRAN2020/month/ocean_mld_mth_"

# 1.2 BRAN2022 Parameters
NAmld        <- -1.00000002004088e+20  # missing values for temperature
Date.Origin   <- "1979-01-01"
Years         <-  seq(1993,2020,1)
Months        <-  sprintf("%02d",seq(1,12,1))

DecadeEarly<-c(1993:2002) # sequence winter to summer with last summer crossing into 2003
DecadeRecent<-c(2010:2019)# sequence winter to summer with last summer crossing into 2020

WinterMths<-c(4,5,6,7,8,9)  # April-September
SummerMths<-c(10,11,12,1,2,3) # October-March

# GIS
Arena<-c(60,90,-60,-45) # West, East, South, North boundaries


outName<-"/perm_storage/home/acon/MEASO/HIMI_MLD_10-S-W_rbrick "
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
names(rTemplate)<-"MLD"

# by period
resPeriod<-lapply(seq(1,length(PeriodNames)),function(p,Name,f){
  Files<-f[[p]]
  # by year
     dYearMLD<-NULL
    for(y in c(1:10)){  
      dSeasonMLD<-NULL
       for(m in c(1:6)){ 
         File_n<-(y-1)*6+m
         print(Files[File_n])
         
         # read file & make into data frame
         FilePointerRaster<-paste(BRANrootRasterReadMonthData,Files[File_n],".nc",sep="")
                   # delete? FilePointerNetCDF<-paste(BRANrootNetCDFreadMonthData,FileDate,".nc",sep="")
         
         r           <- rast(FilePointerRaster)
         rMonthMLD  <- crop(r, ext(Arena))
         
         dMLDdf       <-as.data.frame(rMonthMLD,na.rm=F)  # note that missing values are not dropped with default na.rm=T; come out as NA in vector
         dMLD<-as.matrix(dMLDdf)
         dMLD[dMLD==NAmld]<-NA # make missing values = NA
         rm(dMLDdf)


         dSeasonMLD<-cbind(dSeasonMLD,dMLD)
    } # end m
      
      #       after processing each file, calculate three mean values for each cell for the season to be saved in respective matrices
      #                          i) mean temperature in the top layer
      #                         ii) mean temperature in the middle layer
      #                        iii) for each month calculate the difference between top and middle layers and then take the mean difference for the season
      dYearMLD<-cbind(dYearMLD,apply(dSeasonMLD,1,mean))
      } # end y
  # after processing all years, take thee mean of each of the three metrics in each cell and
  # add them to raster template as values
     rMLDmean<-rTemplate
     values(rMLDmean)<-apply(dYearMLD,1,mean)
     names(rMLDmean)<-paste(Name[p],"_MLD",sep="")
     return(rMLDmean)
########  
    },PeriodNames,perFiles) # end lapply periods
resPeriod<-do.call("c",resPeriod)
print(resPeriod)

writeRaster(resPeriod, filename=paste(outName,Sys.time(),".tiff",sep=""))

print("File Saved and Finished")
