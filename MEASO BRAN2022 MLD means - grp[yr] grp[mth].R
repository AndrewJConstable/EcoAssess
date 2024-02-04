# MEASO BRAN2022 MLD means - grp[yr] grp[mth].R

# function to generate raster of mean mixed layer depths across a group of years and months
# in following sequence
#  1. mean temperature across months as designated in input vector
#. 2. mean temperature across years as designated in input vector

# Important notes for processing BRAN data
# 1. Using Terra package to import netCDF files.
# 2. Rasters are complete x-y matrices, including when subsetting
# 3. When converting rasters to dataframes using as.data.frame, missing values are dropped if the defaults are used.
#          Important to not remove cells with NA by setting na.rm=F i.e. as.data.frame(r, na.rm=F)
#          and manage NAs in code.

#########################################################################################################
# 1. Libraries
##################################

library(RNetCDF)
library(terra)


#########################################################################################################
# 3. Function
##################################

fnBRAN_MLDmean<-function( # return raster of means
    Yrs       # vector of years to average 
    ,Mths      # vector of months to average
){ # begin function

  resYrs<-lapply(Yrs,function(y,Mths,Dpths,dfDepths,readMLD,MLDelement,SplitDpthByMLD){   #  loop year
    resMths<-lapply(Mths,function(m,y,Dpths,dfDepths,readMLD,MLDelement,SplitDpthByMLD){ #     loop month
        rMLD<-rast(paste0(flinkBRAN_MLD,sprintf("%04d",y),"_",sprintf("%02d",m),".nc"))
        crs(rMLD)<-"epsg:4326"
        rMLD  <- crop(rMLD, ext(Arena))
        dfMLD  <-as.data.frame(rMLD,na.rm=F)  # note that missing values are dropped with default na.rm=T; come out as NA in vector
        dMLD   <-as.matrix(dfMLD)
        dimnames(dMLD)<-NULL
        dimnames(dMLD)[[1]]<-row.names(dfMLD) 
        rm(dfMLD)
       return(dMLD)
    },y) # end lapply Mths
    resMths<-do.call(cbind,resMths)
    resMean<-apply(resMths,1,mean,na.rm=TRUE)
    return(resMean)
  },Mths) # end lapply Yrs
  resYrs<-do.call(cbind,resYrs)
  resMean<-apply(resYrs,1,mean,na.rm=TRUE)

  # return raster
  # create template from first month,year MLD
  rTemplate<-rast(paste0(flinkBRAN_MLD,sprintf("%04d",Yrs[1]),"_",sprintf("%02d",Mths[1]),".nc"))
  crs(rTemplate)<-"epsg:4326"
  rTemplate  <- crop(rTemplate, ext(Arena))
  values(rTemplate)<-resMean
  
  
  return(rTemplate)
} # end function

#########################################################################################################
# 3. Input Parameters
##################################


# GIS
Arena<-c(0,360,-80,-40) # West, East, South, North boundaries

outName<-"/perm_storage/home/acon/MEASO/MLD"

#########################################################################################################
# 2. set up BRAN data
##################################

# BRAN File Prefixes ####
flinkBRAN_MLD<-"/vsicurl/https://dapds00.nci.org.au/thredds/fileServer/gb6/BRAN/BRAN2020/month/ocean_mld_mth_"
flinkBRAN_Temp<-"/vsicurl/https://dapds00.nci.org.au/thredds/fileServer/gb6/BRAN/BRAN2020/month/ocean_temp_mth_"

Yrs<-seq(2011,2020,1)
Mths<-seq(1,12,1)
mthNames<-sprintf("%02d",Mths)
rMLD<-NULL
for(m in Mths){
  print(paste0("Start Time - Month ",mthNames[m],": ",Sys.time()))
r<-fnBRAN_MLDmean(Yrs,m)
names(r)<-mthNames[m]
rMLD<-c(rMLD,r)
} # end m
writeRaster(rast(rMLD),file="rMLD.tif",overwrite=TRUE)


