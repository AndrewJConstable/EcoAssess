# MEASO BRAN2022 Temp means - grp[yr] grp[mth] grp[depth].R

# function to generate raster of mean temperatures across a group of years, months and depths
# in following sequence
#. 1. weighted mean temperature between min and max depth (can use mixed layer depth as either minimum or maximum)
#  2. mean temperature across months as designated in input vector
#. 3. mean temperature across years as designated in input vector

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

fnBRAN_TempMean<-function(
    Yrs       # vector of years to average 
    ,Mths      # vector of months to average
    ,Dpths     # depth c(min,max). - note positive values from shallow to deep. NA = use MLD
    ,SplitDpthByMLD # logical to split Dpths by MLD if both depths are entered in Dpths
    ,dfDepths  # dataframe for BRAN depth layers [Depth,Interval,Min,Max]
){ # begin function
  
  readMLD<-FALSE; MLDelement<-0    #  element is for placing MLD in Dpths vector
  if(sum(is.na(Dpths))>0)  {readMLD<-TRUE; MLDelement<-{if(is.na(Dpths[2])) 2 else 1}; SplitDpthByMLD<-FALSE}
  if(SplitDpthByMLD) readMLD<-TRUE
  
  resYrs<-lapply(Yrs,function(y,Mths,Dpths,dfDepths,readMLD,MLDelement,SplitDpthByMLD){   #  loop year
    resMths<-lapply(Mths,function(m,y,Dpths,dfDepths,readMLD,MLDelement,SplitDpthByMLD){ #     loop month
      
      #       read file for temperature and create df
      rTemp<-rast(paste0(flinkBRAN_Temp,sprintf("%04d",y),"_",sprintf("%02d",m),".nc"))
      crs(rTemp)<-"epsg:4326"
      rTemp  <- crop(rTemp, ext(Arena))
      values(rTemp)[values(rTemp)==NAtemp]<-NA   # make missing values = NA
      
      dfTemp  <-as.data.frame(rTemp,na.rm=F)  # note that missing values are dropped with default na.rm=T; come out as NA in vector
      VarNames<-names(dfTemp)
      VarNames<-paste0("Depth_",sapply(sapply(VarNames,function(t){substring(t,15,(unlist(gregexpr(".",t,fixed=TRUE))[1]-1))})
                                       ,function(t){p<-unlist(gregexpr("_",t,fixed=TRUE))[1];if(p<0) p<-nchar(t) else p<-(p-1);
                                       substring(t,1,p)}))
      dTemp   <-as.matrix(dfTemp)
      dimnames(dTemp)<-NULL
      dimnames(dTemp)<-list(row.names(dfTemp),VarNames)
      rm(dfTemp)
      
      #       if readMLD then read file for MLD  
      if(readMLD){
        rMLD<-rast(paste0(flinkBRAN_MLD,sprintf("%04d",y),"_",sprintf("%02d",m),".nc"))
        crs(rMLD)<-"epsg:4326"
        rMLD  <- crop(rMLD, ext(rTemp))
        dfMLD  <-as.data.frame(rMLD,na.rm=F)  # note that missing values are dropped with default na.rm=T; come out as NA in vector
        dMLD   <-as.matrix(dfMLD)
        dimnames(dMLD)<-NULL
        dimnames(dMLD)[[1]]<-row.names(dfMLD) 
        rm(dfMLD)
      } # end readMLD
      
      res<-lapply(as.numeric(row.names(dTemp)),function(Cell,dT,dD,Dpths,dfDepths,readMLD,SplitDpthByMLD,MLDelement){
        
        DpthStrata<-matrix(Dpths,byrow=TRUE,nrow=1)
        if(readMLD){if(SplitDpthByMLD){
          DpthStrata<-rbind(c(Dpths[1],dD[Cell]),c(dD[Cell],Dpths[2]))
        } else {DpthStrata[1,MLDelement]<-dD[Cell]
        }} # end if reaDMLD, SplitDpthByMLD
        # Do by DpthStrata
        Temp<- apply(DpthStrata,1,function(depthRange){
          # return NA if one of the elements of depthRange is NA
          if(sum(is.na(depthRange))>0) return(NA)
          
          # create mask and adjust min-max depth intervals
          
          depthMask<-dfDepths[,"Min"]>=depthRange[1] & dfDepths[,"Min"]<depthRange[2]
          depthWts<-dfDepths[depthMask,"Interval"]
          depthWts[1]<-dfDepths[dfDepths[,"Min"]<=depthRange[1] & dfDepths[,"Max"]>depthRange[1],"Max"]-depthRange[1]
          depthWts[length(depthWts)]<- depthRange[2]-dfDepths[dfDepths[,"Min"]<=depthRange[2] & dfDepths[,"Max"]>depthRange[2],"Min"]
          
          # calculate weighted mean
          return(sum(depthWts*dT[Cell,depthMask],na.rm=TRUE)/sum(depthWts,na.rm=TRUE))
        }) # end apply DpthStrata
        Temp<-matrix(Temp,nrow=1)
        if(ncol(Temp)==1) dimnames(Temp)[[2]]<-"Temp" else dimnames(Temp)[[2]]<-c("Temp_MLD","Temp_Deep")
        Temp<-cbind(as.data.frame(Cell),as.data.frame(Temp))
        return(Temp)  
      },dTemp,dMLD,Dpths,dfDepths,readMLD,SplitDpthByMLD,MLDelement ) # end sapply
      res<-do.call("rbind",res)
      return(res)
      
    },y,Dpths,dfDepths,readMLD,MLDelement,SplitDpthByMLD) # end sapply Mths
    
    # list of df for each month (rows=cells, cols = cells, Temp/Temp_MLD, ?Temp_Deep) - need to account for two strata if split by MLD
    Strata<-seq(1,(if(SplitDpthByMLD) 2 else 1),1)
    res<-do.call(cbind,lapply(Strata,function(s,resMths){
      apply(sapply(seq(1,length(Mths)),function(m,tmp,s){tmp[[m]][,(s+1)]},resMths,s),1,mean) # note (s+1) accounts for Cell column
    },resMths))
    
    return(res)
  },Mths,Dpths,dfDepths,readMLD,MLDelement,SplitDpthByMLD) # end Yrs
  
  # list of df for each year (rows=cells, cols = Temp/Temp_MLD, ?Temp_Deep) - need to account for two strata if split by MLD
  Strata<-seq(1,(if(SplitDpthByMLD) 2 else 1),1)
  resMns<-do.call(cbind,lapply(Strata,function(s,resYrs){
    apply(sapply(seq(1,length(Yrs)),function(y,tmp,s){tmp[[y]][,s]},resYrs,s),1,mean)
  },resYrs))
  
  # return raster brick (Temp/Temp_MLD, ?Temp_Deep)
  # create template from first month,year MLD
  rTemplate<-rast(paste0(flinkBRAN_MLD,sprintf("%04d",Yrs[1]),"_",sprintf("%02d",Mths[1]),".nc"))
  crs(rTemplate)<-"epsg:4326"
  rTemplate  <- crop(rTemplate, ext(Arena))
  rTempMns<-NULL
  for (s in Strata){
    values(rTemplate)<-resMns[,s]
    rTempMns<-c(rTempMns,rTemplate)  
  }
  names(rTempMns)<-if(SplitDpthByMLD) c("Temp_MLD","Temp_Deep") else "Temp"
  return(rast(rTempMns))
} # end function

#########################################################################################################
# 3. Input Parameters
##################################

# 1.1 Directories
 BRANrootRasterReadStatic<-"/vsicurl/https://dapds00.nci.org.au/thredds/fileServer/gb6/BRAN/BRAN2020/static/"
 BRANrootRasterReadMonthData<-"/vsicurl/https://dapds00.nci.org.au/thredds/fileServer/gb6/BRAN/BRAN2020/month/ocean_temp_mth_"
 BRANrootNetCDFreadMonthData<-"https://dapds00.nci.org.au/thredds/dodsC/gb6/BRAN/BRAN2020/month/ocean_temp_mth_"

# 1.2 BRAN2022 Parameters
NAtemp        <- -32768  # missing values for temperature
Date.Origin   <- "1979-01-01"
Years         <-  seq(1993,2020,1)
Months        <-  sprintf("%02d",seq(1,12,1))
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


Dpths<-c(0,1000) # use MLD when NA or when SplitDpthByMLD==TRUE
SplitDpthByMLD<-TRUE

# GIS
Arena<-c(0,360,-80,-40) # West, East, South, North boundaries

outName<-"/perm_storage/home/acon/MEASO/Temp-tmp"

#########################################################################################################
# 2. set up BRAN data
##################################

# 4,1 Cell Attributes dfCellAttr[area_T,ht,kmt] ####

# area_T: Area of T_cell, units: m2
# ht: m depth - note POSITIVE values
# kmt : number of depth layers


# x_T: geographic longitude, units - degrees E
# y_T: geographic latitude, units - degrees N

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

# 4.2 Set dataframe for depth layers - dfDepth[Depth, Interval] ####

# determine depth intervals of each depth bin for weighting calculations of mean etc.
# i. determine intervals between elements in dBins_BRAN
NdBins<-length(dBins_BRAN)
dBins_diff<-c(dBins_BRAN[1],dBins_BRAN[2:NdBins]-dBins_BRAN[1:(NdBins-1)])

# ii. loop through the vector of differences to calculate the interval    
prev<-dBins_diff[1]
depthBins_int<-dBins_diff*0; depthMin<-depthBins_int; depthMax<-depthBins_int
depthBins_int[1]<-prev*2; depthMin[1]<-dBins_BRAN[1]-prev; depthMax[1]<-dBins_BRAN[1]+prev
for (d in seq(2,NdBins,1)) {prev<-dBins_diff[d]-prev; depthBins_int[d]<-prev*2; depthMin[d]<-dBins_BRAN[d]-prev; depthMax[d]<-dBins_BRAN[d]+prev}

# iii. calculate minimum and maximum of depth interval for generating masks

# iv. combine into dataframe and clean up
dfDepths<-data.frame(Depth=dBins_BRAN,Interval=depthBins_int,Min=depthMin,Max=depthMax)
rm(NdBins,dBins_diff,depthBins_int,prev)

# BRAN File Prefixes ####
flinkBRAN_MLD<-"/vsicurl/https://dapds00.nci.org.au/thredds/fileServer/gb6/BRAN/BRAN2020/month/ocean_mld_mth_"
flinkBRAN_Temp<-"/vsicurl/https://dapds00.nci.org.au/thredds/fileServer/gb6/BRAN/BRAN2020/month/ocean_temp_mth_"

Yrs<-seq(2011,2020,1)
Mths<-seq(1,12,1)
mthNames<-sprintf("%02d",Mths)
for(m in Mths){
  print(paste0("Start Time - Month ",mthNames[m],": ",Sys.time()))
rTemp<-fnBRAN_TempMean(Yrs,m,Dpths,SplitDpthByMLD,dfDepths)
writeRaster(rTemp,file=paste0("rTemp_",mthNames[m],".tif"),overwrite=TRUE)
} # end m


