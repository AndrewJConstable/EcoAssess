# Model-Satellite Data functions for Kerguelen Temperature
# BRAN 2022 HIMI Temp Zgroup x Depth Layer x Month PC20-50-80

# Summarise 0.2, 0.5, 0.8 quantiles (cells in space) of monthly mean temperature in zone groups for
#     three depth ranges - epi, meso, bentho pelagic

# Important notes for processing BRAN data
# 1. Using Terra package to import netCDF files.
# 2. Rasters are complete x-y matrices, including when subsetting
# 3. When converting rasters to dataframes using as.data.frame, missing values are dropped using defaults.
#          Important to not remove cells with NA by setting na.rm=F i.e. as.data.frame(r, na.rm=F)
#          and manage NAs in code.



# First update Input Data as needed
InputData$dBRAN<-list(staticRaster = dBRANall$dBRAN_static
                      ,mthRaster   = dBRANall$dBRAN_Temp_mth$Raster
                      ,mthNetCDF   = dBRANall$dBRAN_Temp_mth$NetCDF)
  
generateFileNames<-function(FilesYrMth){
    return(list(Temp = lapply(FilesYrMth,function(p){
          lapply(p,function(y){
            apply(y,1,function(d) paste(d[1],"_",d[2],sep=""))})
})))
} # end function 


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

###############################################
# Analysis function  # note that this requires only one month in a period.
# written as if periods are consecutive years and seasons within years comprise of one month each.

fnAnalysis<-function(p,perFiles,InData){
  depthMask<-InData$depthMask
  dZgrpMasks<-InData$dZgrpMasks
  dfCellAttr<-InData$dfCellAttr
  
  fBRAN      <-perFiles$Temp
  dBRAN      <-InData$dBRAN
  BRAN_lastYear<-InData$BRAN_lastYear

  FilePointerRaster<-paste(dBRAN$mthRaster,fBRAN[p],".nc",sep="")
  FilePointerNetCDF<-paste(dBRAN$mthNetCDF,fBRAN[p],".nc",sep="")

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
  
}  # end fnAnalysis
