# Model-Satellite Data functions for Kerguelen Temperature

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
depthMask<-lapply(InputData$dStrata,function(d){
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

meanNArmTrue<-function(v) mean(v,na.rm=TRUE)


InData<-list(depthMask = depthMask
            ,dBRAN = InputData$dBRAN)
###############################################
# Analysis function
fnAnalysis<-function(p,perFiles,InData){
  fBRAN      <- perFiles$Temp
  depthMask  <- InData$depthMask
  dBRAN      <- InData$dBRAN
 
  resPeriod<-lapply(seq(1,length(fBRAN[[p]])),function(y,yFiles,f){

    resMonthDepths<-lapply(c(1:length(yFiles[[y]])),function(m,mFiles){
           print(mFiles[[m]])
      # Step 1. read file & make into data frame
           FilePointerRaster<-paste(dBRAN$mthRaster,mFiles[[m]],".nc",sep="") 

           r           <- rast(FilePointerRaster)
           rMonthTemp  <- crop(r, AssessmentExtent)
      
           df       <-as.data.frame(rMonthTemp,na.rm=F)  # note that missing values are dropped with default na.rm=T; come out as NA in vector
           matT       <-as.matrix(df)
           matT[matT==NAtemp]<-NA # make missing values = NA
           dimnames(matT)   <-NULL
           dimnames(matT)[[1]]<-row.names(df) 
           rm(df)

      # Step 2. for each cell, mean temperatures in each depth layer
           resDepth<-lapply(depthMask,function(dM,matT){ # for each depth stratum
             
                   if (sum(dM)==0) {  # take temperature from deepest stratum in each cell 
                     matTrows      <- as.numeric(row.names(matT)) # row numbers to identify cells
                     res1        <- rep(NA,length(matTrows))
                     msk         <- !is.na(dfCellAttr[matTrows,"kmt"])
                     
                     matTcols<-dfCellAttr[matTrows,"kmt"]
                     res1   <- sapply(seq(1,nrow(matT),1),function(r,d,c,m){
                       if(m[r] & !is.na(c[r])) return(d[r,c[r]]) else return(NA)
                     },matT,matTcols,msk) # read off bottom layer
                     names(res1)<-row.names(matT)
                     return(res1)
                   } else { # determine mean temp from range of depth cells in depth stratum

                   Tsub<-matT[,dM]
                   DIsub<-depthBins_int[dM]
                   return(apply(Tsub,1,function(dT,dI){
                                 msk<-!is.na(dT)
                                if(sum(msk>0)) return(sum(dT[msk]*dI[msk])/sum(dI[msk]))  else return(NA)
                         },DIsub))
                   } # end if sum(d)
              },matT) # end res
       # carry forward the list - res  
           resDepth<-do.call(cbind,resDepth)
        if(InputData$doDiffs) resDepth<-cbind(resDepth,resDepth[,InputData$diffDepths[1]]-resDepth[,InputData$diffDepths[2]])
       return(resDepth)
    },yFiles[[y]]) # end resYear  lapply years in period
  # Step 3. calculate monthly mean for each depth strata
    
    resDepthYear<-lapply(seq(1,ncol(resMonthDepths[[1]]),1),function(d,resMD){
             res<-do.call(cbind,lapply(resMD,function(Y,d){Y[,d]},d))
             return(apply(res,1,meanNArmTrue))
             },resMonthDepths) # end lapply resDepthYear
    return(do.call(cbind,resDepthYear)) # matrix for year with cols = depth strata
  },fBRAN[[p]]) # lapply end years in period   

  # Step 4 generate list of rasters of means at each depth 

  rNames<-names(InputData$dStrata)
  if(InputData$doDiffs) rNames<-c(rNames,InputData$diffDepthsName)
  rNames<-sapply(rNames,function(d,p,t){paste(p,"-",d,"-",t,sep="")},names(fBRAN)[p],InputData$outVariableName)
  
  resDepthPeriod<-lapply(seq(1,ncol(resPeriod[[1]]),1),function(d,resPD,rNames){
    res<-do.call(cbind,lapply(resPD,function(P,d){P[,d]},d)) # extract all years for depth
    rRes<-rTemplate
    values(rRes)<-apply(res,1,meanNArmTrue)
    names(rRes)<-rNames[d]
    return(rRes)
    },resPeriod,rNames) # end lapply resDepthYear
  resDepthPeriod<-do.call("c",resDepthPeriod)
  if(InputData$outputRasterEachPeriod){ # print rasters individually
    writeRaster(resDepthPeriod, filename=paste(InputData$outFileRoot,"-",names(fBRAN)[p]," ",Sys.time(),".tiff",sep=""))
       } # end print rasters
  
  return(resDepthPeriod)
  ##############
}  # ,PeriodNames,perFiles) # end lapply periods   # end fnAnalysis
