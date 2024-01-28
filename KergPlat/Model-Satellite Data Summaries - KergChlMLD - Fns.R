# Model-Satellite Data functions for Kerguelen Chla with MLD

# MUST HAVE

#        File: "MEASO/HIMI_MLD Recent_Summer.tiff" for use when BRAN does not have most recent year with Chl.


# First update Input Data as needed
InputData$dBRAN<-list(staticRaster = dBRANall$dBRAN_static
                      ,mthRaster   = dBRANall$dBRAN_MLD_mth$Raster
                      ,mthNetCDF   = dBRANall$dBRAN_MLD_mth$NetCDF)

InData<-list(
   AdjustByMLD = InputData$AdjustByMLD
  ,OnlyMLD    = InputData$OnlyMLD
  ,dBRAN      = InputData$dBRAN
  ,BRAN_lastYear = InputData$BRAN_lastYear
  ,BRAN_summerMLDtif = InputData$BRAN_summerMLDtif # with average summer conditions for recent when no BRAN data
) # end list

#     5.2 Chlorophyll file info ####

fChlDates<-seq(fDateStart,fDateEnd,by="1 month")
# vectors for subsetting
YearsChlFiles  <- format(fChlDates,format="%Y")
MonthsChlFiles <- format(fChlDates,format="%m")


# function for generating string for chlorophyll files

returnDate<-function(d){
  which((YearsChlFiles%in%d[1] & MonthsChlFiles%in%d[2]))
} # end return date


#     Period files
generateFileNames<-function(FilesYrMth){
  
  return(list(MLD = lapply(FilesYrMth,function(p){
                      lapply(p,function(y){
                       apply(y,1,function(d) paste(d[1],"_",d[2],sep=""))})
                       })
               ,Chla = lapply(FilesYrMth,function(p){
                             lapply(p,function(y){apply(y,1,returnDate)})})
           ))# end perFiles list
} # end function


  #################################
  Chl_From_Chla_x_MLD<-function(fMLD,fChla,AdjustByMLD,OnlyMLD,dBRAN,BRAN_lastYear){
    if(OnlyMLD) print(paste("MLD : ",fMLD,sep="")) else if (AdjustByMLD) {
      print(paste("MLD : ",fMLD,";   Chla : ",fChla,sep=""))} else print(paste("Chla : ",fChla,sep=""))
    
    if(AdjustByMLD | OnlyMLD) {
      # Step 1 extract MLD
      # read file & make into data frame
      if(as.numeric(substring(fMLD,1,4))>BRAN_lastYear){ 
        if(OnlyMLD) FilePointerRaster<-NULL else FilePointerRaster<-InData$BRAN_summerMLDtif  # mean for the previous 9 years
      } else FilePointerRaster<-paste(dBRAN$mthRaster,fMLD,".nc",sep="")
      
      if(!is.null(FilePointerRaster)){
      r           <- rast(FilePointerRaster)
      crs(r)<-"epsg:4326"
      rMonthMLD  <- crop(r, AssessmentExtent)
      values(rMonthMLD)[values(rMonthMLD)==NAmld]<-NA   # make missing values = NA
      } # end if not is null
      } # end if use MLD
    
    
    if(!OnlyMLD){
      # Step 2 extract Chla
      subChla<-rast(readCHL_month(fChlDates[fChla], xylim = extent(InputData$Domain)))  # raad tools to get chlorophyll  - note it uses package raster so need to convert to spatRaster
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

#################################
meanChl<-function(d){mean(d,na.rm=TRUE)}


#############################
fnAnalysis<-function(p,perFiles,InData){

  fBRAN         <-perFiles$MLD
  fC            <-perFiles$Chla
  AdjustByMLD   <-InData$AdjustByMLD
  OnlyMLD       <-InData$OnlyMLD
  dBRAN         <-InData$dBRAN
  BRAN_lastYear<-InData$BRAN_lastYear
  
  print(paste("Period ",p,sep=""))
  
  if(AdjustByMLD | OnlyMLD)  FilesMLD<-fBRAN[[p]] else FilesMLD<-NULL # list of years with each year have files for the months to be averaged
  FilesChla<-fC[[p]]
 
  # for each year in period 
  dPeriodChl<-lapply(c(1:length(FilesMLD)),function(y,fMLD,fChla,AdjustByMLD,OnlyMLD){
    fMthsMLD<-fMLD[[y]]
    fMthsChla<-fChla[[y]]
    # for each season in year
    dYearChl<-lapply(c(1:length(fMthsMLD)),function(m,fMLD,fChla,AdjustByMLD,OnlyMLD){
      Chl_From_Chla_x_MLD(fMLD[m],fChla[m],AdjustByMLD,OnlyMLD,dBRAN,BRAN_lastYear)  
    },fMthsMLD,fMthsChla,AdjustByMLD,OnlyMLD) # end season
    return(apply(do.call(cbind, dYearChl),1,meanChl))
  },FilesMLD,FilesChla,AdjustByMLD,OnlyMLD) # end year

  if(InputData$doPolygons){ # do dataframe for period
    mChl<-sapply(dPeriodChl,function(Chl){
                  df<-sapply(pZones,function(pZ,d){
                            return(sum((d[mPmasks[,pZ]]*dfCellAttr[mPmasks[,pZ],"area_T"]),na.rm=TRUE))
                         },Chl)
                  })
    dfChl<-as.data.frame(mChl)
    pNum<-rep(p,length(pZones))
    dfChl<-cbind(pNum,pZones,dfChl)
    return(dfChl)
      
  } else { # else do raster for period
       rChlmean<-rTemplate
       values(rChlmean)<-apply(do.call(cbind, dPeriodChl),1,meanChl) # take means across years to generate raster for period
  
       if(OnlyMLD) rName<-paste(names(fBRAN[p]),"_MLD ",sep="") else rName<-paste(names(fBRAN[p]),"_Chla ",sep="")
       names(rChlmean)<-rName
       if(InputData$saveRasterEachPeriod) outfile<-paste(InputData$outFileRoot,rName,Sys.time(),".tiff",sep=""); writeRaster(rChlmean, filename=outfile)
       return(rChlmean)
  } # end else     
  ########  
}   # perFilesMLD,perFilesChla,AdjustByMLD,OnlyMLD) # end lapply periods



