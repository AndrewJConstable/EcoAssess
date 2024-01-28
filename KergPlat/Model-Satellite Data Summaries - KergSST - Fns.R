# Model-Satellite Data functions for Kerguelen SST

InData<-list(NULL)

#     5.2 SST file info ####

fSSTDates<-seq(fDateStart,fDateEnd,by="1 month")
# vectors for subsetting
YearsSSTFiles  <- format(fSSTDates,format="%Y")
MonthsSSTFiles <- format(fSSTDates,format="%m")


# function for generating string for SST files

returnDate<-function(d){
  which((YearsSSTFiles%in%d[1] & MonthsSSTFiles%in%d[2]))
} # end return date


#     Period files
generateFileNames<-function(FilesYrMth){
   fNames<-lapply(FilesYrMth,function(p){
     lapply(p,function(y){apply(y,1,returnDate)})})
  return(list(SST = fNames))# end perFiles list
} # end function


getSST<-function(fSST){
print(paste("SST : ",fSST,sep=""))
    subSST<-rast(readsst(fSSTDates[fSST], time.resolution = "monthly", xylim = extent(InputData$Domain)))  # raad tools to get chlorophyll  - note it uses package raster so need to convert to spatRaster
    crs(subSST)<-"epsg:4326"
   dfSST <-  as.data.frame(subSST,na.rm=F)  # note that missing values are not dropped with default na.rm=T; come out as NA in vector
  return(as.matrix(dfSST))
} # end getSST

get_rTemplate<-function(fSST){
    subSST<-rast(readsst(fSSTDates[fSST], time.resolution = "monthly", xylim = extent(InputData$Domain)))  # raad tools to get chlorophyll  - note it uses package raster so need to convert to spatRaster
    crs(subSST)<-"epsg:4326"
    return(subSST)
  } # end getSST


#################################
meanSST<-function(d){mean(d,na.rm=TRUE)}


#############################
fnAnalysis<-function(p,perFiles,InData){
  fC            <-perFiles$SST

  print(paste("Period ",p,sep=""))
  
  FilesSST<-fC[[p]]
 
  # for each year in period 
  dPeriodSST<-lapply(c(1:length(FilesSST)),function(y,fSST){
    fMthsSST<-fSST[[y]]
    # for each season in year
    dYearSST<-lapply(c(1:length(fMthsSST)),function(m,fSST){
      getSST(fSST[m])  
    },fMthsSST) # end season
    return(apply(do.call(cbind, dYearSST),1,meanSST))
  },FilesSST) # end year

  if(InputData$doPolygons){ # do dataframe for period
    mSST<-sapply(dPeriodSST,function(SST){
                  df<-sapply(pZones,function(pZ,d){
                            return(sum((d[mPmasks[,pZ]]*dfCellAttr[mPmasks[,pZ],"area_T"]),na.rm=TRUE))
                         },SST)
                  })
    dfSST<-as.data.frame(mSST)
    pNum<-rep(p,length(pZones))
    dfSST<-cbind(pNum,pZones,dfSST)
    return(dfSST)
      
  } else { # else do raster for period
       rSSTmean<-rTemplate
       values(rSSTmean)<-apply(do.call(cbind, dPeriodSST),1,meanSST) # take means across years to generate raster for period
  
      rName<-paste(names(fC[p]),"_SST ",sep="")
       names(rSSTmean)<-rName
       if(InputData$saveRasterEachPeriod) {
           print("InsideOutfile")
            outfile<-paste(InputData$outFileRoot,rName,Sys.time(),".tiff",sep="")
            writeRaster(rSSTmean, filename=outfile)
            }
       return(rSSTmean)
  } # end else     
  ########  
}   # perFilesSST



