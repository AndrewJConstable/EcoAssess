# HIMI Chlorophyll a only
# surface density scaled up to cell area (using gebco resolution)  (not multiplied by MLD from BRAN 2020)
# decadal means - summer, winter
# decades - 2000-2009; 2010-2019

#########################################################################################################
# 1. Libraries and Input Parameters
##################################

library(raadtools)
library(raadfiles)
library(blueant)
library(RNetCDF)
library(ggplot2)
library(gplots)
library(terra)

# 1.1 Directories
 FileDir<-"/perm_storage/home/acon/"

refMonths<-c("01","02","03","04","05","06","07","08","09","10","11","12")
 
DecadeEarly<-c(2003:2012) # sequence winter to summer with last summer crossing into 2010
DecadeRecent<-c(2013:2022)# sequence winter to summer with last summer crossing into 2020
DecadeYears<-list(DecadeEarly,DecadeRecent)
nYrsInPeriod<-length(DecadeEarly)

WinterMths<-c(1,2,3)  # April-September
SummerMths<-c(10,11,12)   #  ,1,2,3) # October-March
SeasonMonths<-list(WinterMths,SummerMths)
nMthsInSeason<-length(WinterMths)

###   Chlorophyll input data
FilePrefix<-"Chla"
DaysStart<-as.Date("2002-07-01")  # earliest month available
DaysEnd  <-as.Date("2023-06-30")   
dates<-seq(DaysStart,DaysEnd,by="1 month")
# vectors for subsetting
Years<-format(dates,format="%Y")
Months<-format(dates,format="%m")

# reference years and months
refYears<-unique(Years)

# GIS set domain to be sampled (LongMin,LongMax,LatMin,LatMax)
Domain<-list(Coords=c(63,83,-57,-45) ,Name="HIMI")
AssessmentExtent<-ext(Domain$Coords)

outName<-"/perm_storage/home/acon/MEASO/HIMI_Chla_10-S-W_rbrick "

# 2. Functions ##################################

meanChl_from_rStack<-function(r){
  t<-terra::mean(r)
#  t<-subst(t,NaN,0)
  print(t)
  return(t)
} # end function

meanChl_log_from_rStack<-function(r){
  t<-exp(mean(log(r),na.rm=TRUE))
  t<-subst(t,NaN,0)
  return(t)
} # end function

meanChl<-function(r){
  t<-mean(r,na.rm=TRUE)
  if(is.nan(t)) return(NA) else return(t)
} # end function


# 3. Generate raster for assessment area from a BRAN tiff ##################################
#GEBCO21<-rast(readtopo("gebco_21"))
#GEBCO21_HIMI<-crop(GEBCO21,AssessmentExtent)
#rAreaKm2<-cellSize(GEBCO21_HIMI, unit="km")
#rm(GEBCO21)

rTemplate<-rast("MEASO/rSummer-Top.tif")
names(rTemplate)<-"Chla"

# 4. File Names for each Decade x Season period ######

# Files for generating rasters - decade (early, recent) x season (winter, summer)

dEarly_WinterFiles<-cbind(
   as.character(unlist(lapply(DecadeEarly,function(y,n){rep(y,n)},length(WinterMths))))
  ,refMonths[rep(WinterMths,length(DecadeEarly))])
colnames(dEarly_WinterFiles)<-c("Year","Month")

dEarly_SummerFiles<-cbind(
  as.character(unlist(lapply(DecadeEarly,function(y,m){yrs<-rep(y,length(m)); yrs[m<10]<-y+1; return(yrs)},SummerMths)))
 ,refMonths[rep(SummerMths,length(DecadeEarly))])
colnames(dEarly_SummerFiles)<-c("Year","Month")


dRecent_WinterFiles<-cbind(
  as.character(unlist(lapply(DecadeRecent,function(y,n){rep(y,n)},length(WinterMths))))
  ,refMonths[rep(WinterMths,length(DecadeRecent))])
colnames(dRecent_WinterFiles)<-c("Year","Month")

dRecent_SummerFiles<-cbind(
  as.character(unlist(lapply(DecadeRecent,function(y,m){yrs<-rep(y,length(m)); yrs[m<10]<-y+1; return(yrs)},SummerMths)))
  ,refMonths[rep(SummerMths,length(DecadeRecent))])
colnames(dRecent_SummerFiles)<-c("Year","Month")

perFilesMLD<-list(
  eWinter = apply(dEarly_WinterFiles,1,function(d) paste(d[1],"_",d[2],sep=""))
 ,eSummer = apply(dEarly_SummerFiles,1,function(d) paste(d[1],"_",d[2],sep=""))
 ,rWinter = apply(dRecent_WinterFiles,1,function(d) paste(d[1],"_",d[2],sep=""))
 ,rSummer = apply(dRecent_SummerFiles,1,function(d) paste(d[1],"_",d[2],sep=""))
)

returnDate<-function(d){
     which((Years%in%d[1] & Months%in%d[2]))
     } # end return date

perFilesChla<-list(    # in format of a date
   eWinter = apply(dEarly_WinterFiles,1,returnDate)
  ,eSummer = apply(dEarly_SummerFiles,1,returnDate)
  ,rWinter = apply(dRecent_WinterFiles,1,returnDate)
  ,rSummer = apply(dRecent_SummerFiles,1,returnDate)
) # end list

  
PeriodNames<-c("eWinter","eSummer","rWinter","rSummer")


# 5. Analysis - for each decade x season period,  ######################

# by period
#resPeriod<-lapply(seq(1,length(PeriodNames)),function(p,Name,fM,fC){
resPeriod<-lapply(c(1,3),function(p,Name,fM,fC){
   print(paste("Period ",p,sep=""))
  
    FilesChla<-fC[[p]]
print(FilesChla)

  # by year
    dYearChl<-NULL
    # rYearChl<-NULL
    for(y in c(1:nYrsInPeriod)){  
      dSeasonChl<-NULL
      # rSeasonChl<-NULL
       for(m in c(1:nMthsInSeason)){ 
         File_n<-(y-1)*nMthsInSeason+m
         print(dates[FilesChla[File_n]])
         
         # Step 2 extract Chla
              subChla<-rast(readCHL_month(dates[FilesChla[File_n]], xylim = extent(Domain$Coords)))  # raad tools to get chlorophyll  - note it uses package raster so need to convert to spatRaster
              rMonthChla<-project(subChla,rTemplate)
              rm(subChla)
              rMonthC<-rMonthChla  # here for code with MLD
        # Step 4 convert to matrix and save      
            dChldf <-  as.data.frame(rMonthC,na.rm=F)  # note that missing values are not dropped with default na.rm=T; come out as NA in vector
            dChl   <-  as.matrix(dChldf)
            rm(dChldf,rMonthChla,rMonthC)
            dSeasonChl<-cbind(dSeasonChl,dChl)
             # rSeasonChl<-c(rSeasonChl,rMonthC)      
    } # end m
      
      #       after processing each file, calculate three mean values for each cell for the season to be saved in respective matrices
      #                          i) mean temperature in the top layer
      #                         ii) mean temperature in the middle layer
      #                        iii) for each month calculate the difference between top and middle layers and then take the mean difference for the season
     dYearChl<-cbind(dYearChl,apply(dSeasonChl,1,meanChl))
      #rYearChl<-c(rYearChl,mean(rSeasonChl))
    } # end y
    rChlmean<-rTemplate
    values(rChlmean)<-apply(dYearChl,1,meanChl)
    names(rChlmean)<-paste(Name[p],"_Chla",sep="")
     return(rChlmean)
########  
    },PeriodNames,perFilesMLD,perFilesChla) # end lapply periods


resPeriod<-do.call("c",resPeriod)

writeRaster(resPeriod, filename=paste(outName,Sys.time(),".tiff",sep=""))

print("File Saved and Finished")

######## plot raster #######

getRasterColours<-function(d,breaks,colours){
  col_pal <- colorRampPalette(colours)
  Cmin    <- breaks[1]
  Cmax    <- breaks[length(breaks)]
  Cn      <- 1000

  rCols   <- col_pal(Cn)
  rColour <- rep(NA,length(d))
  Cvalid  <- (!is.nan(d) & !is.na(d))
  rColour[Cvalid] <- rCols[floor((d[Cvalid]-Cmin)/(Cmax-Cmin)*Cn)]
  return(rColour)
} # end get colours

rColours<-col2hex(c("blue","yellow","orange","lightgreen","darkgreen"))
rColBreaks<-c(0,0.1,0.2,0.5,1)

tmp<-resPeriod
values(tmp)[values(tmp)>1.2]<-1.2
hist(values(tmp))

#coltab(tmp)<-getRasterColours(values(tmp),rColBreaks,rColours)
plot(tmp,col=getRasterColours(values(tmp),rColBreaks,rColours))

