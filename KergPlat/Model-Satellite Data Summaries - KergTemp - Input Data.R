# Model-Satellite Data lists for different analyses

# 3. Kerguelen Temperature ####
InputData <-list(
  
   AnalysisTitle          = "Kerguelen Temperature"
  ,useBRAN                = TRUE
  ,useChla                = FALSE
  ,useBRANvar             = "dBRAN_Temp_mth"
  ,BRAN_lastYear          = 2022
  ,doPolygons             = FALSE
  ,saveRasterEachPeriod   = FALSE
  ,doDiffs                = TRUE  # difference between two layers
  ,outVariableName        = "Temp"
  ,outputRasterEachPeriod = TRUE

  ,ReCalculateData        = FALSE # for recalculating all foundation data from BRAN e.g. dfCellAttr
  
  ,outFileRoot          = "/perm_storage/home/acon/MEASO/Kerg_Temp_10-W-S_rbrick "
  
  ,Domain      = c(60,90,-60,-45)  
  ,PeriodYears = list(  Early  = c(1993:2002) # sequence winter to summer with last summer crossing into 2010
                        ,Recent = c(2012:2021)# sequence winter to summer with last summer crossing into 2020
  ) # end periods
  # seasons need to be in sequence for a comparison over a year.  
  #if months are earlier than the first month in the first season then data from the next year will be assembled
  ,SeasonMonths=list(  
     Winter = c(4:9)   
    ,Summer = c(10:12,1:3)
  ) # end seasons
  # Depth layers for assessment
  ,dStrata      = list( epipelagic    = c(0,100)     # Surface conditions (approximately in mixed layer) - surface 120m
                        ,winter_water = c(110,200)   # representative of a subsurface cooler winter water 
                        ,mesopelagic   = c(300,800)  # Mesopelagic conditions - (300-800m)
                        ,benthopelagic = NA          # Bottom cell   NULL indicates whole water column; NA = bottom
                         ) # end list
  ,diffDepths           = c(1,3)  # determine differences between depth strata.  Change order changes sign
  ,diffDepthsName       = "Diff"
) # end InputKerguelenTemp



