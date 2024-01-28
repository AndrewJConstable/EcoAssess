# Model-Satellite Data lists for different analyses

# SST ####
InputData <-list(
   AnalysisTitle        = "Kerguelen SST"
  ,useBRAN              = FALSE
  ,doPolygons           = FALSE
  ,PolygonFile          = "Shapefiles/Zones Merged 2.shp"
  ,saveRasterEachPeriod = FALSE
  
  ,outVariableName = "SST"
  ,outFileRoot     = "/perm_storage/home/acon/MEASO/Kerg_SST_rbrick "

  ,Domain      = c(60,90,-60,-45)  
  ,PeriodYears = list(  Early  = c(1993:2002) 
#                        ,Recent = c(2013:2022)
  ) # end periods
  # seasons need to be in sequence for a comparison over a year. Can cross a split year.  Seasons need to be in sequence across split year 
  #if months are earlier than the first month in the first season then data from the next year will be assembled
  ,SeasonMonths=list(  
     Summer = c(10:12,1:3)   
#     Summer = c(1:3)
  ) # end seasons
) # end InputKerguelenTemp



