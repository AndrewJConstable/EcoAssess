# Model-Satellite Data lists for different analyses

# 3. Kerguelen Chlorophyll with MLD ####
InputData <-list(
   AnalysisTitle        = "Kerguelen Chla with MLD"
  ,useBRAN              = TRUE
  ,useChla              = TRUE
  ,AdjustByMLD          = TRUE
  ,OnlyMLD              = FALSE   # do this for generating rasters for use with chlorophyll when MLD missing
  ,useBRANvar           = "dBRAN_MLD_mth"
  ,BRAN_lastYear        = 2022
  ,BRAN_summerMLDtif    = "MEASO/Kerg_Summer_MLD_mean_2013-2022.tiff"
  ,doPolygons           = TRUE
  ,PolygonFile          = "Shapefiles/Zones Merged 2.shp"
  ,saveRasterEachPeriod = FALSE
  
  ,outVariableName = "Chla"
  ,outFileRoot     = "/perm_storage/home/acon/MEASO/Kerg_Chla_10-S-S_rbrick "

  ,Domain      = c(60,90,-60,-45)  
  ,PeriodYears = list(  Early  = c(2003:2012) 
                        ,Recent = c(2013:2022)
  ) # end periods
  # seasons need to be in sequence for a comparison over a year. Can cross a split year.  Seasons need to be in sequence across split year 
  #if months are earlier than the first month in the first season then data from the next year will be assembled
  ,SeasonMonths=list(  
     Spring = c(10:12)   
    ,Summer = c(1:3)
  ) # end seasons
) # end InputKerguelenTemp



