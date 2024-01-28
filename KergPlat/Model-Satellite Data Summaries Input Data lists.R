# Model-Satellite Data lists for different analyses

# 3. Kerguelen Temperature ####
InputData <-list(
  useBRAN      = TRUE
  ,useChla      = FALSE
  ,useDomain    = "Kerguelen_Plateau"
  ,useBRANvar   = "dBRAN_Temp_mth"
  ,doPolygons   = FALSE
  
  ,outVariableName = "Temp"
  
  ,PeriodYears = list(  Early  = c(1992:2001) # sequence winter to summer with last summer crossing into 2010
                        ,Recent = c(2012:2021)# sequence winter to summer with last summer crossing into 2020
  ) # end periods
  # seasons need to be in sequence for a comparison over a year.  
  #if months are earlier than the first month in the first season then data from the next year will be assembled
  ,SeasonMonths=list(  
    Winter = c(4:9)   
    ,Summer = c(10:3)
  ) # end seasons
  # Depth layers for assessment
  ,dStrata      = list( epipelagic    = c(0,100)     # Surface conditions (approximately in mixed layer) - surface 120m
                        ,winter_water = c(110,200)   # representative of a subsurface cooler winter water 
                        ,mesopelagic   = c(300,800)  # Mesopelagic conditions - (300-800m)
                        ,benthopelagic = NA          # Bottom cell   NULL indicates whole water column; NA = bottom
  ) # end list
  ,dBRAN        = list(staticRaster = dBRANall$dBRAN_static
                       ,mthRaster   = dBRANall[[useBRAN]]$Raster
                       ,mthNetCDF   = dBRANall[[useBRAN]]$NetCDF)
  ,fns  = list(fnInputFiles = 
              ,   )
) # end InputKerguelenTemp



