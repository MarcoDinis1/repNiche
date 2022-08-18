########### Alternative version using RNetCDF package ###############################################################################################

install.packages("RNetCDF")              #! 'RNetCDF' package substitution
library(RNetCDF)
library(raster)

#Enter lat and lon ranges # 
lat.range=c(60.0025109626,30.0941776293)        #! Ranges instead of point values. Order does not matter
lon.range=c(-10.6838986827,52.1494346506)

# enter in variable you want to download see: http://thredds.northwestknowledge.net:8080/thredds/terraclimate_aggregated.html
var="PDSI" #casing matters. Link above takes you to index. CLicking on links within tells you the proper acronym for each variable

baseurlagg <- paste0(paste0("http://thredds.northwestknowledge.net:8080/thredds/dodsC/agg_terraclimate_",var),"_1958_CurrentYear_GLOBE.nc")


nc <- open.nc(baseurlagg)
lon <- var.get.nc(nc, "lon")
lat <- var.get.nc(nc, "lat")
lat.range <- sort(lat.range)                              #!sort user input values from low to high
lon.range <-sort(lon.range)
lat.index <- which(lat>=lat.range[1]&lat<=lat.range[2])    #! index values within specified range
lon.index <- which(lon>=lon.range[1]&lon<=lon.range[2])    
lat.n <- length(lat.index)                                #!value for count
lon.n <- length(lon.index)
start <- c(lon.index[1], lat.index[1], 1)
count <- c(lon.n, lat.n, NA)                            #! parameter change: 'NA' instead of '-1' to signify entire dimension


# read in the full period of record using aggregated files

data <-var.get.nc(nc, variable = var,start = start, count,unpack=TRUE)    


## Calculating rasters of average values across 1958-present)
#object created above contains all values for all cells across all years/months in a single array. The following steps sequentially extracts relevant values for each monthly raster (from 1958-present) and computes the sum. In the end, this is divided by the total number of monthly rasters to produce a mean. Done this way to avoid computational constraints from having >700 rasters in memory at once

ncells <- lat.n*lon.n
pos.start <- 1
pos.finish <- ncells

data.mat <- matrix(data =data[pos.start:pos.finish], nrow = lat.n, ncol = lon.n, byrow = T)
r <- raster(data.mat)
plot(r)

pos.start <- pos.start+ncells
pos.finish <- pos.finish+ncells
while (pos.finish < length(data))
{
  data.mat <- matrix(data =data[pos.start:pos.finish], nrow = lat.n, ncol = lon.n, byrow = T)
  r1 <- raster(data.mat)
  r <- r+r1
  
  pos.start <- pos.start+ncells
  pos.finish <- pos.finish+ncells
  rm(r1)
}
den <- length(data)/ncells
r <- r/den

#Set extent and crs of final raster before saving to disk
extent(r) <- c(-10.6838986827,52.1494346506,30.0941776293, 60.0025109626)
crs(r) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

writeRaster(r, "D:\\PhD\\ModellingViviparity\\FullDist_Sal_Lsal\\egv_prep\\TerraClimate\\vpd.tif")







