## Mosaic and reproject to 1km Copernicus GLC fractional lc data
library(raster)
#1. Mosaic tiles into one file (taking into account their position in multiple folders)
dirs <- list.files("C:\\work_resources\\Features\\landcover\\Copernicus_GLC_2015", full.names = T)
files <- list.files(dirs[1], pattern = "tree-coverfraction-layer_EPSG-4326.tif$", full.names = T)
r1 <- raster(files)
for (i in 2:8)
{
  Sys.sleep(0.1)
  print(i)
  flush.console()
  files <- list.files(dirs[i], pattern = "tree-coverfraction-layer_EPSG-4326.tif$", full.names = T)
  r2 <- raster(files)
  r1 <- merge(r1, r2)
}
writeRaster(r1, filename = "C:\\PhD\\ModellingViviparity\\FullDist_Sal_Lsal\\egv_prep\\roughMask\\forestCover_temp2.tif")

#2. Resample to coarsest resolution and strictest extent
high.res <- raster("C:\\PhD\\ModellingViviparity\\FullDist_Sal_Lsal\\egv_prep\\roughMask\\forestCover_temp.tif")
 


#Version 2 (files resampled in ArcGIS and stored in the same folder)
library(raster)
#1. Mosaic tiles into one file (taking into account their position in multiple folders)
lst <- list.files("C:\\PhD\\ModellingViviparity\\FullDist_Sal_Lsal\\egv_prep\\roughMask\\forestCover_upscale", pattern = "\\.tif$", full.names = T)

r1 <- raster(lst[1])
for (i in 2:8)
{
  Sys.sleep(0.1)
  print(i)
  flush.console()
  r2 <- raster(lst[i])
  r1 <- merge(r1, r2, tolerance = 1)
}
writeRaster(r1, filename = "C:\\PhD\\ModellingViviparity\\FullDist_Sal_Lsal\\egv_prep\\roughMask\\forestCover_upscale\\forestCover.tif")
