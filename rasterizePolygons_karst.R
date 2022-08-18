#Rasterization of WOMAK shapefile data on karst aquifers. The purpose is to create a karst presence/absence raster at 850 x850m resolution which can then be converted to 1 x 1km proportion. I want to keep all values/aquifer clases for now and fill all non-karst cells with '0'.
library(rgdal)
library(raster)

#load shapefile to be rasterized
karst.shape <- readOGR(dsn="D:\\PhD\\ModellingViviparity\\FullDist_Sal_Lsal\\egv_prep\\roughMask\\shapes", layer="whymap_karst")

#load reference raster (template for extent/crs) 
rast <- raster("D:\\PhD\\ModellingViviparity\\FullDist_Sal_Lsal\\egv_prep\\TerraClimate\\aet.tif") 

#create raster template with desired resolution (used reference raster as template for everything except resolution which was manually defined to 2x2km)
rast2 <- raster()
extent(rast2) <- extent(rast)
crs(rast2) <- crs(rast)
origin(rast2) <- origin(rast)
res(rast2) <- c(0.02083334, 0.02083334)

#check proper positioning of shapefile and reference raster
plot(rast)
plot(karst.shape, add = T)

#rasterize. Assigns all values in shapefile to raster cell and treats all empty cells as 0
karst.rast <- rasterize(karst.shape, rast2, field = karst.shape@data[,1], background = 0)

#reclassify into binary karst(0)/non-karst(100)
m <- c(1, 5, 100)
rclmat <- matrix(m, ncol = 3, byrow = T)
karst.rec <- reclassify(karst.rast, rcl = rclmat, right = NA)


#Convert to proportion of karst (not used for final work)
#karst.agg <- aggregate(karst.rec, fact = 2, fun = mean, expand = F)
#origin(karst.agg) <- origin(rast)
#extent(karst.agg) <- extent(rast)
#writeRaster(karst.agg, filename = "D:\\PhD\\ModellingViviparity\\FullDist_Sal_Lsal\\egv_prep\\TerraClimate\\whymap_karst5km.tif")





