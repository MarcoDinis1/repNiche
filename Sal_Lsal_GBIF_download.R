
library(maptools) #geo data handling
#install.packages("rgbif") 
#install.packages("ridigbio")
#install.packages("rvertnet")
library("rgbif") #gbif wrap
library("RCurl") #ftp download wrap
library(dismo) # sdm package

#citation("rgbif")
#citation("dismo")
#citation("maptools")
#citation("RCurl")

#download path
dl.path <- "C:\\PhD\\ModellingViviparity\\FullDist_Sal_Lsal\\presence_prep\\GBIF"
setwd(dl.path)


### Download and plot the data ###

my.sal <- gbif("Salamandra salamandra")

data(wrld_simpl)
plot(wrld_simpl, col = "light yellow", axes = T)
points(my.sal$lon, my.sal$lat, col = "red", cex = 0.5)
#text(-140, -50, "Lyciasalamandra/nspp")


######## checking data quality ####
names(my.sal)  # check all columns for relevant info
sort(my.sal) #same as above but ordered alphabetically
# examples of potentially relevant ones to check:
sort(unique(my.sal$coordinatePrecision))
sort(unique(my.sal$coordinateUncertaintyInMeters))
sort(unique(my.sal$dataGeneralizations))
sort(unique(my.sal$informationWithheld))
sort(unique(my.sal$occurrenceStatus))
sort(unique(my.sal$year))
sort(unique(my.sal$country))
sort(unique(my.sal$species))
sort(unique(my.sal$fieldNotes))
sort(unique(my.sal$georeferenceRemarks))
sort(unique(my.sal2$collectionCode))
sort(unique(my.sal$basisOfRecord))

##Tools for cleaning dataset##
##DON'T RUN THE WHOLE THING!!! Instead use each section as needed and adjust it for the data you want to filter out
#clear records with no coordinates
my.sal <- subset(my.sal, is.na(lat) == FALSE) 
#clear uncertain coordinates
my.sal <- subset(my.sal, coordinateUncertaintyInMeters < 1000 | is.na(coordinateUncertaintyInMeters)) 
arr <- sort(unique(my.sal$dataGeneralizations))
ind <- c(2, 3) #put here stuff you want to exclude
for (i in ind)
{
my.sal <- subset(my.sal, dataGeneralizations != arr[i] | is.na(dataGeneralizations))
}
#clear absence data
my.sal <- subset(my.sal, occurrenceStatus != "absent" | is.na(occurrenceStatus)) 
#clear by country
arr<- sort(unique(my.sal$country))
ind <- c(11, 13, 15, 18, 19, 21, 31) #put here stuff you want to exclude
for (i in ind)
{
my.sal <- subset(my.sal, country != arr[i] | is.na(country)) 
}
#Clear data which have been modified for whatever reason.  (remember to check the actual values in the array to determine if you want to remove them)
arr <- sort(unique(my.sal$informationWithheld))
ind <- c()
for (i in ind) 
{
  my.sal <- subset(my.sal, informationWithheld != arr[i] | is.na(informationWithheld))
}
#clear by georeference remarks (rarely necessary)
arr <- sort(unique(my.sal$georeferenceRemarks))
ind <- c()
my.sal <- subset(my.sal, georeferenceRemarks != arr[i] | is.na(georeferenceRemarks)) 
#clear records collected prior to 1950
my.sal <- subset(my.sal, year >= 1950 | is.na(year))
#clear bad datasets
arr<- sort(unique(my.sal$collectionCode))
ind <- c("amp","Artenfinder Rheinland-Pfalz","BoBO","BOLD","croatia-2011","EDP_SABOR_Herpetology_Baseline","EDP_TUA_Amphibians_PME","EDP_TUA_Quiropteros_Ab_CS","galicia-2012","HERP","KUH","MB04","MCNB-Cord","MNCN_Herpeto","Montseny","naturgucker","NCSM-Herp","Observations","portugal-2003","portugal-2009","romania-2008","slovakia-2005","slovenia-2004","SMNS-Z-herpcoll","user: 1236","user: 1407","user: 1563","user: 1694","user: 2650","user: 3500","user: 4965","user: 7015","user: 7521","VERTEBRATES")
#put here stuff you want to exclude
for (i in ind)
{
  my.sal2 <- subset(my.sal, collectionCode == ind[1] | collectionCode == ind[2] |collectionCode == ind[3] |collectionCode == ind[4] |collectionCode == ind[5] |collectionCode == ind[6] |collectionCode == ind[7] |collectionCode == ind[8] |collectionCode == ind[9] |collectionCode == ind[10] |collectionCode == ind[11] |collectionCode == ind[12] |collectionCode == ind[13] |collectionCode == ind[14] |collectionCode == ind[15] |collectionCode == ind[16] |collectionCode == ind[17] |collectionCode == ind[18] |collectionCode == ind[19] |collectionCode == ind[20] |collectionCode == ind[21] |collectionCode == ind[22] |collectionCode == ind[23] |collectionCode == ind[24] |collectionCode == ind[25] |collectionCode == ind[26] |collectionCode == ind[27] |collectionCode == ind[28] |collectionCode == ind[29] |collectionCode == ind[30] |collectionCode == ind[31] |collectionCode == ind[32] |collectionCode == ind[33] |collectionCode == ind[34] |is.na(collectionCode)) #did this one by hand (adding manually the collections I want to keep). To be improved
}

###write to file###
write.csv2(my.sal,paste(dl.path,"\\Salsal_20190718_inf.csv",sep=""))
