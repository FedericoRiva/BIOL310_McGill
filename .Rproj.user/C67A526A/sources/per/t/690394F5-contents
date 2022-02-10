####
library(tinytex)
library(rmarkdown)

# analysis of diversity
library(vegan)

# to handle raster and vector spatial data
library(raster)
library(sp)

# Open and clean data
library(data.table)

# to get the world map when plotting original data
library(maptools)

# data wrangling
library(dplyr)


#query_gbif = data.table::fread("Data\\raw_gbif.csv", header = TRUE)
query_gbif = data.table::fread("Data\\query_gbif_2013.csv", header = TRUE)


# #query_gbif <- subset(query_gbif, coordinateUncertaintyInMeters < 5000) # & individualCount != 0) potentially remove individual counts
# query_gbif$taxonRank <- as.factor(query_gbif$taxonRank)
# query_gbif <- subset(query_gbif, taxonRank != "GENUS") # & individualCount != 0) potentially remove in
# query_gbif <- subset(query_gbif, identifiedBy == "eButterfly users")

# smaller dataset
#query_gbif <- query_gbif[which(year == 2013),]


# # # was the data download ok? let's check if the cabbage white, a generalist species common in Europe, NA and Oceania, 
# # # occurs only in CA
example <- query_gbif[query_gbif$species == "Pieris rapae"]
data(wrld_simpl)
plot(wrld_simpl, xlim=c(min(example$decimalLongitude)-1,max(example$decimalLatitude)+1), ylim=c(min(example$decimalLongitude)-1,max(example$decimalLatitude)+1), axes=TRUE, col="grey")
points(example$decimalLongitude, example$decimalLatitude, col="blue", pch=20, cex=1.7)
# # # yes, all good. Data restricted to CA, no outliers in the ocean, etc
# note that this plot does not have projected coordinates expressed in metric system, but geographic coordinates expressed in degrees.
# later we will see how to plot datasets in projected systems

####
####
####

# covariates for our models


# euclidead-distance based velocity; https://adaptwest.databasin.org/pages/adaptwest-velocitymed/
VEL_ED <- raster::raster("Data\\velocity\\velocity.ED.tif")

VEL_ED 
# check the properties of this raster; 5 km-resolution raster, projected in WGS84
plot(VEL_ED)
#  up to 25 km/year migration needed to track future climated, simply based on an euclidean distance  (shortest straight line between two cells with analog climates)
# use VLE_ED to set coordinate reference system
crs <- crs(VEL_ED)
crs



# create raster of latitudes
# created a raster having the same extent of VEL_ED, but containing latitude values instead 
LAT <- VEL_ED
raster_coordinates <- coordinates(VEL_ED)
LAT[] <- raster_coordinates[, 2]

# remove pixels around NA for LAT
# create a mask
mask <- VEL_ED
mask[!is.na(mask[])] <- 1
LAT <- LAT*mask
plot(LAT)

# aggregate to assess spatial scales
VEL_ED_50  <- raster::aggregate(VEL_ED, 10)
VEL_ED_500 <- raster::aggregate(VEL_ED_50, 10)

LAT_50  <- raster::aggregate(LAT, 10)
LAT_500 <- raster::aggregate(LAT_50, 10)

# create unique ID cells for every 5, 50 and 500 km cell
ID <- VEL_ED
ID_50 <- VEL_ED_50
ID_500 <- VEL_ED_500

ID[] <- seq(1, ncell(VEL_ED))
ID_50[] <- seq(1, ncell(VEL_ED_50))
ID_500[] <- seq(1, ncell(VEL_ED_500))

# crop NA cells
ID <- ID*mask
ID_50 <- ID_50* (raster::aggregate(mask,10))
ID_500 <- ID_500*(raster::aggregate(mask,100))

# check visually the data
par(mfrow=c(3,3))
plot(VEL_ED)
plot(VEL_ED_50)
plot(VEL_ED_500)
plot(LAT)
plot(LAT_50)
plot(LAT_500)
plot(ID)
plot(ID_50)
plot(ID_500)

## question: how does aggregation affects the covariates?

dev.off() # clean setting for 9 plots
# covariates ready!

# prepare response variable
# project coordinate points



#####
#####
#####

# preparing response data
# always loo kat the data
# project geographic coordinates (from degrees to metric units)

coordinates(query_gbif) <- c("decimalLongitude", "decimalLatitude")
proj4string(query_gbif) <- CRS("+init=epsg:4326") # GBIF geographic coordinates "+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"
query_gbif <- sp::spTransform(query_gbif, crs) # project according to adaptwest raster data
query_gbif <- as.data.frame(query_gbif)

# replace name of columns after projection
names(query_gbif)[names(query_gbif) == 'decimalLongitude'] <- 'x' # changed from degrees to metric
names(query_gbif)[names(query_gbif) == 'decimalLatitude'] <- 'y' # changed from degrees to metric

# prepare coordinates and data for spatial dataframe
coordinates_query_gbif <- query_gbif[ , c("x", "y")]   # coordinates; equivalent to query_gbif[, 49:50]
data_query_gbif   <- query_gbif[ , c(9:10, 31)]   # data you want to keep

# make a SpatialPointsDataFrame object
spdf <- sp::SpatialPointsDataFrame(coords       = coordinates_query_gbif,
                                    data        = data_query_gbif, 
                                    proj4string = crs)


# note that point() does not overlay points properly
# double-check with Laura and Janaina
plot(VEL_ED)
plot(spdf,
     #col = spdf$year,
     add = TRUE)



plot(VEL_ED)
points(query_gbif$x, query_gbif$y)

# same crs; weird
crs(VEL_ED)
crs(spdf)

# super annoying R glitch; does not affect point extraction; see 
# https://stackoverflow.com/questions/28371270/resizing-plot-output-causes-raster-and-points-to-become-misaligned


# look at how many cells have multiple observations
checklists <- raster::rasterize(coordinates_query_gbif, ID, fun='count')

plot(checklists) # difficult to see, 5 km cells are very small at the continental scale
par(mfrow=c(1,3))
plot(log10(checklists))
plot(log10(aggregate(checklists,10))) # 50 km resolution
plot(log10(aggregate(checklists,100))) # 500 km resolution
# clear scaling effects on sampling 

dev.off()


# create an uniqueID per checklist
query_gbif$unique_ID <- paste(query_gbif$eventDate,
                              (paste(query_gbif$x, query_gbif$y))) # a single column combining data and xy coordinates




####
####
####

# inspecting the species occurrence distribution
SOD <- as.data.frame(table(query_gbif$species))
head(SOD)
# as espect the cabbage white (Pieris rapae) is the most common species in these samples (check your backyard if you have never seen one)
# monarch butterflies (Danaus plexippus count for 12.000 observations, but other thant that this is a typical SAD. species abundance distribution)
hist(log10(SOD$Freq))
# most species between 10 and 1000 occurrence data

# inspecting the species occurrence distribution at lower spatial levels (Provinces)
## which provinces we have in the dataframe?
levels(as.factor(query_gbif$stateProvince))

# provinces
query_gbif_subset_SOD_AB <- query_gbif[which(query_gbif$stateProvince == "Alberta"), ]
query_gbif_subset_SOD_QC <- query_gbif[which(query_gbif$stateProvince == "Quebec"), ]
query_gbif_subset_SOD_ON <- query_gbif[which(query_gbif$stateProvince == "Ontario"), ]

SOD_AB <- as.data.frame(table(query_gbif_subset_SOD_AB$species))
SOD_QC <- as.data.frame(table(query_gbif_subset_SOD_QC$species))
SOD_ON <- as.data.frame(table(query_gbif_subset_SOD_ON$species))



##
SOD$Freq <- SOD$Freq/sum(SOD$Freq) # repeat for all other SODs
SOD_AB$Freq <- SOD_AB$Freq/sum(SOD_AB$Freq) # repeat for all other SODs
SOD_QC$Freq <- SOD_QC$Freq/sum(SOD_QC$Freq) # repeat for all other SODs
SOD_ON$Freq <- SOD_ON$Freq/sum(SOD_ON$Freq) # repeat for all other SODs


#par(mfrow=c(1,5))
plot(SOD$Freq[order(SOD$Freq)], xaxt='n', frame.plot=TRUE, 
     main="Species-occurrence ranking plot",
     ylab="Frequency of occurrence", xlab="Species ranked by frequency of occurrence", pch = 16,
     ylim = c(0, 0.10))

par(new=TRUE)
plot(SOD_AB$Freq[order(SOD_AB$Freq)], col = "red", axes=FALSE, frame.plot=TRUE,
     ylab="", xlab="", pch = 16,
     ylim = c(0, 0.10))

par(new=TRUE)
plot(SOD_QC$Freq[order(SOD_QC$Freq)], col = "blue", axes=FALSE, frame.plot=TRUE, 
     ylab="", xlab="", pch = 16,
     ylim = c(0, 0.10))

par(new=TRUE)
plot(SOD_ON$Freq[order(SOD_ON$Freq)], col = "gold", axes=FALSE, frame.plot=TRUE, 
     ylab="", xlab="", pch = 16,
     ylim = c(0, 0.10))

legend(1, 0.09, # position of the legend; the x axis goes from 0 to 0.1 (10% 0f the occurrences observed for a given species)
       legend=c("CA", "AB", "QC", "ON"),
       col=c("black", "red", "blue", "gold"), lty=1, cex=1)

## inference: 
# when considering larger spatial scales, there is a higher dominance on the SOD of very common species (black line always under colored lines)


# is it area? is it # of observations?
SOD <- query_gbif[sample(nrow(query_gbif), 1700), ]
query_gbif_subset_SOD_AB <- query_gbif_subset_SOD_AB[sample(nrow(query_gbif_subset_SOD_AB), 1700), ]
query_gbif_subset_SOD_QC <- query_gbif_subset_SOD_QC[sample(nrow(query_gbif_subset_SOD_QC), 1700), ]
query_gbif_subset_SOD_ON <- query_gbif_subset_SOD_ON[sample(nrow(query_gbif_subset_SOD_ON), 1700), ]

SOD <- as.data.frame(table(query_gbif$species))
SOD_AB <- as.data.frame(table(query_gbif_subset_SOD_AB$species))
SOD_QC <- as.data.frame(table(query_gbif_subset_SOD_QC$species))
SOD_ON <- as.data.frame(table(query_gbif_subset_SOD_ON$species))



##
SOD$Freq <- SOD$Freq/sum(SOD$Freq) # repeat for all other SODs
SOD_AB$Freq <- SOD_AB$Freq/sum(SOD_AB$Freq) # repeat for all other SODs
SOD_QC$Freq <- SOD_QC$Freq/sum(SOD_QC$Freq) # repeat for all other SODs
SOD_ON$Freq <- SOD_ON$Freq/sum(SOD_ON$Freq) # repeat for all other SODs


#par(mfrow=c(1,5))
plot(SOD$Freq[order(SOD$Freq)], xaxt='n', frame.plot=TRUE, 
     main="Species-occurrence ranking plot",
     ylab="Frequency of occurrence", xlab="Species ranked by frequency of occurrence", pch = 16,
     ylim = c(0, 0.10))

par(new=TRUE)
plot(SOD_AB$Freq[order(SOD_AB$Freq)], col = "red", axes=FALSE, frame.plot=TRUE,
     ylab="", xlab="", pch = 16,
     ylim = c(0, 0.10))

par(new=TRUE)
plot(SOD_QC$Freq[order(SOD_QC$Freq)], col = "blue", axes=FALSE, frame.plot=TRUE, 
     ylab="", xlab="", pch = 16,
     ylim = c(0, 0.10))

par(new=TRUE)
plot(SOD_ON$Freq[order(SOD_ON$Freq)], col = "gold", axes=FALSE, frame.plot=TRUE, 
     ylab="", xlab="", pch = 16,
     ylim = c(0, 0.10))

legend(1, 0.09, # position of the legend; the x axis goes from 0 to 0.1 (10% 0f the occurrences observed for a given species)
       legend=c("CA", "AB", "QC", "ON"),
       col=c("black", "red", "blue", "gold"), lty=1, cex=1)

## support for area effect on the species frequency occurrence rank plot

####
####
####


# estract covariates
coordinates_query_gbif$extract_ID <- extract(ID, coordinates_query_gbif[, 1:2])
coordinates_query_gbif$extract_ID_50 <- extract(ID_50, coordinates_query_gbif[, 1:2])
coordinates_query_gbif$extract_ID_500 <- extract(ID_500, coordinates_query_gbif[, 1:2])

coordinates_query_gbif$extract_LAT <- extract(LAT, coordinates_query_gbif[, 1:2])
coordinates_query_gbif$extract_LAT_50 <- extract(LAT_50, coordinates_query_gbif[, 1:2])
coordinates_query_gbif$extract_LAT_500 <- extract(LAT_500, coordinates_query_gbif[, 1:2])

coordinates_query_gbif$extract_VEL_ED <- extract(VEL_ED, coordinates_query_gbif[, 1:2])
coordinates_query_gbif$extract_VEL_ED_50 <- extract(VEL_ED_50, coordinates_query_gbif[, 1:2])
coordinates_query_gbif$extract_VEL_ED_500 <- extract(VEL_ED_500, coordinates_query_gbif[, 1:2])

# # check relationship between LAT and VEL_ED across scales
# # might take a while; select the next three rows and press CTRL + SHIFT + C to run the code
# par(mfrow=c(1,3))
# plot(coordinates_query_gbif$extract_VEL_ED ~ coordinates_query_gbif$extract_LAT)
# plot(coordinates_query_gbif$extract_VEL_ED_50 ~ coordinates_query_gbif$extract_LAT_50)
# plot(coordinates_query_gbif$extract_VEL_ED_500 ~ coordinates_query_gbif$extract_LAT_500)

# create a unique coord field to merge with original data
coordinates_query_gbif$pasted_coord <- paste(coordinates_query_gbif$x, coordinates_query_gbif$y)
query_gbif$pasted_coord <- paste(query_gbif$x, query_gbif$y)

# select instead of a Province, one of the pasted_coord in your data. Check how many entries there are in different cells, e.g., with table()


# 
# ###
# ### data preparation to be removed
# ###
# 
# 
# # clean data
# # split depending on the year
# list_checklists <- split(query_gbif, query_gbif$unique_ID)
# list_checklists <- Filter(function(x) nrow(x) > 1, list_checklists) # we assume that visits with only one observations were not checklists
# 
# ## create parallel list to get the unique ID associated with every row of the final table
# list_checklists2 <- list()
# for(i in 1:length(list_checklists)){
#   list_checklists2[[i]] <- unique(list_checklists[[i]][,c(10,51)])
# }
# list_checklists2 <- Filter(function(x) nrow(x) > 1, list_checklists2) # we assume that visits with only one observations were not checklists
# 
# #keep only uniqueID
# for(i in 1:length(list_checklists2)){
#   list_checklists2[[i]] <- unique(list_checklists2[[i]][,2])
# }
# 
# ##
# 
# list_checklists_species <- list()
# for(i in 1: length(list_checklists)){
#   list_checklists_species[[i]] <- list_checklists[[i]][,10]
# }
# 
# #some lists are repeated
# for(i in 1: length(list_checklists)){
#   list_checklists_species[[i]] <- unique(list_checklists_species[[i]])
# }
# # remove lists that had two observations of the same species
# list_checklists_species <- Filter(function(x) length(x) > 1, list_checklists_species) # we assume that visits with only one observations were not checklists
# 
# 
# 
# list_species <- SOD$Var1
# table_species <- matrix(, ncol = length(list_species), nrow = 1)
# table_species <- as.data.frame(t(rep(0, length(list_species))))
# colnames(table_species) <- list_species
# 
# 
# list_tables <- list()
# for(i in 1:length(list_checklists_species)){
#    list_tables[[i]] <- table_species
# }
# 
# # transform in a table
# for(i in 1:length(list_checklists_species)){
#   list_tables[[i]] <- table_species[,list_checklists_species[[i]]]
# 
# }
# 
# # convert 0s to 1s
# for(i in 1:length(list_checklists_species)){
#   list_tables[[i]][list_tables[[i]] == 0 ] <- 1
# }
# 
# #list_tables[[13]]
# 
# # merge with species table
# for(i in 1:length(list_checklists_species)){
#   list_tables[[i]] <- cbind(list_tables[[i]], table_species[,(setdiff(colnames(table_species), colnames(list_tables[[i]])))])
# }
# 
# #list_tables[[13]]
# 
# #prova33 <- cbind(list_tables[[27]], table_species[,(setdiff(colnames(table_species), colnames(list_tables[[27]])))])
# 
# 
# # order columns
# for(i in 1:length(list_checklists_species)){
#   list_tables[[i]][,order(colnames(list_tables[[i]]))]
# }
# 
# #
# # list_test <- list()
# # for(i in 1:length(list_checklists_species)){
# #   list_test[[i]] <- ncol(list_checklists_species[[i]])
# # }
# #
# # test_table <- do.call(rbind.data.frame, list_test)
# 
# 
# # convert into a table
# checklist_table <- do.call(rbind.data.frame, list_tables)
# id_table <- do.call(rbind.data.frame, list_checklists2)
# colnames(id_table)[1] <- "unique_ID"
# 
# checklist_table$unique_ID <- id_table$unique_ID
# 
# info_needed_for_table <- unique(query_gbif[,c(29,30,31,49,50,51)])
# 
# checklist_table <- left_join(checklist_table, info_needed_for_table)
# checklist_table$pasted_coord <- paste(checklist_table$x, checklist_table$y)
# 
# checklist_table$year <- as.factor(checklist_table$year)
# 
# # relevel to have biennual data (temporal scale #2)
# checklist_table$year2 <- as.factor(checklist_table$year)
# levels(checklist_table$year2)[levels(checklist_table$year2)=="2010"] <-"2010-2011"
# levels(checklist_table$year2)[levels(checklist_table$year2)=="2011"] <-"2010-2011"
# levels(checklist_table$year2)[levels(checklist_table$year2)=="2012"] <-"2012-2013"
# levels(checklist_table$year2)[levels(checklist_table$year2)=="2013"] <-"2012-2013"
# levels(checklist_table$year2)[levels(checklist_table$year2)=="2014"] <-"2014-2015"
# levels(checklist_table$year2)[levels(checklist_table$year2)=="2015"] <-"2014-2015"
# levels(checklist_table$year2)[levels(checklist_table$year2)=="2016"] <-"2016-2017"
# levels(checklist_table$year2)[levels(checklist_table$year2)=="2017"] <-"2016-2017"
# levels(checklist_table$year2)[levels(checklist_table$year2)=="2018"] <-"2018-2019"
# levels(checklist_table$year2)[levels(checklist_table$year2)=="2019"] <-"2018-2019"
# levels(checklist_table$year2)[levels(checklist_table$year2)=="2020"] <-"2020-2021"
# levels(checklist_table$year2)[levels(checklist_table$year2)=="2021"] <-"2020-2021"
# 
# # relevel to have triennial data (temporal scale #3)
# checklist_table$year3 <- as.factor(checklist_table$year)
# levels(checklist_table$year3)[levels(checklist_table$year3)=="2010"] <-"2010-2012"
# levels(checklist_table$year3)[levels(checklist_table$year3)=="2011"] <-"2010-2012"
# levels(checklist_table$year3)[levels(checklist_table$year3)=="2012"] <-"2010-2012"
# levels(checklist_table$year3)[levels(checklist_table$year3)=="2013"] <-"2013-2015"
# levels(checklist_table$year3)[levels(checklist_table$year3)=="2014"] <-"2013-2015"
# levels(checklist_table$year3)[levels(checklist_table$year3)=="2015"] <-"2013-2015"
# levels(checklist_table$year3)[levels(checklist_table$year3)=="2016"] <-"2016-2018"
# levels(checklist_table$year3)[levels(checklist_table$year3)=="2017"] <-"2016-2018"
# levels(checklist_table$year3)[levels(checklist_table$year3)=="2018"] <-"2016-2018"
# levels(checklist_table$year3)[levels(checklist_table$year3)=="2019"] <-"2019-2021"
# levels(checklist_table$year3)[levels(checklist_table$year3)=="2020"] <-"2019-2021"
# levels(checklist_table$year3)[levels(checklist_table$year3)=="2021"] <-"2019-2021"
# 
# 
# ### merge covariates and species table
# sites <- unique(checklist_table$pasted_coord)
# subset_queried_coord <- coordinates_query_gbif[coordinates_query_gbif$pasted_coord %in% sites,]
# subset_queried_coord <- unique(subset_queried_coord)

final_table <- left_join(checklist_table, subset_queried_coord, by = "pasted_coord")
# this just for the class



#write.csv(final_table, "final_table.csv")
#final_table = data.table::fread("Data\\final_table.csv", header = TRUE)


###
###
###
### Analysis at 5 km

## split into the 5-km cell communities
split_5 <- split(final_table, final_table$extract_ID)
# 4674 sites with at least one checklist
# we keep only sites with at least 10 checklists
split_5 <- Filter(function(x) nrow(x) > 10, split_5) 
# 465 sites

# how many transects per year
## split each cell in the 12 years
split_5 <- lapply(split_5, function(x) split(x, x$year))
for(i in 1:length(split_5)){
  split_5[[i]] <- Filter(function(x) nrow(x) > 5, split_5[[i]])
}

# remove sites with less than 2 years
split_5 <- Filter(function(x) length(x) > 2, split_5) # we assume that visits with only one observations were not checklists 


## randomly select in each site an equal number of checklists
for(i in 1:length(split_5)){
  for(j in 1:length(split_5[[i]])){
    split_5[[i]][[j]] <- split_5[[i]][[j]][sample(nrow(split_5[[i]][[j]]), 5, replace = FALSE), ]
  }
}

# remove sites with only one or two years
split_5 <- Filter(function(x) length(x) > 2, split_5) # we assume that visits with only one observations were not checklists 



# data ready


## species richness
split_5_richness <- split_5
for(i in 1:length(split_5_richness)){
  for(j in 1:length(split_5_richness[[i]])){
    split_5_richness[[i]][[j]] <- sum(colSums(split_5[[i]][[j]][,1:262]) != 0)
  }
}

for(i in 1:length(split_5_richness)){
  split_5_richness[[i]] <- do.call(rbind.data.frame, split_5_richness[[i]])
}

# community table
split_5_comm <- split_5
for(i in 1:length(split_5_comm)){
  for(j in 1:length(split_5_comm[[i]])){
    split_5_comm[[i]][[j]] <- as.data.frame(t(colSums(split_5_comm[[i]][[j]][,1:262])))
    split_5_comm[[i]][[j]][split_5_comm[[i]][[j]] > 0] <- 1
  }
}

for(i in 1:length(split_5_comm)){
  split_5_comm[[i]] <- do.call(rbind.data.frame, split_5_comm[[i]])
}


split_5_similarity <- split_5_comm
for(i in 1:length(split_5_similarity)){
  split_5_similarity[[i]] <- mean(vegdist(split_5_similarity[[i]], method = "jaccard"))
}

## years
years <- list()
for(i in 1: length(split_5)){
years[[i]] <-  as.numeric(names(split_5[[i]])) 
}




# trends in richness
models <- list()
for(i in 1:length(years)){
  models[[i]] <- lm(split_5_richness[[i]][,1] ~ years[[i]])
  models[[i]] <- models[[i]]$coefficients[2] # keep only the year effect
}

# for(i in 1:length(years)){
#   models[[i]] <- models[[i]]$coefficients[2]
# }

# velocity
split_5_velocity <- split_5

for(i in 1:length(split_5_velocity)){
    split_5_velocity[[i]] <- split_5_velocity[[i]][[1]]$extract_VEL_ED[[1]]
  }

split_5_latitude <- split_5

for(i in 1:length(split_5_latitude)){
  split_5_latitude[[i]] <- split_5_latitude[[i]][[1]]$extract_LAT[[1]]
}





models <- do.call(rbind.data.frame, models)
colnames(models)[1] <- "slope"
similarity <-do.call(rbind.data.frame, split_5_similarity)
colnames(similarity)[1] <- "similarity"
velocity <- do.call(rbind.data.frame, split_5_velocity)
colnames(velocity)[1] <- "velocity"
latitude <- do.call(rbind.data.frame, split_5_latitude)
colnames(latitude)[1] <- "latitude"

summary(lm(models$slope ~ latitude$latitude))


#Create a function to generate a continuous color palette
rbPal <- colorRampPalette(c('blue','red'))

#This adds a column of color values
# based on the y values
models$Col <- rbPal(10)[as.numeric(cut(models$slope,breaks = 10))]
plot(models$slope ~ latitude$latitude,  
     col = models$Col, 
     pch = 16,
     main="Relationship between biodiversity trends and latitude", 
     xlab="Latitude (m)", ylab="slope of the richness ~ year model")



plot(models$slope ~ velocity$velocity)

# check other relationsips
plot(similarity$similarity ~ latitude$latitude)
plot(models$slope ~ similarity$similarity)



