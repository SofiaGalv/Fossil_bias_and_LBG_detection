# ------------------------------------------------------------------------------
# Publication: Biases can't hide conspicuous biodiversity patterns 
# Last updated: 2024-07-09
# Author: Sofia Galvan
# Email: sofia.galvan@uvigo.es
# Repository: 
# ------------------------------------------------------------------------------

library(devtools)
library(dplyr)
library(ggplot2)
library(paletteer)
library(quantreg)
library(raster)
library(rgdal)
library(rgeos)
library(rnaturalearth)
library(sf)
library(sp)
library(spData)
library(terra)
library(tidyverse)

############################## Load source data ################################

### 1. Load WorldClim climate
setwd("C:/Users/usuario/Desktop/PhD/Paper3 (current_fossils)/GUM_data/Climatic_niche/wc2.1_5m_tavg")
files_dir <- list.files()
tavg_stack <- stack()
for (i in 1:length(files_dir)) {
  tavg_raster <- raster(files_dir[i])
  tavg_stack <- stack(tavg_stack, tavg_raster)
}
x <- raster(resolution = 0.5, ext = extent(-180, 180, -90, 90))
tavg_stack_05 <- resample(tavg_stack, x)
tavg_array <- as.array(tavg_stack_05)
#saveRDS(tavg_array, "tavg_array.rds")

setwd("C:/Users/usuario/Desktop/PhD/Paper3 (current_fossils)/GUM_data/Climatic_niche/wc2.1_5m_prec")
files_dir <- list.files()
files_dir <- files_dir[-1]
prec_stack <- stack()
for (i in 1:length(files_dir)) {
  prec_raster <- raster(files_dir[i])
  prec_stack <- stack(prec_stack, prec_raster)
}
prec_stack_05 <- resample(prec_stack, x)
prec_array <- as.array(prec_stack_05)
#saveRDS(prec_array, "prec_array.rds")

devtools::source_url("https://github.com/MAPASlab/KoppenGeiger_inR/blob/main/kgreclass_Rfunction.R?raw=TRUE")
matrix_reclassified <- kg_reclass(tavg_array, prec_array, type = "broadclass")
wc_raster <- terra::rast(x = matrix_reclassified, 
                         extent = ext(-180, 180, -90, 90))
crs(wc_raster) <- "+proj=longlat +datum=WGS84 +no_defs"
wc_column <- terra::as.data.frame(wc_raster, xy = T, na.rm = F)

wc_moll <- terra::project(wc_raster, 
                          y = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84",
                          method = "near")

land <- ne_download(scale = 10, type = "land", category = "physical")
land_moll <- st_transform(land, "+proj=moll")
plot(land_moll$geometry, border = "grey80", col = NA, lwd=0.08)
raster::plot(wc_moll, 
             col=c("#f9d14a", "#ab3329", "#ed968c", "#7c4b73", "#88a0dc"),
             type = "classes",
             axes = F, 
             levels = c("Tropical", "Arid", "Temperate", "Cold", "Polar"),
             add = T)


### 2. Load Unconsolidated sediments (UcS) file.
setwd("C:/Users/usuario/Desktop/PhD/Paper3 (current_fossils)/GUM_data/GUM_v1.0")
GUM_raster <- rast("gum_v1.0_0point5deg.txt.asc")

compareGeom(GUM_raster, as(x, "SpatRaster"), extent = T)
GUM_raster <- terra::resample(GUM_raster, as(x, "SpatRaster"), method = "near") 
plot(GUM_raster)
matr <- rbind(c(9, NA), c(18, NA), c(19, NA), c(20, NA))
GUM_new <- classify(GUM_raster, matr, include.lowest = T)
matr2 <- matrix(c(1,39,1), ncol = 3, byrow = T)
GUM_new <- classify(GUM_new, matr2, include.lowest = T)
GUM_col <- terra::as.data.frame(GUM_new, xy = T, na.rm = F)
GUM_col[is.nan(GUM_col$gum_v1.0_0point5deg.txt), 3] <- NA
GUM_moll <- terra::project(GUM_new, 
                           y = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84",
                           method = "near")


### 3. Load IUCN mammal data and COMBINE dataset
setwd("C:/Users/usuario/Desktop/PhD/Paper1 (tropics)/Koppen-Geiger biomes/Biomes_new/MAMMALS_TERRESTRIAL_ONLY/Mam_new")
mam_df <- readRDS("mam_final.rds") #Terrestrial mammals
setwd("C:/Users/usuario/Desktop/PhD/Paper3 (current_fossils)/GUM_data/MAMMALS_FRESHWATER")
mam2_df <- readRDS("df_fresman.rds") #Freshwater mammals

mam_df <- cbind(mam_df[,1:5617], mam2_df[,3:126], mam_df[,5618:5619]) #5743 obs
mam_df[sapply(mam_df, is.nan)] <- 0 
mam_names <- colnames(mam_df)[-c(1,2, 5742, 5743)] 

setwd("C:/Users/usuario/Desktop/PhD/Paper3 (current_fossils)/GUM_data/Bergmann")
Combi_data <- read.csv("trait_data_imputed.csv", sep = ",") 
to.drop <- which(duplicated(Combi_data$iucn2020_binomial)) #301 duplicated
Combi_data <- Combi_data[-to.drop,] 

# Keep species from both IUCN and COMBINE
mam_2 <- mam_names[mam_names %in% Combi_data$iucn2020_binomial] 
spp_IUCNtodrop <- setdiff(mam_names, mam_2) #38 IUCN spp are not in COMBINE
mam_df_2 <- mam_df[, !names(mam_df) %in% spp_IUCNtodrop] #5701/5739 = 99.3% from IUCN
Combi_2 <- Combi_data[Combi_data$iucn2020_binomial %in% mam_2,] #261 COMBINE spp are not in IUCN
identical(sort(mam_2), Combi_2$iucn2020_binomial) #TRUE



############################## Join all the data ###############################

mam_df_2[,5706] <- GUM_col[,3] #Add UcS column
colnames(mam_df_2)[colnames(mam_df_2) == "V5706"] <- "GUM_USC" 
mam_df_2[,5707] <- wc_column[,3] #Add WordClim column
colnames(mam_df_2)[colnames(mam_df_2) == "V5707"] <- "WC_zone" 


richness <- rowSums(mam_df_2[,3:5703]) #Add richness richness column
mam_df_2[,5708] <- richness
colnames(mam_df_2)[colnames(mam_df_2) == "V5708"] <- "Richness"
mam_df_2[mam_df_2$Richness == 0, 5708] <- NA


matr <- matrix(c(1,5,1), ncol = 3, byrow = T) #Add land surface column
Land <- classify(wc_raster, matr, include.lowest = T)
Land_col <- terra::as.data.frame(Land, xy = T, na.rm = F)
Land_col[is.nan(Land_col$lyr.1), 3] <- NA
mam_df_2[,5709] <- Land_col[,3]
colnames(mam_df_2)[colnames(mam_df_2) == "V5709"] <- "Land"


area_raster <- terra::cellSize(wc_raster, mask = T, unit = "km") # Add area column
area_column <- terra::as.data.frame(area_raster, xy = T, na.rm = F)
area_column[is.nan(area_column$area), 3] <- NA
identical(sort(colnames(mam_df_2[,3:5703])), sort(Combi_2$iucn2020_binomial)) #T
mam_df_2[,5710] <- area_column[,3] 
colnames(mam_df_2)[5710] <- "Area" 

#For calculating spp range size
mam_forange <- mam_df_2[mam_df_2$Land == 1 & !is.na(mam_df_2$Land),]
mam_forange <- mam_forange[, is.na(colSums(mam_forange)) | colSums(mam_forange !=0, na.rm = T) > 0]
setdiff(colnames(mam_df_2), colnames(mam_forange)) #2 spp lost Land mask

for (i in 3:5701) { 
  print(i)
  mam_forange[,i] <- mam_forange[,i]*mam_forange[,5708]
}
range_size <- colSums(mam_forange[,3:5701], na.rm = T)
rangesize_df <- data.frame(names(range_size), range_size, row.names =  NULL)
colnames(rangesize_df) <- c("iucn2020_binomial", "range_size")
Mam_r_names <- colnames(mam_forange[,3:5701])
Combi_3 <- Combi_2[Combi_2$iucn2020_binomial %in% Mam_r_names,] #5699 rows
identical(sort(colnames(mam_forange[,3:5701])), sort(Combi_3$iucn2020_binomial)) #T
Combi_r <- merge(Combi_3, rangesize_df, by = "iucn2020_binomial")
identical(sort(colnames(mam_forange[,3:5701])), sort(Combi_r$iucn2020_binomial))


#Add Fossil record data
pbdb_down <- read.csv("pbdb_data (1).csv",header = TRUE, 
                      skip = 17, stringsAsFactors = FALSE)
pbdb_down[pbdb_down$class == "Eufolivora", "class"] <- "Mammalia"
Mam_pbdb <-pbdb_down[pbdb_down$class == "Mammalia",] 
Mam_pbdb2 <- Mam_pbdb[,c("lng", "lat")] 
Mam_pbdb2$count <- 1
unique_Mam <- unique(Mam_pbdb2[,1:2]) 
table_Mam <- Mam_pbdb %>% group_by(lat, lng) %>% summarise(count = n()) %>% as.data.frame()
identical(sort(unique_Mam$lng), sort(table_Mam$lng)) #T
identical(sort(unique_Mam$lat), sort(table_Mam$lat)) #T
table_Mam <- data.frame(table_Mam[,2], table_Mam[,1], table_Mam[,3])
colnames(table_Mam) <- c("x", "y", "occ")
Mam_raster <- terra::rasterize(x = as.matrix(table_Mam[,1:2]), y = as(x, "SpatRaster"), 
                               values = table_Mam[,3], update = T)
plot(Mam_raster,col = rev(grDevices::topo.colors(50)))

Mamfos_col <- terra::as.data.frame(Mam_raster, xy = T, na.rm = F)
mam_df_2[,5711] <- Mamfos_col[,3]
colnames(mam_df_2)[colnames(mam_df_2) == "V5711"] <- "Fossils"


########################## Land mask and UcS filtering #########################

mam_df_land <- mam_df_2[mam_df_2$Land == 1 & !is.na(mam_df_2$Land),]
which(colSums(mam_df_land[,3:5703]) == 0) #2 lost spp = 5699 spp
intersect(colnames(mam_df_2), 
          names(which(colSums(mam_df_land[,3:5703]) == 0))) #Nyctophilus howensis and Pteropus rodricensis

Ric_mam_fil <- mam_df_land[mam_df_land$GUM_USC == 1 & !is.na(mam_df_land$GUM_USC),]
length(which(colSums(Ric_mam_fil[,3:5703]) == 0)) 
Mam_fil_names <- names(which(colSums(Ric_mam_fil[,3:5703]) != 0)) #4745 spp
Combi_4 <- Combi_r[Combi_r$iucn2020_binomial %in% Mam_fil_names,] #4745

mam_original <- mam_df_2[,c("X", "Y", "Richness")] #Original
mam_original_raster <- terra::rasterize(x = as.matrix(mam_original[,1:2]), 
                                        y = as(x, "SpatRaster"), 
                                        values = mam_original[,3], update = T)
crs(mam_original_raster) <- "+proj=longlat +datum=WGS84 +no_defs"

mam_df_fr <- mam_df_land[,c("X", "Y", "Richness")] #Data in land 
mam_raster <- terra::rasterize(x = as.matrix(mam_df_fr[,1:2]), 
                               y = as(x, "SpatRaster"), 
                               values = mam_df_fr[,3], update = T)
crs(mam_raster) <- "+proj=longlat +datum=WGS84 +no_defs"

mam_fil_df <- Ric_mam_fil[,c("X", "Y", "Richness")] #UcS filter
mam_fil_raster <- terra::rasterize(x = as.matrix(mam_fil_df[,1:2]), 
                                   y = as(x, "SpatRaster"), 
                                   values = mam_fil_df[,3], update = T)
crs(mam_fil_raster) <- "+proj=longlat +datum=WGS84 +no_defs"


############################## Range size filter ###############################

identical(sort(Mam_fil_names), sort(Combi_4$iucn2020_binomial))
set.seed(288)
rs25 <- sample(seq_len(nrow(Combi_4)), round(0.25*nrow(Combi_4)), 
               prob = Combi_4$range_size)
Combi_25rs <- Combi_4[rs25,] 
rs50 <- sample(seq_len(nrow(Combi_4)), round(0.50*nrow(Combi_4)), 
               prob = Combi_4$range_size)
Combi_50rs <- Combi_4[rs50,] 
rs75 <- sample(seq_len(nrow(Combi_4)), round(0.75*nrow(Combi_4)), 
               prob = Combi_4$range_size)
Combi_75rs <- Combi_4[rs75,] 

Mam_25rs <- Ric_mam_fil[,c("X", "Y",Combi_25rs$iucn2020_binomial, "Land_mass", 
                           "Clim_zone", "GUM_USC", "WC_zone", "Richness", "Land", 
                           "Area", "Fossils")] 
Mam_50rs <- Ric_mam_fil[,c("X", "Y",Combi_50rs$iucn2020_binomial, "Land_mass", 
                           "Clim_zone", "GUM_USC", "WC_zone", "Richness", "Land", 
                           "Area","Fossils")] 
Mam_75rs <- Ric_mam_fil[,c("X", "Y",Combi_75rs$iucn2020_binomial, "Land_mass", 
                           "Clim_zone", "GUM_USC", "WC_zone", "Richness", "Land", 
                           "Area","Fossils")] 

identical(sort(colnames(Mam_25rs[,3:1188])), sort(Combi_25rs$iucn2020_binomial)) #T
identical(sort(colnames(Mam_50rs[,3:2374])), sort(Combi_50rs$iucn2020_binomial)) #T
identical(sort(colnames(Mam_75rs[,3:3561])), sort(Combi_75rs$iucn2020_binomial)) #T

richness <- apply(Mam_25rs[,3:1188], 1, function(x) length(x[x>0])) 
Mam_25rs[,1197] <- richness
colnames(Mam_25rs)[colnames(Mam_25rs) == "V1197"] <- "Richness2" #New column
Mam_25rs[Mam_25rs$Richness2 == 0, 1197] <- NA

richness <- apply(Mam_50rs[,3:2374], 1, function(x) length(x[x>0]))
Mam_50rs[,2383] <- richness
colnames(Mam_50rs)[colnames(Mam_50rs) == "V2383"] <- "Richness2" #New column
Mam_50rs[Mam_50rs$Richness2 == 0, 2383] <- NA

richness <- apply(Mam_75rs[,3:3561], 1, function(x) length(x[x>0]))
Mam_75rs[,3570] <- richness
colnames(Mam_75rs)[colnames(Mam_75rs) == "V3570"] <- "Richness2" #New column
Mam_75rs[Mam_75rs$Richness2 == 0, 3570] <- NA


Mam_25rs_fr <- Mam_25rs[,c("X", "Y", "Richness2")]
Mam_25rs_raster <- terra::rasterize(x = as.matrix(Mam_25rs_fr[,1:2]), 
                                    y = as(x, "SpatRaster"), 
                                    values = Mam_25rs_fr[,3], update = T)
crs(Mam_25rs_raster) <- "+proj=longlat +datum=WGS84 +no_defs"

Mam_50rs_fr <- Mam_50rs[,c("X", "Y", "Richness2")]
Mam_50rs_raster <- terra::rasterize(x = as.matrix(Mam_50rs_fr[,1:2]), 
                                    y = as(x, "SpatRaster"), 
                                    values = Mam_50rs_fr[,3], update = T)
crs(Mam_50rs_raster) <- "+proj=longlat +datum=WGS84 +no_defs"

Mam_75rs_fr <- Mam_75rs[,c("X", "Y", "Richness2")]
Mam_75rs_raster <- terra::rasterize(x = as.matrix(Mam_75rs_fr[,1:2]), 
                                    y = as(x, "SpatRaster"), 
                                    values = Mam_75rs_fr[,3], update = T)
crs(Mam_75rs_raster) <- "+proj=longlat +datum=WGS84 +no_defs"

#Richness before/after
Ric_mam_fil$WC_zone <- as.factor(Ric_mam_fil$WC_zone)
Mam_75rs$WC_zone <- as.factor(Mam_75rs$WC_zone)
Mam_50rs$WC_zone <- as.factor(Mam_50rs$WC_zone)
Mam_25rs$WC_zone <- as.factor(Mam_25rs$WC_zone)

Mam_f1 <- Ric_mam_fil[Ric_mam_fil[,"WC_zone"] == 1,]
Mam_f2 <- Ric_mam_fil[Ric_mam_fil[,"WC_zone"] == 2,]
Mam_f3 <- Ric_mam_fil[Ric_mam_fil[,"WC_zone"] == 3,]
Mam_f4 <- Ric_mam_fil[Ric_mam_fil[,"WC_zone"] == 4,]
Mam_f5 <- Ric_mam_fil[Ric_mam_fil[,"WC_zone"] == 5,]

Mam_75rs_f1 <- Mam_75rs[Mam_75rs[,"WC_zone"] == 1,]
Mam_75rs_f2 <- Mam_75rs[Mam_75rs[,"WC_zone"] == 2,]
Mam_75rs_f3 <- Mam_75rs[Mam_75rs[,"WC_zone"] == 3,]
Mam_75rs_f4 <- Mam_75rs[Mam_75rs[,"WC_zone"] == 4,]
Mam_75rs_f5 <- Mam_75rs[Mam_75rs[,"WC_zone"] == 5,]

Mam_50rs_f1 <- Mam_50rs[Mam_50rs[,"WC_zone"] == 1,]
Mam_50rs_f2 <- Mam_50rs[Mam_50rs[,"WC_zone"] == 2,]
Mam_50rs_f3 <- Mam_50rs[Mam_50rs[,"WC_zone"] == 3,]
Mam_50rs_f4 <- Mam_50rs[Mam_50rs[,"WC_zone"] == 4,]
Mam_50rs_f5 <- Mam_50rs[Mam_50rs[,"WC_zone"] == 5,]

Mam_25rs_f1 <- Mam_25rs[Mam_25rs[,"WC_zone"] == 1,]
Mam_25rs_f2 <- Mam_25rs[Mam_25rs[,"WC_zone"] == 2,]
Mam_25rs_f3 <- Mam_25rs[Mam_25rs[,"WC_zone"] == 3,]
Mam_25rs_f4 <- Mam_25rs[Mam_25rs[,"WC_zone"] == 4,]
Mam_25rs_f5 <- Mam_25rs[Mam_25rs[,"WC_zone"] == 5,]

spp_number <- 0
for (i in 3:5703) {
  if(any(unique(Mam_f1[,i]) == 1, na.rm = T)) {
    spp_number <- spp_number + 1
  } 
}

spp_number <- 0
for (i in 3:3561) {
  if(any(unique(Mam_75rs_f1[,i]) == 1, na.rm = T)) {
    spp_number <- spp_number + 1
  } 
}

spp_number <- 0
for (i in 3:2374) {
  if(any(unique(Mam_50rs_f1[,i]) == 1, na.rm = T)) {
    spp_number <- spp_number + 1
  } 
}

spp_number <- 0
for (i in 3:1188) {
  if(any(unique(Mam_25rs_f1[,i]) == 1, na.rm = T)) {
    spp_number <- spp_number + 1
  } 
}


############################### Body size filter ###############################

Combi_5 <- Combi_75rs[!is.na(Combi_75rs$adult_mass_g),] #3515
set.seed(288)
bs75 <- sample(seq_len(nrow(Combi_5)), round(0.75*nrow(Combi_5)), prob = Combi_5$adult_mass_g)
Combi_75bs <- Combi_5[bs75,] #2636
Mam_75bs <- Mam_75rs[,c("X", "Y",Combi_75bs$iucn2020_binomial, "Land_mass", "Clim_zone",
                        "GUM_USC", "WC_zone", "Richness", "Land", "Area", "Fossils", "Richness2")]#2647  
identical(sort(colnames(Mam_75bs[,3:2638])), sort(Combi_75bs$iucn2020_binomial)) #T
richness <- apply(Mam_75bs[,3:2638], 1, function(x) length(x[x>0]))
Mam_75bs[,2648] <- richness
colnames(Mam_75bs)[colnames(Mam_75bs) == "V2648"] <- "Richness3" #New column
Mam_75bs[Mam_75bs$Richness3 == 0, 2648] <- NA


Combi_5 <- Combi_50rs[!is.na(Combi_50rs$adult_mass_g),] #2350
bs50 <- sample(seq_len(nrow(Combi_5)), round(0.50*nrow(Combi_5)), prob = Combi_5$adult_mass_g)
Combi_50bs <- Combi_5[bs50,] #1175
Mam_50bs <- Mam_50rs[,c("X", "Y",Combi_50bs$iucn2020_binomial, "Land_mass", "Clim_zone",
                        "GUM_USC", "WC_zone", "Richness", "Land", "Area", "Fossils", "Richness2")]#1186 
identical(sort(colnames(Mam_50bs[,3:1177])), sort(Combi_50bs$iucn2020_binomial)) #T
richness <- apply(Mam_50bs[,3:1177], 1, function(x) length(x[x>0]))
Mam_50bs[,1187] <- richness
colnames(Mam_50bs)[colnames(Mam_50bs) == "V1187"] <- "Richness3" #New column
Mam_50bs[Mam_50bs$Richness3 == 0, 1187] <- NA


Combi_5 <- Combi_25rs[!is.na(Combi_25rs$adult_mass_g),] #1179
bs25 <- sample(seq_len(nrow(Combi_5)), round(0.25*nrow(Combi_5)), prob = Combi_5$adult_mass_g)
Combi_25bs <- Combi_5[bs25,] #295
Mam_25bs <- Mam_25rs[,c("X", "Y",Combi_25bs$iucn2020_binomial, "Land_mass", "Clim_zone",
                        "GUM_USC", "WC_zone", "Richness", "Land", "Area", "Fossils", "Richness2")]#306 
identical(sort(colnames(Mam_25bs[,3:297])), sort(Combi_25bs$iucn2020_binomial)) #T
richness <- apply(Mam_25bs[,3:297], 1, function(x) length(x[x>0]))
Mam_25bs[,307] <- richness
colnames(Mam_25bs)[colnames(Mam_25bs) == "V307"] <- "Richness3" #New column
Mam_25bs[Mam_25bs$Richness3 == 0, 307] <- NA

Mam_25bs_fr <- Mam_25bs[,c("X", "Y", "Richness3")]
Mam_25bs_raster <- terra::rasterize(x = as.matrix(Mam_25bs_fr[,1:2]), 
                                    y = as(x, "SpatRaster"), 
                                    values = Mam_25bs_fr[,3], update = T)
crs(Mam_25bs_raster) <- "+proj=longlat +datum=WGS84 +no_defs"

Mam_50bs_fr <- Mam_50bs[,c("X", "Y", "Richness3")]
Mam_50bs_raster <- terra::rasterize(x = as.matrix(Mam_50bs_fr[,1:2]), 
                                    y = as(x, "SpatRaster"), 
                                    values = Mam_50bs_fr[,3], update = T)
crs(Mam_50bs_raster) <- "+proj=longlat +datum=WGS84 +no_defs"

Mam_75bs_fr <- Mam_75bs[,c("X", "Y", "Richness3")]
Mam_75bs_raster <- terra::rasterize(x = as.matrix(Mam_75bs_fr[,1:2]), 
                                    y = as(x, "SpatRaster"), 
                                    values = Mam_75bs_fr[,3], update = T)
crs(Mam_75bs_raster) <- "+proj=longlat +datum=WGS84 +no_defs"

# Richness before/after
Mam_75bs$WC_zone <- as.factor(Mam_75bs$WC_zone)
Mam_50bs$WC_zone <- as.factor(Mam_50bs$WC_zone)
Mam_25bs$WC_zone <- as.factor(Mam_25bs$WC_zone)

Mam_75bs_f1 <- Mam_75bs[Mam_75bs[,"WC_zone"] == 1,]
Mam_75bs_f2 <- Mam_75bs[Mam_75bs[,"WC_zone"] == 2,]
Mam_75bs_f3 <- Mam_75bs[Mam_75bs[,"WC_zone"] == 3,]
Mam_75bs_f4 <- Mam_75bs[Mam_75bs[,"WC_zone"] == 4,]
Mam_75bs_f5 <- Mam_75bs[Mam_75bs[,"WC_zone"] == 5,]

Mam_50bs_f1 <- Mam_50bs[Mam_50bs[,"WC_zone"] == 1,]
Mam_50bs_f2 <- Mam_50bs[Mam_50bs[,"WC_zone"] == 2,]
Mam_50bs_f3 <- Mam_50bs[Mam_50bs[,"WC_zone"] == 3,]
Mam_50bs_f4 <- Mam_50bs[Mam_50bs[,"WC_zone"] == 4,]
Mam_50bs_f5 <- Mam_50bs[Mam_50bs[,"WC_zone"] == 5,]

Mam_25bs_f1 <- Mam_25bs[Mam_25bs[,"WC_zone"] == 1,]
Mam_25bs_f2 <- Mam_25bs[Mam_25bs[,"WC_zone"] == 2,]
Mam_25bs_f3 <- Mam_25bs[Mam_25bs[,"WC_zone"] == 3,]
Mam_25bs_f4 <- Mam_25bs[Mam_25bs[,"WC_zone"] == 4,]
Mam_25bs_f5 <- Mam_25bs[Mam_25bs[,"WC_zone"] == 5,]

spp_number <- 0
for (i in 3:2638) {
  if(any(unique(Mam_75bs_f5[,i]) == 1, na.rm = T)) {
    spp_number <- spp_number + 1
  } 
}

spp_number <- 0
for (i in 3:1177) {
  if(any(unique(Mam_50bs_f5[,i]) == 1, na.rm = T)) {
    spp_number <- spp_number + 1
  } 
}

spp_number <- 0
for (i in 3:297) {
  if(any(unique(Mam_25bs_f4[,i]) == 1, na.rm = T)) {
    spp_number <- spp_number + 1
  } 
}

########################## Human-spatial filter ################################
Mam_75bs$WC_zone <- as.numeric(Mam_75bs$WC_zone)
Mam_75bs[is.na(Mam_75bs$Fossils),2646] <- 0
Mam_75bs[Mam_75bs$Fossils >= 1,2646] <- 1
Mam_75fos <- Mam_75bs[Mam_75bs$Fossils == 1,]
length(which(colSums(Mam_75fos[,3:2638]) == 0)) #614 losing spp, 2636-614 = 2022 spp
Mam_75fos_df <- Mam_75fos[,c("X", "Y", "Richness3")]
Mam_75fos_rast <- rasterize(x = as.matrix(Mam_75fos_df[,1:2]), y = as(x, "SpatRaster"), 
                            values = Mam_75fos_df[,3], update = T)
crs(Mam_75fos_rast) <- "+proj=longlat +datum=WGS84 +no_defs"

Mam_50bs$WC_zone <- as.numeric(Mam_50bs$WC_zone)
Mam_50bs[is.na(Mam_50bs$Fossils),1185] <- 0
Mam_50bs[Mam_50bs$Fossils >= 1,1185] <- 1
Mam_50fos <- Mam_50bs[Mam_50bs$Fossils == 1,]
length(which(colSums(Mam_50fos[,3:1177]) == 0)) #153 losing spp, 1175-153 = 1020 spp
Mam_50fos_df <- Mam_50fos[,c("X", "Y", "Richness3")]
Mam_50fos_rast <- rasterize(x = as.matrix(Mam_50fos_df[,1:2]), y = as(x, "SpatRaster"), 
                            values = Mam_50fos_df[,3], update = T)
crs(Mam_50fos_rast) <- "+proj=longlat +datum=WGS84 +no_defs"

Mam_25bs$WC_zone <- as.numeric(Mam_25bs$WC_zone)
Mam_25bs[is.na(Mam_25bs$Fossils),305] <- 0
Mam_25bs[Mam_25bs$Fossils >= 1,305] <- 1
Mam_25fos <- Mam_25bs[Mam_25bs$Fossils == 1,]
length(which(colSums(Mam_25fos[,3:297]) == 0)) #21 losing spp, 295-21 = 274 spp
Mam_25fos_df <- Mam_25fos[,c("X", "Y", "Richness3")]
Mam_25fos_rast <- rasterize(x = as.matrix(Mam_25fos_df[,1:2]), y = as(x, "SpatRaster"), 
                            values = Mam_25fos_df[,3], update = T)
crs(Mam_25fos_rast) <- "+proj=longlat +datum=WGS84 +no_defs"

# Richness before/after
Mam_75fos$WC_zone <- as.factor(Mam_75fos$WC_zone)
Mam_50fos$WC_zone <- as.factor(Mam_50fos$WC_zone)
Mam_25fos$WC_zone <- as.factor(Mam_25fos$WC_zone)

Mam_75fos_f1 <- Mam_75fos[Mam_75fos[,"WC_zone"] == 1,]
Mam_75fos_f2 <- Mam_75fos[Mam_75fos[,"WC_zone"] == 2,]
Mam_75fos_f3 <- Mam_75fos[Mam_75fos[,"WC_zone"] == 3,]
Mam_75fos_f4 <- Mam_75fos[Mam_75fos[,"WC_zone"] == 4,]
Mam_75fos_f5 <- Mam_75fos[Mam_75fos[,"WC_zone"] == 5,]

Mam_50fos_f1 <- Mam_50fos[Mam_50fos[,"WC_zone"] == 1,]
Mam_50fos_f2 <- Mam_50fos[Mam_50fos[,"WC_zone"] == 2,]
Mam_50fos_f3 <- Mam_50fos[Mam_50fos[,"WC_zone"] == 3,]
Mam_50fos_f4 <- Mam_50fos[Mam_50fos[,"WC_zone"] == 4,]
Mam_50fos_f5 <- Mam_50fos[Mam_50fos[,"WC_zone"] == 5,]

Mam_25fos_f1 <- Mam_25fos[Mam_25fos[,"WC_zone"] == 1,]
Mam_25fos_f2 <- Mam_25fos[Mam_25fos[,"WC_zone"] == 2,]
Mam_25fos_f3 <- Mam_25fos[Mam_25fos[,"WC_zone"] == 3,]
Mam_25fos_f4 <- Mam_25fos[Mam_25fos[,"WC_zone"] == 4,]
Mam_25fos_f5 <- Mam_25fos[Mam_25fos[,"WC_zone"] == 5,]

spp_number <- 0
for (i in 3:2638) {
  if(any(unique(Mam_75fos_f5[,i]) == 1, na.rm = T)) {
    spp_number <- spp_number + 1
  } 
}

spp_number <- 0
for (i in 3:1177) {
  if(any(unique(Mam_50fos_f5[,i]) == 1, na.rm = T)) {
    spp_number <- spp_number + 1
  } 
}

spp_number <- 0
for (i in 3:297) {
  if(any(unique(Mam_25fos_f5[,i]) == 1, na.rm = T)) {
    spp_number <- spp_number + 1
  } 
}

########################### 2 taxonomic filters ################################

#Fixing taxonomic names
Mam_pbdb[Mam_pbdb$genus == "Aotus", "family"] <- "Aotidae"
Mam_pbdb[Mam_pbdb$genus == "Xenothrix", "family"] <- "Aotidae"
Mam_pbdb[Mam_pbdb$genus == "Insulacebus", "family"] <- "Aotidae"

Mam_pbdb[Mam_pbdb$family == "Aplodontidae", "family"] <- "Aplodontiidae"

Mam_pbdb[Mam_pbdb$genus == "Patasola", "family"] <- "Callitrichidae"
Mam_pbdb[Mam_pbdb$genus == "Lagonimico", "family"] <- "Callitrichidae"
Mam_pbdb[Mam_pbdb$genus == "Micodon", "family"] <- "Callitrichidae"
Mam_pbdb[Mam_pbdb$genus == "Cebuella", "family"] <- "Callitrichidae"
Mam_pbdb[Mam_pbdb$genus == "Saguinus", "family"] <- "Callitrichidae"
Mam_pbdb[Mam_pbdb$genus == "Callimico", "family"] <- "Callitrichidae"

#Cistugidae = NO 
Mam_pbdb[Mam_pbdb$family == "Galeopithecidae", "family"] <- "Cynocephalidae"

Mam_pbdb[Mam_pbdb$genus == "Galidia", "family"] <- "Eupleridae"

Mam_pbdb[Mam_pbdb$genus == "Miniopterus", "family"] <- "Miniopteridae"

#Nandiniidae = NO
Mam_pbdb[Mam_pbdb$genus == "Neocometes", "family"] <- "Platacanthomyidae"
Mam_pbdb[Mam_pbdb$genus == "Platacanthomys", "family"] <- "Platacanthomyidae"
Mam_pbdb[Mam_pbdb$genus == "Typhlomys", "family"] <- "Platacanthomyidae"


### 75%
Mam_fos_names <- names(which(colSums(Mam_75fos[,3:2638]) != 0))
Combi_75fos <- Combi_75bs[Combi_75bs$iucn2020_binomial %in% Mam_fos_names,] #2022 spp
spp_Combir <- sort(setdiff(unique(Combi_75fos$family), unique(Mam_pbdb$family)))
#3 families in Combi_75bs are not in Mam_pbdb
sort(table(Combi_75fos[Combi_75fos$family %in% spp_Combir, "family"])) #En total, 4 spp

spp_Mampbdb <- sort(setdiff(unique(Mam_pbdb$family), unique(Combi_75fos$family)))
#488 families in Mam_pbdb are not in Combi_75fos

Combi_75fos <- Combi_75fos[!Combi_75fos$family %in% spp_Combir,] #2018 obs, 118 families, -4 spp 
Mam_pbdb3 <- Mam_pbdb[!Mam_pbdb$family %in% spp_Mampbdb,] #488 lost species = 118 families
Mam_pbdb3 <- as.data.frame(table(Mam_pbdb3$family))
colnames(Mam_pbdb3) <- c("family", "n_obs")
Mam_pbdb3$family <- as.character(Mam_pbdb3$family)
Combi_75fos$Fam_obs <- NA
family_names <- unique(Mam_pbdb3$family)
for (i in 1:length(family_names)) {
  Combi_75fos[Combi_75fos$family == family_names[[i]], 
          "Fam_obs"] <- Mam_pbdb3[Mam_pbdb3$family == family_names[[i]], "n_obs"]
}

## Taxonomic weighted
set.seed(288)
tax75 <- sample(seq_len(nrow(Combi_6)), round(0.75*nrow(Combi_6)), prob = Combi_6$Fam_obs)
Combi_75tax <- Combi_6[tax75,] #1514
Mam_75tax <- Mam_75fos[,c("X", "Y",Combi_75tax$iucn2020_binomial, "Land_mass", "Clim_zone",
                          "GUM_USC", "WC_zone", "Richness", "Land", "Area", "Fossils",
                          "Richness2","Richness3")] #1526
identical(sort(colnames(Mam_75tax[,3:1516])), sort(Combi_75tax$iucn2020_binomial)) #T

richness <- apply(Mam_75tax[,3:1516], 1, function(x) length(x[x>0]))
Mam_75tax[,1527] <- richness
colnames(Mam_75tax)[colnames(Mam_75tax) == "V1527"] <- "Richness4" #New column
Mam_75tax[Mam_75tax$Richness4 == 0, 1527] <- NA

Mam_75tax_fr <- Mam_75tax[,c("X", "Y", "Richness4")]
Mam_75tax_raster <- terra::rasterize(x = as.matrix(Mam_75tax_fr[,1:2]), y = as(x, "SpatRaster"), 
                                     values = Mam_75tax_fr[,3], update = T)
crs(Mam_75tax_raster) <- "+proj=longlat +datum=WGS84 +no_defs"

## Taxonomic "whole-families"
0.75*nrow(Mam_pbdb3) #88.5 = 9 families
sort(Mam_pbdb3$n_obs, decreasing = T)[[89]]
my_families <- Mam_pbdb3[Mam_pbdb3$n_obs >=23,1]
Combi_7 <- Combi_6[Combi_6$family %in% my_families,]
Mam_75tax2 <- Mam_75fos[,c("X", "Y",Combi_7$iucn2020_binomial, "Land_mass", "Clim_zone",
                           "GUM_USC", "WC_zone", "Richness", "Land", "Area", "Fossils", 
                           "Richness2","Richness3")] 
identical(sort(colnames(Mam_75tax2[,3:1902])), sort(Combi_7$iucn2020_binomial)) #T
richness <- apply(Mam_75tax2[,3:1902], 1, function(x) length(x[x>0]))
Mam_75tax2[,1913] <- richness
colnames(Mam_75tax2)[colnames(Mam_75tax2) == "V1913"] <- "Richness4" #New column
Mam_75tax2[Mam_75tax2$Richness4 == 0, 1913] <- NA
Mam_75tax2_fr <- Mam_75tax2[,c("X", "Y", "Richness4")]
Mam_75tax2_raster <- terra::rasterize(x = as.matrix(Mam_75tax2_fr[,1:2]), y = as(x, "SpatRaster"), 
                                      values = Mam_75tax2_fr[,3], update = T)
crs(Mam_75tax2_raster) <- "+proj=longlat +datum=WGS84 +no_defs"


### 50%
Mam_50fos_names <- names(which(colSums(Mam_50fos[,3:1177]) != 0))
Combi_50fos <- Combi_50bs[Combi_50bs$iucn2020_binomial %in% Mam_50fos_names,] #1022 spp
spp_Combi50fos <- sort(setdiff(unique(Combi_50fos$family), unique(Mam_pbdb$family)))
#2 families in Combi_50fos are not in Mam_pbdb
sum(sort(table(Combi_50fos[Combi_50fos$family %in% spp_Combi50fos, "family"]))) 
sort(table(Combi_50fos[Combi_50fos$family %in% spp_Combi50fos, "family"])) # 2 spp
spp_Mampbdb <- sort(setdiff(unique(Mam_pbdb$family), unique(Combi_50fos$family)))
#501 families in Mam_pbdb that are not in Combi_r

Combi_50fos <- Combi_50fos[!Combi_50fos$family %in% spp_Combi50fos,] #1020 obs 105 families, -2 spp
Mam_pbdb50fos <- Mam_pbdb[!Mam_pbdb$family %in% spp_Mampbdb,] #65928 obs, 501 lost families = 105 families

Mam_pbdb50fos <- as.data.frame(table(Mam_pbdb50fos$family))
colnames(Mam_pbdb50fos) <- c("family", "n_obs")
Mam_pbdb50fos$family <- as.character(Mam_pbdb50fos$family)

Combi_50fos$Fam_obs <- NA
family_names <- unique(Mam_pbdb50fos$family)
for (i in 1:length(family_names)) {
  Combi_50fos[Combi_50fos$family == family_names[[i]], 
              "Fam_obs"] <- Mam_pbdb50fos[Mam_pbdb50fos$family == family_names[[i]], "n_obs"]
}

#Taxonomic weighted
set.seed(288)
tax50 <- sample(seq_len(nrow(Combi_50fos)), round(0.50*nrow(Combi_50fos)), prob = Combi_50fos$Fam_obs)
Combi_50tax <- Combi_50fos[tax50,] #510
Mam_50tax <- Mam_50fos[,c("X", "Y",Combi_50tax$iucn2020_binomial, "Land_mass", "Clim_zone",
                          "GUM_USC", "WC_zone", "Richness", "Land", "Area", "Fossils", 
                          "Richness2","Richness3")] 
identical(sort(colnames(Mam_50tax[,3:512])), sort(Combi_50tax$iucn2020_binomial)) #T

richness <- apply(Mam_50tax[,3:512], 1, function(x) length(x[x>0]))
Mam_50tax[,523] <- richness
colnames(Mam_50tax)[colnames(Mam_50tax) == "V523"] <- "Richness4" #New column
Mam_50tax[Mam_50tax$Richness4 == 0, 524] <- NA

Mam_50tax_fr <- Mam_50tax[,c("X", "Y", "Richness4")]
Mam_50tax_raster <- terra::rasterize(x = as.matrix(Mam_50tax_fr[,1:2]), y = as(x, "SpatRaster"), 
                                     values = Mam_50tax_fr[,3], update = T)
crs(Mam_50tax_raster) <- "+proj=longlat +datum=WGS84 +no_defs"

#Taxonomic "whole-families"
round(0.5*nrow(Mam_pbdb50fos)) #52 familias
sort(Mam_pbdb50fos$n_obs, decreasing = T)[[52]]
my_families <- Mam_pbdb50fos[Mam_pbdb50fos$n_obs >=133,1]
Combi_50fos2 <- Combi_50fos[Combi_50fos$family %in% my_families,]

Mam_50tax2 <- Mam_50fos[,c("X", "Y",Combi_50fos2$iucn2020_binomial, "Land_mass", "Clim_zone",
                           "GUM_USC", "WC_zone", "Richness", "Land", "Area", "Fossils", 
                           "Richness2","Richness3")] 
identical(sort(colnames(Mam_50tax2[,3:835])), sort(Combi_50fos2$iucn2020_binomial)) #T

richness <- apply(Mam_50tax2[,3:835], 1, function(x) length(x[x>0]))
Mam_50tax2[,846] <- richness
colnames(Mam_50tax2)[colnames(Mam_50tax2) == "V846"] <- "Richness4" #New column
Mam_50tax2[Mam_50tax2$Richness4 == 0, 846] <- NA

Mam_50tax2_fr <- Mam_50tax2[,c("X", "Y", "Richness4")]
Mam_50tax2_raster <- terra::rasterize(x = as.matrix(Mam_50tax2_fr[,1:2]), y = as(x, "SpatRaster"), 
                                      values = Mam_50tax2_fr[,3], update = T)
crs(Mam_50tax2_raster) <- "+proj=longlat +datum=WGS84 +no_defs"

####################### 25%
Mam_25fos_names <- names(which(colSums(Mam_25fos[,3:297]) != 0))
Combi_25fos <- Combi_25bs[Combi_25bs$iucn2020_binomial %in% Mam_25fos_names,] #274 spp
spp_Combi25fos <- sort(setdiff(unique(Combi_25fos$family), unique(Mam_pbdb$family)))
#1 familie in Combi_25fos is not Mam_pbdb
sort(table(Combi_25fos[Combi_25fos$family %in% spp_Combi25fos, "family"])) #1 spp

spp_Mampbdb <- sort(setdiff(unique(Mam_pbdb$family), unique(Combi_25fos$family)))
#546 families in Mam_pbdb are not in Combi_r

Combi_25fos <- Combi_25fos[!Combi_25fos$family %in% spp_Combi25fos,] #273 obs, 60 families, - 1 spp 
Mam_pbdb25fos <- Mam_pbdb[!Mam_pbdb$family %in% spp_Mampbdb,] #53037 obs, 546 lost families = 60 families

Mam_pbdb25fos <- as.data.frame(table(Mam_pbdb25fos$family))
colnames(Mam_pbdb25fos) <- c("family", "n_obs")
Mam_pbdb25fos$family <- as.character(Mam_pbdb25fos$family)

Combi_25fos$Fam_obs <- NA
family_names <- unique(Mam_pbdb25fos$family)
for (i in 1:length(family_names)) {
  Combi_25fos[Combi_25fos$family == family_names[[i]], "Fam_obs"] <- Mam_pbdb25fos[Mam_pbdb25fos$family == family_names[[i]], "n_obs"]
}


#Taxonomic weighted
set.seed(288)
tax25 <- sample(seq_len(nrow(Combi_25fos)), round(0.25*nrow(Combi_25fos)), prob = Combi_25fos$Fam_obs)
Combi_25tax <- Combi_25fos[tax25,] #68
Mam_25tax <- Mam_25fos[,c("X", "Y",Combi_25tax$iucn2020_binomial, "Land_mass", "Clim_zone",
                          "GUM_USC", "WC_zone", "Richness", "Land", "Area", "Fossils", 
                          "Richness2","Richness3")] 
identical(sort(colnames(Mam_25tax[,3:70])), sort(Combi_25tax$iucn2020_binomial)) #T

richness <- apply(Mam_25tax[,3:70], 1, function(x) length(x[x>0]))
Mam_25tax[,81] <- richness
colnames(Mam_25tax)[colnames(Mam_25tax) == "V81"] <- "Richness4" #New column
Mam_25tax[Mam_25tax$Richness4 == 0, 81] <- NA

Mam_25tax_fr <- Mam_25tax[,c("X", "Y", "Richness4")]
Mam_25tax_raster <- terra::rasterize(x = as.matrix(Mam_25tax_fr[,1:2]), y = as(x, "SpatRaster"), 
                                     values = Mam_25tax_fr[,3], update = T)
crs(Mam_25tax_raster) <- "+proj=longlat +datum=WGS84 +no_defs"

#Taxonomic "whole-families"
round(0.25*nrow(Mam_pbdb25fos)) #15 familias
sort(Mam_pbdb25fos$n_obs, decreasing = T)[[15]]
my_families <- Mam_pbdb25fos[Mam_pbdb25fos$n_obs >=1279,1]
Combi_25fos2 <- Combi_25fos[Combi_25fos$family %in% my_families,]

Mam_25tax2 <- Mam_25fos[,c("X", "Y",Combi_25fos2$iucn2020_binomial, "Land_mass", "Clim_zone",
                           "GUM_USC", "WC_zone", "Richness", "Land", "Area", "Fossils", 
                           "Richness2","Richness3")] 
identical(sort(colnames(Mam_25tax2[,3:161])), sort(Combi_25fos2$iucn2020_binomial)) #T

richness <- apply(Mam_25tax2[,3:161], 1, function(x) length(x[x>0]))
Mam_25tax2[,172] <- richness
colnames(Mam_25tax2)[colnames(Mam_25tax2) == "V172"] <- "Richness4" #New column
Mam_25tax2[Mam_25tax2$Richness4 == 0, 172] <- NA


Mam_25tax2_fr <- Mam_25tax2[,c("X", "Y", "Richness4")]
Mam_25tax2_raster <- terra::rasterize(x = as.matrix(Mam_25tax2_fr[,1:2]), y = as(x, "SpatRaster"), 
                                      values = Mam_25tax2_fr[,3], update = T)
crs(Mam_25tax2_raster) <- "+proj=longlat +datum=WGS84 +no_defs"

### Riqueza antes/despues
Mam_75tax$WC_zone <- as.factor(Mam_75tax$WC_zone)
Mam_50tax$WC_zone <- as.factor(Mam_50tax$WC_zone)
Mam_25tax$WC_zone <- as.factor(Mam_25tax$WC_zone)

Mam_75tax_f1 <- Mam_75tax[Mam_75tax[,"WC_zone"] == 1,]
Mam_75tax_f2 <- Mam_75tax[Mam_75tax[,"WC_zone"] == 2,]
Mam_75tax_f3 <- Mam_75tax[Mam_75tax[,"WC_zone"] == 3,]
Mam_75tax_f4 <- Mam_75tax[Mam_75tax[,"WC_zone"] == 4,]
Mam_75tax_f5 <- Mam_75tax[Mam_75tax[,"WC_zone"] == 5,]
Mam_50tax_f1 <- Mam_50tax[Mam_50tax[,"WC_zone"] == 1,]
Mam_50tax_f2 <- Mam_50tax[Mam_50tax[,"WC_zone"] == 2,]
Mam_50tax_f3 <- Mam_50tax[Mam_50tax[,"WC_zone"] == 3,]
Mam_50tax_f4 <- Mam_50tax[Mam_50tax[,"WC_zone"] == 4,]
Mam_50tax_f5 <- Mam_50tax[Mam_50tax[,"WC_zone"] == 5,]
Mam_25tax_f1 <- Mam_25tax[Mam_25tax[,"WC_zone"] == 1,]
Mam_25tax_f2 <- Mam_25tax[Mam_25tax[,"WC_zone"] == 2,]
Mam_25tax_f3 <- Mam_25tax[Mam_25tax[,"WC_zone"] == 3,]
Mam_25tax_f4 <- Mam_25tax[Mam_25tax[,"WC_zone"] == 4,]
Mam_25tax_f5 <- Mam_25tax[Mam_25tax[,"WC_zone"] == 5,]

spp_number <- 0
for (i in 3:1516) {
  if(any(unique(Mam_75tax_f5[,i]) == 1, na.rm = T)) {
    spp_number <- spp_number + 1
  } 
}

spp_number <- 0
for (i in 3:512) {
  if(any(unique(Mam_50tax_f5[,i]) == 1, na.rm = T)) {
    spp_number <- spp_number + 1
  } 
}

spp_number <- 0
for (i in 3:70) {
  if(any(unique(Mam_25tax_f5[,i]) == 1, na.rm = T)) {
    spp_number <- spp_number + 1
  } 
}


Mam_75tax2$WC_zone <- as.factor(Mam_75tax2$WC_zone)
Mam_50tax2$WC_zone <- as.factor(Mam_50tax2$WC_zone)
Mam_25tax2$WC_zone <- as.factor(Mam_25tax2$WC_zone)

Mam_75tax2_f1 <- Mam_75tax2[Mam_75tax2[,"WC_zone"] == 1,]
Mam_75tax2_f2 <- Mam_75tax2[Mam_75tax2[,"WC_zone"] == 2,]
Mam_75tax2_f3 <- Mam_75tax2[Mam_75tax2[,"WC_zone"] == 3,]
Mam_75tax2_f4 <- Mam_75tax2[Mam_75tax2[,"WC_zone"] == 4,]
Mam_75tax2_f5 <- Mam_75tax2[Mam_75tax2[,"WC_zone"] == 5,]
Mam_50tax2_f1 <- Mam_50tax2[Mam_50tax2[,"WC_zone"] == 1,]
Mam_50tax2_f2 <- Mam_50tax2[Mam_50tax2[,"WC_zone"] == 2,]
Mam_50tax2_f3 <- Mam_50tax2[Mam_50tax2[,"WC_zone"] == 3,]
Mam_50tax2_f4 <- Mam_50tax2[Mam_50tax2[,"WC_zone"] == 4,]
Mam_50tax2_f5 <- Mam_50tax2[Mam_50tax2[,"WC_zone"] == 5,]
Mam_25tax2_f1 <- Mam_25tax2[Mam_25tax2[,"WC_zone"] == 1,]
Mam_25tax2_f2 <- Mam_25tax2[Mam_25tax2[,"WC_zone"] == 2,]
Mam_25tax2_f3 <- Mam_25tax2[Mam_25tax2[,"WC_zone"] == 3,]
Mam_25tax2_f4 <- Mam_25tax2[Mam_25tax2[,"WC_zone"] == 4,]
Mam_25tax2_f5 <- Mam_25tax2[Mam_25tax2[,"WC_zone"] == 5,]

spp_number <- 0
for (i in 3:1902) {
  if(any(unique(Mam_75tax2_f5[,i]) == 1, na.rm = T)) {
    spp_number <- spp_number + 1
  } 
}

spp_number <- 0
for (i in 3:835) {
  if(any(unique(Mam_50tax2_f5[,i]) == 1, na.rm = T)) {
    spp_number <- spp_number + 1
  } 
}

spp_number <- 0
for (i in 3:161) {
  if(any(unique(Mam_25tax2_f5[,i]) == 1, na.rm = T)) {
    spp_number <- spp_number + 1
  } 
}


############################## Plot LBG graphs #################################
mam_r_moll <- terra::project(mam_raster, 
                             y = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84",
                             method = "near",
                             res = 50080)
mam_f_moll <- terra::project(mam_fil_raster, 
                             y = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84",
                             method = "near",
                             res = 50080)
Mam_75rs_moll <- terra::project(Mam_75rs_raster, 
                                y = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84",
                                method = "near",
                                res = 50080)
Mam_50rs_moll <- terra::project(Mam_50rs_raster, 
                                y = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84",
                                method = "near",
                                res = 50080)
Mam_25rs_moll <- terra::project(Mam_25rs_raster, 
                                y = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84",
                                method = "near",
                                res = 50080)
Mam_75bs_moll <- terra::project(Mam_75bs_raster, 
                                y = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84",
                                method = "near",
                                res = 50080)
Mam_50bs_moll <- terra::project(Mam_50bs_raster, 
                                y = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84",
                                method = "near",
                                res = 50080)
Mam_25bs_moll <- terra::project(Mam_25bs_raster, 
                                y = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84",
                                method = "near",
                                res = 50080)
Mam_75fos_moll <- terra::project(Mam_75fos_rast, 
                                 y = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84",
                                 method = "near",
                                 res = 50080)
Mam_50fos_moll <- terra::project(Mam_50fos_rast, 
                                 y = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84",
                                 method = "near",
                                 res = 50080)
Mam_25fos_moll <- terra::project(Mam_25fos_rast, 
                                 y = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84",
                                 method = "near",
                                 res = 50080)
Mam_75tax_moll <- terra::project(Mam_75tax_raster, 
                                 y = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84",
                                 method = "near",
                                 res = 50080)
Mam_50tax_moll <- terra::project(Mam_50tax_raster, 
                                 y = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84",
                                 method = "near",
                                 res = 50080)
Mam_25tax_moll <- terra::project(Mam_25tax_raster, 
                                 y = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84",
                                 method = "near",
                                 res = 50080)
Mam_75tax2_moll <- terra::project(Mam_75tax2_raster, 
                                  y = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84",
                                  method = "near",
                                  res = 50080)
Mam_50tax2_moll <- terra::project(Mam_50tax2_raster, 
                                  y = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84",
                                  method = "near",
                                  res = 50080)
Mam_25tax2_moll <- terra::project(Mam_25tax2_raster, 
                                  y = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84",
                                  method = "near",
                                  res = 50080)

mam_r_moll_df <- terra::as.data.frame(mam_r_moll, xy = T, na.rm = F)
mam_r_moll_df[is.nan(mam_r_moll_df$lyr.1), 3] <- NA
mam_f_moll_df <- terra::as.data.frame(mam_f_moll, xy = T, na.rm = F)
mam_f_moll_df[is.nan(mam_f_moll_df$lyr.1), 3] <- NA
Mam75rs_moll_df <- terra::as.data.frame(Mam_75rs_moll, xy = T, na.rm = F)
Mam75rs_moll_df[is.nan(Mam75rs_moll_df$lyr.1), 3] <- NA
Mam50rs_moll_df <- terra::as.data.frame(Mam_50rs_moll, xy = T, na.rm = F)
Mam50rs_moll_df[is.nan(Mam50rs_moll_df$lyr.1), 3] <- NA
Mam25rs_moll_df <- terra::as.data.frame(Mam_25rs_moll, xy = T, na.rm = F)
Mam25rs_moll_df[is.nan(Mam25rs_moll_df$lyr.1), 3] <- NA
Mam75bs_moll_df <- terra::as.data.frame(Mam_75bs_moll, xy = T, na.rm = F)
Mam75bs_moll_df[is.nan(Mam75bs_moll_df$lyr.1), 3] <- NA
Mam50bs_moll_df <- terra::as.data.frame(Mam_50bs_moll, xy = T, na.rm = F)
Mam50bs_moll_df[is.nan(Mam50bs_moll_df$lyr.1), 3] <- NA
Mam25bs_moll_df <- terra::as.data.frame(Mam_25bs_moll, xy = T, na.rm = F)
Mam25bs_moll_df[is.nan(Mam25bs_moll_df$lyr.1), 3] <- NA
Mam75bsfos_moll_df <- terra::as.data.frame(Mam_75fos_moll, xy = T, na.rm = F)
Mam75bsfos_moll_df[is.nan(Mam75bsfos_moll_df$lyr.1), 3] <- NA
Mam50bsfos_moll_df <- terra::as.data.frame(Mam_50fos_moll, xy = T, na.rm = F)
Mam50bsfos_moll_df[is.nan(Mam50bsfos_moll_df$lyr.1), 3] <- NA
Mam25bsfos_moll_df <- terra::as.data.frame(Mam_25fos_moll, xy = T, na.rm = F)
Mam25bsfos_moll_df[is.nan(Mam25bsfos_moll_df$lyr.1), 3] <- NA
Mam_75tax_moll_df <- terra::as.data.frame(Mam_75tax_moll, xy = T, na.rm = F)
Mam_75tax_moll_df[is.nan(Mam_75tax_moll_df$lyr.1), 3] <- NA
Mam_50tax_moll_df <- terra::as.data.frame(Mam_50tax_moll, xy = T, na.rm = F)
Mam_50tax_moll_df[is.nan(Mam_50tax_moll_df$lyr.1), 3] <- NA
Mam_25tax_moll_df <- terra::as.data.frame(Mam_25tax_moll, xy = T, na.rm = F)
Mam_25tax_moll_df[is.nan(Mam_25tax_moll_df$lyr.1), 3] <- NA
Mam_75tax2_moll_df <- terra::as.data.frame(Mam_75tax2_moll, xy = T, na.rm = F)
Mam_75tax2_moll_df[is.nan(Mam_75tax2_moll_df$lyr.1), 3] <- NA
Mam_50tax2_moll_df <- terra::as.data.frame(Mam_50tax2_moll, xy = T, na.rm = F)
Mam_50tax2_moll_df[is.nan(Mam_50tax2_moll_df$lyr.1), 3] <- NA
Mam_25tax2_moll_df <- terra::as.data.frame(Mam_25tax2_moll, xy = T, na.rm = F)
Mam_25tax2_moll_df[is.nan(Mam_25tax2_moll_df$lyr.1), 3] <- NA


coordinates(mam_r_moll_df) <- c("x", "y")
proj4string(mam_r_moll_df) <- CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84")
coordinates(mam_f_moll_df) <- c("x", "y")
proj4string(mam_f_moll_df) <- CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84")
coordinates(Mam75rs_moll_df) <- c("x", "y")
proj4string(Mam75rs_moll_df) <- CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84")
coordinates(Mam50rs_moll_df) <- c("x", "y")
proj4string(Mam50rs_moll_df) <- CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84")
coordinates(Mam25rs_moll_df) <- c("x", "y")
proj4string(Mam25rs_moll_df) <- CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84")
coordinates(Mam75bs_moll_df) <- c("x", "y")
proj4string(Mam75bs_moll_df) <- CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84")
coordinates(Mam50bs_moll_df) <- c("x", "y")
proj4string(Mam50bs_moll_df) <- CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84")
coordinates(Mam25bs_moll_df) <- c("x", "y")
proj4string(Mam25bs_moll_df) <- CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84")
coordinates(Mam75bsfos_moll_df) <- c("x", "y")
proj4string(Mam75bsfos_moll_df) <- CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84")
coordinates(Mam50bsfos_moll_df) <- c("x", "y")
proj4string(Mam50bsfos_moll_df) <- CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84")
coordinates(Mam25bsfos_moll_df) <- c("x", "y")
proj4string(Mam25bsfos_moll_df) <- CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84")
coordinates(Mam_75tax_moll_df) <- c("x", "y")
proj4string(Mam_75tax_moll_df) <- CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84")
coordinates(Mam_50tax_moll_df) <- c("x", "y")
proj4string(Mam_50tax_moll_df) <- CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84")
coordinates(Mam_25tax_moll_df) <- c("x", "y")
proj4string(Mam_25tax_moll_df) <- CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84")
coordinates(Mam_75tax2_moll_df) <- c("x", "y")
proj4string(Mam_75tax2_moll_df) <- CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84")
coordinates(Mam_50tax2_moll_df) <- c("x", "y")
proj4string(Mam_50tax2_moll_df) <- CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84")
coordinates(Mam_25tax2_moll_df) <- c("x", "y")
proj4string(Mam_25tax2_moll_df) <- CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84")

mam_r_moll_sf <- st_as_sf(mam_r_moll_df)
mam_r_longlat_sf <- st_transform(mam_r_moll_sf, "+proj=longlat +datum=WGS84 +no_defs")
mam_f_moll_sf <- st_as_sf(mam_f_moll_df)
mam_f_longlat_sf <- st_transform(mam_f_moll_sf, "+proj=longlat +datum=WGS84 +no_defs")
Mam75rs_moll_sf <- st_as_sf(Mam75rs_moll_df)
Mam75rs_longlat_sf <- st_transform(Mam75rs_moll_sf, "+proj=longlat +datum=WGS84 +no_defs")
Mam50rs_moll_sf <- st_as_sf(Mam50rs_moll_df)
Mam50rs_longlat_sf <- st_transform(Mam50rs_moll_sf, "+proj=longlat +datum=WGS84 +no_defs")
Mam25rs_moll_sf <- st_as_sf(Mam25rs_moll_df)
Mam25rs_longlat_sf <- st_transform(Mam25rs_moll_sf, "+proj=longlat +datum=WGS84 +no_defs")
Mam75bs_moll_sf <- st_as_sf(Mam75bs_moll_df)
Mam75bs_longlat_sf <- st_transform(Mam75bs_moll_sf, "+proj=longlat +datum=WGS84 +no_defs")
Mam50bs_moll_sf <- st_as_sf(Mam50bs_moll_df)
Mam50bs_longlat_sf <- st_transform(Mam50bs_moll_sf, "+proj=longlat +datum=WGS84 +no_defs")
Mam25bs_moll_sf <- st_as_sf(Mam25bs_moll_df)
Mam25bs_longlat_sf <- st_transform(Mam25bs_moll_sf, "+proj=longlat +datum=WGS84 +no_defs")
Mam75bsfos_moll_sf <- st_as_sf(Mam75bsfos_moll_df)
Mam75bsfos_longlat_sf <- st_transform(Mam75bsfos_moll_sf, "+proj=longlat +datum=WGS84 +no_defs")
Mam50bsfos_moll_sf <- st_as_sf(Mam50bsfos_moll_df)
Mam50bsfos_longlat_sf <- st_transform(Mam50bsfos_moll_sf, "+proj=longlat +datum=WGS84 +no_defs")
Mam25bsfos_moll_sf <- st_as_sf(Mam25bsfos_moll_df)
Mam25bsfos_longlat_sf <- st_transform(Mam25bsfos_moll_sf, "+proj=longlat +datum=WGS84 +no_defs")

Mam_75tax_moll_sf <- st_as_sf(Mam_75tax_moll_df)
Mam75tax_longlat_sf <- st_transform(Mam_75tax_moll_sf, "+proj=longlat +datum=WGS84 +no_defs")
Mam_50tax_moll_sf <- st_as_sf(Mam_50tax_moll_df)
Mam50tax_longlat_sf <- st_transform(Mam_50tax_moll_sf, "+proj=longlat +datum=WGS84 +no_defs")
Mam_25tax_moll_sf <- st_as_sf(Mam_25tax_moll_df)
Mam25tax_longlat_sf <- st_transform(Mam_25tax_moll_sf, "+proj=longlat +datum=WGS84 +no_defs")
Mam_75tax2_moll_sf <- st_as_sf(Mam_75tax2_moll_df)
Mam75tax2_longlat_sf <- st_transform(Mam_75tax2_moll_sf, "+proj=longlat +datum=WGS84 +no_defs")
Mam_50tax2_moll_sf <- st_as_sf(Mam_50tax2_moll_df)
Mam50tax2_longlat_sf <- st_transform(Mam_50tax2_moll_sf, "+proj=longlat +datum=WGS84 +no_defs")
Mam_25tax2_moll_sf <- st_as_sf(Mam_25tax2_moll_df)
Mam25tax2_longlat_sf <- st_transform(Mam_25tax2_moll_sf, "+proj=longlat +datum=WGS84 +no_defs")

identical(mam_f_longlat_sf$lyr.1, mam_f_moll_sf$lyr.1) #TRUE
mam_r_sf_juntas <- cbind(mam_r_moll_sf, mam_r_longlat_sf$geometry)
mam_r_dfforplot <- as.data.frame(mam_r_sf_juntas$lyr.1)
for (i in 1:259200) {
  print(i)
  mam_r_dfforplot[i,2] <- mam_r_sf_juntas$geometry[[i]][2]
  mam_r_dfforplot[i,3] <- mam_r_sf_juntas$geometry.1[[i]][2]
}
colnames(mam_r_dfforplot) <- c("Richness", "Moll", "Longlat")

#write.table(mam_r_dfforplot, "sf_coords.csv", row.names = F)
mam_r_dfforplot <- read.table("sf_coords.csv")

db_forslopes <- data.frame(mam_r_dfforplot$Longlat, mam_r_dfforplot$Moll,
                           mam_r_dfforplot$Richness, mam_f_moll_sf$last, Mam75rs_moll_sf$last,
                           Mam50rs_moll_sf$last, Mam25rs_moll_sf$last, Mam75bs_moll_sf$last,
                           Mam50bs_moll_sf$last, Mam25bs_moll_sf$last, Mam75bsfos_moll_sf$last,
                           Mam50bsfos_moll_sf$last, Mam25bsfos_moll_sf$last, Mam_75tax_moll_sf$last,
                           Mam_50tax_moll_sf$last, Mam_25tax_moll_sf$last, Mam_75tax2_moll_sf$last,
                           Mam_50tax2_moll_sf$last, Mam_25tax2_moll_sf$last)
colnames(db_forslopes) <- c("Longlat", "Moll", "Original", "UcS", "rs75", "rs50",
                            "rs25", "bs75", "bs50", "bs25", "foss75", "foss50",
                            "foss25", "tax75", "tax50", "tax25", "tax2_75", 
                            "tax2_50", "tax2_25")

#setwd("./Figures/Figure_LDGplots")
#pdf("LDG_Mam_taxwholefam25_2.pdf", useDingbats = F, width = 10, height = 7)
ggplot(db_forslopes, aes(Longlat, tax2_25, color = tax2_25)) +
  geom_point(shape = 16,
             size = 3, alpha = 0.9) +
  scale_colour_gradientn(colours = rev(paletteer_c("grDevices::YlOrBr", 400))[70:400], limits = c(1,239)) +
  xlab("") +
  ylim(0,250) +
  ylab("") +
  theme_classic() +
  geom_vline(xintercept = 0, color = "grey70", alpha = 0.5) +
  theme(legend.position = "none")
#dev.off()

########################## Plot Taxonomic richnes map ##########################
