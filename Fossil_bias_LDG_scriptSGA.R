# ------------------------------------------------------------------------------
# Publication: Biases can't hide conspicuous biodiversity patterns 
# Last updated: 2024-07-11
# Author: Sofia Galvan
# Email: sofia.galvan@uvigo.es
# Repository: https://github.com/SofiaGalv/Fossil_bias_and_LBG_detection.git
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
tax75 <- sample(seq_len(nrow(Combi_75fos)), round(0.75*nrow(Combi_75fos)), prob = Combi_75fos$Fam_obs)
Combi_75tax <- Combi_75fos[tax75,] #1514
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
Combi_fos2 <- Combi_75fos[Combi_75fos$family %in% my_families,]
Mam_75tax2 <- Mam_75fos[,c("X", "Y",Combi_fos2$iucn2020_binomial, "Land_mass", "Clim_zone",
                           "GUM_USC", "WC_zone", "Richness", "Land", "Area", "Fossils", 
                           "Richness2","Richness3")] 
identical(sort(colnames(Mam_75tax2[,3:1902])), sort(Combi_fos2$iucn2020_binomial)) #T
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


#All richness values along the filter are stored in "richness_perzone_UcS2.csv"

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

#identical(mam_f_longlat_sf$lyr.1, mam_f_moll_sf$lyr.1) #TRUE
#mam_r_sf_juntas <- cbind(mam_r_moll_sf, mam_r_longlat_sf$geometry)
#mam_r_dfforplot <- as.data.frame(mam_r_sf_juntas$lyr.1)
#for (i in 1:259200) {
#  print(i)
#  mam_r_dfforplot[i,2] <- mam_r_sf_juntas$geometry[[i]][2]
#  mam_r_dfforplot[i,3] <- mam_r_sf_juntas$geometry.1[[i]][2]
#}
#colnames(mam_r_dfforplot) <- c("Richness", "Moll", "Longlat")

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
#setwd("./Figures/Figure_richness")
rast_bc <- rast(nrows = 360, ncols = 720, resolution = 0.5, vals = 1,
                extent = ext(-180, 180, -90, 90))
rast_bc <- terra::project(rast_bc, 
                          y = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84",
                          method = "near",
                          res = 50080)
sp_bc <- as.polygons(rast_bc)

#pdf("Mam_75_taxwholefam.pdf", useDingbats = F)
plot(sp_bc, axes = F)
plot(land_moll$geometry, col = "grey30", border = "grey30", lwd=0.4, add = T) #, bg = "skyblue"
terra::plot(Mam_75bs_moll, col = rev(paletteer_c("grDevices::YlOrBr", 400)), axes = F, 
            range = c(0,250), pax = list(labels = F, tick = F), legend = F, add = T)
plot(land_moll$geometry, border = "grey30", add = T, lwd=0.4)
#dev.off()


########################## Calculating area loss ###############################
wc_moll <- terra::project(wc_raster, 
                          y = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84",
                          method = "near", res = 50080)
wc_molldf <- terra::as.data.frame(wc_moll, xy = T, na.rm = F)
wc_molldf[is.nan(wc_molldf$lyr.1), 3] <- NA
npixel_original <- table(wc_molldf$lyr.1)

GUM_moll <- terra::project(GUM_new, 
                           y = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84",
                           method = "near", res = 50080)
wc_extracted <- terra::mask(wc_moll, GUM_moll)
wc_ext_df <- terra::as.data.frame(wc_extracted, xy = T, na.rm = F)
npixel_extracted <- table(wc_ext_df$lyr.1)

pixels_df <- data.frame(npixel_original, npixel_extracted)[,c(1,2,4)]
colnames(pixels_df) <- c("WC_zone", "npixel_raw", "npixel_fil")
pixels_df$per_loss <- 100-(round((pixels_df$npixel_fil/pixels_df$npixel_raw)*100,2))
pixels_df$UCS_sed <- round((pixels_df$npixel_fil/pixels_df$npixel_raw)*100,2)

Mam_fos_moll <- terra::project(Mam_raster, 
                               y = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84",
                               method = "near", res = 50080)
wc_fos_extracted <- terra::mask(wc_extracted, Mam_fos_moll)
wc_fos_df <- terra::as.data.frame(wc_fos_extracted, xy = T, na.rm = F)
npixel_fos <- table(wc_fos_df$lyr.1)

pixels_df2 <- data.frame(npixel_extracted, npixel_fos)[,c(1,2,4)]
colnames(pixels_df2) <- c("WC_zone", "npixel_UcS", "npixel_fos")
pixels_df2$per_loss <- 100-(round((pixels_df2$npixel_fos/pixels_df2$npixel_UcS)*100,2))
pixels_df2$UCS_sed <- round((pixels_df2$npixel_fos/pixels_df2$npixel_UcS)*100,2)


######################### Calculating richness loss ############################
Mam_rich_total <- read.csv("richness_perzone_UcS2.csv", header = T, sep = ";")
Mam_rich_total <- Mam_rich_total[1:6,]
for (i in 3:19) {
  Mam_rich_total[,i] <- as.numeric(Mam_rich_total[,i])
}

Mam_rich_total[,20] <- round((1-(Mam_rich_total$UcS/Mam_rich_total$Original))*100, 2)
colnames(Mam_rich_total)[colnames(Mam_rich_total) == "V20"] <- "Original_UcS"

Mam_rich_total[,21] <- round((1-(Mam_rich_total$Rangesize_75/Mam_rich_total$UcS))*100, 2)
colnames(Mam_rich_total)[colnames(Mam_rich_total) == "V21"] <- "UcS_rs75"
Mam_rich_total[,22] <- round((1-(Mam_rich_total$Rangesize_50/Mam_rich_total$UcS))*100, 2)
colnames(Mam_rich_total)[colnames(Mam_rich_total) == "V22"] <- "UcS_rs50"
Mam_rich_total[,23] <- round((1-(Mam_rich_total$Rangesize_25/Mam_rich_total$UcS))*100, 2)
colnames(Mam_rich_total)[colnames(Mam_rich_total) == "V23"] <- "UcS_rs25"

Mam_rich_total[,24] <- round((1-(Mam_rich_total$Bodysize_75/Mam_rich_total$Rangesize_75))*100, 2)
colnames(Mam_rich_total)[colnames(Mam_rich_total) == "V24"] <- "rs75_bs75"
Mam_rich_total[,25] <- round((1-(Mam_rich_total$Bodysize_50/Mam_rich_total$Rangesize_50))*100, 2)
colnames(Mam_rich_total)[colnames(Mam_rich_total) == "V25"] <- "rs50_bs50"
Mam_rich_total[,26] <- round((1-(Mam_rich_total$Bodysize_25/Mam_rich_total$Rangesize_25))*100, 2)
colnames(Mam_rich_total)[colnames(Mam_rich_total) == "V26"] <- "rs25_bs25"

Mam_rich_total[,27] <- round((1-(Mam_rich_total$Fossil_75/Mam_rich_total$Bodysize_75))*100, 2)
colnames(Mam_rich_total)[colnames(Mam_rich_total) == "V27"] <- "bs75_fos75"
Mam_rich_total[,28] <- round((1-(Mam_rich_total$Fossil_50/Mam_rich_total$Bodysize_50))*100, 2)
colnames(Mam_rich_total)[colnames(Mam_rich_total) == "V28"] <- "bs50_fos50"
Mam_rich_total[,29] <- round((1-(Mam_rich_total$Fossil_25/Mam_rich_total$Bodysize_25))*100, 2)
colnames(Mam_rich_total)[colnames(Mam_rich_total) == "V29"] <- "bs25_fos25"

Mam_rich_total[,30] <- round((1-(Mam_rich_total$Taxweight_75/Mam_rich_total$Fossil_75))*100, 2)
colnames(Mam_rich_total)[colnames(Mam_rich_total) == "V30"] <- "fos75_taxw75"
Mam_rich_total[,31] <- round((1-(Mam_rich_total$Taxweight_50/Mam_rich_total$Fossil_50))*100, 2)
colnames(Mam_rich_total)[colnames(Mam_rich_total) == "V31"] <- "fos50_taxw50"
Mam_rich_total[,32] <- round((1-(Mam_rich_total$Taxweight_25/Mam_rich_total$Fossil_25))*100, 2)
colnames(Mam_rich_total)[colnames(Mam_rich_total) == "V32"] <- "fos25_taxw25"

Mam_rich_total[,33] <- round((1-(Mam_rich_total$Taxfamily_75/Mam_rich_total$Fossil_75))*100, 2)
colnames(Mam_rich_total)[colnames(Mam_rich_total) == "V33"] <- "fos75_taxf75"
Mam_rich_total[,34] <- round((1-(Mam_rich_total$Taxfamily_50/Mam_rich_total$Fossil_50))*100, 2)
colnames(Mam_rich_total)[colnames(Mam_rich_total) == "V34"] <- "fos50_taxf50"
Mam_rich_total[,35] <- round((1-(Mam_rich_total$Taxfamily_25/Mam_rich_total$Fossil_25))*100, 2)
colnames(Mam_rich_total)[colnames(Mam_rich_total) == "V35"] <- "fos25_tafw25"

var_toplot <- colnames(Mam_rich_total[,20:35])
labels_toplot <- list(c("Tropical (21.8%)", "Arid (14.0%)", "Temperate (27.0%)", "Cold (22.9%)", "Polar (42.4%)"),
                      c("Tropical (22.4%)", "Arid (13.1%)", "Temperate (17.1%)", "Cold (10.3%)", "Polar (6.8%)"),
                      c("Tropical (47.5%)", "Arid (32.6%)", "Temperate (38.8%)", "Cold (27.9%)", "Polar (20.7%)"),
                      c("Tropical (73.2%)", "Arid (60.7%)", "Temperate (66.5%)", "Cold (57.0%)", "Polar (44.9%)"),
                      c("Tropical (25.7%)", "Arid (27.8%)", "Temperate (27.5%)", "Cold (26.9%)", "Polar (25.4%)"),
                      c("Tropical (49.0%)", "Arid (51.4%)", "Temperate (51.6%)", "Cold (49.5%)", "Polar (50.4%)"),
                      c("Tropical (73.4%)", "Arid (74.1%)", "Temperate (74.3%)", "Cold (74.0%)", "Polar (71.7%)"),
                      c("Tropical (33.7%)", "Arid (24.7%)", "Temperate (32.9%)", "Cold (22.4%)", "Polar (73.7%)"),
                      c("Tropical (21.8%)", "Arid (21.3%)", "Temperate (29.5%)", "Cold (17.9%)", "Polar (78.6%)"),
                      c("Tropical (15.9%)", "Arid (12.7%)", "Temperate (20.2%)", "Cold (14.3%)", "Polar (71.9%)"),
                      c("Tropical (33.5%)", "Arid (20.9%)", "Temperate (21.4%)", "Cold (8.1%)", "Polar (23.1%)"),
                      c("Tropical (57.9%)", "Arid (43.5%)", "Temperate (47.4%)", "Cold (32.4%)", "Polar (33.3%)"),
                      c("Tropical (81.1%)", "Arid (74.5%)", "Temperate (76.3%)", "Cold (67.9%)", "Polar (80.0%)"),
                      c("Tropical (9.4%)", "Arid (2.6%)", "Temperate (2.2%)", "Cold (0.2%)", "Polar (0.0%)"),
                      c("Tropical (24.9%)", "Arid (13.5%)", "Temperate (14.2%)", "Cold (4.2%)", "Polar (10.4%)"),
                      c("Tropical (52.9%)", "Arid (38.5%)", "Temperate (42.5%)", "Cold (17.9%)", "Polar (24.0%)"))
names(labels_toplot) <- colnames(Mam_rich_total[,20:35])

richloss_plots <- list()
for (i in 1:16) {
  mam_df <- Mam_rich_total[,var_toplot[i]][2:6]
  df_forplot <- data.frame(mam_df, labels_toplot[[i]])
  colnames(df_forplot) <- c("value1", "labels_names")
  df_forplot$labels_names <- factor(df_forplot$labels_names, levels = df_forplot$labels_names)
  pl <- ggplot(df_forplot, 
               aes(x = labels_names, y = value1)) +
    geom_col(color = c("#f9d14a", "#ab3329", "#ed968c", "#7c4b73", "#88a0dc"), 
             lwd = 1, fill = c("#f9d14a", "#ab3329", "#ed968c", "#7c4b73", "#88a0dc"), 
             width = 0.6) +
    labs(x = "", y = "Richness loss (%)") +
    theme_minimal() +
    ggtitle(var_toplot[i]) +
    coord_flip()
  richloss_plots[[i]] <- pl
}

#setwd("./Figures/Figure_speciesloss")
#pdf("fos75_taxf75_spploss.pdf", useDingbats = F)
richloss_plots[[14]]
#dev.off()

#################### Correlation between area/richness loss ####################
df_forcor <- as.data.frame(cbind(pixels_df$WC_zone, pixels_df$per_loss, Mam_rich_total$Original_UcS[2:6]))
colnames(df_forcor) <- c("WC_zone", "Area_loss", "Sp_loss")
#setwd("./Figures/Figure_arealoss")
#pdf("Cor_UcS.pdf", useDingbats = F, width = 10, height = 7)
ggplot(data = df_forcor, aes(x = Area_loss, y = Sp_loss)) +
  geom_point(size = 4, alpha = 1, color = c("#f9d14a", "#ab3329", "#ed968c", "#7c4b73", "#88a0dc"))  +
  stat_smooth(method = 'lm', se = TRUE, alpha = 0.1, color = "black") +
  xlim(c(45,100)) + ylim(c(-75, 120)) +
  #scale_x_continuous(limits = c(130, 1350), breaks = seq(100,1400, 300)) +
  #scale_y_continuous(limits = c(-1500,6000)) +
  labs(x = "Area loss (%)", y = "Species loss (%)") +
  theme_minimal() +
  theme(axis.line = element_line(colour = "grey50", size = 0.6), legend.position = 'none')
#dev.off()
mamcor_ric_area <- cor.test(x = df_forcor$Area_loss, 
                            y = df_forcor$Sp_loss, 
                            method = 'pearson')
mamcor_ric_area$estimate
mamcor_ric_area$p.value

df_forcor2 <- as.data.frame(cbind(pixels_df2$WC_zone, pixels_df2$per_loss, Mam_rich_total$bs75_fos75[2:6]))
df_forcor2 <- as.data.frame(cbind(pixels_df2$WC_zone, pixels_df2$per_loss, Mam_rich_total$bs50_fos50[2:6]))
df_forcor2 <- as.data.frame(cbind(pixels_df2$WC_zone, pixels_df2$per_loss, Mam_rich_total$bs25_fos25[2:6]))
colnames(df_forcor2) <- c("WC_zone", "Area_loss", "Sp_loss")
#pdf("Cor_Foss25.pdf", useDingbats = F, width = 10, height = 7)
ggplot(data = df_forcor2, aes(x = Area_loss, y = Sp_loss)) +
  geom_point(size = 4, alpha = 1, color = c("#f9d14a", "#ab3329", "#ed968c", "#7c4b73", "#88a0dc"))  +
  stat_smooth(method = 'lm', se = TRUE, alpha = 0.1, color = "black") +
  xlim(c(45,100)) + ylim(c(-75, 120)) +
  #scale_x_continuous(limits = c(130, 1350), breaks = seq(100,1400, 300)) +
  #scale_y_continuous(limits = c(-1500,6000)) +
  labs(x = "Area loss (%)", y = "Species loss (%)") +
  theme_minimal() +
  theme(axis.line = element_line(colour = "grey50", size = 0.6), legend.position = 'none')
#dev.off()
mamcor_ric_area <- cor.test(x = df_forcor2$Area_loss, 
                            y = df_forcor2$Sp_loss, 
                            method = 'pearson')
mamcor_ric_area$estimate
mamcor_ric_area$p.value



########################### 0.9 quantile regressions ############################
set.seed(134)
seeds <- sample(1:2000, 100, replace = F)

tell <- 10
iterative_qr <- function(bbdd_entrada, var_name) {
  qr_models <- list()
  qr_iterations <- NULL
  my_subdb <- bbdd_entrada[!is.na(bbdd_entrada[,var_name]),]
  for (i in 1:100) {
    set.seed(seeds[i])
    sample1 <- my_subdb[sample(seq_len(nrow(my_subdb)), round(0.50*nrow(my_subdb))),]
    qr1 <- rq(sample1[,var_name] ~ abs(Longlat), data=sample1, tau = 0.9)
    qr_models[[i]] <- summary(qr1, se = "nid")
    qr_iterations <- as.data.frame(rbind(qr_iterations, c(seeds[i], qr_models[[i]]$coefficients[1,]),
                                         c(seeds[i], qr_models[[i]]$coefficients[2,])))
    if (i %% tell ==0) cat(paste("iteration", i, "completed\n"))
  }
  return(qr_iterations)
}

vari_names <- colnames(db_forslopes)[3:19]
qr_all <- list()
for (i in 1:length(vari_names)) {
  print(paste("Running variable called", vari_names[[i]], sep = " "))
  qr_1var <- iterative_qr(db_forslopes, vari_names[[i]])
  qr_all[[i]] <- qr_1var
}

qr_means_slope <- function(ddbb_input) {
  medias <- c()
  CI <- c()
  mins <- c()
  maxs <- c()
  pval_sig <- c()
  pares <- seq(2, 200, by = 2)
  media_value <- mean(ddbb_input[pares,"Value"]) 
  medias <- c(medias, media_value)
  CI_values <- quantile(ddbb_input[pares,"Value"], probs=c(0.025, 0.5, 0.975))
  CI <- rbind(CI, CI_values)
  min_values <- min(ddbb_input[pares,"Value"])
  mins <- c(mins, min_values)
  max_values <- max(ddbb_input[pares,"Value"])
  maxs <- c(maxs, max_values)
  pval_value <- length(which(ddbb_input[pares, "Pr(>|t|)"]< 0.05))
  pval_sig <- c(pval_sig, pval_value)
  qr_summaries <- c(vari_names[[i]], medias, mins, maxs, CI[,1], CI[,2], CI[,3], pval_sig)
  return(qr_summaries)
}

qrmean_all <- NULL
for (i in 1:length(qr_all)) {
  qrmeans <- qr_means_slope(qr_all[[i]])
  qrmean_all <- rbind(qrmean_all, qrmeans)
}
qrmean_all <- as.data.frame(qrmean_all)
colnames(qrmean_all) <- c("Name", "Mean", "Min", "Max", "Lower95", "Median", "Upper95", "N_pvalsig")
for (i in 2:7) {
  qrmean_all[,i] <- as.numeric(qrmean_all[,i])
}

### Slope decay plot
#setwd("./Figures/Figure_09regression")
#pdf("Slope_decress.pdf", useDingbats = F, width = 10, height = 7)
ggplot(qrmean_all[c(1,2,3,6,9,15),], aes(1:6, -Median)) +
  xlab("") +
  #ylim(0,2.5) +
  ylab("") +
  theme_classic() +
  geom_line() +
  geom_line(data = qrmean_all[c(1,2,3,6,9,12),], color = "black") +
  geom_line(data = qrmean_all[c(1,2,4,7,10,16),], color = "honeydew4") +
  geom_line(data = qrmean_all[c(1,2,4,7,10,13),], color = "honeydew4") +
  geom_line(data = qrmean_all[c(1,2,5,8,11,17),], color = "honeydew3") +
  geom_line(data = qrmean_all[c(1,2,5,8,11,14),], color = "honeydew3") +
  geom_ribbon(data = qrmean_all[c(1,2,3,6,9,15),], aes(ymin = -Lower95, ymax=-Upper95), fill = "black", alpha = 0.2) +
  geom_ribbon(data = qrmean_all[c(1,2,3,6,9,12),], aes(ymin = -Lower95, ymax=-Upper95), fill = "black", alpha = 0.2) +
  geom_ribbon(data = qrmean_all[c(1,2,4,7,10,16),], aes(ymin = -Lower95, ymax=-Upper95), fill = "honeydew4", alpha = 0.2) +
  geom_ribbon(data = qrmean_all[c(1,2,4,7,10,13),], aes(ymin = -Lower95, ymax=-Upper95), fill = "honeydew4", alpha = 0.2) +
  geom_ribbon(data = qrmean_all[c(1,2,5,8,11,17),], aes(ymin = -Lower95, ymax=-Upper95), fill = "honeydew3", alpha = 0.2) +
  geom_ribbon(data = qrmean_all[c(1,2,5,8,11,14),], aes(ymin = -Lower95, ymax=-Upper95), fill = "honeydew3", alpha = 0.2)
#dev.off()

qr_means_inter <- function(ddbb_input) {
  medias <- c()
  CI <- c()
  mins <- c()
  maxs <- c()
  impares <- seq(1, 60, by = 2)
  media_value <- mean(ddbb_input[impares,"Value"]) 
  medias <- c(medias, media_value)
  CI_values <- quantile(ddbb_input[impares,"Value"], probs=c(0.025, 0.5, 0.975))
  CI <- rbind(CI, CI_values)
  min_values <- min(ddbb_input[impares,"Value"])
  mins <- c(mins, min_values)
  max_values <- max(ddbb_input[impares,"Value"])
  maxs <- c(maxs, max_values)
  qr_summaries <- c(vari_names[[i]], medias, mins, maxs, CI[,1], CI[,2], CI[,3])
  return(qr_summaries)
}


qrmean_all_inter <- NULL
for (i in 1:length(qr_all)) {
  qrmeans <- qr_means_inter(qr_all[[i]])
  qrmean_all_inter <- rbind(qrmean_all_inter, qrmeans)
}
qrmean_all_inter <- as.data.frame(qrmean_all_inter)
colnames(qrmean_all_inter) <- c("Name", "Mean", "Min", "Max", "Lower95", "Median", "Upper95")
for (i in 2:7) {
  qrmean_all_inter[,i] <- as.numeric(qrmean_all_inter[,i])
}


### 0.9 quantile regressions plots
ribbons_plot <- list()
xlims_plot <- list()
ylims_plot <- list()
for (i in 1:17) {
  print(i)
  val_intcpt <- c(qrmean_all_inter[i,5], qrmean_all_inter[i,7])
  val_ss <- c(qrmean_all[i,5], qrmean_all[i,7])
  ss <- c(val_ss, rev(val_ss))
  p.prueba <- ggplot(mam_r_dfforplot, aes(abs(Longlat), Richness)) + geom_blank() +
    xlim(0, 90) +
    xlab("") +
    ylim(0,200) +
    ylab("") +
    theme_classic() +
    geom_abline(aes(intercept=qrmean_all_inter[i,6], slope=qrmean_all[i,6], color = "1.Original"), lwd=1) +
    geom_abline(aes(intercept=qrmean_all_inter[i,7], slope=qrmean_all[i,7]), linetype = 'dotted', lwd=1) +
    geom_abline(aes(intercept=qrmean_all_inter[i,5], slope=qrmean_all[i,5]), linetype = 'dotted', lwd=1) +
    scale_colour_manual(values = rev(my_colours))
  p_x <- layer_scales(p.prueba)$x$get_limits()
  p_y <- layer_scales(p.prueba)$y$get_limits()
  df <- data.frame(
    x = rep(c(p_x[1] - (p_x[2] - p_x[1]),
              p_x[2] + (p_x[2] - p_x[1])), each = 2),
    intcpt = c(val_intcpt, rev(val_intcpt))
  ) %>%
    mutate(y = intcpt + ss * x)
  ribbons_plot[[i]] <- df
  xlims_plot[[i]] <- p_x
  ylims_plot[[i]] <- p_y
}

#Y ploteamos
p <- ggplot(mam_r_dfforplot, aes(abs(Longlat), Richness)) + geom_blank() +
  #xlim(0, 200) +
  xlab("") +
  #ylim(0,200) +
  ylab("") +
  theme_classic() +
  geom_abline(aes(intercept=qrmean_all_inter[1,6], slope=qrmean_all[1,6], color = "1.Original"), lwd=1) +
  geom_abline(aes(intercept=qrmean_all_inter[1,7], slope=qrmean_all[1,7]), linetype = 'blank', lwd=1) +
  geom_abline(aes(intercept=qrmean_all_inter[1,5], slope=qrmean_all[1,5]), linetype = 'blank', lwd=1) +
  geom_abline(aes(intercept=qrmean_all_inter[2,6], slope=qrmean_all[2,6],color = "2.UcS"), lwd=1) +
  geom_abline(aes(intercept=qrmean_all_inter[2,7], slope=qrmean_all[2,7]), linetype = 'blank', lwd=1) +
  geom_abline(aes(intercept=qrmean_all_inter[2,5], slope=qrmean_all[2,5]), linetype = 'blank', lwd=1) +
  geom_abline(aes(intercept=qrmean_all_inter[3,6], slope=qrmean_all[3,6], color = "3.Range size"), lwd=1) +
  geom_abline(aes(intercept=qrmean_all_inter[3,7], slope=qrmean_all[3,7]), linetype = 'blank', lwd=1) +
  geom_abline(aes(intercept=qrmean_all_inter[3,5], slope=qrmean_all[3,5]), linetype = 'blank', lwd=1) +
  geom_abline(aes(intercept=qrmean_all_inter[6,6], slope=qrmean_all[6,6], color = "4.Body size"), lwd=1) +
  geom_abline(aes(intercept=qrmean_all_inter[6,7], slope=qrmean_all[6,7]), linetype = 'blank', lwd=1) +
  geom_abline(aes(intercept=qrmean_all_inter[6,5], slope=qrmean_all[6,5]), linetype = 'blank', lwd=1) +
  geom_abline(aes(intercept=qrmean_all_inter[9,6], slope=qrmean_all[9,6],color = "5. Fossil"), lwd=1) +
  geom_abline(aes(intercept=qrmean_all_inter[9,7], slope=qrmean_all[9,7]), linetype = 'blank', lwd=1) +
  geom_abline(aes(intercept=qrmean_all_inter[9,5], slope=qrmean_all[9,5]), linetype = 'blank', lwd=1) +
  geom_abline(aes(intercept=qrmean_all_inter[12,6], slope=qrmean_all[12,6],color = "6. Taxonomy1"), lwd=1) +
  geom_abline(aes(intercept=qrmean_all_inter[12,7], slope=qrmean_all[12,7]), linetype = 'blank', lwd=1) +
  geom_abline(aes(intercept=qrmean_all_inter[12,5], slope=qrmean_all[12,5]), linetype = 'blank', lwd=1) +
  geom_abline(aes(intercept=qrmean_all_inter[15,6], slope=qrmean_all[15,6],color = "7. Taxonomy2"), lwd=1) +
  geom_abline(aes(intercept=qrmean_all_inter[15,7], slope=qrmean_all[15,7]), linetype = 'blank', lwd=1) +
  geom_abline(aes(intercept=qrmean_all_inter[15,5], slope=qrmean_all[15,5]), linetype = 'blank', lwd=1) +
  annotate(geom = "polygon", x = ribbons_plot[[1]]$x, y = ribbons_plot[[1]]$y, alpha = 0.15, 
           fill = rev(my_colours)[1]) +
  coord_cartesian(xlim = xlims_plot[[1]], ylim = ylims_plot[[1]]) +
  annotate(geom = "polygon", x = ribbons_plot[[2]]$x, y = ribbons_plot[[2]]$y, alpha = 0.15, 
           fill = rev(my_colours)[2]) +
  coord_cartesian(xlim = xlims_plot[[2]], ylim = ylims_plot[[2]]) +
  annotate(geom = "polygon", x = ribbons_plot[[3]]$x, y = ribbons_plot[[3]]$y, alpha = 0.15, 
           fill = rev(my_colours)[3]) +
  coord_cartesian(xlim = xlims_plot[[3]], ylim = ylims_plot[[3]]) +
  annotate(geom = "polygon", x = ribbons_plot[[6]]$x, y = ribbons_plot[[6]]$y, alpha = 0.15, 
           fill = rev(my_colours)[4]) +
  coord_cartesian(xlim = xlims_plot[[6]], ylim = ylims_plot[[6]]) +
  annotate(geom = "polygon", x = ribbons_plot[[9]]$x, y = ribbons_plot[[9]]$y, alpha = 0.15, 
           fill = rev(my_colours)[5]) +
  coord_cartesian(xlim = xlims_plot[[9]], ylim = ylims_plot[[9]]) +
  annotate(geom = "polygon", x = ribbons_plot[[12]]$x, y = ribbons_plot[[12]]$y, alpha = 0.15, 
           fill = rev(my_colours)[6]) +
  coord_cartesian(xlim = xlims_plot[[12]], ylim = ylims_plot[[12]]) +
  annotate(geom = "polygon", x = ribbons_plot[[15]]$x, y = ribbons_plot[[15]]$y, alpha = 0.15,
           fill = rev(my_colours)[7]) +
  coord_cartesian(xlim = xlims_plot[[15]], ylim = ylims_plot[[15]]) +
  scale_colour_manual(values = rev(my_colours)) +
  scale_x_continuous(breaks = c(0, 20, 40, 60, 80, 90), expand = expansion(mult = c(0, 0))) + 
  scale_y_continuous(breaks = c(0, 25, 50, 75, 100, 125, 150, 175, 200), expand = expansion(mult = c(0, 0)))

#pdf("09reg_75.pdf", useDingbats = F, width = 10, height = 7)
p
#dev.off()


p2 <- ggplot(mam_r_dfforplot, aes(abs(Longlat), Richness)) + geom_blank() +
  #xlim(0, 200) +
  xlab("") +
  #ylim(0,200) +
  ylab("") +
  theme_classic() +
  geom_abline(aes(intercept=qrmean_all_inter[1,6], slope=qrmean_all[1,6], color = "1.Original"), lwd=1) +
  geom_abline(aes(intercept=qrmean_all_inter[1,7], slope=qrmean_all[1,7]), linetype = 'blank', lwd=1) +
  geom_abline(aes(intercept=qrmean_all_inter[1,5], slope=qrmean_all[1,5]), linetype = 'blank', lwd=1) +
  geom_abline(aes(intercept=qrmean_all_inter[2,6], slope=qrmean_all[2,6],color = "2.UcS"), lwd=1) +
  geom_abline(aes(intercept=qrmean_all_inter[2,7], slope=qrmean_all[2,7]), linetype = 'blank', lwd=1) +
  geom_abline(aes(intercept=qrmean_all_inter[2,5], slope=qrmean_all[2,5]), linetype = 'blank', lwd=1) +
  geom_abline(aes(intercept=qrmean_all_inter[4,6], slope=qrmean_all[4,6], color = "3.Range size"), lwd=1) +
  geom_abline(aes(intercept=qrmean_all_inter[4,7], slope=qrmean_all[4,7]), linetype = 'blank', lwd=1) +
  geom_abline(aes(intercept=qrmean_all_inter[4,5], slope=qrmean_all[4,5]), linetype = 'blank', lwd=1) +
  geom_abline(aes(intercept=qrmean_all_inter[7,6], slope=qrmean_all[7,6], color = "4.Body size"), lwd=1) +
  geom_abline(aes(intercept=qrmean_all_inter[7,7], slope=qrmean_all[7,7]), linetype = 'blank', lwd=1) +
  geom_abline(aes(intercept=qrmean_all_inter[7,5], slope=qrmean_all[7,5]), linetype = 'blank', lwd=1) +
  geom_abline(aes(intercept=qrmean_all_inter[10,6], slope=qrmean_all[10,6],color = "5. Fossil"), lwd=1) +
  geom_abline(aes(intercept=qrmean_all_inter[10,7], slope=qrmean_all[10,7]), linetype = 'blank', lwd=1) +
  geom_abline(aes(intercept=qrmean_all_inter[10,5], slope=qrmean_all[10,5]), linetype = 'blank', lwd=1) +
  geom_abline(aes(intercept=qrmean_all_inter[13,6], slope=qrmean_all[13,6],color = "6. Taxonomy1"), lwd=1) +
  geom_abline(aes(intercept=qrmean_all_inter[13,7], slope=qrmean_all[13,7]), linetype = 'blank', lwd=1) +
  geom_abline(aes(intercept=qrmean_all_inter[13,5], slope=qrmean_all[13,5]), linetype = 'blank', lwd=1) +
  geom_abline(aes(intercept=qrmean_all_inter[16,6], slope=qrmean_all[16,6],color = "7. Taxonomy2"), lwd=1) +
  geom_abline(aes(intercept=qrmean_all_inter[16,7], slope=qrmean_all[16,7]), linetype = 'blank', lwd=1) +
  geom_abline(aes(intercept=qrmean_all_inter[16,5], slope=qrmean_all[16,5]), linetype = 'blank', lwd=1) +
  annotate(geom = "polygon", x = ribbons_plot[[1]]$x, y = ribbons_plot[[1]]$y, alpha = 0.15, 
           fill = rev(my_colours)[1]) +
  coord_cartesian(xlim = xlims_plot[[1]], ylim = ylims_plot[[1]]) +
  annotate(geom = "polygon", x = ribbons_plot[[2]]$x, y = ribbons_plot[[2]]$y, alpha = 0.15, 
           fill = rev(my_colours)[2]) +
  coord_cartesian(xlim = xlims_plot[[2]], ylim = ylims_plot[[2]]) +
  annotate(geom = "polygon", x = ribbons_plot[[4]]$x, y = ribbons_plot[[4]]$y, alpha = 0.15, 
           fill = rev(my_colours)[3]) +
  coord_cartesian(xlim = xlims_plot[[4]], ylim = ylims_plot[[4]]) +
  annotate(geom = "polygon", x = ribbons_plot[[7]]$x, y = ribbons_plot[[7]]$y, alpha = 0.15, 
           fill = rev(my_colours)[4]) +
  coord_cartesian(xlim = xlims_plot[[7]], ylim = ylims_plot[[7]]) +
  annotate(geom = "polygon", x = ribbons_plot[[10]]$x, y = ribbons_plot[[10]]$y, alpha = 0.15, 
           fill = rev(my_colours)[5]) +
  coord_cartesian(xlim = xlims_plot[[10]], ylim = ylims_plot[[10]]) +
  annotate(geom = "polygon", x = ribbons_plot[[13]]$x, y = ribbons_plot[[13]]$y, alpha = 0.15, 
           fill = rev(my_colours)[6]) +
  coord_cartesian(xlim = xlims_plot[[13]], ylim = ylims_plot[[13]]) +
  annotate(geom = "polygon", x = ribbons_plot[[16]]$x, y = ribbons_plot[[16]]$y, alpha = 0.15,
           fill = rev(my_colours)[7]) +
  coord_cartesian(xlim = xlims_plot[[16]], ylim = ylims_plot[[16]]) +
  scale_colour_manual(values = rev(my_colours)) +
  scale_x_continuous(breaks = c(0, 20, 40, 60, 80, 90), expand = expansion(mult = c(0, 0))) + 
  scale_y_continuous(breaks = c(0, 25, 50, 75, 100, 125, 150, 175, 200), expand = expansion(mult = c(0, 0)))

#pdf("09reg_50.pdf", useDingbats = F, width = 10, height = 7)
p2
#dev.off()


p3 <- ggplot(mam_r_dfforplot, aes(abs(Longlat), Richness)) + geom_blank() +
  #xlim(0, 200) +
  xlab("") +
  #ylim(0,200) +
  ylab("") +
  theme_classic() +
  geom_abline(aes(intercept=qrmean_all_inter[1,6], slope=qrmean_all[1,6], color = "1.Original"), lwd=1) +
  geom_abline(aes(intercept=qrmean_all_inter[1,7], slope=qrmean_all[1,7]), linetype = 'blank', lwd=1) +
  geom_abline(aes(intercept=qrmean_all_inter[1,5], slope=qrmean_all[1,5]), linetype = 'blank', lwd=1) +
  geom_abline(aes(intercept=qrmean_all_inter[2,6], slope=qrmean_all[2,6],color = "2.UcS"), lwd=1) +
  geom_abline(aes(intercept=qrmean_all_inter[2,7], slope=qrmean_all[2,7]), linetype = 'blank', lwd=1) +
  geom_abline(aes(intercept=qrmean_all_inter[2,5], slope=qrmean_all[2,5]), linetype = 'blank', lwd=1) +
  geom_abline(aes(intercept=qrmean_all_inter[5,6], slope=qrmean_all[5,6], color = "3.Range size"), lwd=1) +
  geom_abline(aes(intercept=qrmean_all_inter[5,7], slope=qrmean_all[5,7]), linetype = 'blank', lwd=1) +
  geom_abline(aes(intercept=qrmean_all_inter[5,5], slope=qrmean_all[5,5]), linetype = 'blank', lwd=1) +
  geom_abline(aes(intercept=qrmean_all_inter[8,6], slope=qrmean_all[8,6], color = "4.Body size"), lwd=1) +
  geom_abline(aes(intercept=qrmean_all_inter[8,7], slope=qrmean_all[8,7]), linetype = 'blank', lwd=1) +
  geom_abline(aes(intercept=qrmean_all_inter[8,5], slope=qrmean_all[8,5]), linetype = 'blank', lwd=1) +
  geom_abline(aes(intercept=qrmean_all_inter[11,6], slope=qrmean_all[11,6],color = "5. Fossil"), lwd=1) +
  geom_abline(aes(intercept=qrmean_all_inter[11,7], slope=qrmean_all[11,7]), linetype = 'blank', lwd=1) +
  geom_abline(aes(intercept=qrmean_all_inter[11,5], slope=qrmean_all[11,5]), linetype = 'blank', lwd=1) +
  geom_abline(aes(intercept=qrmean_all_inter[14,6], slope=qrmean_all[14,6],color = "6. Taxonomy1"), lwd=1) +
  geom_abline(aes(intercept=qrmean_all_inter[14,7], slope=qrmean_all[14,7]), linetype = 'blank', lwd=1) +
  geom_abline(aes(intercept=qrmean_all_inter[14,5], slope=qrmean_all[14,5]), linetype = 'blank', lwd=1) +
  geom_abline(aes(intercept=qrmean_all_inter[17,6], slope=qrmean_all[17,6],color = "7. Taxonomy2"), lwd=1) +
  geom_abline(aes(intercept=qrmean_all_inter[17,7], slope=qrmean_all[17,7]), linetype = 'blank', lwd=1) +
  geom_abline(aes(intercept=qrmean_all_inter[17,5], slope=qrmean_all[17,5]), linetype = 'blank', lwd=1) +
  annotate(geom = "polygon", x = ribbons_plot[[1]]$x, y = ribbons_plot[[1]]$y, alpha = 0.15, 
           fill = rev(my_colours)[1]) +
  coord_cartesian(xlim = xlims_plot[[1]], ylim = ylims_plot[[1]]) +
  annotate(geom = "polygon", x = ribbons_plot[[2]]$x, y = ribbons_plot[[2]]$y, alpha = 0.15, 
           fill = rev(my_colours)[2]) +
  coord_cartesian(xlim = xlims_plot[[2]], ylim = ylims_plot[[2]]) +
  annotate(geom = "polygon", x = ribbons_plot[[5]]$x, y = ribbons_plot[[5]]$y, alpha = 0.15, 
           fill = rev(my_colours)[3]) +
  coord_cartesian(xlim = xlims_plot[[5]], ylim = ylims_plot[[5]]) +
  annotate(geom = "polygon", x = ribbons_plot[[8]]$x, y = ribbons_plot[[8]]$y, alpha = 0.15, 
           fill = rev(my_colours)[4]) +
  coord_cartesian(xlim = xlims_plot[[8]], ylim = ylims_plot[[8]]) +
  annotate(geom = "polygon", x = ribbons_plot[[11]]$x, y = ribbons_plot[[11]]$y, alpha = 0.15, 
           fill = rev(my_colours)[5]) +
  coord_cartesian(xlim = xlims_plot[[11]], ylim = ylims_plot[[11]]) +
  annotate(geom = "polygon", x = ribbons_plot[[14]]$x, y = ribbons_plot[[14]]$y, alpha = 0.15, 
           fill = rev(my_colours)[6]) +
  coord_cartesian(xlim = xlims_plot[[14]], ylim = ylims_plot[[14]]) +
  annotate(geom = "polygon", x = ribbons_plot[[17]]$x, y = ribbons_plot[[17]]$y, alpha = 0.15,
           fill = rev(my_colours)[7]) +
  coord_cartesian(xlim = xlims_plot[[17]], ylim = ylims_plot[[17]]) +
  scale_colour_manual(values = rev(my_colours)) +
  scale_x_continuous(breaks = c(0, 20, 40, 60, 80, 90), expand = expansion(mult = c(0, 0))) + 
  scale_y_continuous(breaks = c(0, 25, 50, 75, 100, 125, 150, 175, 200), expand = expansion(mult = c(0, 0)))

#pdf("09reg_25.pdf", useDingbats = F, width = 10, height = 7)
p3
#dev.off()

#pdf("09reg_75_zoom.pdf", useDingbats = F, width = 10, height = 7)
p + coord_cartesian(xlim = c(0,5), ylim= c(160,180))
#dev.off()
#pdf("09reg_50_zoom.pdf", useDingbats = F, width = 10, height = 7)
p2 + coord_cartesian(xlim = c(0,5), ylim= c(160,180))
#dev.off()
#pdf("09reg_25_zoom.pdf", useDingbats = F, width = 10, height = 7)
p3 + coord_cartesian(xlim = c(0,5), ylim= c(160,180))
#dev.off()


#75%
sloss_75 <- c(1-(qrmean_all[2,6]/qrmean_all[1,6]),
              1-(qrmean_all[3,6]/qrmean_all[2,6]),
              1-(qrmean_all[6,6]/qrmean_all[3,6]),
              1-(qrmean_all[9,6]/qrmean_all[6,6]),
              1-(qrmean_all[15,6]/qrmean_all[9,6]),
              1-(qrmean_all[12,6]/qrmean_all[9,6]))
#50%
sloss_50 <- c(1-(qrmean_all[2,6]/qrmean_all[1,6]),
              1-(qrmean_all[4,6]/qrmean_all[2,6]),
              1-(qrmean_all[7,6]/qrmean_all[4,6]),
              1-(qrmean_all[10,6]/qrmean_all[7,6]),
              1-(qrmean_all[16,6]/qrmean_all[10,6]),
              1-(qrmean_all[13,6]/qrmean_all[10,6]))
#50%
sloss_25 <- c(1-(qrmean_all[2,6]/qrmean_all[1,6]),
              1-(qrmean_all[5,6]/qrmean_all[2,6]),
              1-(qrmean_all[8,6]/qrmean_all[5,6]),
              1-(qrmean_all[11,6]/qrmean_all[8,6]),
              1-(qrmean_all[17,6]/qrmean_all[11,6]),
              1-(qrmean_all[14,6]/qrmean_all[11,6]))
