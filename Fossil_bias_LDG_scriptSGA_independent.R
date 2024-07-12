# ------------------------------------------------------------------------------
# Publication: Biases can't hide conspicuous biodiversity patterns 
# Last updated: 2024-07-12
# Author: Sofia Galvan
# Email: sofia.galvan@uvigo.es
# Repository: https://github.com/SofiaGalv/Fossil_bias_and_LBG_detection.git
# ------------------------------------------------------------------------------

################################################################################
################### Independent bias-simulating workflow #######################
################################################################################

#Re-run L.1 - L.213 of "Fossil_bias_LDG_scriptSGA.R"

### Aplicar el filtro de range size sobre el original (mam_df_land)
set.seed(288)
library(scales)
rs25 <- sample(seq_len(nrow(Combi_r)), round(0.25*nrow(Combi_r)), prob = Combi_r$range_size)
Combi_25rs <- Combi_r[rs25,] #1425
rs50 <- sample(seq_len(nrow(Combi_r)), round(0.50*nrow(Combi_r)), prob = Combi_r$range_size)
Combi_50rs <- Combi_r[rs50,] #2850
rs75 <- sample(seq_len(nrow(Combi_r)), round(0.75*nrow(Combi_r)), prob = Combi_r$range_size)
Combi_75rs <- Combi_r[rs75,] #4274

Mam_25rs <- mam_df_land[,c("X", "Y",Combi_25rs$iucn2020_binomial, "Land_mass", "Clim_zone",
                           "GUM_USC", "WC_zone", "Richness", "Land", "Area", "Fossils")] #1435
Mam_50rs <- mam_df_land[,c("X", "Y",Combi_50rs$iucn2020_binomial, "Land_mass", "Clim_zone",
                           "GUM_USC", "WC_zone", "Richness", "Land", "Area","Fossils")] #2860
Mam_75rs <- mam_df_land[,c("X", "Y",Combi_75rs$iucn2020_binomial, "Land_mass", "Clim_zone",
                           "GUM_USC", "WC_zone", "Richness", "Land", "Area","Fossils")] #4284

identical(sort(colnames(Mam_25rs[,3:1427])), sort(Combi_25rs$iucn2020_binomial)) #T
identical(sort(colnames(Mam_50rs[,3:2852])), sort(Combi_50rs$iucn2020_binomial)) #T
identical(sort(colnames(Mam_75rs[,3:4276])), sort(Combi_75rs$iucn2020_binomial)) #T

richness <- apply(Mam_25rs[,3:1427], 1, function(x) length(x[x>0])) #rowSums(Mam_25s[,3:1189])
Mam_25rs[,1436] <- richness
colnames(Mam_25rs)[colnames(Mam_25rs) == "V1436"] <- "Richness2" #New column
Mam_25rs[Mam_25rs$Richness2 == 0, 1436] <- NA

richness <- apply(Mam_50rs[,3:2852], 1, function(x) length(x[x>0]))
Mam_50rs[,2861] <- richness
colnames(Mam_50rs)[colnames(Mam_50rs) == "V2861"] <- "Richness2" #New column
Mam_50rs[Mam_50rs$Richness2 == 0, 2861] <- NA

richness <- apply(Mam_75rs[,3:4276], 1, function(x) length(x[x>0]))
Mam_75rs[,4285] <- richness
colnames(Mam_75rs)[colnames(Mam_75rs) == "V4285"] <- "Richness2" #New column
Mam_75rs[Mam_75rs$Richness2 == 0, 4285] <- NA

#Reproject
Mam_25rs_fr <- Mam_25rs[,c("X", "Y", "Richness2")]
Mam_25rs_raster <- terra::rasterize(x = as.matrix(Mam_25rs_fr[,1:2]), y = as(x, "SpatRaster"), 
                                    values = Mam_25rs_fr[,3], update = T)
crs(Mam_25rs_raster) <- "+proj=longlat +datum=WGS84 +no_defs"

Mam_50rs_fr <- Mam_50rs[,c("X", "Y", "Richness2")]
Mam_50rs_raster <- terra::rasterize(x = as.matrix(Mam_50rs_fr[,1:2]), y = as(x, "SpatRaster"), 
                                    values = Mam_50rs_fr[,3], update = T)
crs(Mam_50rs_raster) <- "+proj=longlat +datum=WGS84 +no_defs"

Mam_75rs_fr <- Mam_75rs[,c("X", "Y", "Richness2")]
Mam_75rs_raster <- terra::rasterize(x = as.matrix(Mam_75rs_fr[,1:2]), y = as(x, "SpatRaster"), 
                                    values = Mam_75rs_fr[,3], update = T)
crs(Mam_75rs_raster) <- "+proj=longlat +datum=WGS84 +no_defs"

### Riqueza antes/despues
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
  if(any(unique(Mam_f5[,i]) == 1, na.rm = T)) {
    spp_number <- spp_number + 1
    print(spp_number)
  } 
}

spp_number <- 0
for (i in 3:4276) {
  if(any(unique(Mam_75rs_f1[,i]) == 1, na.rm = T)) {
    spp_number <- spp_number + 1
    print(spp_number)
  } 
}

spp_number <- 0
for (i in 3:2852) {
  if(any(unique(Mam_50rs_f1[,i]) == 1, na.rm = T)) {
    spp_number <- spp_number + 1
    print(spp_number)
  } 
}

spp_number <- 0
for (i in 3:1427) {
  if(any(unique(Mam_25rs_f1[,i]) == 1, na.rm = T)) {
    spp_number <- spp_number + 1
    print(spp_number)
  } 
}


### Aplicar el filtro de body size 
Combi_r2 <- Combi_r[!is.na(Combi_r$adult_mass_g),] 

set.seed(288)
bs75 <- sample(seq_len(nrow(Combi_r2)), round(0.75*nrow(Combi_r2)), prob = Combi_r2$adult_mass_g)
Combi_75bs <- Combi_r2[bs75,] #4124
Mam_75bs <- mam_df_land[,c("X", "Y",Combi_75bs$iucn2020_binomial, "Land_mass", "Clim_zone",
                           "GUM_USC", "WC_zone", "Richness", "Land", "Area", "Fossils")] #4134
identical(sort(colnames(Mam_75bs[,3:4126])), sort(Combi_75bs$iucn2020_binomial)) #T

richness <- apply(Mam_75bs[,3:4126], 1, function(x) length(x[x>0]))
Mam_75bs[,4135] <- richness
colnames(Mam_75bs)[colnames(Mam_75bs) == "V4135"] <- "Richness3" #New column
Mam_75bs[Mam_75bs$Richness3 == 0, 4135] <- NA


bs50 <- sample(seq_len(nrow(Combi_r2)), round(0.50*nrow(Combi_r2)), prob = Combi_r2$adult_mass_g)
Combi_50bs <- Combi_r2[bs50,] #2750
Mam_50bs <- mam_df_land[,c("X", "Y",Combi_50bs$iucn2020_binomial, "Land_mass", "Clim_zone",
                           "GUM_USC", "WC_zone", "Richness", "Land", "Area", "Fossils")] #2760
identical(sort(colnames(Mam_50bs[,3:2752])), sort(Combi_50bs$iucn2020_binomial)) #T

richness <- apply(Mam_50bs[,3:2752], 1, function(x) length(x[x>0]))
Mam_50bs[,2761] <- richness
colnames(Mam_50bs)[colnames(Mam_50bs) == "V2761"] <- "Richness3" #New column
Mam_50bs[Mam_50bs$Richness3 == 0, 2761] <- NA


bs25 <- sample(seq_len(nrow(Combi_r2)), round(0.25*nrow(Combi_r2)), prob = Combi_r2$adult_mass_g)
Combi_25bs <- Combi_r2[bs25,] #1375
Mam_25bs <- mam_df_land[,c("X", "Y",Combi_25bs$iucn2020_binomial, "Land_mass", "Clim_zone",
                           "GUM_USC", "WC_zone", "Richness", "Land", "Area", "Fossils")] #1385
identical(sort(colnames(Mam_25bs[,3:1377])), sort(Combi_25bs$iucn2020_binomial)) #T

richness <- apply(Mam_25bs[,3:1377], 1, function(x) length(x[x>0]))
Mam_25bs[,1386] <- richness
colnames(Mam_25bs)[colnames(Mam_25bs) == "V1386"] <- "Richness3" #New column
Mam_25bs[Mam_25bs$Richness3 == 0, 1386] <- NA

#Reproject
Mam_25bs_fr <- Mam_25bs[,c("X", "Y", "Richness3")]
Mam_25bs_raster <- terra::rasterize(x = as.matrix(Mam_25bs_fr[,1:2]), y = as(x, "SpatRaster"), 
                                    values = Mam_25bs_fr[,3], update = T)
crs(Mam_25bs_raster) <- "+proj=longlat +datum=WGS84 +no_defs"

Mam_50bs_fr <- Mam_50bs[,c("X", "Y", "Richness3")]
Mam_50bs_raster <- terra::rasterize(x = as.matrix(Mam_50bs_fr[,1:2]), y = as(x, "SpatRaster"), 
                                    values = Mam_50bs_fr[,3], update = T)
crs(Mam_50bs_raster) <- "+proj=longlat +datum=WGS84 +no_defs"

Mam_75bs_fr <- Mam_75bs[,c("X", "Y", "Richness3")]
Mam_75bs_raster <- terra::rasterize(x = as.matrix(Mam_75bs_fr[,1:2]), y = as(x, "SpatRaster"), 
                                    values = Mam_75bs_fr[,3], update = T)
crs(Mam_75bs_raster) <- "+proj=longlat +datum=WGS84 +no_defs"

#Riqueza antes/despues
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
for (i in 3:4126) {
  if(any(unique(Mam_75bs_f1[,i]) == 1, na.rm = T)) {
    spp_number <- spp_number + 1
    print(spp_number)
  } 
}

spp_number <- 0
for (i in 3:2752) {
  if(any(unique(Mam_50bs_f1[,i]) == 1, na.rm = T)) {
    spp_number <- spp_number + 1
    print(spp_number)
  } 
}

spp_number <- 0
for (i in 3:1377) {
  if(any(unique(Mam_25bs_f1[,i]) == 1, na.rm = T)) {
    spp_number <- spp_number + 1
    print(spp_number)
  } 
}


### Aplicar el filtro de fossil 
mam_df_land$WC_zone <- as.numeric(mam_df_land$WC_zone)
mam_df_land[is.na(mam_df_land$Fossils),5711] <- 0
mam_df_land[mam_df_land$Fossils >= 1,5711] <- 1
Mam_fos <- mam_df_land[mam_df_land$Fossils == 1,]
length(which(colSums(Mam_fos[,3:5703]) == 0)) #Se pierden 1267 (mas las 2 que se perdieron antes) 5701-1267 = 4434 spp

Mam_fos_df <- Mam_fos[,c("X", "Y", "Richness")]
Mam_fos_rast <- rasterize(x = as.matrix(Mam_fos_df[,1:2]), y = as(x, "SpatRaster"), 
                          values = Mam_fos_df[,3], update = T)
crs(Mam_fos_rast) <- "+proj=longlat +datum=WGS84 +no_defs"
Mam_fos_moll <- terra::project(Mam_fos_rast, 
                               y = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84",
                               method = "near", res = 50080)
Mamfos_moll_df <- terra::as.data.frame(Mam_fos_moll, xy = T, na.rm = F)
colnames(Mamfos_moll_df)[colnames(Mamfos_moll_df) == "last"] <- "Richness3"
Mamfos_moll_df[is.nan(Mamfos_moll_df$Richness3), 3] <- NA


#Riqueza por zona climatica
Mam_fos$WC_zone <- as.factor(Mam_fos$WC_zone)

Mam_fos_f1 <- Mam_fos[Mam_fos[,"WC_zone"] == 1,]
Mam_fos_f2 <- Mam_fos[Mam_fos[,"WC_zone"] == 2,]
Mam_fos_f3 <- Mam_fos[Mam_fos[,"WC_zone"] == 3,]
Mam_fos_f4 <- Mam_fos[Mam_fos[,"WC_zone"] == 4,]
Mam_fos_f5 <- Mam_fos[Mam_fos[,"WC_zone"] == 5,]

spp_number <- 0
for (i in 3:5703) {
  if(any(unique(Mam_fos_f2[,i]) == 1, na.rm = T)) {
    spp_number <- spp_number + 1
    print(spp_number)
  } 
}


### Aplicar el filtro de fossil taxonomico
Mam_pbdb[Mam_pbdb$genus == "Ptilocercus", "family"] <- "Ptilocercidae"
######################### 75%
Mam_fos_names <- names(which(colSums(mam_df_land[,3:5703]) != 0))
Combi_2 <- Combi_r[Combi_r$iucn2020_binomial %in% Mam_fos_names,] #5699 spp
spp_Combir <- sort(setdiff(unique(Combi_2$family), unique(Mam_pbdb$family)))
#Ahora hay 6 familias en Combi_75bs que no estan en Mam_pbdb (Ptilocercidae no tiene ocurrencias)
sort(table(Combi_2[Combi_2$family %in% spp_Combir, "family"])) #En total, 8 spp

spp_Mampbdb <- sort(setdiff(unique(Mam_pbdb$family), unique(Combi_2$family)))
#Hay 474 familias en Mam_pbdb que no estan en Combi_r

Combi_2 <- Combi_2[!Combi_2$family %in% spp_Combir,] #5691 obs, 133 familias, 8 spp menos
Mam_pbdb3 <- Mam_pbdb[!Mam_pbdb$family %in% spp_Mampbdb,] #67150 obs, 474 familias menos = 133 familias

Mam_pbdb3 <- as.data.frame(table(Mam_pbdb3$family))
colnames(Mam_pbdb3) <- c("family", "n_obs")
Mam_pbdb3$family <- as.character(Mam_pbdb3$family)

Combi_2$Fam_obs <- NA
family_names <- unique(Mam_pbdb3$family)
for (i in 1:length(family_names)) {
  Combi_2[Combi_2$family == family_names[[i]], "Fam_obs"] <- Mam_pbdb3[Mam_pbdb3$family == family_names[[i]], "n_obs"]
}

#Submuestreo "weighted"
set.seed(288)
tax75 <- sample(seq_len(nrow(Combi_2)), round(0.75*nrow(Combi_2)), prob = Combi_2$Fam_obs)
Combi_75tax <- Combi_2[tax75,] #4268
Mam_75tax <- mam_df_land[,c("X", "Y",Combi_75tax$iucn2020_binomial, "Land_mass", "Clim_zone",
                            "GUM_USC", "WC_zone", "Richness", "Land", "Area", "Fossils")] #4278
identical(sort(colnames(Mam_75tax[,3:4270])), sort(Combi_75tax$iucn2020_binomial)) #T

richness <- apply(Mam_75tax[,3:4270], 1, function(x) length(x[x>0]))
Mam_75tax[,4279] <- richness
colnames(Mam_75tax)[colnames(Mam_75tax) == "V4279"] <- "Richness2" #New column
Mam_75tax[Mam_75tax$Richness2 == 0, 4279] <- NA

Mam_75tax_fr <- Mam_75tax[,c("X", "Y", "Richness2")]
Mam_75tax_raster <- terra::rasterize(x = as.matrix(Mam_75tax_fr[,1:2]), y = as(x, "SpatRaster"), 
                                     values = Mam_75tax_fr[,3], update = T)
crs(Mam_75tax_raster) <- "+proj=longlat +datum=WGS84 +no_defs"

#Quedandonos con el 75% de las familias con mas observaciones
round(0.75*nrow(Mam_pbdb3)) #100 familias
sort(Mam_pbdb3$n_obs, decreasing = T)[[100]]
my_families <- Mam_pbdb3[Mam_pbdb3$n_obs >=17,1]
Combi_3 <- Combi_2[Combi_2$family %in% my_families,]

Mam_75tax2 <- mam_df_land[,c("X", "Y",Combi_3$iucn2020_binomial, "Land_mass", "Clim_zone",
                             "GUM_USC", "WC_zone", "Richness", "Land", "Area", "Fossils")] 
identical(sort(colnames(Mam_75tax2[,3:5400])), sort(Combi_3$iucn2020_binomial)) #T

richness <- apply(Mam_75tax2[,3:5400], 1, function(x) length(x[x>0]))
Mam_75tax2[,5409] <- richness
colnames(Mam_75tax2)[colnames(Mam_75tax2) == "V5409"] <- "Richness2" #New column
Mam_75tax2[Mam_75tax2$Richness2 == 0, 5409] <- NA

Mam_75tax2_fr <- Mam_75tax2[,c("X", "Y", "Richness2")]
Mam_75tax2_raster <- terra::rasterize(x = as.matrix(Mam_75tax2_fr[,1:2]), y = as(x, "SpatRaster"), 
                                      values = Mam_75tax2_fr[,3], update = T)
crs(Mam_75tax2_raster) <- "+proj=longlat +datum=WGS84 +no_defs"


####################### 50%
#Submuestreo "weighted"
set.seed(288)
tax50 <- sample(seq_len(nrow(Combi_2)), round(0.50*nrow(Combi_2)), prob = Combi_2$Fam_obs)
Combi_50tax <- Combi_2[tax50,] #2846
Mam_50tax <- mam_df_land[,c("X", "Y",Combi_50tax$iucn2020_binomial, "Land_mass", "Clim_zone",
                            "GUM_USC", "WC_zone", "Richness", "Land", "Area", "Fossils")] 
identical(sort(colnames(Mam_50tax[,3:2848])), sort(Combi_50tax$iucn2020_binomial)) #T

richness <- apply(Mam_50tax[,3:2848], 1, function(x) length(x[x>0]))
Mam_50tax[,2858] <- richness
colnames(Mam_50tax)[colnames(Mam_50tax) == "V2858"] <- "Richness2" #New column
Mam_50tax[Mam_50tax$Richness2 == 0, 2858] <- NA

Mam_50tax_fr <- Mam_50tax[,c("X", "Y", "Richness2")]
Mam_50tax_raster <- terra::rasterize(x = as.matrix(Mam_50tax_fr[,1:2]), y = as(x, "SpatRaster"), 
                                     values = Mam_50tax_fr[,3], update = T)
crs(Mam_50tax_raster) <- "+proj=longlat +datum=WGS84 +no_defs"


#Quedandonos con el 50% de las familias con mas observaciones
round(0.5*nrow(Mam_pbdb3)) #66 familias
sort(Mam_pbdb3$n_obs, decreasing = T)[[66]]
my_families <- Mam_pbdb3[Mam_pbdb3$n_obs >=82,1]
Combi_3 <- Combi_2[Combi_2$family %in% my_families,]

Mam_50tax2 <- mam_df_land[,c("X", "Y",Combi_3$iucn2020_binomial, "Land_mass", "Clim_zone",
                             "GUM_USC", "WC_zone", "Richness", "Land", "Area", "Fossils")] 
identical(sort(colnames(Mam_50tax2[,3:4669])), sort(Combi_3$iucn2020_binomial)) #T

richness <- apply(Mam_50tax2[,3:4669], 1, function(x) length(x[x>0]))
Mam_50tax2[,4678] <- richness
colnames(Mam_50tax2)[colnames(Mam_50tax2) == "V4678"] <- "Richness2" #New column
Mam_50tax2[Mam_50tax2$Richness2 == 0, 4678] <- NA


Mam_50tax2_fr <- Mam_50tax2[,c("X", "Y", "Richness2")]
Mam_50tax2_raster <- terra::rasterize(x = as.matrix(Mam_50tax2_fr[,1:2]), y = as(x, "SpatRaster"), 
                                      values = Mam_50tax2_fr[,3], update = T)
crs(Mam_50tax2_raster) <- "+proj=longlat +datum=WGS84 +no_defs"



####################### 25%
#Submuestreo "weighted"
set.seed(288)
tax25 <- sample(seq_len(nrow(Combi_2)), round(0.25*nrow(Combi_2)), prob = Combi_2$Fam_obs)
Combi_25tax <- Combi_2[tax25,] #1423
Mam_25tax <- mam_df_land[,c("X", "Y",Combi_25tax$iucn2020_binomial, "Land_mass", "Clim_zone",
                            "GUM_USC", "WC_zone", "Richness", "Land", "Area", "Fossils")] 
identical(sort(colnames(Mam_25tax[,3:1425])), sort(Combi_25tax$iucn2020_binomial)) #T

richness <- apply(Mam_25tax[,3:1425], 1, function(x) length(x[x>0]))
Mam_25tax[,1434] <- richness
colnames(Mam_25tax)[colnames(Mam_25tax) == "V1434"] <- "Richness2" #New column
Mam_25tax[Mam_25tax$Richness2 == 0, 1434] <- NA

Mam_25tax_fr <- Mam_25tax[,c("X", "Y", "Richness2")]
Mam_25tax_raster <- terra::rasterize(x = as.matrix(Mam_25tax_fr[,1:2]), y = as(x, "SpatRaster"), 
                                     values = Mam_25tax_fr[,3], update = T)
crs(Mam_25tax_raster) <- "+proj=longlat +datum=WGS84 +no_defs"


#Quedandonos con el 25% de las familias con mas observaciones
round(0.25*nrow(Mam_pbdb3)) #33 familias
sort(Mam_pbdb3$n_obs, decreasing = T)[[33]]
my_families <- Mam_pbdb3[Mam_pbdb3$n_obs >=504,1]
Combi_3 <- Combi_2[Combi_2$family %in% my_families,]

Mam_25tax2 <- mam_df_land[,c("X", "Y",Combi_3$iucn2020_binomial, "Land_mass", "Clim_zone",
                             "GUM_USC", "WC_zone", "Richness", "Land", "Area", "Fossils")] 
identical(sort(colnames(Mam_25tax2[,3:3466])), sort(Combi_3$iucn2020_binomial)) #T

richness <- apply(Mam_25tax2[,3:3466], 1, function(x) length(x[x>0]))
Mam_25tax2[,3475] <- richness
colnames(Mam_25tax2)[colnames(Mam_25tax2) == "V3475"] <- "Richness2" #New column
Mam_25tax2[Mam_25tax2$Richness2 == 0, 3475] <- NA

Mam_25tax2_fr <- Mam_25tax2[,c("X", "Y", "Richness2")]
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
for (i in 3:4270) {
  if(any(unique(Mam_75tax_f5[,i]) == 1, na.rm = T)) {
    spp_number <- spp_number + 1
    print(spp_number)
  } 
}

spp_number <- 0
for (i in 3:2847) {
  if(any(unique(Mam_50tax_f5[,i]) == 1, na.rm = T)) {
    spp_number <- spp_number + 1
    print(spp_number)
  } 
}

spp_number <- 0
for (i in 3:1424) {
  if(any(unique(Mam_25tax_f5[,i]) == 1, na.rm = T)) {
    spp_number <- spp_number + 1
    print(spp_number)
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
for (i in 3:5400) {
  if(any(unique(Mam_75tax2_f5[,i]) == 1, na.rm = T)) {
    spp_number <- spp_number + 1
    print(spp_number)
  } 
}

spp_number <- 0
for (i in 3:4669) {
  if(any(unique(Mam_50tax2_f5[,i]) == 1, na.rm = T)) {
    spp_number <- spp_number + 1
    print(spp_number)
  } 
}

spp_number <- 0
for (i in 3:3466) {
  if(any(unique(Mam_25tax2_f5[,i]) == 1, na.rm = T)) {
    spp_number <- spp_number + 1
    print(spp_number)
  } 
}


################################# Richness maps ################################
#setwd("./Figures/Figure_Indepcase/")
rast_bc <- rast(nrows = 360, ncols = 720, resolution = 0.5, vals = 1,
                extent = ext(-180, 180, -90, 90))
rast_bc <- terra::project(rast_bc, 
                          y = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84",
                          method = "near",
                          res = 50080)
sp_bc <- as.polygons(rast_bc)
#pdf("Indep_Mam_25taxwholefam.pdf", useDingbats = F)
plot(sp_bc, axes = F)
plot(land_moll$geometry, col = "grey30", border = "grey30", lwd=0.4, add = T) #, bg = "skyblue"
terra::plot(Mam_25tax2_moll, col = rev(paletteer_c("grDevices::YlOrBr", 400)), axes = F, 
            range = c(0,250), pax = list(labels = F, tick = F), legend = F, add = T)
plot(land_moll$geometry, border = "grey30", add = T, lwd=0.4)
#dev.off()

################################# Area loss ####################################
wc_moll <- terra::project(wc_raster, 
                          y = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84",
                          method = "near", res = 50080)
wc_molldf <- terra::as.data.frame(wc_moll, xy = T, na.rm = F)
wc_molldf[is.nan(wc_molldf$lyr.1), 3] <- NA
npixel_original <- table(wc_molldf$lyr.1)


wc_fos_extracted <- wc_extracted <- terra::mask(wc_moll, Mam_fos_moll)
wc_fos_df <- terra::as.data.frame(wc_fos_extracted, xy = T, na.rm = F)
npixel_fos <- table(wc_fos_df$lyr.1)

pixels_df2 <- data.frame(npixel_original, npixel_fos)[,c(1,2,4)]
colnames(pixels_df2) <- c("WC_zone", "npixel_raw", "npixel_fos")
pixels_df2$per_loss <- 100-(round((pixels_df2$npixel_fos/pixels_df2$npixel_raw)*100,2))
pixels_df2$UCS_sed <- round((pixels_df2$npixel_fos/pixels_df2$npixel_raw)*100,2)

#pdf("Indep_CZ_Humanspatial.pdf", useDingbats = F)
plot(sp_bc, axes = F)
plot(land_moll$geometry, col = "grey30", border = "grey30", lwd=0.4, add = T) #, bg = "skyblue"
terra::plot(wc_fos_extracted, col = c("#f9d14a", "#ab3329", "#ed968c", "#7c4b73", "#88a0dc"),
            axes = F, add = T)
plot(land_moll$geometry, border = "grey30", add = T, lwd=0.4)
#dev.off()

##########################  Richness loss estimates ############################
Mam_rich_total <- read.csv("Filters_alone_richness.csv", header = T, sep = ";")
for (i in 3:17) {
  Mam_rich_total[,i] <- as.numeric(Mam_rich_total[,i])
}

Mam_rich_total[,18] <- round((1-(Mam_rich_total$UcS/Mam_rich_total$Original))*100, 2)
colnames(Mam_rich_total)[colnames(Mam_rich_total) == "V18"] <- "Orig_UcS"

Mam_rich_total[,19] <- round((1-(Mam_rich_total$Rangesize_75/Mam_rich_total$Original))*100, 2)
colnames(Mam_rich_total)[colnames(Mam_rich_total) == "V19"] <- "Orig_rs75"
Mam_rich_total[,20] <- round((1-(Mam_rich_total$Rangesize_50/Mam_rich_total$Original))*100, 2)
colnames(Mam_rich_total)[colnames(Mam_rich_total) == "V20"] <- "Orig_rs50"
Mam_rich_total[,21] <- round((1-(Mam_rich_total$Rangesize_25/Mam_rich_total$Original))*100, 2)
colnames(Mam_rich_total)[colnames(Mam_rich_total) == "V21"] <- "Orig_rs25"

Mam_rich_total[,22] <- round((1-(Mam_rich_total$Bodysize_75/Mam_rich_total$Original))*100, 2)
colnames(Mam_rich_total)[colnames(Mam_rich_total) == "V22"] <- "Orig_bs75"
Mam_rich_total[,23] <- round((1-(Mam_rich_total$Bodysize_50/Mam_rich_total$Original))*100, 2)
colnames(Mam_rich_total)[colnames(Mam_rich_total) == "V23"] <- "Orig_bs50"
Mam_rich_total[,24] <- round((1-(Mam_rich_total$Bodysize_25/Mam_rich_total$Original))*100, 2)
colnames(Mam_rich_total)[colnames(Mam_rich_total) == "V24"] <- "Orig_bs25"

Mam_rich_total[,25] <- round((1-(Mam_rich_total$Fossil/Mam_rich_total$Original))*100, 2)
colnames(Mam_rich_total)[colnames(Mam_rich_total) == "V25"] <- "Orig_fos7"

Mam_rich_total[,26] <- round((1-(Mam_rich_total$Taxweight_75/Mam_rich_total$Original))*100, 2)
colnames(Mam_rich_total)[colnames(Mam_rich_total) == "V26"] <- "Orig_taxw75"
Mam_rich_total[,27] <- round((1-(Mam_rich_total$Taxweight_50/Mam_rich_total$Original))*100, 2)
colnames(Mam_rich_total)[colnames(Mam_rich_total) == "V27"] <- "Orig_taxw50"
Mam_rich_total[,28] <- round((1-(Mam_rich_total$Taxweight_25/Mam_rich_total$Original))*100, 2)
colnames(Mam_rich_total)[colnames(Mam_rich_total) == "V28"] <- "Orig_taxw25"

Mam_rich_total[,29] <- round((1-(Mam_rich_total$Taxfamily_75/Mam_rich_total$Original))*100, 2)
colnames(Mam_rich_total)[colnames(Mam_rich_total) == "V29"] <- "Orig_taxf75"
Mam_rich_total[,30] <- round((1-(Mam_rich_total$Taxfamily_50/Mam_rich_total$Original))*100, 2)
colnames(Mam_rich_total)[colnames(Mam_rich_total) == "V30"] <- "Orig_taxf50"
Mam_rich_total[,31] <- round((1-(Mam_rich_total$Taxfamily_25/Mam_rich_total$Original))*100, 2)
colnames(Mam_rich_total)[colnames(Mam_rich_total) == "V31"] <- "Orig_tafw25"

var_toplot <- colnames(Mam_rich_total[,18:31])
labels_toplot <- list(c("Tropical (21.8%)", "Arid (14.0%)", "Temperate (27.0%)", "Cold (22.9%)", "Polar (42.4%)"),
                      c("Tropical (24.3%)", "Arid (9.7%)", "Temperate (17.1%)", "Cold (8.7%)", "Polar (6.2%)"),
                      c("Tropical (49.5%)", "Arid (28.0%)", "Temperate (40.1%)", "Cold (25.6%)", "Polar (21.1%)"),
                      c("Tropical (74.5%)", "Arid (59.0%)", "Temperate (67.3%)", "Cold (55.2%)", "Polar (47.4%)"),
                      c("Tropical (27.5%)", "Arid (28.7%)", "Temperate (28.6%)", "Cold (27.7%)", "Polar (27.6%)"),
                      c("Tropical (50.8%)", "Arid (52.6%)", "Temperate (52.5%)", "Cold (51.5%)", "Polar (51.7%)"),
                      c("Tropical (74.4%)", "Arid (75.5%)", "Temperate (76.1%)", "Cold (75.1%)", "Polar (73.8%)"),
                      c("Tropical (30.9%)", "Arid (19.1%)", "Temperate (23.3%)", "Cold (20.7%)", "Polar (47.5%)"),
                      c("Tropical (29.9%)", "Arid (21.5%)", "Temperate (23.3%)", "Cold (9.2%)", "Polar (16.3%)"),
                      c("Tropical (54.9%)", "Arid (47.6%)", "Temperate (48.4%)", "Cold (33.9%)", "Polar (41.0%)"),
                      c("Tropical (78.3%)", "Arid (73.1%)", "Temperate (73.4%)", "Cold (67.1%)", "Polar (69.3%)"),
                      c("Tropical (6.7%)", "Arid (3.0%)", "Temperate (4.1%)", "Cold (1.1%)", "Polar (1.2%)"),
                      c("Tropical (22.9%)", "Arid (11.9%)", "Temperate (14.8%)", "Cold (3.2%)", "Polar (5.9%)"),
                      c("Tropical (45.5%)", "Arid (36.0%)", "Temperate (37.4%)", "Cold (16.3%)", "Polar (29.4%)"))
names(labels_toplot) <- colnames(Mam_rich_total[,18:31])


richloss_plots <- list()
for (i in 1:14) {
  mam_df <- Mam_rich_total[,var_toplot[i]][2:6]
  df_forplot <- data.frame(mam_df, labels_toplot[[i]])
  colnames(df_forplot) <- c("value1", "labels_names")
  df_forplot$labels_names <- factor(df_forplot$labels_names, levels = df_forplot$labels_names)
  pl <- ggplot(df_forplot, #El reordering no funciona
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

richloss_plots[[14]]

#setwd("./Figures/Figure_Indepcase")
#pdf("Indep_Orig_taxf25_spploss.pdf", useDingbats = F)
richloss_plots[[14]]
#dev.off()


################## Correlation between area/richness loss Original-UcS
df_forcor2 <- as.data.frame(cbind(pixels_df2$WC_zone, pixels_df2$per_loss, Mam_rich_total$Orig_fos7[2:6]))
colnames(df_forcor2) <- c("WC_zone", "Area_loss", "Sp_loss")
df_forcor2$WC_zone <- as.factor(df_forcor2$WC_zone)

#pdf("Indep_Cor_OrigHumanSpatial.pdf", useDingbats = F, width = 10, height = 7)
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


############################ 0.9 quantile regression ###########################
### Slopes
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
Mam_fos_moll <- terra::project(Mam_fos_rast, 
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
Mamfos_moll_df <- terra::as.data.frame(Mam_fos_moll, xy = T, na.rm = F)
Mamfos_moll_df[is.nan(Mamfos_moll_df$lyr.1), 3] <- NA
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
coordinates(Mamfos_moll_df) <- c("x", "y")
proj4string(Mamfos_moll_df) <- CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84")
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
Mamfos_moll_sf <- st_as_sf(Mamfos_moll_df)
Mamfos_longlat_sf <- st_transform(Mamfos_moll_sf, "+proj=longlat +datum=WGS84 +no_defs")
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
mam_r_dfforplot <- read.table("sf_coords.csv")

db_forslopes <- data.frame(mam_r_dfforplot$Longlat, mam_r_dfforplot$Moll,
                           mam_r_dfforplot$Richness, mam_f_moll_sf$last, Mam75rs_moll_sf$last,
                           Mam50rs_moll_sf$last, Mam25rs_moll_sf$last, Mam75bs_moll_sf$last, Mam50bs_moll_sf$last,
                           Mam25bs_moll_sf$last, Mamfos_moll_sf$last, Mam_75tax_moll_sf$last, Mam_50tax_moll_sf$last,
                           Mam_25tax_moll_sf$last, Mam_75tax2_moll_sf$last, Mam_50tax2_moll_sf$last,
                           Mam_25tax2_moll_sf$last)
colnames(db_forslopes) <- c("Longlat", "Moll", "Original", "UcS", "rs75", "rs50",
                            "rs25", "bs75", "bs50", "bs25", "foss",
                            "tax75", "tax50", "tax25", "tax2_75", 
                            "tax2_50", "tax2_25")
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

vari_names <- colnames(db_forslopes)[3:17]
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


qr_means_inter <- function(ddbb_input) {
  medias <- c()
  CI <- c()
  mins <- c()
  maxs <- c()
  impares <- seq(1, 200, by = 2)
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


### Loop para aplicarlo a todas las variables
ribbons_plot <- list()
xlims_plot <- list()
ylims_plot <- list()
for (i in 1:15) {
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
  geom_abline(aes(intercept=qrmean_all_inter[10,6], slope=qrmean_all[10,6],color = "6. Taxonomy1"), lwd=1) +
  geom_abline(aes(intercept=qrmean_all_inter[10,7], slope=qrmean_all[10,7]), linetype = 'blank', lwd=1) +
  geom_abline(aes(intercept=qrmean_all_inter[10,5], slope=qrmean_all[10,5]), linetype = 'blank', lwd=1) +
  geom_abline(aes(intercept=qrmean_all_inter[13,6], slope=qrmean_all[13,6],color = "7. Taxonomy2"), lwd=1) +
  geom_abline(aes(intercept=qrmean_all_inter[13,7], slope=qrmean_all[13,7]), linetype = 'blank', lwd=1) +
  geom_abline(aes(intercept=qrmean_all_inter[13,5], slope=qrmean_all[13,5]), linetype = 'blank', lwd=1) +
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
  annotate(geom = "polygon", x = ribbons_plot[[10]]$x, y = ribbons_plot[[10]]$y, alpha = 0.15, 
           fill = rev(my_colours)[6]) +
  coord_cartesian(xlim = xlims_plot[[10]], ylim = ylims_plot[[10]]) +
  annotate(geom = "polygon", x = ribbons_plot[[13]]$x, y = ribbons_plot[[13]]$y, alpha = 0.15,
           fill = rev(my_colours)[7]) +
  coord_cartesian(xlim = xlims_plot[[13]], ylim = ylims_plot[[13]]) +
  scale_colour_manual(values = rev(my_colours)) +
  scale_x_continuous(breaks = c(0, 20, 40, 60, 80, 90), expand = expansion(mult = c(0, 0))) + 
  scale_y_continuous(breaks = c(0, 25, 50, 75, 100, 125, 150, 175, 200), expand = expansion(mult = c(0, 0)))

#pdf("09reg_75ind.pdf", useDingbats = F, width = 10, height = 7)
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
  geom_abline(aes(intercept=qrmean_all_inter[9,6], slope=qrmean_all[9,6],color = "5. Fossil"), lwd=1) +
  geom_abline(aes(intercept=qrmean_all_inter[9,7], slope=qrmean_all[9,7]), linetype = 'blank', lwd=1) +
  geom_abline(aes(intercept=qrmean_all_inter[9,5], slope=qrmean_all[9,5]), linetype = 'blank', lwd=1) +
  geom_abline(aes(intercept=qrmean_all_inter[11,6], slope=qrmean_all[11,6],color = "6. Taxonomy1"), lwd=1) +
  geom_abline(aes(intercept=qrmean_all_inter[11,7], slope=qrmean_all[11,7]), linetype = 'blank', lwd=1) +
  geom_abline(aes(intercept=qrmean_all_inter[11,5], slope=qrmean_all[11,5]), linetype = 'blank', lwd=1) +
  geom_abline(aes(intercept=qrmean_all_inter[14,6], slope=qrmean_all[14,6],color = "7. Taxonomy2"), lwd=1) +
  geom_abline(aes(intercept=qrmean_all_inter[14,7], slope=qrmean_all[14,7]), linetype = 'blank', lwd=1) +
  geom_abline(aes(intercept=qrmean_all_inter[14,5], slope=qrmean_all[14,5]), linetype = 'blank', lwd=1) +
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
  annotate(geom = "polygon", x = ribbons_plot[[9]]$x, y = ribbons_plot[[9]]$y, alpha = 0.15, 
           fill = rev(my_colours)[5]) +
  coord_cartesian(xlim = xlims_plot[[9]], ylim = ylims_plot[[9]]) +
  annotate(geom = "polygon", x = ribbons_plot[[11]]$x, y = ribbons_plot[[11]]$y, alpha = 0.15, 
           fill = rev(my_colours)[6]) +
  coord_cartesian(xlim = xlims_plot[[11]], ylim = ylims_plot[[11]]) +
  annotate(geom = "polygon", x = ribbons_plot[[14]]$x, y = ribbons_plot[[14]]$y, alpha = 0.15,
           fill = rev(my_colours)[7]) +
  coord_cartesian(xlim = xlims_plot[[14]], ylim = ylims_plot[[14]]) +
  scale_colour_manual(values = rev(my_colours)) +
  scale_x_continuous(breaks = c(0, 20, 40, 60, 80, 90), expand = expansion(mult = c(0, 0))) + 
  scale_y_continuous(breaks = c(0, 25, 50, 75, 100, 125, 150, 175, 200), expand = expansion(mult = c(0, 0)))

#pdf("09reg_50ind.pdf", useDingbats = F, width = 10, height = 7)
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
  annotate(geom = "polygon", x = ribbons_plot[[5]]$x, y = ribbons_plot[[5]]$y, alpha = 0.15, 
           fill = rev(my_colours)[3]) +
  coord_cartesian(xlim = xlims_plot[[5]], ylim = ylims_plot[[5]]) +
  annotate(geom = "polygon", x = ribbons_plot[[8]]$x, y = ribbons_plot[[8]]$y, alpha = 0.15, 
           fill = rev(my_colours)[4]) +
  coord_cartesian(xlim = xlims_plot[[8]], ylim = ylims_plot[[8]]) +
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

#pdf("09reg_25ind.pdf", useDingbats = F, width = 10, height = 7)
p3
#dev.off()

#pdf("09reg_75_zoomind.pdf", useDingbats = F, width = 10, height = 7)
p + coord_cartesian(xlim = c(0,5), ylim= c(160,180))
#dev.off()
#pdf("09reg_50_zoomind.pdf", useDingbats = F, width = 10, height = 7)
p2 + coord_cartesian(xlim = c(0,5), ylim= c(160,180))
#dev.off()
#pdf("09reg_25_zoomind.pdf", useDingbats = F, width = 10, height = 7)
p3 + coord_cartesian(xlim = c(0,5), ylim= c(160,180))
#dev.off()

#75%
sloss_75 <- c(1-(qrmean_all[2,6]/qrmean_all[1,6]),
              1-(qrmean_all[3,6]/qrmean_all[1,6]),
              1-(qrmean_all[6,6]/qrmean_all[1,6]),
              1-(qrmean_all[9,6]/qrmean_all[1,6]),
              1-(qrmean_all[10,6]/qrmean_all[1,6]),
              1-(qrmean_all[13,6]/qrmean_all[1,6]))
#50%
sloss_50 <- c(1-(qrmean_all[2,6]/qrmean_all[1,6]),
              1-(qrmean_all[4,6]/qrmean_all[1,6]),
              1-(qrmean_all[7,6]/qrmean_all[1,6]),
              1-(qrmean_all[9,6]/qrmean_all[1,6]),
              1-(qrmean_all[11,6]/qrmean_all[1,6]),
              1-(qrmean_all[14,6]/qrmean_all[1,6]))
#50%
sloss_25 <- c(1-(qrmean_all[2,6]/qrmean_all[1,6]),
              1-(qrmean_all[5,6]/qrmean_all[1,6]),
              1-(qrmean_all[8,6]/qrmean_all[1,6]),
              1-(qrmean_all[9,6]/qrmean_all[1,6]),
              1-(qrmean_all[12,6]/qrmean_all[1,6]),
              1-(qrmean_all[15,6]/qrmean_all[1,6]))
