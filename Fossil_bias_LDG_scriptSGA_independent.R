# ------------------------------------------------------------------------------
# Publication: Biases can't hide conspicuous biodiversity patterns 
# Last updated: 2024-07-11
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