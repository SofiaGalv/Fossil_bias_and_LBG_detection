# ------------------------------------------------------------------------------
# Publication: Biases can't hide conspicuous biodiversity patterns 
# Last updated: 2024-07-11
# Author: Sofia Galvan
# Email: sofia.galvan@uvigo.es
# Repository: https://github.com/SofiaGalv/Fossil_bias_and_LBG_detection.git
# ------------------------------------------------------------------------------

################################################################################
####################### Workflow with 9 replicates #############################
################################################################################

set.seed(13)
seeds <- sample(1:2000, 9, replace = F)

#Re-run L.1 - L.213 of "Fossil_bias_LDG_scriptSGA.R"

################################# seed = 1496 ##################################
identical(sort(Mam_fil_names), sort(Combi_4$iucn2020_binomial))

### Range size filter

set.seed(1496)
library(scales)
rs25 <- sample(seq_len(nrow(Combi_4)), round(0.25*nrow(Combi_4)), prob = Combi_4$range_size)
Combi_25rs <- Combi_4[rs25,] #1186
rs50 <- sample(seq_len(nrow(Combi_4)), round(0.50*nrow(Combi_4)), prob = Combi_4$range_size)
Combi_50rs <- Combi_4[rs50,] #2372
rs75 <- sample(seq_len(nrow(Combi_4)), round(0.75*nrow(Combi_4)), prob = Combi_4$range_size)
Combi_75rs <- Combi_4[rs75,] #3559

Mam_25rs <- Ric_mam_fil[,c("X", "Y",Combi_25rs$iucn2020_binomial, "Land_mass", "Clim_zone",
                           "GUM_USC", "WC_zone", "Richness", "Land", "Area", "Fossils")] #1196
Mam_50rs <- Ric_mam_fil[,c("X", "Y",Combi_50rs$iucn2020_binomial, "Land_mass", "Clim_zone",
                           "GUM_USC", "WC_zone", "Richness", "Land", "Area","Fossils")] #2382
Mam_75rs <- Ric_mam_fil[,c("X", "Y",Combi_75rs$iucn2020_binomial, "Land_mass", "Clim_zone",
                           "GUM_USC", "WC_zone", "Richness", "Land", "Area","Fossils")] #3569

identical(sort(colnames(Mam_25rs[,3:1188])), sort(Combi_25rs$iucn2020_binomial)) #T
identical(sort(colnames(Mam_50rs[,3:2374])), sort(Combi_50rs$iucn2020_binomial)) #T
identical(sort(colnames(Mam_75rs[,3:3561])), sort(Combi_75rs$iucn2020_binomial)) #T

richness <- apply(Mam_25rs[,3:1188], 1, function(x) length(x[x>0])) #rowSums(Mam_25s[,3:1189])
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
for (i in 3:(ncol(Mam_25rs_f5)-10)) {
  if(any(unique(Mam_25rs_f5[,i]) == 1, na.rm = T)) {
    spp_number <- spp_number + 1
    print(spp_number)
  } 
}

#Store it in "richness_replicates.xlsx"

### Body size filter
Combi_5 <- Combi_75rs[!is.na(Combi_75rs$adult_mass_g),] #3516

set.seed(1496)
bs75 <- sample(seq_len(nrow(Combi_5)), round(0.75*nrow(Combi_5)), prob = Combi_5$adult_mass_g)
Combi_75bs <- Combi_5[bs75,] #2637
Mam_75bs <- Mam_75rs[,c("X", "Y",Combi_75bs$iucn2020_binomial, "Land_mass", "Clim_zone",
                        "GUM_USC", "WC_zone", "Richness", "Land", "Area", "Fossils", "Richness2")]#2649  
identical(sort(colnames(Mam_75bs[,3:2639])), sort(Combi_75bs$iucn2020_binomial)) #T

richness <- apply(Mam_75bs[,3:2639], 1, function(x) length(x[x>0]))
Mam_75bs[,2649] <- richness
colnames(Mam_75bs)[colnames(Mam_75bs) == "V2649"] <- "Richness3" #New column
Mam_75bs[Mam_75bs$Richness3 == 0, 2649] <- NA


Combi_5 <- Combi_50rs[!is.na(Combi_50rs$adult_mass_g),] #2347
bs50 <- sample(seq_len(nrow(Combi_5)), round(0.50*nrow(Combi_5)), prob = Combi_5$adult_mass_g)
Combi_50bs <- Combi_5[bs50,] #1174
Mam_50bs <- Mam_50rs[,c("X", "Y",Combi_50bs$iucn2020_binomial, "Land_mass", "Clim_zone",
                        "GUM_USC", "WC_zone", "Richness", "Land", "Area", "Fossils", "Richness2")]#1186 
identical(sort(colnames(Mam_50bs[,3:1176])), sort(Combi_50bs$iucn2020_binomial)) #T

richness <- apply(Mam_50bs[,3:1176], 1, function(x) length(x[x>0]))
Mam_50bs[,1186] <- richness
colnames(Mam_50bs)[colnames(Mam_50bs) == "V1186"] <- "Richness3" #New column
Mam_50bs[Mam_50bs$Richness3 == 0, 1186] <- NA


Combi_5 <- Combi_25rs[!is.na(Combi_25rs$adult_mass_g),] #1182
bs25 <- sample(seq_len(nrow(Combi_5)), round(0.25*nrow(Combi_5)), prob = Combi_5$adult_mass_g)
Combi_25bs <- Combi_5[bs25,] #296
Mam_25bs <- Mam_25rs[,c("X", "Y",Combi_25bs$iucn2020_binomial, "Land_mass", "Clim_zone",
                        "GUM_USC", "WC_zone", "Richness", "Land", "Area", "Fossils", "Richness2")]#308 
identical(sort(colnames(Mam_25bs[,3:298])), sort(Combi_25bs$iucn2020_binomial)) #T

richness <- apply(Mam_25bs[,3:298], 1, function(x) length(x[x>0]))
Mam_25bs[,308] <- richness
colnames(Mam_25bs)[colnames(Mam_25bs) == "V308"] <- "Richness3" #New column
Mam_25bs[Mam_25bs$Richness3 == 0, 308] <- NA

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

### Riqueza antes/despues
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
for (i in 3:(ncol(Mam_25bs_f5)-11)) {
  if(any(unique(Mam_25bs_f5[,i]) == 1, na.rm = T)) {
    spp_number <- spp_number + 1
    print(spp_number)
  } 
}

### Fossil sites filter
Mam_75bs$WC_zone <- as.numeric(Mam_75bs$WC_zone)
Mam_75bs[is.na(Mam_75bs$Fossils),2647] <- 0
Mam_75bs[Mam_75bs$Fossils >= 1,2647] <- 1
Mam_75fos <- Mam_75bs[Mam_75bs$Fossils == 1,]
length(which(colSums(Mam_75fos[,3:2639]) == 0)) #Se pierden 597, 2637-597 = 2040 spp

Mam_75fos_df <- Mam_75fos[,c("X", "Y", "Richness3")]
Mam_75fos_rast <- rasterize(x = as.matrix(Mam_75fos_df[,1:2]), y = as(x, "SpatRaster"), 
                            values = Mam_75fos_df[,3], update = T)
crs(Mam_75fos_rast) <- "+proj=longlat +datum=WGS84 +no_defs"

Mam_50bs$WC_zone <- as.numeric(Mam_50bs$WC_zone)
Mam_50bs[is.na(Mam_50bs$Fossils),1184] <- 0
Mam_50bs[Mam_50bs$Fossils >= 1,1184] <- 1
Mam_50fos <- Mam_50bs[Mam_50bs$Fossils == 1,]
length(which(colSums(Mam_50fos[,3:1176]) == 0)) #Se pierden 171, 1174-171 = 1003 spp

Mam_50fos_df <- Mam_50fos[,c("X", "Y", "Richness3")]
Mam_50fos_rast <- rasterize(x = as.matrix(Mam_50fos_df[,1:2]), y = as(x, "SpatRaster"), 
                            values = Mam_50fos_df[,3], update = T)
crs(Mam_50fos_rast) <- "+proj=longlat +datum=WGS84 +no_defs"

Mam_25bs$WC_zone <- as.numeric(Mam_25bs$WC_zone)
Mam_25bs[is.na(Mam_25bs$Fossils),306] <- 0
Mam_25bs[Mam_25bs$Fossils >= 1,306] <- 1
Mam_25fos <- Mam_25bs[Mam_25bs$Fossils == 1,]
length(which(colSums(Mam_25fos[,3:298]) == 0)) #Se pierden 22, 296-22 = 274 spp

Mam_25fos_df <- Mam_25fos[,c("X", "Y", "Richness3")]
Mam_25fos_rast <- rasterize(x = as.matrix(Mam_25fos_df[,1:2]), y = as(x, "SpatRaster"), 
                            values = Mam_25fos_df[,3], update = T)
crs(Mam_25fos_rast) <- "+proj=longlat +datum=WGS84 +no_defs"

### Riqueza antes/despues
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
for (i in 3:(ncol(Mam_25fos_f5)-11)) {
  if(any(unique(Mam_25fos_f5[,i]) == 1, na.rm = T)) {
    spp_number <- spp_number + 1
    print(spp_number)
  } 
}


######################### TAXONOMIC FILTER #####################################
Mam_pbdb[Mam_pbdb$genus == "Ptilocercus", "family"] <- "Ptilocercidae"
######################### 75%
Mam_fos_names <- names(which(colSums(Mam_75fos[,3:2639]) != 0))
Combi_6 <- Combi_75bs[Combi_75bs$iucn2020_binomial %in% Mam_fos_names,] #2040 spp
spp_Combir <- sort(setdiff(unique(Combi_6$family), unique(Mam_pbdb$family)))
#Ahora hay 2 familias en Combi_75bs que no estan en Mam_pbdb (Nandiniidae Prionodontidae)
sum(sort(table(Combi_6[Combi_6$family %in% spp_Combir, "family"]))) #2 familias, 3spp
sort(table(Combi_6[Combi_6$family %in% spp_Combir, "family"])) #En total, 3 spp

spp_Mampbdb <- sort(setdiff(unique(Mam_pbdb$family), unique(Combi_6$family)))
#Hay 488 familias en Mam_pbdb que no estan en Combi_r

Combi_6 <- Combi_6[!Combi_6$family %in% spp_Combir,] #2037 obs, 119 familias, 3 spp menos
Mam_pbdb3 <- Mam_pbdb[!Mam_pbdb$family %in% spp_Mampbdb,] #66614 obs, 488 familias menos = 119 familias
Mam_pbdb3 <- as.data.frame(table(Mam_pbdb3$family))
colnames(Mam_pbdb3) <- c("family", "n_obs")
Mam_pbdb3$family <- as.character(Mam_pbdb3$family)

Combi_6$Fam_obs <- NA
family_names <- unique(Mam_pbdb3$family)
for (i in 1:length(family_names)) {
  Combi_6[Combi_6$family == family_names[[i]], "Fam_obs"] <- Mam_pbdb3[Mam_pbdb3$family == family_names[[i]], "n_obs"]
}


#Submuestreo "weighted"
set.seed(1496)
tax75 <- sample(seq_len(nrow(Combi_6)), round(0.75*nrow(Combi_6)), prob = Combi_6$Fam_obs)
Combi_75tax <- Combi_6[tax75,] #1528
Mam_75tax <- Mam_75fos[,c("X", "Y",Combi_75tax$iucn2020_binomial, "Land_mass", "Clim_zone",
                          "GUM_USC", "WC_zone", "Richness", "Land", "Area", "Fossils",
                          "Richness2","Richness3")] #1540
identical(sort(colnames(Mam_75tax[,3:1530])), sort(Combi_75tax$iucn2020_binomial)) #T

richness <- apply(Mam_75tax[,3:1530], 1, function(x) length(x[x>0]))
Mam_75tax[,1541] <- richness
colnames(Mam_75tax)[colnames(Mam_75tax) == "V1541"] <- "Richness4" #New column
Mam_75tax[Mam_75tax$Richness4 == 0, 1541] <- NA

Mam_75tax_fr <- Mam_75tax[,c("X", "Y", "Richness4")]
Mam_75tax_raster <- terra::rasterize(x = as.matrix(Mam_75tax_fr[,1:2]), y = as(x, "SpatRaster"), 
                                     values = Mam_75tax_fr[,3], update = T)
crs(Mam_75tax_raster) <- "+proj=longlat +datum=WGS84 +no_defs"

#Quedandonos con el 75% de las familias con mas observaciones
round(0.75*nrow(Mam_pbdb3)) #89 familias
sort(Mam_pbdb3$n_obs, decreasing = T)[[89]]
my_families <- Mam_pbdb3[Mam_pbdb3$n_obs >=22,1]
Combi_7 <- Combi_6[Combi_6$family %in% my_families,]

Mam_75tax2 <- Mam_75fos[,c("X", "Y",Combi_7$iucn2020_binomial, "Land_mass", "Clim_zone",
                           "GUM_USC", "WC_zone", "Richness", "Land", "Area", "Fossils", 
                           "Richness2","Richness3")] 
identical(sort(colnames(Mam_75tax2[,3:1929])), sort(Combi_7$iucn2020_binomial)) #T

richness <- apply(Mam_75tax2[,3:1929], 1, function(x) length(x[x>0]))
Mam_75tax2[,1940] <- richness
colnames(Mam_75tax2)[colnames(Mam_75tax2) == "V1940"] <- "Richness4" #New column
Mam_75tax2[Mam_75tax2$Richness4 == 0, 1940] <- NA

Mam_75tax2_fr <- Mam_75tax2[,c("X", "Y", "Richness4")]
Mam_75tax2_raster <- terra::rasterize(x = as.matrix(Mam_75tax2_fr[,1:2]), y = as(x, "SpatRaster"), 
                                      values = Mam_75tax2_fr[,3], update = T)
crs(Mam_75tax2_raster) <- "+proj=longlat +datum=WGS84 +no_defs"

####################### 50%
Mam_50fos_names <- names(which(colSums(Mam_50fos[,3:1176]) != 0))
Combi_50fos <- Combi_50bs[Combi_50bs$iucn2020_binomial %in% Mam_50fos_names,] #1003 spp
spp_Combi50fos <- sort(setdiff(unique(Combi_50fos$family), unique(Mam_pbdb$family)))
#Hay 2 familias en Combi_50fos que no estan en Mam_pbdb
sort(table(Combi_50fos[Combi_50fos$family %in% spp_Combi50fos, "family"])) #En total, 3 spp
spp_Mampbdb <- sort(setdiff(unique(Mam_pbdb$family), unique(Combi_50fos$family)))
#Hay 503 familias en Mam_pbdb que no estan en Combi_r

Combi_50fos <- Combi_50fos[!Combi_50fos$family %in% spp_Combi50fos,] #1000 obs, 104 familias, 2 spp menos
Mam_pbdb50fos <- Mam_pbdb[!Mam_pbdb$family %in% spp_Mampbdb,] #66089 obs, 501 familias menos = 104 familias

Mam_pbdb50fos <- as.data.frame(table(Mam_pbdb50fos$family))
colnames(Mam_pbdb50fos) <- c("family", "n_obs")
Mam_pbdb50fos$family <- as.character(Mam_pbdb50fos$family)

Combi_50fos$Fam_obs <- NA
family_names <- unique(Mam_pbdb50fos$family)
for (i in 1:length(family_names)) {
  Combi_50fos[Combi_50fos$family == family_names[[i]], "Fam_obs"] <- Mam_pbdb50fos[Mam_pbdb50fos$family == family_names[[i]], "n_obs"]
}

#Submuestreo "weighted"
set.seed(1496)
tax50 <- sample(seq_len(nrow(Combi_50fos)), round(0.50*nrow(Combi_50fos)), prob = Combi_50fos$Fam_obs)
Combi_50tax <- Combi_50fos[tax50,] #500
Mam_50tax <- Mam_50fos[,c("X", "Y",Combi_50tax$iucn2020_binomial, "Land_mass", "Clim_zone",
                          "GUM_USC", "WC_zone", "Richness", "Land", "Area", "Fossils", 
                          "Richness2","Richness3")] 
identical(sort(colnames(Mam_50tax[,3:502])), sort(Combi_50tax$iucn2020_binomial)) #T

richness <- apply(Mam_50tax[,3:502], 1, function(x) length(x[x>0]))
Mam_50tax[,513] <- richness
colnames(Mam_50tax)[colnames(Mam_50tax) == "V513"] <- "Richness4" #New column
Mam_50tax[Mam_50tax$Richness4 == 0, 513] <- NA

Mam_50tax_fr <- Mam_50tax[,c("X", "Y", "Richness4")]
Mam_50tax_raster <- terra::rasterize(x = as.matrix(Mam_50tax_fr[,1:2]), y = as(x, "SpatRaster"), 
                                     values = Mam_50tax_fr[,3], update = T)
crs(Mam_50tax_raster) <- "+proj=longlat +datum=WGS84 +no_defs"

#Quedandonos con el 50% de las familias con mas observaciones
round(0.5*nrow(Mam_pbdb50fos)) #52 familias
sort(Mam_pbdb50fos$n_obs, decreasing = T)[[52]]
my_families <- Mam_pbdb50fos[Mam_pbdb50fos$n_obs >=138,1]
Combi_50fos2 <- Combi_50fos[Combi_50fos$family %in% my_families,]

Mam_50tax2 <- Mam_50fos[,c("X", "Y",Combi_50fos2$iucn2020_binomial, "Land_mass", "Clim_zone",
                           "GUM_USC", "WC_zone", "Richness", "Land", "Area", "Fossils", "Richness2","Richness3")] 
identical(sort(colnames(Mam_50tax2[,3:797])), sort(Combi_50fos2$iucn2020_binomial)) #T

richness <- apply(Mam_50tax2[,3:797], 1, function(x) length(x[x>0]))
Mam_50tax2[,808] <- richness
colnames(Mam_50tax2)[colnames(Mam_50tax2) == "V808"] <- "Richness4" #New column
Mam_50tax2[Mam_50tax2$Richness4 == 0, 808] <- NA

Mam_50tax2_fr <- Mam_50tax2[,c("X", "Y", "Richness4")]
Mam_50tax2_raster <- terra::rasterize(x = as.matrix(Mam_50tax2_fr[,1:2]), y = as(x, "SpatRaster"), 
                                      values = Mam_50tax2_fr[,3], update = T)
crs(Mam_50tax2_raster) <- "+proj=longlat +datum=WGS84 +no_defs"

####################### 25%
Mam_25fos_names <- names(which(colSums(Mam_25fos[,3:298]) != 0))
Combi_25fos <- Combi_25bs[Combi_25bs$iucn2020_binomial %in% Mam_25fos_names,] #274 spp
spp_Combi25fos <- sort(setdiff(unique(Combi_25fos$family), unique(Mam_pbdb$family)))
#Hay 1 familia en Combi_25fos que no estan en Mam_pbdb
sort(table(Combi_25fos[Combi_25fos$family %in% spp_Combi25fos, "family"])) #En total, 2 spp

spp_Mampbdb <- sort(setdiff(unique(Mam_pbdb$family), unique(Combi_25fos$family)))
#Hay 546 familias en Mam_pbdb que no estan en Combi_r

Combi_25fos <- Combi_25fos[!Combi_25fos$family %in% spp_Combi25fos,] #272 obs, 61 familias, 1 spp menos
Mam_pbdb25fos <- Mam_pbdb[!Mam_pbdb$family %in% spp_Mampbdb,] #55949 obs, 546 familias menos = 61 familias

Mam_pbdb25fos <- as.data.frame(table(Mam_pbdb25fos$family))
colnames(Mam_pbdb25fos) <- c("family", "n_obs")
Mam_pbdb25fos$family <- as.character(Mam_pbdb25fos$family)

Combi_25fos$Fam_obs <- NA
family_names <- unique(Mam_pbdb25fos$family)
for (i in 1:length(family_names)) {
  Combi_25fos[Combi_25fos$family == family_names[[i]], "Fam_obs"] <- Mam_pbdb25fos[Mam_pbdb25fos$family == family_names[[i]], "n_obs"]
}

#Submuestreo "weighted"
set.seed(1496)
tax25 <- sample(seq_len(nrow(Combi_25fos)), round(0.25*nrow(Combi_25fos)), prob = Combi_25fos$Fam_obs)
Combi_25tax <- Combi_25fos[tax25,] #62
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

#Quedandonos con el 25% de las familias con mas observaciones
round(0.25*nrow(Mam_pbdb25fos)) #15 familias
sort(Mam_pbdb25fos$n_obs, decreasing = T)[[15]]
my_families <- Mam_pbdb25fos[Mam_pbdb25fos$n_obs >=1076,1]
Combi_25fos2 <- Combi_25fos[Combi_25fos$family %in% my_families,]

Mam_25tax2 <- Mam_25fos[,c("X", "Y",Combi_25fos2$iucn2020_binomial, "Land_mass", "Clim_zone",
                           "GUM_USC", "WC_zone", "Richness", "Land", "Area", "Fossils", "Richness2","Richness3")] 
identical(sort(colnames(Mam_25tax2[,3:156])), sort(Combi_25fos2$iucn2020_binomial)) #T

richness <- apply(Mam_25tax2[,3:156], 1, function(x) length(x[x>0]))
Mam_25tax2[,167] <- richness
colnames(Mam_25tax2)[colnames(Mam_25tax2) == "V167"] <- "Richness4" #New column
Mam_25tax2[Mam_25tax2$Richness4 == 0, 167] <- NA


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
for (i in 3:(ncol(Mam_25tax_f5)-12)) {
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
for (i in 3:(ncol(Mam_25tax2_f5)-12)) {
  if(any(unique(Mam_25tax2_f5[,i]) == 1, na.rm = T)) {
    spp_number <- spp_number + 1
    print(spp_number)
  } 
}



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
mam_f_moll_df <- terra::as.data.frame(mam_f_moll, xy = T, na.rm = F)
Mam75rs_moll_df <- terra::as.data.frame(Mam_75rs_moll, xy = T, na.rm = F)
Mam50rs_moll_df <- terra::as.data.frame(Mam_50rs_moll, xy = T, na.rm = F)
Mam25rs_moll_df <- terra::as.data.frame(Mam_25rs_moll, xy = T, na.rm = F)
Mam75bs_moll_df <- terra::as.data.frame(Mam_75bs_moll, xy = T, na.rm = F)
Mam50bs_moll_df <- terra::as.data.frame(Mam_50bs_moll, xy = T, na.rm = F)
Mam25bs_moll_df <- terra::as.data.frame(Mam_25bs_moll, xy = T, na.rm = F)
Mam75bsfos_moll_df <- terra::as.data.frame(Mam_75fos_moll, xy = T, na.rm = F)
Mam50bsfos_moll_df <- terra::as.data.frame(Mam_50fos_moll, xy = T, na.rm = F)
Mam25bsfos_moll_df <- terra::as.data.frame(Mam_25fos_moll, xy = T, na.rm = F)
Mam_75tax_moll_df <- terra::as.data.frame(Mam_75tax_moll, xy = T, na.rm = F)
Mam_50tax_moll_df <- terra::as.data.frame(Mam_50tax_moll, xy = T, na.rm = F)
Mam_25tax_moll_df <- terra::as.data.frame(Mam_25tax_moll, xy = T, na.rm = F)
Mam_75tax2_moll_df <- terra::as.data.frame(Mam_75tax2_moll, xy = T, na.rm = F)
Mam_50tax2_moll_df <- terra::as.data.frame(Mam_50tax2_moll, xy = T, na.rm = F)
Mam_25tax2_moll_df <- terra::as.data.frame(Mam_25tax2_moll, xy = T, na.rm = F)


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

identical(mam_f_longlat_sf$last, mam_f_moll_sf$last) #TRUE
#mam_r_sf_juntas <- cbind(mam_r_moll_sf, mam_r_longlat_sf$geometry)
#mam_r_dfforplot <- as.data.frame(mam_r_sf_juntas$last)
#for (i in 1:259200) {
#  print(i)
#  mam_r_dfforplot[i,2] <- mam_r_sf_juntas$geometry[[i]][2]
#  mam_r_dfforplot[i,3] <- mam_r_sf_juntas$geometry.1[[i]][2]
#}
#colnames(mam_r_dfforplot) <- c("Richness", "Moll", "Longlat")
#write.table(mam_r_dfforplot, "sf_coords.csv", row.names = F)
mam_r_dfforplot <- read.table("sf_coords2.csv")

identical(mam_r_moll_sf$geometry, Mam_25tax_moll_sf$geometry) #TRUE

db_forslopes <- data.frame(mam_r_dfforplot$Longlat, mam_r_dfforplot$Moll,
                           mam_r_dfforplot$Richness, mam_f_moll_sf$last, Mam75rs_moll_sf$last,
                           Mam50rs_moll_sf$last, Mam25rs_moll_sf$last, Mam75bs_moll_sf$last, Mam50bs_moll_sf$last,
                           Mam25bs_moll_sf$last, Mam75bsfos_moll_sf$last, Mam50bsfos_moll_sf$last,
                           Mam25bsfos_moll_sf$last, Mam_75tax_moll_sf$last, Mam_50tax_moll_sf$last,
                           Mam_25tax_moll_sf$last, Mam_75tax2_moll_sf$last, Mam_50tax2_moll_sf$last,
                           Mam_25tax2_moll_sf$last)
colnames(db_forslopes) <- c("Longlat", "Moll", "Original", "UcS", "rs75", "rs50",
                            "rs25", "bs75", "bs50", "bs25", "foss75", "foss50",
                            "foss25", "tax75", "tax50", "tax25", "tax2_75", 
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
