##################################################################################:
##################################################################################:
## Part 2 Stacked Species Data and Modelling
##################################################################################:
##################################################################################:
## requires Part 1 results: covariates and shapeAOI_UTM
##################################################################################:
## author Daniel Gruschwitz
##################################################################################:


#--------------------------------------------------------------------------------:)
### Environment setup #########################################################
#--------------------------------------------------------------------------------:)

# required packages  
reqp <- c("caret",  "devtools", "doParallel","dplyr",    
          "foreach", "fuzzySim","geodata", "raster", 
          "pracma","rgbif", "RStoolbox", "SSDM", "sf",  
          "stringr", 'terra', "tidyverse", "wellknown") 

# check if already installed
ip <- reqp[!(reqp %in% installed.packages()[, "Package"])] 



if(length(ip)) install.packages(ip, dependencies = TRUE) # install packages
lapply(reqp, require, character.only = TRUE) # load packages

rm(reqp, ip) # remove vector


reqp <- c("scrubr", "flexsdm" )
reqp2 <- c("ropensci/scrubr", "sjevelazco/flexsdm")
ip <- reqp[!(reqp %in% installed.packages()[, "Package"])] # subset packages


if(length(ip)) devtools::install_github(reqp2, dependencies = TRUE, upgrade = "ask") # install packages


lapply(reqp, require, character.only = TRUE) # load packages

rm(reqp, ip, reqp2) # remove vector

#--------------------------------------------------------------------------------:)
### define directory path #######################################################
#same directory as script 1
directory <- "your_directory"
setwd(directory)
#--------------------------------------------------------------------------------:)
### load data of script 1 #######################################################
covariates <- rast("./covariates/covariates_new.grd")

plot(covariates)
names(covariates)

shapeAOI <- st_read("Species_modelling_env.gpkg", "shapeAOI")
shapeAOI_UTM <- st_read("Species_modelling_env.gpkg", "shapeAOI_UTM")
crs_UTM <- terra::crs(shapeAOI_UTM, describe= F)


#--------------------------------------------------------------------------------:)
### define species ##############################################################
#--------------------------------------------------------------------------------:)
  
### latin scientific species name

#tree species of Austria e.g. based on Baumartenatlas https://bfw.ac.at/700/2092_1.html

species_list <- c("Pinus nigra", 'Pinus cembra L.', 'Pinus sylvestris', 'Picea abies', 'Abies alba', 'Larix decidua', 'Fagus sylvatica', 
                  'Pyrus pyraster', 'Quercus robur', 'Quercus petraea', 'Quercus cerris', 'Quercus pubescens', 'Carpinus betulus', 'Fraxinus excelsior', 
                  'Acer pseudoplatanus', 'Acer campestre', 'Acer platanoides','Ulmus glabra', 'Ulmus minor', 'Castanea sativa', 'Robinia pseudoacacia', 
                  'Prunus avium', 'Sorbus aucuparia', 'Sorbus torminalis', 'Sorbus aria (L.) Crantz', 'Malus sylvestris', 'Betula pendula', 'Betula pubescens', 
                  'Alnus glutinosa (L.) Gaertn.', 'Alnus incana', 'Tilia cordata','Tilia platyphyllos', 'Populus tremula', 'Populus alba', 'Populus nigra',
                  'Salix caprea', 'Salix alba L.', 'Taxus baccata', 'Pseudotsuga menziesii', 'Juglans regia') 


#  ---------------------------------------------------------------------------------.
#---------------------------------------------------------------------------------:)
### GBIF species data ###########################################################
#---------------------------------------------------------------------------------:)
#  ---------------------------------------------------------------------------------.
# https://www.r-bloggers.com/2021/03/downloading-and-cleaning-gbif-data-with-r/
# presence data preparation 

### function for downloading and data cleaning of Gbif data
### if error occurs likely that there are not enough data entries at one point

species_data_preparation <- function (species_name) {
  
  ### download GBIF data (species presence data with lat/lon)
  ### geometry in wkt format, limit set to maximum
  ### filter for year and coordinate uncerntainty
  ### speeds up data fetching 
  gbif_data <- occ_data(scientificName = species_name, hasCoordinate = TRUE, geometry = sf_convert(shapeAOI), 
                        geom_big = "bbox", coordinateUncertaintyInMeters = '0,250', year = '1995,2022', limit = 100000)
  
  # exclude data with no entries from final dataset
  if (is.null(dim(gbif_data$data))) {next
    print(paste0("no entries ", species_list[i]))
    } 
   
    #keep only relevant columns
    species_coords <- gbif_data$data[ , c("decimalLongitude", "decimalLatitude", "occurrenceStatus", "coordinateUncertaintyInMeters", "institutionCode",  
                                          "year", "month", "day")] %>% date_create(year, month, day) 
    species_coords <- species_coords[,!names(species_coords) %in% c("year","month", "day")]  #year/month/day united in one column
    print(paste0('all gbif data entries: ', nrow(species_coords))) 
    if (nrow(species_coords) < 20) {next
      print(paste0("not enough occurences for ", species_list[i])) 
    } 
      
      ### cleaning data
      
      
      # all absence data will be removed (often just stratified sampled data)
      species_coords <- species_coords[!(species_coords$occurrenceStatus 
                                         %in% c("absent", "Absent", "ABSENT", "ausente", "Ausente", "AUSENTE")),]
      # all records prior to 1995 will be removed
      
      species_coords$date <- as.Date(species_coords$date)
      species_coords <- species_coords[(!is.na(species_coords$date)),]
      species_coords <- species_coords[species_coords$date > '1995-01-01',]
      
      vague_rows <- which(is.na(species_coords$coordinateUncertaintyInMeters) 
                          | species_coords$coordinateUncertaintyInMeters >= 250)
      if (length(vague_rows) > 0){   # if statement is necessary because null vector would eliminate all entries 'loiseleuria procumbens'
      species_coords <- species_coords[-vague_rows,]}   
      print(paste0('gbif data entries with location uncertainty < 1km: ', nrow(species_coords)))
      
      #   removal of points with erroneous location
      species_coords <- coord_incomplete(coord_imprecise(coord_impossible(coord_unlikely(species_coords))))
      print(paste0('gbif data entries without erroneous location: ', nrow(species_coords)))    
      
      if (nrow(species_coords) < 10){next
        print(paste0("not enough occurences for ", species_list[i])) 
      } 
        
          # create column with species name
          species <- rep(species_name, nrow(species_coords))
          species_coords <- cbind(species, species_coords)
          
          # creation of point data
          species_points <- st_as_sf(species_coords, coords = c("decimalLongitude", "decimalLatitude"), crs=crs_WGS84) %>% st_transform(crs_UTM)
          
          
          ### sample bias correction: only 1 sample per raster cell (250x250m)
          species_points_sampled <- gridRecords(rst = covariates,
                                                pres.coords = data.frame(st_coordinates(species_points)), 
                                                absences = F,
                                                na.rm = T)
          
          print(paste0('gbif data entries after sample bias correction: ', nrow(species_points_sampled))) 
          
          species_points_sampled <- na.omit(species_points_sampled)
              
          print(paste0('gbif data entries after NA removal: ', nrow(species_points_sampled)))
          
          # cannot be joined back via coordinates because centroid coordinates of cells are kept
          # create again spatial points data
          species_points_sampled <-   st_as_sf(species_points_sampled, coords = c('x', 'y'), crs=crs_UTM, remove= F)
          
          
          # join back gbif data via nearest point
          # potentially false joins
          
          species_data_cleaned <- data.frame()
          for (j in 0:nrow(species_points_sampled)) {
            nearest <- st_nearest_feature(species_points_sampled[j,], species_points)
            species_points_sampled2 <- cbind(species_coords[nearest,], species_points_sampled[j,])
            if (j > 0) {
              species_data_cleaned <- rbind(species_data_cleaned, species_points_sampled2)
            }
          } 
          
          return(species_data_cleaned)
  }

### setup for parallel processing

# Use the detectCores() function to find the number of cores in system
no_cores <- detectCores()-1

# This function has to be used if something goes wrong in the first try
 unregister_dopar <- function() {
   env <- foreach:::.foreachGlobals
   rm(list=ls(name=env), pos=env)
 }
 unregister_dopar()
 
# Setup cluster
clust <- makeCluster(no_cores) 



registerDoParallel(clust)

# export the loaded data to each worker
# maybe this step is redundant because it is automatically done by foreach function
clusterExport(clust, list("species_data_preparation", "species_list", "covariates", "shapeAOI"))

# loop through different species resulting in one dataset for all species

tic('Presence Points preparation: ')
species_data <- foreach(i = 1:length(species_list), 
                           .combine=rbind,
                           .packages = c("scrubr", "fuzzySim", "sf", "rgbif", "wellknown", "tcltk"),
                           .inorder = F,
                           .verbose= T) %dopar% {
          temporary <- species_data_preparation(species_name = species_list[i])
          if (class(temporary) == "data.frame"){       # this is to check if something was downloaded otherwise merge would fail
            result <- data.frame(temporary)
            #rm(temporary)
            }
             
                          }                
toc()                             
stopCluster(clust)                             

species_data <- data.frame()
tic('Presence Points preparation: ')                           
for(i in 1:length(species_list)){
  print(species_list[i])
  tic(paste("time for: ", species_list[i]))
      temporary <- species_data_preparation(species_name = species_list[i])
    # make sure no empty data is joined which would cause error
    #if (class(temporary) == "data.frame"){
      species_data <- rbind(species_data, temporary)
      #}
      rm(temporary)
  toc()
}
toc()

# remove columns that are no longer needed
species_data <- subset(species_data, select = -c(cells, decimalLongitude, decimalLatitude, occurrenceStatus, geometry))

# reorder alphabetically
species_data <- species_data[order(species_data$species),]
rownames(species_data) <- seq(1:nrow(species_data))

# save dataset
write.csv2(species_data, file = "tree_species.csv", row.names=F)

#load dataset
#species_data <- read.csv2("tree_species.csv", header= T)

#---------------------------------------------------------------------------------:)
###  occurence count per species ################################################
#---------------------------------------------------------------------------------:)
  
occurence_number <- count(species_data, species)

# which species have 0 valid presence data?
no_occurences <- species_list[-which(species_list %in% occurence_number[,1])]
no_occurences <- data.frame(cbind(no_occurences, rep(0, length(no_occurences))))
names(no_occurences) <- c("species", "n")

# merge into one dataset and save
occurence_number <- rbind(occurence_number, no_occurences)
occurence_number <- occurence_number[order(occurence_number$species),]
rm(no_occurences)
write.csv2(occurence_number, file = "tree_occurence_numbers.csv")

# species with too few entries could be excluded
# remove sorbus e.g. torminalis because of too few occurence entries
# species_data <- species_data[species_data$species != "Sorbus torminalis",]

#---------------------------------------------------------------------------------:)
### show presence data distribution #############################################
#---------------------------------------------------------------------------------:)
  
### create an empty raster
raster_blank <- matrix(NA, nrow = nrow(covariates), ncol = ncol(covariates)) %>% 
  rast(crs= crs(covariates),
       extent= ext(covariates)) %>% 
  mask(shapeAOI_UTM)

### aggregate to 5km raster cells
raster_blank <- terra::aggregate(raster_blank, fact=20)

point_count <- terra::rasterize(st_as_sf(species_data[,c("x", "y")], coords = c('x','y'), crs=crs_UTM), 
                                y= raster_blank,
                                fun='length')

### save plot as png image
png("species occurrences.png", width= 300, height= 150, unit="mm", res=200)

plot(point_count, col=heat.colors(50, rev=T), colNA= c("lightgrey"),
     #main= "Number of species occurences per 5km raster cell"
) 
plot(vect(shapeAOI_states), lwd=2, add=T)

dev.off()

### calculate point density per federal state
# federal states of nation used for visualization and statistics
shapeAOI_states <- geodata::gadm(country= country, level = 1, path=getwd())   # download shapes from GADM
shapeAOI_states <-    shapeAOI_states[, "NAME_1"] %>% 
  st_as_sf() %>% 
  st_transform(crs =crs_UTM)  

# state area calculation
shapeAOI_states$area <- st_area(shapeAOI_states)/1000000

# number of presence points within one state
points_contained <- st_contains(shapeAOI_states, 
                                st_as_sf(species_data[,c("x", "y")], coords = c('x','y'), crs=crs_UTM))
names(points_contained) <- shapeAOI_states$NAME_1
points_contained <- sapply(points_contained, function (state){length(state)} )

# calculation of point density and export of dataframe
point_density <- data.frame(round(rbind(shapeAOI_states$area, points_contained,  points_contained/shapeAOI_states$area), 2),
                            row.names = c("Area_in_sqkm", "Point number", "Point density"))
write.csv2(point_density, "density.csv")



  #---------------------------------------------------------------------------------.
#---------------------------------------------------------------------------------:)
### create pseudo-absence data ###################################################
#---------------------------------------------------------------------------------:)
  #---------------------------------------------------------------------------------.
  
  
# non NA values raster mask to create points only in non NA area

covariates_mask <- !is.na(sum(mask(covariates, shapeAOI_UTM)))
covariates_mask <- terra::classify(covariates_mask, matrix(c(-Inf, 0, NA), ncol = 3, byrow=T))

terra::writeRaster(covariates_mask, "./Cov_Mask.tif", gdal="GTiff", overwrite=T)
plot(covariates_mask)

# Generate background data using flexSDM package
# applying a 2km buffer around presence points
# placing absence data only in low fit areas of BIOCLIM model 
# create as many absence points as species points

registerDoParallel(clust)

# potentially trouble with foreach and terra package (SpatRaster class)
set.seed(10)
background_data <- foreach(i = 1:length(unique(species_data$species)), 
        .combine=rbind,
        .packages = c("terra", "flexsdm", "sf", "magrittr"),
        #.errorhandling = "pass",
        .noexport= c("background_points","background_points_spatial"),
        .inorder = F,
        .verbose= T) %dopar% { 
  background_points <- data.frame(sample_pseudoabs(species_data[species_data$species == unique(species_data$species)[i], c("species","x","y", "presence") ],
                                                   x= "x",
                                                   y= "y",
                                                   n= nrow(species_data[species_data$species == unique(species_data$species)[i],]),
                                                   method = c("geo_env_const", width= "2000", env=rast("./covariates/predictors.grd")),
                                                   rlayer = rast("./Cov_Mask.tif"),
                                                   maskval = 1,
                                                   sp_name = unique(species_data$species)[i])) 
  names(background_points) <- c('species', 'x', 'y', 'presence')
  background_points_spatial <- st_as_sf(background_points, coords = c('x','y'), crs=crs_UTM)
  background2 <- data.frame(background_points, terra::extract(rast("./covariates/predictors.grd"), background_points_spatial, method = "bilinear"))
    }

stopCluster(clust)

## create further columns as in presence data and save
background <- data.frame(background_data[,1:4], date = Sys.Date(), coordinateUncertaintyInMeters = NA, institutionCode = NA, background_data[,6:ncol(background_data)] )
write.csv2(background, "background_env_const.csv")

## combine presence and absence data
species_data2 <- rbind(background, species_data) 
species_data2 <- species_data2[order(species_data2$species),]

# save points data in geopackage
for (i in 1:length(unique(species_data2$species))){
  species_points_spatial <- st_as_sf(species_data2[species_data2$species == unique(species_data2$species)[i],],
                                     coords = c('x','y'), remove=F, crs=crs_UTM)
  st_write(species_points_spatial, "Species_modelling_env.gpkg", unique(species_data2$species)[i], append=T)
}
rm(species_points_spatial)


#---------------------------------------------------------------------------------:)
### partition in training and testing data #####################################
#---------------------------------------------------------------------------------:)
  
set.seed(0)
split <- createDataPartition(species_data$species, p=0.75, list=F)
pres_train <- species_data[split,]
pres_test <- species_data[!(rownames(species_data) %in% split) ,]

set.seed(0)
split <- createDataPartition(background$species, p=0.75, list=F)
background_train <- background[split,]
background_test <- background[!(rownames(species_data) %in% split) ,]

traindata <- data.frame(rbind(pres_train, background_train))
traindata <- traindata[order(traindata$species),]
rownames(traindata) <- seq(1:nrow(traindata))
testdata <- data.frame(rbind(pres_test, background_test))
testdata <- testdata[order(testdata$species),]
rownames(testdata) <- seq(1:nrow(testdata))


  #---------------------------------------------------------------------------------:)
#---------------------------------------------------------------------------------.
### Ensemble Modelling #########################################################
#---------------------------------------------------------------------------------.
  #---------------------------------------------------------------------------------:)
  
  
registerDoParallel(clust)

covariates2 <- stack("./covariates/predictors.grd") # SSDM does only accept raster package stacks

dir.create("Modells")
# create for each species an ensemble model 
# each ensemble model is based on an Random Forest Model (RF), Maxent and Support Vector Machine (SVM)
# those 3 models show usually best performance
# cross validation is here based on k-fold subsampling
# model accuracy for ensembling is 0.5

tic('ensemble time complete:')
foreach(i = 1:length(unique(traindata$species)), 
              .packages = c("SSDM", "terra"),
              .inorder = F,
              .verbose= T) %dopar% {
  
    pracma::tic(paste0('ensemble time: ', unique(traindata$species)[i]))
        ensemble_modelling(algorithm = c('RF', 'MAXENT', 'SVM'), 
              Occurrences= traindata[traindata$species == unique(traindata$species)[i],] ,
              Env= covariates2, # does only accept raster package stacks
              Xcol = "x",
              Ycol = "y",
              cv = "k-fold",
              cv.param = c(5,10),
              Pcol = "presence",
              rep = 1,
              name= unique(traindata$species)[i],
              path = paste0(getwd(), "/Modells"),
              save=T,
              #tmp = T,
              ensemble.metric = 'Kappa',
              ensemble.thresh = 0.5,
              weight = T,
              uncertainty = T,
              ntree = 3000,
              nodesize= 0.8,
              verbose = T,
              cores = 7)
    pracma::toc()}
toc()

stopCluster(clust)


### load all the models as single variables (required by stacking operator)
### check if models exists because some might be dropped because of too low kappa

for (i in 1:length(unique(species_data$species))) {
  if (file.exists(paste0(getwd(), "/Modells/", unique(species_data$species)[i]))) {
  assign(x = gsub(pattern= " ", x= unique(species_data$species)[i], fixed = T, replacement = "_"), 
         value = load_esdm(unique(species_data$species)[i], paste0(getwd(), "/Modells")))
  }
}

plot(Acer_campestre)




#---------------------------------------------------------------------------------:)
# Ensemble validation ###########################################################
#---------------------------------------------------------------------------------:)
  
# load all models in a list
load_esdm_func <- function(tree) {if (file.exists(paste0(getwd(), "/Modells/", unique(species_data$species)[i]))) {
  ESDM <- load_esdm(tree, paste0(getwd(), "/Modells"))
  return(ESDM)}
}
Ensemple_list <- lapply(unique(species_data$species), FUN = load_esdm_func)

# validate binary output map based on separated testdata
validate_ESDM <- function(modell) {validateMap(modell@binary, 
                                         st_as_sf(testdata[unique(testdata$species == str_sub(modell@name, end= -14) ),],
                                                  coords = c('x','y'), crs=crs_UTM), 
                                         responseCol = "presence")}


registerDoParallel(clust)

clusterExport(clust, c("testdata", "crs_UTM"))
clusterEvalQ(clust, library(RStoolbox))
clusterEvalQ(clust, library(sf))
clusterEvalQ(clust, library(stringr))

Validation_list <- parLapply(clust, Ensemple_list, validate_ESDM)
stopCluster(clust)

# combine data and export dataframe
Ensemble_validation <- mapply(c, Ensemple_list,Validation_list, SIMPLIFY=F)
names(Ensemble_validation) <- unique(species_data$species)


validation <- sapply(Ensemble_validation, function(vali) {round(vali$performance$overall,4) })
write.csv2(validation, "tree_species_validation.csv")




  #---------------------------------------------------------------------------------:)
#---------------------------------------------------------------------------------.
### Stacking of all the ensemble models #########################################
#---------------------------------------------------------------------------------.
  #---------------------------------------------------------------------------------:)
  
# ### create a character string of all valid model names for stacking input 
# not working yet
# Esdm_list <- character()
# for (i in 1:length(unique(species_data$species))) {
#   if (i == 1) {
#     if (file.exists(paste0(getwd(), "/Modells/", unique(species_data$species)[i]))) {
#       Esdm_list <- gsub(pattern= " ", x= unique(species_data$species)[i], fixed = T, replacement = "_")}
#   } 
#   if(i >1) {
#   if (file.exists(paste0(getwd(), "/Modells/", unique(species_data$species)[i]))) {
#   Esdm_list <- append(x = Esdm_list, values = gsub(pattern= " ", x= unique(species_data$species)[i], fixed = T, replacement = "_"))}}
# }

tic('stacking time complete:')
Tree_SSDM <- SSDM::stacking(Abies_alba, Acer_campestre, Acer_platanoides, 
                            Acer_pseudoplatanus, `Alnus_glutinosa_(L.)_Gaertn.`, 
                            Alnus_incana, Betula_pendula, Betula_pubescens, 
                            Carpinus_betulus, Castanea_sativa,  
                            Fagus_sylvatica, Fraxinus_excelsior, Juglans_regia, 
                            Larix_decidua, Picea_abies, Pinus_cembra_L., 
                            Pinus_nigra, Pinus_sylvestris, Populus_alba, 
                            Populus_nigra, Populus_tremula, Prunus_avium, Pseudotsuga_menziesii, 
                            Pyrus_pyraster, Quercus_cerris, Quercus_petraea, 
                            Quercus_pubescens, Quercus_robur, Robinia_pseudoacacia, 
                            Salix_alba_L., Salix_caprea, `Sorbus_aria_(L.)_Crantz`, 
                            Sorbus_aucuparia, Taxus_baccata, Tilia_cordata, 
                            Tilia_platyphyllos, Ulmus_glabra, Ulmus_minor,
                      name= 'Tree_SSDM_new',
                      method= 'pSSDM',
                      endemism = c('WEI', 'NbOcc'),
                      verbose= T)
toc()


save.stack(Tree_SSDM, name = 'Tree_SSDM', path = getwd())
plot(Tree_SSDM)
Tree_SSDM@evaluation

png("Tree_richness.png", width= 700, height= 385, unit="mm", res=500)


terra::plot(Tree_SSDM@diversity.map, mar=c(2, 2, 2, 4),
            colNA="lightgrey") 
terra::plot(urban_mask, col="#EC7C66", border="#F79F9D", lwd=0.1, add=T)
plot(water_mask, col= "#9799FD", border="#9799FD", lwd=0.1, add =T) 
plot(glacier_mask, col="#9EE2E2", border="#9EE2E2", lwd=0.1, add=T)
plot(vect(shapeAOI_UTM), lwd= 3, add=T)
legend("bottomright", c("urban", "water", "glacier"), 
       pch = 19, 
       col = c("#EC7C66","#9799FD","#9EE2E2"),
       title= "CORINE\n landcover masks",
       bty = "n",
       xjust = 1,
       cex=2)
  
dev.off()









