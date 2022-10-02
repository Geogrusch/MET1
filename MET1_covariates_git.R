###############################################################################:
###  MT1 Spatial Modeling and Prediction final project Daniel Gruschwitz
###  stacked species modeling
###  exemplified for tree species richness in Austria
###  Part 1 environmental covarariates layer
###############################################################################:
###  platform: x86_64-w64-mingw32
###  RVersion: R version 4.2.1 (2022-06-23 ucrt)
###############################################################################:
### written by Daniel Gruschwitz (EAGLE student)
### date: 01.10.2022
### https://github.com/Geogrusch/MET1
###############################################################################:


#--------------------------------------------------------------------------------:)
### Software and account requirements #########################################
### requires SAGA download for terrain analysis: https://sourceforge.net/projects/saga-gis/
# path to SAGA folder with saga.exe
path_SAGA <- "your_SAGA_path"

### requires GDAL: if not installed then automatic R download
gdal_directory = "your_GDAL_path"

### reguires NASA Earthdata account https://urs.earthdata.nasa.gov/home: 
# NASA Earthdata server user name
NASA_user <- "your_username"

### requires Google earth engine account and python installation for rgee package 
# create your own python GEE environment before named ee
python_env <- "your_python_environment/envs/ee"
GEE_emailaddress <- "your_email"

### non automatic data download (save in working directory): 
### Corine Landcover 100m tif https://land.copernicus.eu/pan-european/corine-land-cover/clc2018?tab=download
### Corine landcover classes description: https://www.eea.europa.eu/data-and-maps/data/corine-land-cover-2/corine-land-cover-classes-and/clc_legend.csv
### Biomass data CEDA: choose tiles according to Lat/Lon https://data.ceda.ac.uk/neodc/esacci/biomass/data/agb/maps/v2.0/geotiff/2018

#--------------------------------------------------------------------------------:)
### Directory setup ############################################################
#--------------------------------------------------------------------------------:)
### define directory path
# better create an empty folder where all the files will be saved
directory <- "your_directory"
setwd(directory)


#--------------------------------------------------------------------------------:)
### define AOI ##################################################################
#--------------------------------------------------------------------------------:)

# processing for entire country takes a while
country <- "Austria"
level <- 0                   # if entire country is wished set level to 0, for states within a country to 1
#  states <- c("Vorarlberg")    # define states

# That's it for the manual setup!

#----------------------------------------------------------------------------------:)
### Environment setup #########################################################
#----------------------------------------------------------------------------------:)

  
rp <- c("caret", "corrplot", "devtools", "dplyr", "fuzzySim","geodata", "getPass",   
        "MODIStsp", "raster", "reticulate", 
        "rfUtilities","rgdal", "Rsagacmd", "RStoolbox", "sf",  
        "stringr", 'terra', "tidyverse") # required packages

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


#----------------------------------------------------------------------------------:)
# GEE setup ######################################################################
# via Google Earth Engine Landsat data will be processed and downloaded
# alternatively processing can also be done via GEE directly, exported and imported to R later on
 
# Google Earth Engine Connection setup via rgee and reticulate
# rgee is based on python 


# python environment
conda_list()
reticulate::py_version()
use_condaenv(python_env)

# installation only has to be run once
# requires Python installation, numpy and ee python packages, Earth Engine Python API setup and Google Cloud installation
# gcloud installation: https://cloud.google.com/sdk/docs/install#deb
# credentials have to be saved
# rgee::ee_install_set_pyenv(py_path = python_env,
#                            py_env = "ee",
#                            Renviron = paste0(python_env, "/RGEE) )


ee_check() 
ee_check_credentials()
# if issues appear check for recent changes on https://github.com/r-spatial/rgee
# e.g. concerning gcloud https://github.com/r-spatial/rgee/issues/267


# initialize setup has to be done every session
ee_Initialize(user = GEE_emailaddress, drive=T)

ee_user_info()


#----------------------------------------------------------------------------------:)

# GDAL setup
# optional if more GDAL versions installed
gdal_directory = "C:/OSGeo4W64/bin"
gdal_setInstallation(search_path = gdal_directory)


# SAGA setup (required for terrain analysis)

saga <- saga_gis(saga_bin = paste0(path_SAGA,"/saga_cmd.exe"), raster_format = "GeoTIFF")



#----------------------------------------------------------------------------------:)
### get GADM data ##############################################################
#----------------------------------------------------------------------------------:)


shapeAOI <- geodata::gadm(country= country, level = level, path=getwd())   # download outline shape from GADM

if (level == 0){  # processing for entire country
  shapeAOI <- shapeAOI[,  -(1:ncol(shapeAOI))] %>% 
    st_as_sf()
}
if (level == 1){ # processing for selected states
  shapeAOI@data[["NAME_1"]] # print overview of available states names
  shapeAOI <- shapeAOI[shapeAOI$NAME_1 %in% states,  -(1:ncol(shapeAOI))] %>% 
    st_as_sf() %>% 
    st_union()  %>%   # dissolve states' boundaries
    sf::st_sf()       # convert sfc geometry collection (always the result of st_union) back to single sf object
}

### alternatively load your own shapefile 

#shapeAOI <- st_read(path.shp")


### reproject to UTM and setting up variables

  if (!st_crs(shapeAOI)== 4326){
    shapeAOI <- st_transform(shapeAOI,4326)  # project to a common crs first
  } 

st_write(shapeAOI, "Species_modelling_env.gpkg", "shapeAOI", append=T)

coords <- st_centroid(shapeAOI) %>%          # center coordinates to derive UTM zone
                st_coordinates()
UTM_zone <- floor((coords[1,] + 180) / 6) + 1


crs_UTM  <-  paste0( "+proj=utm +zone=", UTM_zone[1], " +datum=WGS84 +units=m +no_defs")  # different crs strings definition
crs_WGS84 <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
crs_igh <- '+proj=igh +lat_0=0 +lon_0=0 +datum=WGS84 +units=m +no_defs'

shapeAOI_UTM <- st_transform(shapeAOI, crs =crs_UTM)                                      # UTM reprojection
st_write(shapeAOI_UTM, "Species_modelling_env.gpkg", "shapeAOI_UTM", append=T)

# different bbox definitions:   
extent_ll <- ext(shapeAOI)                                                               
extent_AOI_UTM <- ext(st_buffer(shapeAOI_UTM, dist = 5000))
extent_AOI_igh <- ext(st_transform(st_buffer(shapeAOI_UTM, dist = 5000), crs=crs_igh))


################################################################################:)
### create Covariates ###########################################################
################################################################################:)

################################################################################:)
### ISRIC SOIL Data #############################################################
# potentially better download function in geodata package soil_world function

# https://www.isric.org/explore/soilgrids
# https://git.wur.nl/isric/soilgrids/soilgrids.notebooks/-/blob/master/markdown/wcs_from_R.md
# https://git.wur.nl/isric/soilgrids/soilgrids.notebooks/-/blob/master/markdown/webdav_from_R.md
# https://rpubs.com/ials2un/soilgrids_webdav # this is the one implemented here


dir.create("soildata")

### bounding box in interrupted Goode homolosine projection, it is better handled by the webservice

bb=c(round(as.numeric(extent_AOI_igh@xmin), 0), round(as.numeric(extent_AOI_igh@ymax), 0),round(as.numeric(extent_AOI_igh@xmax), 0),round(as.numeric(extent_AOI_igh@ymin), 0)) 

# define soil metrics
soilgriddata <- c("bdod", "cec", "cfvo", "clay", "nitrogen", "phh2o", "sand", "silt", "soc", "ocd")
# define soil depths
soildepth <- c("_0-5","_5-15", "_15-30", "_30-60", "_60-100", "_100-200")

### download different soil layers in 250m res for six different depths
sg_url <- "/vsicurl/https://files.isric.org/soilgrids/latest/data/"

for (j in 1:length(soildepth)) {
  for (i in 1:length(soilgriddata)){
  
  soildata  <-  paste0(soilgriddata[i], "/", soilgriddata[i],  soildepth[j], "cm_mean.vrt")
  soiloutput <- paste0("./soildata/", soilgriddata[i], soildepth[j], ".tif")
  gdal_translate(paste0(sg_url,soildata), soiloutput ,
                 tr=c(250,250),
                 projwin=bb,
                 projwin_srs =crs_igh,
                 verbose=TRUE)
  }
}  
  
rm(soildata, soiloutput, sg_url, i, bb)

# list all the downloaded files
soil_files <- list.files("./soildata/", full.names = T, pattern = "*.tif")

# reproject and create a raster stack
# the soil raster is the projection reference for the other covariates!!
soil_raster_tmp <- rast()
for (i in 1:length(soil_files)){
      terra::add(soil_raster_tmp) <- rast(soil_files[i]) %>% 
      terra::project(y=crs_UTM, method="bilinear", res=c(250,250))
}

names(soil_raster_tmp)

### soilgrids stores data in integer format therefore conversion in conventional units
# this is skipped because scaling is applied later on 
# for (i in 1:nlyr(soil_raster_tmp)) {
#   if (names(soil_raster_tmp[[i]]) %in% grep("bdod*", names(soil_raster_tmp), value =T) | names(soil_raster_tmp[[i]]) %in% grep("nitrogen*", names(soil_raster_tmp), value =T) ) {
#     soil_raster_tmp[[i]] <- soil_raster_tmp[[i]]/100
#   } else {
#     soil_raster_tmp[[i]] <- soil_raster_tmp[[i]]/10}
# }




# pH to H3O+ concentration (metric scaled instead of ordinal)
# afterwards very small values
 
phh2o <- grep("phh2o*", names(soil_raster_tmp), value =T)
for (i in 1:nlyr(soil_raster_tmp)) {
  if (names(soil_raster_tmp[[i]]) %in% phh2o) {
      new_H3O_layer <- soil_raster_tmp[[i]]/10 
      new_H3O_layer <- '^'(10,(-soil_raster_tmp[[i]]))
      names(new_H3O_layer) <- paste0("H3O_", str_sub(names(soil_raster_tmp[[i]]), start= 7))
      terra::add(soil_raster_tmp) <- new_H3O_layer
      rm(new_H3O_layer)
  }
}


soil_raster_tmp <- terra::subset(soil_raster_tmp, phh2o, negate=T)  
soilgriddata <- c("bdod", "cec", "cfvo", "clay", "nitrogen", "sand", "silt", "soc", "ocd", "H3O")

# compute PCAs for each soil metric over the different soil depths 
# reduces collinearity and absolute number of input numbers

soil_raster <- rast()

for(p in 1:length(soilgriddata)) {
                          soil_subset <- terra::subset(soil_raster_tmp ,grep(soilgriddata[p], names(soil_raster_tmp)))
                          soil_pca <- rasterPCA(soil_subset, nComp=1, spca = T)
                          terra::add(soil_raster) <- rast(soil_pca$map) %>% setNames(soilgriddata[p])
                          rm(soil_subset, soil_pca)
                          }

plot(soil_raster)
res(soil_raster)
rm(soil_raster_tmp, soil_files, soilgriddata, soildepth)


################################################################################:)
### Terrain Metrics #############################################################
################################################################################:)

## SRTM DEM Download

### define download SRTM extent (might not work for large countries)
### trying to ensure that entire AOI is downloaded as SRTM tiles using the corner points and centroid
### alternatively adaption working with SRTM Tiles .shp feasible

SRTM_Download_for_cornerpoints <- function (boundingbox, centroid_coords) {
  SW <- st_point(c(boundingbox@xmin, boundingbox@ymin))
  SE <- st_point(c(boundingbox@xmax, boundingbox@ymin))
  NW <- st_point(c(boundingbox@xmin, boundingbox@ymax))
  NE <- st_point(c(boundingbox@xmax, boundingbox@ymax))
  SRTM <- elevation_3s(lon = centroid_coords[1], lat = centroid_coords[2], path = getwd() )
  SRTM_extent <- as.polygons(ext(SRTM))
  # check if any corner points are not contained and if so download another tile  
  
  if (!st_contains(SW, st_as_sf(SRTM_extent)  , sparse = F)) { 
    SRTM2 <- elevation_3s(lon = boundingbox@xmin, lat = boundingbox@ymin, path = getwd())
    SRTM <- mosaic(SRTM, SRTM2, fun="mean")
    SRTM_extent <- as.polygons(ext(SRTM))
    rm(SRTM2)
  }
  if (!st_contains(SE,st_as_sf(SRTM_extent), sparse = F)) { 
    SRTM2 <- elevation_3s(lon = boundingbox@xmax, lat = boundingbox@ymin, path = getwd())
    SRTM <- mosaic(SRTM, SRTM2, fun="mean")
    SRTM_extent <- as.polygons(ext(SRTM))
    rm(SRTM2)
  } 
  if (!st_contains(NW, st_as_sf(SRTM_extent), sparse = F)) { 
    SRTM2 <- elevation_3s(lon = boundingbox@xmin, lat = boundingbox@ymax, path = getwd())
    SRTM <- mosaic(SRTM, SRTM2, fun="mean")
    SRTM_extent <- as.polygons(ext(SRTM))
    rm(SRTM2)
  }
  if (!st_contains(NE, st_as_sf(SRTM_extent), sparse = F)) { 
    SRTM2 <- elevation_3s(lon = boundingbox@xmax, lat = boundingbox@ymax, path = getwd())
    SRTM <- mosaic(SRTM, SRTM2, fun="mean")
    rm(SRTM2, SRTM_extent)
  }
  return(SRTM)
}

SRTM <- SRTM_Download_for_cornerpoints(boundingbox = extent_ll, centroid_coords = coords)

plot(SRTM)
plot(shapeAOI, add=T)


dir.create("terrain")

# reprojection of DEM and saving
SRTM <- terra::project(SRTM, 
                       soil_raster, 
                       method = "bilinear", 
                       filename = "./terrain/Dem.tif", 
                       overwrite=T ) 
plot(SRTM) 


### calculation of terrain parameters with SAGA #################################

setwd("./terrain")

## This compound function calculates various metrics at once
## For SDM the slope, Topographic Wetness Index (TWI) and the vertical distance to channels is chosen
saga$ta_compound$basic_terrain_analysis(elevation = SRTM,
                                        slope = "slope.tif",
                                        #aspect = "aspect.tif",
                                        #hcurv = "plan_curv.tif",
                                        #vcurv = "profile_curv.tif",
                                        #convergence = "convergence.tif",
                                        #flow = "catchment_area.tif",
                                        wetness = "TWI.tif",
                                        #lsfactor = "LSfactor.tif",
                                        #chnl_base = "channel_base.tif",
                                        #vall_depth = "valley_depth.tif",
                                        #rsp = "relative_slope_position.tif",
                                        chnl_dist = "channel_distance.tif") #vertical distance to channel base


## calculates the Topographic Position Index (TPI)
saga$ta_morphometry$topographic_position_index_tpi(dem= SRTM,
                                                   tpi = "TPI.tif",
                                                   radius_min = 0,
                                                   radius_max = 1000) # alternatives Wind exposition index or topographic openness

## calculates the Topographic Ruggedness Index (TRI)
saga$ta_morphometry$terrain_ruggedness_index_tri(dem = SRTM,
                                                 tri= "TRI.tif")

## caclculates the Sky View Factor (intermediate product)
saga$ta_lighting$sky_view_factor(dem = SRTM,
                                 svf = "Sky_view_factor.tif",
                                 radius = 10000)

## calculates the potential incoming solar radiation based on DEM and Sky view factor for one year
## results in a direct and diffus insolation layer
## processing intensive!!
saga$ta_lighting$potential_incoming_solar_radiation(grd_dem = SRTM,
                                                    grd_svf = "Sky_view_factor.tif",
                                                    grd_direct = "direct_insolation.tif",
                                                    grd_diffus = "diffus_insolation.tif",
                                                    location = 1,
                                                    period = 2,
                                                    day = "2021-01-01",
                                                    day_stop = "2021-12-31",
                                                    days_step = 30
)


## list all the terrain rasters
saga_raster_files <- list.files("./", full.names = T,  pattern="*.tif$") 

### building stack of rasters
terrain_raster <- rast()

for (i in 1:length(saga_raster_files)){
  terra::add(terrain_raster) <- rast(saga_raster_files[i])
}

plot(terrain_raster$direct_insolation)
plot(shapeAOI_UTM, add=T)


#adjust slope to degree  
terrain_raster$slope <- (terrain_raster$slope)*180/pi 

names(terrain_raster)

## drop Sky view factor
terrain_raster <- terra::subset(terrain_raster, "Sky_view_factor", negate=T)

# DEM might be named differently based on the tile
#names(terrain_raster$srtm_39_03) <- c("DEM", names(terrain_raster[[2:8,]]))


setwd("..")
rm(SRTM, SRTM_Download_for_cornerpoints)



#################################################################################:)
### Tree Biomass data ###########################################################
#################################################################################:)
### Download data from CEDA does not work automatically yet. Required is a certificate that can be obtained via a python script. 
### However incorporation in R script did not work out yet. 
### Data available at: https://data.ceda.ac.uk/neodc/esacci/biomass/data/agb/maps/v2.0/geotiff/2018
### citation: Santoro, M.; Cartus, O. (2021): ESA Biomass Climate Change Initiative (Biomass_cci): Global datasets of forest above-ground biomass for the years 2010, 2017 and 2018, v3. NERC EDS Centre for Environmental Data Analysis, 26 November 2021. doi:10.5285/5f331c418e9f4935b8eb1b836f8a91b8

# read in previously (manually) downloaded raster files
biomass_raw <- list.files("./Biomass", full.names = T,  pattern="*.tif$", recursive=T) 

biomass_list <- lapply(1:length(biomass_raw), function(x) {rast(biomass_raw[x]) })

# if required mosaic different tiles
biomass_mosaic <- do.call(terra::mosaic,biomass_list) %>% 
  terra::project(soil_raster, method = "bilinear")
names(biomass_mosaic) <- "biomass"

plot(biomass_mosaic)
plot(shapeAOI_UTM, add=T)
rm(biomass_raw,biomass_list)


#################################################################################:)
### World Clim Data #############################################################
#################################################################################:)

Worldclim <- geodata::worldclim_tile("worldclim", var = "bio", res=0.5, lon=round(coords[1]), lat= round(coords[2])) %>% 
  terra::project(soil_raster, method = "bilinear")

World_clim_names <- c("BIO1_mean_temp", "BIO2_meanDiurnal_temp", "BIO3_iso_Temp", 
                      "BIO4_Seasonality_Temp", "BIO5_max_temp", "BIO6_min_temp", 
                      "BIO7_range_temp", "BIO8_mean_temp_wet", "BIO9_mean_temp_dry", 
                      "BIO10_mean_temp_warm", "BIO11_mean_temp_cold", "BIO12_annPrec", 
                      "BIO13_prec_wetmonth", "BIO14_prec_drymonth", "BIO15_prec_seasonality", 
                      "BIO16_prec_wetquart", "BIO17_prec_dryquart", "BIO18_prec_warmquart", 
                      "BIO19_prec_coldquart")
names(Worldclim) <- World_clim_names

plot(Worldclim$BIO1_mean_temp)
plot(shapeAOI_UTM, add=T)

# conversion of temperature (due to online integer storage convention)

for (i in 1:nlyr(Worldclim)) {
  if (names(Worldclim[[i]]) %in% grep("temp*", World_clim_names, value =T)) {
    Worldclim[[i]] <- Worldclim[[i]]/10
  } else{}
  if (names(Worldclim[[i]]) %in% grep("Temp*", World_clim_names, value =T)) {
    Worldclim[[i]] <- Worldclim[[i]]/100
  }
}
rm(World_clim_names)


# further potential input
### CORINE Landcover ############################################################
### extent only Europe
### load previously downloaded CORINE raster 

Corine <- rast("./U2018_CLC2018_V2020_20u1.tif") 

#reproject
Corine <- terra::project(Corine, soil_raster, method = "near" ) %>% terra::crop(soil_raster)

### categorical raster data
Corine <- as.factor(classify(Corine, matrix(c(44, Inf, NA), ncol = 3, byrow = T)))

names(Corine) <- "landcover2018"



### explanation to classes, apart from that csv incorporated
clc_classes <- read.csv("./clc_legend.csv", header = T)
print(clc_classes)

# convert the RGB color code to Hex color code
colors <- grDevices::rgb(red= as.numeric(str_sub(clc_classes$RGB[1:44], start= 0, end= -9)),
              green= as.numeric(str_sub(clc_classes$RGB[1:44], start= 5, end= -5)),
              blue = as.numeric(str_sub(clc_classes$RGB[1:44], start= 9)), 
              maxColorValue = 255)

# new more compact data frame
cls <- data.frame(id=1:44, cover= clc_classes$LABEL3[1:44], colors)
rm(clc_classes, colors)

#adjust levels of Corine layer for visualization 
levels(Corine) <- cls

# this shows you the correct color code but without labels
terra::plot(Corine, 
            y="cover",
     type="classes", 
     #levels = cls[cls$cover %in% uniqueCorine$cover, 2],
     col= data.frame(value=cls$id, color = cls$colors),
     all_levels=F)

# shows you labels but without color code
terra::plot(Corine,  
            all_levels=F)
terra::polys(vect(shapeAOI_UTM))

writeRaster(Corine, "./Corine_clip.tif")



### masks for displaying
water_mask <- as.polygons(mask(classify(Corine, matrix(c(0, 39, NA,  39, 44, 1,  45, Inf, NA), ncol = 3, byrow=T)), shapeAOI_UTM), dissolve = T)
urban_mask <- as.polygons(mask(classify(Corine, matrix(c(1, 6, 1,  6, 45, NA), ncol = 3, byrow=T)), shapeAOI_UTM), dissolve = T)
glacier_mask <- as.polygons(mask(classify(Corine, matrix(c(0, 33, NA, 33, 34, 1, 34, Inf, NA), ncol = 3, byrow=T)), shapeAOI_UTM), dissolve = T)


plot(soil_raster$bulk_density_0_5)
plot(water_mask, col= 'blue', add= T)
plot(urban_mask, col = 'red', add= T)
plot(glacier_mask, col= "#A6E6CC", add=T)
plot(shapeAOI_UTM, add=T)



################################################################################:)
### MODIS Global Landcover data ################################################
################################################################################:)
# Download from NASA earthdata server (maintenance on wednesday)
dir.create("MODIS_landcover")

bb_utm=c(round(as.numeric(extent_AOI_UTM@xmin), 0), round(as.numeric(extent_AOI_UTM@ymin), 0), 
     round(as.numeric(extent_AOI_UTM@xmax), 0), round(as.numeric(extent_AOI_UTM@ymax), 0)) 

MODIStsp_get_prodnames()
MODIStsp_get_prodlayers("LandCover_Type_Yearly_500m (MCD12Q1)")
MODIStsp(gui = FALSE,
         out_folder = "./MODIS_landcover",
         out_folder_mod = "./MODIS_landcover", 
         selprod = "LandCover_Type_Yearly_500m (MCD12Q1)",
         bandsel = "LC1",
         user = NASA_user,
         password = getPass(),
         start_date  = "2019.01.01", 
         end_date = "2019.12.31",
         spatmeth = "bbox",
         bbox = bb_utm,
         output_proj = crs_UTM,
         resampling = "near",
         out_format = "GTiff")

MODIS_LC <- rast("./MODIS_landcover/LandCover_Type_Yearly_500m_v6/LC1/MCD12Q1_LC1_2019_001.tif") %>% 
  terra::as.factor() %>% 
  terra::project(covariates, method = "near" ) %>% 
    setNames("LC_MODIS")
plot(MODIS_LC)


### adding Landcover classes of the scene as level
level_orig <- levels(MODIS_LC)[[1]]

# creation of temporary dataframe (all possible IDs and landcover names)
ID <- 1:17
landcover <- c( "Evergreen needleleaf forests",
                "Evergreen broadleaf forests",
                "Deciduous needleleaf forests",
                "Deciduous broadleaf forests",
                "Mixed forests",
                "Closed shrublands",
                "Open shrublands",
                "Woody savannas",
                "Savannas",
                "Grasslands",
                "Permanent wetlands",
                "Croplands",
                "Urban and built-up lands",
                "Cropland/natural vegetation mosaics",
                "Snow and ice",
                "Barren",
                "Water bodies")

MODIS_class <- data.frame(ID, landcover)
# check which landcovers are covered by scene
MODIS_class <- subset.data.frame(MODIS_class, MODIS_class$ID %in% level_orig$ID)

# rename to selected landcovers
level_orig$landcover <- MODIS_class$landcover

#adjust MODIS levels
levels(MODIS_LC) <- level_orig
MODIS_LC
rm(ID, landcover, MODIS_class, level_orig)
levels(MODIS_LC)

# create landcover masks
water_mask_MODIS <-  as.polygons(mask(MODIS_LC[[1]] == 17, shapeAOI_UTM), dissolve = T)
plot(water_mask_MODIS)

urban_mask_MODIS <- as.polygons(mask(MODIS_LC[[1]] == 13, shapeAOI_UTM), dissolve = T)
plot(urban_mask_MODIS)


#################################################################################:)
### Terra Evopotranspiration yearly #############################################

dir.create("./Evapotranspiration")
MODIStsp(gui = FALSE,
         out_folder = "./Evapotranspiration",
         out_folder_mod = "./Evapotranspiration", 
         selprod = "Net_ETgf_Yearly_500m (M*D16A3GF)"  ,
         bandsel = "ET_500m",
         sensor = "Terra",
         user = NASA_user,
         password = getPass(),
         start_date  = "2000.01.01", 
         end_date = "2021.12.31",
         spatmeth = "bbox",
         bbox = bb_utm,
         output_proj = crs_UTM,
         resampling = "bilinear",
         out_format = "GTiff")

# if only one year:
# MODIS_Evp <- rast("./Evapotranspiration/Net_ETgf_Yearly_500m_v6/ET_500m/MOD16A3GF_ET_500m_2019_001.tif") %>% 
#   classify(matrix(c(32761, Inf, NA), ncol = 3, byrow = T)) %>%
#   terra::project(soil_raster, method = "bilinear" ) %>% 
#   setNames("Evapotranspiration")

# for multiple years:
# does not account for mosaicing
MODIS_EVP_files <- list.files("./Evapotranspiration/Net_ETgf_Yearly_500m_v6/ET_500m/", pattern = "*.tif$", full.names = T)

for (i in 1:length(MODIS_EVP_files)) {
  if (i == 1){
    MODIS_EVP_stack <- rast(MODIS_EVP_files[i]) %>% 
      terra::project(soil_raster, method="bilinear") 
    
  }
  else {
    add(MODIS_EVP_stack) <- rast(MODIS_EVP_files[i]) %>% 
      terra::project(soil_raster, method="bilinear") 
   }
}

#calculate the mean over all the years
EVP <- terra::app(MODIS_EVP_stack, fun= "mean", filename="EVP_MODIS.tif", overwrite=T ) %>% setNames("EVP")
plot(EVP)
rm(MODIS_EVP_files, MODIS_EVP_stack)

#################################################################################:)
### Terra NPP yearly ############################################################
### however large gaps due to missing LAI values in processing

dir.create("./NPP")
MODIStsp(gui = F,
         out_folder = "./NPP",
         out_folder_mod = "./NPP", 
         selprod = "Net_PP_GapFil_Yearly_500m (M*D17A3HGF)"  ,
         bandsel = "Npp",
         sensor = "Terra",
         user = NASA_user,
         password = getPass(),
         start_date  = "2000.01.01", 
         end_date = "2021.12.31",
         spatmeth = "bbox",
         bbox = bb_utm,
         output_proj = crs_UTM,
         resampling = "near",
         nodata_change = F,
         out_format = "GTiff")


# if only one year
# MODIS_Npp <- rast("./NPP/Net_PP_Yearly_500m_v6/Npp/MOD17A3HGF_Npp_2019_001.tif") %>% 
#   project(soil_raster, method = "bilinear" ) %>% 
#   terra::classify(matrix(c(32761, Inf, NA), ncol = 3, byrow = T)) %>% 
#   setNames("NPP")

MODIS_NPP_files <- list.files("./NPP/Net_PP_Yearly_500m_v6/Npp/", pattern = "*.tif$", full.names = T)

for (i in 1:length(MODIS_NPP_files)) {
  if (i == 1){
    MODIS_raster <- rast(MODIS_NPP_files[i]) %>% 
      terra::classify(matrix(c(32700, Inf, NA), ncol = 3, byrow = T)) %>% 
      terra::project(covariates2, method="bilinear") 
      
        }
  else {
    add(MODIS_raster) <- rast(MODIS_NPP_files[i]) %>% 
     terra::classify(matrix(c(32700, Inf, NA), ncol = 3, byrow = T)) %>% 
     terra::project(covariates2, method="bilinear") 
      
  }
}
# https://lpdaac.usgs.gov/products/mod17a3hv006/ fill value classes

## calculate the mean over all years
NPP <- terra::app(MODIS_raster, fun= "mean", filename="NPP_MODIS.tif", overwrite=T ) %>% setNames("NPP")

plot(NPP)

plot(MODIS_raster$MOD17A3HGF_Npp_2001_001)
rm(MODIS_NPP_files, MODIS_raster)

#################################################################################:)
### NDVI Landsat 8 via Google Earth Engine ######################################
#################################################################################:)

# This function gets NDVI from Landsat 8 imagery.
addNDVI <- function(image) {
  return(image$addBands(image$normalizedDifference(c("B5", "B4"))$rename("NDVI")$float()))
}

# Landsat 8 cloud, cloud shadow and snow mask function
cloudMaskL8 <- function(image){
  qa_band <- image$select("pixel_qa")
  cloud <-qa_band$bitwiseAnd(56)$eq(0)
  return(image$updateMask(cloud))
}


## GEE Asset creation, only once needed

# shape_GEE <- sf_as_ee(st_buffer(shapeAOI, 10000), 
#                       via = "getInfo_to_asset",
#                       assetId = sprintf("%s/%s", ee_get_assethome(), 'rgee_SDM_buffer'))

## load Asset
shape_GEE <- ee$FeatureCollection(sprintf("%s/%s", ee_get_assethome(), 'rgee_SDM_buffer')) 



## Landsat 8 ImageCollection with NDVI, cloud and cloud shadow filtered
L8Image <- ee$ImageCollection('LANDSAT/LC08/C01/T1_SR')$
  filterBounds(shape_GEE)$
  filterDate("2014-01-01", "2022-01-01")$
  filter(ee$Filter$calendarRange(5,9, 'month'))$
  map(addNDVI)$
  map(cloudMaskL8)$
  map(function(image) {return(image$clip(shape_GEE))})$
  select("NDVI")

## create median and standard deviation composites
L8_NDVI_median <- L8Image$reduce(ee$Reducer$median())

L8_NDVI_stdv <- L8Image$reduce(ee$Reducer$stdDev())

#interpolation of masked out gaps
interpolation_gaps <-  function(image){
  focal_min <-  image$focalMin(radius= 11, kernelType= 'square', units= 'pixels')
  gaps <-  image$unmask()$clip(shape_GEE)$eq(0)
  gaps_values <-  focal_min$updateMask(gaps)
  image_interpolated  <-  image$blend(gaps_values)
  return(image_interpolated) 
}

L8_NDVI_median_interpolated <- interpolation_gaps(L8_NDVI_median)
L8_NDVI_stdv_interpolated <- interpolation_gaps(L8_NDVI_stdv)

## add as leaflet map viewer
Map$centerObject(shape_GEE)
Map$addLayer(
  eeObject = L8_NDVI_median_interpolated,
  visParams = list(min = -1, max= 1, palette = c('brown', 'khaki', 'green')),
  name = "median_NDVI"
)

Map$addLayer(
  eeObject = L8_NDVI_stdv_interpolated,
  visParams = list(min = -0.5, max= 0.5, palette = c('brown', 'grey', 'green')),
  name = "stdv_NDVI"
)

# export and conversion into R raster, takes some time
# create a google drive folder named GEE
NDVI_median <- ee_as_raster(L8_NDVI_median_interpolated,
                            via = "drive",
                            region= shape_GEE$first()$geometry(),
                            dsn = "L8_median_NDVI.tif",
                            container = "GEE",
                            scale = 250,
                            maxPixels = 1e+12)

NDVI_stdv <- ee_as_raster(L8_NDVI_stdv_interpolated,
                            via = "drive",
                            region= shape_GEE$first()$geometry(),
                            dsn = "L8_stdv_NDVI",
                            container = "GEE",
                            scale = 250,
                            maxPixels = 1e+12)

rm(L8_NDVI_median_interpolated, L8_NDVI_stdv_interpolated, interpolation_gaps, addNDVI, 
   L8_NDVI_median, L8_NDVI_stdv, L8Image, cloudMaskL8, shape_GEE)

# projection
NDVI_median <- rast(NDVI_median) %>% terra::project(soil_raster, method="bilinear")
NDVI_stdv <- rast(NDVI_stdv) %>% terra::project(soil_raster, method="bilinear")

# manual data loading
#NDVI_median <- rast("./VI/L8_median_NDVI_2022_09_27_11_26_05.tif") %>% terra::project(soil_raster, method="bilinear")
#NDVI_stdv <- rast("./VI/L8_stdv_NDVI_2022_09_27_11_36_57.tif") %>% terra::project(soil_raster, method="bilinear")


plot(NDVI_stdv)
plot(NDVI_median)
plot(shapeAOI_UTM, add=T)

#----------------------------------------------------------------------------------:)
### multicollinearity ##########################################################


### correlation between grouped rasters

corrplot_calculation <- function (raster_stack) {
  cor <- layerCor(raster_stack, 'pearson', na.rm=T)
  cor <- cor[[1]]
  cor[which(cor>1)] <- 1
  cor[which(cor<(-1))] <- -1
  stackcor <- corrplot(cor, method = "circle", type = "lower", order = "FPC", col = rev.default(COL2('RdBu', 200)), tl.col = "#010101", tl.srt = 45)
  return(stackcor)
}


terrain_corplot <- corrplot_calculation(terrain_raster)
soil_corplot <- corrplot_calculation(soil_raster)
climate_corplot <- corrplot_calculation(Worldclim)

# exclude MODIS NPP, EVP and Landcover

#----------------------------------------------------------------------------------:)
################################################################################:)
# stacking and scaling ##########################################################

covariates <- terra::rast(list(soil_raster, terrain_raster, biomass_mosaic, Worldclim, NDVI_median, NDVI_stdv))
plot(covariates$bdod)

# center and scale
covariates_scaled <- terra::scale(covariates, center=T, scale = T)
plot(covariates_scaled$H3O)
summary(covariates_scaled$H3O)

### calculate variance inflation factor with flexsdm package
### to reduce multicollinearity
vif <- correct_colinvar(env_layer = covariates_scaled, method = c("vif", th= "10"))

vif$env_layer
vif$removed_variables
view(vif$vif_table)    

### subset covariates to layers with vif value < 10
covariates <- vif$env_layer
covariates <- terra::rast(list(covariates)) # necessary conversion

### crop and mask to a buffer mask
buffer_shape <- st_buffer(shapeAOI_UTM, 10000)
covariates <- terra::crop(covariates2, y=terra::ext(vect(buffer_shape)), snap='near', mask=T)
covariates  <- terra::mask(covariates, mask=vect(buffer_shape))

### save covariates stack as .grd
dir.create("covariates")
terra::writeRaster(covariates,"./covariates/predictors.grd",  gdal="ENVI")

