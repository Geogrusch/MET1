# MET1

###Stacked species distribution modeling to model species richness <br />

Here exemplified for 40 tree species in Austria.<br />


The first script downloads and generates many environmental raster layers as covariates/predictors
Including: <br />

  ISRIC soildata<br />
  SRTM DEM and Saga based terrain options (Rsagacmd package)<br />
  Biomass (ESA CCI product has to be downloaded at https://data.ceda.ac.uk/neodc/esacci/biomass/data/agb/maps/v2.0/geotiff/2018)<br />
  Worldclim data<br />
  Landsat8 NDVI via rgee package<br />
  MODIS data via MODIStsp package<br />
  Corine landcover (download at: https://land.copernicus.eu/pan-european/corine-land-cover/clc2018?tab=download)<br />

Raster layers are stacked &
filtered to reduce collinearity 


Script 2: <br />

  downloads species presence data from gbif.org<br />
  cleans data <br />
  applies sampling bias reduction<br />
  creates pseudo absence data via flexsdm package with distance buffer and BIOCLIM modell<br />
  generates ensemble models for every species<br />
  validates models<br />
  stacks ensemble models for species richness map
