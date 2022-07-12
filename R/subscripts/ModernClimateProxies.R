# Load libraries ----------------------------------------------------------
library(raster)
library(sf)
# Coral Reefs -------------------------------------------------------------
# Mean Annual Sea Surface Temperature
# Source:
# https://www.bio-oracle.org/index.php
sst <- raster("./data/raw/climate/BioOracle_20220711/Present.Surface.Temperature.Mean.asc")
# Generate shallow marine mask for coral reefs
# ETOPO5 world digital elevation model
# Source:
# https://www.eea.europa.eu/data-and-maps/data/world-digital-elevation-model-etopo5
dem <- raster("./data/raw/DEM/alwdgg.tif")
dem[dem > 0] <- NA
dem <- abs(dem)
dem[dem > 200] <- NA
#mask sst layer by available shallow marine habitat
sst <- mask(x = sst, mask = dem)
vals_sst <- getValues(sst)
hist(vals_sst, main = "Available sea surface temperature (continental shelf)",
     xlab = "Mean sea surface temperature (\u00b0C)",
     col=rgb(1,0,0,0.5))

# Load coral reef data 
# Source: https://data.unep-wcmc.org/datasets/1
# Extract coordinates and unique cell values within raster
reef <- sf::read_sf(
  dsn = "./data/raw/occurrences/coral_20220711/14_001_WCMC008_CoralReefs2018_v4_1/01_Data/WCMC008_CoralReef2018_Py_v4_1.shp")
reef <- sf::st_coordinates(reef)
cells <- unique(raster::cellFromXY(object = sst, xy = reef))
reef <- xyFromCell(object = sst, cell = cells)
vals_reef <- sst[cells]
hist(vals_reef, main = "Warm-water coral reefs",
     xlab = "Mean sea surface temperature (\u00b0C)",
     col=rgb(0,0,1,0.5))
reef <- cbind(reef, vals_reef)
colnames(reef) <- c("lng", "lat", "sst")
reef <- na.omit(reef)
saveRDS(reef, "./data/processed/modern_reef_temp.RDS")

# Mangroves ---------------------------------------------------------------
# Mean Annual Sea Surface Temperature
# Source:
# https://www.bio-oracle.org/index.php
sst <- raster("./data/raw/climate/BioOracle_20220711/Present.Surface.Temperature.Mean.asc")
# Generate shallow marine mask for coral reefs
# ETOPO5 world digital elevation model
# Source:
# https://www.eea.europa.eu/data-and-maps/data/world-digital-elevation-model-etopo5
dem <- raster("./data/raw/DEM/alwdgg.tif")
dem[dem > 10] <- NA
dem <- abs(dem)
dem[dem > 10] <- NA
#mask sst layer by available shallow marine habitat
sst <- mask(x = sst, mask = dem)
vals_sst <- getValues(sst)
hist(vals_sst, main = "Available sea surface temperature (coastlines)",
     xlab = "Mean sea surface temperature (\u00b0C)",
     col=rgb(1,0,0,0.5))

# Load mangrove data 
# Source: https://data.unep-wcmc.org/datasets/5
# Extract coordinates and unique cell values within raster
mangrove <- sf::read_sf(
  dsn = "./data/raw/occurrences/mangrove_20220711/WCMC011_AtlasMangrove2010_v3_1/01_Data/WCMC011_AtlasMangrove2010_Py_v3_1.shp")
mangrove <- sf::st_coordinates(mangrove)
cells <- unique(raster::cellFromXY(object = sst, xy = mangrove))
mangrove <- xyFromCell(object = sst, cell = cells)
vals_mangrove <- sst[cells]
hist(vals_mangrove, main = "Mangroves",
     xlab = "Mean sea surface temperature (\u00b0C)",
     col=rgb(0,1,0,0.5))
mangrove <- cbind(mangrove, vals_mangrove)
colnames(mangrove) <- c("lng", "lat", "sst")
mangrove <- na.omit(mangrove)
saveRDS(mangrove, "./data/processed/modern_mangrove_temp.RDS")
