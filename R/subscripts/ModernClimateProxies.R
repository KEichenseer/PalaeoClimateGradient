library(raster)
library(sf)

# Mean Annual Sea Surface Temperature
# Source: https://www.bio-oracle.org/index.php
sst <- raster("./data/raw/climate/BioOracle_20220711/Present.Surface.Temperature.Mean.asc")
vals <- getValues(sst)
hist(vals, main = "Available sea surface temperature",
     xlab = "Mean sea surface temperature (\u00b0C)")

# Load coral reef data 
# Source: https://data.unep-wcmc.org/datasets/1
# Extract coordinates and unique cell values within raster
reef <- sf::read_sf(
  dsn = "./data/raw/occurrences/coral_20220711/14_001_WCMC008_CoralReefs2018_v4_1/01_Data/WCMC008_CoralReef2018_Py_v4_1.shp")
reef <- sf::st_coordinates(reef)
cells <- unique(raster::cellFromXY(object = sst, xy = reef))
reef <- xyFromCell(object = sst, cell = cells)
vals <- sst[cells]
hist(vals, main = "Warm-water coral reefs",
     xlab = "Mean sea surface temperature (\u00b0C)")
reef <- cbind(reef, vals)
colnames(reef) <- c("lng", "lat", "sst")
reef <- na.omit(reef)
saveRDS(reef, "./data/processed/modern_reef_temp.RDS")

# Load mangrove data 
# Source: https://data.unep-wcmc.org/datasets/5
# Extract coordinates and unique cell values within raster
mangrove <- sf::read_sf(
  dsn = "./data/raw/occurrences/mangrove_20220711/WCMC011_AtlasMangrove2010_v3_1/01_Data/WCMC011_AtlasMangrove2010_Py_v3_1.shp")
mangrove <- sf::st_coordinates(mangrove)
cells <- unique(raster::cellFromXY(object = sst, xy = mangrove))
mangrove <- xyFromCell(object = sst, cell = cells)
vals <- sst[cells]
hist(vals, main = "Mangroves",
     xlab = "Mean sea surface temperature (\u00b0C)")
mangrove <- cbind(mangrove, vals)
colnames(mangrove) <- c("lng", "lat", "sst")
mangrove <- na.omit(mangrove)
saveRDS(mangrove, "./data/processed/modern_mangrove_temp.RDS")
