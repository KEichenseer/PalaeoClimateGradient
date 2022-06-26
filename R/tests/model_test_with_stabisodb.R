### A test with StabisoDB data
### 26/06/2022
###

dat <- readRDS("data/processed/StabisoDB_processed_26_06_2022.rds")

tith <- subset(dat,stage_2020 == "Tithonian")#table(dat$stage_2020)
plot(tith$latitude,tith$temperature)

# wait a second - how do we deal with situationsn with 1 value per locality?