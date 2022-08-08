#run-analyses

data <- read.csv("./data/raw/Eocene/bio_proxies.csv")
data

saveRDS("./data/processed/bio_proxies_2022_08_08.RDS")