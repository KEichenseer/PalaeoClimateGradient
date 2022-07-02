setwd("C://Users/keichenseer/OneDrive - University of Plymouth/data/reefzones")
setwd("D://OneDrive - Durham University/Plymouth/data/reefzones")

reef <- read.csv("pared_2021_06_12.csv") # reed the PARED data
reef_text <- read.csv("pared_2021_06_12_text.csv") # reed the PARED data

timebins <- read.csv(file="clean_timebins.csv", header=TRUE, sep=",") # read Phanerozoic international stratigraphy (Gradstein 2016)
other_intervals <- read.csv(file="other_intervals.csv", header=TRUE, sep = ";") # read list of other stratigraphic intervals
trim.leading <- function (x)  sub("^\\s+", "", x) # get rid of leading spaces
other_intervals$early_stage <- trim.leading(other_intervals$early_stage)
other_intervals$late_stage <- trim.leading(other_intervals$late_stage)

length(table(reef$intervall))

reef$intervall <- gsub("\\?","",reef$intervall) # remove the question marks from the intervall column
reef$intervall <- gsub("\\(","",reef$intervall) # remove the parentheses from the intervall column
reef$intervall <- gsub("\\)","",reef$intervall) # # remove the parentheses from the intervall column

reef$intervall <- gsub("\\/","-",reef$intervall) # transform "/" into "-"
reef$intervall <- gsub("\\ to ","-",reef$intervall) # transform " to " into "-"
reef$intervall <- gsub("\\-"," ",reef$intervall) # transform "-" to " "
reef$intervall <- gsub("\\  "," ",reef$intervall) # transform "  " to " "
reef$intervall <- gsub("\\  "," ",reef$intervall) # transform "   " to " "
reef$intervall <- gsub("\\Upper","upper",reef$intervall) # standardise capitalisation
reef$intervall <- gsub("\\Middle","middle",reef$intervall) # standardise capitalisation
reef$intervall <- gsub("\\Lower","lower",reef$intervall) # standardise capitalisation
reef$intervall <- gsub("\\Late","late",reef$intervall) # standardise capitalisation
reef$intervall <- gsub("\\Early","early",reef$intervall) # standardise capitalisation

reef$intervall[which(reef$intervall == "middle upper Viséan")] <- "middle upper Visean" # correct error in "Visean"
reef$intervall[which(reef$intervall == "upper Viséan")] <- "upper Visean"
length(table(reef$intervall))

# Seperate out multiple intervals / words in the Intervall column
timetab <- matrix(NA,nrow = nrow(reef), ncol = 12)
for (i in 1:nrow(reef)) {
  temp <- unlist(strsplit(as.character(reef$intervall[i]), " "))
  timetab[i,1:length(temp)] <- temp
} 

timetab2 <- as.data.frame(matrix(NA,nrow = nrow(reef), ncol = 13))
colnames(timetab2) <- c("specif1_1","specif1_2","early","specif2_1", "specif2_2" ,"late", "early_int", "late_int", "early_int2", "late_int2", "early_stage", "late_stage", "bin")

# Function that makes first letter upper case:
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
specifiers <- c("early", "middle", "late", "upper", "lower", "latest", "earliest","topmost", "mid", "lowermost", "pre")
for (i in 1:nrow(reef)) {
  if (firstup(timetab[i,1]) == timetab[i,1]) { timetab2$early[i] <- timetab[i,1] # e
  if (is.na(timetab[i,2])) next
  if (firstup(timetab[i,2]) == timetab[i,2])  timetab2$late[i] <- timetab[i,2] # e l 
  if (timetab[i,2] %in% specifiers) { timetab2$specif2_1[i] <- timetab[i,2] # e s
  if (is.na(timetab[i,3])) next
  if (firstup(timetab[i,3]) == timetab[i,3])  timetab2$late[i] <- timetab[i,3] # e s l 
  if (timetab[i,3] %in% specifiers) { timetab2$specif2_2[i] <- timetab[i,3] # e s s
  if (is.na(timetab[i,4])) next
  if (firstup(timetab[i,4]) == timetab[i,4])  timetab2$late[i] <- timetab[i,4] # e s s l
  }
  }
  }
  if (timetab[i,1] %in% specifiers) { timetab2$specif1_1[i] <- timetab[i,1] # s
  if (timetab[i,1] %in% specifiers & firstup(timetab[i,2]) == timetab[i,2]) { timetab2$early[i] <- timetab[i,2] # s e
  if (is.na(timetab[i,3])) next
  if (firstup(timetab[i,3]) == timetab[i,3])  timetab2$late[i] <- timetab[i,3] # s e l
  if (timetab[i,3] %in% specifiers) { timetab2$specif2_1[i] <-  timetab[i,3] # s e s
  if (is.na(timetab[i,4])) next
  if (firstup(timetab[i,4]) == timetab[i,4])  timetab2$late[i] <- timetab[i,4] # s e s l 
  if (timetab[i,4] %in% specifiers) { timetab2$specif2_2[i] <- timetab[i,4] # s e s s
  if (firstupp(timetab[i,5]) == timetab[i,5])  timetab2$late[i] <- timetab[i,5] # s e s s l
  }
  
  }
  }
  if (timetab[i,1] %in% specifiers & timetab[i,2] %in% specifiers) { timetab2$specif1_2[i] <- timetab[i,2] # s s
  if (firstup(timetab[i,3]) == timetab[i,3])  {timetab2$early[i] <- timetab[i,3] # s s e
  if (is.na(timetab[i,4])) next
  if (firstup(timetab[i,4]) == timetab[i,4])  timetab2$late[i] <- timetab[i,4] # s s e l
  if (is.na(timetab[i,5])) next
  if (timetab[i,4] %in% specifiers) { timetab2$specif2_1[i] <- timetab[i,4] # s s e s
  if (firstup(timetab[i,5]) == timetab[i,5])  timetab2$late[i] <- timetab[i,5] # s s e s l
  if (is.na(timetab[i,6])) next
  if (timetab[i,5] %in% specifiers) { timetab2$specif2_2[i] <- timetab[i,5] # s s e s s
  if (firstupp(timetab[i,6]) == timetab[i,6])  timetab2$late[i] <- timetab[i,6] # s e s s l
  }
  }
  }
  }
  }
}


# Corrections to timetab2 
timetab2$early[which(timetab2$early == "Westphal")] <- "Westphalian"
timetab2$late[which(timetab2$late %in% c(1:100,"A","B", "C", "D", "E", "F", "G", "H", "I", "J"))] <- timetab2$early[
  which(timetab2$late %in% c(1:100,"A","B", "C", "D", "E", "F", "G", "H", "I", "J"))] # correct assignments of units with numbers or letters
timetab2$early[which(timetab2$early == "Llandoverian")] <- "Llandovery"
timetab2$early[which(timetab2$early == "Ludlov")] <- "Ludlow"
timetab2$early[which(timetab2$early == "Ludlovian")] <- "Ludlow"
timetab2$late[which(timetab2$late == "Ludlovian")] <- "Ludlow"
timetab2$early[which(timetab2$early == "Moskovian")] <- "Moscovian"
timetab2$early[which(timetab2$early == "Pridolian")] <- "Pridoli"
timetab2$late[which(timetab2$late == "Pridolian")] <- "Pridoli"
timetab2$early[which(timetab2$early == "Titonian")] <- "Tithonian"
timetab2$early[which(timetab2$early == "Tournasian")] <- "Tournaisian"
timetab2$early[which(timetab2$early == "Wenlockian")] <- "Wenlock"
timetab2$late[which(timetab2$late == "Tn3")] <- "Osagean"
timetab2$late[which(timetab2$late == "Portlandian")] <- "Tithonian"
timetab2$early[which(timetab2$early == "Wenlockian")] <- "Wenlock"
timetab2$early[which(timetab2$early == "VisÃ©an")] <- "Visean"
timetab2$late[which(timetab2$late == "VisÃ©an")] <- "Visean"

timetab2$specif2_1[which(is.na(timetab2$late) & is.na(timetab2$specif2_1))] <- timetab2$specif1_1[which(is.na(timetab2$late) & is.na(timetab2$specif2_1))]
timetab2$specif2_2[which(is.na(timetab2$late) & is.na(timetab2$specif2_2))] <- timetab2$specif1_2[which(is.na(timetab2$late) & is.na(timetab2$specif2_2))]
timetab2$late[which(is.na(timetab2$late))] <- timetab2$early[which(is.na(timetab2$late))]

head(timetab)

timetab2[timetab2=="lower"] <- "early"
timetab2[timetab2=="upper"] <- "late"

# recompose to match other intervals

for(i in 1:nrow(reef)) {
  
  if(is.na(timetab2$specif1_1[i]) & is.na(timetab2$specif1_2[i])) timetab2$early_int[i] <- timetab2$early[i]
  if(!(is.na(timetab2$specif1_1[i])) & is.na(timetab2$specif1_2[i])) timetab2$early_int[i] <- paste((timetab2$specif1_1[i]), timetab2$early[i], sep = " ")
  if(!(is.na(timetab2$specif1_1[i])) & !(is.na(timetab2$specif1_2[i]))) timetab2$early_int[i] <- paste(timetab2$specif1_1[i], timetab2$specif1_2[i], timetab2$early[i], sep = " ")

  if(is.na(timetab2$specif2_1[i]) & is.na(timetab2$specif2_2[i])) timetab2$late_int[i] <- timetab2$late[i]
  if(!(is.na(timetab2$specif2_1[i])) & is.na(timetab2$specif2_2[i])) timetab2$late_int[i] <- paste((timetab2$specif2_1[i]), timetab2$late[i], sep = " ")
  if(!(is.na(timetab2$specif2_1[i])) & !(is.na(timetab2$specif2_2[i]))) timetab2$late_int[i] <- paste(timetab2$specif2_1[i], timetab2$specif2_2[i], timetab2$late[i], sep = " ")

}


for (i in 1: nrow(reef)) {
  # match and assign stages
  if (timetab2$early[i] %in% timebins$stage) timetab2$early_stage[i] <- as.character(timebins$stage[which(timetab2$early[i] == timebins$stage)])
  if (timetab2$late[i] %in% timebins$stage) timetab2$late_stage[i] <- as.character(timebins$stage[which(timetab2$late[i] == timebins$stage)])
}


for (i in 1: nrow(reef)) {
  # match and assign stages
  print(i)
  if (timetab2$early_int[i] %in% other_intervals$ori_name & is.na(timetab2$early_stage[i])) timetab2$early_stage[i] <- 
      unique(as.character(other_intervals$early_stage[which(timetab2$early_int[i] == other_intervals$ori_name)]))
  
  if (timetab2$late_int[i] %in% other_intervals$ori_name& is.na(timetab2$late_stage[i])) timetab2$late_stage[i] <- 
      unique(as.character(other_intervals$late_stage[which(timetab2$late_int[i] == other_intervals$ori_name)]))
}

table(timetab2$early_int[which(is.na(timetab2$early_stage))])
table(timetab2$late_int[which(is.na(timetab2$late_int))])

reef <- cbind(timetab2$early_stage, timetab2$late_stage, reef)
colnames(reef)[1:2] <- c("early_stage", "late_stage")


reef$early_stage[which(reef$early_stage == "Early Darriwilian")] <- "Darriwilian"
reef$early_stage[which(reef$early_stage == "Late Darriwilian")] <- "Darriwilian"
reef$early_stage[which(reef$early_stage == "Early Katian")] <- "Katian"
reef$early_stage[which(reef$early_stage == "Late Katian")] <- "Katian"
reef$early_stage[which(reef$early_stage == "Early Visean")] <- "Visean"
reef$early_stage[which(reef$early_stage == "Late Visean")] <- "Visean"

reef$late_stage[which(reef$late_stage == "Early Darriwilian")] <- "Darriwilian"
reef$late_stage[which(reef$late_stage == "Late Darriwilian")] <- "Darriwilian"
reef$late_stage[which(reef$late_stage == "Early Katian")] <- "Katian"
reef$late_stage[which(reef$late_stage == "Late Katian")] <- "Katian"
reef$late_stage[which(reef$late_stage == "Early Visean")] <- "Visean"
reef$late_stage[which(reef$late_stage == "Late Visean")] <- "Visean"


reef$early_stage[which(reef$early_stage=="Stage 5")] <- "Wuliuan"

allstages <- timebins$stage[-c(14,17,36)]
allstages[5] <- "Wuliuan"
allstages[14] <- "Darriwilian"
allstages[16] <- "Katian"
allstages[34] <- "Visean"

table(reef$late_stage[which(!(reef$late_stage %in% allstages))])
table(reef$early_stage[which(!(reef$early_stage %in% allstages))])

table(reef$early_stage)


# Check if everything is ok
stageage_early <- rep(NA,nrow(reef))
for (i in 1:nrow(reef)) if(reef$early_stage[i] %in% timebins$stage) stageage_early[i] <-  timebins$age[which(timebins$stage == reef$early_stage[i])]

stageage_late <- rep(NA,nrow(reef))
for (i in 1:nrow(reef)) if(reef$late_stage[i] %in% timebins$stage) stageage_late[i] <-  timebins$age[which(timebins$stage == reef$late_stage[i])+1]


stageage_early_diff <- sapply(reef$early_stage, function(x) timebins$age[which(timebins$stage == x)])

length(stageage_early)

plot(reef$max_ma,stageage_early, col = "black") #, ylim = c(-25,25))
points(reef$min_ma,stageage_late, col = "red")

View(cbind(reef$intervall,reef$max_ma,  reef$early_stage, stageage_early,  reef$min_ma,reef$late_stage,stageage_late))

points(reef$min_ma,stageage_late, col = "red")

View(reef[which(log(reef$max_ma)-log(stageage_early)>1),])
View(reef[which((reef$min_ma)-stageage_late>10),])

View(reef[which(is.na(reef$late_stage)),])
View(cbind(stageage_late[which((reef$min_ma)-stageage_late < -10 )], reef[which((reef$min_ma)-stageage_late < -10),]))

table(reef$biota_main)

tropcoral <- subset(reef, biota_main == 1 & tropical == 0 & bathymetry == 1 & early_stage == late_stage)

stages <- c("Pre-Cambrian boundary", "Fortunian", "Stage 2", "Stage 3", "Stage 4", "Wuliuan", "Drumian", 
            "Guzhangian", "Paibian", "Jiangshanian", "Stage 10", "Tremadocian", 
            "Floian", "Dapingian", "Darriwilian", "Sandbian", "Katian", "Hirnantian", 
            "Rhuddanian", "Aeronian", "Telychian", "Sheinwoodian", "Homerian", 
            "Gorstian", "Ludfordian", "Pridoli", "Lochkovian", "Pragian", 
            "Emsian", "Eifelian", "Givetian", "Frasnian", "Famennian", "Tournaisian", 
            "Visean", "Serpukhovian", "Bashkirian", "Moscovian", "Kasimovian", 
            "Gzhelian", "Asselian", "Sakmarian", "Artinskian", "Kungurian", 
            "Roadian", "Wordian", "Capitanian", "Wuchiapingian", "Changhsingian", 
            "Induan", "Olenekian", "Anisian", "Ladinian", "Carnian", "Norian", 
            "Rhaetian", "Hettangian", "Sinemurian", "Pliensbachian", "Toarcian", 
            "Aalenian", "Bajocian", "Bathonian", "Callovian", "Oxfordian", 
            "Kimmeridgian", "Tithonian", "Berriasian", "Valanginian", "Hauterivian", 
            "Barremian", "Aptian", "Albian", "Cenomanian", "Turonian", "Coniacian", 
            "Santonian", "Campanian", "Maastrichtian", "Danian", "Selandian", 
            "Thanetian", "Ypresian", "Lutetian", "Bartonian", "Priabonian", 
            "Rupelian", "Chattian", "Aquitanian", "Burdigalian", "Langhian", 
            "Serravallian", "Tortonian", "Messinian", "Zanclean", "Piacenzian", 
            "Gelasian", "Calabrian", "Chibanian", "Upper Pleistocene", "Holocene")
unistagecor <- unique(tropcoral$early_stage)
unistagecor[which(!(unistagecor %in% stages))] # Ludlow Devonian Modern - not use
tropcoral$early_stage[which(tropcoral$early_stage == "Late Pleistocene")] <- "Upper Pleistocene"
tropcoral$late_stage[which(tropcoral$late_stage == "Late Pleistocene")] <- "Upper Pleistocene"
tropcoral$early_stage[which(tropcoral$early_stage == "Middle Pleistocene")] <- "Chibanian"
tropcoral$late_stage[which(tropcoral$late_stage == "Middle Pleistocene")] <- "Chibanian"
saveRDS(tropcoral,"data/processed/PARED_coralreefs_processed_02_07_2022.rds")
