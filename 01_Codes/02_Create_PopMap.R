# Info --------------------------------------------------------------------
# 
# Authors: Audrey Bourret, Luca Montana
# Affiliation: Fisheries and Oceans Canada (DFO)
# Group: Genomic laboratory, Demersal and Benthic Sciences Branch 
# Location: Maurice Lamontagne Institute
# Date: 2022-03-09
# 
# Overview: Creation of POP INFO files + barcodes for RAD-seq pipeline with STACKS
# 
#


# Library -----------------------------------------------------------------

library(readxl)
library(tidyverse)

# Internal functions
for(i in 1:length( list.files("./01_Codes/Functions") )){
  source(file.path("./01_Codes/Functions",  list.files("./01_Codes/Functions")[i]))  
}


# Data --------------------------------------------------------------------

list.files(get.value("info.path"))

IBIS.barcode <-  read_excel(file.path(get.value("info.path"),"BarcodesIBIS.xlsx"))
IBIS.barcode

# Upload specimens sent to Novaseq ddRAD analyses 2019
GQ.data <- read.csv(file.path(get.value("access.path"),"MOBELS_ddRAD811_230717.csv"), stringsAsFactors = F)
GQ.data <- GQ.data %>% mutate(Cat_sample = ifelse(duplicated(Numero_unique_specimen ), "Duplicate", "Sample"),
                              ID_GQ = ifelse(Cat_sample == "Duplicate", paste0(Numero_unique_specimen, "_rep"), Numero_unique_specimen))


# Create a new dataset
data <- GQ.data %>% select(ID_GQ, Numero_unique_specimen, Numero_reception_specimen, Numero_unique_extrait, No_plaque_envoi, No_puits_envoi, Cat_sample, Genre, Age, 
                           Sexe_visuel, Sexe_laboratoire, Region_echantillonnage, Lieu_echantillonnage, Latitude_echantillonnage_DD, Longitude_echantillonnage_DD,
                           Annee_echantillonnage,Mois_echantillonnage,Jour_echantillonnage, Communaute, Nom_chasseur, Haplotype) %>%
  left_join(IBIS.barcode %>% select(No_puits_envoi = Cell2, Barcode))
data
data <- data %>% mutate(Cat_sample = ifelse(duplicated(ID_GQ), "Duplicate", Cat_sample),  # only if some samples are duplicated multiple times
                        ID_GQ = ifelse(duplicated(ID_GQ), paste0(ID_GQ, "_1"), ID_GQ))

data %>% View()
colnames(data)[8] <- "Espece"

## Format Sampling region (Pop)
data %>% filter(Cat_sample %in% "Sample") %>%
         group_by(Region_echantillonnage) %>%
         summarise(N = n()) %>% View()

data <- data %>% arrange(No_plaque_envoi, No_puits_envoi)

write_csv(data, file = file.path(get.value("info.path"),"Project_Infos.csv"))


# Create the first popmap  

pop.info <- data %>% select(ID_GQ, Region_echantillonnage) %>% arrange(ID_GQ)
pop.info

write.table(pop.info, 
            file = file.path(get.value("info.path"), "popmap.txt"),
            quote = FALSE, sep = "\t",
            row.names = F, col.names = F)


# Create barcode files
## Check that plate ID are the same as the files WHAT IF THEY ARE NOT?
data %>% pull(No_plaque_envoi) %>% unique()


## Check that every individual as an unique ID
data %>% pull(ID_GQ) %>% length()
data %>% pull(ID_GQ) %>% unique() %>% length()
data %>% filter(duplicated(Numero_unique_specimen))  # pull(ID_GQ[duplicated(ID_GQ)])# %>% duplicated()

# Loop to create all the files relating individuals to barcodes
for(x in unique(data$No_plaque_envoi)){
  
  write.table(data %>% filter(No_plaque_envoi == x) %>%
                select(Barcode, ID_GQ),
              file = file.path(get.value("info.path"), paste(x,"barcodes.txt", sep = "_")) %>%
                str_remove("[A-Z]*[_][:digit:][:digit:][_][:digit:][:digit:][:digit:][:digit:]"),  # included this because names of raw sequence do not have the form
                                                                                                   # PE_20_00001, so I'm just keeping the last digit (otherwise issues
                                                                                                   # when demultiplexing as the when looking for barcode 2 it will pick
                                                                                                   # up all barcodes because of the _20_)
              quote = FALSE, sep = "\t",
              row.names = F, col.names = F)
  
}

