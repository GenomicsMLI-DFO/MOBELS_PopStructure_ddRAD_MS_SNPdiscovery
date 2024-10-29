# Info --------------------------------------------------------------------
# 
# Authors: Luca Montana
# Affiliation: Fisheries and Oceans Canada (DFO)
# Group: Genomic laboratory, Demersal and Benthic Sciences Branch 
# Location: Maurice Lamontagne Institute
# Date: 2022-04-01
# 
# Overview: Create directories and move scripts in appropriate '01_Codes directory'
# 
#


# Housekeeping ------------------------------------------------------------

# Verify if you're in the right directory
getwd()  # the home directory should be where the Rproj is (/media/genyoda/Extra_Storage/Projets/Plies/Flatfish_ddRAD_Novaseq)

# Clear workspace
rm(list = ls())

# Libraries

# Functions



# Create directories ------------------------------------------------------

# 00_Data
if(dir.exists("./00_Data") == F){
  dir.create("./00_Data")
}

## Functions 00_FileInfos in 00_Data
if(dir.exists("./00_Data/00_FileInfos") == F){
  dir.create("./00_Data/00_FileInfos")
}

## Functions 01_Raw in 00_Data
if(dir.exists("./00_Data/01_Raw") == F){
  dir.create("./00_Data/01_Raw")
}

# 01_Codes
if(dir.exists("./01_Codes") == F){
  dir.create("./01_Codes")
}

## Functions folder in 01_Codes
if(dir.exists("./01_Codes/Functions") == F){
  dir.create("./01_Codes/Functions")
}

# 02_Results
if(dir.exists("./02_Results") == F){
  dir.create("./02_Results")
}

## 00_Stacks folder in 02_Results
if(dir.exists("./02_Results/00_Stacks") == F){
  dir.create("./02_Results/00_Stacks")
}



# Move scripts in 01_Codes ------------------------------------------------

# Select files to move
my_files <- list.files(pattern = "\\.R$")
my_files

# Copy files from home directory to 01_Codes
file.copy(from = paste0("./", my_files),
          to = paste0("./01_Codes/", my_files))

# Remove files from home directory
file.remove(from = paste0("./", my_files))



# Unzip and move fasta sequences in 00_Data/01_Raw (Terminal) -------------

# cd /media/genyoda/Extra_Storage/Projets/MOBELS/
# unzip P1_P9_NS.zip -d ./Beluga_ddRAD_Novaseq_2019/00_Data/01_Raw/

# cd ../../media/genyoda/Extra_Storage/Projets/MOBELS/Beluga_ddRAD_Novaseq_2019/00_Data/01_Raw/P1_P9_NS/
# mv *.fastq.gz ../
# mv *.fastq.gz.md5