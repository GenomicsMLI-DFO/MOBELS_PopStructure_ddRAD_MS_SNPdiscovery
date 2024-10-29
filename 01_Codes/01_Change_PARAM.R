# Info --------------------------------------------------------------------
# 
# Authors: Audrey Bourret, Luca Montana
# Affiliation: Fisheries and Oceans Canada (DFO)
# Group: Genomic laboratory, Demersal and Benthic Sciences Branch 
# Location: Maurice Lamontagne Institute
# Date: 2021-01-12
# 
# Overview: This file allowed to set/modify folder structure
# 
#


# Library -----------------------------------------------------------------

# Internal functions
if(length(list.files("./01_Codes/Functions"))>=1){ 
  for(i in 1:length( list.files("./01_Codes/Functions") )){
    source(file.path("./01_Codes/Functions",  list.files("./01_Codes/Functions")[i]))  
  }
}

PARAM <- data.frame(x = character(), value = character())


# Specific functions ------------------------------------------------------

"add.param<-" <- function(PARAM, ..., value){
  PARAM <- rbind(PARAM, value)
  names(PARAM) <- c("x", "value")
  PARAM$x     <- as.character(PARAM$x)
  PARAM$value <- as.character(PARAM$value)
  PARAM
  
}


# Path --------------------------------------------------------------------

# 00_Data
add.param(PARAM) <- c("info.path", "./00_Data/00_FileInfos")
add.param(PARAM) <- c("access.path", "../ACCESS")
add.param(PARAM) <- c("raw.path", "./00_Data/01_Raw") 
#add.param(PARAM) <- c("cutadapt.path", "./00_Data/02_Cutadapt") 
add.param(PARAM) <- c("trimmo.path", "./00_Data/02_Trimmomatic") 

add.param(PARAM) <- c("demulti.path", "./00_Data/03a_Demultiplex") 
add.param(PARAM) <- c("align.path", "./00_Data/03b_Align") # External disk
#add.param(PARAM) <- c("align.mtDNA.path", "./00_Data/03c_Align_mtDNA") # External disk

#add.param(PARAM) <- c("samples.path", "./00_Data/04_Samples") # External disk

add.param(PARAM) <- c("test.denovo.path", "./00_Data/04a_Test.denovo") # External disk
add.param(PARAM) <- c("test.ref.path", "./00_Data/04b_Test.ref") # External disk
add.param(PARAM) <- c("stacks.ref.path", "./00_Data/05b_Stacks.ref") # External disk
add.param(PARAM) <- c("filter.ref.path", "./00_Data/06b_Filtering.ref")

add.param(PARAM) <- c("ref.genome.path", "./00_Data/99_REF_Genome")

# Results - Stacks

add.param(PARAM) <- c("fastqc.path", "./02_Results/00_Stacks/01_FastQC")

add.param(PARAM) <- c("fastqc.raw.path", "./02_Results/00_Stacks/01_FastQC/01_Raw")
add.param(PARAM) <- c("fastqc.trim.path", "./02_Results/00_Stacks/01_FastQC/02_Trimmomatic")
add.param(PARAM) <- c("fastqc.demulti.path", "./02_Results/00_Stacks/01_FastQC/03_Demultiplex")


add.param(PARAM) <- c("cutadapt.log", "./02_Results/00_Stacks/02_Cutadapt")
add.param(PARAM) <- c("trimmo.log", "./02_Results/00_Stacks/02_Trimmomatic")
add.param(PARAM) <- c("demulti.log", "./02_Results/00_Stacks/03_Demultiplex")
add.param(PARAM) <- c("stacks.log", "./02_Results/00_Stacks/05a_Stacks.denovo")
add.param(PARAM) <- c("stacks.ref.log", "./02_Results/00_Stacks/05b_Stacks.ref")

add.param(PARAM) <- c("stacks.ref.path.bringloe", "./02_Results/00_Stacks/05b_Stacks.ref")
stacks.ref.path.bringloe


# Update file before continue to ensure next section will do OK!
if(dir.exists("./00_Data/00_FileInfos") == F){
  dir.create("./00_Data/00_FileInfos")
}

write.csv2(PARAM, file = file.path("./00_Data/00_FileInfos", "Options.csv"), row.names=F)


# Files -------------------------------------------------------------------

# Adapter file

add.param(PARAM) <- c("adapters", file.path(get.value("info.path"),"adapters.fasta"))


# External programs -------------------------------------------------------


# Data --------------------------------------------------------------------


# Save Parameters ---------------------------------------------------------

PARAM

write.csv2(PARAM, file = file.path("./00_Data/00_FileInfos", "Options.csv"), row.names=F)

get.value("raw.path")



