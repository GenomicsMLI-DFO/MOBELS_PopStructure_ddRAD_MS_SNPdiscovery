# Info --------------------------------------------------------------------
# 
# Authors: Audrey Bourret, Luca Montana
# Affiliation: Fisheries and Oceans Canada (DFO)
# Group: Genomic laboratory, Demersal and Benthic Sciences Branch 
# Location: Maurice Lamontagne Institute
# Date: 2022-03-10
# 
# Overview: RAD-seq pipeline with STACKS version >2, up-to Gstacks
# NS data
# Trimmomatics + no rad check on one side (cut 5 pb) MAYBE
# With and without ref genome
#



# Library -----------------------------------------------------------------

library(parallel)
library(tidyverse)
library(readxl)


# Internal functions
for(i in 1:length( list.files("./01_Codes/Functions") )){
  source(file.path("./01_Codes/Functions",  list.files("./01_Codes/Functions")[i]))  
}

"%nin%" <- Negate(`%in%`)

# Add python env to this specific project
Sys.setenv(PATH = paste(c("/home/genyoda/Documents/PythonVenv/GenoBaseEnv/bin",
                          Sys.getenv("PATH")),
                        collapse = .Platform$path.sep))

system2("multiqc", "--help")

stack2.55_path <- "/home/genyoda/Documents/Programs/stacks-2.55/"



# Check folders -----------------------------------------------------------
# To run if you need the folder skeleton 

#auto.folder()




# Data --------------------------------------------------------------------

pop.info <- read.table(file.path(get.value("info.path"), "popmap.txt"))
names(pop.info) <- c("Sample", "POP")
pop.info


pop.data <- read_csv(file.path(get.value("info.path"),"Project_Infos.csv"))
pop.data 
colnames(pop.data)[3] <- "No_plaque"

# Define initial working directory (just in case something goes wrong)
current.wd <- getwd()

numCores <- if(get_os() %in% c("os","linux")){
  detectCores() # Utilise le max de coeurs si sur linux
} else 1

cat("There's", numCores, "cores available in this computer", sep = " ")




# FastQC ------------------------------------------------------------------

# Took around 30 min by plate (NS and HI)
# Now with an overwrite function

fastqc <- function(folder.in, folder.out, overwrite = F, nthread) {
  #@folder.in : where to find the data
  #@folder.out : where to put the data
  #@overwrite : if T, will remove everything within the folder  
  
  if(get_os() %in% c("os","linux")){ # to run only on os and not windows
    
    cat("Performing a FastQC analysis\n")
    
    # Remove old files or not (addition June 2020)
    if(isTRUE(overwrite)){
      file.remove(list.files(folder.out, full.name = T, pattern ="fastqc"))
      files.to.use <- list.files(folder.in, full.names = T, pattern = "fastq") %>% 
        str_remove(".md5") %>% 
        str_subset("R1.fastq|R2.fastq") %>% unique()
    }
    
    if(isTRUE(overwrite == F)){
      previous.analysis <- list.files(folder.out, pattern = ".html") %>% str_replace("_fastqc.html", ".fastq.gz")
      
      files.to.use <- list.files(folder.in, full.names = T, pattern = "fastq") %>% 
        str_remove(".md5") %>% unique() %>% 
        str_subset(paste(previous.analysis, collapse = "|"), negate = T) %>% 
        str_subset("R1.fastq|R2.fastq")
      
    }
    
    cat("Results could be find here:", folder.out ,"\n")
    
    mclapply(files.to.use,
             FUN = function(x){
               
               cmd <- paste("--outdir", folder.out, x, 
                            "-t", 1)
               system2("fastqc", cmd) 
               
             } ,
             mc.cores = nthread
    )
    
    
  } else {cat("Cannot perform FastQC on windows yet -- sorry!!")}
} # End of my function


# Test a multiqc function

multiqc <- function(folder.in){
  # Multi QC aggregation - run pretty fast ...
  for(s in c("R1", "R2")){
    print(s)  
    cmd <- paste(paste(list.files(folder.in, full.names = T) %>% 
                         #str_subset(l) %>% 
                         str_subset(paste0("_",s)) %>% 
                         str_subset(".zip"), collapse = " "),
                 "--outdir", file.path(folder.in, "MultiQC_report"),
                 "--filename", paste0("multiqc_report_", s, ".html"),
                 "-f" # to rewrite on previous data
    )
    
    system2("multiqc", cmd)
    
  } 
  
}

if(dir.exists("./02_Results/00_Stacks/01_FastQC") == F){
  dir.create("./02_Results/00_Stacks/01_FastQC")
}

if(dir.exists("./02_Results/00_Stacks/01_FastQC/01_Raw") == F){
  dir.create("./02_Results/00_Stacks/01_FastQC/01_Raw")
}


# Run FastQC - change for overwrite T for the first time to or it will not work
fastqc(folder.in = get.value("raw.path"), folder.out = get.value("fastqc.raw.path"), overwrite = F, nthread = 40)


# Multi QC aggregation - runs pretty fast ...
multiqc(folder.in = get.value("fastqc.raw.path"))

# CHECK on the R2 side that the adapter "CGG" is present. 
# IF not, the PART2 of trimmomatic allows to remove 5 pb on this side




# Trimmomatics ------------------------------------------------------------

if(dir.exists("./00_Data/02_Trimmomatic") == F){  # trimmo.path
  dir.create("./00_Data/02_Trimmomatic")
}

if(dir.exists("./02_Results/00_Stacks/02_Trimmomatic") == F){  # trimmo.log
  dir.create("./02_Results/00_Stacks/02_Trimmomatic")
}

# In paired-end
# To remove the Illumina adapter and to cut 3pb in R2

# Check that you can reach the program
trimmomatic.path <- "/home/genyoda/Documents/Programs/Trimmomatic-0.39/trimmomatic-0.39.jar"
system2("java", paste("-jar", trimmomatic.path, "PE", "-version"), stdout=T, stderr=T) 

# Check which files you want to use

files.to.use <- list.files(get.value("raw.path"), full.names = T) %>% 
  str_subset("R1") %>% 
  str_subset(".md5", negate = T)

files.to.use

#start_time3 <- Sys.time()

### BENCHMARK ###
# On 4 files of 1000 reads
# 80 cores in Trimmomatics: 1.249837 secs
# 8 cores in Trimmomatics: 1.179123 secs 
# 4 cores in Trimmomatics / parallele in R 4 cores: 0.3734035 secs
#  1 core in Trimmomatics / parallele in R 4 cores: 0.3527803 secs

# With 60 cores, memory is full, so 20 seems a better compromise
# checked in Terminal with "free -m" or "top"

# STEP 1 - REMOVE 3 PB on R2 (must do it first because of misalignment)

mclapply(files.to.use,
         FUN = function(x){
           
           # Names of the important files
           file1 <- x
           file2 <- file1 %>% str_replace("_R1", "_R2")
           
           file2.out <- file2 %>% str_replace(get.value("raw.path"), get.value("trimmo.path")) %>% 
             str_replace("_R2.fastq.gz", "_R2_HC3.fastq.gz")
           
           logout <- file1 %>% str_replace(get.value("raw.path"), get.value("trimmo.log")) %>% 
             str_replace("_R1.fastq.gz", ".log")
           
           # The command
           cmd1 <- paste("-jar",
                         trimmomatic.path, "SE",
                         "-threads", 2,
                         #"-trimlog", logout %>% str_replace(".log", "_HC.log"),
                         file2, file2.out,
                         #"ILLUMINACLIP:00_Data/00_FileInfos/adapters/TruSeq3-PE-2.fa:2:30:10",
                         "HEADCROP:3",
                         sep = " ")
           
           A1 <- system2("java", cmd1, stdout=T, stderr=T) # from R console to Linux bash terminal
           A1
           
           # save a file log
           cat(file = logout,
               #"STEP1 - REMOVE ILLUMINA ADAPTORS", cmd1, A1,
               "STEP1 - CUT 3 pb ON R2 PAIRED because of quality drop on Novaseq",cmd1, A1, # what to put in my file
               append= F, sep = "\n\n")
           
         } ,
         mc.cores = 15
)


# STEP 2 - REMOVE Illumina adaptors

mclapply(files.to.use,
         FUN = function(x){
           
           # Names of the important files
           file1 <- x
           file2 <- file1 %>% str_replace(get.value("raw.path"), get.value("trimmo.path")) %>% 
             str_replace("_R1.fastq.gz", "_R2_HC3.fastq.gz")
           
           fileout <- file1 %>%  str_replace(get.value("raw.path"), get.value("trimmo.path")) %>% 
             str_replace("_R1.fastq.gz", ".fastq.gz")   
           
           logout <- file1 %>% str_replace(get.value("raw.path"), get.value("trimmo.log")) %>% 
             str_replace("_R1.fastq.gz", ".log")
           
           # The command
           
           cmd1 <- paste("-jar",
                         trimmomatic.path, "PE",
                         "-threads", 2,
                         #"-trimlog", logout,
                         file1, file2,
                         "-baseout", fileout,
                         "ILLUMINACLIP:00_Data/00_FileInfos/adapters/TruSeq3-PE-2.fa:2:30:10",
                         #"HEADCROP:1-99",
                         sep = " ")
           
           A1 <- system2("java", cmd1, stdout=T, stderr=T) # from R console to Linux bash terminal
           A1
           
           # save a file log
           cat(file = logout,
               "STEP2 - REMOVE ILLUMINA ADAPTORS", cmd1, A1,
               #"STEP1 - CUT 3 pb ON R2 PAIRED because of quality drop on Novaseq",cmd1, A1, # what to put in my file
               append= T, sep = "\n\n")
           
         } ,
         mc.cores = 15
)


# SOME POST ANALYSIS FILE MANIPULATION

# Rename the files

old.name <-list.files(get.value("trimmo.path"), full.names = T)
old.name

new.name <- old.name %>% 
  str_replace("_1P.fastq", "_R1.fastq") %>% 
  str_replace("_2P.fastq", "_R2.fastq") 

new.name

file.rename(from = old.name,
            to = new.name)

# Compute stats by plate/RunSeq

# Will work for paired reads

trim.data <- data.frame(ID = character(),
                        RunSeq = character(),
                        No_plaque = character(),
                        Ntotal = numeric(),
                        Nsurvival = numeric(),
                        stringsAsFactors = F)

for(x in  list.files(get.value("trimmo.log"), pattern = ".log")){
  
  ID <- x %>% str_remove(".log")
  ID.int = ID %>% str_replace("[.][A-Z][:digit:][:digit:][:digit:][.][A-Z]([a-z]*)[-]", "___")
  RunSeq = sapply(str_split(ID.int, "___"), `[`, 1)
  No_plaque = sapply(str_split(ID.int, "___"), `[`, 2)
  temp <- readLines(file.path(get.value("trimmo.log"), x ))
  temp <- temp %>% str_subset(pattern = "Input Read Pairs")
  Ntotal    <- sapply(str_split(temp, " "), `[`, 4)
  Nsurvival <- sapply(str_split(temp, " "), `[`, 7)
  
  trim.data <- bind_rows(trim.data,
                         data.frame(ID = ID,
                                    RunSeq = RunSeq,
                                    No_plaque = No_plaque,
                                    Ntotal = as.numeric(Ntotal),
                                    Nsurvival = as.numeric(Nsurvival),
                                    stringsAsFactors = F))
  
}

trim.data <- trim.data %>% filter(Ntotal != Nsurvival)

write_csv(x = trim.data, file = file.path(get.value("trimmo.log"), "Summary_Nreads.csv"))

gg.trim <- trim.data %>% ggplot(aes(x = No_plaque, y = Nsurvival, fill = RunSeq)) +
  geom_bar(stat = "identity") +
  labs(title = "N reads by plates post trimmomatic") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) 

gg.trim

gg.trim <- trim.data %>% mutate(perc_surv = Nsurvival/Ntotal) %>% 
  ggplot(aes(x = perc_surv)) +
  geom_histogram() +
  labs(title = "N reads by plates post trimmomatic") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) 

gg.trim

gg.trim <- trim.data %>% mutate(perc_surv = Nsurvival/Ntotal) %>% 
  ggplot(aes(x = Nsurvival, fill = No_plaque)) +
  geom_histogram() +
  labs(title = "N reads by plates post trimmomatic") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) 

gg.trim

gg.trim <- trim.data %>% mutate(perc_surv = Nsurvival/Ntotal) %>% 
  ggplot(aes(x = Ntotal, y = Nsurvival, col = No_plaque)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  labs(title = "N reads by plates post trimmomatic") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) 

gg.trim

# ggsave(filename = file.path(get.value("trimmo.log"), "Summary_Nreads.png"),
#        plot = gg.trim,
#        width = 5, height = 4, units = "in")


# If you want the log file to be ignored, run the following :

cat("*.log", "!.gitignore", sep = "\n",
    file = file.path(get.value("trimmo.log"), ".gitignore"))


# Run FastQC - change for overwrite T first the first time to or it will not work
# THIS PART SHOULD BE BENCHMARKED TOO

if(dir.exists("./02_Results/00_Stacks/01_FastQC/02_Trimmomatic") == F){
  dir.create("./02_Results/00_Stacks/01_FastQC/02_Trimmomatic")
}

fastqc(folder.in = get.value("trimmo.path"), folder.out = get.value("fastqc.trim.path"), 
       overwrite = F, nthread = 20)

multiqc(folder.in = get.value("fastqc.trim.path"))




# Demultiplex -------------------------------------------------------------

# Barcode files are important for thise step
# They are created with the 01_Create_PopMap.R

# TAKE CARE, all individuals with similar names will be collapsed
length(pop.info$Sample)
length(pop.info$Sample %>% unique())

# If you specify the same output directory for two different
# process_radtags runs, the second run will overwrite identical filenames
# from the first run. Output the data into separate directories, and then
# concatenate the shared samples together after the runs complete

sub.dir <- c("NS.1272.002", "NS.1272.003")
sub.dir

# Before running, check that the right barcode file is found - here an example
# not in a loop

files.to.use <-  list.files(get.value("trimmo.path"), full.names = T, pattern = "_R1") %>% 
  str_subset(sub.dir[1])
files.to.use

barcode <- list.files(get.value("info.path"), full.names = T, pattern = "barcodes.txt") %>% 
  str_subset(files.to.use[1] %>% str_remove(get.value("trimmo.path")) %>% 
               str_remove("_R1.fastq.gz") %>% 
               str_remove(sub.dir[1]) %>%
               str_remove("[.][A-Z][:digit:][:digit:][:digit:][.][A-Z]([a-z]*)[-][A-Z]") %>% 
               str_remove("/")) 
barcode


if(dir.exists("./00_Data/03a_Demultiplex") == F){
  dir.create("./00_Data/03a_Demultiplex")
}

for(i in sub.dir){
  
  files.to.use <-  list.files(get.value("trimmo.path"), full.names = T, pattern = "_R1") %>% 
    str_subset(i) #%>% str_subset("P02|P03")
  
  # files.to.use
  
  cat(paste("\nWorking with the run:", i),
      paste("There are" , length(files.to.use) , "plates within this run"),
      sep= "\n")
  
  # Create a sub-directory if needed (i.e. it doesn't exists already)
  if(file.exists(file.path(get.value("demulti.path"), i)) == F) {
    
    cat(paste("\nCreating a new directory:", file.path(get.value("demulti.path"), i)),
        sep= "\n")
    
    dir.create(file.path(get.value("demulti.path"), i), recursive = T)}
  
  
  # Parallel version of process_radtag 
  
  stack2.55_path <- "/home/genyoda/Documents/Programs/stacks-2.55/" 
  
  mclapply(files.to.use,
           FUN = function(x){
             
             # Files
             file1 <- x
             file2 <- file1 %>% str_replace("_R1", "_R2")
             
             barcode <- list.files(get.value("info.path"), full.names = T, pattern = "barcodes.txt") %>% 
               str_subset(x %>% str_remove(get.value("trimmo.path")) %>% 
                            str_remove("_R1.fastq.gz") %>% 
                            str_remove(i) %>%
                            str_remove("[.][A-Z][:digit:][:digit:][:digit:][.][A-Z]([a-z]*)[-][A-Z]") %>%  # adjust string
                            str_remove("/")) 
             
             plate <- barcode %>% str_remove(get.value("info.path")) %>% 
               str_remove("_barcodes.txt") %>% 
               str_remove("/")
             
             # Create a sub-directory if it doesn't exists already
             if(file.exists(file.path(get.value("demulti.path"), i, plate)) == F) {
               
               #cat(paste("\nCreating a new directory:", file.path(get.value("demulti.path"), i)),
               #     sep= "\n")
               
               dir.create(file.path(get.value("demulti.path"), i, plate), recursive = T)}
             
             # Command
             cmd <- paste("--paired",
                          "-1", file1,
                          "-2", file2,
                          "-o", file.path(get.value("demulti.path"), i, plate),
                          "--inline_null",   
                          "-b", barcode,
                          "--renz_1", "pstI", # Check RAD site on R1 
                          #"--renz_2", "mspI", # CGG on R2 
                          "-E", "phred33",
                          "--filter-illumina",
                          "-c", # clean
                          "-r", # rescue
                          "-q", # check quality
                          "-t", 135,  # truncate at 150 - 6 (restriction site) - 8 (max barcode) = 136 (all reads within sample must be the same length)
                          "-i", "gzfastq"
             )
             
             A <- system2(paste0(stack2.55_path, "process_radtags"), cmd, stdout = T, stderr = T)
             A
             # save a log file 
             log.file <- file1 %>% str_replace(get.value("trimmo.path"), get.value("demulti.log")) %>% 
               str_replace(".fastq.gz","_summary.log") %>% 
               str_remove("_R1")
             
             cat(file = log.file,
                 cmd, "\n\n",
                 A, # what to put in my file
                 append= F, sep = "\n")
             
             # Detailed log file
             file.rename(from = list.files(file.path(get.value("demulti.path"), i, plate), full.names = T, pattern = ".log"),
                         to = log.file %>% str_replace("summary", "detailed")
             )
             
           } ,
           mc.cores = 20
  )
  
  gc()
  
}

# If you want the log/csv file to be ignored, run the following :

cat("*", "!.gitignore", "!*.png", "!AllIndividuals_Nreads.csv", sep = "\n",
    file = file.path(get.value("demulti.log"), ".gitignore"))



# Extract data for each summary_detailed
# You must change the str_subset code that is the index for each project

for(x in list.files(get.value("demulti.log"), pattern = "detailed", full.names = T)){
  
  data <- readLines(x) %>% 
    str_subset("S_") %>% 
    str_split(pattern = "\t")
  
  cat("\n",length(data), "samples retrieved in", x %>% str_remove(get.value("demulti.log")) %>% str_remove("/"))
  
  data <-  data.frame(matrix(unlist(data), nrow= length(data), byrow = T))
  # If there is a RUN column, it's because some individuals come from more than one sample
  names(data) <- c("Barcode", "Filename", "Total", "NoRadTag", "LowQuality", "Retained")
  
  data <- data %>% mutate(Total = as.numeric(as.character(Total)),
                          NoRadTag = as.numeric(as.character(NoRadTag)),
                          LowQuality = as.numeric(as.character(LowQuality)),
                          Retained = as.numeric(as.character(Retained)),
                          Run = x %>% str_remove(get.value("demulti.log")) %>% 
                            str_remove("_log_detailed.txt")  %>% 
                            str_remove("/")
  )
  write_csv(data, x %>% str_replace("detailed.log", "Nreads.csv"))
}

# Create one big log file 

Nreads.data <- data.frame()

for(x in list.files(get.value("demulti.log"), pattern = "Nreads", full.names = T) %>% str_subset(".csv") %>% str_subset("NS.")){
  data.int <- read_csv(x)
  
  Nreads.data <- bind_rows(Nreads.data, data.int)
  
}

# Compute N read removed
Nreads.data$Removed <- Nreads.data$Total - Nreads.data$Retained  

# Add a column for the run name
Nreads.data$RunSeq <- paste(sapply(str_split(Nreads.data$Run, "[.]"),`[`,1),
                            sapply(str_split(Nreads.data$Run, "[.]"),`[`,2),
                            sapply(str_split(Nreads.data$Run, "[.]"),`[`,3),
                            sep = ".")


# Save the result
write_csv(Nreads.data, file.path(get.value("demulti.log"),"AllIndividuals_Nreads.csv"))

Nreads.data <- read_csv(file.path(get.value("demulti.log"),"AllIndividuals_Nreads.csv"))

trim.data <- read_csv(file = file.path(get.value("trimmo.log"), "Summary_Nreads.csv"))
trim.data$No_plaque <- paste0(str_sub(trim.data$No_plaque, 1, 1), "E_20_0000", str_sub(trim.data$No_plaque, 2, 2))  # update name plate to match that in Nreads.data

Nreads.data <- Nreads.data %>% left_join(pop.data, by = c("Filename" = "ID_GQ")) 

# Check that we have all the data - all with unique data

head(Nreads.data)
Nreads.data %>% nrow() / length(sub.dir)
Nreads.data$Filename %>% unique() %>% length()

# Stats by plate
Nreads.data %>% group_by(No_plaque) %>% summarise(Retained = sum(Retained)) %>%
  group_by() %>% 
  summarise(mean = mean(Retained)/2,  # divided by 2 because each sequence has R1 and R2
            sd = sd(Retained)/2,
            min = min(Retained)/2,
            max = max(Retained)/2)

# Verify content (reads retained) of each well

gg.retained.well <- Nreads.data %>%
  group_by(No_puits_envoi, No_plaque) %>%
  mutate(Column = as.integer(str_remove(No_puits_envoi, "[A-Z]")),
         Row = str_remove(No_puits_envoi, "[:digit:][:digit:]")) %>%
  ggplot(aes(x = Row, y = factor(Column), fill = log(Retained))) +
  geom_bin2d() +
  stat_bin2d(geom = "text", aes(label = Filename), binwidth = 1, color = "white") +
  facet_grid(RunSeq ~ No_plaque) + theme_bw()

gg.retained.well

# Impact of overall process

gg.radtag <- Nreads.data %>% group_by(No_plaque, RunSeq) %>% 
  #mutate(RunSeq = paste()) 
  summarise(Total = sum(Total)/2,
            Retained = sum(Retained)/2) %>% 
  left_join(trim.data %>% select(RunSeq, No_plaque, Nsurvival)) %>% 
  mutate(perc_with_barcode = Total/Nsurvival,
         perc_keep = Retained / Nsurvival,
         perc_keep_over_barcode = Retained / Total) %>% 
  select(RunSeq, No_plaque, perc_with_barcode, perc_keep,  perc_keep_over_barcode) %>% 
  pivot_longer(c(perc_with_barcode, perc_keep,  perc_keep_over_barcode), names_to = "Cat", values_to = "Perc") %>% 
  ggplot(aes(x = No_plaque, y = Perc, col = RunSeq, group = RunSeq)) +
  geom_point(position= position_dodge(width = 1/4)) +
  labs(title = "% reads retained post demultiplexing") +
  facet_grid(Cat ~ . ) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) 

gg.radtag

# ggsave(filename = file.path(get.value("demulti.log"), "Summary_Preads_byPlates.png"),
#       plot = gg.radtag,
#       width = 8, height = 6, units = "in")


gg.retained.plate <- Nreads.data %>% #gather(Removed, Retained, key = "Cat", value = "N") %>% #head()
  group_by(RunSeq, No_plaque) %>% summarise(Retained = sum(Retained)/2,
                                            Removed = sum(Removed)/2) %>% 
  ggplot(aes(y = Retained, x = No_plaque, fill = RunSeq)) + 
  geom_bar(stat = "identity")+
  theme_bw() + 
  labs(title = "N reads retained post demultiplexing") +
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "bottom")

gg.retained.plate

# ggsave(filename = file.path(get.value("demulti.log"), "Summary_Nreads_byPlates.png"),
#       plot = gg.retained.plate,
#       width = 5, height = 4, units = "in")


# Check that we have all the data - all with unique data

head(Nreads.data)
Nreads.data %>% nrow() / length(sub.dir)
Nreads.data$Filename %>% unique() %>% length()

# Stats by IND

Nreads.data %>% group_by(Filename) %>% summarise(Retained = sum(Retained)) %>%
  group_by() %>% 
  summarise(mean = mean(Retained)/2,
            sd = sd(Retained)/2,
            min = min(Retained)/2,
            max = max(Retained)/2)


gg.retained.ind <- Nreads.data %>% gather(Removed, Retained, key = "Cat", value = "N") %>% #head()
  group_by(Filename, No_plaque, Cat) %>% summarise(N = sum(N)/2) %>% 
  ggplot(aes(x = Filename, y = N, fill= Cat)) + 
  geom_bar(stat= "identity") +
  #scale_y_continuous(breaks = c(1000000, 2000000, 3000000, 4000000, 5000000, 6000000, 7000000),
  #                   labels = c(1:7))+
  labs(y = "N reads", x = "Samples", title = "N reads by ind") +
  facet_wrap(~ No_plaque, scale = "free_x") +
  theme_minimal() + 
  theme(axis.text.x = element_blank(),
        legend.position = "bottom")

gg.retained.ind

# ggsave(filename = file.path(get.value("demulti.log"), "Summary_Nreads_byInd.png"),
#       plot = gg.retained.ind,
#       width = 8, height = 8, units = "in")


gg.retained.ind2.a <- Nreads.data %>% 
  mutate(barcode_length = nchar(Barcode.y),
         first.nuc = str_sub(Barcode.y, 1,1)) %>% 
  group_by(Filename, barcode_length, first.nuc) %>% summarise(Retained = sum(Retained)/2) %>%
  #  arrange(Retained) %>% 
  ggplot(aes(x = reorder(Filename, Retained), y = Retained)) +
  scale_y_continuous(trans = "log10") +
  geom_bar(stat = "identity") +
  #facet_grid(barcode_length ~ first.nuc, scale = "free_x", space = "free") +
  theme_bw() +
  theme(axis.text.x = element_blank())

gg.retained.ind2.a


gg.retained.ind2.b <- Nreads.data %>% 
  mutate(barcode_length = nchar(Barcode.y),
         first.nuc = str_sub(Barcode.y, 1,1)) %>% 
  group_by(Filename, barcode_length, first.nuc) %>% summarise(Retained = sum(Retained)/2) %>%
  #  arrange(Retained) %>% 
  ggplot(aes(x = reorder(Filename, Retained), y = Retained)) +
  scale_y_continuous(trans = "log10") +
  geom_bar(stat = "identity") +
  facet_grid(barcode_length ~ first.nuc, scale = "free_x", space = "free") +
  theme_bw() +
  theme(axis.text.x = element_blank())

gg.retained.ind2.b

gg.retained.ind2 <- ggpubr::ggarrange(gg.retained.ind2.a, gg.retained.ind2.b,
                                      nrow = 2, heights = c(1:3))
gg.retained.ind2

# ggsave(filename = file.path(get.value("demulti.log"), "Summary_Nreads_byInd_withBarcodes.png"),
#       plot = gg.retained.ind2,
#       width = 8, height = 8, units = "in")

gg.retained.pop <- Nreads.data %>% 
  group_by(Filename, Region_echantillonnage, No_plaque) %>% summarise(Retained = sum(Retained)/2) %>% 
  ggplot(aes(x = Region_echantillonnage, y = Retained)) + 
  geom_hline(yintercept = 0) +
  geom_boxplot() +
  geom_jitter(aes(col= factor(No_plaque)), cex = 1, alpha = 0.5) +
  labs(y = "N reads") +
  scale_y_continuous(trans = "log10") +
  facet_grid(.~Region_echantillonnage, space = "free", scales = "free") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "none")

gg.retained.pop

# ggsave(filename = file.path(get.value("demulti.log"), "Summary_Nreads_byPop.png"),
#       plot = gg.retained.pop,
#       width =7, height = 5, units = "in")

gg.barcodes <- Nreads.data %>% 
  group_by(No_plaque, Barcode.y, Filename,RunSeq) %>% summarise(Total = sum(Total)/2) %>% 
  ggplot(aes(x = Barcode.y, y = Total)) + 
  geom_boxplot() +
  geom_jitter(aes(col =RunSeq), cex = 1, alpha = 0.5) +
  # facet_wrap(~Barcode.y)
  # scale_y_continuous(breaks = c(1000000, 2000000, 3000000, 4000000, 5000000, 6000000, 7000000, 8000000, 9000000, 10000000),
  #                   labels = c(1:10))+
  #  scale_y_continuous(trans = "log10") +
  labs(y = "N total reads (log)", x = "Barcodes") +
  # facet_grid(. ~ Espece, scale = "free_x", space = "free") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1 ),
        legend.position = "bottom")

gg.barcodes

# ggsave(filename = file.path(get.value("demulti.log"), "Summary_Nreads_byBarcode.png"),
#       plot = gg.barcodes,
#       width =8, height = 5, units = "in")


# This is THE graph
gg.barcodes.summary <- Nreads.data %>% 
  group_by(Barcode.y, Filename) %>% summarise(Total =mean(Total)/2, 
                                              Retained =mean(Retained)/2,
                                              NoRadTag = mean(NoRadTag)/2,
                                              LowQuality = mean(LowQuality)/2) %>%
  pivot_longer(c(Total, Retained, NoRadTag, LowQuality), names_to = "Cat", values_to = "N") %>% 
  mutate(barcode_length = nchar(Barcode.y),
         first.nuc = str_sub(Barcode.y, 1,1)) %>% 
  ggplot(aes(x = barcode_length, y =N, group = barcode_length) )+ 
  geom_violin() +
  geom_jitter(alpha = 0.5, aes(col = first.nuc), cex = 1) +
  # geom_violin()+
  # facet_wrap(~Barcode.y)
  scale_y_continuous(trans = "log10") +
  labs(y = "N reads (log)", x = "Barcodes length (4-8)") +
  facet_grid(Cat ~ ., scale = "free_y") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1 ),
        legend.position = "bottom")

gg.barcodes.summary

# ggsave(filename = file.path(get.value("demulti.log"), "Summary_Nreads_byBarcode_Summary.png"),
#       plot = gg.barcodes.summary,
#       width =6, height = 6, units = "in")

gg.barcodes.summary2 <- Nreads.data %>% 
  group_by(Barcode.y, Filename) %>% summarise(Total =mean(Total)/2, 
                                              Retained =mean(Retained)/2,
                                              NoRadTag = mean(NoRadTag)/2,
                                              LowQuality = mean(LowQuality)/2) %>%
  pivot_longer(c(Total, Retained, NoRadTag, LowQuality), names_to = "Cat", values_to = "N") %>% 
  mutate(barcode_length = nchar(Barcode.y),
         first.nuc = str_sub(Barcode.y, 1,1)) %>% 
  group_by(first.nuc, barcode_length) %>% summarise(Nretained = mean(N[Cat == "Retained"])) %>% 
  arrange(Nretained) %>% 
  ggplot(aes(x = factor(barcode_length), y = first.nuc, fill = Nretained)) +
  geom_bin2d() + theme_bw()

gg.barcodes.summary2

# ggsave(filename = file.path(get.value("demulti.log"), "Summary_Nreads_byBarcode_Summary2.png"),
#       plot = gg.barcodes.summary2,
#       width =5, height = 4, units = "in")


## Create one file with up to 3 files by individual
## UPDATE: should work with up to 3 subdir
## IF only one run, just copy-paste in 'main' directory - code added 220506

sub.dir <- c("NS.1272.002", "NS.1272.003")
sub.dir
plates <- file.path(get.value("demulti.path"), sub.dir[1]) %>% list.files()
plates

for(p in plates){
  
  files.to.use <- list.files(file.path(get.value("demulti.path"), sub.dir[1], p), full.names = T)
  
  cat(paste("\nWorking with the plate:", p),
      paste("There are" , length(files.to.use) , "files to process"),
      sep= "\n")
  
  mclapply(files.to.use,
           FUN = function(x){
             
             # Files
             file1 <-  x
             
             if(length(sub.dir) >= 2){
               file2 <- file1 %>% str_replace(sub.dir[1], sub.dir[2])  
             }
             
             if(length(sub.dir) >= 3){
               file3 <- file1 %>% str_replace(sub.dir[1], sub.dir[3])  
             }
             
             file.join <- file1 %>% str_remove(file.path(sub.dir[1], p))
             
             if(length(sub.dir) == 1){
               fq.gz <- list.files(file.path(get.value("demulti.path"), sub.dir[1], p), pattern = ".fq.gz")
               path.from <- file.path(get.value("demulti.path"), sub.dir[1], p)
               path.to <- file.path(get.value("demulti.path"))
               file.copy(from = file.path(path.from, fq.gz),
                         to = file.path(path.to, fq.gz))
             }
             
             if(length(sub.dir) == 2){
               # Cat command - this is so simple !
               cmd <- paste(file1, file2, ">", file.join)
               system2("cat", cmd, stdout=T, stderr=T)
             }
             
             if(length(sub.dir) == 3){
               # Cat command - this is so simple !
               cmd <- paste(file1, file2, file3, ">", file.join)
               system2("cat", cmd, stdout=T, stderr=T)
             }
             
             
           } ,
           mc.cores = 20
  )
  
  gc()
  
}

# Number of files observed
list.files(file.path(get.value("demulti.path")), pattern = ".fq.gz") %>% length()

# Number of files expected
384 * length(plates)  # normally we want full plates as otherwise the number of reads results biased - it doesn't seem the case here
# the number obtained in the row above is 811*4 (4 files per sequence)

# FastQC - fastqc and multiqc functions work with multiple files at once, but don't with one sample at a time
# 
# fastqc(folder.in = get.value("demulti.path"), folder.out = get.value("fastqc.demulti.path"), 
#        overwrite = T, nthread = 20)
# 
# multiqc(get.value("fastqc.demulti.path"))




# Select individuals for tests --------------------------------------------

# Select between 10-20 representative individuals
# Will be used both for testing alignment and stacks

Nreads.data <- read_csv(file.path(get.value("demulti.log"), "AllIndividuals_Nreads.csv"))
Nreads.data <- Nreads.data %>% left_join(pop.data, by = c("Filename" = "ID_GQ"))

# Create the list of individuals for testing parameters :

hist(Nreads.data$Retained)
#summary(Nreads.data$Retained.by.sample)
summary(Nreads.data)

test.ID <- Nreads.data %>% filter(Espece %in% "Delphinapterus") %>%  # filter for belugas only as test.ID - also possibly remove Narwhals from alignment later on
  group_by(Filename) %>% 
  summarise(Retained = sum(Retained)) %>% 
  filter(Retained >= quantile(Retained, probs = 0.25),
         Retained <= quantile(Retained, probs = 0.75)) %>%
  sample_n(10) %>% pull(Filename)  

test.ID  # although each individual is duplicated in AllIndividuals_Nreads.csv it has one reads file list.files(get.value("demulti.path"), pattern = ".1.fq.gz", full.names = T)

write.table(pop.info %>% select(Sample, POP) %>% 
              filter(Sample %in% test.ID) %>% 
              mutate(POP = "noPOP"), 
            file = file.path(get.value("info.path"), "popmap.test_samples.txt"),
            quote = FALSE, sep = "\t",
            row.names = F, col.names = F)




# BWA ---------------------------------------------------------------------

# Indexing reference genome from Bringloe and Parent 2023. DOI: 10.1186/s12864-023-09779-3
# S_20_00703: https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_029941455.3/

## Indexing ---------------------------------------------------------------

list.files(get.value("ref.genome.path"))

cmd2 <- paste("index",
              "-p",  file.path(get.value("ref.genome.path"),"Beluga.align.TB.final.23Jan20"),
              "-a", "bwtsw",
              file.path(get.value("ref.genome.path"),"S_20_00703_MOBELS_0202.22xii22.FINAL_scaffolds-001.fasta")
)


A2 <- system2("bwa", cmd2, stdout=T, stderr=T)
A2

# save a log file
cat(file = file.path(get.value("ref.genome.path"), "Beluga.align.TB.final.23Jan20" ),
    cmd2, "\n\n",
    A2, # what to put in my file
    append= T, sep = "\n")



## Alignment --------------------------------------------------------------

# This part takes time, can start on a subset (and test Stack at the same times)
# 2-3 min by file

TEST.ID <- read.table(file.path(get.value("info.path"), "popmap.test_samples.txt")) %>% pull(V1) %>% paste(collapse = "|")

# ATTENTION - HERE FOR PAIRED 
# Test version
# demulti.files <- list.files(get.value("demulti.path"), pattern = ".1.fq.gz", full.names = T) %>%
#   str_subset(".rem.", negate = T) %>%
#   str_subset(TEST.ID)

# Complete version
demulti.files <- list.files(get.value("demulti.path"), pattern = ".1.fq.gz", full.names = T) %>%
  str_subset(".rem.", negate = T)

demulti.files %>% length()

#demulti.files[1:10]

mclapply(demulti.files,
         FUN = function(x){
           # How I will rename all this : 
           file.R1 <- x
           file.R2 <- x %>% str_replace(".1.fq.gz", ".2.fq.gz")
           file.bam <- x %>% str_replace(get.value("demulti.path"), file.path(get.value("align.path"), "01_Bringloe_22xii22")) %>% 
             str_replace(".1.fq.gz", ".bam")
           file.sort.bam <- file.bam %>% str_replace(".bam", ".sorted.bam")
           stat.tsv <- file.bam %>% str_replace(".bam", ".stat.tsv")
           
           # DO THE ALIGMENT  
           cmd1 <- paste("mem",
                         "-t", 32,
                         "-M",
                         file.path(get.value("ref.genome.path"),"Beluga.align.TB.final.23Jan20"), # the index ref genome
                         file.R1,
                         file.R2,
                         "2> /dev/null",
                         "| samtools", "view", "-Sb", 
                         "--threads", 15, # Number of additional threads
                         "-o", file.bam#,
           )# the file
           
           A <- system2("bwa", cmd1, stdout=T, stderr=T)
           
           # Sort the file 
           cmd2 <- paste("sort",
                         "--threads", 15,
                         #"-n",
                         "-O", "BAM",
                         # the new file (sorted)
                         "-o", file.sort.bam,
                         # the bam file to sort
                         file.bam
                         
           )
           
           A <- system2("samtools", cmd2, stdout=T, stderr=T)
           
           # Compute stats
           
           cmd3 <- paste("flagstat",
                         "--threads", 15,
                         "-O", "tsv",
                         file.sort.bam,
                         ">",
                         stat.tsv
           )
           
           A <- system2("samtools", cmd3, stdout=T, stderr=T)
           
         },
         mc.cores = 2
) 

cat("*.bam", "*.tsv", "!.gitignore", sep = "\n",
    file = file.path(get.value("align.path"), ".gitignore"))




# Compute the aligned reads -----------------------------------------------

map.res <- data.frame(ID = character(),
                      total = numeric(),
                      secondary = numeric(),
                      supplementary = numeric(),
                      duplicates = numeric(),
                      mapped = numeric(),
                      mapped_perc = numeric(),
                      paired = numeric(),
                      read1 = numeric(),
                      read2 = numeric(),
                      properly_paired = numeric(),
                      properly_paired_perc = numeric(),
                      twith_itself_mate_mapped = numeric(),
                      singletons = numeric(),
                      singletons_perc = numeric(),
                      twith_mate_mapped_diff_chr = numeric(),
                      twith_mate_mapped_diff_chr_HMQ = numeric(),
                      #Nmappedprim = numeric(),
                      stringsAsFactors = F)

files.to.use <- list.files(file.path(get.value("align.path"),"01_Bringloe_22xii22"), pattern = ".stat.tsv", full.names = T)  # specify file path to requested alignment

for(x in seq_along(files.to.use)){
  ID <-   files.to.use[x] %>% str_remove(file.path(get.value("align.path"),"01_Bringloe_22xii22")) %>%  # specify file path to requested alignment
    str_remove(".stat.tsv") %>% 
    str_remove("/")
  temp <-  sapply(str_split(readLines(files.to.use[x]), "\t"), `[`,1)  
  map.res[x,] <- c(ID, temp %>% str_remove("%"))
}  

for(x in 2:ncol(map.res)){
  map.res[,x] <- as.numeric(as.character(map.res[,x]))
}  

nrow(map.res)
head(map.res)  
summary(map.res) 

map.res <- map.res %>%  left_join(pop.data, by = c("ID" = "ID_GQ"))

graph1.0 <- map.res %>% #filter(ID != "Dp_2066") %>% 
  ggplot(aes(x = mapped_perc)) + 
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = 96, lty = "dashed") +
  #facet_grid(Espece ~ . , scales = "free")
  labs(x = "Percentage of reads mapped", y = "Density") +
  theme_bw()  
graph1.0

ggsave(filename = file.path(get.value("demulti.log"), "Alignment_density_TB_final_230123.png"),
       plot = graph1.0,
       width = 5, height = 4, units = "in")


graph1.1 <-  map.res %>%  ggplot(aes(y = mapped_perc, x = Region_echantillonnage)) + 
  geom_boxplot() +
  geom_jitter(aes(col = factor(Annee_echantillonnage)), height = 0, alpha = 0.5) +
  #facet_grid(Espece ~ . , scales = "free")
  labs(y = "Percentage of reads mapped") +
  geom_hline(yintercept = 96, lty = "dashed") +
  theme_bw() + 
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) 
graph1.1 

ggsave(filename = file.path(get.value("demulti.log"), "Alignment_boxplot_TB_final_230123.png"),
       plot = graph1.1,
       width = 5, height = 4, units = "in")


graph1.2 <- map.res %>% #filter(Espece_revision %in% c("Pb", "Pm")) %>% 
  ggplot(aes(x = total, y = mapped_perc, col = Region_echantillonnage)) +
  geom_point(alpha = 0.5)+
  geom_hline(yintercept = 96, lty = "dashed") +
  #scale_x_continuous(trans = "log10") + 
  labs(y = "Percentage of reads mapped", x = "N read total") +
  theme_bw() 
graph1.2


graph1.3 <- map.res %>% #filter(Espece_revision %in% c("Pb"))  %>% 
  ggplot(aes(y = mapped_perc, x = Region_echantillonnage)) + 
  geom_boxplot(alpha = 0.5) +
  geom_hline(yintercept = 96, lty = "dashed") +
  geom_jitter(aes(col = total),height = 0, alpha = 0.75) +
  scale_color_distiller(palette = "Spectral", trans = "log10") +
  #facet_wrap(~ Gen_ZONE)
  labs(y = "Percentage of reads mapped") +
  theme_bw()  +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
graph1.3 


graph1.4 <- map.res %>% #filter(Espece_revision %in% c("Pb"))  %>% 
  ggplot(aes(x = mapped, col = No_plaque)) +
  geom_density() +
  #geom_vline(xintercept = map.res %>% filter(Espece_revision %in% c("Pb")) %>% pull(mapped) %>% median, 
  #           lty = "dashed" )+
  #facet_wrap(~ Gen_ZONE)+
  labs(x = "N reads mapped") +
  theme_bw()
graph1.4

graph1.5 <- map.res %>% #filter(Espece_revision %in% c("Pb"))  %>% 
  ggplot(aes(x = mapped)) +
  geom_density() +
  #geom_vline(xintercept = map.res %>% filter(Espece_revision %in% c("Pb")) %>% pull(mapped) %>% median, 
  #           lty = "dashed" )+
  labs(x = "N reads mapped") +
  #facet_wrap(~ Pop)+
  theme_bw()
graph1.5


## Overall popmap

id.rem <- read.csv("../Beluga_PopStructure_Global/00_Data/01_Filtering.ref/duplicated.ID.to.remove.Bringloe.csv", stringsAsFactors = F)
write.table(map.res %>% filter(mapped_perc >= 96) %>%  # remove ID with mapped_perc < 96
              filter(Espece %in% "Delphinapterus") %>%   # remove narwhals
              filter(Lieu_echantillonnage %nin% "Riviere_Marralik") %>%  # remove Marralik belugas
              filter(ID %nin% id.rem$ID_GQ) %>%  # remove duplicated IDs
              # you should also remove duplicated IDs - those with lower read coverage (slice_max(as.numeric(as.character(mean_cov_ns))))
              select(ID) %>% mutate(Pop  = "NoPop"),
            file = file.path(get.value("info.path"), "popmap_good_align_popgen_Bringloe_230201.txt"),
            quote = FALSE, sep = "\t",
            row.names = F, col.names = F)

# Used for CSB ResDoc, includes ALL samples (narwhals, duplicates, Marralik)
# write.table(map.res %>% filter(mapped_perc >= 96) %>%  # remove ID with mapped_perc < 96
#               select(ID) %>% mutate(Pop  = "NoPop"),
#             file = file.path(get.value("info.path"), "popmap_good_align_all_TB_final_230123.txt"),
#             quote = FALSE, sep = "\t",
#             row.names = F, col.names = F)

map.res %>% #group_by(Espece_revision) %>% 
  summarise(MedianNread = median(total),
            MedianPerc = median(mapped_perc),
            MinPerc = min(mapped_perc),
            MaxPerc = max(mapped_perc))

plot(map.res[,c("total", "secondary", "mapped", "paired", "properly_paired", "singletons", "twith_itself_mate_mapped", "twith_mate_mapped_diff_chr")])

map.res %>% #filter(Espece_revision %in% c("Pb", "Pm")) %>%
  ggplot(aes(x = total, y = secondary, color = No_plaque)) +
  geom_smooth(method = "lm")+
  geom_point() +
  #facet_wrap(~Espece_revision)+
  theme_bw()

map.res %>% #filter(Espece_revision %in% c("Pb", "Pm")) %>% 
  ggplot(aes(x = secondary, y = singletons, col = No_plaque)) +
  #geom_smooth()+
  geom_point() +
  # facet_wrap(~Espece_revision)+
  theme_bw()

map.res %>% #filter(Espece_revision %in% c("Pb", "Pm")) %>% 
  ggplot(aes(x = total, y = singletons, col = No_plaque)) +
  geom_point()

map.res %>% #filter(Espece_revision %in% c("Pb", "Pm")) %>% 
  ggplot(aes(x = twith_itself_mate_mapped, y = singletons, col = No_plaque)) +
  geom_point()

# map.res %>% filter(Cat_Sample == "Sample",  # use this after performing alignment with ALL individuals
#                    between(Permapped, 0.95, 0.97)) %>% 
#   group_by(Espece, Gen_ZONE) %>% 
#   summarise(N = n()) %>% arrange(desc(N)) #%>% pull(N) %>% sum()

hist(map.res$total) 

quantile(map.res$mapped, c(0.01,0.05,0.1))# %>% summary()

# map.res %>%  %>% summarise(Mean = mean(Permapped))

# write_csv(map.res, file = file.path(get.value("demulti.log"), "Mapping_Draft_TB_final_230123.csv"))




# Stacks - Testing parameters - Ref-genome --------------------------------

# This part will run the ref_map.pl pipeline with a subset of samples,
# just to be sure that everything is alright

cmd <- paste("--samples", file.path(get.value("align.path"),"01_Bringloe_22xii22"),
             "--popmap", file.path(get.value("info.path"), "popmap.test_samples.txt"),
             "-o", file.path(get.value("test.ref.path"),"01_Bringloe_22xii22"),
             "-T", 32,
             "-X", "\"populations:-r 0.80 --genepop\"",
             "-X", "\"gstacks:-S .sorted.bam\"" # en espérant que ça fonctionne ...
)

stacks.pipelines <- "/media/genyoda/Extra_Storage/Projets/MOBELS/stacks2.55-scripts/"  # specify path for ref_map_new.pl (version with path for stacks 2 specified)

A <- system2(paste0(stacks.pipelines ,"ref_map_new.pl"), cmd, stdout=T, stderr=T)
A

# Stats

# Extract statistics

cmd1 <- paste(file.path(get.value("test.ref.path"), "01_Bringloe_22xii22", "populations.log.distribs"),
              "samples_per_loc_postfilters")

res1 <- system2(paste0(stack2.55_path, "/scripts/", "stacks-dist-extract"), cmd1, stdout = T)  # modify path for stacks 2.55 because default stacks in use is an older version

data.int <- res1[c(-1, -2)] %>% str_split(pattern = "\t")

data.int <-  data.frame(matrix(unlist(data.int), nrow= length(data.int), byrow = T))
names(data.int) <- c("n_samples", "n_loci")

n_loci <- as.numeric(as.character(data.int$n_loci)) %>%  sum()  # 206785 - 199895 loci TB final

# N SNP / locus

cmd2 <- paste(file.path(get.value("test.ref.path"), "01_Bringloe_22xii22", "populations.log.distribs"),
              "snps_per_loc_postfilters")

res2 <- system2(paste0(stack2.55_path, "/scripts/", "stacks-dist-extract"), cmd2, stdout = T)  # modify path for stacks 2.55 because default stacks in use is an older version

data.int <- res2[c(-1, -2)] %>% str_split(pattern = "\t")

data.int <-  data.frame(matrix(unlist(data.int), nrow= length(data.int), byrow = T))
names(data.int) <- c("n_snps", "n_loci")

data.int$n_snps <- as.numeric(as.character(data.int$n_snps))
data.int$n_loci <- as.numeric(as.character(data.int$n_loci))

no_div <- data.int %>% filter(n_snps == 0) %>% pull(n_loci)  %>% as.numeric()
n_loci_poly <- n_loci - no_div

cat("\nThere is", n_loci_poly, "polymorphic loci (r80) out of", n_loci, "loci")  # There is 74861 polymorphic loci (r80) out of 206785 loci
# TB final: There is 72627 polymorphic loci (r80) out of 199895 loci




# Gstacks - Ref - discovered SNPs -----------------------------------------

if(dir.exists("./00_Data/05b_Stacks.ref/01_Bringloe_22xii22/01_popgen") == F){
  dir.create("./00_Data/05b_Stacks.ref/01_Bringloe_22xii22/01_popgen")
}


cmd <- paste("-I", file.path(get.value("align.path"),"01_Bringloe_22xii22"),
             "-M", file.path(get.value("info.path"), "popmap_good_align_popgen_Bringloe_230201.txt"),  # changed from popmap_good_align_all_TB_final_230123
             "-O", file.path(get.value("stacks.ref.path"),"01_Bringloe_22xii22","01_popgen"),
             "-S", ".sorted.bam",
             "-t", 12)  # was 16

A <- system2(paste0(stack2.55_path,"gstacks"), cmd, stdout=T, stderr=T)  # default filter for PHRED-scaled mapping quality to consider a read (Q = 10)
A

cat(file = file.path(get.value("stacks.ref.log"), "01_Bringloe_22xii22", "01_popgen", "gstacks.ref.Bringloe.log"),
    "\n", cmd, "\n",
    A, # what to put in my file
    append= F, sep = "\n")

# If you want the big files to be ignored, run the following :

cat("catalog.*", "*.calls", "!.gitignore", sep = "\n",
    file = file.path(get.value("stacks.ref.path"), ".gitignore"))



# Check the distribution of ...

# N SNP / locus

cmd <- paste(file.path(get.value("stacks.ref.path"),"01_Bringloe_22xii22","01_popgen","gstacks.log.distribs"),  # check if gstacks.log.distribs is empty: if so, this command won't work
             "effective_coverages_per_sample")

res <- system2(paste0(stack2.55_path, "/scripts/", "stacks-dist-extract"), cmd, stdout = T)  # runs very fast
# Providing the number of loci in each sample, the number of forward reads used in the final loci for that sample (recall that gstacks takes the forward 
# reads assembled by ustacks, and finds the corresponding paired-end reads from which to build the paired-end contig), the mean depth of coverage, and a 
# weighted mean depth of coverage (weighting loci found in more samples higher)

data <- res[c(-1, -2, -3)] %>% str_split(pattern = "\t")

data <-  data.frame(matrix(unlist(data), nrow= length(data), byrow = T))
names(data) <- c("sample", "n_loci", "n_used_fw_reads", "mean_cov", "mean_cov_ns")

data <- data %>% left_join(pop.data, by = c("sample" = "ID_GQ"))
data %>% View()

data %>% group_by(Annee_echantillonnage) %>% summarise(N = n())
data %>% group_by(Region_echantillonnage) %>% summarise(N = n())
# pop.data %>% group_by(Pop) %>% summarise(N = n())

# find out which specimen is in map.res but not in data
which(map.res$ID %nin% data$sample)  # 26 samples removed
map.res$ID[which(map.res$ID %nin% data$sample)]
# "S_20_00626"     "S_20_00646"     "S_20_00647"     "S_20_00649_rep" "S_20_00650_rep" "S_20_00651_rep" "S_20_00669"     "S_20_00698"     "S_20_00700"     "S_20_00708"    
# "S_20_00726_rep" "S_20_00786"     "S_20_00835_rep" "S_20_00835"     "S_20_00849"     "S_20_00851"     "S_20_00852"     "S_20_00853_rep" "S_20_00854"     "S_20_00855_rep"
# "S_20_00856_rep" "S_20_00857"     "S_20_00860_rep" "S_20_00861"     "S_20_00862"     "S_20_00916"     "S_20_00917_rep" "S_20_00918"     "S_20_00921"     "S_20_01012"    
# "S_20_01013"     "S_20_01018"     "S_20_01019_rep" "S_20_01021_rep" "S_20_01032"     "S_20_01095"     "S_20_01096"     "S_20_01097_rep" "S_20_01098_rep" "S_20_01099_rep"
# "S_20_01213"     "S_20_01215"     "S_20_01216"     "S_20_01217"     "S_20_01218"     "S_20_01219"     "S_20_01220"     "S_20_01221"     "S_20_01222"     "S_20_01223"    
# "S_20_01773"     "S_20_01783"     "S_20_01786"     "S_20_01795"     "S_20_01893"     "S_20_01944_rep" "S_20_03229_rep" "S_20_03475"     "S_20_03476"     "S_20_03505"    
# "S_20_03542"     "S_20_03661"     "S_20_03662"     "S_20_03663"     "S_20_03665"
# Mapped perc < 96% 
# "S_20_00698" "S_20_00786" "S_20_01773" "S_20_01786" "S_20_01795" "S_20_01893" "S_20_03475" "S_20_03505" "S_20_03542"
# Duplicates (in id.rem)
# "S_20_00626"     "S_20_00646"     "S_20_00647"     "S_20_00649_rep" "S_20_00650_rep" "S_20_00651_rep" "S_20_00726_rep" "S_20_00835_rep" "S_20_00835"    
# "S_20_00849"     "S_20_00851"     "S_20_00852"     "S_20_00853_rep" "S_20_00854"     "S_20_00855_rep" "S_20_00856_rep" "S_20_00857"     "S_20_00860_rep"
# "S_20_00861"     "S_20_00862"     "S_20_00916"     "S_20_00917_rep" "S_20_00918"     "S_20_00921"     "S_20_01012"     "S_20_01013"     "S_20_01018"    
# "S_20_01019_rep" "S_20_01021_rep" "S_20_01095"     "S_20_01096"     "S_20_01097_rep" "S_20_01098_rep" "S_20_01099_rep" "S_20_01783"     "S_20_01944_rep"
# "S_20_03229_rep" "S_20_03476"     "S_20_03665"
# Marralik river
# "S_20_03661" "S_20_03662" "S_20_03663"
# Narwhals
# "S_20_00669" "S_20_00700" "S_20_00708" "S_20_01032" "S_20_01213" "S_20_01215" "S_20_01216" "S_20_01217" "S_20_01218" "S_20_01219" "S_20_01220" "S_20_01221"
# "S_20_01222" "S_20_01223"


library(ggridges)
graph2.1 <- data %>% 
  subset(Espece %in% "Delphinapterus") %>% 
  ggplot(aes(x = as.numeric(as.character(mean_cov_ns)), y  = Region_echantillonnage, fill = ..x..)) +
  ggridges::geom_density_ridges_gradient(alpha = 1/2, scale = 3) +
  scale_fill_viridis_c(name = "Mean coverage") +
  geom_vline(xintercept = c(5,10), lty = "dashed", col = "black") +
  scale_x_continuous(breaks = c(0,5,10,20,40,60,80, 100))  + #facet_wrap(~ Objective) +
  labs(x = "Mean coverage", title = "Mean coverage post gstacks in Belugas") +
  # theme_ridges() +
  theme_bw() +
  theme(legend.position = "right")

graph2.1

ggsave(filename = file.path(get.value("stacks.ref.log"),"01_Bringloe_22xii22","01_popgen","CoverageByPop_TB_final_230203.png"),
       plot = graph2.1,
       width = 15, height = 10, units = "in")

graph2.2 <- data %>%
  subset(Espece %in% "Delphinapterus") %>% 
  ggplot(aes(y = as.numeric(as.character(mean_cov_ns)), x = Region_echantillonnage, col = as.factor(Region_echantillonnage))) +
  geom_violin(col = "black") +
  geom_jitter(height = 0, alpha = 0.5) +
  #scale_y_continuous(breaks = c(0, 60000, 120000, 180000, 240000, 300000, 360000, 420000, 480000)) +
  #geom_vline(xintercept = 5, lty = "dashed", col = "darkgray") +
  geom_hline(yintercept = c(5,10), lty = "dashed", col = "darkgray") +
  #facet_wrap(~Gen_ZONE, nrow = 2) +
  labs(y = "Mean coverage", x = "Sampling region") +
  theme_bw() +
  theme(legend.position = "none")

graph2.2

ggsave(filename = file.path(get.value("stacks.ref.log"),"01_Bringloe_22xii22","01_popgen","CoverageByPop2_TB_final_230203.png"),
       plot = graph2.2,
       width = 14, height = 8, units = "in")

graph2.3 <- data %>%
  subset(Espece %in% "Delphinapterus") %>% 
  ggplot(aes(x = as.numeric(as.character(mean_cov_ns)), y = as.numeric(as.character(n_loci)), col = as.factor(Region_echantillonnage))) +
  geom_point(alpha = 0.5) +
  #scale_y_continuous(breaks = c(0, 60000, 120000, 180000, 240000, 300000, 360000, 420000, 480000)) +
  geom_vline(xintercept = c(5,10), lty = "dashed", col = "darkgray") +
  geom_hline(yintercept = 50000, lty = "dashed", col = "darkgray") +
  #facet_wrap(~Gen_ZONE, nrow = 2) +
  labs(x = "Mean coverage", y = "N loci", colour = "Sampling region") +
  theme_bw() +
  theme(legend.position = "right")

graph2.3

ggsave(filename = file.path(get.value("stacks.ref.log"),"01_Bringloe_22xii22","01_popgen","CoverageVSloci_TB_final_230203.png"),
       plot = graph2.3,
       width = 14, height = 9, units = "in")

graph2.4 <- data %>%
  subset(Espece %in% "Delphinapterus") %>% 
  ggplot(aes(y = as.numeric(as.character(mean_cov_ns)), x = as.factor(Annee_echantillonnage), col = as.factor(Region_echantillonnage))) +
  geom_boxplot(col = "black") +  
  geom_jitter(height = 0, alpha = 0.5) +
  #scale_y_continuous(breaks = c(0, 60000, 120000, 180000, 240000, 300000, 360000, 420000, 480000)) +
  #geom_vline(xintercept = 5, lty = "dashed", col = "darkgray") +
  geom_hline(yintercept = 5, lty = "dashed", col = "darkgray") +
  #facet_wrap(~Gen_ZONE, nrow = 2) +
  labs(y = "Mean coverage", x = "Sampling year", colour = "Sampling region") + 
  theme_bw() +
  theme(legend.position = "none")

graph2.4

ggsave(filename = file.path(get.value("stacks.ref.log"),"01_Bringloe_22xii22","01_popgen","CoverageVsYear_TB_final_230203.png"),
       plot = graph2.4,
       width = 14, height = 8, units = "in")

graph2.5 <- data %>% left_join(map.res %>% select(ID, total), by = c("sample" = "ID")) %>% 
  ggplot(aes(y = as.numeric(as.character(mean_cov_ns)), x = as.numeric(as.character(total)), col = as.factor(Espece))) +
  geom_point(alpha = 0.5) +
  scale_y_continuous(trans = "log10") +
  scale_x_continuous(trans = "log10") +
  geom_hline(yintercept = 5, lty = "dashed", col = "darkgray") +
  #geom_hline(yintercept = 60000, lty = "dashed", col = "darkgray") +
  #  facet_wrap(~Gen_ZONE, nrow = 2) +
  labs(x = "Nreads obtained", y = "Coverage") + 
  theme_bw() +
  theme(legend.position = "none")

graph2.5

ggsave(filename = file.path(get.value("stacks.ref.log"),"01_Bringloe_22xii22","01_popgen","CoverageVsNreads_TB_final_230203.png"),
       plot = graph2.5,
       width = 14, height = 8, units = "in")


data %>% mutate(cat_cov = ifelse(as.numeric(as.character(mean_cov_ns)) <3, "LOW", "OK")) %>% group_by(cat_cov) %>% summarise(N = n())
# TB popgenL LOW 10 vs OK 736
# TB final: LOW 13 vs OK 789

# Save the results

readr::write_csv(data, file.path(get.value("stacks.ref.log"),"01_Bringloe_22xii22","01_popgen","BelugaID_unique_noMarralik_RefGen_Bringloe_NreadsNloci.csv"))
