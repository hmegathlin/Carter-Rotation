setwd("~/Desktop/R Working Directory/NanoString")
library(tidyverse)


############################################Load and organize LEV Data################################################
meta.LEV <- read_csv("Raw Data/metadata_NSdata_LEV_samples.csv") %>% 
  select(-c(`Duplex ID`, `Cohort ID`)) %>%
  mutate(Age=6, Treatment="LEV")
colnames(meta.LEV) <- c("Sample", "Sex", "Age", "Genotype", "Dose", "Treatment")
LEV <- read_csv("Modified Data/Transpose_Nanostringdata_logNormalized_MergedTg_Levi_129Samples.csv")
colnames <- colnames(LEV[,2:ncol(LEV)])
colnames(LEV) <- c("Sample", colnames)
meta.LEV$Dose[meta.LEV$Dose==0] <- "VEH"
meta.LEV$Dose[meta.LEV$Dose==10] <- "LOW"
meta.LEV$Dose[meta.LEV$Dose==30] <- "MED"
meta.LEV$Dose[meta.LEV$Dose==56] <- "HIGH"
data.LEV <- full_join(meta.LEV, LEV)
write_csv(data.LEV, "Modified Data/LEV_meta_data_merged.csv")

############################################Load and organize VER Data################################################
meta.VER <- read_csv("Raw Data/metadata_NSdata_VER_samples.csv",col_types = cols(...1 = col_skip())) %>% 
  select(-c(Cohort, PEN.ID, MNBF.ID)) %>% 
  mutate(Age=6, Treatment="VER") %>%
  relocate(Age, .before="Genotype")
meta.VER$Genotype[meta.VER$Genotype=="TG"] <- "Tg"
VER <- read_csv("Modified Data/Transpose_Nanostringdata_logNormalized_MergedTg_VER_60Samples.csv")
colnames <- colnames(VER[,2:ncol(VER)])
colnames(VER) <- c("Sample", colnames)
data.VER <- full_join(meta.VER, VER)
write_csv(data.VER, "Modified Data/VER_meta_data_merged.csv")


############################################Load and organize No Treatment Data################################################
NT <- read_csv("Modified Data/Transpose_Nanostringdata_logNormalized_BatchCorrected_MergedTg_All_937Samples_02Feb23.csv")
colnames <- colnames(NT[,2:ncol(NT)])
colnames(NT) <- c("Sample",colnames)
meta.NT <- read_csv("Raw Data/metadata_NSdata_937samples.csv", col_types = cols(...1 = col_skip(), BindingDensity = col_skip(), RLF = col_skip(), Batch = col_skip(), 
                                          Cohort = col_skip()))
colnames(meta.NT) <- c("Sample", "Sex", "Age", "Genotype")
meta.NT <- mutate(meta.NT, Treatment="NO")

data.NT <- full_join(meta.NT, NT)

########## Merge all data #################################################################################################
data.merge <- full_join(data.VER, data.NT) %>% full_join(data.LEV)
data.merge <- data.merge %>% replace_na(list(Dose="NO"))

apply(is.na(data.merge), 2, which) #gives where NAs are 
data.merge <- filter(data.merge, Sample!="GT19-03096") #Remove the sample that had no data. 

na.col <- data.merge %>% select_if(~ any(is.na(.))) #Finds columns that have NAs in them. 
na.col <- colnames(na.col) #pull the names of the columns that have NAs
data.merge <- data.merge %>% select(!all_of(na.col)) #Selects only the genes that every sample is measured for. 
write_csv(data.merge, "Modified Data/data.merged.csv")

######################################## Normalize Genotype Names ########################
unique(data.merge$Genotype)

data.merge$Genotype <- gsub("APOE4Trem2", data.merge$Genotype, replacement="A/T")
data.merge$Genotype <- gsub("hATA", data.merge$Genotype, replacement="A/T")
data.merge$Genotype <- gsub("LOAD1", data.merge$Genotype, replacement="A/T")
data.merge$Genotype <- gsub("Tg", data.merge$Genotype, replacement="5xFAD")
write_csv(data.merge, "Modified Data/data.merged.csv")
