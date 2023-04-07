setwd("~/Desktop/R Working Directory/NanoString")
library("tidyverse")
library(car)
library(pandsas)
data.merge <- read_csv("Modified Data/data.merged.csv")
gene.list <- (colnames(data.merge[7:ncol(data.merge)]))

genotypes <- unique(data.merge$Genotype)
ages <- unique(data.merge$Age)
sex <- unique(data.merge$Sex)
dose <- unique(data.merge$Dose)
treatments <- c("VER", "LEV")
dd2 <- append(genotypes, c(ages, sex, dose, treatments))

test <- data.merge 
dd <-data.frame(matrix(nrow=nrow(data.merge), ncol=length(dd2)))
colnames(dd) <- dd2
dd <- dd  %>% mutate(across(everything(), replace_na, 0))
dd <- mutate(dd, Sample=data.merge$Sample, .before = '5xFAD')
test <- pivot_longer(test, cols=7:ncol(test), names_to="Gene", values_to="Expression")
test <- left_join(test, dd, by="Sample", relationship='many-to-many')
test <- test %>% group_by(across(2:6)) %>% nest()

for(i in 1:nrow(test)){
  test[[6]][[i]]
}





