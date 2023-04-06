setwd("~/Desktop/R Working Directory/NanoString")
library("tidyverse")
library(car)
library(broom)
data.merge <- read_csv("Modified Data/data.merged.csv")
treatment <- filter(data.merge, Treatment!="NO")
gene.list <- (colnames(data.merge[7:ncol(data.merge)]))
################################ Find the significant genes between treatments ###################################
treatment <- treatment %>% filter(Dose != "MED" | Dose !="LOW")
VER <- treatment %>% 
  filter(Treatment=="VER" & Genotype=="5xFAD") %>% 
  pivot_longer(cols=7:ncol(treatment), names_to="Gene", values_to="Expression")
VER$Dose[VER$Dose != "VEH"] <- c("Treated")
VER$Dose <- factor(VER$Dose,levels=c("VEH", "Treated"))

LEV <- treatment %>% 
  filter(Treatment=="LEV" & Genotype=="5xFAD") %>% 
  pivot_longer(cols=7:ncol(treatment), names_to="Gene", values_to="Expression")
LEV$Dose[LEV$Dose != "VEH"] <- c("Treated")
LEV$Dose <- factor(LEV$Dose,levels=c("VEH", "Treated"))

ver.model <- data.frame()
for (i in 1:length(gene.list)) {
  NSdata.Gene <- filter(VER, Gene==gene.list[i])
  fit <- lm(Expression ~ Dose, data=NSdata.Gene)
  beta <- summary(fit)$coeff[,1]
  pval <- summary(fit)$coeff[,4]
  X <- cbind(beta,pval,gene.list[i])
  ver.model <- rbind(ver.model, X)
}
colnames(ver.model) <- c("Beta", "pval","Gene")
ver.model$covariates <- gsub("^\\d+|\\d+$", "", rownames(ver.model))
ver.model$Beta <- as.numeric(ver.model$Beta)
ver.model$pval <- as.numeric(ver.model$pval)
ver.model <- filter(ver.model, covariates!="(Intercept)")
ver.sig.gene <- ver.model %>% filter(pval<=0.05)


lev.model <- data.frame()
for (i in 1:length(gene.list)) {
  NSdata.Gene <- filter(LEV, Gene==gene.list[i])
  fit <- glm(Dose ~ Expression, data=NSdata.Gene, family="binomial")
  beta <- summary(fit)$coeff[,1]
  pval <-  summary(fit)$coeff[,4]
  X <- cbind(beta,pval,gene.list[i])
  lev.model <- rbind(lev.model, X)
}
rm(X, pval, beta, NSdata.Gene)
colnames(lev.model) <- c("Beta", "pval","Gene")
lev.model$covariates <- gsub("^\\d+|\\d+$", "", rownames(lev.model))
lev.model$Beta <- as.numeric(lev.model$Beta)
lev.model$pval <- as.numeric(lev.model$pval)
lev.model <- filter(lev.model, covariates!="(Intercept)")
lev.sig.gene <- lev.model %>% filter(pval<=0.05)

write_csv(lev.sig.gene, "Results/LEV_Significant_Genes.csv")
write_csv(ver.sig.gene, "Results/VER_Significant_Genes.csv")

ggplot(data=ver.model) + 
  geom_point(aes(x=Beta, y=-log10(pval))) + 
  geom_hline(aes(yintercept=1.30103), color="red") +
  labs(title="Treatment Linear Regression VER", caption="N=36")
ggsave("Results/VER_Treatment_Regression.png")

######## Anovas of LEV and VER #####
VER <- treatment %>% 
  filter(Treatment=="VER" & Genotype=="5xFAD" & Dose=="VEH" |
           Treatment=="VER" & Genotype=="5xFAD" & Dose=="HIGH") %>% 
  pivot_longer(cols=7:ncol(treatment), names_to="Gene", values_to="Expression")

LEV <- treatment %>% filter(Treatment=="LEV" & Genotype=="5xFAD" & Dose=="VEH" |
                              Treatment=="LEV" & Genotype=="5xFAD" & Dose=="HIGH") %>% 
  pivot_longer(cols=7:ncol(treatment), names_to="Gene", values_to="Expression")

VER.Anova <- data.frame()
for(i in 1:length(gene.list)){
  temp <- filter(VER, Gene==gene.list[i])
  temp.aov <- aov(Expression ~ Dose, data=temp)
  pval <- tidy(temp.aov)$p.value[1]
  X <- cbind(gene.list[i], pval)
  VER.Anova <- rbind(VER.Anova, X)
}
colnames(VER.Anova) <- c("Gene", "VER_pval")
sig.VER.Anova <- filter(VER.Anova, VER_pval <0.05)


LEV.Anova <- data.frame()
for(i in 1:length(gene.list)){
  temp <- filter(LEV, Gene==gene.list[i])
  temp.aov <- aov(Expression ~ Dose, data=temp)
  pval <- tidy(temp.aov)$p.value[1]
  X <- cbind(gene.list[i], pval)
  LEV.Anova <- rbind(LEV.Anova, X)
}
colnames(LEV.Anova) <- c("Gene", "LEV_pval")
sig.LEV.Anova <- filter(LEV.Anova, LEV_pval <0.05)

compare <- full_join(LEV.Anova, VER.Anova, by="Gene")
compare.sig <- filter(compare, LEV_pval<0.05 & VER_pval <0.05)
LEV.sig <- filter(compare, LEV_pval<0.05)
VER.sig <- filter(compare, VER_pval <0.05)

##########################################Comparing Expression of 3 sig genes ##################################################################
treat.Lyrm9 <- treat %>% filter(Gene=="Lyrm9" & Treatment=="VER" & Genotype!="WT" & Dose!="LOW" & Dose!="MED" |
                                  Gene=="Lyrm9" & Treatment=="LEV" & Genotype!="WT" & Dose!="LOW" & Dose!="MED" )

treat.Lyrm9 %>% ggplot() + 
  geom_boxplot(aes(Dose, Expression, color=Treatment)) + 
  labs(title="Lyrm9 Expression in 5xFAD mice") +
  xlim("VEH", "HIGH")
ggsave("Results/Lyrm9_Treatment_5xFAD.png")

treat.Trp53bp2 <- treat %>% filter(Gene=="Trp53bp2" & Treatment=="VER" & Genotype!="WT" & Dose!="LOW" & Dose!="MED" |
                                     Gene=="Trp53bp2" & Treatment=="LEV" & Genotype!="WT" & Dose!="LOW" & Dose!="MED" )

treat.Trp53bp2 %>% ggplot() + 
  geom_boxplot(aes(Dose, Expression, color=Treatment)) + 
  labs(title="Trp53bp2 Expression in 5xFAD mice")+
  xlim("VEH", "HIGH")
ggsave("Results/Trp53bp2_Treatment_5xFAD.png")


treat.Rras <- treat %>% filter(Gene=="Rras" & Treatment=="VER" & Genotype!="WT" & Dose!="LOW" & Dose!="MED" |
                                 Gene=="Rras" & Treatment=="LEV" & Genotype!="WT" & Dose!="LOW" & Dose!="MED" )

treat.Rras %>% ggplot() + 
  geom_boxplot(aes(Dose, Expression, color=Treatment)) + 
  labs(title="Rras Expression in 5xFAD mice")+
  xlim("VEH", "HIGH")
ggsave("Results/Rras_Treatment_5xFAD.png")
####################### Look at pathways for each treatment ###################
library(gprofiler2)
ver.sig.gene <- arrange(ver.sig.gene, beta)
gost.ver <- gost(ver.sig.gene$Gene, organism="mmusculus", ordered_query=TRUE)
gost.ver.res <- gost.ver$result

lev.sig.gene <- arrange(lev.sig.gene, beta)
gost.lev <- gost(lev.sig.gene$Gene, organism="mmusculus", ordered_query=TRUE)
gost.lev.res <- gost.lev$result



##################################   Genes significantly effected by LEV in other genotypes ####
sig.LEV.gene <- data.merge %>% select(1:7, all_of(lev.gene))
sig.LEV.gene <- sig.LEV.gene[str_detect(sig.LEV.gene$Genotype, "HFD", negate=TRUE),]
sig.LEV.gene <- sig.LEV.gene %>% pivot_longer(cols=7:ncol(sig.LEV.gene), names_to="Gene", values_to="Expression")
sig.LEV.gene <- sig.LEV.gene %>% filter(Treatment=="NO" & Genotype!="WT")
genotypes <- unique(sig.LEV.gene$Genotype)

sig.LEV.gene.model <- data.frame()
for(i in 1:length(genotypes)){
  for(j in 1:length(lev.gene)){
  temp <- filter(sig.LEV.gene, Genotype==genotypes[i] | Genotype=="A/T.MIXED-WT")
  temp <- filter(temp, Gene==lev.gene[j])
  temp.model <- lm(Expression ~ Age * Genotype, data=temp)
  m<- summary(temp.model)
  beta <- summary(temp.model)$coefficients[,1]
  pval <- summary(temp.model)$coefficients[,4]
  X <- cbind(genotypes[i], lev.gene[j], beta,pval)
  sig.LEV.gene.model <- rbind(sig.LEV.gene.model, X)
}
}

colnames(sig.LEV.gene.model) <- c("Genotype","Gene", "Beta", "pval")
sig.LEV.gene.model$covariates <- gsub(pattern="^\\d+|\\d+$",  rownames(sig.LEV.gene.model), replacement="")
sig.LEV.gene.model$Beta <- as.numeric(sig.LEV.gene.model$Beta)
sig.LEV.gene.model$pval <- as.numeric(sig.LEV.gene.model$pval)
sig.LEV.gene.model <-filter(sig.LEV.gene.model, covariates!="(Intercept)")
sig.LEV <- filter(sig.LEV.gene.model, pval<=0.05)

lev.genotypes.counts <- data.frame()
for (i in 1:length(genotypes)){
  temp <- filter(sig.LEV, Genotype==genotypes[i])
  up <- temp %>% filter(Beta>0) %>% summarise(n())
  down <- temp %>% filter(Beta<0) %>% summarise(n())
  total <- temp %>% summarise(n())
  X <- cbind(genotypes[i], up, down, total)
  lev.genotypes.counts <- rbind(lev.genotypes.counts, X)
}
colnames(lev.genotypes.counts) <- c("Genotype", "num.Up", "num.Down", "Total")


lev.gene.counts <- data.frame()
for (i in 1:length(lev.gene)){
  temp <- filter(sig.LEV, Gene==lev.gene[i])
  up <- temp %>% filter(Beta>0) %>% summarise(n())
  down <- temp %>% filter(Beta<0) %>% summarise(n())
  total <- temp %>% summarise(n())
  X <- cbind(lev.gene[i], up, down, total)
  lev.gene.counts <- rbind(lev.gene.counts, X)
}
colnames(lev.gene.counts) <- c("Gene", "num.Up", "num.Down", "Total")
lev.sig.dose.gene <- filter(lev.sig.gene, covariates=="DoseTreated")
lev.gene.counts <- left_join(lev.sig.dose.gene, lev.gene.counts, by="Gene") %>% select(-covariates)
write_csv(lev.gene.counts, "Results/LEV_Significant_Gene_Genotype_Counts.csv")
####

###### Does the dose affect expression? ####
VER.dose  <- treatment %>% 
  filter(Treatment=="VER" & Genotype=="5xFAD")
VER.dose <- pivot_longer(VER.dose, cols=7:ncol(treatment), names_to="Gene", values_to="Expression")
VER.dose$Dose <- factor(VER.dose$Dose, levels=c("VEH", "LOW", "MED", "HIGH"))

VER.dose <- filter(VER.dose, Dose != "VEH")
VER.dose$Dose <- factor(VER.dose$Dose, levels=c("LOW", "MED", "HIGH"))
VER.aov <- data.frame()
for (i in 1:length(gene.list)) { #This should give all the genes where dose significantly affects the expression among mice that were given treatment. 
  NSdata.Gene <- filter(VER.dose, Gene==gene.list[i])
  aov <- Anova(lm(Expression ~ Dose, data=NSdata.Gene))
  pval <-  tidy(aov)$p.value[1]
  X <- cbind(gene.list[i], pval)
  VER.aov <- rbind(VER.aov, X)
}

colnames(VER.aov) <- c("Gene", "pval")
test1 <- filter(VER.aov, pval < 0.05)  #Shows that there are 47 genes where the dose is a significant factor in the expression level

VER.dose  <- treatment %>% 
  filter(Treatment=="VER" & Genotype=="5xFAD")
VER.dose <- pivot_longer(VER.dose, cols=7:ncol(treatment), names_to="Gene", values_to="Expression")
VER.dose$Dose <- factor(VER.dose$Dose, levels=c("VEH", "LOW", "MED", "HIGH"))
VER.aov <- data.frame()
for (i in 1:length(gene.list)) { #This should give all the genes where dose significantly affects the expression among mice that were given treatment. 
  NSdata.Gene <- filter(VER.dose, Gene==gene.list[i])
  aov <- Anova(lm(Expression ~ Dose, data=NSdata.Gene))
  pval <-  tidy(aov)$p.value[1]
  X <- cbind(gene.list[i], pval)
  VER.aov <- rbind(VER.aov, X)
}

colnames(VER.aov) <- c("Gene", "pval")
test2 <- filter(VER.aov, pval < 0.05)  #Shows that there are 47 genes where the dose is a significant factor in the expression level

LEV.dose  <- treatment %>% 
  filter(Treatment=="LEV" & Genotype=="5xFAD")
LEV.dose <- pivot_longer(LEV.dose, cols=7:ncol(treatment), names_to="Gene", values_to="Expression")
LEV.dose$Dose <- factor(LEV.dose$Dose, levels=c("VEH", "LOW", "MED", "HIGH"))

LEV.dose <- filter(LEV.dose, Dose != "VEH")
LEV.dose$Dose <- factor(LEV.dose$Dose, levels=c("LOW", "MED", "HIGH"))
LEV.aov <- data.frame()
for (i in 1:length(gene.list)) { #This should give all the genes where dose significantly affects the expression among mice that were given treatment. 
  NSdata.Gene <- filter(LEV.dose, Gene==gene.list[i])
  aov <- Anova(lm(Expression ~ Dose, data=NSdata.Gene))
  pval <-  tidy(aov)$p.value[1]
  X <- cbind(gene.list[i], pval)
  LEV.aov <- rbind(LEV.aov, X)
}

colnames(LEV.aov) <- c("Gene", "pval")
test1 <- filter(LEV.aov, pval < 0.05)  #Shows that there are 47 genes where the dose is a significant factor in the expression level

LEV.dose  <- treatment %>% 
  filter(Treatment=="LEV" & Genotype=="5xFAD")
LEV.dose <- pivot_longer(LEV.dose, cols=7:ncol(treatment), names_to="Gene", values_to="Expression")
LEV.dose$Dose <- factor(LEV.dose$Dose, levels=c("VEH", "LOW", "MED", "HIGH"))
LEV.aov <- data.frame()
for (i in 1:length(gene.list)) { #This should give all the genes where dose significantly affects the expression among mice that were given treatment. 
  NSdata.Gene <- filter(LEV.dose, Gene==gene.list[i])
  aov <- Anova(lm(Expression ~ Dose, data=NSdata.Gene))
  pval <-  tidy(aov)$p.value[1]
  X <- cbind(gene.list[i], pval)
  LEV.aov <- rbind(LEV.aov, X)
}

colnames(LEV.aov) <- c("Gene", "pval")
test2 <- filter(LEV.aov, pval < 0.05)


##### Keeping Dose in mind, which genes are changed by the treatment? ###
VER <- treatment %>% 
  filter(Treatment=="VER" & Genotype=="5xFAD") %>%
  pivot_longer(cols=7:ncol(treatment), names_to="Gene", values_to="Expression")
VER$Dose <- factor(VER$Dose,levels=c("VEH", "LOW", "MED", "HIGH"))

LEV <- treatment %>% filter(Treatment=="LEV" & Genotype=="5xFAD") %>% 
  pivot_longer(cols=7:ncol(treatment), names_to="Gene", values_to="Expression")
LEV$Dose <- factor(LEV$Dose, levels=c("VEH", "LOW", "MED", "HIGH"))

ver.model <- data.frame()
for (i in 1:length(gene.list)) {
  NSdata.Gene <- filter(VER, Gene==gene.list[i])
  fit <- lm(Expression ~ Dose + Sex, data=NSdata.Gene)
  beta <- summary(fit)$coeff[,1]
  pval <- summary(fit)$coeff[,4]
  X <- cbind(beta,pval,gene.list[i])
  ver.model <- rbind(ver.model, X)
}
colnames(ver.model) <- c("Beta", "pval","Gene")
ver.model$covariates <- gsub("^\\d+|\\d+$", "", rownames(ver.model))
ver.model$Beta <- as.numeric(ver.model$Beta)
ver.model$pval <- as.numeric(ver.model$pval)
ver.model <- filter(ver.model, covariates!="(Intercept)")
ver.sig.gene <- ver.model %>% filter(pval<=0.05)
ver.model.wide <- ver.model %>% pivot_wider(values_from=c("Beta", "pval"), names_from=covariates)
write_csv(ver.model.wide, "Results/VER_Dose_Model.csv")


lev.model <- data.frame()
for (i in 1:length(gene.list)) {
  NSdata.Gene <- filter(LEV, Gene==gene.list[i])
  fit <- lm(Expression ~ Dose, data=NSdata.Gene)
  beta <- summary(fit)$coeff[,1]
  pval <-  summary(fit)$coeff[,4]
  X <- cbind(beta,pval,gene.list[i])
  lev.model <- rbind(lev.model, X)
}
rm(X, pval, beta, NSdata.Gene)
colnames(lev.model) <- c("Beta", "pval","Gene")
lev.model$covariates <- gsub("^\\d+|\\d+$", "", rownames(lev.model))
lev.model$Beta <- as.numeric(lev.model$Beta)
lev.model$pval <- as.numeric(lev.model$pval)
lev.model <- filter(lev.model, covariates!="(Intercept)")
lev.sig.gene <- lev.model %>% filter(pval<=0.05)
lev.genes <- unique(lev.sig.gene $Gene)

summary(lm(Expression ~ Sample, data=VER))
