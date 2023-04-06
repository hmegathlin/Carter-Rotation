library(tidyverse)
library(corrr)
library(ggcorrplot)
library(factoextra)
library(FactoMineR)

data.merge <- read_csv("Modified Data/data.merged.csv")
treatment <- data.merge[str_detect(data.merge$Genotype, "HFD", negate=TRUE),]
treatment <- filter(data.merge, Treatment!="NO")
treatment <- data.frame(treatment)
rownames(treatment) <- treatment$Sample
meta <- treatment %>% select(1:6)
treatment <- select(treatment, 2:ncol(treatment), -Age)


pca <- PCA(treatment, quali.sup=1:4)

test <- data.frame(pca$ind$coord)
test$Sample <- rownames(test)
test <- left_join(test, meta)

ggplot(test) + geom_point(aes(x=Dim.1, y=Dim.2, color=Treatment, shape=Genotype))
 ############ VER PCA #########
treatment <- data.merge[str_detect(data.merge$Genotype, "HFD", negate=TRUE),]
VER <- filter(data.merge, Treatment=="VER")
VER <- data.frame(VER)
rownames(VER) <- VER$Sample
meta <- VER %>% select(1:6)
VER <- select(VER, 2:ncol(VER), -Age, -Treatment)


pca <- PCA(VER, quali.sup=1:3)

test <- data.frame(pca$ind$coord)
test$Sample <- rownames(test)
test <- left_join(test, meta)

test$Dose <- factor(test$Dose, levels=c("VEH", "LOW", "MED", "HIGH"))
ggplot(test) + geom_point(aes(x=Dim.1, y=Dim.2, alpha=Dose, color=Sex, shape=Genotype))


##########################
data <- data.frame(data.merge)
data <- filter(data, Treatment == "NO")
data <- data[str_detect(data$Genotype, "HFD", negate=TRUE),]

rownames(data) <- data$Sample
meta <- data %>% select(1:4)                          
data <- data %>% select(2:ncol(data), -Dose, -Treatment)
pca <- PCA(data, quali.sup=1:3)
fviz_pca_ind(pca)

pca.data <- data.frame(pca$ind$coord)
pca.data$Sample <- rownames(pca.data)
pca.data <- left_join(pca.data, meta)

pca.data.males <- filter(pca.data, Sex=="M")
pca.data.females <- filter(pca.data, Sex=="F")

ggplot() + 
  geom_point(mapping=aes(Dim.1, Dim.2, color=Genotype), data=pca.data.males) + scale_color_viridis_d(option="D") +
  new_scale_color() +
  geom_point(aes(Dim.1, Dim.2), data=pca.data.females) + scale_colour_viridis_d(option="A")



ggplot(pca.data, aes(group=Sex)) + geom_point(aes(x=Dim.1, y=Dim.2, color=Genotype, shape=Sex, alpha=Age)) + 
  xlab("Dim 1 (21.9%)") + ylab("Dim 2 (8.9%)")
ggsave("Results/Genotype_PCA.png", width=8, height=5)
