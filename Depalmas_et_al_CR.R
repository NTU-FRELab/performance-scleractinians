#' ---
#' title: "Coral_Strategies"
#' author: "Stéphane De Palmas"
#' date: "2025-01-23"
#' output:
#'   html_document:
#'     toc: true
#'     toc_depth: 2
#'   pdf_document:
#'     toc: true
#'     toc_depth: '2'
#' ---
#' 
#' Data and script for replicating analyses from De Palmas et al. 2025 (https://doi.org/10.1007/s00338-025-02617-w).
#' 
#' If used in full or in part, please cite the original publication: Stéphane De Palmas, Qi Chen, Arnaud Guerbet, Tsai-Hsuan Tony Hsu, Yunli Eric Hsieh, Yuting Vicky Lin, Pei-Ling Wang, Nicolas Sturaro, Vianney Denis (2025) Energy acquisition and allocation strategies in scleractinian corals: insights from intraspecific trait variability. Coral Reefs. https://doi.org/10.1007/s00338-025-02617-w
#' 
#' The full dataset is also available on [https://doi.org/10.5061/dryad.xwdbrv1hw](https://doi.org/10.5061/dryad.xwdbrv1hw). This is a study from [FRElab](https://www.dipintothereef.com/).
#' 
#' ## **Generate R Script**
#' 

#' 
#' ## **Packages**
## ----eval = T, echo = T, message = F, warning = F-----------------------------
library(car)
library(cluster)
library(corrplot)
library(dplyr)
library(factoextra)
library(ggplot2)
library(ggpubr)
library(multcompView)
library(multcomp)
library(NbClust)
library(devtools)
library(pairwiseAdonis)
library(PMCMRplus)
library(PupillometryR)
library(rcompanion)
library(rlang)
library(tidyr)
library(vegan)
library(tidyverse)
library(ggraph)
library(rstatix)
library(lme4)
library(lmerTest)
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

#' 
#' ## **Importing the raw dataset**
## ----eval = T, echo = T, message = F, warning = F, results = "hide"-----------
#Import the raw data
data1 <- read.csv ("Data/DePalmas_et_al_CR_raw.csv", header = TRUE, sep=",", dec=".")
#Set the date format
data1$date <- as.Date(data1$date, "%d/%m/%Y")
#Display "delta" symbol for isotopic analyses
Sys.setlocale("LC_ALL","Greek") 

#' 
#' ## **Trait calculations from the raw data** {.tabset}
#' ### **Skeletal measurements (S)**
## ----eval = T, echo = T, message = F, warning = F, results = "hide"-----------
# Skeleton density 
data1$dens <- data1$Weight/data1$Volume 

#' 
#' ### **Isotopic measurements (I)**
## ----eval = T, echo = T, message = F, warning = F, results = "hide"-----------
# Host (H) and Symbiodiniaceae (S) carbon and nitrogen Isotopes 
carbon_title = expression(paste(delta^{13}, "C (\u2030)"))
nitrogen_title = expression(paste(delta^{15}, "N (\u2030)"))
data1$H15N <- data1$Hdelta15N 
data1$H13C <- data1$Hdelta13C 
data1$S15N <- data1$Zdelta15N 
data1$S13C <- data1$Zdelta13C 
data1$HCN <- data1$H_C/data1$H_N
data1$SCN <- data1$Z_C/data1$Z_N

#' 
#' ### **Other organismal measurements (O)**
## ----eval = T, echo = T, message = F, warning = F, results = "hide"-----------
# Tissue biomass (AFDW, mg) - used to standardize all the organismal traits
# w.fil (weight of the filter), w.fil.tis (weight of the filter with tissue), w.fil.tis.bur (weight of the filter with tissue after burning), Vi (total volume of slurry), Vbio (volume used for biomass) 
tot.AFDW <- (data1$W.fil.tis-data1$W.fil.tis.bur)*(data1$Vi/data1$V.bio)*1000 

# Symbiodiniaceae count (Sc)
col.countS <- paste("countS.",rep(1:5),sep="") # column with Symbiodiniaceae measurement replicates
average.countS <- rowMeans(subset(data1, select = col.countS), na.rm = TRUE) # add new column with average Symbiodiniaceae count
tot.countS <- average.countS*10000*data1$Vi # total count in our Vi
data1$Sc <- tot.countS/tot.AFDW  # standardization by ash free dry weight (counts per mgAFDW)

# cnidocyte count (Cc)
col.countC <- paste("countC.",rep(1:5),sep="") # column with cnidocyte measurement
average.countC <- rowMeans(subset(data1, select = col.countC), na.rm = TRUE) # add new column with average cnidocyte count
tot.countC <- average.countC*10000*data1$Vi # total count in our Vi
data1$Cc <-tot.countC/tot.AFDW  # standardization by ash free dry weight

# chlorophyll concentrations from the literature
# Use trichromatic equations to convert absorbances into ug mL-1
# chla = 11.85 * (A664 - A750) - 1.54 * (A647 - A750) - 0.08 * (A630 - A750),
# chlb = 21.03 * (A647 - A750) - 5.43 * (A664 - A750) - 2.66 * (A630 - A750),
# chlc = 24.52 * (A630 - A750) - 1.67 * (A664 - A750) - 7.60 * (A647 - A750),
# chlorophyll in all tissue
chlaT <- (11.85*(data1$A664T-data1$A750T)-1.54*(data1$A647T-data1$A750T)-0.08*(data1$A630T-data1$A750T)) *data1$V.extraction.T* data1$Vi /  data1$V.chl.slurry
chlbT <- (21.03*(data1$A647T-data1$A750T)-5.43*(data1$A664T-data1$A750T)-2.66*(data1$A630T-data1$A750T)) *data1$V.extraction.T* data1$Vi /  data1$V.chl.slurry
chlcT <- (24.52*(data1$A630T-data1$A750T)-1.67*(data1$A664T-data1$A750T)-7.60*(data1$A647T-data1$A750T)) *data1$V.extraction.T* data1$Vi /  data1$V.chl.slurry
# chlorophyll in all skeleton
chlaS <- (11.85*(data1$A664S-data1$A750S)-1.54*(data1$A647S-data1$A750S)-0.08*(data1$A630S-data1$A750S)) *data1$V.extraction.S
chlbS <- (21.03*(data1$A647S-data1$A750S)-5.43*(data1$A664S-data1$A750S)-2.66*(data1$A630S-data1$A750S)) *data1$V.extraction.S
chlcS <- (24.52*(data1$A630S-data1$A750S)-1.67*(data1$A664S-data1$A750S)-7.60*(data1$A647S-data1$A750S)) *data1$V.extraction.S
# adding chlorophylls in tissue and skeleton
chla <- chlaT+chlaS
chlb <- chlbT+chlbS
chlc <- chlcT+chlcS
# chlorophyll concentrations standardization
data1$chla <- chla/tot.AFDW # by afdw
data1$chlb <- chlb/tot.AFDW # by afdw
# negative values of chlorophyll b should be 0
data1$chlb[data1$chlb<0] <- 0
data1$chlc <- chlc/tot.AFDW # by afdw
# ratio chl c / chl a
data1$chlca = data1$chlc/data1$chla 

# protein concentrations 
abs.prot.zoox1.corr <- data1$abs.prot.zoox1-data1$abs.prot.blank # correction absorbance with blank
abs.prot.zoox2.corr <-data1$abs.prot.zoox2-data1$abs.prot.blank # correction absorbance with blank
abs.prot.host1.corr <- data1$abs.prot.host1-data1$abs.prot.blank # correction absorbance with blank
abs.prot.host2.corr <- data1$abs.prot.host2-data1$abs.prot.blank # correction absorbance with blank
prot.zoox.1 <- data1$STCX4*abs.prot.zoox1.corr^4+data1$STCX3*abs.prot.zoox1.corr^3+data1$STCX2*abs.prot.zoox1.corr^2+data1$STCX1*abs.prot.zoox1.corr+data1$STCE
prot.zoox.2 <- data1$STCX4*abs.prot.zoox2.corr^4+data1$STCX3*abs.prot.zoox2.corr^3+data1$STCX2*abs.prot.zoox2.corr^2+data1$STCX1*abs.prot.zoox2.corr+data1$STCE
prot.host.1 <- data1$STCX4*abs.prot.host1.corr^4+data1$STCX3*abs.prot.host1.corr^3+data1$STCX2*abs.prot.host1.corr^2+data1$STCX1*abs.prot.host1.corr+data1$STCE
prot.host.2 <- data1$STCX4*abs.prot.host2.corr^4+data1$STCX3*abs.prot.host2.corr^3+data1$STCX2*abs.prot.host2.corr^2+data1$STCX1*abs.prot.host2.corr+data1$STCE
# Host protein was measured with a different curve and blank (host and Symbiodiniaceae on different plates)
# corrected blank
abs.prot.host1.corr1 <- data1$abs.prot.host1-data1$abs.prot.blank.b # correction absorbance with blank
abs.prot.host2.corr1 <- data1$abs.prot.host2-data1$abs.prot.blank.b # correction absorbance with blank
# corrected absorbance value
prot.host.1b <- data1$STCX4b*abs.prot.host1.corr1^4+data1$STCX3b*abs.prot.host1.corr1^3+data1$STCX2b*abs.prot.host1.corr1^2+data1$STCX1b*abs.prot.host1.corr1+data1$STCEb
prot.host.2b <- data1$STCX4b*abs.prot.host2.corr1^4+data1$STCX3b*abs.prot.host2.corr1^3+data1$STCX2b*abs.prot.host2.corr1^2+data1$STCX1b*abs.prot.host2.corr1+data1$STCEb
# Replace the 13 first values with the corrected ones, and average on 2 values 
prot.host.1[1:13] <- prot.host.1b[1:13]
prot.host.2[1:13] <- prot.host.2b[1:13]
# total protein
prot.zoox <- apply(rbind(prot.zoox.1,prot.zoox.2),2,mean,na.rm = TRUE)*data1$SDS1V*data1$Vi/data1$V.prot.slurry
prot.host <- apply(rbind(prot.host.1,prot.host.2),2,mean,na.rm = TRUE)*data1$Vi
# concentration proteins
data1$Spro <- prot.zoox/tot.AFDW # by afdw
data1$Hpro <- prot.host/tot.AFDW # by afdw

#' 
#' ## **Preparation for multivariate and univariate analyses**
## ----eval = T, echo = T, message = F, warning = F, results = "hide"-----------
# preparation for statistics and multivariate analyses. 
# we create a  data frame and fill it with the information we need
resdata <- data.frame(matrix(ncol = 0, nrow = 196))
resdata$spe = data1$spe
resdata$dens = data1$dens
resdata$Sc = data1$Sc
resdata$Cc = data1$Cc
resdata$chla = data1$chla
resdata$chlb = data1$chlb
resdata$chlc = data1$chlc
resdata$chlca = data1$chlca
resdata$Spro = data1$Spro
resdata$Hpro = data1$Hpro
resdata$H13C = data1$H13C
resdata$H15N = data1$H15N
resdata$S13C = data1$S13C
resdata$S15N = data1$S15N
resdata$HCN = data1$HCN
resdata$SCN = data1$SCN

#we merge environmental data into it for now (only to remove the incomplete data, but eventually we will separate traits and environmental data later)
resdata$SST = data1$SST_mean
resdata$SST_sd = data1$SST_sd
resdata$light = data1$light_mean
resdata$light_sd = data1$light_sd
resdata$spe <- as.factor(resdata$spe) # to specify that species is a factor

#dataset with complete data only
resdata.full = resdata[complete.cases(resdata),] # keep only the complete data with no NA

# only keep numeric data (without the species factor) for later analyses and subset into traits and environmental data
resdata.fully = subset(resdata.full, select = c(dens, Sc, Cc, chla, chlb, chlc, chlca, Spro, Hpro, H13C,H15N, S13C, S15N, HCN, SCN))
envdata.fully = subset(resdata.full, select = c(SST, SST_sd, light, light_sd))

### Important notice: note that we have the same dataset with NA (resdata) and without NA (resdata.full). For analyses that were sensitive to the presence of NA, we use the later in our analyses. 

#' 
#' ## **Statistical analyses** {.tabset}
#' ### **Trait correlations**
## ----eval = T, echo = T, message = F, warning = F, results = "hide", fig.show='hide'----
# VIF check
model <- lm(resdata$dens ~ resdata$Sc + resdata$Cc + resdata$chla + resdata$chlb + resdata$chlc + resdata$chlca + resdata$Spro + resdata$Hpro + resdata$H13C + resdata$H15N + resdata$S13C + resdata$S15N + resdata$HCN + resdata$SCN, data = resdata)
summary(model)
vif(model)
# correlation matrix
matrice = cor(resdata.fully)
res1 = cor.mtest(resdata.fully, conf.level = .95) # significance test
p_value = res1$p
row.names(p_value) = row.names(matrice)
colnames(p_value) = colnames(matrice)
p_value[which(p_value > 0.05)] = "N.S."

col = colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(matrice, method = "color", col = col(200),
         type = "upper",tl.cex=1, order = "hclust",
         addCoef.col = "black",
         tl.col = "black", tl.srt = 90,
         p.mat=res1$p, sig.level = 0.05,
         diag = FALSE, number.cex=0.6)
# Notes : matrix shows 2 sets of correlated traits: H13C/S13C and chla/chlc while VIF shows which trait has higher correlation values
# Subset is done by removing trait with highest VIF one by one 
subsetdata1 = subset(resdata.fully, select=-c(chla))
### VIF check
model1 <- lm(resdata$dens ~ resdata$Sc + resdata$Cc + resdata$chlb + resdata$chlc + resdata$chlca + resdata$Spro + resdata$Hpro + resdata$H13C + resdata$S13C +resdata$H15N + resdata$S15N + resdata$HCN + resdata$SCN, data = resdata)
summary(model1)
vif(model1)
# correlation matrix on subset data1
matrice = cor(subsetdata1)
res1 = cor.mtest(subsetdata1, conf.level = .95) # significance test
p_value = res1$p
row.names(p_value) = row.names(matrice)
colnames(p_value) = colnames(matrice)
p_value[which(p_value > 0.05)] = "N.S."
col = colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(matrice, method = "color", col = col(200),
         type = "upper",tl.cex=1, order = "hclust",
         addCoef.col = "black",
         tl.col = "black", tl.srt = 90,
         p.mat=res1$p,sig.level = 0.05,
         diag = FALSE, number.cex=0.6)
# subset the second trait
subsetdata2 = subset(resdata.fully, select = -c(chla, H13C))
### VIF check
model2 <- lm(resdata$dens ~ resdata$Sc + resdata$Cc + resdata$chlb + resdata$chlc + resdata$chlca + resdata$Spro + resdata$Hpro + resdata$S13C + resdata$H15N + resdata$S15N + resdata$HCN + resdata$SCN, data = resdata)
summary(model2)
vif(model2)
# correlation matrix on subset data1
matrice = cor(subsetdata2)
res1 = cor.mtest(subsetdata2, conf.level = .95) # significance test
p_value = res1$p
row.names(p_value) = row.names(matrice)
colnames(p_value) = colnames(matrice)
p_value[which(p_value > 0.05)] = "N.S."

col = colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(matrice, method = "color", col = col(200),
         type = "upper",tl.cex=1, order = "hclust",
         addCoef.col = "black",
         tl.col = "black", tl.srt = 90,
         p.mat=res1$p,sig.level = 0.05,
         diag = FALSE, number.cex=0.6)

#' 
#' ### **PERMANOVA ~ species**
## ----eval = T, echo = T, message = F, warning = F, results = "hide", fig.show='hide'----
#  PERMANOVA
# We first need to transform the data 
stand.try<-decostand((subsetdata2), method= "range")
distance <- vegdist(as.matrix(stand.try))
# betadisper
mod <- betadisper(distance, resdata.full$spe)
anova(mod)
boxplot(mod) # visualization of the dispersion
TukeyHSD(mod)
# PERMANOVA and post hoc to find where are the differences
permanova.species<-adonis2(subsetdata2~resdata.full$spe, permutations = 999, method = 'euclidean', na.rm = TRUE)
permanova.species
pairwise.adonis(distance,
                resdata.full$spe, 
                p.adjust.m = "bonferroni",
                reduce = NULL,
                perm = 999)

#' 
#' ### **Cluster analyses**
## ----eval = T, echo = T, message = F, warning = F, results = "hide", fig.show='hide'----
# testing a wide range of algorithms using NbClust
set.seed(1234)
t <- NbClust(stand.try, distance = "euclidean", min.nc = 2, max.nc = 15, method = 'kmeans', index = "all", alphaBeale = 0.90)
part <- as.factor(t$Best.partition) # we can set the best partition (here k=3) as factor to replace the species factor

# Sub-setting by species and by partition 
# by species
resdata.full.amur<- resdata.full[resdata.full$spe == 'A. muricata', ]
resdata.full.cmic<- resdata.full[resdata.full$spe == 'C. microphthalma', ]
resdata.full.ipal<- resdata.full[resdata.full$spe == 'I. palifera', ]
resdata.full.plut<- resdata.full[resdata.full$spe == 'P. lutea', ]
resdata.full.ppro<- resdata.full[resdata.full$spe == 'P. profundacella', ]
resdata.full.pspe<- resdata.full[resdata.full$spe == 'P. speciosa', ]
resdata.full.spis<- resdata.full[resdata.full$spe == 'S. pistillata', ]

# by partition
partition.table<- as_factor(c(part)) 
resdata.full$partition<-partition.table
resdata.full.partition1<- resdata.full[resdata.full$partition == '1', ]
resdata.full.partition2<- resdata.full[resdata.full$partition == '2', ]
resdata.full.partition3<- resdata.full[resdata.full$partition == '3', ]
write.csv(resdata.full)

# Summarizing the data for additional analyses (mean, variance and median for all traits according to SPECIES and PARTITION)
# by species
resdata.mv = as_tibble(resdata[-1]) #just need to remove the first column with factor
Trait_mean = resdata.mv %>% group_by(resdata$spe) %>% summarise_all (mean, na.rm = T)
Trait_var = resdata.mv %>% group_by(resdata$spe) %>% summarise_all (var, na.rm = T)
Trait_median = resdata.mv %>% group_by(resdata$spe) %>% summarise_all (median, na.rm = T)

# by partition
resdata.mv.partition = as_tibble(resdata.full[,-c(1, 21)])#just need to remove the first and last column with factors
Trait_mean_partition = resdata.mv.partition %>% group_by(resdata.full$partition) %>% summarise_all (mean, na.rm = T)
Trait_var_partition = resdata.mv.partition %>% group_by(resdata.full$partition) %>% summarise_all (var, na.rm = T)
Trait_median_partition = resdata.mv.partition %>% group_by(resdata.full$partition) %>% summarise_all (median, na.rm = T)

#' 
#' ### **PERMANOVA ~ clusters**
## ----eval = T, echo = T, message = F, warning = F, results = "hide", fig.show='hide'----

permanova.partition<-adonis2(subsetdata2~part, permutations = 999, method = 'euclidean', na.rm = TRUE)
permanova.partition
pairwise.adonis(distance,
                part, 
                p.adjust.m = "bonferroni",
                reduce = NULL,
                perm = 999)

#' 
#' ### **Trait differences among clusters**
## ----eval = T, echo = T, message = F, warning = F, results = "hide"-----------
# UNIVARIATE STATISTICS: Our goal is to test variance and mean differences between traits for partitions (the following code could also be adapted to test variance and mean differences among species)
# First we tested homogeneity of variance (Levene's test) and normality of data (Shapiro's test).
# Whenever homogeneity of variance AND normality of data were verified we proceed to an ANOVA followed by post hoc test (Tukey's test)
# Whenever homogeneity of variance was violated but the data were normal, we performed a Games-Howell test.
# Whenever homogeneity of variance and normality were violated we proceed to mean comparisons using a Kruskal test followed by a Wilcoxon post hoc
# If homogeneity of variance was violated, we choose to perform an ANOVA on the residuals to compare differences of variance among partition
# Then we proceed to test the  
# Please not that we also use Compact Letter Display to visualize relative differences of VARIANCE between partitions (lowercase, in red) and differences of MEAN between partitions (uppercase, in black)

### Skeleton density (dens) 
# test on variance difference between partitions
leveneTest(dens ~ partition, resdata.full) # non significant => variance are equal 
# test on mean difference between partitions
# normality of the distribution
shapiro.test(resdata.full$dens) # significant => data not normal
kruskal.test(dens ~ partition, resdata.full)
wilcox.dens.part <- pairwise.wilcox.test(resdata.full$dens, resdata.full$partition, p.adjust.method = "bonferroni")
## post hoc display
cld.dens.mean.part <- as.data.frame.list(multcompLetters(fullPTable(wilcox.dens.part$p.value), compare="<", threshold=0.05,Letters= LETTERS, reversed = FALSE))

### Symbiodiniaceae count (Sc) 
# test on variance difference between partitions
leveneTest(Sc ~ partition, resdata.full) # significant => variance are not equal
# test on mean difference between partitions
# normality of the distribution
shapiro.test(resdata.full$Sc) # significant => data not normal
kruskal.test(Sc ~ partition, resdata.full)
wilcox.Sc.part <- pairwise.wilcox.test(resdata.full$Sc, resdata.full$partition, p.adjust.method = "bonferroni")
# post hoc display
cld.Sc.mean.part <- as.data.frame.list(multcompLetters(fullPTable(wilcox.Sc.part$p.value), compare="<", threshold=0.05,Letters= LETTERS, reversed = FALSE))

### Cnidocyte count (Cc) 
# test on variance difference between partitions
leveneTest(Cc ~ partition, resdata.full) # significant => variance not equal
# test on mean difference between partitions
# normality of the distribution
shapiro.test(resdata.full$Cc) # significant => data not normal
kruskal.test(Cc ~ partition, resdata.full)
wilcox.Cc.part <- pairwise.wilcox.test(resdata.full$Cc, resdata.full$partition, p.adjust.method = "bonferroni")
# post hoc display
cld.Cc.mean.part <- as.data.frame.list(multcompLetters(fullPTable(wilcox.Cc.part$p.value), compare="<", threshold=0.05,Letters= LETTERS, reversed = FALSE))

### chlorophyll a (chla) 
# test on variance difference between partitions
leveneTest(chla ~ partition, resdata.full) # non significant => variance equal 
# test on mean difference between partitions
# normality of the distribution
shapiro.test(resdata.full$chla) # non significant => data normal, so technically we could directly compare species using anova
anova.chla.part <- aov(resdata.full$chla ~ partition, resdata.full)
summary(anova.chla.part)
tukey.chla.part <- TukeyHSD(anova.chla.part)
# post hoc display 
cld.chlameanpart <- multcompLetters4(anova.chla.part, tukey.chla.part, Letters= LETTERS)
cld.chlameanpart.display <- as.data.frame.list(cld.chlameanpart$partition)
t4 <- cld.chlameanpart.display[c('1', '2', '3'),]
cld.chlameanpart.display <- as.data.frame.list(t4)
Trait_median_partition$cld.chlameanpart <- cld.chlameanpart.display$Letters
Trait_median_partition$cld.chlameanpart.up <- lapply(Trait_median_partition$cld.chlameanpart, toupper)

### chlorophyll b (Chlb) 
# test on variance difference between partitions
leveneTest(chlb ~ partition, resdata.full) # significant => variance not equal
# test on mean difference between partitions
# normality of the distribution
shapiro.test(resdata.full$chlb) # significant => data not normal
kruskal.test(chlb ~ partition, resdata.full)
wilcox.chlb.part <- pairwise.wilcox.test(resdata.full$chlb, resdata.full$partition, p.adjust.method = "bonferroni")
# post hoc display
cld.chlb.mean.part <- as.data.frame.list(multcompLetters(fullPTable(wilcox.chlb.part$p.value), compare="<", threshold=0.05,Letters= LETTERS, reversed = FALSE))

### chlorophyll c (chlc)
# test on variance difference between partitions
leveneTest(chlc ~ resdata.full$partition, resdata.full) # significant => variance not equal
# test on mean difference between partitions
# normality of the distribution
shapiro.test(resdata.full$chlc) # significant => data not normal
kruskal.test(chlc ~ partition, resdata.full)
wilcox.chlc.part <- pairwise.wilcox.test(resdata.full$chlc, resdata.full$partition, p.adjust.method = "bonferroni")
# post hoc display
cld.chlc.mean.part <- as.data.frame.list(multcompLetters(fullPTable(wilcox.chlc.part$p.value), compare="<", threshold=0.05,Letters= LETTERS, reversed = FALSE))

### chlorophyll ratio c/a (chlca) 
# test on variance difference between partitions
leveneTest(chlca ~ partition, resdata.full) # significant => variance non equal 
# test on mean difference between partitions
# normality of the distribution
shapiro.test(resdata.full$chlca) # non significant => data normal
kruskal.test(chlca ~ partition, resdata.full)
wilcox.chlca.part <- pairwise.wilcox.test(resdata.full$chlca, resdata.full$partition, p.adjust.method = "bonferroni")
# post hoc display
cld.chlca.mean.part <- as.data.frame.list(multcompLetters(fullPTable(wilcox.chlca.part$p.value), compare="<", threshold=0.05,Letters= LETTERS, reversed = FALSE))

### Symbiodiniaceae protein (Spro) 
# test on variance difference between partition
leveneTest(Spro ~ partition, resdata.full) # significant => variance not equal 
# test on mean difference between partition
# normality of the distribution
shapiro.test(resdata.full$Spro) # significant => data not normal
kruskal.test(Spro ~ partition, resdata.full)
wilcox.Spro.part <- pairwise.wilcox.test(resdata.full$Spro, resdata.full$partition, p.adjust.method = "bonferroni")
#post hoc display
cld.Spro.mean.part <- as.data.frame.list(multcompLetters(fullPTable(wilcox.Spro.part$p.value), compare="<", threshold=0.05,Letters= LETTERS, reversed = FALSE))

### Host protein (Hpro) 
# test on variance difference between partition
leveneTest(Hpro ~ partition, resdata.full) # non significant => variance are equal 
# test on mean difference between partition
# normality of the distribution
shapiro.test(resdata.full$Hpro) # non significant => data normal, so technically we could directly compare species using anova
anova.Hpro.part <- aov(resdata.full$Hpro ~ partition, resdata.full)
summary(anova.Hpro.part)
tukey.Hpro.part <- TukeyHSD(anova.Hpro.part)
#post hoc display
cld.Hpromeanpart <- multcompLetters4(anova.Hpro.part, tukey.Hpro.part, Letters= LETTERS)
cld.Hpromeanpart.display <- as.data.frame.list(cld.Hpromeanpart$partition)
t9 <- cld.Hpromeanpart.display[c('1', '2', '3'),]
cld.Hpromeanpart.display <- as.data.frame.list(t9)
Trait_median_partition$cld.Hpromeanpart <- cld.chlameanpart.display$Letters
Trait_median_partition$cld.Hpromeanpart.up <- lapply(Trait_median_partition$cld.Hpromeanpart, toupper)

### H13C
# test on variance difference between partition
leveneTest(H13C ~ partition, resdata.full) # significant => variance are not equal 
# test on mean difference between partition
# normality of the distribution
shapiro.test(resdata.full$H13C) # significant => data not normal
kruskal.test(H13C ~ partition, resdata.full)
wilcox.H13C.part <- pairwise.wilcox.test(resdata.full$H13C, resdata.full$partition, p.adjust.method = "bonferroni")
#post hoc display
cld.H13C.mean.part <- as.data.frame.list(multcompLetters(fullPTable(wilcox.H13C.part$p.value), compare="<", threshold=0.05,Letters= LETTERS, reversed = FALSE))

### H15N
# test on variance difference between partition
leveneTest(H15N ~ partition, resdata.full) # significant => variance are not equal 
# test on mean difference between partition
# normality of the distribution
shapiro.test(resdata.full$H15N) # significant => data not normal
kruskal.test(H15N ~ partition, resdata.full)
wilcox.H15N.part <- pairwise.wilcox.test(resdata.full$H15N, resdata.full$partition, p.adjust.method = "bonferroni")
#post hoc display
cld.H15N.mean.part <- as.data.frame.list(multcompLetters(fullPTable(wilcox.H15N.part$p.value), compare="<", threshold=0.05,Letters= LETTERS, reversed = FALSE))

### S13C 
# test on variance difference between partition
leveneTest(S13C ~ partition, resdata.full) # significant => variance are not equal 
# test on mean difference between partition
# normality of the distribution
shapiro.test(resdata.full$S13C) # significant => data not normal
kruskal.test(S13C ~ partition, resdata.full)
wilcox.S13C.part <- pairwise.wilcox.test(resdata.full$S13C, resdata.full$partition, p.adjust.method = "bonferroni")
#post hoc display
cld.S13C.mean.part <- as.data.frame.list(multcompLetters(fullPTable(wilcox.S13C.part$p.value), compare="<", threshold=0.05,Letters= LETTERS, reversed = FALSE))

### S15N
# test on variance difference between partition
leveneTest(S15N ~ partition, resdata.full) # non significant => variance are equal 
# test on mean difference between partition
# normality of the distribution
shapiro.test(resdata.full$S15N) # significant => data not normal
kruskal.test(S15N ~ partition, resdata.full)
wilcox.S15N.part <- pairwise.wilcox.test(resdata.full$S15N, resdata.full$partition, p.adjust.method = "bonferroni")
#post hoc display
cld.S15N.mean.part <- as.data.frame.list(multcompLetters(fullPTable(wilcox.S13C.part$p.value), compare="<", threshold=0.05,Letters= LETTERS, reversed = FALSE))

### HCN
# test on variance difference between partition
leveneTest(HCN ~ partition, resdata.full) #  significant => variance are not equal 
# test on mean difference between partition
# normality of the distribution
kruskal.test(HCN ~ partition, resdata.full)
wilcox.HCN.part <- pairwise.wilcox.test(resdata.full$HCN, resdata.full$partition, p.adjust.method = "bonferroni")
#post hoc display
cld.HCN.mean.part <- as.data.frame.list(multcompLetters(fullPTable(wilcox.HCN.part$p.value), compare="<", threshold=0.05,Letters= LETTERS, reversed = FALSE))

### SCN 
# test on variance difference between species
leveneTest(SCN ~ partition, resdata.full) # significant => variance not equal 
# test on mean difference between partition
# normality of the distribution
shapiro.test(resdata.full$SCN) # significant => data not normal
kruskal.test(SCN ~ partition, resdata.full)
wilcox.SCN.part <- pairwise.wilcox.test(resdata.full$SCN, resdata.full$partition, p.adjust.method = "bonferroni")
#post hoc display
cld.SCN.mean.part <- as.data.frame.list(multcompLetters(fullPTable(wilcox.SCN.part$p.value), compare="<", threshold=0.05,Letters= LETTERS, reversed = FALSE))

#' 
#' ### **Mantel tests**
## ----eval = T, echo = T, message = F, warning = F, results = "hide"-----------
### Mantel tests (overall)
matrice.tr = daisy(resdata.fully, metric = c("euclidean"), stand = F)
matrice.env = daisy(envdata.fully, metric = c("euclidean"), stand = F)
test_env = mantel(matrice.tr, matrice.env, method = "pearson", permutations = 999, na.rm = TRUE)

### Mantel tests (partitions)
#Partition 1
tradata.full.partition1 = subset(resdata.full.partition1, select = c(dens, Sc, Cc, chla, chlb, chlc, chlca, Spro, Hpro, H13C,H15N, S13C, S15N, HCN, SCN))
envdata.full.partition1 = subset(resdata.full.partition1, select = c(SST, SST_sd, light, light_sd))
test_env_part1 = mantel(dist(tradata.full.partition1), dist(envdata.full.partition1), method = "pearson", permutations = 999, na.rm = TRUE)

#Partition 2
tradata.full.partition2 = subset(resdata.full.partition1, select = c(dens, Sc, Cc, chla, chlb, chlc, chlca, Spro, Hpro, H13C,H15N, S13C, S15N, HCN, SCN))
envdata.full.partition2 = subset(resdata.full.partition1, select = c(SST, SST_sd, light, light_sd))
test_env_part2 = mantel(dist(tradata.full.partition2), dist(envdata.full.partition2), method = "pearson", permutations = 999, na.rm = TRUE)

#Partition 3
tradata.full.partition3 = subset(resdata.full.partition3, select = c(dens, Sc, Cc, chla, chlb, chlc, chlca, Spro, Hpro, H13C,H15N, S13C, S15N, HCN, SCN))
envdata.full.partition3 = subset(resdata.full.partition3, select = c(SST, SST_sd, light, light_sd))
test_env_part3 = mantel(dist(tradata.full.partition3), dist(envdata.full.partition3), method = "pearson", permutations = 999, na.rm = TRUE)

## Mantel tests (per species)
#A. muricata
tradata.fully.amur = subset(resdata.full.amur, select = c(dens, Sc, Cc, chla, chlb, chlc, chlca, Spro, Hpro, H13C,H15N, S13C, S15N, HCN, SCN))
envdata.full.amur = subset(resdata.full.amur, select = c(SST, SST_sd, light, light_sd))
test_env_amur = mantel(dist(tradata.fully.amur), dist(envdata.full.amur), method = "pearson", permutations = 999, na.rm = TRUE)

#C. microphthalma
tradata.fully.cmic = subset(resdata.full.cmic, select = c(dens, Sc, Cc, chla, chlb, chlc, chlca, Spro, Hpro, H13C,H15N, S13C, S15N, HCN, SCN))
envdata.fully.cmic = subset(resdata.full.cmic, select = c(SST, SST_sd, light, light_sd))
test_env_cmic = mantel(dist(tradata.fully.cmic), dist(envdata.fully.cmic), method = "pearson", permutations = 999, na.rm = TRUE)

#I. palifera
tradata.fully.ipal = subset(resdata.full.ipal, select = c(dens, Sc, Cc, chla, chlb, chlc, chlca, Spro, Hpro, H13C,H15N, S13C, S15N, HCN, SCN))
envdata1.fully.ipal = subset(resdata.full.ipal, select = c(SST, SST_sd, light, light_sd))
test_env_ipal = mantel(dist(tradata.fully.ipal), dist(envdata1.fully.ipal), method = "pearson", permutations = 999, na.rm = TRUE)

#P. lutea
tradata.fully.plut = subset(resdata.full.plut, select = c(dens, Sc, Cc, chla, chlb, chlc, chlca, Spro, Hpro, H13C,H15N, S13C, S15N, HCN, SCN))
envdata1.fully.plut = subset(resdata.full.plut, select = c(SST, SST_sd, light, light_sd))
test_env_plut = mantel(dist(tradata.fully.plut), dist(envdata1.fully.plut), method = "pearson", permutations = 999, na.rm = TRUE)

#P. profundacella
tradata.fully.ppro = subset(resdata.full.ppro, select = c(dens, Sc, Cc, chla, chlb, chlc, chlca, Spro, Hpro, H13C,H15N, S13C, S15N, HCN, SCN))
envdata1.fully.ppro = subset(resdata.full.ppro, select = c(SST, SST_sd, light, light_sd))
test_env_ppro = mantel(dist(tradata.fully.ppro), dist(envdata1.fully.ppro), method = "pearson", permutations = 999, na.rm = TRUE)
count(resdata.full.spis)

#P. speciosa
tradata.fully.pspe = subset(resdata.full.pspe, select = c(dens, Sc, Cc, chla, chlb, chlc, chlca, Spro, Hpro, H13C,H15N, S13C, S15N, HCN, SCN))
envdata1.fully.pspe = subset(resdata.full.pspe, select = c(SST, SST_sd, light, light_sd))
test_env_pspe = mantel(dist(tradata.fully.pspe), dist(envdata1.fully.pspe), method = "pearson", permutations = 999, na.rm = TRUE)

#S. pistillata
tradata.fully.spis = subset(resdata.full.spis, select = c(dens, Sc, Cc, chla, chlb, chlc, chlca, Spro, Hpro, H13C,H15N, S13C, S15N, HCN, SCN))
envdata1.fully.spis = subset(resdata.full.spis, select = c(SST, SST_sd, light, light_sd))
test_env_spis = mantel(dist(tradata.fully.spis), dist(envdata1.fully.spis), method = "pearson", permutations = 999, na.rm = TRUE)

#' 
#' ### **GLMM and GLM**
## ----eval = T, echo = T, message = F, warning = F, results = "hide"-----------
# standardization for trait data
# import data

write.csv(resdata.full,'Data/DePalmas_et_al_CR_res.csv')
resdata <- read.csv ("Data/DePalmas_et_al_CR_res.csv", header = TRUE, sep=",", dec=".")

# standardization for trait data
trait_data <- resdata[,c(3:17)]
trait_data_trans_raw <- sapply(trait_data, scale) 
trait_data_trans_raw <- as.data.frame(trait_data_trans_raw)
trait_data_trans_raw <- as.data.frame(cbind(trait_data_trans_raw, resdata[,c(2,18:21)]))
trait_data_trans <- trait_data_trans_raw 

# visually check the normality of trait data
hist(trait_data_trans$dens)
hist(trait_data_trans$Sc)
hist(trait_data_trans$Cc)
hist(log(trait_data_trans$Cc+1)) # Cc is log(x+1) transformed 
hist(trait_data_trans$chla)
hist(trait_data_trans$chlb)
hist(log(trait_data_trans$chlb+1)) # Chlb is log(x+1) transformed 
hist(trait_data_trans$chlc)
hist(trait_data_trans$chlca)
hist(trait_data_trans$Spro)
hist(trait_data_trans$Hpro)
hist(trait_data_trans$H13C)
hist(trait_data_trans$H15N)
hist(trait_data_trans$S13C)
hist(trait_data_trans$S15N)
hist(trait_data_trans$HCN)
hist(log(trait_data_trans$HCN+1)) # HCN is log(x+1) transformed 
hist(trait_data_trans$SCN)

trait_data_trans <-trait_data_trans_raw
CcLog <- log(trait_data_trans$Cc+1) # log(x+1) transformed Cc
chlbLog <- log(trait_data_trans$chlb+1) # log(x+1) transformed chlb
HCNLog <- log(trait_data_trans$HCN+1) # log(x+1) transformed HCN
trait_data_trans$Cc <- NULL
trait_data_trans$chlb <- NULL
trait_data_trans$HCN <- NULL
trait_data_trans_log <- data.frame(CcLog, chlbLog, HCNLog, trait_data_trans )
traits_list_log <- colnames(trait_data_trans_log)[1:15]

# visually check the normality of environmental data
hist(trait_data_trans$SST)
hist(trait_data_trans$SST_sd)
hist(trait_data_trans$light)
hist(trait_data_trans$light_sd)

# linear mixed effect models were built for each trait, with environmental factors 
# as predictors, species as a random factor, and trait as responsive factors (15 LMMs)
fix_est_mat <- as.data.frame(matrix(NA, nrow = length(traits_list_log)*4, 6))
for (i in 1:length(traits_list_log)){
  trait <- trait_data_trans_log[, colnames(trait_data_trans_log) == traits_list_log[i]]
  lmm_all <- lmer(trait ~ SST + SST_sd + light + light_sd + (1 | spe), data = trait_data_trans_log)
  lmm_all_sum <- summary(lmm_all) 
  
  fix_est_mat[(4*i-3),] <- c(rownames(lmm_all_sum$coefficients)[2], as.numeric(round(lmm_all_sum$coefficients[2,1],2)), 
                             round((lmm_all_sum$coefficients[2,1] + lmm_all_sum$coefficients[2,2]),2), 
                             round((lmm_all_sum$coefficients[2,1] - lmm_all_sum$coefficients[2,2]),2),
                             round(lmm_all_sum$coefficients[2,5],2), ifelse(round(lmm_all_sum$coefficients[2,5],2)<0.05, '*','') )
  fix_est_mat[(4*i-2),] <- c(rownames(lmm_all_sum$coefficients)[3], as.numeric(round(lmm_all_sum$coefficients[3,1],2)),
                             round((lmm_all_sum$coefficients[3,1] + lmm_all_sum$coefficients[3,2]),2), 
                             round((lmm_all_sum$coefficients[3,1] - lmm_all_sum$coefficients[3,2]),2),
                             round(lmm_all_sum$coefficients[3,5],2), ifelse(round(lmm_all_sum$coefficients[3,5],2)<0.05, '*',''))
  fix_est_mat[(4*i-1),] <- c(rownames(lmm_all_sum$coefficients)[4], as.numeric(round(lmm_all_sum$coefficients[4,1],3)), 
                             round((lmm_all_sum$coefficients[4,1] + lmm_all_sum$coefficients[4,2]),3), 
                             round((lmm_all_sum$coefficients[4,1] - lmm_all_sum$coefficients[4,2]),3),
                             round(lmm_all_sum$coefficients[4,5],2), ifelse(round(lmm_all_sum$coefficients[4,5],2)<0.05, '*',''))
  fix_est_mat[(4*i),] <- c(rownames(lmm_all_sum$coefficients)[5], as.numeric(round(lmm_all_sum$coefficients[5,1],2)), 
                           round((lmm_all_sum$coefficients[5,1] + lmm_all_sum$coefficients[5,2]),2), 
                           round((lmm_all_sum$coefficients[5,1] - lmm_all_sum$coefficients[5,2]),2),
                           round(lmm_all_sum$coefficients[5,5],2), ifelse(round(lmm_all_sum$coefficients[5,5],2)<0.05, '*',''))
}

trait <- c(rep('dens',4),rep('sc',4),rep('cc',4),rep('chla',4),rep('chlb',4),
           rep('chlc',4),rep('chlca',4),rep('spro',4),rep('hpro',4),rep('h13c',4),
           rep('h15n',4),rep('s13c',4),rep('s15n',4),rep('HCN',4),rep('SCN',4))
fix_est_mat <- cbind(trait, fix_est_mat)
colnames(fix_est_mat) <- c('trait','environment', 'estimate','upperestimate', 'lowerestimate','pvalue','sig')
head(fix_est_mat)

fix_est_mat_wide <- fix_est_mat%>%
  select(trait, environment, estimate)%>%
  spread(environment, estimate)
write.csv(fix_est_mat_wide, 'Results/15LMM.csv')

fix_est_mat_sig_wide <- fix_est_mat%>%
  select(trait, environment, sig)%>%
  spread(environment, sig)
write.csv(fix_est_mat_sig_wide, 'Results/15LMM.sig.csv')


# linear models were built for each trait and species, with environmental factors 
# as predictors, and trait as responsive factors (105 LMs)
species_list <- unique(trait_data_trans_raw$spe)
traits_list <- colnames(trait_data_trans_raw)[1:15]
species_lm_mat <- as.data.frame(matrix(NA, nrow = 4 * 7 , ncol = 8))
species_tra_lm_raw <- as.data.frame(matrix(NA, nrow = 4 * 7 * 15 , ncol = 8))
species_tra_lm_mat <- as.data.frame(matrix(NA, nrow = 4 * 7 * 15 , ncol = 8))

for (k in 1:length(traits_list)){
  for (i in 1:length(species_list)){
    tra <- trait_data_trans_raw[, colnames(trait_data_trans_raw) == traits_list[k]]
    SST <- trait_data_trans_raw$SST
    SST_sd <- trait_data_trans_raw$SST_sd
    light <- trait_data_trans_raw$light
    light_sd <- trait_data_trans_raw$light_sd
    spe <- trait_data_trans_raw$spe
    trait_spe_mat <- data.frame(tra, spe, SST, SST_sd, light, light_sd)
    
    species_mat <- subset(trait_spe_mat, spe == species_list[i])
    species_lm <- glm(tra-1 ~ SST + SST_sd + light + light_sd, species_mat, family = gaussian)
    species_lm_sum <- summary(species_lm)
    
    if (nrow(species_lm_sum$coefficients) == 5){
      var <- rownames(species_lm_sum$coefficients)[2:5]
      species <- rep(species_list[i],4)
      trait <- rep(traits_list[k],4)
      sig <- ifelse(species_lm_sum$coefficients[2:5,4]<0.05,'*','')
      species_lm_mat[(4*i-3):(4*i),] <- data.frame(trait, species,var, round(species_lm_sum$coefficients[2:5,],2),sig)
      
    }else {
      var <- rownames(species_lm_sum$coefficients)[2:4]
      species <- rep(species_list[i],3)
      trait <- rep(traits_list[k],3)
      sig <- ifelse(species_lm_sum$coefficients[2:4,4]<0.05,'*','')
      species_lm_mat[(4*i-3):(4*i-1),] <- data.frame(trait, species,var, round(species_lm_sum$coefficients[2:4,],2),sig)
    }
    
  }
  species_tra_lm_raw[(28*k-27):(28*k),] <- species_lm_mat
  colnames(species_tra_lm_raw) <- c('trait','species','environment','estimate','Std. Error','t-value','p-value','significance')
  species_tra_lm_raw
  species_tra_lm_mat <- na.omit(species_tra_lm_raw)
  
}

head(species_tra_lm_mat)

species_est <-  species_tra_lm_mat %>%
  select(species, trait, environment, estimate) %>%
  spread(environment, estimate)
write.csv(species_est, 'Results/105LM.est.csv')

species_sig <-  species_tra_lm_mat %>%
  select(species, trait, environment, significance) %>%
  spread(environment, significance)
write.csv(species_sig, 'Results/105LM.sig.csv')


#' 
#' ## **Figure preparation**
## ----eval = T, echo = T, message = F, warning = F, results = "hide"-----------
### skeleton density 
density.plot.part = resdata.full[,c('partition','dens')]
S = ggplot(data = density.plot.part, aes(x = partition, y = dens, fill = partition))
density_main_title = expression(paste("Skeleton density"))
density_title = expression(paste( g.cm ^ ' -3')) 
p1<-S + geom_flat_violin(position = position_nudge(x = .2, y = 0),adjust =2)+
  geom_jitter(position = position_jitter(width = .15), shape =16)+
  geom_boxplot(width = .1, alpha = 0.5)+
  labs(title = density_main_title, y = density_title)+
  theme_bw()+
  theme(legend.position = 'none')+
  theme(axis.text.x= element_blank())+
  theme(plot.title = element_text(hjust = "0.5"))+
  theme(axis.title.x = element_blank())+
  scale_y_continuous(limits=c(0.0000, 0.0040))+
  stat_summary(geom = 'text', label = cld.dens.mean.part$Letters, fun = min, vjust = 1.4, color = "black")+ 
  coord_cartesian(expand=T)

### Symbiodiniaceae content (Sc) 
symbio.plot.part = resdata.full[,c('partition','Sc')]
S = ggplot(data = symbio.plot.part, aes(x=partition, y= Sc, fill=partition))
symbio_main_title = expression(paste("Symbiodiniaceae content"))
symbio_title = expression(paste(cell.mgAFDW^ ' -1')) 
p2<-S + geom_flat_violin(position = position_nudge(x = .2, y = 0),adjust =2)+
  geom_jitter(position = position_jitter(width = .15), shape =16)+
  geom_boxplot(width = .1, alpha = 0.5)+
  labs(title= symbio_main_title, y=symbio_title)+
  theme_bw()+
  theme(legend.position='none')+
  theme(axis.text.x= element_blank())+
  theme(plot.title= element_text(hjust="0.5"))+
  theme(axis.title.x = element_blank())+
  scale_y_continuous(limits=c(-180000, 1600000))+
  stat_summary(geom = 'text', label = cld.Sc.mean.part$Letters, fun = min, vjust = 1.4, color = "black")+ 
  coord_cartesian(expand=T)

### Cnidocyte content (Cc) 
cnido.plot.part = resdata.full[,c('partition','Cc')]
S = ggplot(data = cnido.plot.part, aes(x=partition, y= Cc, fill=partition))
cnido_main_title = expression(paste("Cnidocyte content"))
cnido_title = expression(paste(cnidocyte.mgAFDW^ ' -1')) 
p3<-S + geom_flat_violin(position = position_nudge(x = .2, y = 0),adjust =2)+
  geom_jitter(position = position_jitter(width = .15), shape =16)+
  geom_boxplot(width = .1, alpha = 0.5)+
  labs(title= cnido_main_title, y=cnido_title)+
  theme_bw()+
  theme(legend.position='none')+
  theme(axis.text.x= element_blank())+
  theme(plot.title= element_text(hjust="0.5"))+
  theme(axis.title.x = element_blank())+
  scale_y_continuous(limits=c(-100000, 760000))+
  stat_summary(geom = 'text', label = cld.Cc.mean.part$Letters, fun = min, vjust = 1.4, color = "black")+ 
  coord_cartesian(expand=T)

### Chlorophyll a (chl a) 
chla.plot.part = resdata.full[,c('partition','chla')]
S = ggplot(data = chla.plot.part, aes(x=partition, y= chla, fill=partition))
chla_main_title = expression(paste("Chlorophyll ", italic ("a")))
chla_title = expression(paste(μg.mgAFDW^ ' -1')) 
p4<-S + geom_flat_violin(position = position_nudge(x = .2, y = 0),adjust =2)+
  geom_jitter(position = position_jitter(width = .15), shape =16)+
  geom_boxplot(width = .1, alpha = 0.5)+
  labs(title= chla_main_title, y=chla_title)+
  theme_bw()+
  theme(legend.position='none')+
  theme(axis.text.x= element_blank())+
  theme(plot.title= element_text(hjust="0.5"))+
  theme(axis.title.x = element_blank())+
  scale_y_continuous(limits=c(-0.5, 8))+
  stat_summary(geom = 'text', label = cld.chlameanpart.display$Letters, fun = min, vjust = 1.4, color = "black")+ 
  coord_cartesian(expand=T)

### Chlorophyll b (chl b) 
chlb.plot.part = resdata.full[,c('partition','chlb')]
S = ggplot(data = chlb.plot.part, aes(x=partition, y= chlb, fill=partition))
chlb_main_title = expression(paste("Chlorophyll ", italic ("b")))
chlb_title = expression(paste(μg.mgAFDW^ ' -1')) 
p5<-S + geom_flat_violin(position = position_nudge(x = .2, y = 0),adjust =2)+
  geom_jitter(position = position_jitter(width = .15), shape =16)+
  geom_boxplot(width = .1, alpha = 0.5)+
  labs(title= chlb_main_title, y=chlb_title)+
  theme_bw()+
  theme(legend.position='none')+
  theme(axis.text.x= element_blank())+
  theme(plot.title= element_text(hjust="0.5"))+
  theme(axis.title.x = element_blank())+
  scale_y_continuous(limits=c(-0.5, 2.9))+
  stat_summary(geom = 'text', label = cld.chlb.mean.part$Letters, fun = min, vjust = 1.4, color = "black")+ 
  coord_cartesian(expand=T)

### Chlorophyll c (chl c) 
chlc.plot.part = resdata.full[,c('partition','chlc')]
S = ggplot(data = chlc.plot.part, aes(x=partition, y= chlc, fill=partition))
chlc_main_title = expression(paste("Chlorophyll ", italic ("c")))
chlc_title = expression(paste(μg.mgAFDW^ ' -1')) 
p6<-S + geom_flat_violin(position = position_nudge(x = .2, y = 0),adjust =2)+
  geom_jitter(position = position_jitter(width = .15), shape =16)+
  geom_boxplot(width = .1, alpha = 0.5)+
  labs(title= chlc_main_title, y=chlc_title)+
  theme_bw()+
  theme(legend.position='none')+
  theme(axis.text.x= element_blank())+
  theme(plot.title= element_text(hjust="0.5"))+
  theme(axis.title.x = element_blank())+
  scale_y_continuous(limits=c(-0.5, 2.6))+
  stat_summary(geom = 'text', label = cld.chlc.mean.part$Letters, fun = min, vjust = 1.4, color = "black")+ 
  coord_cartesian(expand=T)

### ratio Chlorophyll c/a 
ratio.plot.part = resdata.full[,c('partition','chlca')]
S = ggplot(data = ratio.plot.part, aes(x=partition, y= chlca, fill=partition))
ratio_main_title = expression(paste("Chlorophyll ", italic ("c/a"), " ratio"))
ratio_title = expression(paste("ratio chlorophyll ",italic ("c/a"))) 
p7<-S + geom_flat_violin(position = position_nudge(x = .2, y = 0),adjust =2)+
  geom_jitter(position = position_jitter(width = .15), shape =16)+
  geom_boxplot(width = .1, alpha = 0.5)+
  labs(title= ratio_main_title, y=ratio_title)+
  theme_bw()+
  theme(legend.position='none')+
  theme(axis.text.x= element_blank())+
  theme(plot.title= element_text(hjust="0.5"))+
  theme(axis.title.x = element_blank())+
  scale_y_continuous(limits=c(-0.05, 0.7))+
  stat_summary(geom = 'text', label = cld.chlca.mean.part$Letters, fun = min, vjust = 1.4, color = "black")+ 
  coord_cartesian(expand=T)

### Symbiodiniaceae protein (Spro)
prozoox.plot.part = resdata.full[,c('partition','Spro')]
S = ggplot(data = prozoox.plot.part, aes(x=partition, y= Spro, fill=partition))
prozoox_main_title = expression(paste("Total Symbiodiniaceae protein"))
prozoox_title = expression(paste(μg.mgAFDW^ ' -1')) 
p8<-S + geom_flat_violin(position = position_nudge(x = .2, y = 0),adjust =2)+
  geom_jitter(position = position_jitter(width = .15), shape =16)+
  geom_boxplot(width = .1, alpha = 0.5)+
  labs(title= prozoox_main_title, y=prozoox_title)+
  theme_bw()+
  theme(legend.position='none')+
  theme(axis.text.x= element_blank())+
  theme(plot.title= element_text(hjust="0.5"))+
  theme(axis.title.x = element_blank())+
  scale_y_continuous(limits=c(-50, 250))+
  stat_summary(geom = 'text', label = cld.Spro.mean.part$Letters, fun = min, vjust = 1.4, color = "black")+ 
  coord_cartesian(expand=T)

### Coral protein (Hpro)
prohost.plot.part = resdata.full[,c('partition','Hpro')]
S = ggplot(data = prohost.plot.part, aes(x=partition, y= Hpro, fill=partition))
prohost_main_title = expression(paste("Total host protein"))
prohost_title = expression(paste(μg.mgAFDW^ ' -1')) 
p9<-S + geom_flat_violin(position = position_nudge(x = .2, y = 0),adjust =2)+
  geom_jitter(position = position_jitter(width = .15), shape =16)+
  geom_boxplot(width = .1, alpha = 0.5)+
  labs(title= prohost_main_title, y=prohost_title)+
  theme_bw()+
  theme(legend.position='none')+
  theme(axis.text.x= element_blank())+
  theme(plot.title= element_text(hjust="0.5"))+
  theme(axis.title.x = element_blank())+
  scale_y_continuous(limits=c(0, 390))+
  stat_summary(geom = 'text', label = Trait_median_partition$cld.Hpromeanpart.up, fun = min, vjust = 1.4, color = "black")+ 
  coord_cartesian(expand=T)

### H13C
Hdelta13C.plot.part = resdata.full[,c('partition','H13C')]
S = ggplot(data = Hdelta13C.plot.part, aes(x=partition, y= H13C, fill=partition))
delta13C_main_title = expression(paste("Host ", delta^13, " C"))
delta13C_title = expression(paste(delta^13, " C (\u2030)")) 
p10<-S + geom_flat_violin(position = position_nudge(x = .2, y = 0),adjust =2)+
  geom_jitter(position = position_jitter(width = .15), shape =16)+
  geom_boxplot(width = .1, alpha = 0.5)+
  labs(title= delta13C_main_title, y=delta13C_title)+
  theme_bw()+
  theme(legend.position='none')+
  theme(axis.text.x= element_blank())+
  theme(plot.title= element_text(hjust="0.5"))+
  theme(axis.title.x = element_blank())+
  scale_y_continuous(limits=c(-26, -10))+
  stat_summary(geom = 'text', label = cld.H13C.mean.part$Letters, fun = min, vjust = 1.4, color = "black")+ 
  coord_cartesian(expand=T)

### H15N
Hdelta15N.plot.part = resdata.full[,c('partition','H15N')]
S = ggplot(data = Hdelta15N.plot.part, aes(x=partition, y= H15N, fill=partition))
delta15N_main_title = expression(paste("Host ", delta^15, " N"))
delta15N_title = expression(paste(delta^15, " N (\u2030)")) 
p11<-S + geom_flat_violin(position = position_nudge(x = .2, y = 0),adjust =2)+
  geom_jitter(position = position_jitter(width = .15), shape =16)+
  geom_boxplot(width = .1, alpha = 0.5)+
  labs(title= delta15N_main_title, y=delta15N_title)+
  theme_bw()+
  theme(legend.position='none')+
  theme(axis.text.x= element_blank())+
  theme(plot.title= element_text(hjust="0.5"))+
  theme(axis.title.x = element_blank())+
  scale_y_continuous(limits=c(-8.5, 9.5))+
  stat_summary(geom = 'text', label = cld.H15N.mean.part$Letters, fun = min, vjust = 1.4, color = "black")+ 
  coord_cartesian(expand=T)

### S13C
Zdelta13C.plot.part = resdata.full[,c('partition','S13C')]
S = ggplot(data = Zdelta13C.plot.part, aes(x=partition, y= S13C, fill=partition))
Zdelta13C_main_title = expression(paste("Symbiodiniaceae ", delta^13, " C"))
Zdelta13C_title = expression(paste(delta^13, " C (\u2030)")) 
p12<-S + geom_flat_violin(position = position_nudge(x = .2, y = 0),adjust =2)+
  geom_jitter(position = position_jitter(width = .15), shape =16)+
  geom_boxplot(width = .1, alpha = 0.5)+
  labs(title= Zdelta13C_main_title, y=Zdelta13C_title)+
  theme_bw()+
  theme(legend.position='none')+
  theme(axis.text.x= element_blank())+
  theme(plot.title= element_text(hjust="0.5"))+
  theme(axis.title.x = element_blank())+
  scale_y_continuous(limits=c(-24, -10))+
  stat_summary(geom = 'text', label = cld.S13C.mean.part$Letters, fun = min, vjust = 1.4, color = "black")+ 
  coord_cartesian(expand=T)

### S15N
Zdelta15N.plot.part = resdata.full[,c('partition','S15N')]
S = ggplot(data = Zdelta15N.plot.part, aes(x=partition, y= S15N, fill=partition))
Zdelta15N_main_title = expression(paste("Symbiodiniaceae ", delta^15, " N"))
Zdelta15N_title = expression(paste(delta^15, " N (\u2030)")) 
p13<-S + geom_flat_violin(position = position_nudge(x = .2, y = 0),adjust =2)+
  geom_jitter(position = position_jitter(width = .15), shape =16)+
  geom_boxplot(width = .1, alpha = 0.5)+
  labs(title= Zdelta15N_main_title, y=Zdelta15N_title)+
  theme_bw()+
  theme(legend.position='none')+
  theme(axis.text.x= element_blank())+
  theme(plot.title= element_text(hjust="0.5"))+
  theme(axis.title.x = element_blank())+
  scale_y_continuous(limits=c(0.6, 6))+
  stat_summary(geom = 'text', label = cld.S15N.mean.part$Letters, fun = min, vjust = 1.4, color = "black")+ 
  coord_cartesian(expand=T)

### HCN
H.CN.plot.part = resdata.full[,c('partition','HCN')]
S = ggplot(data = H.CN.plot.part, aes(x=partition, y=HCN, fill=partition))
H.CN_main_title = expression(paste("Host C:N ratio"))
H.CN_y_title = expression(paste("Ratio")) 
p14<-S + geom_flat_violin(position = position_nudge(x = .2, y = 0),adjust =2)+
  geom_jitter(position = position_jitter(width = .15), shape =16)+
  geom_boxplot(width = .1, alpha = 0.5)+
  labs(title= H.CN_main_title, y=H.CN_y_title)+
  theme_bw()+
  theme(legend.position='none')+
  theme(axis.text.x= element_blank())+
  theme(plot.title= element_text(hjust="0.5"))+
  theme(axis.title.x = element_blank())+
  scale_y_continuous(limits=c(2, 18))+
  stat_summary(geom = 'text', label = cld.HCN.mean.part$Letters, fun = min, vjust = 1.4, color = "black")+ 
  coord_cartesian(expand=T)

### SCN
Z.CN.plot.part = resdata.full[,c('partition','SCN')]
S = ggplot(data = Z.CN.plot.part, aes(x= partition, y=SCN, fill=partition))
Z.CN_main_title = expression(paste("Symbiodiniaceae C:N ratio"))
Z.CN_y_title = expression(paste("Ratio")) 
p15<-S + geom_flat_violin(position = position_nudge(x = .2, y = 0),adjust =2)+
  geom_jitter(position = position_jitter(width = .15), shape =16)+
  geom_boxplot(width = .1, alpha = 0.5)+
  labs(title= Z.CN_main_title, y=Z.CN_y_title)+
  theme_bw()+
  theme(legend.position='none')+
  theme(axis.text.x= element_blank())+
  theme(plot.title= element_text(hjust="0.5"))+
  theme(axis.title.x = element_blank())+
  scale_y_continuous(limits=c(3.8, 16.8))+
  stat_summary(geom = 'text', label = cld.SCN.mean.part$Letters, fun = min, vjust = 1.4, color = "black")+ 
  coord_cartesian(expand=T)

# EMPTY plot for legend capture
emp <- ggplot(resdata.full, aes(x = partition, y = dens, fill=partition))+
  geom_boxplot()+
  theme_bw()+
  theme(legend.position = c(0.5,0.5),
        legend.text = element_text(size =  12),
        legend.title = element_blank())
leg <- get_legend(emp)
p16 <- as_ggplot(leg)

### Visualisation of species niche in a PCA
pca = prcomp(subsetdata2, scale=T) # (scale=TRUE) : scales the variables so that all have unit variance . This is necessary if the data1 has different units
species_pca <- fviz_pca_biplot(pca, 
               habillage = resdata.full$spe,
               label = "var", 
               geom = "point",
               alpha = 0.3,
               invisible = "quali",
               addEllipses = TRUE,
               ellipse.type = "convex",
               ellipse.alpha = 0.1, 
               legend.position = "none",
               legend.title = "Species") +
                geom_point(aes(shape = resdata.full$spe), size = 2) +
  scale_y_continuous(limits=c(-6.9, 4.5))

### Visualisation of the partition nich in the same PCA
cluster_pca <- fviz_pca_biplot(pca, 
              habillage = part,
              label = "var", 
              geom = "point",
              alpha = 0.3,
              invisible = "quali",
              addEllipses = TRUE,
              ellipse.level = 0.8,
              ellipse.alpha = 0.1, 
              legend.position = "none",
              legend.title = "Partition") +
              geom_point(shape = resdata.full$spe, size = 2) +
  scale_y_continuous(limits=c(-6.9, 4.5))

#' 
#' ## **Figures**
## ----eval = T, echo = T, message = F, warning = F-----------------------------
#### generating figure1 #### 
#figure 1 
ggarrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, 
          p10, p11, p12, p13, p14, p15, p16, ncol = 4, nrow = 4, 
          font.label = list(size = 18, color = "black", face = "bold"), 
          labels = "auto", hjust = -1.5, vjust= 1.2)

#supplementary_figure S1 
corrplot_final <- corrplot(matrice, method = "color", col = col(200),
             type = "upper",tl.cex=0.6, order = "hclust",
             addCoef.col = "black",
             tl.col = "black", tl.srt = 90,
             p.mat=res1$p,sig.level = 0.05,
             diag = FALSE, number.cex=0.6)

#figure 2
ggarrange(species_pca, cluster_pca, ncol = 1, nrow = 2, 
          font.label = list(size = 25, color = "black", face = "bold"), 
          labels = "auto", hjust = -8, vjust= 0.9)

#supplementary_figure S2 
hist(t$Best.nc[1,], breaks = max(na.omit(t$Best.nc[1,])), 
                     xlab = "Number of clusters k",
                     ylab = "Frequencies among all indices", 
                     main = "")
#figure3
plot(part~resdata.full$spe, 
                        xlab = "",
                        ylab = "Behavior frequencies among conspecifics", 
                        main = "",
                        col = c(c("#619CFF", "#00BA38", "#F98880")))

