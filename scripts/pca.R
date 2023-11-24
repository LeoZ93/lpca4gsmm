## Performing PCA

library(readr)
library(dplyr)
library(tidyr)
library(utils)


###############################
### Cancer and Tissue cells ###
###############################

# import dataset
df_1 <- read_csv("/home/users/lzehetner/data/logPCA/cancer_and_tissue.csv")

# prepare dataframe
gene_names <- df_1$Gene
df_1 <- df_1[,-c(1,2)]
df_1[is.na(df_1)] <- 0.0

# z-normalization transcriptomic data by tissue
df_1 <- as.data.frame(lapply(df_1, function(x) ((x - mean(x)) / sd(x))))

# perform PCA on transposed dataset
df_1 <- t(df_1)
data.pca <- prcomp(df_1, center = TRUE, scale = TRUE)
summary(data.pca)

# add labels for healthy and cancerous tissues
class.labels = vector()
class.labels[c(1:50)] = "#0000FF" # 
class.labels[c(51:80)] = "#FF0000" #

# create the pca plot
plot(data.pca$x[, 1:2],  # x and y data
     pch=21,           # point shape
     col=class.labels, # 
     bg=class.labels,  #
     cex=1,          # point size
     main="",     # title of plot
     xlab = "PC1 (60%)", # PC1 = 60 w/o Adrenal Cancer
     ylab = "PC2 (8%)",
     cex.axis = 1,
     cex.lab = 1
)

colors <- c("#0000FF", "#FF0000")
legend("bottomleft", legend=c("Normal Tissue", "Cancer"), fill=colors, cex=1)



######################################################
### All metabolic Genes in Cancer and Tissue cells ###
######################################################

# import the genes dataset from human1 repository
human_gem_genes <- read_delim("data/human1/genes.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)

# prepare dataset
gene_names <- human_gem_genes$genes

# import cancer and tissue transcriptomic data
df <- read_csv("/home/users/lzehetner/data/logPCA/cancer_and_tissue.csv")

# rename genes column
human_gem_genes <- human_gem_genes %>%
  rename(Gene = genes)

# extract all metabolic genes from transcriptomic data
metab_gene_expr <- human_gem_genes %>%
  left_join(df, by = "Gene" )

# remove unnecessary columns and prepare dataframe
df_1 <- metab_gene_expr[, -c(1:11)]

df_1[is.na(df_1)] <- 0.0

# z-normalization of transcriptomic data
df_1 <- as.data.frame(lapply(df_1, function(x) ((x - mean(x)) / sd(x))))

# perform PCA on transposed dataset
df_1 <- t(df_1)

data.pca <- prcomp(df_1, center = TRUE, scale = TRUE)

summary(data.pca)

# add labels for healthy and cancer tissues
class.labels = vector()
class.labels[c(1:50)] = "#0000FF" # 
class.labels[c(51:80)] = "#FF0000" #

# create pca plot for metabolic genes from human1
plot(data.pca$x[, 1:2],  # x and y data
     pch=21,           # point shape
     col=class.labels, # 
     bg=class.labels,  #
     cex=1.5,          # point size
     main="",     # title of plot
     xlab = "PC1 (56%)", # PC1 = 53 w/o adrenal cancer and 56 w. adrenal cancer
     ylab = "PC2 (7%)",
     cex.axis = 1.5,
     cex.lab = 1.5
)

colors <- c("#0000FF", "#FF0000")
legend("bottomleft", legend=c("Normal Tissue", "Cancer"), fill=colors, cex=1)

####################################################################
### Genes from differential reactions in Cancer and Tissue cells ###
####################################################################

# import differential genes dataset from extracted reactions
df_1 <- read_csv("/home/users/lzehetner/data/logPCA/differential_genes_tissue_and_cancer.csv")

# prepare dataset
gene_names <- df_1$Gene
df_1 <- df_1[,-c(1,2)]

df_1[is.na(df_1)] <- 0.0

# z-normalization of transcriptomic data
df_1 <- as.data.frame(lapply(df_1, function(x) ((x - mean(x)) / sd(x))))

# perform pca on normalized transcriptomic data
df_1 <- t(df_1)

data.pca <- prcomp(df_1, center = TRUE, scale = TRUE)

summary(data.pca)

# add labels for healthy and cancerous tissue
class.labels = vector()
class.labels[c(1:50)] = "#0000FF" # 
class.labels[c(51:80)] = "#FF0000" #

# create plot using differential genes across tissue specific gsmm
plot(data.pca$x[, 1:2],  # x and y data
     pch=21,           # point shape
     col=class.labels, # 
     bg=class.labels,  #
     cex=1.5,          # point size
     main="",     # title of plot
     xlab = "PC1 (47%)", # PC1 = 46 w/o adrenal cancer
     ylab = "PC2 (7%)",
     cex.axis = 1.5,
     cex.lab = 1.5
)


##################################################################
### PCA on simulated growth rates from 222 Escherichia strains ###
##################################################################

# import dataset, obtained from Monk, 2022
df <- read_excel("data/logPCA/rstb20210236_si_003.xlsx", sheet = "growth predictions")

# prepare dataframe
df_1 <- df[, -c(1:4)]

# perform pca using a correlation-based approach to obtain a similar clustering as in Monk, 2022
df_1 <- t(df_1)

corr_matrix <- cor(df_1, use = "everything", method = "spearman")
df_1[is.na(df_1)] <- 0.0
data.pca <- princomp(corr_matrix)

summary(data.pca)

# add labels for each clade. Same colors are taken as from the original publication in Monk, 2022
class.labels = vector()
class.labels[c(154:207)] = "#127db9" # E.albertii
class.labels[c(128:135)] = "#2c9f28" # Clade III
class.labels[c(90:92)] = "#ffe4ca" # S.flexneri
class.labels[c(1:87, 93, 94)] = "#fd7d08" # E.coli
class.labels[c(88:89)] = "#feb166" # S.dysenteriae
class.labels[c(95:135)] = "#abcaff" # Clade V
class.labels[c(136:153)] = "#e46864" # E.fergusonii
class.labels[c(222)] = "#915347" # Clade VIII
class.labels[c(216:219)] = "#c5aedb" # Clade II
class.labels[c(208:215)] = "#90e585" # Clade IV
class.labels[c(221)] = "#feb0d2" # Clade VI
class.labels[c(220)] = "#9bdeef" # Clade VII

# create the reconstructed pca plot, compare to Monk, 2022
plot(data.pca$scores[, 1:2],  # x and y data
     pch=21,           # point shape
     col=class.labels, # 
     bg=class.labels,  #
     cex=1,          # point size
     main="",     # title of plot
     xlab = "PC1 (49%)", # 
     ylab = "PC2 (20%)",
     cex.axis = 1.3,
     cex.lab = 1.3
)

colors <- c("#fd7d08", "#e46864", "#127db9","#feb166","#ffe4ca","#c5aedb","#2c9f28","#90e585","#abcaff","#feb0d2","#9bdeef", "#915347")
legend("bottomleft", legend = c("E.coli", "E.fergusonii", "E.albertii", "S.dysenteriae", "S.flexneri", "Clade II", "Clade III", "Clade IV", "Clade V", "Clade VI", "Clade VII", "Clade VIII"), fill=colors, cex=0.7)



