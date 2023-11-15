## Performing PCA

library(readr)
library(dplyr)
library(tidyr)
library(utils)


###############################
### Cancer and Tissue cells ###
###############################

df_1 <- read_csv("/home/users/lzehetner/data/logPCA/cancer_and_tissue.csv")

gene_names <- df_1$Gene

df_1 <- df_1[,-c(1,2)]
#df_1 <- df_1[, -c(51)] # to remove Adreno cancer since it was and outlier in logPCA

df_1[is.na(df_1)] <- 0.0

# df_1 <- t(df_1)

df_1 <- as.data.frame(lapply(df_1, function(x) ((x - mean(x)) / sd(x))))

df_1 <- t(df_1)

# df_1 <- df_1 / sapply(df_1, sum) * 100

#df_1 <- t(df_1)

data.pca <- prcomp(df_1, center = TRUE, scale = TRUE)

summary(data.pca)

loadings <- as.data.frame(data.pca$rotation[, c(1,2)])

loadings$Gene <- gene_names

Ensembl2Reactome <- read.delim("~/data/logPCA/Ensembl2Reactome.txt", header=FALSE)

filtered_df <- Ensembl2Reactome %>%
  filter(grepl("^ENSG", V1) & V6 == "Homo sapiens")

filtered_df <- filtered_df %>%
  rename(Gene = V1)

result_df <- loadings %>%
  left_join(filtered_df, by = "Gene" ) %>%
  select(-V2, -V3, -V5, -V6) %>%  # adjust this to remove unwanted columns
  unnest(V4, keep_empty = TRUE) 

result_df <- result_df %>%
  rename(subsystems = V4)


loading_results <- result_df %>%
  group_by(subsystems) %>%
  summarize(
    loading_1 = mean(PC1),
    loading_2 = mean(PC2)
  )

loading_results$length <- sqrt(loading_results$loading_1^2 + loading_results$loading_2^2)

# Sort the dataframe by 'length' in descending order
sorted_result <- loading_results %>% arrange(desc(length))

# Take the top 10 longest vectors
top_longest_vectors <- head(sorted_result, 5)

class.labels = vector()
class.labels[c(1:50)] = "#0000FF" # 
class.labels[c(51:80)] = "#FF0000" #

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

arrows(x0 = rep(0, nrow(top_longest_vectors)), 
       y0 = rep(0, nrow(top_longest_vectors)), 
       x1 = top_longest_vectors$loading_1*1000, 
       y1 = top_longest_vectors$loading_2*1000, 
       col = "black", # You can change the color
       angle = 25, # Angle of the arrow head
       length = 0.1) # Length of the arrow head

for(i in 1:nrow(top_longest_vectors)) {
  jitter_factor <- 30  # You can adjust this value based on how much jitter you want.
  
  # Calculate jittered positions for the labels
  x_jittered <- jitter(top_longest_vectors$loading_1[i]*1000, factor=jitter_factor)
  y_jittered <- jitter(top_longest_vectors$loading_2[i]*1000, factor=jitter_factor)
  
  if (top_longest_vectors$loading_1[i] < 0) {
    text(x = x_jittered, 
         y = y_jittered, 
         labels = top_longest_vectors$subsystems[i],
         pos = 2,  # Places the text above the point
         col = "black",  # Text color
         cex = 1  # Text size
    )
  }
  if (top_longest_vectors$loading_1[i] > 0) {
    text(x = x_jittered, 
         y = y_jittered, 
         labels = top_longest_vectors$subsystems[i],
         pos = 4,  # Places the text above the point
         col = "black",  # Text color
         cex = 1  # Text size
    )
  }
}

colors <- c("#0000FF", "#FF0000")
legend("bottomleft", legend=c("Normal Tissue", "Cancer"), fill=colors, cex=1)



######################################################
### All metabolic Genes in Cancer and Tissue cells ###
######################################################

human_gem_genes <- read_delim("data/human1/genes.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)

gene_names <- human_gem_genes$genes

df <- read_csv("/home/users/lzehetner/data/logPCA/cancer_and_tissue.csv")

#df <- df[, -c(53)] # here adrenocortical cancer is in col 53, since the genes are still in the df

human_gem_genes <- human_gem_genes %>%
  rename(Gene = genes)

metab_gene_expr <- human_gem_genes %>%
  left_join(df, by = "Gene" )

df_1 <- metab_gene_expr[, -c(1:11)]

df_1[is.na(df_1)] <- 0.0

df_1 <- as.data.frame(lapply(df_1, function(x) ((x - mean(x)) / sd(x))))

df_1 <- t(df_1)

# df_1 <- t(df_1)

data.pca <- prcomp(df_1, center = TRUE, scale = TRUE)

summary(data.pca)

loadings <- as.data.frame(data.pca$rotation[, c(1,2)])

class.labels = vector()
class.labels[c(1:50)] = "#0000FF" # 
class.labels[c(51:80)] = "#FF0000" #

loadings$Gene <- gene_names

Ensembl2Reactome <- read.delim("~/data/logPCA/Ensembl2Reactome.txt", header=FALSE)

filtered_df <- Ensembl2Reactome %>%
  filter(grepl("^ENSG", V1) & V6 == "Homo sapiens")

filtered_df <- filtered_df %>%
  rename(Gene = V1)

result_df <- loadings %>%
  left_join(filtered_df, by = "Gene" ) %>%
  select(-V2, -V3, -V5, -V6) %>%  # adjust this to remove unwanted columns
  unnest(V4, keep_empty = TRUE) 

result_df <- result_df %>%
  rename(subsystems = V4)


loading_results <- result_df %>%
  group_by(subsystems) %>%
  summarize(
    loading_1 = mean(PC1),
    loading_2 = mean(PC2)
  )

loading_results$length <- sqrt(loading_results$loading_1^2 + loading_results$loading_2^2)

# Sort the dataframe by 'length' in descending order
sorted_result <- loading_results %>% arrange(desc(length))

# Take the top 10 longest vectors
top_longest_vectors <- head(sorted_result, 5)

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

arrows(x0 = rep(0, nrow(top_longest_vectors)), 
       y0 = rep(0, nrow(top_longest_vectors)), 
       x1 = top_longest_vectors$loading_1*200, 
       y1 = top_longest_vectors$loading_2*200, 
       col = "black", # You can change the color
       angle = 25, # Angle of the arrow head
       length = 0.1) # Length of the arrow head

for(i in 1:nrow(top_longest_vectors)) {
  jitter_factor <- 20  # You can adjust this value based on how much jitter you want.
  
  # Calculate jittered positions for the labels
  x_jittered <- jitter(top_longest_vectors$loading_1[i]*200, factor=jitter_factor)
  y_jittered <- jitter(top_longest_vectors$loading_2[i]*200, factor=jitter_factor)
  
  if (top_longest_vectors$loading_1[i] < 0) {
    text(x = x_jittered, 
         y = y_jittered, 
         labels = top_longest_vectors$subsystems[i],
         pos = 2,  # Places the text above the point
         col = "black",  # Text color
         cex = 0.8  # Text size
    )
  }
  if (top_longest_vectors$loading_1[i] > 0) {
    text(x = x_jittered, 
         y = y_jittered, 
         labels = top_longest_vectors$subsystems[i],
         pos = 4,  # Places the text above the point
         col = "black",  # Text color
         cex = 0.8  # Text size
    )
  }
}

colors <- c("#0000FF", "#FF0000")
legend("bottomleft", legend=c("Normal Tissue", "Cancer"), fill=colors, cex=1)

####################################################################
### Genes from differential reactions in Cancer and Tissue cells ###
####################################################################

df_1 <- read_csv("/home/users/lzehetner/data/logPCA/differential_genes_tissue_and_cancer.csv")

gene_names <- df_1$Gene

df_1 <- df_1[,-c(1,2)]
# df_1 <- df_1[, -c(51)] # to remove Adreno cancer since it was and outlier in logPCA

df_1[is.na(df_1)] <- 0.0

df_1 <- as.data.frame(lapply(df_1, function(x) ((x - mean(x)) / sd(x))))

df_1 <- t(df_1)

# df_1 <- t(df_1)

data.pca <- prcomp(df_1, center = TRUE, scale = TRUE)

summary(data.pca)

loadings <- as.data.frame(data.pca$rotation[, c(1,2)])

class.labels = vector()
class.labels[c(1:50)] = "#0000FF" # 
class.labels[c(51:80)] = "#FF0000" #

loadings$Gene <- gene_names

Ensembl2Reactome <- read.delim("~/data/logPCA/Ensembl2Reactome.txt", header=FALSE)

filtered_df <- Ensembl2Reactome %>%
  filter(grepl("^ENSG", V1) & V6 == "Homo sapiens")

filtered_df <- filtered_df %>%
  rename(Gene = V1)

result_df <- loadings %>%
  left_join(filtered_df, by = "Gene" ) %>%
  select(-V2, -V3, -V5, -V6) %>%  # adjust this to remove unwanted columns
  unnest(V4, keep_empty = TRUE) 

result_df <- result_df %>%
  rename(subsystems = V4)


loading_results <- result_df %>%
  group_by(subsystems) %>%
  summarize(
    loading_1 = mean(PC1),
    loading_2 = mean(PC2)
  )

loading_results$length <- sqrt(loading_results$loading_1^2 + loading_results$loading_2^2)

# Sort the dataframe by 'length' in descending order
sorted_result <- loading_results %>% arrange(desc(length))

# Take the top 10 longest vectors
top_longest_vectors <- head(sorted_result, 5)

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

arrows(x0 = rep(0, nrow(top_longest_vectors)), 
       y0 = rep(0, nrow(top_longest_vectors)), 
       x1 = top_longest_vectors$loading_1*50, 
       y1 = top_longest_vectors$loading_2*50, 
       col = "black", # You can change the color
       angle = 25, # Angle of the arrow head
       length = 0.1) # Length of the arrow head

for(i in 1:nrow(top_longest_vectors)) {
  jitter_factor <- 20  # You can adjust this value based on how much jitter you want.
  
  # Calculate jittered positions for the labels
  x_jittered <- jitter(top_longest_vectors$loading_1[i]*50, factor=jitter_factor)
  y_jittered <- jitter(top_longest_vectors$loading_2[i]*50, factor=jitter_factor)
  
  if (top_longest_vectors$loading_1[i] < 0) {
    text(x = x_jittered, 
         y = y_jittered, 
         labels = top_longest_vectors$subsystems[i],
         pos = 2,  # Places the text above the point
         col = "black",  # Text color
         cex = 0.8  # Text size
    )
  }
  if (top_longest_vectors$loading_1[i] > 0) {
    text(x = x_jittered, 
         y = y_jittered, 
         labels = top_longest_vectors$subsystems[i],
         pos = 4,  # Places the text above the point
         col = "black",  # Text color
         cex = 0.8  # Text size
    )
  }
}

colors <- c("#0000FF", "#FF0000")
legend("bottomleft", legend=c("Normal Tissue", "Cancer"), fill=colors, cex=1)


##################################################################
### PCA on simulated growth rates from 222 Escherichia strains ###
##################################################################

df <- read_excel("data/logPCA/rstb20210236_si_003.xlsx", sheet = "growth predictions")

df_1 <- df[, -c(1:4)]

#df_1 <- as.data.frame(lapply(df_1, function(x) ((x - mean(x)) / sd(x))))

#df_1[is.na(df_1)] <- 0.0

# df_1 <- t(df_1)

# df_1 <- df_1 / sapply(df_1, sum) * 100

df_1 <- t(df_1)

# data.pca <- prcomp(df_1, center = FALSE, scale = FALSE)

corr_matrix <- cor(df_1, use = "everything", method = "spearman")
df_1[is.na(df_1)] <- 0.0
data.pca <- princomp(corr_matrix)

summary(data.pca)

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



