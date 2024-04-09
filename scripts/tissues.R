## Script to analyze the tissue dataset 
## This script includes data preprocessing, LPCA, t-sne, Jaccard similarity, phylogeny and cophenetic correlation


library(logisticPCA)
library(readr)
library(rARPACK)
library(readxl)
library(ggplot2)
library(reshape2)
library(ComplexHeatmap)
library(ape) 
library(Rtsne)
library(dplyr)
library(nnet)

###############################
### Cancer and Tissue cells ###
###############################

subsystems <- read_csv("/home/users/lzehetner/data/logPCA/subsystems.csv", col_names = TRUE)
subsystems <- subsystems[, -c(1)]

spec_rxns <- read_csv("/home/users/lzehetner/data/logPCA/differential_rxns_tissue_and_cancer.csv", col_names = TRUE)

spec_cell_rxns <- spec_rxns[, -c(1,2)]

# to remove Adrenocortical cancer, since it clusters far away
# spec_cell_rxns <- spec_cell_rxns[, -c(51)]


spec_cell_rxns_t <- t(spec_cell_rxns)

class.labels = vector()
class.labels[c(1:50)] = "#0000FF" # Normal tissue
class.labels[c(51:80)] = "#FF0000" # Cancer

shape.labels = vector("numeric")
shape.labels[c(1:50)] = 19      # Solid circle
shape.labels[c(51:80)] = 17 # Triangle point up

K = 1

## logistic PCA model
logpca.model = logisticPCA(spec_cell_rxns_t, # binary data
                           k=K, # number of PCs
                           m=0, # approximation of natural parameter
                           main_effects = TRUE,
                           partial_decomp = TRUE) # including offset term

PC1 <- logpca.model$prop_deviance_expl

K = 2

## logistic PCA model
logpca.model = logisticPCA(spec_cell_rxns_t, # binary data
                           k=K, # number of PCs
                           m=0, # approximation of natural parameter
                           main_effects = TRUE,
                           partial_decomp = TRUE) # including offset term
# here PC1 <- 0.3586717 - PC2
PC2 <- logpca.model$prop_deviance_expl - PC1

cancer.vs.tissue.model <- logpca.model

logpca.scores = cancer.vs.tissue.model$PCs # extract score matrix
logpca.loadings = cancer.vs.tissue.model$U # extract loading matrix

x <- as.data.frame(logpca.loadings)
x$rxns <- spec_rxns$Reaction
x <- merge(x, subsystems, by = "rxns", all.x = TRUE)

loading_results <- x %>%
  group_by(subsystems) %>%
  summarize(
    loading_1 = mean(V1),
    loading_2 = mean(V2)
  )

loading_results$length <- sqrt(loading_results$loading_1^2 + loading_results$loading_2^2)

# Sort the dataframe by 'length' in descending order
sorted_result <- loading_results %>% arrange(desc(length))

pathway_importance_lpca_tissue <- sorted_result[, c(1,4)]

# Take the top 10 longest vectors
top_longest_vectors <- head(sorted_result, 5)

png("tissue_diff_rxns_lpca.png", width = 1200, height = 700)

plot(logpca.scores,  # x and y data
     pch=shape.labels,           # point shape
     col=class.labels, # 
     bg=class.labels,  #
     cex=1.8,          # point size
     main="",     # title of plot
     xlab = "PC1 (20%)", # for K = 1, m = 0 the PC = 0.163   sum(PC1, PC2) = 0.273
     ylab = "PC2 (15%)", # for K = 2, m = 0, the PC = 0.11
     cex.axis = 1.5,
     cex.lab = 1.5
)

arrows(x0 = rep(0, nrow(top_longest_vectors)), 
       y0 = rep(0, nrow(top_longest_vectors)), 
       x1 = top_longest_vectors$loading_1*2000, 
       y1 = top_longest_vectors$loading_2*2000, 
       col = "black", # You can change the color
       angle = 25, # Angle of the arrow head
       length = 0.1) # Length of the arrow head

text(x = top_longest_vectors$loading_1[1]*2000, 
     y = top_longest_vectors$loading_2[1]*3000, 
     labels = top_longest_vectors$subsystems[1],
     pos = 2,  # Places the text above the point
     col = "black",  # Text color
     cex = 1.3  # Text size
)

text(x = top_longest_vectors$loading_1[2]*2500, 
     y = top_longest_vectors$loading_2[2]*2500, 
     labels = top_longest_vectors$subsystems[2],
     pos = 2,  # Places the text above the point
     col = "black",  # Text color
     cex = 1.3  # Text size
)

text(x = top_longest_vectors$loading_1[3]*2000, 
     y = top_longest_vectors$loading_2[3]*2500, 
     labels = top_longest_vectors$subsystems[3],
     pos = 2,  # Places the text above the point
     col = "black",  # Text color
     cex = 1.3  # Text size
)

text(x = top_longest_vectors$loading_1[4]*2500, 
     y = top_longest_vectors$loading_2[4]*2300, 
     labels = top_longest_vectors$subsystems[4],
     pos = 2,  # Places the text above the point
     col = "black",  # Text color
     cex = 1.3  # Text size
)

text(x = top_longest_vectors$loading_1[5]*2500, 
     y = top_longest_vectors$loading_2[5]*3000, 
     labels = top_longest_vectors$subsystems[5],
     pos = 2,  # Places the text above the point
     col = "black",  # Text color
     cex = 1.3  # Text size
)

dev.off()

selected_rxns <- x %>%
  filter(subsystems == "Sulfur metabolism")


###########################################
###### Jaccard cancer/healthy tissue ######


# Calculate Jaccard similarity matrix for the binary matrix
df <- calculate_jaccard_similarity(spec_cell_rxns_t)

group_assignments <- c(rep(1, 50), rep(2, 30))

group_factors <- factor(group_assignments, levels = c(1, 2), labels = c("Healthy", "Cancer"))

group_colors <- c("Healthy" = "blue", "Cancer" = "red")

custom_colors <- colorRampPalette(c("white", "yellow", "blue", "black"))(100)

annotation_col <- HeatmapAnnotation(
  Group = group_factors,
  col = list(Group = group_colors),
  which = "col" 
)

png("tissue_diff_rxns_jacc.png", width = 1200, height = 1200)

Heatmap(df, 
        top_annotation = annotation_col, 
        name = "Similarity", 
        column_title = "",
        col = custom_colors,
        column_names_side = "bottom",
        column_names_gp = gpar(col = NA))

dev.off()

#########################################
###### tsne cancer/healthy tissues ######

tsne <- Rtsne(spec_cell_rxns_t, dims = 2, distance = 'hamming', check_duplicates = FALSE, perplexity = 10)

tsne_data <- as.data.frame(tsne$Y)
colnames(tsne_data) <- c('tSNE1', 'tSNE2')


tsne_data$color <- class.labels

# if all reactions -> change to "escherichia_all_rxns_tsne.png"
png(filename = "tissue_diff_rxns_tsne.png", width = 1200, height = 700)

plot(tsne_data$tSNE1,
     tsne_data$tSNE2,
     pch=shape.labels,           # point shape
     col=class.labels, # 
     bg=class.labels,  #
     cex=1.8,          # point size
     main="",     # title of plot
     xlab = "tSNE1", # 
     ylab = "tSNE2", # 
     cex.axis = 1.5,
     cex.lab = 1.5
)


dev.off()

########################################
###------------ legend --------------###

group_names <- c("Healthy Tissue", 
                 "Cancer Tissue" 
)

# Vector with color assignments for each group
class.labels <- c("Healthy Tissue" = "#0000FF", 
                  "Cancer Tissue" = "#FF0000"
)

# Vector with shape assignments for each group
shape.labels <- c("Healthy Tissue" = 19, 
                  "Cancer Tissue" = 17
)


png(filename = "tissue_legend.png", width = 1000, height = 300)

plot(1, type="n", xlab="", ylab="", xlim=c(0, 4), ylim=c(1, 2), xaxt='n', yaxt='n', bty='n', main="Legend")

# Loop through each group to add symbols (shapes) and colors to the legend, with 4 per row
for(i in 1:length(group_names)) {
  # Calculate x and y positions for 4 items per row
  x_pos <- 0.5 + ((i-1) %% 4)  # Shifts right for each new item, resets every 4 items
  y_pos <- 2 - 0.1 * floor((i-1) / 4)  # Adjusted to decrease space between rows
  
  # Plot symbol
  points(x_pos, y_pos, pch=shape.labels[group_names[i]], col=class.labels[group_names[i]], cex=1.5)
  
  # Add text, closer to the symbol
  text(x_pos + 0.2, y_pos, labels=group_names[i], pos=4, cex=1, offset=0.1)
}

dev.off()

#########################################################
#####-----subsystem analysis using logistic GLM-----#####

subsystems <- read_csv("/home/users/lzehetner/data/logPCA/subsystems.csv", col_names = TRUE)
subsystems <- subsystems[, -c(1)]

spec_rxns <- read_csv("/home/users/lzehetner/data/logPCA/differential_rxns_tissue_and_cancer.csv", col_names = TRUE)

differential.rxns <- spec_rxns$Reaction

spec_rxns <- spec_rxns[, -c(1,2)]

spec_rxns_t <- as.data.frame(t(spec_rxns))

spec_rxns_t$Tissue <- NA

# Assign "Healthy Tissue" to the first 50 rows in the new column
spec_rxns_t$Tissue[1:50] <- "Healthy"

# Assign "Cancerous Tissue" to rows 51 to 80 in the new column
spec_rxns_t$Tissue[51:80] <- "Cancer"

# Fit a multinomial logistic regression model
multinom_model_tissue <- multinom(Tissue ~ ., data = spec_rxns_t, MaxNWts = 30000)

log_model_smry_tissue <- summary(multinom_model_tissue)

log_model_coeff <- as.data.frame(log_model_smry_tissue[["coefficients"]])

log_model_coeff <- as.data.frame(log_model_coeff[-c(1),])

log_model_coeff$rxns <- differential.rxns
log_model_coeff <- merge(log_model_coeff, subsystems, by = "rxns", all.x = TRUE)

# Then, aggregate these squared sum importance scores by pathway -> use the mean instead
pathway_importance_mlr <- log_model_coeff %>%
  group_by(subsystems) %>%
  summarise(total_mean_importance = mean(`log_model_coeff[-c(1), ]`, na.rm = TRUE))

#### Compare MLR and LPCA subsystem enrichment ####

merged_rankings <- merge(pathway_importance_mlr, pathway_importance_lpca_tissue, by = "subsystems")

colnames(merged_rankings)[2] = "MLR"
colnames(merged_rankings)[3] = "LPCA"

merged_rankings$MLR <- merged_rankings$MLR / max(merged_rankings$MLR)
merged_rankings$LPCA <- merged_rankings$LPCA / max(merged_rankings$LPCA)

heatmap_data <- as.matrix(merged_rankings[, 2:3])
rownames(heatmap_data) <- merged_rankings$subsystems


color_palette <- colorRampPalette(c("white", "yellow", "blue", "darkblue", "black"))(100)

png("comparison_lpca_mlr_tissues.png", width = 2500, height = 3200)

heatmap.2(heatmap_data, scale = "none", trace = "none", col = color_palette, 
          margin = c(30, 50), cexRow = 2, cexCol = 2.5,
          key = FALSE, density.info = "none",
          main = "", xlab = "", ylab = "",
          Colv = FALSE)

dev.off()

####################################
### legend mlr and lpca analysis ###


data <- data.frame(Value = seq(0, 1, by = 0.01))

png("legend_mlr_lpca.png", width = 1300, height = 200)

ggplot(data, aes(x = Value, y = 1, fill = Value)) +
  geom_tile(stat = "identity") +
  scale_fill_gradientn(colors = color_palette) +
  theme_minimal() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(size = 70)) +
  labs(fill = "Val")

dev.off()


###############################################################
### All metabolic Genes in cancer and healthy tissues cells ###

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

shape.labels = vector("numeric")
shape.labels[c(1:50)] = 19 # plus
shape.labels[c(51:80)] = 17  # Triangle point up

png("tissue_metab_genes_pca.png", width = 1200, height = 700)

plot(data.pca$x[, 1:2],  # x and y data
     pch=shape.labels,           # point shape
     col=class.labels, # 
     bg=class.labels,  #
     cex=1.5,          # point size
     main="",     # title of plot
     xlab = "PC1 (56%)", # PC1 = 53 w/o adrenal cancer and 56 w. adrenal cancer
     ylab = "PC2 (7%)",
     cex.axis = 1.4,
     cex.lab = 1.4
)

dev.off()

arrows(x0 = rep(0, nrow(top_longest_vectors)), 
       y0 = rep(0, nrow(top_longest_vectors)), 
       x1 = top_longest_vectors$loading_1*200, 
       y1 = top_longest_vectors$loading_2*200, 
       col = "black", # You can change the color
       angle = 25, # Angle of the arrow head
       length = 0.1) # Length of the arrow head

for(i in 1:nrow(top_longest_vectors)) {
  jitter_factor <- 40  # You can adjust this value based on how much jitter you want.
  
  # Calculate jittered positions for the labels
  x_jittered <- jitter(top_longest_vectors$loading_1[i]*200, factor=jitter_factor)
  y_jittered <- jitter(top_longest_vectors$loading_2[i]*200, factor=jitter_factor)
  
  if (top_longest_vectors$loading_1[i] < 0) {
    text(x = x_jittered, 
         y = y_jittered, 
         labels = top_longest_vectors$subsystems[i],
         pos = 4,  # Places the text above the point
         col = "black",  # Text color
         cex = 1.3  # Text size
    )
  }
  if (top_longest_vectors$loading_1[i] > 0) {
    text(x = x_jittered, 
         y = y_jittered, 
         labels = top_longest_vectors$subsystems[i],
         pos = 4,  # Places the text above the point
         col = "black",  # Text color
         cex = 1.3  # Text size
    )
  }
}

colors <- c("#0000FF", "#FF0000")
legend("bottomleft", legend=c("Normal Tissue", "Cancer"), fill=colors, cex=1.3)

#######################################################################
### Genes from differential reactions in Cancer and healthy tissues ###

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

shape.labels = vector("numeric")
shape.labels[c(1:50)] = 19 # plus
shape.labels[c(51:80)] = 17  # Triangle point up

png("tissue_metab_genes_from_gsmms_pca.png", width = 1200, height = 700)

plot(data.pca$x[, 1:2],  # x and y data
     pch=shape.labels,           # point shape
     col=class.labels, # 
     bg=class.labels,  #
     cex=1.5,          # point size
     main="",     # title of plot
     xlab = "PC1 (56%)", # PC1 = 53 w/o adrenal cancer and 56 w. adrenal cancer
     ylab = "PC2 (7%)",
     cex.axis = 1.4,
     cex.lab = 1.4
)

dev.off()


arrows(x0 = rep(0, nrow(top_longest_vectors)), 
       y0 = rep(0, nrow(top_longest_vectors)), 
       x1 = top_longest_vectors$loading_1*50, 
       y1 = top_longest_vectors$loading_2*50, 
       col = "black", # You can change the color
       angle = 25, # Angle of the arrow head
       length = 0.1) # Length of the arrow head

for(i in 1:nrow(top_longest_vectors)) {
  jitter_factor <- 40  # You can adjust this value based on how much jitter you want.
  
  # Calculate jittered positions for the labels
  x_jittered <- jitter(top_longest_vectors$loading_1[i]*50, factor=jitter_factor)
  y_jittered <- jitter(top_longest_vectors$loading_2[i]*50, factor=jitter_factor)
  
  if (top_longest_vectors$loading_1[i] < 0) {
    text(x = x_jittered, 
         y = y_jittered, 
         labels = top_longest_vectors$subsystems[i],
         pos = 2,  # Places the text above the point
         col = "black",  # Text color
         cex = 1.3  # Text size
    )
  }
  if (top_longest_vectors$loading_1[i] > 0) {
    text(x = x_jittered, 
         y = y_jittered, 
         labels = top_longest_vectors$subsystems[i],
         pos = 4,  # Places the text above the point
         col = "black",  # Text color
         cex = 1.3  # Text size
    )
  }
}

colors <- c("#0000FF", "#FF0000")
legend("bottomleft", legend=c("Normal Tissue", "Cancer"), fill=colors, cex=1)


###########################################################
### Whole transcriptome from cancer and healthy tissues ###


df_1 <- read_csv("/home/users/lzehetner/data/logPCA/cancer_and_tissue.csv")

gene_names <- df_1$Gene

df_1 <- df_1[,-c(1,2)]
#df_1 <- df_1[, -c(51)] # to remove Adreno cancer since it was and outlier in logPCA

df_1[is.na(df_1)] <- 0.0

df_1 <- as.data.frame(lapply(df_1, function(x) ((x - mean(x)) / sd(x))))

df_1 <- t(df_1)

data.pca <- prcomp(df_1, center = TRUE, scale = TRUE)

summary(data.pca)

class.labels = vector()
class.labels[c(1:50)] = "#0000FF" # 
class.labels[c(51:80)] = "#FF0000" #

shape.labels = vector("numeric")
shape.labels[c(1:50)] = 19 # plus
shape.labels[c(51:80)] = 17  # Triangle point up

png("tissue_pca.png", width = 1200, height = 700)

plot(data.pca$x[, 1:2],  # x and y data
     pch=shape.labels,           # point shape
     col=class.labels, # 
     bg=class.labels,  #
     cex=1.5,          # point size
     main="",     # title of plot
     xlab = "PC1 (56%)", # PC1 = 53 w/o adrenal cancer and 56 w. adrenal cancer
     ylab = "PC2 (7%)",
     cex.axis = 1.4,
     cex.lab = 1.4
)

dev.off()

arrows(x0 = rep(0, nrow(top_longest_vectors)), 
       y0 = rep(0, nrow(top_longest_vectors)), 
       x1 = top_longest_vectors$loading_1*1000, 
       y1 = top_longest_vectors$loading_2*1000, 
       col = "black", # You can change the color
       angle = 25, # Angle of the arrow head
       length = 0.1) # Length of the arrow head

for(i in 1:nrow(top_longest_vectors)) {
  jitter_factor <- 40  # You can adjust this value based on how much jitter you want.
  
  # Calculate jittered positions for the labels
  x_jittered <- jitter(top_longest_vectors$loading_1[i]*1000, factor=jitter_factor)
  y_jittered <- jitter(top_longest_vectors$loading_2[i]*1000, factor=jitter_factor)
  
  if (top_longest_vectors$loading_1[i] < 0) {
    text(x = x_jittered, 
         y = y_jittered, 
         labels = top_longest_vectors$subsystems[i],
         pos = 2,  # Places the text above the point
         col = "black",  # Text color
         cex = 1.3  # Text size
    )
  }
  if (top_longest_vectors$loading_1[i] > 0) {
    text(x = x_jittered, 
         y = y_jittered, 
         labels = top_longest_vectors$subsystems[i],
         pos = 4,  # Places the text above the point
         col = "black",  # Text color
         cex = 1.3  # Text size
    )
  }
}

colors <- c("#0000FF", "#FF0000")
legend("bottomleft", legend=c("Normal Tissue", "Cancer"), fill=colors, cex=1)

