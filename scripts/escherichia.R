## Script to analyze the Escherichia dataset derived from Monk (2022)
## This script includes data preprocessing, LPCA, t-sne, Jaccard similarity, phylogenty, pca, cophenetic correlation, and mlr analyses


library(logisticPCA) # for logistic PCA ans logistic SVD
library(readr)
library(rARPACK)
library(readxl)
library(ggplot2)
library(reshape2)
library(ComplexHeatmap)
library(ape) 
library(Rtsne)
library(dplyr)
library(ellipse)

#install.packages("nnet", quiet = TRUE)
library(nnet)

set.seed(42)

###############################
##### PanGEMs Escherichia #####
###############################

escherichia <- read_excel("data/logPCA/rstb20210236_si_002.xlsx", sheet = "Reaction Info")
clades <- read_excel("data/logPCA/rstb20210236_si_002.xlsx", sheet = "Model reaction counts")

escherichia <- escherichia[, -c(2:5)]

escherichia[is.na(escherichia)] <- 0.0

# skip this command if analysis of all reactions is of interest
escherichia <- escherichia[!(apply(escherichia[,3:224], 1, function(row) all(row == 1) | all(row == 0))), ]

differential.rxns <- escherichia$...1
subsystems <- escherichia$Subsystem

escherichia <- escherichia[, -c(1,2)]

matching_cols <- match(colnames(escherichia), clades$strain_id)

unique_groups <- unique(clades$Species)
grouped_uid_lists <- lapply(unique_groups, function(group) {
  return(clades$strain_id[clades$Species == group])
})

index_lists <- lapply(grouped_uid_lists, function(uid_list) {
  return(which(colnames(escherichia) %in% uid_list))
})

E.albertii <- index_lists[[1]]
Clade.III <- index_lists[[2]]
S.flexneri <- index_lists[[3]]
E.coli <- index_lists[[4]]
S.dysenteriae <- index_lists[[5]]
Clade.V <- index_lists[[6]]
E.fergusonii <- c(index_lists[[7]], index_lists[[10]])
Clade.VIII <- index_lists[[8]]
Clade.II <- index_lists[[9]]
Clade.IV <- index_lists[[11]]
Clade.VI <- index_lists[[12]]
Clade.VII <- index_lists[[13]]

rem.clades.cluster <- c(index_lists[[3]], index_lists[[4]], index_lists[[5]], index_lists[[6]], index_lists[[8]], index_lists[[12]], index_lists[[13]])

escherichia_t <- t(escherichia)

class.labels = vector()
class.labels[c(E.coli)] = "#fd7d08" # 
class.labels[c(E.fergusonii)] = "#e46864" #
class.labels[c(E.albertii)] = "#127db9" # 
class.labels[c(S.dysenteriae)] = "#616161" #
class.labels[c(S.flexneri)] = "#FBC02D" #
class.labels[c(Clade.II)] = "#c5aedb" #
class.labels[c(Clade.III)] = "#7B1FA2" #  
class.labels[c(Clade.IV)] = "#90e585" # 
class.labels[c(Clade.V)] = "#388E3C" #  
class.labels[c(Clade.VI)] = "#feb0d2" # 
class.labels[c(Clade.VII)] = "#9bdeef" # 
class.labels[c(Clade.VIII)] = "#915347" # 

shape.labels = vector("numeric")
shape.labels[c(E.coli)] = 19      # Solid circle
shape.labels[c(E.fergusonii)] = 17 # Triangle point up
shape.labels[c(E.albertii)] = 15   # Filled square
shape.labels[c(S.dysenteriae)] = 3 # Plus
shape.labels[c(S.flexneri)] = 4    # Cross
shape.labels[c(Clade.II)] = 8      # Star
shape.labels[c(Clade.III)] = 9     # Diamond
shape.labels[c(Clade.IV)] = 0      # Square
shape.labels[c(Clade.V)] = 5       # Triangle point down
shape.labels[c(Clade.VI)] = 1      # Circle
shape.labels[c(Clade.VII)] = 2     # Triangle point up filled
shape.labels[c(Clade.VIII)] = 6    # Asterisk

##########################################
#### tsne plot on escherichia dataset ####

tsne <- Rtsne(escherichia_t, dims = 2, distance = 'hamming', check_duplicates = FALSE)

tsne_data <- as.data.frame(tsne$Y)
colnames(tsne_data) <- c('tSNE1', 'tSNE2')


tsne_data$color <- class.labels

# if all reactions -> change to "escherichia_all_rxns_tsne.png"
png(filename = "escherichia_diff_rxns_tsne.png", width = 1200, height = 700)

p <- plot(tsne_data$tSNE1,
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

group_names <- c("E.coli", 
                 "E.fergusonii", 
                 "E.albertii", 
                 "S.dysenteriae", 
                 "S.flexneri", 
                 "Clade II", 
                 "Clade III", 
                 "Clade IV", 
                 "Clade V", 
                 "Clade VI", 
                 "Clade VII", 
                 "Clade VIII" 
)

# Vector with color assignments for each group
class.labels <- c("E.coli" = "#fd7d08", 
                  "E.fergusonii" = "#e46864", 
                  "E.albertii" = "#127db9", 
                  "S.dysenteriae" = "#616161", 
                  "S.flexneri" = "#FBC02D", 
                  "Clade II" = "#c5aedb", 
                  "Clade III" = "#7B1FA2", 
                  "Clade IV" = "#90e585",
                  "Clade V" = "#388E3C",
                  "Clade VI" = "#feb0d2",
                  "Clade VII" = "#9bdeef",
                  "Clade VIII" = "#915347"
)

# Vector with shape assignments for each group
shape.labels <- c("E.coli" = 19, 
                  "E.fergusonii" = 17, 
                  "E.albertii" = 15, 
                  "S.dysenteriae" = 3, 
                  "S.flexneri" = 4, 
                  "Clade II" = 8, 
                  "Clade III" = 9, 
                  "Clade IV" = 0, 
                  "Clade V" = 5,
                  "Clade VI" = 1,
                  "Clade VII" = 2,
                  "Clade VIII" = 6
)


png(filename = "escherichia_legend.png", width = 1000, height = 300)

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

#####################################
### LPCA analysis for Escherichia ###
K = 1

## logistic PCA model
logpca.model = logisticPCA(escherichia_t, # binary data
                           k=K, # number of PCs
                           m=0, # approximation of natural parameter
                           main_effects = TRUE,
                           partial_decomp = TRUE) # including offset term

PC1 <- logpca.model$prop_deviance_expl

K = 2

## logistic PCA model
logpca.model = logisticPCA(escherichia_t, # binary data
                           k=K, # number of PCs
                           m=0, # approximation of natural parameter
                           main_effects = TRUE,
                           partial_decomp = TRUE) # including offset term

# escherichia.model.all if you use all reactions
escherichia.model.diff <- logpca.model
#escherichia.model.all <- logpca.model

escherichia.model <- escherichia.model.diff
#escherichia.model <- escherichia.model.all

PC2 <- escherichia.model$prop_deviance_expl - PC1

logpca.scores = escherichia.model$PCs # extract score matrix
logpca.loadings = escherichia.model$U # extract loading matrix

x <- as.data.frame(logpca.loadings)
x$rxns <- differential.rxns
x$subsystems <- subsystems

# add new subsystem cholesterol degradation
chol_degr <- c("CHOL1",
               "CYP1",
               "CYP2",
               "CYP3",
               "FADD17",
               "FADE281",
               "ECHA191",
               "FADB2",
               "FADE51",
               "CHOLCOA",
               "ECHA192",
               "FADB3",
               "FADE52",
               "FADE282",
               "ECHA193",
               "PRCOA1",
               "KSTD",
               "KSH",
               "3HSA",
               "HSA",
               "HSAC",
               "HDAD",
               "FADD3",
               "HIPB",
               "FADE30",
               "ECHA20",
               "CHCD",
               "HIAA",
               "HIAt"
)

x$subsystems[x$rxns %in% chol_degr] <- "Cholesterol degradation"

x$length <- sqrt(x$V1^2 + x$V2^2)

filtered_groups <- x %>%
  group_by(length) %>%
  filter(n() >= 4)

loading_results <- x %>%
  group_by(subsystems) %>%
  summarize(
    loading_1 = mean(V1),
    loading_2 = mean(V2)
  )

loading_results$length <- sqrt(loading_results$loading_1^2 + loading_results$loading_2^2)

# Sort the dataframe by 'length' in descending order
sorted_result <- loading_results %>% arrange(desc(length))
pathway_importance_lpca <- sorted_result[, c(1,4)]

# Take the top 10 longest vectors
top_longest_vectors <- head(sorted_result, 5)

rem.cluster_ellipse <- ellipse(cov(logpca.scores[rem.clades.cluster,]), center=colMeans(logpca.scores[rem.clades.cluster,]))
ealbertii_ellipse <- ellipse(cov(logpca.scores[E.albertii,]), center=colMeans(logpca.scores[E.albertii,]))
efergusonii_ellipse <- ellipse(cov(logpca.scores[E.fergusonii,]), center=colMeans(logpca.scores[E.fergusonii,]))

png("escherichia_diff_rxns_lpca.png", width = 1200, height = 700)
# change to "escherichia_all_rxns_lpca.png" if all reactions should be considered


plot(logpca.scores, 
     pch=shape.labels,
     col=class.labels, 
     bg=class.labels,
     cex=1.8,
     main="",
     xlab = "PC1 (19%)", 
     ylab = "PC2 (15%)", 
     cex.axis = 1.5,
     cex.lab = 1.5
)

arrows(x0 = rep(0, nrow(top_longest_vectors)), 
       y0 = rep(0, nrow(top_longest_vectors)), 
       x1 = top_longest_vectors$loading_1*2000, 
       y1 = top_longest_vectors$loading_2*2000, 
       col = "black", 
       angle = 25,
       length = 0.1)

for(i in 1:nrow(top_longest_vectors)) {
  if (top_longest_vectors$loading_1[i] < 0) {
    text(x = top_longest_vectors$loading_1[i]*1500, 
         y = top_longest_vectors$loading_2[i]*2500, 
         labels = top_longest_vectors$subsystems[i],
         pos = 2, 
         col = "black",
         cex = 1.3
    )
  }
  if (top_longest_vectors$loading_1[i] > 0) {
    text(x = top_longest_vectors$loading_1[i]*2500, 
         y = top_longest_vectors$loading_2[i]*2500, 
         labels = top_longest_vectors$subsystems[i],
         pos = 4, 
         col = "black",
         cex = 1.3  
    )
  }
}


lines(rem.cluster_ellipse, col="#fd7d08", lwd=2)
lines(ealbertii_ellipse, col="#127db9", lwd=2)
lines(efergusonii_ellipse, col="#e46864", lwd=2)

dev.off()

# select the corresponding subsystem for analysis
selected_rxns <- x %>%
  filter(subsystems == "Murein Biosynthesis")

#####################################################
#### loadings plot differential reactions E.coli ####

escherichia.model <- escherichia.model.diff
#escherichia.model <- escherichia.model.all

PC2 <- escherichia.model$prop_deviance_expl - PC1

logpca.scores = escherichia.model$PCs # extract score matrix
logpca.loadings = escherichia.model$U # extract loading matrix

x <- as.data.frame(logpca.loadings)
x$rxns <- differential.rxns
x$subsystems <- subsystems

x$subsystems[x$rxns %in% chol_degr] <- "Cholesterol degradation"

alaalad = x[202,]
prcoa1 = x[1344,]

png("escherichia_diff_rxns_loadings.png", width = 1200, height = 700)

# Plotting
ggplot(x, aes(x = V1, y = V2, color = subsystems)) + # Assigning aesthetics; x and y coordinates, color by group
  geom_point() + # Use point geom for scatter plot
  geom_text(data = alaalad, aes(label = "ALAALAD", x = V1, y = V2), hjust = 0.5, vjust = -2) +
  geom_text(data = prcoa1, aes(label = "PRCOA1", x = V1, y = V2), hjust = 0.5, vjust = -2) +
  theme_minimal() + # Optional: use a minimal theme for a nice look
  labs(title = "", x = "PC1 (19%)", y = "PC2 (15%)", color = "Subsystem") + # Adding labels
  theme(text = element_text(size = 16), # Increase base text size
        axis.title = element_text(size = 18), # Increase axis titles size
        axis.text = element_text(size = 14))

dev.off()

######################################
#### Jaccard distance and heatmap ####

## define the jaccard similary function ##

calculate_jaccard_similarity <- function(matrix) {
  num_rows <- nrow(matrix)
  jaccard_matrix <- matrix(0, nrow = num_rows, ncol = num_rows)
  
  for (i in 1:num_rows) {
    for (j in i:num_rows) { 
      intersection <- sum(matrix[i,] & matrix[j,])
      union <- sum(matrix[i,] | matrix[j,])
      jaccard_similarity <- intersection / union
      jaccard_matrix[i, j] <- jaccard_similarity
      jaccard_matrix[j, i] <- jaccard_similarity
    }
  }
  return(jaccard_matrix)
}

df <- calculate_jaccard_similarity(escherichia_t)

max_index <- max(c(E.albertii, E.coli, E.fergusonii, S.dysenteriae, S.flexneri, Clade.II, Clade.III, Clade.IV, Clade.V, Clade.VI, Clade.VII, Clade.VIII))
assignments <- vector("character", max_index)

# Function to fill the assignments vector
fill_assignments <- function(index_list, id, assignment_vector) {
  assignment_vector[index_list] <- as.character(id)
  return(assignment_vector)
}

# Fill in the assignments
assignments <- fill_assignments(E.coli, 1, assignments)
assignments <- fill_assignments(E.albertii, 2, assignments)
assignments <- fill_assignments(E.fergusonii, 3, assignments)
assignments <- fill_assignments(S.flexneri, 4, assignments)
assignments <- fill_assignments(S.dysenteriae, 5, assignments)
assignments <- fill_assignments(Clade.II, 6, assignments)
assignments <- fill_assignments(Clade.III, 7, assignments)
assignments <- fill_assignments(Clade.IV, 8, assignments)
assignments <- fill_assignments(Clade.V, 9, assignments)
assignments <- fill_assignments(Clade.VI, 10, assignments)
assignments <- fill_assignments(Clade.VII, 11, assignments)
assignments <- fill_assignments(Clade.VIII, 12, assignments)

group_factors <- factor(assignments, levels = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), labels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"))

group_colors <- c("1" = "#fd7d08", "2" = "#127db9", "3" = "#e46864", "4" = "#FBC02D", "5" = "#616161", "6" = "#c5aedb", "7" = "#7B1FA2", "8" = "#90e585", "9" = "#388E3C", "10" = "#feb0d2", "11" = "#9bdeef", "12" = "#915347")

custom_colors <- colorRampPalette(c("yellow", "lightblue", "blue", "black"))(100)

annotation_col <- HeatmapAnnotation(
  Group = group_factors,
  col = list(Group = group_colors),
  which = "col" 
)

png("escherichia_diff_rxns_jacc.png", width = 1200, height = 1200)
# change the png title to "escherichia_all_rxns_jacc.png" in case the pan-reactions are considered

Heatmap(df, 
        top_annotation = annotation_col, 
        name = "Similarity", 
        column_title = "",
        col = custom_colors,
        column_names_side = "bottom",
        column_names_gp = gpar(col = NA)
)

dev.off()


### legend for Jaccard similarity based heatmap ###

custom_colors <- colorRampPalette(c("white", "lightyellow", "yellow", "lightblue", "blue", "black"))(100)

data <- data.frame(Value = seq(0, 1, by = 0.01))

png("legend_heatmaps.png", width = 1300, height = 200)

ggplot(data, aes(x = Value, y = 1, fill = Value)) +
  geom_tile(stat = "identity") +
  scale_fill_gradientn(colors = custom_colors) +
  theme_minimal() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(size = 70)) +
  labs(fill = "Val")

dev.off()


#########################################################
#####-----subsystem analysis using logistic GLM-----#####

escherichia <- read_excel("data/logPCA/rstb20210236_si_002.xlsx", sheet = "Reaction Info")
clades <- read_excel("data/logPCA/rstb20210236_si_002.xlsx", sheet = "Model reaction counts")

escherichia <- escherichia[, -c(2:5)]

escherichia[is.na(escherichia)] <- 0.0

# once do it with this command for differential reactions and once without for full reactions
escherichia <- escherichia[!(apply(escherichia[,3:224], 1, function(row) all(row == 1) | all(row == 0))), ]

differential.rxns <- escherichia$...1
subsystems <- escherichia$Subsystem

escherichia <- escherichia[, -c(1,2)]

escherichia_t <- as.data.frame(t(escherichia))

escherichia_t$strain_id <- rownames(escherichia_t)

clades <- clades[, -c(2,4,5,6,7)]

merged_df <- merge(escherichia_t, clades, by.x = "strain_id", by.y = "strain_id", all.x = TRUE)

merged_df <- merged_df[,c(-1)]

merged_df$Species <- gsub("Escherichia  fergusonii", "Escherichia fergusonii", merged_df$Species)



# Fit a multinomial logistic regression model
multinom_model <- multinom(Species ~ ., data = merged_df, MaxNWts = 30000)

log_model_smry <- summary(multinom_model)

log_model_coeff <- log_model_smry[["coefficients"]]

log_model_coeff <- log_model_coeff[,-c(1)]

log_model_coeff_t <- t(log_model_coeff)

log_model_coeff_t <- as.data.frame(log_model_coeff_t)

log_model_coeff_t$rxns <- differential.rxns
log_model_coeff_t$subsystems <- subsystems

log_model_coeff_t <- log_model_coeff_t %>%
  rowwise() %>%
  mutate(squared_sum_importance = sqrt(sum(c_across(1:12)^2)))

pathway_importance_mlr_2 <- log_model_coeff_t %>%
  group_by(subsystems) %>%
  summarise(total_sqrt_sum_importance = mean(squared_sum_importance, na.rm = TRUE))

#### Compare RF, MLR and LPCA subsystem enrichment ####

merged_rankings <- merge(pathway_importance_mlr_2, pathway_importance_lpca, by = "subsystems")

colnames(merged_rankings)[2] = "MLR"
colnames(merged_rankings)[3] = "LPCA"

## heatmap

merged_rankings$MLR <- merged_rankings$MLR / max(merged_rankings$MLR)
merged_rankings$LPCA <- merged_rankings$LPCA / max(merged_rankings$LPCA)

heatmap_data <- as.matrix(merged_rankings[, 2:3])
rownames(heatmap_data) <- merged_rankings$subsystems


color_palette <- colorRampPalette(c("white", "yellow", "blue", "darkblue", "black"))(100)

png("comparison_lpca_mlr_2.png", width = 1500, height = 2200)

heatmap.2(heatmap_data, scale = "none", trace = "none", col = color_palette, 
          margin = c(15, 40), cexRow = 2, cexCol = 2.5,
          key = FALSE, density.info = "none",
          main = "", xlab = "", ylab = "",
          Colv = FALSE)

dev.off()


##################################
### legend for mlr & lpca plot ###

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


#####################################
### PCA on simulated growth rates ###

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

shape.labels = vector("numeric")
shape.labels[c(1:87, 93, 94)] = 19      # Solid circle
shape.labels[c(136:153)] = 17 # Triangle point up
shape.labels[c(154:207)] = 15   # Filled square
shape.labels[c(88:89)] = 3 # Plus
shape.labels[c(90:92)] = 4    # Cross
shape.labels[c(216:219)] = 8      # Star
shape.labels[c(128:135)] = 9     # Diamond
shape.labels[c(208:215)] = 0      # Square
shape.labels[c(95:135)] = 5       # Triangle point down
shape.labels[c(221)] = 1      # Circle
shape.labels[c(220)] = 2     # Triangle point up filled
shape.labels[c(222)] = 6    # Asterisk


png("escherichia_sim_growth_rates_pca.png", width = 1200, height = 700)

plot(data.pca$scores[, 1:2],  # x and y data
     pch=shape.labels,           # point shape
     col=class.labels, # 
     bg=class.labels,  #
     cex=1.4,          # point size
     main="",     # title of plot
     xlab = "PC1 (49%)", # 
     ylab = "PC2 (20%)",
     cex.axis = 1.5,
     cex.lab = 1.5
)

dev.off()

#####################################
### Phylogenetic tree Escherichia ###


## Whole genome based phylogeny

tree <- read.tree("/home/users/lzehetner/data/logPCA/SpeciesTree_rooted.txt")

data <- read_excel("data/rstb20210236_si_003.xlsx", sheet = "growth predictions")

data <- data[-223,]

# Split the data into a list by species
split_data <- split(data$strain_id, data$Species)

# Convert list elements to vectors and name them
species_vectors <- lapply(split_data, as.vector)

e.coli <- species_vectors[["Escherichia coli"]]
e.albertii <- species_vectors[["Escherichia albertii"]]
e.fergusonii <- species_vectors[["Escherichia fergusonii"]]
s.flexneri <- species_vectors[["Shigella flexneri"]]
s.dysenteriae <- species_vectors[["Shigella dysenteriae"]]
clade.II <- species_vectors[["Clade II"]]
clade.III <- species_vectors[["Clade III"]]
clade.IV <- species_vectors[["Clade IV"]]
clade.V <- species_vectors[["Clade V"]]
clade.VI <- species_vectors[["Clade VI"]]
clade.VII <- species_vectors[["Clade VII"]]
clade.VIII <- species_vectors[["Clade VIII"]]


groups <- list(e.coli,
               e.albertii,
               e.fergusonii,
               s.flexneri,
               s.dysenteriae,
               clade.II,
               clade.III,
               clade.IV,
               clade.V,
               clade.VI,
               clade.VII,
               clade.VIII
)


colors <- c("#fd7d08",
            "#127db9" ,
            "#e46864",
            "#ffe4ca",
            "#feb166",
            "#c5aedb",
            "#2c9f28",
            "#90e585",
            "#abcaff",
            "#feb0d2",
            "#9bdeef",
            "#915347"
)

assign_edge_colors <- function(tree, groups, colors) {
  edge_colors <- rep("grey", length(tree$edge[,1])) # Default color for edges
  for (i in seq_along(groups)) {
    for (species in groups[[i]]) {
      tip_index <- which(tree$tip.label == species)
      if (length(tip_index) > 0) {
        # Find the edge that corresponds to this tip
        edge_index <- which(tree$edge[,2] == tip_index)
        edge_colors[edge_index] <- colors[i]
      }
    }
  }
  return(edge_colors)
}

edge_colors <- assign_edge_colors(tree, groups, colors)


plot.phylo(tree, type = "phylogram", show.tip.label = FALSE, edge.color = edge_colors)

legend_labels <- c("E. coli", "E. albertii", "E. fergusonii", 
                   "S. flexneri", "S. dysenteriae", "Clade II", 
                   "Clade III", "Clade IV", "Clade V", 
                   "Clade VI", "Clade VII", "Clade VIII")

# Add the legend to the plot
legend("topright", # position of the legend
       legend = legend_labels, 
       col = colors, 
       lty = 1, # line type
       cex = 1.1) # font size of the legend text


##################################################################
###### Cophenetic correlation of all escherichias ################

tree <- read.tree("/home/users/lzehetner/data/logPCA/SpeciesTree_rooted.txt")
tree_ultrametric <- chronos(tree, lambda = 0.00000001)
tree_binary <- multi2di(tree_ultrametric)
dend <- as.dendrogram(tree_binary)
dend_cophenetic <- cophenetic(dend)

distance_matrix_logpca <- dist(logpca.scores, method = "euclidian")
hc <- hclust(distance_matrix_logpca)
binary_dend <- as.dendrogram(hc)
binary_cophenetic_logpca <- cophenetic(binary_dend)

correlation_coefficient <- cor(dend_cophenetic, binary_cophenetic_logpca)
print(correlation_coefficient)


## differential reactions
# coph corr. LPCA vs Phylo (manhattan): 0.61  (euclidian): 0.37

## all reactions
# coph corr. LPCA vs Phylo (manhattan): 0.63  (euclidian): 0.36

