## Script to analyze the Firmicutes from the Agora2 dataset derived from Heineken (2023)
## This script includes data preprocessing, LPCA, t-sne, Jaccard similarity, phylogeny and cophenetic correlation

library(logisticPCA) # for logistic PCA ans logistic SVD
library(readr)
library(rARPACK)
library(readxl)
library(ggplot2)
library(reshape2)
library(ComplexHeatmap)
library(Rtsne)
library(dplyr)

set.seed(42)

##################
### Firmicutes ###

# import binary reaction matrix

bacteria <- read_csv("/home/users/lzehetner/data/agora2/reactions_presence_matrix.csv")

bacteria <- bacteria[!(apply(bacteria[,2:5557], 1, function(row) all(row == 1) | all(row == 0))), ]
reactions <- bacteria[, 1]
colnames(reactions)[1] <- "rxns"

bacteria <- bacteria[, -c(1)]

# import metadata from Agora2

df <- read_excel("/home/users/lzehetner/data/agora2/41587_2022_1628_MOESM3_ESM.xlsx", sheet = "Table_S1")
df <- df[-c(1,2),]

df <- df[, c(1,8,9,10)]

names(df)[1] <- "species"

filtered_df <- df[df$species %in% colnames(bacteria), ]

df2 <- filtered_df %>% 
  filter(filtered_df$...10 == "Firmicutes")

column_names_to_select <- df2[[1]]

# Step 2: Subset df2 by these column names
df3 <- bacteria[, column_names_to_select, drop = FALSE]

df3$rxns <- reactions$rxns

df3 <- df3[!(apply(df3[,1:2943], 1, function(row) all(row == 1) | all(row == 0))), ]
differential.rxns.fir <- df3$rxns

firmicutes <- df3[, -c(2944)]
firmicutes_t <- t(firmicutes)


species_names <- colnames(firmicutes)

by_class <- split(df2$species, df2$...8)

index_lists <- lapply(by_class, function(species_names) {
  match(species_names, colnames(firmicutes))
})

Acidaminococcales <- index_lists$Acidaminococcales
Bacillales <- index_lists$Bacillales
C.Borkfalkiales <- index_lists$`Candidatus Borkfalkiales`
Clostridiales <- index_lists$Clostridiales
Erysipelotrichales <- index_lists$Erysipelotrichales
Eubacteriales <- index_lists$Eubacteriales
Lactobacillales <- index_lists$Lactobacillales
Selenomonadales <- index_lists$Selenomonadales
Thermoanaerobacterales <- index_lists$Thermoanaerobacterales
Tissierellales <- index_lists$Tissierellales
Tissierellia <- index_lists$Tissierellia
U.Tissierellia <- index_lists$`unclassified Tissierellia`
Veillonellales <- index_lists$Veillonellales

class.labels = vector()
class.labels[c(Acidaminococcales)] = "#d32f2f" # red
class.labels[c(Bacillales)] = "#fbc02d" # yellow
class.labels[c(C.Borkfalkiales)] = "#7b1fa2" # purple
class.labels[c(Clostridiales)] = "#5d4037" # brown
class.labels[c(Erysipelotrichales)] = "#303f9f" # indigo
class.labels[c(Eubacteriales)] = "#1976d2" # blue
class.labels[c(Lactobacillales)] = "#616161" # grey
class.labels[c(Selenomonadales)] = "#0097a7" # cyan
class.labels[c(Thermoanaerobacterales)] = "#f57c00" # orange
class.labels[c(Tissierellales)] = "#388e3c" # green
class.labels[c(Tissierellia)] = "#000000" # black
class.labels[c(U.Tissierellia)] = "darkgreen" # 
class.labels[c(Veillonellales)] = "pink" # pink

shape.labels = vector("numeric")
shape.labels[c(Acidaminococcales)] = 9      # Solid circle
shape.labels[c(Bacillales)] = 17 # Triangle point up
shape.labels[c(C.Borkfalkiales)] = 15   # Filled square
shape.labels[c(Clostridiales)] = 3 # Plus
shape.labels[c(Erysipelotrichales)] = 4    # Cross
shape.labels[c(Eubacteriales)] = 8      # Star
shape.labels[c(Lactobacillales)] = 19     # Diamond
shape.labels[c(Selenomonadales)] = 0      # Square
shape.labels[c(Thermoanaerobacterales)] = 5       # Triangle point down
shape.labels[c(Tissierellales)] = 1      # Circle
shape.labels[c(Tissierellia)] = 2     # Triangle point up filled
shape.labels[c(U.Tissierellia)] = 6    # Asterisk
shape.labels[c(Veillonellales)] = 7 

## LPCA
K = 1

## logistic PCA model
logpca.model = logisticPCA(firmicutes_t, # binary data
                           k=K, # number of PCs
                           m=0, # approximation of natural parameter
                           main_effects = TRUE,
                           partial_decomp = TRUE) # including offset term

PC1.fir <- logpca.model$prop_deviance_expl

K = 2

## logistic PCA model
logpca.model = logisticPCA(firmicutes_t, # binary data
                           k=K, # number of PCs
                           m=0, # approximation of natural parameter
                           main_effects = TRUE,
                           partial_decomp = TRUE) # including offset term

fir.model <- logpca.model

logpca.scores <- fir.model$PCs
logpca.loadings <- fir.model$U

y <- as.data.frame(logpca.scores)

x <- as.data.frame(logpca.loadings)
x$rxns <- differential.rxns.fir

x$length <- sqrt(x$V1^2 + x$V2^2)

# Sort the dataframe by 'length' in descending order
sorted_result <- x %>% arrange(desc(length))

# Take the top 10 longest vectors
top_longest_vectors <- head(sorted_result, 5)

filtered_groups <- x %>%
  group_by(length) %>%
  filter(n() >= 4)


png("firmicutes_lpca.png", width = 1200, height = 700)

plot(logpca.scores,  # x and y data
     pch=shape.labels,           # point shape
     col=class.labels, # 
     bg=class.labels,  #
     cex=1.5,          # point size
     main="",     # title of plot
     xlab = "PC1 (19%)", # for "all rxns": PC1 = 19%, PC2 = 15%
     ylab = "PC2 (10%)", # for "differential rxns": PC1 = 19%, PC2 = 15%
     cex.axis = 1.5,
     cex.lab = 1.5
)

dev.off()

selected_rxns <- x %>%
  filter(subsystems == "Pyruvate Metabolism")

########################################
### Jaccard similarity on Firmicutes ###

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

df <- calculate_jaccard_similarity(firmicutes_t)

max_index <- max(c(Acidaminococcales, Bacillales, C.Borkfalkiales,Clostridiales,Erysipelotrichales,Eubacteriales,
                   Lactobacillales,Selenomonadales,Thermoanaerobacterales,Tissierellales,Tissierellia,U.Tissierellia,Veillonellales))
assignments <- vector("character", max_index)

# Function to fill the assignments vector
fill_assignments <- function(index_list, id, assignment_vector) {
  assignment_vector[index_list] <- as.character(id)
  return(assignment_vector)
}

# Fill in the assignments
assignments <- fill_assignments(Acidaminococcales, 1, assignments)
assignments <- fill_assignments(Bacillales, 2, assignments)
assignments <- fill_assignments(C.Borkfalkiales, 3, assignments)
assignments <- fill_assignments(Clostridiales, 4, assignments)
assignments <- fill_assignments(Erysipelotrichales, 5, assignments)
assignments <- fill_assignments(Eubacteriales, 6, assignments)
assignments <- fill_assignments(Lactobacillales, 7, assignments)
assignments <- fill_assignments(Selenomonadales, 8, assignments)
assignments <- fill_assignments(Thermoanaerobacterales, 9, assignments)
assignments <- fill_assignments(Tissierellales, 10, assignments)
assignments <- fill_assignments(Tissierellia, 11, assignments)
assignments <- fill_assignments(U.Tissierellia, 12, assignments)
assignments <- fill_assignments(Veillonellales, 13, assignments)


group_factors <- factor(assignments, levels = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13), 
                        labels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13"))

group_colors <- c("1" = "#d32f2f", "2" = "#fbc02d", "3" = "#7b1fa2", "4" = "#5d4037", "5" = "#303f9f", 
                  "6" = "#1976d2", "7" = "#616161", "8" = "#0097a7", "9" = "#f57c00", "10" = "#388e3c", 
                  "11" = "#000000", "12" = "darkgreen", "13" = "pink")

custom_colors <- colorRampPalette(c("yellow", "lightblue", "blue", "black"))(100)

annotation_col <- HeatmapAnnotation(
  Group = group_factors,
  col = list(Group = group_colors),
  which = "col" 
)

png("firmicutes_jacc.png", width = 2000, height = 2000)

Heatmap(df, 
        top_annotation = annotation_col, 
        name = "Similarity", 
        column_title = "",
        col = custom_colors,
        column_names_side = "bottom",
        column_names_gp = gpar(col = NA)
)

dev.off()

############################################
### t-SNE analysis on Firmicutes dataset ###

tsne <- Rtsne(firmicutes_t, dims = 2, distance = 'hamming', check_duplicates = FALSE)

tsne_data <- as.data.frame(tsne$Y)
colnames(tsne_data) <- c('tSNE1', 'tSNE2')


tsne_data$color <- class.labels

# if all reactions -> change to "escherichia_all_rxns_tsne.png"
png(filename = "firmicutes_tsne.png", width = 1200, height = 700)

p <- plot(tsne_data$tSNE1,
          tsne_data$tSNE2,
          pch=shape.labels,           # point shape
          col=class.labels, # 
          bg=class.labels,  #
          cex=1.5,          # point size
          main="",     # title of plot
          xlab = "tSNE1", # 
          ylab = "tSNE2", # 
          cex.axis = 1.3,
          cex.lab = 1.3
)

#p <- colors <- c("#d32f2f", "#fbc02d", "#7b1fa2", "#5d4037", "#303f9f", "#1976d2", 
#                 "#616161", "#0097a7", "#f57c00", "#388e3c", "#000000", "darkgreen", "pink")
#p <- legend("topleft", 
#            legend = c("Acidaminococcales", "Bacillales", "C.Borkfalkiales", "Clostridiales", 
#                       "Erysipelotrichales", "Eubacteriales", "Lactobacillales", "Selenomonadales", 
#                       "Thermoanaerobacterales", "Tissierellales", "Tissierellia", "U.Tissierellia", "Veillonellales"),
#            fill=colors, 
#            cex=1.3)

dev.off()

########################################
###------------ legend --------------###

group_names <- c("Acidaminococcales", 
                 "Bacillales", 
                 "C.Borkfalkiales", 
                 "Clostridiales", 
                 "Erysipelotrichales", 
                 "Eubacteriales", 
                 "Lactobacillales", 
                 "Selenomonadales", 
                 "Thermoanaerobacterales", 
                 "Tissierellales", 
                 "Tissierellia", 
                 "U.Tissierellia",
                 "Veillonellales"
)

# Vector with color assignments for each group
class.labels <- c("Acidaminococcales" = "#d32f2f",
                  "Bacillales" = "#fbc02d",
                  "C.Borkfalkiales" = "#7b1fa2",
                  "Clostridiales" = "#5d4037",
                  "Erysipelotrichales" = "#303f9f",
                  "Eubacteriales" = "#1976d2",
                  "Lactobacillales" = "#616161",
                  "Selenomonadales" = "#0097a7",
                  "Thermoanaerobacterales" = "#f57c00",
                  "Tissierellales" = "#388e3c",
                  "Tissierellia" = "#000000" ,
                  "U.Tissierellia" = "darkgreen",
                  "Veillonellales" = "pink"
)

# Vector with shape assignments for each group
shape.labels <- c("Acidaminococcales" = 9,
                  "Bacillales" = 17,
                  "C.Borkfalkiales" = 15,
                  "Clostridiales" = 3,
                  "Erysipelotrichales" = 4,
                  "Eubacteriales" = 8,
                  "Lactobacillales" = 19,
                  "Selenomonadales" = 0,
                  "Thermoanaerobacterales" = 5,
                  "Tissierellales" = 1,
                  "Tissierellia" = 2,
                  "U.Tissierellia" = 6,
                  "Veillonellales" = 7
)

png(filename = "firmicutes_legend.png", width = 1000, height = 300)

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
