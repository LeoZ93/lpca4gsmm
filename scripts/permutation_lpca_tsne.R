library(Rtsne)
library(ggplot2)
library(dplyr)

### Data ###

escherichia <- read_excel("data/logPCA/rstb20210236_si_002.xlsx", sheet = "Reaction Info")
clades <- read_excel("data/logPCA/rstb20210236_si_002.xlsx", sheet = "Model reaction counts")

escherichia <- escherichia[, -c(2:5)]

escherichia[is.na(escherichia)] <- 0.0

# once do it with this command for differential reactions and once without for full reactions
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


########## Logistic PCA #############

library(logisticPCA)
library(ggplot2)
library(dplyr)


# Function to perform logistic PCA and collect results
run_logisticPCA <- function(data) {
  logpca.model <- logisticPCA(data, k = 2, m = 0, main_effects = TRUE, partial_decomp = TRUE)
  scores <- logpca.model$PCs  # Extracting the first two principal components
  return(data.frame(PC1 = scores[,1], PC2 = scores[,2]))
}

# Permute and run logistic PCA 10 times
K <- 2  # Number of principal components
results <- lapply(1:10, function(i) {
  permuted_data <- escherichia_t[sample(nrow(escherichia_t)), ]
  pca_results <- run_logisticPCA(permuted_data)
  pca_results$run <- i
  return(pca_results)
})

# Combine results and calculate statistics
all_results <- do.call(rbind, results)
stats <- all_results %>%
  group_by(run) %>%
  summarise(mean_PC1 = mean(PC1), mean_PC2 = mean(PC2),
            sd_PC1 = sd(PC1), sd_PC2 = sd(PC2))

png(filename = "lpca_permutation.png", width = 1200, height = 700)

# Plotting the results
ggplot(all_results, aes(x = PC1, y = PC2, color = factor(run))) +
  geom_point(alpha = 0.4) +
  theme_minimal() +
  labs(title = "", color = "Permutation", x = "PC1", y = "PC2") +
  theme(plot.title = element_text(size = 15),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 15),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 15) )

dev.off()


#print(p1)

print(stats)



########################################

set.seed(42)

# Function to perform t-SNE and collect results
run_tsne <- function(input_data) {  # Use a different name for the function argument
  tsne <- Rtsne(input_data, dims = 2, perplexity = 30, verbose = FALSE)
  return(data.frame(tsne1 = tsne$Y[,1], tsne2 = tsne$Y[,2]))
}

# Permute and run t-SNE 10 times
results <- lapply(1:10, function(i) {
  permuted_data <- escherichia_t[sample(nrow(escherichia_t)), ]  # Refer to the global data matrix
  tsne_results <- run_tsne(permuted_data)
  tsne_results$run <- i
  return(tsne_results)
})

# Combine results and calculate statistics
all_results_tsne <- do.call(rbind, results)
stats <- all_results_tsne %>%
  group_by(run) %>%
  summarise(mean_tsne1 = mean(tsne1), mean_tsne2 = mean(tsne2),
            sd_tsne1 = sd(tsne1), sd_tsne2 = sd(tsne2))


png(filename = "tsne_permutation.png", width = 1200, height = 700)
# Plotting the results
ggplot(all_results_tsne, aes(x = tsne1, y = tsne2, color = factor(run))) +
  geom_point(alpha = 0.4) +
  theme_minimal() +
  labs(title = "", color = "Permutation", x = "tSNE1", y = "tSNE2") +
  theme(plot.title = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15) )

dev.off()

print(stats)


