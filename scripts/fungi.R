## Script to analyze the Fungi dataset derived from Li (2021)
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


fungi_df <- read_csv("/home/users/lzehetner/data/logPCA/differential_rxns_per_fungus.csv", col_names = TRUE)

fungi <- fungi_df[, -c(1,2,3)]

fungi_t <- t(fungi)

## fungi classes ##
alloascoideaceae = c('Alloascoidea_hylecoeti')
ascomycota = c('Arthrobotrys_oligospora', 'Aspergillus_nidulans', 'Botrytis_cinerea', 'Coccidioides_immitis', 'Fusarium_graminearum', 'Neurospora_crassa', 'Saitoella_complicata', 'Schizosaccharomyces_pombe', 'Sclerotinia_sclerotiorum', 'Stagonospora_nodorum', 'Xylona_heveae')
cug_ala = c('Nakazawaea_holstii', 'Nakazawaea_peltata', 'Pachysolen_tannophilus', 'Peterozyma_toletana', 'Peterozyma_xylosa')
cug_ser1 = c("Aciculoconidium_aculeatum", "Babjeviella_inositovora", "Candida_albicans", "Candida_ascalaphidarum", "Candida_athensensis", "Candida_auris", "Candida_blattae", "Candida_carpophila", "Candida_corydali", "Candida_dubliniensis", "Candida_fragi", "Candida_fructus", "Candida_golubevii", "Candida_gorgasii", "Candida_gotoi", "Candida_hawaiiana", "Candida_heveicola", "Candida_intermedia", "Candida_oregonensis", "Candida_orthopsilosis", "Candida_parapsilosis", "Candida_restingae", "Candida_rhagii", "Candida_schatavii", "Candida_sojae", "Candida_tammaniensis", "Candida_wancherniae", "Cephaloascus_albidus", "Cephaloascus_fragrans", "Clavispora_lusitaniae", "Danielozyma_ontarioensis", "Debaryomyces_fabryi", "Debaryomyces_hansenii", "Debaryomyces_maramus", "Debaryomyces_nepalensis", "Debaryomyces_prosopidis", "Debaryomyces_subglobosus", "Hyphopichia_burtonii", "Hyphopichia_heimii", "Hyphopichia_homilentoma", "Kodamaea_laetipori", "Kodamaea_ohmeri", "Kurtzmaniella_cleridarum", "Lodderomyces_elongisporus", "Metschnikowia_aberdeeniae", "Metschnikowia_arizonensis", "Metschnikowia_bicuspidata_var._bicuspidata", "Metschnikowia_borealis", "Metschnikowia_bowlesiae", "Metschnikowia_cerradonensis", "Metschnikowia_continentalis", "Metschnikowia_dekortorum", "Metschnikowia_drakensbergensis", "Metschnikowia_hamakuensis", "Metschnikowia_hawaiiensis", "Metschnikowia_hibisci", "Metschnikowia_ipomoeae", "Metschnikowia_kamakouana", "Metschnikowia_kipukae", "Metschnikowia_lochheadii", "Metschnikowia_matae_var._maris", "Metschnikowia_matae_var._matae", "Metschnikowia_mauinuiana", "Metschnikowia_proteae", "Metschnikowia_santaceciliae", "Metschnikowia_shivogae", "Metschnikowia_similis", "Meyerozyma_caribbica", "Meyerozyma_guilliermondii", "Millerozyma_acaciae", "Priceomyces_carsonii", "Priceomyces_castillae", "Priceomyces_haplophilus", "Priceomyces_medius", "Scheffersomyces_lignosus", "Scheffersomyces_stipitis", "Spathaspora_arborariae", "Spathaspora_girioi", "Spathaspora_gorwiae", "Spathaspora_hagerdaliae", "Spathaspora_passalidarum", "Suhomyces_canberraensis", "Suhomyces_emberorum", "Suhomyces_pyralidae", "Suhomyces_tanzawaensis", "Teunomyces_cretensis", "Teunomyces_gatunensis", "Teunomyces_kruisii", "Wickerhamia_fluorescens", "Yamadazyma_nakazawae", "Yamadazyma_philogaea", "Yamadazyma_scolyti", "Yamadazyma_tenuis")
cug_ser2 = c('Ascoidea_asiatica', 'Ascoidea_rubescens', 'Saccharomycopsis_capsularis', 'Saccharomycopsis_malanga')
dipodascaceae = c('Blastobotrys_adeninivorans', 'Blastobotrys_americana', 'Blastobotrys_mokoenaii', 'Blastobotrys_muscicola', 'Blastobotrys_nivea', 'Blastobotrys_peoriensis', 'Blastobotrys_proliferans', 'Blastobotrys_raffinosifermentans', 'Blastobotrys_serpentis', 'Candida_hispaniensis', 'Candida_incommunis', 'Deakozyma_indianensis', 'Diddensiella_caesifluorescens', 'Dipodascus_albidus', 'Dipodascus_geniculatus', 'Geotrichum_candidum', 'Groenewaldozyma_salmanticensis', 'Magnusiomyces_tetraspermus', 'Middelhovenomyces_tepae', 'Nadsonia_fulvescens_var._elongata', 'Nadsonia_fulvescens_var._fulvescens', 'Saprochaete_clavata', 'Spencermartinsiella_europaea', 'Starmerella_apicola', 'Starmerella_bombicola', 'Sugiyamaella_lignohabitans', 'Wickerhamiella_cacticola', 'Wickerhamiella_domercqiae', 'Wickerhamiella_infanticola', 'Wickerhamiella_versatilis', 'Yarrowia_bubula', 'Yarrowia_deformans', 'Yarrowia_divulgata', 'Yarrowia_keelungensis', 'Yarrowia_lipolytica', 'Zygoascus_meyerae', 'Zygoascus_ofunaensis')
lipomycetaceae = c('Lipomyces_arxii', 'Lipomyces_doorenjongii', 'Lipomyces_japonicus', 'Lipomyces_kononenkoae', 'Lipomyces_lipofer', 'Lipomyces_mesembrius', 'Lipomyces_oligophaga', 'Lipomyces_starkeyi', 'Lipomyces_suomiensis')
phaffomycetaceae = c('Barnettozyma_californica', 'Barnettozyma_hawaiiensis', 'Barnettozyma_populi', 'Barnettozyma_pratensis', 'Barnettozyma_salicaria', 'Candida_freyschussii', 'Candida_montana', 'Candida_mycetangii', 'Candida_orba', 'Candida_ponderosae', 'Candida_stellimalicola', 'Candida_vartiovaarae', 'Cyberlindnera_americana', 'Cyberlindnera_fabianii', 'Cyberlindnera_jadinii', 'Cyberlindnera_maclurae', 'Cyberlindnera_misumaiensis', 'Cyberlindnera_mrakii', 'Cyberlindnera_petersonii', 'Cyberlindnera_saturnus', 'Cyberlindnera_suaveolens', 'Cyberlindnera_xylosilytica', 'Phaffomyces_antillensis', 'Phaffomyces_opuntiae', 'Phaffomyces_thermotolerans', 'Starmera_amethionina', 'Starmera_quercuum', 'Wickerhamomyces_alni', 'Wickerhamomyces_anomalus', 'Wickerhamomyces_bovis', 'Wickerhamomyces_canadensis', 'Wickerhamomyces_ciferrii', 'Wickerhamomyces_hampshirensis', 'Wickerhamomyces_sp.')
pichiaceae = c('Ambrosiozyma_ambrosiae', 'Ambrosiozyma_kashinagacola', 'Ambrosiozyma_maleeae', 'Ambrosiozyma_monospora', 'Ambrosiozyma_oregonensis', 'Ambrosiozyma_philentoma', 'Ambrosiozyma_pseudovanderkliftii', 'Ambrosiozyma_vanderkliftii', 'Brettanomyces_anomalus', 'Brettanomyces_bruxellensis', 'Brettanomyces_custersianus', 'Candida_arabinofermentans', 'Candida_boidinii', 'Candida_sorboxylosa', 'Candida_succiphila', 'Citeromyces_hawaiiensis', 'Citeromyces_matritensis', 'Citeromyces_siamensis', 'Komagataella_pastoris', 'Komagataella_populi', 'Komagataella_pseudopastoris', 'Kregervanrija_delftensis', 'Kregervanrija_fluxuum', 'Kuraishia_capsulata', 'Kuraishia_molischiana', 'Kuraishia_ogatae', 'Martiniozyma_abiesophila', 'Ogataea_glucozyma', 'Ogataea_henricii', 'Ogataea_kodamae', 'Ogataea_methanolica', 'Ogataea_methylivora', 'Ogataea_minuta', 'Ogataea_naganishii', 'Ogataea_nitratoaversa', 'Ogataea_nonfermentans', 'Ogataea_parapolymorpha', 'Ogataea_philodendri', 'Ogataea_pilisensis', 'Ogataea_pini', 'Ogataea_polymorpha', 'Ogataea_populialbae', 'Ogataea_ramenticola', 'Ogataea_trehaloabstinens', 'Ogataea_trehalophila', 'Ogataea_zsoltii', 'Pichia_exigua', 'Pichia_heedii', 'Pichia_kudriavzevii', 'Pichia_membranifaciens', 'Pichia_nakasei', 'Pichia_norvegensis', 'Pichia_occidentalis', 'Pichia_terricola', 'Saturnispora_dispora', 'Saturnispora_hagleri', 'Saturnispora_mendoncae', 'Saturnispora_saitoi', 'Saturnispora_serradocipensis', 'Saturnispora_silvae', 'Saturnispora_zaruensis')
saccharomycetaceae = c('Ashbya_aceri', 'Candida_bracarensis', 'Candida_castellii', 'Candida_glabrata', 'Candida_nivariensis', 'Eremothecium_coryli', 'Eremothecium_cymbalariae', 'Eremothecium_gossypii', 'Eremothecium_sinecaudum', 'Kazachstania_aerobia', 'Kazachstania_africana', 'Kazachstania_bromeliacearum', 'Kazachstania_intestinalis', 'Kazachstania_kunashirensis', 'Kazachstania_martiniae', 'Kazachstania_naganishii', 'Kazachstania_rosinii', 'Kazachstania_siamensis', 'Kazachstania_solicola', 'Kazachstania_spencerorum', 'Kazachstania_taianensis', 'Kazachstania_transvaalensis', 'Kazachstania_turicensis', 'Kazachstania_unispora', 'Kazachstania_viticola', 'Kazachstania_yakushimaensis', 'Kluyveromyces_aestuarii', 'Kluyveromyces_dobzhanskii', 'Kluyveromyces_lactis', 'Kluyveromyces_marxianus', 'Lachancea_cidri', 'Lachancea_dasiensis', 'Lachancea_fantastica_nom._nud.', 'Lachancea_fermentati', 'Lachancea_kluyveri', 'Lachancea_lanzarotensis', 'Lachancea_meyersii', 'Lachancea_mirantina', 'Lachancea_nothofagi', 'Lachancea_quebecensis', 'Lachancea_thermotolerans', 'Lachancea_waltii', 'Nakaseomyces_bacillisporus', 'Nakaseomyces_delphensis', 'Naumovozyma_castellii', 'Naumovozyma_dairenensis', 'Saccharomyces_arboricola', 'Saccharomyces_cerevisiae', 'Saccharomyces_eubayanus', 'Saccharomyces_kudriavzevii', 'Saccharomyces_mikatae', 'Saccharomyces_paradoxus', 'Saccharomyces_uvarum', 'Tetrapisispora_blattae', 'Tetrapisispora_fleetii', 'Tetrapisispora_iriomotensis', 'Tetrapisispora_namnaonensis', 'Tetrapisispora_phaffii', 'Torulaspora_delbrueckii', 'Torulaspora_franciscae', 'Torulaspora_maleeae', 'Torulaspora_microellipsoides', 'Torulaspora_pretoriensis', 'Vanderwaltozyma_polyspora', 'Yueomyces_sinensis', 'Zygosaccharomyces_bailii', 'Zygosaccharomyces_bisporus', 'Zygosaccharomyces_kombuchaensis', 'Zygosaccharomyces_rouxii', 'Zygotorulaspora_florentina', 'Zygotorulaspora_mrakii')
saccharomycodaceae = c('Hanseniaspora_clermontiae', 'Hanseniaspora_osmophila', 'Hanseniaspora_pseudoguilliermondii', 'Hanseniaspora_singularis', 'Hanseniaspora_uvarum', 'Hanseniaspora_valbyensis', 'Hanseniaspora_vineae', 'Kloeckera_hatyaiensis')
sporopachydermia = c('Sporopachydermia_lactativora', 'Sporopachydermia_quercuum')
trigonopsidaceae = c('Botryozyma_nematodophila', 'Tortispora_caseinolytica', 'Tortispora_ganteri', 'Tortispora_starmeri', 'Trigonopsis_variabilis', 'Trigonopsis_vinaria')

index_1 <- match(alloascoideaceae, names(fungi))
index_2 <- match(ascomycota, names(fungi))
index_3 <- match(cug_ala, names(fungi))
index_4 <- match(cug_ser1, names(fungi))
index_5 <- match(cug_ser2, names(fungi))
index_6 <- match(dipodascaceae, names(fungi))
index_7 <- match(lipomycetaceae, names(fungi))
index_8 <- match(phaffomycetaceae, names(fungi))
index_9 <- match(pichiaceae, names(fungi))
index_10 <- match(saccharomycetaceae, names(fungi))
index_11 <- match(saccharomycodaceae, names(fungi))
index_12 <- match(sporopachydermia, names(fungi))
index_13 <- match(trigonopsidaceae, names(fungi))

# for cophenetic correlation, ascomycota need to be removed, since they were not part of the original analysis by Shen.
#fungi_t <- fungi_t[!(row.names(fungi_t) %in% ascomycota), ]

class.labels = vector()
class.labels[c(index_1)] = "#559668" # alloascoideaceae
class.labels[c(index_2)] = "#000000" # outgrouped
class.labels[c(index_3)] = "#633717" # cug_ala
class.labels[c(index_4)] = "#fac34e" # cug_ser1
class.labels[c(index_5)] = "#a7cd57" # cug_ser2
class.labels[c(index_6)] = "#f83e32" # dipodascaceae
class.labels[c(index_7)] = "#942d8c" # lipomycetaceae
class.labels[c(index_8)] = "#2cc5da" # phaffomycetaceae
class.labels[c(index_9)] = "#fc803e" # pichiaceae
class.labels[c(index_10)] = "#0a5aa2" # saccharomycetaceae
class.labels[c(index_11)] = "#04988d" # saccharomycodaceae
class.labels[c(index_12)] = "#58539e" # sporopachydermia
class.labels[c(index_13)] = "#f5348f" # trigonopsidaceae

shape.labels = vector("numeric")
shape.labels[c(index_1)] = 4 # plus
shape.labels[c(index_2)] = 17 # Triangle point up
shape.labels[c(index_3)] = 5   # Filled square
shape.labels[c(index_4)] = 19 # solid circle
shape.labels[c(index_5)] = 4    # Cross
shape.labels[c(index_6)] = 8      # Star
shape.labels[c(index_7)] = 9     # Diamond
shape.labels[c(index_8)] = 0      # Square
shape.labels[c(index_9)] = 15       # Triangle point down
shape.labels[c(index_10)] = 2      # Circle
shape.labels[c(index_11)] = 1     # Triangle point up filled
shape.labels[c(index_12)] = 6    # Asterisk
shape.labels[c(index_13)] = 7  # 

K = 1

## logistic PCA model
logpca.model = logisticPCA(fungi_t, # binary data
                           k=K, # number of PCs
                           m=0, # approximation of natural parameter
                           main_effects = TRUE,
                           partial_decomp = TRUE) # including offset term

fungi.K1 <- logpca.model

K = 2

## logistic PCA model
logpca.model = logisticPCA(fungi_t, # binary data
                           k=K, # number of PCs
                           m=0, # approximation of natural parameter
                           main_effects = TRUE,
                           partial_decomp = TRUE) # including offset term

fungi.K2 <- logpca.model
#fungi.K2.wo.asc <- logpca.model

logpca.scores = fungi.K2$PCs  
logpca.loadings = fungi.K2$U 

x <- as.data.frame(logpca.loadings)
x$rxns <- fungi_df$rxns
x$subsystems <- substring(fungi_df$subsystems, 11)


loading_results <- x %>%
  group_by(subsystems) %>%
  summarize(
    loading_1 = mean(V1),
    loading_2 = mean(V2)
  )

loading_results$length <- sqrt(loading_results$loading_1^2 + loading_results$loading_2^2)

# Sort the dataframe by 'length' in descending order
sorted_result <- loading_results %>% arrange(desc(length))

# Take the top 10 longest vectors
top_longest_vectors <- head(sorted_result, 5)

lipo_ellipse <- ellipse(cov(logpca.scores[lipomycetaceae,]), center=colMeans(logpca.scores[lipomycetaceae,]))
sacc_ellipse <- ellipse(cov(logpca.scores[saccharomycodaceae,]), center=colMeans(logpca.scores[saccharomycodaceae,]))

png(filename = "fungi_diff_rxns_lpca.png", width = 1200, height = 700)

plot(logpca.scores,  # x and y data
     pch=shape.labels,           # point shape
     col=class.labels, # 
     bg=class.labels,  #
     cex=1.5,          # point size
     main="",     # title of plot
     xlab = "PC1 (16%)", # for K = 1, m = 0 the PC = 0.163   sum(PC1, PC2) = 0.273
     ylab = "PC2 (11%)", # for K = 2, m = 0, the PC = 0.11
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

for(i in 1:nrow(top_longest_vectors)) {
  if (top_longest_vectors$loading_1[i] < -0.02) {
    text(x = top_longest_vectors$loading_1[i]*2500, 
         y = top_longest_vectors$loading_2[i]*2500, 
         labels = top_longest_vectors$subsystems[i],
         pos = 2,  # Places the text above the point
         col = "black",  # Text color
         cex = 1.3  # Text size
    )
  }
  if (top_longest_vectors$loading_1[i] > -0.02) {
    text(x = top_longest_vectors$loading_1[i]*2500, 
         y = top_longest_vectors$loading_2[i]*2500, 
         labels = top_longest_vectors$subsystems[i],
         pos = 4,  # Places the text above the point
         col = "black",  # Text color
         cex = 1.3  # Text size
    )
  }
}

lines(lipo_ellipse, col="#942d8c", lwd=2)
lines(sacc_ellipse, col="#04988d", lwd=2)

dev.off()


####################################
#### tsne plot on fungi dataset ####

tsne <- Rtsne(fungi_t, dims = 2, distance = 'hamming', check_duplicates = FALSE)

tsne_data <- as.data.frame(tsne$Y)
colnames(tsne_data) <- c('tSNE1', 'tSNE2')


tsne_data$color <- class.labels

png(filename = "fungi_diff_rxns_tsne.png", width = 1200, height = 700)

p <- plot(tsne_data$tSNE1,
          tsne_data$tSNE2,
          pch=shape.labels,           # point shape
          col=class.labels, # 
          bg=class.labels,  #
          cex=1.5,          # point size
          main="",     # title of plot
          xlab = "tSNE1", # 
          ylab = "tSNE2", # 
          cex.axis = 1.5,
          cex.lab = 1.5
)

dev.off()

##################################
#### Jaccard distance heatmap ####

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

df <- calculate_jaccard_similarity(fungi_t)

max_index <- max(c(index_1, index_2, index_3, index_4, index_5, index_6, index_7, index_8, index_9, index_10, index_11, index_12, index_13))
assignments <- vector("character", max_index)

# Function to fill the assignments vector
fill_assignments <- function(index_list, id, assignment_vector) {
  assignment_vector[index_list] <- as.character(id)
  return(assignment_vector)
}

# Fill in the assignments
assignments <- fill_assignments(index_1, 1, assignments)
assignments <- fill_assignments(index_2, 2, assignments)
assignments <- fill_assignments(index_3, 3, assignments)
assignments <- fill_assignments(index_4, 4, assignments)
assignments <- fill_assignments(index_5, 5, assignments)
assignments <- fill_assignments(index_6, 6, assignments)
assignments <- fill_assignments(index_7, 7, assignments)
assignments <- fill_assignments(index_8, 8, assignments)
assignments <- fill_assignments(index_9, 9, assignments)
assignments <- fill_assignments(index_10, 10, assignments)
assignments <- fill_assignments(index_11, 11, assignments)
assignments <- fill_assignments(index_12, 12, assignments)
assignments <- fill_assignments(index_13, 13, assignments)

group_factors <- factor(assignments, levels = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13), labels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13"))

group_colors <- c("1" = "#559668", "2" = "#000000", "3" = "#633717", "4" = "#fac34e", 
                  "5" = "#a7cd57", "6" = "#f83e32", "7" = "#942d8c", "8" = "#2cc5da", 
                  "9" = "#fc803e", "10" = "#0a5aa2", "11" = "#04988d", "12" = "#58539e", "13" = "#f5348f")

custom_colors <- colorRampPalette(c("yellow", "lightblue", "blue", "black"))(100)

annotation_col <- HeatmapAnnotation(
  Group = group_factors,
  col = list(Group = group_colors),
  which = "col" 
)

png("fungi_diff_rxns_jacc.png", width = 1200, height = 1200)

Heatmap(df, 
        top_annotation = annotation_col, 
        name = "Similarity", 
        column_title = "",
        col = custom_colors,
        column_names_side = "bottom",
        column_names_gp = gpar(col = NA)
)

dev.off()

###------------ legend --------------###

group_names <- c("Alloascoideaceae",
                 "Outgrouped",
                 "Cug_Ala",
                 "Cug_Ser1",
                 "Cug_Ser2",
                 "Dipodascaceae",
                 "Lipomycetaceae",
                 "Phaffomycetaceae",
                 "Pichiaceae",
                 "Saccharomycetaceae",
                 "Saccharomycodaceae",
                 "Sporopachydermia",
                 "Trigonopsidaceae"
)

# Vector with color assignments for each group
class.labels <- c("Alloascoideaceae" = "#559668",
                  "Outgrouped" = "#000000",
                  "Cug_Ala" = "#633717",
                  "Cug_Ser1" = "#fac34e",
                  "Cug_Ser2" = "#a7cd57",
                  "Dipodascaceae" = "#f83e32",
                  "Lipomycetaceae" = "#942d8c",
                  "Phaffomycetaceae" = "#2cc5da",
                  "Pichiaceae" = "#fc803e",
                  "Saccharomycetaceae" = "#0a5aa2",
                  "Saccharomycodaceae" = "#04988d",
                  "Sporopachydermia" = "#58539e",
                  "Trigonopsidaceae" = "#f5348f"
)

# Vector with shape assignments for each group
shape.labels <- c("Alloascoideaceae" = 4,
                  "Outgrouped" = 17,
                  "Cug_Ala" = 5,
                  "Cug_Ser1" = 19,
                  "Cug_Ser2" = 4,
                  "Dipodascaceae" = 8,
                  "Lipomycetaceae" = 9,
                  "Phaffomycetaceae" = 0,
                  "Pichiaceae" = 15,
                  "Saccharomycetaceae" = 2,
                  "Saccharomycodaceae" = 1,
                  "Sporopachydermia" = 6,
                  "Trigonopsidaceae" = 7
)

png(filename = "fungi_legend.png", width = 1500, height = 300)

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

####################################################
###### Cophenetic correlation fungi ################

tree <- read.tree("/home/users/lzehetner/data/logPCA/332_2408OGs_time-calibrated_phylogeny_species-names_updated.newick")
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