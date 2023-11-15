## Performing logpca on binary matrices from GSMMs

library(logisticPCA) # for logistic PCA ans logistic SVD
library(readr)
library(rARPACK)
library(readxl)

###############################
##### PanGSMM Escherichia #####
###############################

escherichia <- read_excel("data/logPCA/rstb20210236_si_002.xlsx", sheet = "Reaction Info") # file obtained from Monk, 2022, 
clades <- read_excel("data/logPCA/rstb20210236_si_002.xlsx", sheet = "Model reaction counts")

# data cleaning and preparation
escherichia <- escherichia[, -c(2:5)]

escherichia[is.na(escherichia)] <- 0.0

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

rem.clades.cluster <- c(index_lists[[2]], index_lists[[3]], index_lists[[4]], index_lists[[5]], index_lists[[6]], index_lists[[8]], index_lists[[9]], index_lists[[11]], index_lists[[12]], index_lists[[13]])

escherichia_t <- t(escherichia)

class.labels = vector()
class.labels[c(E.albertii)] = "#127db9" # 
class.labels[c(Clade.III)] = "#2c9f28" # 
class.labels[c(S.flexneri)] = "#ffe4ca" # 
class.labels[c(E.coli)] = "#fd7d08" # 
class.labels[c(S.dysenteriae)] = "#feb166" # 
class.labels[c(Clade.V)] = "#abcaff" # 
class.labels[c(E.fergusonii)] = "#e46864" # 
class.labels[c(Clade.VIII)] = "#915347" # 
class.labels[c(Clade.II)] = "#c5aedb" # 
class.labels[c(Clade.IV)] = "#90e585" # 
class.labels[c(Clade.VI)] = "#feb0d2" # 
class.labels[c(Clade.VII)] = "#9bdeef" # 

# Running logpca on differential reactions
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

escherichia.model <- logpca.model

PC2 <- escherichia.model$prop_deviance_expl - PC1

# analyzing results from logpca
logpca.scores = escherichia.model$PCs # extract score matrix
logpca.loadings = escherichia.model$U # extract loading matrix

x <- as.data.frame(logpca.loadings)
x$rxns <- differential.rxns
x$subsystems <- subsystems

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

rem.cluster_ellipse <- ellipse(cov(logpca.scores[rem.clades.cluster,]), center=colMeans(logpca.scores[rem.clades.cluster,]))
ealbertii_ellipse <- ellipse(cov(logpca.scores[E.albertii,]), center=colMeans(logpca.scores[E.albertii,]))
efergusonii_ellipse <- ellipse(cov(logpca.scores[E.fergusonii,]), center=colMeans(logpca.scores[E.fergusonii,]))

# plotting results from logpca
plot(logpca.scores,  # x and y data
     pch=21,           # point shape
     col=class.labels, # 
     bg=class.labels,  #
     cex=1,          # point size
     main="",     # title of plot
     xlab = "PC1 (19%)", # 
     ylab = "PC2 (15%)", #
     cex.axis = 1.3,
     cex.lab = 1.3
)

arrows(x0 = rep(0, nrow(top_longest_vectors)), 
       y0 = rep(0, nrow(top_longest_vectors)), 
       x1 = top_longest_vectors$loading_1*2000, 
       y1 = top_longest_vectors$loading_2*2000, 
       col = "black", # You can change the color
       angle = 25, # Angle of the arrow head
       length = 0.1) # Length of the arrow head

for(i in 1:nrow(top_longest_vectors)) {
  if (top_longest_vectors$loading_1[i] < 0) {
    text(x = top_longest_vectors$loading_1[i]*2500, 
         y = top_longest_vectors$loading_2[i]*2500, 
         labels = top_longest_vectors$subsystems[i],
         pos = 2,  # Places the text above the point
         col = "black",  # Text color
         cex = 1.3  # Text size
    )
  }
  if (top_longest_vectors$loading_1[i] > 0) {
    text(x = top_longest_vectors$loading_1[i]*2500, 
         y = top_longest_vectors$loading_2[i]*2500, 
         labels = top_longest_vectors$subsystems[i],
         pos = 4,  # Places the text above the point
         col = "black",  # Text color
         cex = 1.3  # Text size
    )
  }
}

colors <- c("#fd7d08", "#e46864", "#127db9","#feb166","#ffe4ca","#c5aedb","#2c9f28","#90e585","#abcaff","#feb0d2","#9bdeef", "#915347")
legend("bottomleft", legend = c("E.coli", "E.fergusonii", "E.albertii", "S.dysenteriae", "S.flexneri", "Clade II", "Clade III", "Clade IV", "Clade V", "Clade VI", "Clade VII", "Clade VIII"), fill=colors, cex=0.8)

lines(rem.cluster_ellipse, col="#fd7d08", lwd=2)
lines(ealbertii_ellipse, col="#127db9", lwd=2)
lines(efergusonii_ellipse, col="#e46864", lwd=2)


#### different approach of calculating loadings: first take length of vectors, then mean by subsystem
#### Cave: visualization is not possible

logpca.scores = escherichia.model$PCs # extract score matrix
logpca.loadings = escherichia.model$U # extract loading matrix

x <- as.data.frame(logpca.loadings)
x$rxns <- differential.rxns
x$subsystems <- subsystems

x$length <- sqrt(x$V1^2 + x$V2^2)

loading_results <- x %>%
  group_by(subsystems) %>%
  summarize(
    loading = mean(length)
  )

sorted_result <- loading_results %>% arrange(desc(loading))

# Take the top 10 longest vectors
top_longest_vectors.2 <- head(sorted_result, 5)

# top_longest_vectors and top_longest_vectors.2 correspond to system-loadings 1 and system-loadings 2 in table 1 of the paper


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

# Take the top 10 longest vectors
top_longest_vectors <- head(sorted_result, 5)
#top_longest_vectors <- top_longest_vectors[c(3:6),]

plot(logpca.scores,  # x and y data
     pch=21,           # point shape
     col=class.labels, # 
     bg=class.labels,  #
     cex=1.5,          # point size
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

for(i in 1:nrow(top_longest_vectors)) {
  if (top_longest_vectors$loading_1[i] < 0) {
    text(x = top_longest_vectors$loading_1[i]*2500, 
         y = top_longest_vectors$loading_2[i]*2500, 
         labels = top_longest_vectors$subsystems[i],
         pos = 2,  # Places the text above the point
         col = "black",  # Text color
         cex = 0.8  # Text size
    )
  }
  if (top_longest_vectors$loading_1[i] > 0) {
    text(x = top_longest_vectors$loading_1[i]*2500, 
         y = top_longest_vectors$loading_2[i]*2500, 
         labels = top_longest_vectors$subsystems[i],
         pos = 4,  # Places the text above the point
         col = "black",  # Text color
         cex = 0.8  # Text size
    )
  }
}

colors <- c("#0000FF", "#FF0000")
legend("bottomright", legend=c("Normal Tissue", "Cancer"), fill=colors, cex=1)


selected_rxns <- x %>%
  filter(subsystems == "Linoleate metabolism")



###################################
############# Fungi ###############
###################################

fungi_df <- read_csv("/home/users/lzehetner/data/logPCA/differential_rxns_per_fungus.csv", col_names = TRUE)

fungi <- fungi_df[, -c(1,2,3)]

# spec_cell_rxns <- spec_cell_rxns[, -c(21,69)]

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

class.labels = vector()
class.labels[c(index_1)] = "#559668" # alloascoideaceae
class.labels[c(index_2)] = "#000000" # ascomycota
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

logpca.scores = fungi.K2$PCs # extract score matrix
logpca.loadings = fungi.K2$U # extract loading matrix

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
#top_longest_vectors$subsystems <- substring(top_longest_vectors$subsystems, 9)

# y <- x[x$V1 < -0.02, ]

#png(file="/home/users/lzehetner/data/logPCA/logpca_spec_cells_test.png",
#    width=1059, height=660)

plot(logpca.scores,  # x and y data
     pch=21,           # point shape
     col=class.labels, # 
     bg=class.labels,  #
     cex=1,          # point size
     main="",     # title of plot
     xlab = "PC1 (16%)", # for K = 1, m = 0 the PC = 0.163   sum(PC1, PC2) = 0.273
     ylab = "PC2 (11%)", # for K = 2, m = 0, the PC = 0.11
     cex.axis = 1.3,
     cex.lab = 1.3
)

arrows(x0 = rep(0, nrow(top_longest_vectors)), 
       y0 = rep(0, nrow(top_longest_vectors)), 
       x1 = top_longest_vectors$loading_1*2000, 
       y1 = top_longest_vectors$loading_2*2000, 
       col = "black", # You can change the color
       angle = 25, # Angle of the arrow head
       length = 0.1) # Length of the arrow head

for(i in 1:nrow(top_longest_vectors)) {
  if (top_longest_vectors$loading_1[i] < 0) {
    text(x = top_longest_vectors$loading_1[i]*2500, 
         y = top_longest_vectors$loading_2[i]*2500, 
         labels = top_longest_vectors$subsystems[i],
         pos = 2,  # Places the text above the point
         col = "black",  # Text color
         cex = 1.3  # Text size
    )
  }
  if (top_longest_vectors$loading_1[i] > 0) {
    text(x = top_longest_vectors$loading_1[i]*2500, 
         y = top_longest_vectors$loading_2[i]*2500, 
         labels = top_longest_vectors$subsystems[i],
         pos = 4,  # Places the text above the point
         col = "black",  # Text color
         cex = 1.3  # Text size
    )
  }
}

colors <- c("#559668", "#000000", "#633717", "#fac34e", "#a7cd57", "#f83e32", "#942d8c", "#2cc5da", "#fc803e", "#0a5aa2", "#04988d", "#58539e", "#f5348f")

legend("topleft", legend=c("Alloascoideaceae", "Outgroup", "Cug_Ala", "Cug_Ser1", "Cug_Ser2", "Dipodascaceae", "Lipomycetaceae", "Phaffomycetaceae", "Pichiaceae", "Saccharomycetaceae", "Saccharomycodaceae", "Sporopachydermia", "Trigonopsidaceae"), fill=colors, cex=0.8)

