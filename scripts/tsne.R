## Performing t-sne on fungal GSMMs

library(readr)
library(Rtsne)
library(plotly)

# import dataset containing pan-reactions for fungi
df_1 <- read_csv("data/logPCA/extracted_rxns_per_fungus.csv")

df_1 <- df_1[, -c(1:3)]

# assign species to clades and label them
alloascoideaceae = c('Alloascoidea_hylecoeti')
outgrouped = c('Arthrobotrys_oligospora', 'Aspergillus_nidulans', 'Botrytis_cinerea', 'Coccidioides_immitis', 'Fusarium_graminearum', 'Neurospora_crassa', 'Saitoella_complicata', 'Schizosaccharomyces_pombe', 'Sclerotinia_sclerotiorum', 'Stagonospora_nodorum', 'Xylona_heveae')
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

index_1 <- match(alloascoideaceae, names(df_1))
index_2 <- match(cug_ala, names(df_1))
index_3 <- match(cug_ser1, names(df_1))
index_4 <- match(cug_ser2, names(df_1))
index_5 <- match(dipodascaceae, names(df_1))
index_6 <- match(lipomycetaceae, names(df_1))
index_7 <- match(outgrouped, names(df_1))
index_8 <- match(phaffomycetaceae, names(df_1))
index_9 <- match(pichiaceae, names(df_1))
index_10 <- match(saccharomycetaceae, names(df_1))
index_11 <- match(saccharomycodaceae, names(df_1))
index_12 <- match(sporopachydermia, names(df_1))
index_13 <- match(trigonopsidaceae, names(df_1))

class.labels = vector()
class.labels[c(index_1)] = "#559668" # alloascoideaceae
class.labels[c(index_2)] = "#633717" # cug_ala
class.labels[c(index_3)] = "#fac34e" # cug_ser1
class.labels[c(index_4)] = "#a7cd57" # cug_ser2
class.labels[c(index_5)] = "#f83e32" # dipodascaceae
class.labels[c(index_6)] = "#942d8c" # lipomycetaceae
class.labels[c(index_7)] = "#000000" # outgrouped
class.labels[c(index_8)] = "#2cc5da" # phaffomycetaceae
class.labels[c(index_9)] = "#fc803e" # pichiaceae
class.labels[c(index_10)] = "#0a5aa2" # saccharomycetaceae
class.labels[c(index_11)] = "#04988d" # saccharomycodaceae
class.labels[c(index_12)] = "#58539e" # sporopachydermia
class.labels[c(index_13)] = "#f5348f" # trigonopsidaceae

df_1 <- t(df_1)

# perform tsne analysis with 3 dimensions
tsne <- Rtsne(df_1, dims = 3, distance = 'hamming', check_duplicates = FALSE)

tsne_data <- as.data.frame(tsne$Y)
colnames(tsne_data) <- c('tSNE1', 'tSNE2', 'tSNE3')

tsne_data$color <- class.labels

# Ensure the 'color' column is a factor and map it correctly in the plot
tsne_data$color <- as.factor(tsne_data$color)

# Plot using plot_ly with corrected color referencing
plot_ly(tsne_data, x = ~tSNE1, y = ~tSNE2, z = ~tSNE3, type = 'scatter3d', mode = 'markers',
        marker = list(color = ~color, colorscale = 'Viridis', size = 3)) %>%
  layout(title = '',
         scene = list(xaxis = list(title = 'tSNE1'),
                      yaxis = list(title = 'tSNE2'),
                      zaxis = list(title = 'tSNE3')))

### Create 2D tsne plot

tsne <- Rtsne(df_1, dims = 2, distance = 'hamming', check_duplicates = FALSE)

tsne_data <- as.data.frame(tsne$Y)
colnames(tsne_data) <- c('tSNE1', 'tSNE2')


tsne_data$color <- class.labels


# Ensure the 'color' column is a factor and map it correctly in the plot
#tsne_data$color <- as.factor(tsne_data$color)

plot(tsne_data$tSNE1,
     tsne_data$tSNE2,
     pch=19,           # point shape
     col=class.labels, # 
     bg=class.labels,  #
     cex=1,          # point size
     main="",     # title of plot
     xlab = "tSNE1", # 
     ylab = "tSNE2", # 
     cex.axis = 1.4,
     cex.lab = 1.4
)

