### Reconstructing Phylogenetic tree from fungal species

library(ape)



#####################
#### Escherichia ####
#####################


## Whole genome based phylogeny

tree <- read.tree("/home/users/lzehetner/data/prodigal/proteomes/OrthoFinder/Results_Nov14/Species_Tree/SpeciesTree_rooted.txt")

data <- read_excel("data/rstb20210236_si_003.xlsx", sheet = "growth predictions")

data <- data[-223,]

# Split the data into a list by species
split_data <- split(data$strain_id, data$Species)

# Convert list elements to vectors and name them
species_vectors <- lapply(split_data, as.vector)
#names(species_vectors) <- unique(data$Species)

# Print the vectors
# print(species_vectors)

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



###############
#### Fungi ####
###############

tree <- read.tree("/home/users/lzehetner/data/logPCA/332_2408OGs_time-calibrated_phylogeny_species-names_updated.newick")

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


groups <- list(alloascoideaceae, outgrouped, cug_ala, cug_ser1, cug_ser2, dipodascaceae, lipomycetaceae, phaffomycetaceae, pichiaceae, saccharomycetaceae, saccharomycodaceae, sporopachydermia, trigonopsidaceae)

#groups <- list(ascomycota) # Add more groups as necessary
colors <- c("#559668", "#000000", "#633717", "#fac34e", "#a7cd57", "#f83e32", "#942d8c", "#2cc5da", "#fc803e", "#0a5aa2", "#04988d", "#58539e", "#f5348f") # Add more colors corresponding to your groups

assign_edge_colors <- function(tree, groups, colors) {
  edge_colors <- rep("black", length(tree$edge[,1])) # Default color for edges
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

plot(tree, edge.color = edge_colors, show.tip.label = FALSE)

legend("bottomleft", legend=c("Alloascoideaceae", "Outgrouped", "Cug_Ala", "Cug_Ser1", "Cug_Ser2", "Dipodascaceae", "Lipomycetaceae", "Phaffomycetaceae", "Pichiaceae", "Saccharomycetaceae", "Saccharomycodaceae", "Sporopachydermia", "Trigonopsidaceae"), fill=colors, cex=0.8)
