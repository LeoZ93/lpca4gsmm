# Differential logistic PCA for large-scale comparison of GSMMs and pathway analysis

Leopold Zehetner, Diana Szeliova, Barbara Kraus, Juan A. Hernandez Bort, and Jürgen Zanghellini

## Abstract
Genome-scale metabolic models (gsmms) provide a platform for metabolic analyses at the systems level by mathematically depicting biochemical reactions. Recent advances led to a large amount of gsmms for different species or tissue-specific reconstructions from global gsmms. Currently, the comparison of gsmms mostly depends on methods like dimension reduction algorithms, distance matrices, or predefined environmental conditions. However, these approaches are either hard to interpret or rely on subjective assumptions.
Here we introduce a novel use of logistic principal component analysis (lpca) to simplify the comparison of gsmms. Lpca is especially suited for binary data, derived from the presence or absence of reactions across gsmms. Its advantage is its ability to cluster gsmms while pinpointing the reactions that influence the categorization. When we applied lpca to three distinct datasets, it effectively highlighted phylogenetic relationships and distinctly separated the gsmms of healthy tissues from those of cancer tissues. Furthermore, the identified subsystems and their associated reactions that influenced the clustering align with existing knowledge.

## Data and Code

The scripts folder contains all scripts for data processing and to reconstruct results from the paper.  
1) extraction_cancer_healthy_tissue.ipynb: jupyter notebook to extract reactions and genes specific to transcriptomic datasets of cancer and healthy tissues from the Human Protein Atlas [1]
2) logpca.R: logistic pca performed on three different datasets, including 222 Escherichia-specific gsmms [2], 343 yeast-specific gsmms [3], and 80 healthy and cancer tissue specific reaction sets.
3) pca.R: pca performed on transcriptomic data from healthy and cancer tissues [1], as well as simulated growth rates from 222 Escherichia-specific gsmms [2]
4) phylogeny.R: phylogenetic tree for 222 Escherichia strains, based on their genomes, downloaded from Enterobase [4] and phylogenetic tree reconstruction based on 332 sequenced yeast genomes [5].
5) reaction_extraction_yeast.ipynb: extraction of pan-reactions from yeast-specific gsmms [3]
6) subsystem_extraction_human1.ipynb: assignment of subsystems to reactions, extracted from human1 gsmm [6]
7) tsne.R: reconstruction of tsne analysis based on yeast-specific gsmms [3].

Datasets newly generated for this study are stored in the data folder, which includes:  
1) SpeciesTree_rooted.txt: newick file obtained from 222 Escherichia strains, obtained from orthofinder [7].
2) cancer_and_tissue.csv: combined transcriptomic datasets containing healthy and cancer tissues, obtained from Human Protein Atlas [1]
3) differential_genes_tissue_and_cancer.csv: binary matrix containing extracted genes from human1 [6] assigned to reactions from tissue-specific reaction sets.  
4) differential_rxns_per_fungus.csv: binary matrix containing differential reactions from 343 yeast-specific gsmms [3].
5) differential_rxns_tissue_and_cancer.csv: binary matrix containing extracted reactions from human1 [6] based on transcriptomic data for healthy and cancer tissues, obtained from Human Protein Atlas [1].
6) extracted_rxns_per_fungus.csv: binary matrix containing pan-reactions from 343 yeast-specific gsmms [3].
7) subsystems.csv: assigned subsystems to reactions from human1 [6]
8) gpr_human1.csv: containing gene-protein reaction rules of human1, obtained from [6]

Additional datasets that need to be downloaded from other resources:  
1) [human1 gsmm](https://github.com/SysBioChalmers/Human-GEM/blob/main/model/Human-GEM.xml)  
2) [genes.csv from human1 gsmm](https://github.com/SysBioChalmers/Human-GEM/blob/main/model/genes.tsv)  
3) [binary reaction matrix of 222 Escherichia-specific gsmms](https://rs.figshare.com/articles/dataset/Supplementary_Data_File_2_from_Genome-scale_metabolic_network_reconstructions_of_diverse_i_Escherichia_i_strains_reveal_strain-specific_adaptations/20236554?backTo=/collections/Supplementary_material_from_Genome-scale_metabolic_network_reconstructions_of_diverse_i_Escherichia_i_strains_reveal_strain-specific_adaptations_/6080730)  
4) [simulated growth rates of 222 Escherichia-specific gsmms](https://rs.figshare.com/articles/dataset/Supplementary_Data_File_3_from_Genome-scale_metabolic_network_reconstructions_of_diverse_i_Escherichia_i_strains_reveal_strain-specific_adaptations/20236560?backTo=/collections/Supplementary_material_from_Genome-scale_metabolic_network_reconstructions_of_diverse_i_Escherichia_i_strains_reveal_strain-specific_adaptations_/6080730)  
5) [Phylogenetic results from 332 yeast-specific genomes](https://figshare.com/articles/dataset/Tempo_and_mode_of_genome_evolution_in_the_budding_yeast_subphylum/5854692?file=12977468)  
6) [Yeast-specific gsmms](https://www.ebi.ac.uk/biomodels/search?query=Lu2021&domain=biomodels)

## References
<a id="1">[1]</a> 
Uhlen M, Fagerberg L, Hallstrom BM, Lindskog C, Oksvold P, Mardinoglu A, et~al.
Tissue-based map of the human proteome.
Science. 2015;347(6220):1260419.  

<a id="2">[2]</a> 
Monk JM.
Genome-scale metabolic network reconstructions of diverse Escherichia strains reveal strain-specific adaptations.
Philosophical Transactions of the Royal Society B. 2022;377(1861):20210236.  

<a id="3">[3]</a> 
Lu H, Li F, Yuan L, Domenzain I, Yu R, Wang H, et al.
Yeast metabolic innovations emerged via expanded metabolic network and gene positive selection.
Molecular Systems Biology. 2021;17(10):e10427.  

<a id="4">[4]</a> 
Zhou, Zhemin, et al. 
The EnteroBase user's guide, with case studies on Salmonella transmissions, Yersinia pestis phylogeny, and Escherichia core genomic diversity. 
Genome research 30.1 (2020): 138-152.  

<a id="5">[5]</a> 
Shen, Xing-Xing, et al. 
Tempo and mode of genome evolution in the budding yeast subphylum. 
Cell 175.6 (2018): 1533-1545.

<a id="6">[6]</a> 
Robinson, Jonathan L., et al. 
An atlas of human metabolism.
Science signaling 13.624 (2020): eaaz1482.

<a id="7">[7]</a> 
Emms, David M., and Steven Kelly.
OrthoFinder: phylogenetic orthology inference for comparative genomics
Genome biology 20 (2019): 1-14.
