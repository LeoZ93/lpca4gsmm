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
1) 
