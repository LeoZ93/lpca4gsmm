# lpca4gsmm

Differential logistic PCA for large-scale comparison of GSMMs and pathway analysis

Leopold Zehetner, Diana Szeliova, Barbara Kraus, Juan A. Hernandez Bort, and JÜrgen Zanghellini

Genome-scale metabolic models (gsmms) provide a platform for metabolic analyses at the systems level by mathematically depicting biochemical reactions. Recent advances led to a large amount of gsmms for different species or tissue-specific reconstructions from global \gsmms. Currently, the comparison of gsmms mostly depends on methods like dimension reduction algorithms, distance matrices, or predefined environmental conditions. However, these approaches are either hard to interpret or rely on subjective assumptions.
Here we introduce a novel use of logistic principal component analysis (lpca) to simplify the comparison of gsmms. Lpca is especially suited for binary data, derived from the presence or absence of reactions across gsmms. Its advantage is its ability to cluster gsmms while pinpointing the reactions that influence the categorization. When we applied lpca to three distinct datasets, it effectively highlighted phylogenetic relationships and distinctly separated the gsmms of healthy tissues from those of cancer tissues. Furthermore, the identified subsystems and their associated reactions that influenced the clustering align with existing knowledge.

