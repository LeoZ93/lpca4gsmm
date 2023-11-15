import cobra
import pandas as pd
# import numpy as np
import re

## preparation of gpr dataset from human-gsmm
gpr_ori = pd.read_csv("/home/users/lzehetner/gpr_human1.csv", sep = ";")
mod = cobra.io.read_sbml_model("/home/users/lzehetner/data/human1/human1.xml")

gpr = gpr_ori.iloc[:, [0,2]]
rxns_wo_genes = gpr[gpr.iloc[:, 1].isna()]
gpr.dropna(inplace=True)
gpr_parent = gpr[gpr['Gene-reaction association'].str.contains('\(')]
gpr = gpr[~gpr.isin(gpr_parent)].dropna()
gpr_and = gpr[gpr['Gene-reaction association'].str.contains('and')]
gpr = gpr[~gpr.isin(gpr_and)].dropna()
gpr_or = gpr[gpr['Gene-reaction association'].str.contains('or')]
gpr_one = gpr[~gpr.isin(gpr_or)].dropna()
gpr_parent_and = gpr_parent[gpr_parent['Gene-reaction association'].str.contains(' and \(')]
gpr_parent_or = gpr_parent[~gpr_parent.isin(gpr_parent_and)].dropna()
gpr_parent_and = gpr_parent_and.iloc[:-1 , :]


## define a function for reaction extraction from transcriptomic data

def rxns_extraction(model, df, gene_column, column):
    genes = pd.DataFrame(df[[gene_column, column]])
    genes = df[df[column] > 0.2]

    rxns = []
    for r in model.reactions:
        a = model.reactions.get_by_id(r.id).gene_reaction_rule
        if any(a in x for x in gpr_or.iloc[:, 1]):
            b = a.split(' or ')
            relevant_genes = genes[genes[gene_column].isin(b)]
            if len(relevant_genes) != 0:
                rxns.append(r.id)
        elif any(a in x for x in gpr_and.iloc[:, 1]):
            b = a.split(' and ')
            relevant_genes = genes[genes[gene_column].isin(b)]
            if len(relevant_genes) == len(b):
                rxns.append(r.id)
        elif any(a in x for x in gpr_one.iloc[:, 1]):
            if any(a in x for x in genes.iloc[:, 0]):
                rxns.append(r.id)      
        elif any(a in x for x in gpr_parent_or.iloc[:, 1]):
            parentheses = a.split(') or (')
            parentheses = [group.strip(' () ') for group in parentheses]   
            for group in parentheses:
                b = group.split(' and ') 
                relevant_genes = genes[genes[gene_column].isin(b)]
                if len(relevant_genes) == len(b):
                    rxns.append(r.id)
        elif any(a in x for x in gpr_parent_and.iloc[:, 1]):
            parentheses_groups = [group.strip(' () ') for group in a.split(' and ')]
            for group in parentheses_groups:
                b = group.split(' or ')
                relevant_genes = genes[genes[gene_column].isin(b)]
                if len(relevant_genes) != 0:
                    rxns.append(r.id)
    unique_list = []
    for item in rxns:
        if item not in unique_list:
            unique_list.append(item)
    return unique_list


#### Extrqaction of tissue and cancer specific reactions
## extraction of cancer tissue specific reactions 
data = pd.read_csv("/home/users/lzehetner/data/logPCA/rna_celline_cancer.tsv", sep = "\t")

# prepare transcriptomic dataset for reaction extraction
df = pd.DataFrame(data)
reshaped_df = df.pivot(index='Gene', columns='Cancer', values='nTPM')
cell_types = reshaped_df.columns.to_list()
df_reset = reshaped_df.reset_index()
mod_rxns = []
for r in mod.reactions:
    mod_rxns.append(r.id)

# extract cancer specific reactions
extracted_rxns_per_cell_1 = pd.DataFrame()
extracted_rxns_per_cell_1["Reaction"] = mod_rxns
for cell in cell_types:
    cell_rxns = rxns_extraction(mod, df_reset, "Gene", cell)
    extracted_rxns_per_cell_1[cell] = extracted_rxns_per_cell_1["Reaction"].isin(cell_rxns).astype(int)

extracted_rxns_per_cell_1 = extracted_rxns_per_cell_1.drop(extracted_rxns_per_cell_1.columns[0], axis=1)



## extraction of tissue specific reaction
data = pd.read_csv("/home/users/lzehetner/data/logPCA/rna_tissue_consensus.tsv", sep = "\t")

# prepare transcriptomic dataset for reaction extraction
df = pd.DataFrame(data)
reshaped_df = df.pivot(index='Gene', columns='Tissue', values='nTPM')
cell_types = reshaped_df.columns.to_list()
df_reset = reshaped_df.reset_index()
mod_rxns = []
for r in mod.reactions:
    mod_rxns.append(r.id)

# extract tissue specific reactions
extracted_rxns_per_cell = pd.DataFrame()
extracted_rxns_per_cell["Reaction"] = mod_rxns
for cell in cell_types:
    cell_rxns = rxns_extraction(mod, df_reset, "Gene", cell)
    extracted_rxns_per_cell[cell] = extracted_rxns_per_cell["Reaction"].isin(cell_rxns).astype(int)

# combine tissue and cancer specific reaction datasets
result = pd.concat([extracted_rxns_per_cell, extracted_rxns_per_cell_1], axis=1)
rows_to_remove = result.iloc[:, 1:].apply(lambda row: row.sum() == 0 or row.sum() == len(result.columns) - 1, axis=1)

differential_rxns_per_cell = result.loc[~rows_to_remove]
differential_rxns_per_cell.to_csv('/home/users/lzehetner/data/logPCA/differential_rxns_tissue_and_cancer.csv', )




# combine tissue and cancer transcriptomes

data = pd.read_csv("/home/users/lzehetner/data/logPCA/rna_celline_cancer.tsv", sep = "\t")
df = pd.DataFrame(data)
reshaped_df = df.pivot(index='Gene', columns='Cancer', values='nTPM')
df_reset_1 = reshaped_df.reset_index()
data = pd.read_csv("/home/users/lzehetner/data/logPCA/rna_tissue_consensus.tsv", sep = "\t")
df = pd.DataFrame(data)
reshaped_df = df.pivot(index='Gene', columns='Tissue', values='nTPM')
df_reset = reshaped_df.reset_index()
merged_df = df_reset.merge(df_reset_1, on='Gene', how='outer')
merged_df.to_csv('/home/users/lzehetner/data/logPCA/cancer_and_tissue.csv', )

# extraction of gsmm specific genes

diff_rxns = differential_rxns_per_cell["Reaction"]
gpr_ori = pd.read_csv("/home/users/lzehetner/gpr_human1.csv", sep = ";")
df1 = gpr_ori[gpr_ori.iloc[:, 0].isin(diff_rxns)]
gene_matches = df1.iloc[:, 2].str.extractall(r'(ENSG\d+)')
all_genes = gene_matches[0].tolist()
unique_genes = list(set(all_genes))
data = pd.read_csv("/home/users/lzehetner/data/logPCA/rna_celline_cancer.tsv", sep = "\t")
df = pd.DataFrame(data)
reshaped_df = df.pivot(index='Gene', columns='Cancer', values='nTPM')
df_reset_1 = reshaped_df.reset_index()
data = pd.read_csv("/home/users/lzehetner/data/logPCA/rna_tissue_consensus.tsv", sep = "\t")
df = pd.DataFrame(data)
reshaped_df = df.pivot(index='Gene', columns='Tissue', values='nTPM')
df_reset = reshaped_df.reset_index()
merged_df = pd.merge(df_reset, df_reset_1, on=df_reset.columns[0], how='inner')
filtered_merged_df = merged_df[merged_df.iloc[:, 0].isin(unique_genes)]
filtered_merged_df.to_csv('/home/users/lzehetner/data/logPCA/differential_genes_tissue_and_cancer.csv', )


## extract subsystems for every reaction from human1 gsmm

mod = cobra.io.read_sbml_model("/home/users/lzehetner/data/human1/human1.xml")
rxn_subsystem_dict = {}
    
for reaction in mod.reactions:
    subsystem = reaction.subsystem if reaction.subsystem else "None"  # Use "None" if subsystem is empty
    rxn_subsystem_dict[reaction.id] = subsystem

subsystems = pd.DataFrame({
    'rxns': list(rxn_subsystem_dict.keys()),
    'subsystems': list(rxn_subsystem_dict.values())
})

subsystems.to_csv('/home/users/lzehetner/data/logPCA/subsystems.csv')
## files necessary

# subsystems.csv
# differential_rxns_tissue_and_cancer.csv
# cancer_and_tissue.csv
# genes.tsv
# Ensembl2Reactome.txt
# differential_genes_tissue_and_cancer.csv
# cancer_and_tissue.csv



# rstb20210236_si_003.xlsx", sheet = "growth predictions"