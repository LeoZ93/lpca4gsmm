import cobra
import pandas as pd
import os

loaded_models = []

# Specify the directory where the XML files are located
folder_path = '/home/users/lzehetner/fungi/fungal_models_Lu_2021/fungal_models/'

rxns_df = pd.DataFrame(columns=['rxns', 'subsystems'])

for filename in os.listdir(folder_path):
    if filename.endswith('.xml'):  # Ensure the file is an XML file
        full_path = os.path.join(folder_path, filename)
        model_name = filename[:-4]
        try:
            model = cobra.io.read_sbml_model(full_path) # Read the model into a COBRA object
            model.id = model_name
            loaded_models.append(model)  # Add the model to the list
            rxn_subsystem_dict = {}
    
            for reaction in model.reactions:
                subsystem = reaction.subsystem if reaction.subsystem else "None"  # Use "None" if subsystem is empty
                rxn_subsystem_dict[reaction.id] = subsystem

            # Create a temporary DataFrame for this model
            temp_df = pd.DataFrame({
#                'model_name': [model_name] * len(rxn_subsystem_dict),
                'rxns': list(rxn_subsystem_dict.keys()),
                'subsystems': list(rxn_subsystem_dict.values())
            })

            # Append this temporary DataFrame to the main DataFrame
            rxns_df = pd.concat([rxns_df, temp_df], ignore_index=True)

#            print(f'Successfully loaded {filename}')
        except Exception as e:
            print(f'Could not load {filename}: {e}')

extracted_rxns_per_fungus = rxns_df.drop_duplicates(subset=['rxns'])

for model in loaded_models:
    mod_rxns = []
    for r in model.reactions:
        mod_rxns.append(r.id)
    extracted_rxns_per_fungus[model.id] = extracted_rxns_per_fungus["rxns"].isin(mod_rxns).astype(int)

rows_to_remove = extracted_rxns_per_fungus.iloc[:, 2:].apply(lambda row: row.sum() == 0 or row.sum() == len(loaded_models), axis=1)

# Drop those rows
differential_rxns_per_fungus = extracted_rxns_per_fungus.loc[~rows_to_remove]

differential_rxns_per_fungus.to_csv('/home/users/lzehetner/data/logPCA/differential_rxns_per_fungus.csv')