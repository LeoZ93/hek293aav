# # flux sampling file

from cobra.sampling import sample
from cobra.sampling import OptGPSampler, ACHRSampler
import pandas as pd
import numpy as np
import cobra
from cobra import Model, Reaction, Metabolite
import matplotlib.pyplot as plt
import seaborn as sns
import time
from Bio import SeqIO, Seq
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from cobra.flux_analysis.loopless import loopless_solution

reaction_matrix = pd.read_csv("/home/users/lzehetner/data/hek/reaction_matrix.csv")

# I_T_04 = ["lp_v1_tp2", "lp_v2_tp2", "lp_v3_tp2", "lp_v4_tp2"]
# I_T_24 = ["lp_v1_tp3", "lp_v2_tp3", "lp_v3_tp3", "lp_v4_tp3"]
# I_T_48 = ["lp_v1_tp4", "lp_v2_tp4", "lp_v3_tp4", "lp_v4_tp4"]
# I_T_72 = ["lp_v1_tp5", "lp_v2_tp5", "lp_v3_tp5", "lp_v4_tp5"]

# I_M_04 = ["lp_v5_tp2", "lp_v6_tp2", "lp_v7_tp2", "lp_v8_tp2"]
# I_M_24 = ["lp_v5_tp3", "lp_v6_tp3", "lp_v7_tp3", "lp_v8_tp3"]
# I_M_48 = ["lp_v5_tp4", "lp_v6_tp4", "lp_v7_tp4", "lp_v8_tp4"]
# I_M_72 = ["lp_v5_tp5", "lp_v6_tp5", "lp_v7_tp5", "lp_v8_tp5"]

#P_T_04 = ["hp_v1_tp2", "hp_v2_tp2", "hp_v3_tp2", "hp_v4_tp2"]
#P_T_24 = ["hp_v1_tp3", "hp_v2_tp3", "hp_v3_tp3", "hp_v4_tp3"]
#P_T_48 = ["hp_v1_tp4", "hp_v2_tp4", "hp_v3_tp4", "hp_v4_tp4"]
#P_T_72 = ["hp_v1_tp5", "hp_v2_tp5", "hp_v3_tp5", "hp_v4_tp5"]

P_M_04 = ["hp_v5_tp2", "hp_v6_tp2", "hp_v7_tp2", "hp_v8_tp2"]
P_M_24 = ["hp_v5_tp3", "hp_v6_tp3", "hp_v7_tp3", "hp_v8_tp3"]
P_M_48 = ["hp_v5_tp4", "hp_v6_tp4", "hp_v7_tp4", "hp_v8_tp4"]
P_M_72 = ["hp_v5_tp5", "hp_v6_tp5", "hp_v8_tp5"]

input_rates = "/home/users/lzehetner/data/hek/specific_exchange_rates.csv"
df = pd.read_csv(input_rates)

filtered_dfs = {}

states = ['HP_TR', 'HP_MO', 'LP_TR', 'LP_MO']
time_points = [4, 24, 48, 72]

for state in states:
    for time_point in time_points:
        key = f"{state}_{time_point}"
        filtered_dfs[key] = df[(df['state'] == state) & (df['time_point'] == time_point)]


#df_hp_tr_4 = filtered_dfs['HP_TR_4']
#df_hp_tr_24 = filtered_dfs['HP_TR_24']
#df_hp_tr_48 = filtered_dfs['HP_TR_48']
#df_hp_tr_72 = filtered_dfs['HP_TR_72']

df_hp_mo_4 = filtered_dfs['HP_MO_4']
df_hp_mo_24 = filtered_dfs['HP_MO_24']
df_hp_mo_48 = filtered_dfs['HP_MO_48']
df_hp_mo_72 = filtered_dfs['HP_MO_72']

# df_lp_tr_4 = filtered_dfs['LP_TR_4']
# df_lp_tr_24 = filtered_dfs['LP_TR_24']
# df_lp_tr_48 = filtered_dfs['LP_TR_48']
# df_lp_tr_72 = filtered_dfs['LP_TR_72']

# df_lp_mo_4 = filtered_dfs['LP_MO_4']
# df_lp_mo_24 = filtered_dfs['LP_MO_24']
# df_lp_mo_48 = filtered_dfs['LP_MO_48']
# df_lp_mo_72 = filtered_dfs['LP_MO_72']

# remove last row with ammonia, since there are no exchange rate fitted
#df_hp_tr_72 = df_hp_tr_72[:-1]

metabolite_to_reaction = {'Alanine': 'MAR09061',
                         'Arginine': 'MAR09066',
                         'Asparagine': 'MAR09062',
                         'Aspartic acid': 'MAR09070',
                         'Glutamate': 'MAR09071',
                         'Glutamine': 'MAR09063',
                         'Gly': 'MAR09067',
                         'Histidine': 'MAR09038',
                         'Isoleucine': 'MAR09039',
                         'Leucine': 'MAR09040',
                         'Lysine': 'MAR09041',
                         'Methionine': 'MAR09042',
                         'Phenylalanine': 'MAR09043',
                         'Pro': 'MAR09068',
                         'Serine': 'MAR09069',
                         'Threonine': 'MAR09044',
                         'Tryptophan': 'MAR09045',
                         'Tyrosine': 'MAR09064',
                         'Valine': 'MAR09046',
                         'Glucose': 'MAR09034',
                         'Nh3': 'MAR11420',
                         'Lac': 'MAR09135'
                         }

input_model = "/home/users/lzehetner/data/human1/human1.xml"
mod = cobra.io.read_sbml_model(input_model)
mod.solver = "cplex"

min_media_comp = [
    "MAR09061",
    "MAR09034",
    "MAR09035",
    "MAR09036",
    "MAR09038",
    "MAR09039",
    "MAR09040",
    "MAR09041",
    "MAR09042",
    "MAR09043",
    "MAR09044",
    "MAR09045",
    "MAR09046",
    "MAR09047",
    "MAR09048",
    "MAR09062",
    "MAR09063",
    "MAR09064",
    "MAR09065",
    "MAR09066",
    "MAR09068",
    "MAR09069",
    "MAR09076",
    "MAR09146",
    "MAR09109",
    "MAR09143",
    "MAR09144",
    "MAR09147",
    "MAR09150",
    "MAR09151",
    "MAR09153",
    "MAR09158",
    "MAR09159",
    "MAR09167",
    "MAR09269",
    "MAR09072",
    "MAR09145",
    "MAR09070",
    "MAR09071",
    "MAR09067",
    "MAR11420",
    "MAR09135"
]


def growth_reconstruction(column_list, df_state, binary_matrix):
    
    mod = cobra.io.read_sbml_model(input_model)
    mod.solver = "cplex"
    
    filtered_matrix = binary_matrix[column_list]
    mask = filtered_matrix.any(axis=1)
    rxns = binary_matrix.loc[mask, "Rxn name"].tolist()
    
    rxns.extend(["MAR07861", "MAR07862", "MAR07857"]) # add reactions for A, G, T into the nucleus
    
    b = []
    for r in mod.reactions:
        b.append(mod.reactions.get_by_id(r.id))
#    len(b)
    
    c = [item.id for item in b if item.id not in rxns]
    
    mod.remove_reactions(c)
    len(mod.reactions)
    
    a = mod.exchanges
    b = []
    for r in a:
        b.append(r.id)
    # print(len(b))
    for r in b:
        mod.reactions.get_by_id(r).upper_bound = 1000
        mod.reactions.get_by_id(r).lower_bound = 0

    for r in min_media_comp:
        mod.reactions.get_by_id(r).lower_bound = -1000
    
    for index, row in df_state.iterrows():
        metabolite = row['metabolite']
    
        if metabolite in metabolite_to_reaction:
            reaction_id = metabolite_to_reaction[metabolite]
            try:
                reaction = mod.reactions.get_by_id(reaction_id)
                reaction.lower_bound = row['lb']
                reaction.upper_bound = row['ub']
#                print(f"Set bounds for {reaction_id}: lower = {row['lb']}, upper = {row['ub']}")
            except KeyError:
                print(f"Reaction {reaction_id} not found in the model.")
        else:
            continue
            print(f"Metabolite {metabolite} not found in the dictionary.")
    
    mod.objective = "MAR13082"
    growth_rate = mod.optimize().objective_value
    
    orphan_metabolites = [m for m in mod.metabolites if len(m.reactions) == 0]
    mod.remove_metabolites(orphan_metabolites)
    
    # Remove orphan genes (not associated with any reactions)
    orphan_genes = [g for g in mod.genes if len(g.reactions) == 0]
    mod.genes = [g for g in mod.genes if g not in orphan_genes]
    
    return growth_rate, mod

#hp_tr_4 = growth_reconstruction(P_T_04, df_hp_tr_4, reaction_matrix)
#hp_tr_24 = growth_reconstruction(P_T_24, df_hp_tr_24, reaction_matrix)
#hp_tr_48 = growth_reconstruction(P_T_48, df_hp_tr_48, reaction_matrix)
#hp_tr_72 = growth_reconstruction(P_T_72, df_hp_tr_72, reaction_matrix)

hp_mo_4 = growth_reconstruction(P_M_04, df_hp_mo_4, reaction_matrix)
hp_mo_24 = growth_reconstruction(P_M_24, df_hp_mo_24, reaction_matrix)
hp_mo_48 = growth_reconstruction(P_M_48, df_hp_mo_48, reaction_matrix)
hp_mo_72 = growth_reconstruction(P_M_72, df_hp_mo_72, reaction_matrix)

# lp_tr_4 = growth_reconstruction(I_T_04, df_lp_tr_4, reaction_matrix)
# lp_tr_24 = growth_reconstruction(I_T_24, df_lp_tr_24, reaction_matrix)
# lp_tr_48 = growth_reconstruction(I_T_48, df_lp_tr_48, reaction_matrix)
# lp_tr_72 = growth_reconstruction(I_T_72, df_lp_tr_72, reaction_matrix)

# lp_mo_4 = growth_reconstruction(I_M_04, df_lp_mo_4, reaction_matrix)
# lp_mo_24 = growth_reconstruction(I_M_24, df_lp_mo_24, reaction_matrix)
# lp_mo_48 = growth_reconstruction(I_M_48, df_lp_mo_48, reaction_matrix)
# lp_mo_72 = growth_reconstruction(I_M_72, df_lp_mo_72, reaction_matrix)

#hp_tr_4_model = hp_tr_4[1]
#hp_tr_24_model = hp_tr_24[1]
#hp_tr_48_model = hp_tr_48[1]
#hp_tr_72_model = hp_tr_72[1]

hp_mo_4_model = hp_mo_4[1]
hp_mo_24_model = hp_mo_24[1]
hp_mo_48_model = hp_mo_48[1]
hp_mo_72_model = hp_mo_72[1]

# lp_tr_4_model = lp_tr_4[1]
# lp_tr_24_model = lp_tr_24[1]
# lp_tr_48_model = lp_tr_48[1]
# lp_tr_72_model = lp_tr_72[1]

# lp_mo_4_model = lp_mo_4[1]
# lp_mo_24_model = lp_mo_24[1]
# lp_mo_48_model = lp_mo_48[1]
# lp_mo_72_model = lp_mo_72[1]


print("Reconstruction done.")

# ## constrain models with mus and qps

df = pd.read_csv("/home/users/lzehetner/data/paper4_aav/HEK_AAV/mus_qps.csv", index_col=0)

reactions = {
    "mu": "MAR13082"
#    "qp_full": "MAR_SINK",
#    "qp_empty": "MAR_SINK_EMPTY"
}

def mu_qp_constraining(model, column):
    """
    Apply reaction constraints to a metabolic model based on a column from the DataFrame.

    Parameters:
        model (cobra.Model): The metabolic model to constrain.
        column (str): The column name from the DataFrame to use for constraints.
    """
    if column not in df.columns:
        raise ValueError(f"Column '{column}' not found in DataFrame.")

    model.reactions.MAR13082.upper_bound = 1000
    model.reactions.MAR13082.lower_bound = 0
    
#    for index, reaction_id in reactions.items():
#        if reaction_id in model.reactions:
#            value = df.loc[index, column]
#            reaction = model.reactions.get_by_id(reaction_id)
#            reaction.lower_bound = value
#            reaction.upper_bound = value
#            print(f"Set bounds for reaction {reaction_id} in model '{column}': {value}")

    mu_pred = model.optimize().objective_value
    print(mu_pred)

    model.reactions.MAR13082.upper_bound = mu_pred
    model.reactions.MAR13082.lower_bound = mu_pred

    return model

## Correct for loops in flux distributions
def calculate_loopless_fluxes(model, flux_distributions_df):
    # Initialize an empty list to store loopless flux distributions
    loopless_fluxes = []

    solver_status = []
    
    # Iterate through each row (flux distribution) in the dataframe
    for _, row in flux_distributions_df.iterrows():
        # Convert the row into a dictionary (reaction_id: flux_value)
        flux_distribution = row.to_dict()

        # Calculate the loopless flux distribution
        loopless_flux_distribution = loopless_solution(model, flux_distribution)

        # Append the loopless fluxes as a new row
        loopless_fluxes.append(loopless_flux_distribution.fluxes)

        solver_status.append(loopless_flux_distribution)

    # Convert the list of loopless flux dictionaries back to a DataFrame
    loopless_fluxes_df = pd.DataFrame(loopless_fluxes)

    solver_status_df = pd.DataFrame(solver_status)
    
    return loopless_fluxes_df, solver_status_df

# HP MO 48
selected_column = "hp_mo_48"
cstr_hp_mo_48_model = mu_qp_constraining(hp_mo_48_model, selected_column)
cstr_hp_mo_48_model.solver = "cplex"
cstr_hp_mo_48_model.solver.configuration.tolerances.feasibility = 1e-8

print("Flux sampling (HP-MO-48) started...")

optgp = OptGPSampler(cstr_hp_mo_48_model, processes=10)

flux_samples = pd.DataFrame()

target_samples = 10000

start_time = time.time()
#optgp = OptGPSampler(cstr_hp_tr_4_model, processes=50)
while len(flux_samples) < target_samples:
    try:
        # Number of samples to collect in each iteration
        remaining_samples = target_samples - len(flux_samples)
        sample_count = min(remaining_samples, 10)  # Collect up to 10 samples per iteration
        s2 = optgp.sample(sample_count)

        # Convert the samples to a DataFrame
        sample_df = pd.DataFrame(s2, columns=cstr_hp_mo_48_model.reactions.list_attr("id"))

        # Append the valid samples to the main DataFrame
        flux_samples = pd.concat([flux_samples, sample_df], ignore_index=True)

    except RuntimeError as e:
        # Handle the numerically unstable error without breaking the loop
        print(f"Encountered an issue: {e}. Skipping this batch of samples.")
    except Exception as e:
        # Handle any other unexpected errors
        print(f"Unexpected error: {e}.")
        break  # Exit the loop on critical errors
end_time = time.time()
elapsed_time = end_time - start_time
print(elapsed_time)

print("Flux sampling done.")
print("Loop correction starts ...")

start_time = time.time()
results = calculate_loopless_fluxes(cstr_hp_mo_48_model, flux_samples)
end_time = time.time()
elapsed_time = end_time - start_time
print(elapsed_time)
flux_distributions_df = results[0]
solver_results_df = results[1]

print("Loop correction done.")

flux_samples.to_csv("/home/users/lzehetner/data/paper4_aav/HEK_AAV/flux_sampling/flux_samples_hp_mo_48.csv", index=False)
flux_distributions_df.to_csv("/home/users/lzehetner/data/paper4_aav/HEK_AAV/flux_sampling/loopless_flux_samples_hp_mo_48.csv", index = False)
solver_results_df.to_csv("/home/users/lzehetner/data/paper4_aav/HEK_AAV/flux_sampling/solver_status_hp_mo_48.csv", index = False)

# HP MO 72
selected_column = "hp_mo_72"
cstr_hp_mo_72_model = mu_qp_constraining(hp_mo_72_model, selected_column)
cstr_hp_mo_72_model.solver = "cplex"
cstr_hp_mo_72_model.solver.configuration.tolerances.feasibility = 1e-8

print("Flux sampling (HP-MO-72) started...")

optgp = OptGPSampler(cstr_hp_mo_72_model, processes=10)

flux_samples = pd.DataFrame()

target_samples = 10000

start_time = time.time()
#optgp = OptGPSampler(cstr_hp_tr_4_model, processes=50)
while len(flux_samples) < target_samples:
    try:
        # Number of samples to collect in each iteration
        remaining_samples = target_samples - len(flux_samples)
        sample_count = min(remaining_samples, 10)  # Collect up to 10 samples per iteration
        s2 = optgp.sample(sample_count)

        # Convert the samples to a DataFrame
        sample_df = pd.DataFrame(s2, columns=cstr_hp_mo_72_model.reactions.list_attr("id"))

        # Append the valid samples to the main DataFrame
        flux_samples = pd.concat([flux_samples, sample_df], ignore_index=True)

    except RuntimeError as e:
        # Handle the numerically unstable error without breaking the loop
        print(f"Encountered an issue: {e}. Skipping this batch of samples.")
    except Exception as e:
        # Handle any other unexpected errors
        print(f"Unexpected error: {e}.")
        break  # Exit the loop on critical errors
end_time = time.time()
elapsed_time = end_time - start_time
print(elapsed_time)

print("Flux sampling done.")
print("Loop correction starts ...")

start_time = time.time()
results = calculate_loopless_fluxes(cstr_hp_mo_72_model, flux_samples)
end_time = time.time()
elapsed_time = end_time - start_time
print(elapsed_time)
flux_distributions_df = results[0]
solver_results_df = results[1]

print("Loop correction done.")

flux_samples.to_csv("/home/users/lzehetner/data/paper4_aav/HEK_AAV/flux_sampling/flux_samples_hp_mo_72.csv", index=False)
flux_distributions_df.to_csv("/home/users/lzehetner/data/paper4_aav/HEK_AAV/flux_sampling/loopless_flux_samples_hp_mo_72.csv", index = False)
solver_results_df.to_csv("/home/users/lzehetner/data/paper4_aav/HEK_AAV/flux_sampling/solver_status_hp_mo_72.csv", index = False)

# HP MO 24
selected_column = "hp_mo_24"
cstr_hp_mo_24_model = mu_qp_constraining(hp_mo_24_model, selected_column)
cstr_hp_mo_24_model.solver = "cplex"
cstr_hp_mo_24_model.solver.configuration.tolerances.feasibility = 1e-8

print("Flux sampling (HP-MO-24) started...")

optgp = OptGPSampler(cstr_hp_mo_24_model, processes=10)

flux_samples = pd.DataFrame()

target_samples = 10000

start_time = time.time()
#optgp = OptGPSampler(cstr_hp_tr_4_model, processes=50)
while len(flux_samples) < target_samples:
    try:
        # Number of samples to collect in each iteration
        remaining_samples = target_samples - len(flux_samples)
        sample_count = min(remaining_samples, 10)  # Collect up to 10 samples per iteration
        s2 = optgp.sample(sample_count)

        # Convert the samples to a DataFrame
        sample_df = pd.DataFrame(s2, columns=cstr_hp_mo_24_model.reactions.list_attr("id"))

        # Append the valid samples to the main DataFrame
        flux_samples = pd.concat([flux_samples, sample_df], ignore_index=True)

    except RuntimeError as e:
        # Handle the numerically unstable error without breaking the loop
        print(f"Encountered an issue: {e}. Skipping this batch of samples.")
    except Exception as e:
        # Handle any other unexpected errors
        print(f"Unexpected error: {e}.")
        break  # Exit the loop on critical errors
end_time = time.time()
elapsed_time = end_time - start_time
print(elapsed_time)

print("Flux sampling done.")
print("Loop correction starts ...")

start_time = time.time()
results = calculate_loopless_fluxes(cstr_hp_mo_24_model, flux_samples)
end_time = time.time()
elapsed_time = end_time - start_time
print(elapsed_time)
flux_distributions_df = results[0]
solver_results_df = results[1]

print("Loop correction done.")

flux_samples.to_csv("/home/users/lzehetner/data/paper4_aav/HEK_AAV/flux_sampling/flux_samples_hp_mo_24.csv", index=False)
flux_distributions_df.to_csv("/home/users/lzehetner/data/paper4_aav/HEK_AAV/flux_sampling/loopless_flux_samples_hp_mo_24.csv", index = False)
solver_results_df.to_csv("/home/users/lzehetner/data/paper4_aav/HEK_AAV/flux_sampling/solver_status_hp_mo_24.csv", index = False)


# HP MO 4
selected_column = "hp_mo_4"
cstr_hp_mo_4_model = mu_qp_constraining(hp_mo_4_model, selected_column)
cstr_hp_mo_4_model.solver = "cplex"
cstr_hp_mo_4_model.solver.configuration.tolerances.feasibility = 1e-8

print("Flux sampling (HP-MO-4) started...")

optgp = OptGPSampler(cstr_hp_mo_4_model, processes=10)

flux_samples = pd.DataFrame()

target_samples = 10000

start_time = time.time()
#optgp = OptGPSampler(cstr_hp_tr_4_model, processes=50)
while len(flux_samples) < target_samples:
    try:
        # Number of samples to collect in each iteration
        remaining_samples = target_samples - len(flux_samples)
        sample_count = min(remaining_samples, 10)  # Collect up to 10 samples per iteration
        s2 = optgp.sample(sample_count)

        # Convert the samples to a DataFrame
        sample_df = pd.DataFrame(s2, columns=cstr_hp_mo_4_model.reactions.list_attr("id"))

        # Append the valid samples to the main DataFrame
        flux_samples = pd.concat([flux_samples, sample_df], ignore_index=True)

    except RuntimeError as e:
        # Handle the numerically unstable error without breaking the loop
        print(f"Encountered an issue: {e}. Skipping this batch of samples.")
    except Exception as e:
        # Handle any other unexpected errors
        print(f"Unexpected error: {e}.")
        break  # Exit the loop on critical errors
end_time = time.time()
elapsed_time = end_time - start_time
print(elapsed_time)

print("Flux sampling done.")
print("Loop correction starts ...")
start_time = time.time()
results = calculate_loopless_fluxes(cstr_hp_mo_4_model, flux_samples)
end_time = time.time()
elapsed_time = end_time - start_time
print(elapsed_time)
flux_distributions_df = results[0]
solver_results_df = results[1]

print("Loop correction done.")

flux_samples.to_csv("/home/users/lzehetner/data/paper4_aav/HEK_AAV/flux_sampling/flux_samples_hp_mo_4.csv", index=False)
flux_distributions_df.to_csv("/home/users/lzehetner/data/paper4_aav/HEK_AAV/flux_sampling/loopless_flux_samples_hp_mo_4.csv", index = False)
solver_results_df.to_csv("/home/users/lzehetner/data/paper4_aav/HEK_AAV/flux_sampling/solver_status_hp_mo_4.csv", index = False)



