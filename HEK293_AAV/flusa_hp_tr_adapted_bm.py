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

P_T_04 = ["hp_v1_tp2", "hp_v2_tp2", "hp_v3_tp2", "hp_v4_tp2"]
P_T_24 = ["hp_v1_tp3", "hp_v2_tp3", "hp_v3_tp3", "hp_v4_tp3"]
P_T_48 = ["hp_v1_tp4", "hp_v2_tp4", "hp_v3_tp4", "hp_v4_tp4"]
P_T_72 = ["hp_v1_tp5", "hp_v2_tp5", "hp_v3_tp5", "hp_v4_tp5"]

# P_M_04 = ["hp_v5_tp2", "hp_v6_tp2", "hp_v7_tp2", "hp_v8_tp2"]
# P_M_24 = ["hp_v5_tp3", "hp_v6_tp3", "hp_v7_tp3", "hp_v8_tp3"]
# P_M_48 = ["hp_v5_tp4", "hp_v6_tp4", "hp_v7_tp4", "hp_v8_tp4"]
# P_M_72 = ["hp_v5_tp5", "hp_v6_tp5", "hp_v8_tp5"]

input_rates = "/home/users/lzehetner/data/hek/specific_exchange_rates.csv"
df = pd.read_csv(input_rates)

filtered_dfs = {}

states = ['HP_TR', 'HP_MO', 'LP_TR', 'LP_MO']
time_points = [4, 24, 48, 72]

for state in states:
    for time_point in time_points:
        key = f"{state}_{time_point}"
        filtered_dfs[key] = df[(df['state'] == state) & (df['time_point'] == time_point)]


df_hp_tr_4 = filtered_dfs['HP_TR_4']
df_hp_tr_24 = filtered_dfs['HP_TR_24']
df_hp_tr_48 = filtered_dfs['HP_TR_48']
df_hp_tr_72 = filtered_dfs['HP_TR_72']

# df_hp_mo_4 = filtered_dfs['HP_MO_4']
# df_hp_mo_24 = filtered_dfs['HP_MO_24']
# df_hp_mo_48 = filtered_dfs['HP_MO_48']
# df_hp_mo_72 = filtered_dfs['HP_MO_72']

# df_lp_tr_4 = filtered_dfs['LP_TR_4']
# df_lp_tr_24 = filtered_dfs['LP_TR_24']
# df_lp_tr_48 = filtered_dfs['LP_TR_48']
# df_lp_tr_72 = filtered_dfs['LP_TR_72']

# df_lp_mo_4 = filtered_dfs['LP_MO_4']
# df_lp_mo_24 = filtered_dfs['LP_MO_24']
# df_lp_mo_48 = filtered_dfs['LP_MO_48']
# df_lp_mo_72 = filtered_dfs['LP_MO_72']

# remove last row with ammonia, since there are no exchange rate fitted
df_hp_tr_72 = df_hp_tr_72[:-1]

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

# # AAV production reaction
# #### Alanine reaction

ala_reaction = Reaction('MAR_ALA_nuc')
ala_reaction.name = 'Exchange reaction of alanine between cytosole and nucleus'
ala_reaction.subsystem = ''
ala_reaction.lower_bound = 0 # but needs to be adopted from paper
ala_reaction.upper_bound = 1000 # but needs also to be adopted

MAM01307n = Metabolite(
    'MAM01307n',
    formula= mod.metabolites.get_by_id("MAM01307c").formula,
    name='alanine',
    compartment='n')

ala_reaction.add_metabolites({
    mod.metabolites.get_by_id("MAM01307c"): -1, # Alanine
    MAM01307n: 1
})

# #### Arginine reaction

arg_reaction = Reaction('MAR_ARG_nuc')
arg_reaction.name = 'Exchange reaction of arginine between cytosole and nucleus'
arg_reaction.subsystem = ''
arg_reaction.lower_bound = 0 # but needs to be adopted from paper
arg_reaction.upper_bound = 1000 # but needs also to be adopted

MAM01365n = Metabolite(
    'MAM01365n',
    formula= mod.metabolites.get_by_id("MAM01365c").formula,
    name='arginine',
    compartment='n')

arg_reaction.add_metabolites({
    mod.metabolites.get_by_id("MAM01365c"): -1, # Alanine
    MAM01365n: 1
})

# #### Asparagine reaction

asn_reaction = Reaction('MAR_ASN_nuc')
asn_reaction.name = 'Exchange reaction of asparagine between cytosole and nucleus'
asn_reaction.subsystem = ''
asn_reaction.lower_bound = 0 # but needs to be adopted from paper
asn_reaction.upper_bound = 1000 # but needs also to be adopted

MAM01369n = Metabolite(
    'MAM01369n',
    formula= mod.metabolites.get_by_id("MAM01369c").formula,
    name='asparagine',
    compartment='n')

asn_reaction.add_metabolites({
    mod.metabolites.get_by_id("MAM01369c"): -1, # Alanine
    MAM01369n: 1
})

# #### Aspartate reaction

asp_reaction = Reaction('MAR_ASP_nuc')
asp_reaction.name = 'Exchange reaction of aspartate between cytosole and nucleus'
asp_reaction.subsystem = ''
asp_reaction.lower_bound = 0 # but needs to be adopted from paper
asp_reaction.upper_bound = 1000 # but needs also to be adopted

MAM01370n = Metabolite(
    'MAM01370n',
    formula= mod.metabolites.get_by_id("MAM01370c").formula,
    name='aspartate',
    compartment='n')

asp_reaction.add_metabolites({
    mod.metabolites.get_by_id("MAM01370c"): -1, # Alanine
    MAM01370n: 1
})

# #### Cysteine reaction

cys_reaction = Reaction('MAR_CYS_nuc')
cys_reaction.name = 'Exchange reaction of cysteine between cytosole and nucleus'
cys_reaction.subsystem = ''
cys_reaction.lower_bound = 0 # but needs to be adopted from paper
cys_reaction.upper_bound = 1000 # but needs also to be adopted

MAM01628n = Metabolite(
    'MAM01628n',
    formula= mod.metabolites.get_by_id("MAM01628c").formula,
    name='cysteine',
    compartment='n')

cys_reaction.add_metabolites({
    mod.metabolites.get_by_id("MAM01628c"): -1, # Alanine
    MAM01628n: 1
})

# #### Glutamine reaction

gln_reaction = Reaction('MAR_GLN_nuc')
gln_reaction.name = 'Exchange reaction of glutamine between cytosole and nucleus'
gln_reaction.subsystem = ''
gln_reaction.lower_bound = 0 # but needs to be adopted from paper
gln_reaction.upper_bound = 1000 # but needs also to be adopted

MAM01975n = Metabolite(
    'MAM01975n',
    formula= mod.metabolites.get_by_id("MAM01975c").formula,
    name='glutamine',
    compartment='n')

gln_reaction.add_metabolites({
    mod.metabolites.get_by_id("MAM01975c"): -1, # Alanine
    MAM01975n: 1
})

# #### Glutamate reaction

glu_reaction = Reaction('MAR_GLU_nuc')
glu_reaction.name = 'Exchange reaction of glutamate between cytosole and nucleus'
glu_reaction.subsystem = ''
glu_reaction.lower_bound = 0 # but needs to be adopted from paper
glu_reaction.upper_bound = 1000 # but needs also to be adopted

MAM01974n = Metabolite(
    'MAM01974n',
    formula= mod.metabolites.get_by_id("MAM01974c").formula,
    name='glutamate',
    compartment='n')

glu_reaction.add_metabolites({
    mod.metabolites.get_by_id("MAM01974c"): -1, # Alanine
    MAM01974n: 1
})

# #### Glycine reaction

gly_reaction = Reaction('MAR_GLY_nuc')
gly_reaction.name = 'Exchange reaction of glycine between cytosole and nucleus'
gly_reaction.subsystem = ''
gly_reaction.lower_bound = 0 # but needs to be adopted from paper
gly_reaction.upper_bound = 1000 # but needs also to be adopted

MAM01986n = Metabolite(
    'MAM01986n',
    formula= mod.metabolites.get_by_id("MAM01986c").formula,
    name='glycine',
    compartment='n')

gly_reaction.add_metabolites({
    mod.metabolites.get_by_id("MAM01986c"): -1, # Alanine
    MAM01986n: 1
})

# #### Histidine reaction

his_reaction = Reaction('MAR_HIS_nuc')
his_reaction.name = 'Exchange reaction of histidine between cytosole and nucleus'
his_reaction.subsystem = ''
his_reaction.lower_bound = 0 # but needs to be adopted from paper
his_reaction.upper_bound = 1000 # but needs also to be adopted

MAM02125n = Metabolite(
    'MAM02125n',
    formula= mod.metabolites.get_by_id("MAM02125c").formula,
    name='histidine',
    compartment='n')

his_reaction.add_metabolites({
    mod.metabolites.get_by_id("MAM02125c"): -1, # Alanine
    MAM02125n: 1
})

# #### Isoleucine reaction

ile_reaction = Reaction('MAR_ILE_nuc')
ile_reaction.name = 'Exchange reaction of isoleucine between cytosole and nucleus'
ile_reaction.subsystem = ''
ile_reaction.lower_bound = 0 # but needs to be adopted from paper
ile_reaction.upper_bound = 1000 # but needs also to be adopted

MAM02184n = Metabolite(
    'MAM02184n',
    formula= mod.metabolites.get_by_id("MAM02184c").formula,
    name='isoleucine',
    compartment='n')

ile_reaction.add_metabolites({
    mod.metabolites.get_by_id("MAM02184c"): -1, # Alanine
    MAM02184n: 1
})

# #### Leucine reaction

leu_reaction = Reaction('MAR_LEU_nuc')
leu_reaction.name = 'Exchange reaction of leucine between cytosole and nucleus'
leu_reaction.subsystem = ''
leu_reaction.lower_bound = 0 # but needs to be adopted from paper
leu_reaction.upper_bound = 1000 # but needs also to be adopted

MAM02360n = Metabolite(
    'MAM02360n',
    formula= mod.metabolites.get_by_id("MAM02360c").formula,
    name='leucine',
    compartment='n')

leu_reaction.add_metabolites({
    mod.metabolites.get_by_id("MAM02360c"): -1, # Alanine
    MAM02360n: 1
})

# #### Lysine reaction

lys_reaction = Reaction('MAR_LYS_nuc')
lys_reaction.name = 'Exchange reaction of lysine between cytosole and nucleus'
lys_reaction.subsystem = ''
lys_reaction.lower_bound = 0 # but needs to be adopted from paper
lys_reaction.upper_bound = 1000 # but needs also to be adopted

MAM02426n = Metabolite(
    'MAM02426n',
    formula= mod.metabolites.get_by_id("MAM02426c").formula,
    name='lysine',
    compartment='n')

lys_reaction.add_metabolites({
    mod.metabolites.get_by_id("MAM02426c"): -1,
    MAM02426n: 1
})

# #### Methionine reaction

met_reaction = Reaction('MAR_MET_nuc')
met_reaction.name = 'Exchange reaction of methionine between cytosole and nucleus'
met_reaction.subsystem = ''
met_reaction.lower_bound = 0 # but needs to be adopted from paper
met_reaction.upper_bound = 1000 # but needs also to be adopted

MAM02471n = Metabolite(
    'MAM02471n',
    formula= mod.metabolites.get_by_id("MAM02471c").formula,
    name='methionine',
    compartment='n')

met_reaction.add_metabolites({
    mod.metabolites.get_by_id("MAM02471c"): -1,
    MAM02471n: 1
})

# #### Phenylalanine reaction

phe_reaction = Reaction('MAR_PHE_nuc')
phe_reaction.name = 'Exchange reaction of phenylalanine between cytosole and nucleus'
phe_reaction.subsystem = ''
phe_reaction.lower_bound = 0 # but needs to be adopted from paper
phe_reaction.upper_bound = 1000 # but needs also to be adopted

MAM02724n = Metabolite(
    'MAM02724n',
    formula= mod.metabolites.get_by_id("MAM02724c").formula,
    name='phenylalanine',
    compartment='n')

phe_reaction.add_metabolites({
    mod.metabolites.get_by_id("MAM02724c"): -1,
    MAM02724n: 1
})

# #### Proline reaction

pro_reaction = Reaction('MAR_PRO_nuc')
pro_reaction.name = 'Exchange reaction of proline between cytosole and nucleus'
pro_reaction.subsystem = ''
pro_reaction.lower_bound = 0 # but needs to be adopted from paper
pro_reaction.upper_bound = 1000 # but needs also to be adopted

MAM02770n = Metabolite(
    'MAM02770n',
    formula= mod.metabolites.get_by_id("MAM02770c").formula,
    name='proline',
    compartment='n')

pro_reaction.add_metabolites({
    mod.metabolites.get_by_id("MAM02770c"): -1,
    MAM02770n: 1
})

# #### Serine reaction

ser_reaction = Reaction('MAR_SER_nuc')
ser_reaction.name = 'Exchange reaction of serine between cytosole and nucleus'
ser_reaction.subsystem = ''
ser_reaction.lower_bound = 0 # but needs to be adopted from paper
ser_reaction.upper_bound = 1000 # but needs also to be adopted

MAM02896n = Metabolite(
    'MAM02896n',
    formula= mod.metabolites.get_by_id("MAM02896c").formula,
    name='serine',
    compartment='n')

ser_reaction.add_metabolites({
    mod.metabolites.get_by_id("MAM02896c"): -1,
    MAM02896n: 1
})

# #### Threonine reaction

thr_reaction = Reaction('MAR_THR_nuc')
thr_reaction.name = 'Exchange reaction of threonine between cytosole and nucleus'
thr_reaction.subsystem = ''
thr_reaction.lower_bound = 0 # but needs to be adopted from paper
thr_reaction.upper_bound = 1000 # but needs also to be adopted

MAM02993n = Metabolite(
    'MAM02993n',
    formula= mod.metabolites.get_by_id("MAM02993c").formula,
    name='threonine',
    compartment='n')

thr_reaction.add_metabolites({
    mod.metabolites.get_by_id("MAM02993c"): -1,
    MAM02993n: 1
})

# #### Tryptophane reaction

trp_reaction = Reaction('MAR_TRP_nuc')
trp_reaction.name = 'Exchange reaction of tryptophane between cytosole and nucleus'
trp_reaction.subsystem = ''
trp_reaction.lower_bound = 0 # but needs to be adopted from paper
trp_reaction.upper_bound = 1000 # but needs also to be adopted

MAM03089n = Metabolite(
    'MAM03089n',
    formula= mod.metabolites.get_by_id("MAM03089c").formula,
    name='tryptophane',
    compartment='n')

trp_reaction.add_metabolites({
    mod.metabolites.get_by_id("MAM03089c"): -1,
    MAM03089n: 1
})

# #### Tyrosine reaction

tyr_reaction = Reaction('MAR_TYR_nuc')
tyr_reaction.name = 'Exchange reaction of tyrosine between cytosole and nucleus'
tyr_reaction.subsystem = ''
tyr_reaction.lower_bound = 0 # but needs to be adopted from paper
tyr_reaction.upper_bound = 1000 # but needs also to be adopted

MAM03101n = Metabolite(
    'MAM03101n',
    formula= mod.metabolites.get_by_id("MAM03101c").formula,
    name='tyrosine',
    compartment='n')

tyr_reaction.add_metabolites({
    mod.metabolites.get_by_id("MAM03101c"): -1,
    MAM03101n: 1
})

# #### Valine reaction

val_reaction = Reaction('MAR_VAL_nuc')
val_reaction.name = 'Exchange reaction of valine between cytosole and nucleus'
val_reaction.subsystem = ''
val_reaction.lower_bound = 0 # but needs to be adopted from paper
val_reaction.upper_bound = 1000 # but needs also to be adopted

MAM03135n = Metabolite(
    'MAM03135n',
    formula= mod.metabolites.get_by_id("MAM03135c").formula,
    name='valine',
    compartment='n')

val_reaction.add_metabolites({
    mod.metabolites.get_by_id("MAM03135c"): -1,
    MAM03135n: 1
})

# the file was downloaded from ncbi db
infile = open("/home/users/lzehetner/sequence_1.fasta")

aav_dict = SeqIO.to_dict(SeqIO.parse(infile, "fasta"))
# print(aav_dict)
vp1_seq = aav_dict["YP_077180.1"].seq
vp2_seq = aav_dict["YP_077180.2"].seq
vp3_seq = aav_dict["YP_077180.3"].seq

vp1_counts = {
    "A": vp1_seq.count('A'),
    "R": vp1_seq.count('R'),
    "N": vp1_seq.count('N'),
    "D": vp1_seq.count('D'),
    "C": vp1_seq.count('C'),
    "Q": vp1_seq.count('Q'),
    "E": vp1_seq.count('E'),
    "G": vp1_seq.count('G'),
    "H": vp1_seq.count('H'),
    "I": vp1_seq.count('I'),
    "L": vp1_seq.count('L'),
    "K": vp1_seq.count('K'),
    "M": vp1_seq.count('M'),
    "F": vp1_seq.count('F'),
    "P": vp1_seq.count('P'),
    "S": vp1_seq.count('S'),
    "T": vp1_seq.count('T'),
    "W": vp1_seq.count('W'),
    "Y": vp1_seq.count('Y'),
    "V": vp1_seq.count('V'),
}

vp2_counts = {
    "A": vp2_seq.count('A'),
    "R": vp2_seq.count('R'),
    "N": vp2_seq.count('N'),
    "D": vp2_seq.count('D'),
    "C": vp2_seq.count('C'),
    "Q": vp2_seq.count('Q'),
    "E": vp2_seq.count('E'),
    "G": vp2_seq.count('G'),
    "H": vp2_seq.count('H'),
    "I": vp2_seq.count('I'),
    "L": vp2_seq.count('L'),
    "K": vp2_seq.count('K'),
    "M": vp2_seq.count('M'),
    "F": vp2_seq.count('F'),
    "P": vp2_seq.count('P'),
    "S": vp2_seq.count('S'),
    "T": vp2_seq.count('T'),
    "W": vp2_seq.count('W'),
    "Y": vp2_seq.count('Y'),
    "V": vp2_seq.count('V'),
}

vp3_counts = {
    "A": vp3_seq.count('A'),
    "R": vp3_seq.count('R'),
    "N": vp3_seq.count('N'),
    "D": vp3_seq.count('D'),
    "C": vp3_seq.count('C'),
    "Q": vp3_seq.count('Q'),
    "E": vp3_seq.count('E'),
    "G": vp3_seq.count('G'),
    "H": vp3_seq.count('H'),
    "I": vp3_seq.count('I'),
    "L": vp3_seq.count('L'),
    "K": vp3_seq.count('K'),
    "M": vp3_seq.count('M'),
    "F": vp3_seq.count('F'),
    "P": vp3_seq.count('P'),
    "S": vp3_seq.count('S'),
    "T": vp3_seq.count('T'),
    "W": vp3_seq.count('W'),
    "Y": vp3_seq.count('Y'),
    "V": vp3_seq.count('V'),
}

print("VP1: ", vp1_counts)
print("VP2: ", vp2_counts)
print("VP3: ", vp3_counts)

aa_compositions = {
    'A': {'C': 3, 'H': 7, 'O': 2, 'N': 1, 'S': 0},
    'C': {'C': 3, 'H': 7, 'O': 2, 'N': 1, 'S': 1},
    'D': {'C': 4, 'H': 7, 'O': 4, 'N': 1, 'S': 0},
    'E': {'C': 5, 'H': 9, 'O': 4, 'N': 1, 'S': 0},
    'F': {'C': 9, 'H': 11, 'O': 2, 'N': 1, 'S': 0},
    'G': {'C': 2, 'H': 5, 'O': 2, 'N': 1, 'S': 0},
    'H': {'C': 6, 'H': 9, 'O': 2, 'N': 3, 'S': 0},
    'I': {'C': 6, 'H': 13,'O': 2, 'N': 1, 'S': 0},
    'K': {'C': 6, 'H': 14,'O': 2, 'N': 2, 'S': 0},
    'L': {'C': 6, 'H': 13,'O': 2, 'N': 1, 'S': 0},
    'M': {'C': 5, 'H': 11, 'O': 2, 'N': 1, 'S': 1},
    'N': {'C': 4, 'H': 8, 'O': 3, 'N': 2, 'S': 0},
    'P': {'C': 5, 'H': 9, 'O': 2, 'N': 1, 'S': 0},
    'Q': {'C': 5, 'H': 10, 'O': 3, 'N': 2, 'S': 0},
    'R': {'C': 6, 'H': 14,'O': 2, 'N': 4, 'S': 0},
    'S': {'C': 3, 'H': 7, 'O': 3, 'N': 1, 'S': 0},
    'T': {'C': 4, 'H': 9, 'O': 3, 'N': 1, 'S': 0},
    'V': {'C': 5, 'H': 11, 'O': 2, 'N': 1, 'S': 0},
    'W': {'C': 11,'H': 12,'O': 2, 'N': 2, 'S': 0},
    'Y': {'C': 9, 'H': 11, 'O': 3, 'N': 1, 'S': 0},
}

# VP1
protein_seq = str(vp1_seq)
analysis = ProteinAnalysis(protein_seq)

# Count amino acids in the protein sequence
aa_counts = analysis.count_amino_acids()

# Calculate the sum formula based on atomic composition
sum_formula_counts = {'C': 0, 'H': 0, 'O': 0, 'N': 0, 'S': 0}

for aa, count in aa_counts.items():
    for atom, atom_count in aa_compositions[aa].items():
        sum_formula_counts[atom] += atom_count * count

# Adjust for the released water molecules during peptide bond formation
num_peptide_bonds = len(protein_seq) - 1
sum_formula_counts['H'] -= 2 * num_peptide_bonds
sum_formula_counts['O'] -= num_peptide_bonds

# Create the sum formula string
sum_formula_vp1 = "".join([f"{atom}{count}" for atom, count in sum_formula_counts.items()])
print(f"Sum Formula: {sum_formula_vp1}")

# VP2
protein_seq = str(vp2_seq)
analysis = ProteinAnalysis(protein_seq)

# Count amino acids in the protein sequence
aa_counts = analysis.count_amino_acids()

# Calculate the sum formula based on atomic composition
sum_formula_counts = {'C': 0, 'H': 0, 'O': 0, 'N': 0, 'S': 0}

for aa, count in aa_counts.items():
    for atom, atom_count in aa_compositions[aa].items():
        sum_formula_counts[atom] += atom_count * count

# Adjust for the released water molecules during peptide bond formation
num_peptide_bonds = len(protein_seq) - 1
sum_formula_counts['H'] -= 2 * num_peptide_bonds
sum_formula_counts['O'] -= num_peptide_bonds

# Create the sum formula string
sum_formula_vp2 = "".join([f"{atom}{count}" for atom, count in sum_formula_counts.items()])
print(f"Sum Formula: {sum_formula_vp2}")

# VP3
protein_seq = str(vp3_seq)
analysis = ProteinAnalysis(protein_seq)

# Count amino acids in the protein sequence
aa_counts = analysis.count_amino_acids()

# Calculate the sum formula based on atomic composition
sum_formula_counts = {'C': 0, 'H': 0, 'O': 0, 'N': 0, 'S': 0}

for aa, count in aa_counts.items():
    for atom, atom_count in aa_compositions[aa].items():
        sum_formula_counts[atom] += atom_count * count

# Adjust for the released water molecules during peptide bond formation
num_peptide_bonds = len(protein_seq) - 1
sum_formula_counts['H'] -= 2 * num_peptide_bonds
sum_formula_counts['O'] -= num_peptide_bonds

# Create the sum formula string
sum_formula_vp3 = "".join([f"{atom}{count}" for atom, count in sum_formula_counts.items()])
print(f"Sum Formula: {sum_formula_vp3}")




protein_seq = 5 * str(vp1_seq) + 5 * str(vp2_seq) + 50 * str(vp3_seq)
analysis = ProteinAnalysis(protein_seq)

# Count amino acids in the protein sequence
aa_counts = analysis.count_amino_acids()

# Calculate the sum formula based on atomic composition
sum_formula_counts = {'C': 0, 'H': 0, 'O': 0, 'N': 0, 'S': 0}

for aa, count in aa_counts.items():
    for atom, atom_count in aa_compositions[aa].items():
        sum_formula_counts[atom] += atom_count * count

# Adjust for the released water molecules during peptide bond formation
num_peptide_bonds = len(protein_seq) - 1
sum_formula_counts['H'] -= 2 * num_peptide_bonds
sum_formula_counts['O'] -= num_peptide_bonds

# Create the sum formula string
sum_formula_capsid = "".join([f"{atom}{count}" for atom, count in sum_formula_counts.items()])
print(f"Sum Formula of Capsid: {sum_formula_capsid}")

# #### vp1 reaction

vp1_reaction = Reaction('MAR_VP1')
vp1_reaction.name = 'Synthesis of viral protein 1'
vp1_reaction.subsystem = ''
vp1_reaction.lower_bound = 0 # but needs to be adopted from paper
vp1_reaction.upper_bound = 1000 # but needs also to be adopted

VP1_n = Metabolite(
    'VP1_n',
    formula= {sum_formula_vp1},
    name='Viral Protein 1',
    compartment='n')

vp1_reaction.add_metabolites({
    MAM01307n: -vp1_counts["A"], # Alanine
    MAM01365n: -vp1_counts["R"], # Arginine
    MAM01369n: -vp1_counts["N"], # Asparagine
    MAM01370n: -vp1_counts["D"], # Aspartate
    MAM01628n: -vp1_counts["C"], # Cysteine
    MAM01975n: -vp1_counts["Q"], # Glutamine
    MAM01974n: -vp1_counts["E"], # Glutamate
    MAM01986n: -vp1_counts["G"], # Glycine
    MAM02125n: -vp1_counts["H"], # Histidine
    MAM02184n: -vp1_counts["I"], # Ile
    MAM02360n: -vp1_counts["L"], # Leu
    MAM02426n: -vp1_counts["K"], # Lysin
    MAM02471n: -vp1_counts["M"], # Met
    MAM02724n: -vp1_counts["F"], # Phe
    MAM02770n: -vp1_counts["P"], # Pro
    MAM02896n: -vp1_counts["S"], # Serin
    MAM02993n: -vp1_counts["T"], # Thr
    MAM03089n: -vp1_counts["W"], # Trp
    MAM03101n: -vp1_counts["Y"], # Tyr
    MAM03135n: -vp1_counts["V"], # Val
    mod.metabolites.get_by_id("MAM02034n"): - 3 * sum(vp1_counts.values()), # GTP
    mod.metabolites.get_by_id("MAM02040n"): -( 3 * sum(vp1_counts.values()) + 1), # Water
    mod.metabolites.get_by_id("MAM01948n"): 3 * sum(vp1_counts.values()), # GDP
    mod.metabolites.get_by_id("MAM02751n"): 3 * sum(vp1_counts.values()), # phosphate
    mod.metabolites.get_by_id("MAM02039n"): 0.75 * 3 * sum(vp1_counts.values()), # proton
    VP1_n: 1.0
})

##### vp2 reaction
vp2_reaction = Reaction('MAR_VP2')
vp2_reaction.name = 'Synthesis of viral protein 2'
vp2_reaction.subsystem = ''
vp2_reaction.lower_bound = 0 # but needs to be adopted from paper
vp2_reaction.upper_bound = 1000 # but needs also to be adopted

VP2_n = Metabolite(
    'VP2_n',
    formula= {sum_formula_vp2},
    name='Viral Protein 2',
    compartment='n')

vp2_reaction.add_metabolites({
    MAM01307n: -vp2_counts["A"], # Alanine
    MAM01365n: -vp2_counts["R"], # Arginine
    MAM01369n: -vp2_counts["N"], # Asparagine
    MAM01370n: -vp2_counts["D"], # Aspartate
    MAM01628n: -vp2_counts["C"], # Cysteine
    MAM01975n: -vp2_counts["Q"], # Glutamine
    MAM01974n: -vp2_counts["E"], # Glutamate
    MAM01986n: -vp2_counts["G"], # Glycine
    MAM02125n: -vp2_counts["H"], # Histidine
    MAM02184n: -vp2_counts["I"], # Ile
    MAM02360n: -vp2_counts["L"], # Leu
    MAM02426n: -vp2_counts["K"], # Lysin
    MAM02471n: -vp2_counts["M"], # Met
    MAM02724n: -vp2_counts["F"], # Phe
    MAM02770n: -vp2_counts["P"], # Pro
    MAM02896n: -vp2_counts["S"], # Serin
    MAM02993n: -vp2_counts["T"], # Thr
    MAM03089n: -vp2_counts["W"], # Trp
    MAM03101n: -vp2_counts["Y"], # Tyr
    MAM03135n: -vp2_counts["V"], # Val
    mod.metabolites.get_by_id("MAM02034n"): - 3 * sum(vp2_counts.values()), # GTP
    mod.metabolites.get_by_id("MAM02040n"): -( 3 * sum(vp2_counts.values()) + 1), # Water
    mod.metabolites.get_by_id("MAM01948n"): 3 * sum(vp2_counts.values()), # GDP
    mod.metabolites.get_by_id("MAM02751n"): 3 * sum(vp2_counts.values()), # phosphate
    mod.metabolites.get_by_id("MAM02039n"): 0.75 * 3 * sum(vp2_counts.values()), # proton
    VP2_n: 1.0
})

##### vp3 reaction
vp3_reaction = Reaction('MAR_VP3')
vp3_reaction.name = 'Synthesis of viral protein 3'
vp3_reaction.subsystem = ''
vp3_reaction.lower_bound = 0 # but needs to be adopted from paper
vp3_reaction.upper_bound = 1000 # but needs also to be adopted

VP3_n = Metabolite(
    'VP3_n',
    formula= {sum_formula_vp3},
    name='Viral Protein 3',
    compartment='n')

vp3_reaction.add_metabolites({
    MAM01307n: -vp3_counts["A"], # Alanine
    MAM01365n: -vp3_counts["R"], # Arginine
    MAM01369n: -vp3_counts["N"], # Asparagine
    MAM01370n: -vp3_counts["D"], # Aspartate
    MAM01628n: -vp3_counts["C"], # Cysteine
    MAM01975n: -vp3_counts["Q"], # Glutamine
    MAM01974n: -vp3_counts["E"], # Glutamate
    MAM01986n: -vp3_counts["G"], # Glycine
    MAM02125n: -vp3_counts["H"], # Histidine
    MAM02184n: -vp3_counts["I"], # Ile
    MAM02360n: -vp3_counts["L"], # Leu
    MAM02426n: -vp3_counts["K"], # Lysin
    MAM02471n: -vp3_counts["M"], # Met
    MAM02724n: -vp3_counts["F"], # Phe
    MAM02770n: -vp3_counts["P"], # Pro
    MAM02896n: -vp3_counts["S"], # Serin
    MAM02993n: -vp3_counts["T"], # Thr
    MAM03089n: -vp3_counts["W"], # Trp
    MAM03101n: -vp3_counts["Y"], # Tyr
    MAM03135n: -vp3_counts["V"], # Val
    mod.metabolites.get_by_id("MAM02034n"): - 3 * sum(vp3_counts.values()), # GTP
    mod.metabolites.get_by_id("MAM02040n"): -( 3 * sum(vp3_counts.values()) + 1), # Water
    mod.metabolites.get_by_id("MAM01948n"): 3 * sum(vp3_counts.values()), # GDP
    mod.metabolites.get_by_id("MAM02751n"): 3 * sum(vp3_counts.values()), # phosphate
    mod.metabolites.get_by_id("MAM02039n"): 0.75 * 3 * sum(vp3_counts.values()), # proton
    VP3_n: 1.0
})

# #### capsid reaction

capsid_reaction = Reaction('MAR_CAPSID')
capsid_reaction.name = 'Synthesis of viral Capsid'
capsid_reaction.subsystem = ''
capsid_reaction.lower_bound = 0 # but needs to be adopted from paper
capsid_reaction.upper_bound = 1000 # but needs also to be adopted

capsid_n = Metabolite(
    'Capsid',
    formula = {sum_formula_capsid},
    name = 'Capsid',
    compartment = 'n')

capsid_reaction.add_metabolites({
    VP1_n : -5,
    VP2_n : -5,
    VP3_n : -50,
    capsid_n : 1
})

##### product reaction
nucleotide_reaction_1 = Reaction('NT_SYNTH_1')
nucleotide_reaction_1.name = 'Synthesis of 47 nucleotide sequence'
nucleotide_reaction_1.subsystem = ''
nucleotide_reaction_1.lower_bound = 0
nucleotide_reaction_1.upper_bound = 1000

nucleotide_seq_1 = Metabolite(
    '47_Nucleotide_Sequence',
    formula = 'C458H716N176O411P47',
    name = '47-NT Building Block',
    compartment = 'n'
)

nucleotide_reaction_1.add_metabolites({
    mod.metabolites.get_by_id("MAM01688n") : -11.75,  # G
    mod.metabolites.get_by_id("MAM01642n") : -11.75,  # A
    mod.metabolites.get_by_id("MAM01753n") : -11.75,  # T
    mod.metabolites.get_by_id("MAM01645n") : -11.75,  # C
    nucleotide_seq_1 : 1,
    mod.metabolites.get_by_id("MAM02759n") : 46 # PPi
})

nucleotide_reaction_2 = Reaction('NT_SYNTH_2')
nucleotide_reaction_2.name = 'Synthesis of 470 nucleotide sequence'
nucleotide_reaction_2.subsystem = ''
nucleotide_reaction_2.lower_bound = 0
nucleotide_reaction_2.upper_bound = 1000

nucleotide_seq_2 = Metabolite(
    '470_Nucleotide_Sequence',
    formula = 'C4580H7160N1760O4110P470',
    name = '470-NT Building Block',
    compartment = 'n'
)

nucleotide_reaction_2.add_metabolites({
    nucleotide_seq_1 : -10,
    nucleotide_seq_2 : 1,
    mod.metabolites.get_by_id("MAM02759n") : 10 # PPi
})


target_gene_reaction = Reaction('MAR_TARGET_GENE')
target_gene_reaction.name = 'Synthesis of target Gene'
target_gene_reaction.subsystem = ''
target_gene_reaction.lower_bound = 0 # but needs to be adopted from paper
target_gene_reaction.upper_bound = 1000 # but needs also to be adopted

target_gene_n = Metabolite(
    'Target_Gene_Nucleus',
    formula = 'C45825H71600O41100N17600P4699',
    name = 'Target Gene',
    compartment = 'n')

target_gene_reaction.add_metabolites({
    nucleotide_seq_2 : -10, 
    target_gene_n : 1,
    mod.metabolites.get_by_id("MAM02759n") : 10 # PPi
})

# #### assembly reaction

assembly_reaction = Reaction('MAR_ASSEMBLY')
assembly_reaction.name = 'Assembly of Product'
assembly_reaction.subsystem = ''
assembly_reaction.lower_bound = 0 # but needs to be adopted from paper
assembly_reaction.upper_bound = 1000 # but needs also to be adopted

product_n = Metabolite(
    'Product_Nucleus',
    formula = 'C211920H306299O79399N63505P4699S905',
    name = 'Product Nucleus',
    compartment = 'n')

assembly_reaction.add_metabolites({
    target_gene_n : -1,
    capsid_n : -1,
    product_n : 1
})

# #### export of aav from nucleus to cytosole

prod_transport_reaction = Reaction('MAR_TRANSPORT')
prod_transport_reaction.name = 'Export of AAV to Cytosole'
prod_transport_reaction.subsystem = ''
prod_transport_reaction.lower_bound = 0 # but needs to be adopted from paper
prod_transport_reaction.upper_bound = 1000 # but needs also to be adopted

product_c = Metabolite(
    'Product_Cytosole',
    formula = 'C211920H306299O79399N63505P4699S905',
    name = 'Product Cytosole',
    compartment = 'c')

prod_transport_reaction.add_metabolites({
    product_n : -1,
    product_c : 1
})

# #### export of product

prod_export_reaction = Reaction('MAR_EXPORT')
prod_export_reaction.name = 'Export of Product'
prod_export_reaction.subsystem = ''
prod_export_reaction.lower_bound = 0 # but needs to be adopted from paper
prod_export_reaction.upper_bound = 1000 # but needs also to be adopted

product_s = Metabolite(
    'Product_Extracellular',
    formula = 'C211920H306299O79399N63505P4699S905',
    name = 'Product Extracellular',
    compartment = 's')

prod_export_reaction.add_metabolites({
    product_c : -1,
    product_s : 1
})

##### sink of product
prod_sink_reaction = Reaction('MAR_SINK')
prod_sink_reaction.name = 'Sink of Product'
prod_sink_reaction.subsystem = ''
prod_sink_reaction.lower_bound = 0 # but needs to be adopted from paper
prod_sink_reaction.upper_bound = 1000 # but needs also to be adopted

prod_sink_reaction.add_metabolites({
    product_s : -1
})

# #### export of empty aav from nucleus to cytosole

empty_capsid_transport_reaction = Reaction('MAR_TRANSPORT_EMPTY')
empty_capsid_transport_reaction.name = 'Transport of empty AAV to Cytosole'
empty_capsid_transport_reaction.subsystem = ''
empty_capsid_transport_reaction.lower_bound = 0 # but needs to be adopted from paper
empty_capsid_transport_reaction.upper_bound = 1000 # but needs also to be adopted

capsid_c = Metabolite(
    'Empty_AAV_Cytosole',
    formula = {sum_formula_capsid},
    name = 'Empty AAV Cytosole',
    compartment = 'c')

empty_capsid_transport_reaction.add_metabolites({
    capsid_n : -1,
    capsid_c : 1
})

# #### export of empty aav from cytosole to extracellular space

empty_export_reaction = Reaction('MAR_EXPORT_EMPTY')
empty_export_reaction.name = 'Export of empty AAV to Cytosole'
empty_export_reaction.subsystem = ''
empty_export_reaction.lower_bound = 0 # but needs to be adopted from paper
empty_export_reaction.upper_bound = 1000 # but needs also to be adopted

capsid_s = Metabolite(
    'Empty_AAV_Extracellular',
    formula = {sum_formula_capsid},
    name = 'Empty AAV Cytosole',
    compartment = 's')

empty_export_reaction.add_metabolites({
    capsid_c : -1,
    capsid_s : 1
})

##### sink of empty aav
empty_prod_sink_reaction = Reaction('MAR_SINK_EMPTY')
empty_prod_sink_reaction.name = 'Sink of empty AAV'
empty_prod_sink_reaction.subsystem = ''
empty_prod_sink_reaction.lower_bound = 0 # but needs to be adopted from paper
empty_prod_sink_reaction.upper_bound = 1000 # but needs also to be adopted

empty_prod_sink_reaction.add_metabolites({
    capsid_s : -1
})

print("Reactions done.")


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
    
    mod.add_reactions([vp1_reaction,
                   vp2_reaction,
                   vp3_reaction,
                   capsid_reaction,
                   nucleotide_reaction_1,
                   nucleotide_reaction_2,
                   target_gene_reaction,
                   assembly_reaction,
                   prod_transport_reaction,
                   prod_export_reaction,
                   prod_sink_reaction,
                   empty_export_reaction,
                   empty_capsid_transport_reaction,
                   empty_prod_sink_reaction ])
    
    mod.add_reactions([ala_reaction,
                  arg_reaction,
                  asn_reaction,
                  asp_reaction,
                  cys_reaction,
                  gln_reaction,
                  glu_reaction,
                  gly_reaction,
                  his_reaction,
                  ile_reaction,
                  leu_reaction,
                  lys_reaction,
                  met_reaction,
                  phe_reaction,
                  pro_reaction,
                  ser_reaction,
                  thr_reaction,
                  trp_reaction,
                  tyr_reaction,
                  val_reaction])
    
    orphan_metabolites = [m for m in mod.metabolites if len(m.reactions) == 0]
    mod.remove_metabolites(orphan_metabolites)
    
    # Remove orphan genes (not associated with any reactions)
    orphan_genes = [g for g in mod.genes if len(g.reactions) == 0]
    mod.genes = [g for g in mod.genes if g not in orphan_genes]
    
    return growth_rate, mod

hp_tr_4 = growth_reconstruction(P_T_04, df_hp_tr_4, reaction_matrix)
hp_tr_24 = growth_reconstruction(P_T_24, df_hp_tr_24, reaction_matrix)
hp_tr_48 = growth_reconstruction(P_T_48, df_hp_tr_48, reaction_matrix)
hp_tr_72 = growth_reconstruction(P_T_72, df_hp_tr_72, reaction_matrix)

# hp_mo_4 = growth_reconstruction(P_M_04, df_hp_mo_4, reaction_matrix)
# hp_mo_24 = growth_reconstruction(P_M_24, df_hp_mo_24, reaction_matrix)
# hp_mo_48 = growth_reconstruction(P_M_48, df_hp_mo_48, reaction_matrix)
# hp_mo_72 = growth_reconstruction(P_M_72, df_hp_mo_72, reaction_matrix)

# lp_tr_4 = growth_reconstruction(I_T_04, df_lp_tr_4, reaction_matrix)
# lp_tr_24 = growth_reconstruction(I_T_24, df_lp_tr_24, reaction_matrix)
# lp_tr_48 = growth_reconstruction(I_T_48, df_lp_tr_48, reaction_matrix)
# lp_tr_72 = growth_reconstruction(I_T_72, df_lp_tr_72, reaction_matrix)

# lp_mo_4 = growth_reconstruction(I_M_04, df_lp_mo_4, reaction_matrix)
# lp_mo_24 = growth_reconstruction(I_M_24, df_lp_mo_24, reaction_matrix)
# lp_mo_48 = growth_reconstruction(I_M_48, df_lp_mo_48, reaction_matrix)
# lp_mo_72 = growth_reconstruction(I_M_72, df_lp_mo_72, reaction_matrix)

hp_tr_4_model = hp_tr_4[1]
hp_tr_24_model = hp_tr_24[1]
hp_tr_48_model = hp_tr_48[1]
hp_tr_72_model = hp_tr_72[1]

# hp_mo_4_model = hp_mo_4[1]
# hp_mo_24_model = hp_mo_24[1]
# hp_mo_48_model = hp_mo_48[1]
# hp_mo_72_model = hp_mo_72[1]

# lp_tr_4_model = lp_tr_4[1]
# lp_tr_24_model = lp_tr_24[1]
# lp_tr_48_model = lp_tr_48[1]
# lp_tr_72_model = lp_tr_72[1]

# lp_mo_4_model = lp_mo_4[1]
# lp_mo_24_model = lp_mo_24[1]
# lp_mo_48_model = lp_mo_48[1]
# lp_mo_72_model = lp_mo_72[1]

print("Reconstruction done.")

print("Biomass adaption starts...")

models = [hp_tr_4_model, hp_tr_24_model, hp_tr_48_model, hp_tr_72_model]

bm_coeffs = pd.read_csv("/home/users/lzehetner/data/paper4_aav/HEK_AAV/bm_coefficients.csv", index_col = "Unnamed: 0")

# biomass metabolite ids
bm_metabolite_ids = ["MAM10014c", "MAM10013c", "MAM01721n", "MAM02847c", "MAM03161c"]

bm_coeffs.iloc[:, 0] = bm_metabolite_ids

# Set the first column as the index
bm_coeffs.set_index(bm_coeffs.columns[0], inplace=True)

# Split the DataFrame into groups of 4 columns and calculate the mean
means = [
    bm_coeffs.iloc[:, i:i+4].mean(axis=1)  # Calculate mean of every 4 columns
    for i in range(0, bm_coeffs.shape[1], 4)
]

# Combine the means into a new DataFrame
means_df = pd.concat(means, axis=1)

# Rename the columns for clarity (e.g., "Group 1", "Group 2", etc.)
means_df.columns = [f"Group {i+1}" for i in range(means_df.shape[1])]

# Reset the index to add the identifier column back as the first column
means_df.reset_index(inplace=True)

means_df = means_df.drop(columns=means_df.columns[[1, 2, 11, 12]])

model_names = ["hp_tr_4_model", "hp_tr_24_model", "hp_tr_48_model",  
               "hp_tr_72_model"]

columns_to_rename = means_df.columns[1:]  # Select columns from index 1 onwards
rename_mapping = dict(zip(columns_to_rename, model_names))
means_df.rename(columns=rename_mapping, inplace=True)

biomass_df = means_df

for i, model in enumerate(models):
    # Get the coefficients for the current model from the dataframe
    coefficients = biomass_df.iloc[:, i+1]  # Assuming column index 1 to 16 correspond to models
    metabolite_ids = biomass_df.iloc[:, 0]  # First column contains the metabolite IDs
    
    # Find the biomass reaction (adjust ID if needed)
    biomass_reaction = model.reactions.get_by_id("MAR13082")
    
    for metabolite_id, coefficient in zip(metabolite_ids, coefficients):
        metabolite = model.metabolites.get_by_id(metabolite_id)
        
        # Remove the current coefficient
        current_coefficient = biomass_reaction.metabolites.get(metabolite, 0)  # Use .get() to avoid KeyError
        biomass_reaction.add_metabolites({metabolite: -current_coefficient})
        
        # Add the new coefficient
        biomass_reaction.add_metabolites({metabolite: -coefficient})  # Add the new coefficient (negative as per your logic)
#        print(metabolite, "Removed:", current_coefficient, "Added:", -coefficient)

print("Biomass coefficients updated for all models.")

# ## constrain models with mus and qps

df = pd.read_csv("/home/users/lzehetner/data/paper4_aav/HEK_AAV/mus_qps.csv", index_col=0)

reactions = {
    "mu": "MAR13082",
    "qp_full": "MAR_SINK",
    "qp_empty": "MAR_SINK_EMPTY"
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

    model.reactions.MAR_SINK.upper_bound = 1000
    model.reactions.MAR_SINK.lower_bound = 0

    model.reactions.MAR_SINK_EMPTY.upper_bound = 1000
    model.reactions.MAR_SINK_EMPTY.lower_bound = 0

    model.reactions.MAR13082.upper_bound = 1000
    model.reactions.MAR13082.lower_bound = 0
    
    for index, reaction_id in reactions.items():
        if reaction_id in model.reactions:
            value = df.loc[index, column]
            reaction = model.reactions.get_by_id(reaction_id)
            reaction.lower_bound = value
            reaction.upper_bound = value
            print(f"Set bounds for reaction {reaction_id} in model '{column}': {value}")

    model.reactions.MAR13082.upper_bound = 1000
    model.reactions.MAR13082.lower_bound = 0
    
    mu_pred = model.optimize().objective_value
    print(mu_pred)

    model.reactions.MAR13082.upper_bound = mu_pred * 1.05
    model.reactions.MAR13082.lower_bound = mu_pred * 0.95

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

# HP TR 4
#selected_column = "hp_tr_4"
#cstr_hp_tr_4_model = mu_qp_constraining(hp_tr_4_model, selected_column)
#cstr_hp_tr_4_model.solver = "cplex"
#cstr_hp_tr_4_model.solver.configuration.tolerances.feasibility = 1e-8

#file_path = "/home/users/lzehetner/data/paper4_aav/HEK_AAV/flux_sampling/adapted_bm/flux_samples_hp_tr_4.csv"
#flux_samples = pd.read_csv(file_path)

#print("Flux sampling (HP-TR-4) started...")

#optgp = OptGPSampler(cstr_hp_tr_4_model, processes=10)

#flux_samples = pd.DataFrame()

#target_samples = 1000

#start_time = time.time()
#while len(flux_samples) < target_samples:
#    try:
        # Number of samples to collect in each iteration
#        remaining_samples = target_samples - len(flux_samples)
#        sample_count = min(remaining_samples, 10)  # Collect up to 10 samples per iteration
#        s2 = optgp.sample(sample_count)

        # Convert the samples to a DataFrame
#        sample_df = pd.DataFrame(s2, columns=cstr_hp_tr_4_model.reactions.list_attr("id"))

        # Append the valid samples to the main DataFrame
#        flux_samples = pd.concat([flux_samples, sample_df], ignore_index=True)

#    except RuntimeError as e:
        # Handle the numerically unstable error without breaking the loop
#        print(f"Encountered an issue: {e}. Skipping this batch of samples.")
#    except Exception as e:
        # Handle any other unexpected errors
#        print(f"Unexpected error: {e}.")
#        break  # Exit the loop on critical errors
#end_time = time.time()
#elapsed_time = end_time - start_time
#print(elapsed_time)

#print("Flux sampling done.")
#print("Loop correction starts ...")

#cstr_hp_tr_4_model.reactions.MAR13082.upper_bound = 1000
#cstr_hp_tr_4_model.reactions.MAR13082.lower_bound = 0


#start_time = time.time()
#results = calculate_loopless_fluxes(cstr_hp_tr_4_model, flux_samples)
#end_time = time.time()
#elapsed_time = end_time - start_time
#print(elapsed_time)
#flux_distributions_df = results[0]
#solver_results_df = results[1]

#print("Loop correction done.")

#flux_samples.to_csv("/home/users/lzehetner/data/paper4_aav/HEK_AAV/flux_sampling/adapted_bm/flux_samples_hp_tr_4.csv", index=False)
#flux_distributions_df.to_csv("/home/users/lzehetner/data/paper4_aav/HEK_AAV/flux_sampling/adapted_bm/loopless_flux_samples_hp_tr_4.csv", index = False)
#solver_results_df.to_csv("/home/users/lzehetner/data/paper4_aav/HEK_AAV/flux_sampling/adapted_bm/solver_status_hp_tr_4.csv", index = False)


# HP TR 24
#selected_column = "hp_tr_24"
#cstr_hp_tr_24_model = mu_qp_constraining(hp_tr_24_model, selected_column)
#cstr_hp_tr_24_model.solver = "cplex"
#cstr_hp_tr_24_model.solver.configuration.tolerances.feasibility = 1e-8

#file_path = "/home/users/lzehetner/data/paper4_aav/HEK_AAV/flux_sampling/adapted_bm/flux_samples_hp_tr_24.csv"
#flux_samples = pd.read_csv(file_path)

#print("Flux sampling (HP-TR-24) started...")

#optgp = OptGPSampler(cstr_hp_tr_24_model, processes=10)

#flux_samples = pd.DataFrame()

#target_samples = 1000

#start_time = time.time()
#while len(flux_samples) < target_samples:
#    try:
        # Number of samples to collect in each iteration
#        remaining_samples = target_samples - len(flux_samples)
#        sample_count = min(remaining_samples, 10)  # Collect up to 10 samples per iteration
#        s2 = optgp.sample(sample_count)

        # Convert the samples to a DataFrame
#        sample_df = pd.DataFrame(s2, columns=cstr_hp_tr_24_model.reactions.list_attr("id"))

        # Append the valid samples to the main DataFrame
#        flux_samples = pd.concat([flux_samples, sample_df], ignore_index=True)

#    except RuntimeError as e:
        # Handle the numerically unstable error without breaking the loop
#        print(f"Encountered an issue: {e}. Skipping this batch of samples.")
#    except Exception as e:
        # Handle any other unexpected errors
#        print(f"Unexpected error: {e}.")
#        break  # Exit the loop on critical errors
#end_time = time.time()
#elapsed_time = end_time - start_time
#print(elapsed_time)

#print("Flux sampling done.")
#print("Loop correction starts ...")

#cstr_hp_tr_24_model.reactions.MAR13082.upper_bound = 1000
#cstr_hp_tr_24_model.reactions.MAR13082.lower_bound = 0

#start_time = time.time()
#results = calculate_loopless_fluxes(cstr_hp_tr_24_model, flux_samples)
#end_time = time.time()
#elapsed_time = end_time - start_time
#print(elapsed_time)
#flux_distributions_df = results[0]
#solver_results_df = results[1]

#print("Loop correction done.")

#flux_samples.to_csv("/home/users/lzehetner/data/paper4_aav/HEK_AAV/flux_sampling/adapted_bm/flux_samples_hp_tr_24.csv", index=False)
#flux_distributions_df.to_csv("/home/users/lzehetner/data/paper4_aav/HEK_AAV/flux_sampling/adapted_bm/loopless_flux_samples_hp_tr_24.csv", index = False)
#solver_results_df.to_csv("/home/users/lzehetner/data/paper4_aav/HEK_AAV/flux_sampling/adapted_bm/solver_status_hp_tr_24.csv", index = False)

# HP TR 48
#selected_column = "hp_tr_48"
#cstr_hp_tr_48_model = mu_qp_constraining(hp_tr_48_model, selected_column)
#cstr_hp_tr_48_model.solver = "cplex"
#cstr_hp_tr_48_model.solver.configuration.tolerances.feasibility = 1e-8

#file_path = "/home/users/lzehetner/data/paper4_aav/HEK_AAV/flux_sampling/adapted_bm/flux_samples_hp_tr_48.csv"
#flux_samples = pd.read_csv(file_path)

#print("Flux sampling (HP-TR-48) started...")

#optgp = OptGPSampler(cstr_hp_tr_48_model, processes=10)

#flux_samples = pd.DataFrame()

#target_samples = 1000

#start_time = time.time()
#while len(flux_samples) < target_samples:
#    try:
        # Number of samples to collect in each iteration
#        remaining_samples = target_samples - len(flux_samples)
#        sample_count = min(remaining_samples, 10)  # Collect up to 10 samples per iteration
#        s2 = optgp.sample(sample_count)

        # Convert the samples to a DataFrame
#        sample_df = pd.DataFrame(s2, columns=cstr_hp_tr_48_model.reactions.list_attr("id"))

        # Append the valid samples to the main DataFrame
#        flux_samples = pd.concat([flux_samples, sample_df], ignore_index=True)

#    except RuntimeError as e:
        # Handle the numerically unstable error without breaking the loop
#        print(f"Encountered an issue: {e}. Skipping this batch of samples.")
#    except Exception as e:
        # Handle any other unexpected errors
#        print(f"Unexpected error: {e}.")
#        break  # Exit the loop on critical errors
#end_time = time.time()
#elapsed_time = end_time - start_time
#print(elapsed_time)

#print("Flux sampling done.")
#print("Loop correction starts ...")

#cstr_hp_tr_48_model.reactions.MAR13082.upper_bound = 1000
#cstr_hp_tr_48_model.reactions.MAR13082.lower_bound = 0

#start_time = time.time()
#results = calculate_loopless_fluxes(cstr_hp_tr_48_model, flux_samples)
#end_time = time.time()
#elapsed_time = end_time - start_time
#print(elapsed_time)
#flux_distributions_df = results[0]
#solver_results_df = results[1]

#print("Loop correction done.")

#flux_samples.to_csv("/home/users/lzehetner/data/paper4_aav/HEK_AAV/flux_sampling/adapted_bm/flux_samples_hp_tr_48.csv", index=False)
#flux_distributions_df.to_csv("/home/users/lzehetner/data/paper4_aav/HEK_AAV/flux_sampling/adapted_bm/loopless_flux_samples_hp_tr_48.csv", index = False)
#solver_results_df.to_csv("/home/users/lzehetner/data/paper4_aav/HEK_AAV/flux_sampling/adapted_bm/solver_status_hp_tr_48.csv", index = False)

# HP TR 72
selected_column = "hp_tr_72"
cstr_hp_tr_72_model = mu_qp_constraining(hp_tr_72_model, selected_column)
cstr_hp_tr_72_model.solver = "cplex"
cstr_hp_tr_72_model.solver.configuration.tolerances.feasibility = 1e-8

file_path = "/home/users/lzehetner/data/paper4_aav/HEK_AAV/flux_sampling/adapted_bm/flux_samples_hp_tr_72.csv"
flux_samples = pd.read_csv(file_path)

#print("Flux sampling (HP-TR-72) started...")

#optgp = OptGPSampler(cstr_hp_tr_72_model, processes=10)

#flux_samples = pd.DataFrame()

#target_samples = 1000

#start_time = time.time()
#while len(flux_samples) < target_samples:
#    try:
        # Number of samples to collect in each iteration
#        remaining_samples = target_samples - len(flux_samples)
#        sample_count = min(remaining_samples, 10)  # Collect up to 10 samples per iteration
#        s2 = optgp.sample(sample_count)

        # Convert the samples to a DataFrame
#        sample_df = pd.DataFrame(s2, columns=cstr_hp_tr_72_model.reactions.list_attr("id"))

        # Append the valid samples to the main DataFrame
#        flux_samples = pd.concat([flux_samples, sample_df], ignore_index=True)

#    except RuntimeError as e:
        # Handle the numerically unstable error without breaking the loop
#        print(f"Encountered an issue: {e}. Skipping this batch of samples.")
#    except Exception as e:
        # Handle any other unexpected errors
#        print(f"Unexpected error: {e}.")
#        break  # Exit the loop on critical errors
#end_time = time.time()
#elapsed_time = end_time - start_time
#print(elapsed_time)

#print("Flux sampling done.")

print("Loop correction starts ...")

#cstr_hp_tr_72_model.reactions.MAR13082.upper_bound = 1000
#cstr_hp_tr_72_model.reactions.MAR13082.lower_bound = 0

start_time = time.time()
results = calculate_loopless_fluxes(cstr_hp_tr_72_model, flux_samples)
end_time = time.time()
elapsed_time = end_time - start_time
print(elapsed_time)
flux_distributions_df = results[0]
solver_results_df = results[1]

print("Loop correction done.")

#flux_samples.to_csv("/home/users/lzehetner/data/paper4_aav/HEK_AAV/flux_sampling/adapted_bm/flux_samples_hp_tr_72.csv", index=False)
flux_distributions_df.to_csv("/home/users/lzehetner/data/paper4_aav/HEK_AAV/flux_sampling/adapted_bm/loopless_flux_samples_hp_tr_72.csv", index = False)
solver_results_df.to_csv("/home/users/lzehetner/data/paper4_aav/HEK_AAV/flux_sampling/adapted_bm/solver_status_hp_tr_72.csv", index = False)



print("Done with everything")











