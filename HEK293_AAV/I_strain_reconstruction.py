#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import cobra
from cobra.io import read_sbml_model
from corda import CORDA
from corda import reaction_confidence
from cobra import Model, Reaction, Metabolite


output_rxns = "/home/users/lzehetner/data/paper4_aav/reaction_matrix_LP.csv"
input_transcriptome = "/home/users/lzehetner/data/hek/transcriptome/tpm_I.csv"
input_gpr = "/home/users/lzehetner/gpr_human1.csv"
input_model = "/home/users/lzehetner/data/human1/human1.xml"
mod = cobra.io.read_sbml_model(input_model)

reactions_ids = [reaction.id for reaction in mod.reactions]
reactions_df = pd.DataFrame(reactions_ids, columns=['reaction_id'])

mod.objective = "MAR13082"
print(mod.optimize().objective_value)

# Extract data from transcriptome
tr = pd.read_csv(input_transcriptome)

tr = tr.drop(tr.columns[0], axis=1)

# function for reaction extraction:

gpr_ori = pd.read_csv(input_gpr, sep = ";")

def rxn_extraction(gene_list):
    gpr = gpr_ori.iloc[:, [0,2]]
    rxn_wo_genes = gpr[gpr.iloc[:, 1].isna()]
    gpr.dropna(inplace=True)
    gpr_and = gpr[gpr['Gene-reaction association'].str.contains('and')]
    gpr = gpr[~gpr.isin(gpr_and)].dropna()
    gpr_or = gpr[gpr['Gene-reaction association'].str.contains('or')]
    gpr_one = gpr[~gpr.isin(gpr_or)].dropna()
    
    tr_list = gene_list.values.tolist()
    
    dict_gpr_one = gpr_one.set_index("Rxn name")["Gene-reaction association"].to_dict()
    rxns_w_one = [key for (key,value) in dict_gpr_one.items() if value in tr_list]

    gpr_and_1 = gpr_and[gpr_and['Gene-reaction association'].str.contains("\(")]

    gpr_and_2 = gpr_and_1[gpr_and_1['Gene-reaction association'].str.contains("\) or")]
    x = gpr_and_2["Gene-reaction association"].str.split(' or ').tolist()
    y = gpr_and_2["Rxn name"].tolist()
    gpr_and_dict_2 = {y[i]: x[i] for i in range(len(y))}
    df = pd.DataFrame(list(gpr_and_dict_2.items()))
    df = df.rename(columns={1: 'Genes'})
    df = df.explode('Genes')
    df_and = df[df["Genes"].str.contains("\(")]
    df_and['Genes'] = df_and['Genes'].str.replace("\( ", "")
    df_and['Genes'] = df_and['Genes'].str.replace(" \)", "")
    df_or = df[df["Genes"].str.contains(r'^E')]

    x = gpr_or["Gene-reaction association"].str.split(' or ').tolist()
    y = gpr_or["Rxn name"].tolist()
    gpr_or_dict = {y[i]: x[i] for i in range(len(y))}
    df = pd.DataFrame(list(gpr_or_dict.items()))
    df = df.rename(columns={1: 'Genes'})
    df1 = df.explode('Genes')
    df = pd.concat([df1, df_or])
    rxns_w_iso = df[df["Genes"].isin(tr_list)]
    rxns_w_iso_unique = rxns_w_iso[0].unique().tolist()

    x = gpr_and[~gpr_and.isin(gpr_and_1)].dropna()
    x1 = x["Gene-reaction association"].str.split(' and ').tolist()
    y = x["Rxn name"].tolist()
    gpr_and_dict = {y[i]: x1[i] for i in range(len(y))}
    df = pd.DataFrame(list(gpr_and_dict.items()))
    df1 = df.rename(columns={1: 'Genes'})
    df = pd.concat([df1, df_and])
    df['new'] = df['Genes'].explode().isin(tr_list).groupby(level=0).any()
    z = df[df[["new"]].all(axis=1)]
    z = z.iloc[:, [0,1]]
    rxns_w_compl = z.iloc[:, 0].values.tolist()

    nested_rxns = ["MAR04137", "MAR07161", "MAR07162"]
    rxns = list(set(rxns_w_one + rxns_w_iso_unique + rxns_w_compl + nested_rxns))
    
    return rxns

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

secr_comp = ["MAR09422", # Thymine
             "MAR09142", # Nicotinate
             "MAR09437", # Uracil
             "MAR09378", # Nicotinamide
             "MAR11428", # Xanthine
             "MAR09358", # Hypoxanthine
             "MAR09021", # 5-Methyladenosine
             "MAR09253", # Adenine
             "MAR09341", # Betaine
             "MAR09353", # Guanine
             "MAR09848", # Dimethylglycine
             "MAR09290", # Creatine
             "MAR09131", # Sarcosine
             "MAR10428", # 5-Aminoevulinic acid
             "MAR04849", # Citrulline
             "MAR09087", # Ornithine
             "MAR11400", # Fumarate
             "MAR09415", # Succinate
             "MAR09025", # Pyroglutamate / 5-oxoproline
             "MAR11404", # Malate
             "MAR09286", # Citrate
             "MAR09861", # Xanthosine
             "MAR10431", # Allantoine
             "MAR09285" # Cholesterole
            ]

# select strain dependent transcriptome, mean them and remove all genes below 0.2 (as usual for TPM) to get expressed genes

for column in tr.columns:
    if column != 'gene':
        
        genes = pd.DataFrame(tr[["gene", column]])
        genes['mean'] = genes.mean(axis=1)
        genes = genes.loc[genes["mean"] >= 0.2]

        q = genes.quantile([0.00, 0.10, 0.25, 0.50, 0.75, 0.90, 1.00])
        col = 'mean'
        a = genes[((genes[col]>=q[col][0.10]) & (genes[col]<=q[col][1.00]))]
        b = a.iloc[:, 0]
        ntile = b.to_frame()

        # ## Rxn Extraction

        rxns = pd.DataFrame(rxn_extraction(ntile["gene"]), columns = ['Rxns'])

        print("Rxns are extracted")
        print(len(rxns))

        for r in mod.exchanges:
            mod.reactions.get_by_id(r.id).lower_bound = 0.0

        for r in min_media_comp:
            mod.reactions.get_by_id(r).lower_bound = -1000.0

        mod.solver = 'cplex'

        conf = {}
        for r in mod.reactions: conf[r.id] = -1
        for r in rxns["Rxns"]: conf[r] = 3    
        a = mod.exchanges
        b = []
        for r in a:
            b.append(r.id)
        for r in b: conf[r] = -1

        for r in min_media_comp: conf[r] = 3
        for r in secr_comp: conf[r] = 3
        conf["MAR13082"] = 3 # Biomass
        conf["MAR10023"] = 3 # Biomass export
        conf["MAR10024"] = 3 # Biomass export
        conf["MAR07160"] = 3 # DNA Pool reaction
        conf["MAR13086"] = 3 # RNA Pool reaction
        conf["MAR10065"] = 3 # Cofactor Pool reaction
        conf["MAR09727"] = 3 # Glycogen Pool reaction
        conf["MAR10063"] = 3 # Lipid Pool reaction
        conf["MAR10064"] = 3 # Metabolite Pool reaction
        conf["MAR10062"] = 3 # Protein Pool reaction
        conf.pop("Rxns", None)

        print(f"Model reconstruction starts for {column}")

        model = CORDA(mod, conf)
        model.build()
        model_reconstr = model.cobra_model()

        model_reconstr.objective = "MAR13082"
        print(model_reconstr.optimize().objective_value)

        reactions = []
        for r in model_reconstr.reactions:
            reactions.append(r.id)

        # >25% reconstruction

        mod = cobra.io.read_sbml_model(input_model)    

        a = genes[((genes[col]>=q[col][0.25]) & (genes[col]<=q[col][1.00]))]
        b = a.iloc[:, 0]
        ntile = b.to_frame()

        rxns = pd.DataFrame(rxn_extraction(ntile["gene"]), columns = ['Rxns'])

        for r in mod.exchanges:
            mod.reactions.get_by_id(r.id).lower_bound = 0.0

        for r in min_media_comp:
            mod.reactions.get_by_id(r).lower_bound = -1000.0

        mod.solver = 'cplex'

        conf = {}
        for r in mod.reactions: conf[r.id] = -1
        for r in rxns["Rxns"]: conf[r] = 3    
        a = mod.exchanges
        b = []
        for r in a:
            b.append(r.id)
        for r in b: conf[r] = -1

        for r in min_media_comp: conf[r] = 3
        for r in secr_comp: conf[r] = 3
        conf["MAR13082"] = 3 # Biomass
        conf["MAR10023"] = 3 # Biomass export
        conf["MAR10024"] = 3 # Biomass export
        conf["MAR07160"] = 3 # DNA Pool reaction
        conf["MAR13086"] = 3 # RNA Pool reaction
        conf["MAR10065"] = 3 # Cofactor Pool reaction
        conf["MAR09727"] = 3 # Glycogen Pool reaction
        conf["MAR10063"] = 3 # Lipid Pool reaction
        conf["MAR10064"] = 3 # Metabolite Pool reaction
        conf["MAR10062"] = 3 # Protein Pool reaction
        conf.pop("Rxns", None)

        print(f"Model reconstruction starts for {column}")

        model = CORDA(mod, conf)
        model.build()
        model_reconstr = model.cobra_model()

        for r in model_reconstr.reactions:
            reactions.append(r.id)

        # # >50% reconstruction

        mod = cobra.io.read_sbml_model(input_model)    

        a = genes[((genes[col]>=q[col][0.50]) & (genes[col]<q[col][1.00]))]
        b = a.iloc[:, 0]
        ntile = b.to_frame()

        rxns = pd.DataFrame(rxn_extraction(ntile["gene"]), columns = ['Rxns'])

        for r in mod.exchanges:
            mod.reactions.get_by_id(r.id).lower_bound = 0.0

        for r in min_media_comp:
            mod.reactions.get_by_id(r).lower_bound = -1000.0

        mod.solver = 'cplex'

        conf = {}
        for r in mod.reactions: conf[r.id] = -1
        for r in rxns["Rxns"]: conf[r] = 3    
        a = mod.exchanges
        b = []
        for r in a:
            b.append(r.id)
        for r in b: conf[r] = -1

        for r in min_media_comp: conf[r] = 3
        for r in secr_comp: conf[r] = 3
        conf["MAR13082"] = 3 # Biomass
        conf["MAR10023"] = 3 # Biomass export
        conf["MAR10024"] = 3 # Biomass export
        conf["MAR07160"] = 3 # DNA Pool reaction
        conf["MAR13086"] = 3 # RNA Pool reaction
        conf["MAR10065"] = 3 # Cofactor Pool reaction
        conf["MAR09727"] = 3 # Glycogen Pool reaction
        conf["MAR10063"] = 3 # Lipid Pool reaction
        conf["MAR10064"] = 3 # Metabolite Pool reaction
        conf["MAR10062"] = 3 # Protein Pool reaction
        conf.pop("Rxns", None)

        print(f"Model reconstruction starts for {column}")

        model = CORDA(mod, conf)
        model.build()
        model_reconstr = model.cobra_model()

        for r in model_reconstr.reactions:
            reactions.append(r.id)

        # # >75% reconstruction

        mod = cobra.io.read_sbml_model(input_model)    

        a = genes[((genes[col]>=q[col][0.75]) & (genes[col]<q[col][1.00]))]
        b = a.iloc[:, 0]
        ntile = b.to_frame()

        rxns = pd.DataFrame(rxn_extraction(ntile["gene"]), columns = ['Rxns'])

        for r in mod.exchanges:
            mod.reactions.get_by_id(r.id).lower_bound = 0.0

        for r in min_media_comp:
            mod.reactions.get_by_id(r).lower_bound = -1000.0

        mod.solver = 'cplex'

        conf = {}
        for r in mod.reactions: conf[r.id] = -1
        for r in rxns["Rxns"]: conf[r] = 3    
        a = mod.exchanges
        b = []
        for r in a:
            b.append(r.id)
        for r in b: conf[r] = -1

        for r in min_media_comp: conf[r] = 3
        for r in secr_comp: conf[r] = 3
        conf["MAR13082"] = 3 # Biomass
        conf["MAR10023"] = 3 # Biomass export
        conf["MAR10024"] = 3 # Biomass export
        conf["MAR07160"] = 3 # DNA Pool reaction
        conf["MAR13086"] = 3 # RNA Pool reaction
        conf["MAR10065"] = 3 # Cofactor Pool reaction
        conf["MAR09727"] = 3 # Glycogen Pool reaction
        conf["MAR10063"] = 3 # Lipid Pool reaction
        conf["MAR10064"] = 3 # Metabolite Pool reaction
        conf["MAR10062"] = 3 # Protein Pool reaction
        conf.pop("Rxns", None)

        print(f"Model reconstruction starts for {column}")

        model = CORDA(mod, conf)
        model.build()
        model_reconstr = model.cobra_model()

        for r in model_reconstr.reactions:
            reactions.append(r.id)

        # # >90% reconstruction

        mod = cobra.io.read_sbml_model(input_model)    

        a = genes[((genes[col]>=q[col][0.90]) & (genes[col]<q[col][1.00]))]
        b = a.iloc[:, 0]
        ntile = b.to_frame()

        rxns = pd.DataFrame(rxn_extraction(ntile["gene"]), columns = ['Rxns'])

        for r in mod.exchanges:
            mod.reactions.get_by_id(r.id).lower_bound = 0.0

        for r in min_media_comp:
            mod.reactions.get_by_id(r).lower_bound = -1000.0

        mod.solver = 'cplex'

        conf = {}
        for r in mod.reactions: conf[r.id] = -1
        for r in rxns["Rxns"]: conf[r] = 3    
        a = mod.exchanges
        b = []
        for r in a:
            b.append(r.id)
        for r in b: conf[r] = -1

        for r in min_media_comp: conf[r] = 3
        for r in secr_comp: conf[r] = 3
        conf["MAR13082"] = 3 # Biomass
        conf["MAR10023"] = 3 # Biomass export
        conf["MAR10024"] = 3 # Biomass export
        conf["MAR07160"] = 3 # DNA Pool reaction
        conf["MAR13086"] = 3 # RNA Pool reaction
        conf["MAR10065"] = 3 # Cofactor Pool reaction
        conf["MAR09727"] = 3 # Glycogen Pool reaction
        conf["MAR10063"] = 3 # Lipid Pool reaction
        conf["MAR10064"] = 3 # Metabolite Pool reaction
        conf["MAR10062"] = 3 # Protein Pool reaction
        conf.pop("Rxns", None)

        print(f"Model reconstruction starts for {column}")

        model = CORDA(mod, conf)
        model.build()
        model_reconstr = model.cobra_model()

        for r in model_reconstr.reactions:
            reactions.append(r.id)

        # print("10% reconstruction led to:", len(reactions))
        
        reactions = list(set(reactions))

        reactions_df[column] = reactions_df['reaction_id'].apply(lambda x: 1 if x in reactions else 0)

#reactions_df = pd.DataFrame(reactions, columns = ['Rxns']).drop_duplicates()
# the export file needs also to be changed every time
reactions_df.to_csv(output_rxns, index=False, header=True)

