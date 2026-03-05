# %% Imports
import os
from pathlib import Path
from datetime import datetime
import tempfile
import os
import subprocess
import sys

import pandas as pd
import cobra
import memote

# Import the utility function
from python.utils.model_utils import load_models, get_models_dir

# %% Set-up models dirs
MODELS_DIR = get_models_dir()

# %% Load all models using utility function
models, summary_df = load_models()

print(summary_df[["model_name", "objective_rxn"]])

# %% Check biomass metabolites for each model
for model_name, model in models.items():
    # Get the objective reaction expression as string
    obj_expr = str(model.objective.expression)

    # Heuristic: get first reaction ID from expression
    if '*' in obj_expr:
        # Example: '1.0*biomass_cho - 1.0*biomass_cho_reverse'
        rxn_id = obj_expr.split('*')[1].split()[0]
    else:
        rxn_id = obj_expr  # fallback if expression is just a single ID

    # Safely get reaction and count metabolites
    try:
        biomass_rxn = model.reactions.get_by_id(rxn_id)
        n_mets = len(biomass_rxn.metabolites)
        print(f"{model_name}: Biomass '{rxn_id}' has {n_mets} metabolites")
    except KeyError:
        print(f"{model_name}: Reaction '{rxn_id}' not found in model")

# %% Extract exchange reactions for each model
# Model-specific exchange reactions for experimental metabolites
model_uptake_names = {
    "iCHOv1": {
        "Ala": "EX_ala__L_e",
        "Arg": "EX_arg__L_e",
        "Asn": "EX_asn__L_e",
        "Asp": "EX_asp__L_e",
        #"Cys": "EX_cys__L_e",
        "Glu": "EX_glu__L_e",
        "Gln": "EX_gln__L_e",
        "Gly": "EX_gly_e",
        "His": "EX_his__L_e",
        "Ile": "EX_ile__L_e",
        "Leu": "EX_leu__L_e",
        "Lys": "EX_lys__L_e",
        "Met": "EX_met__L_e",
        "Phe": "EX_phe__L_e",
        "Pro": "EX_pro__L_e",
        "Ser": "EX_ser__L_e",
        "Thr": "EX_thr__L_e",
        "Trp": "EX_trp__L_e",
        "Tyr": "EX_tyr__L_e",
        "Val": "EX_val__L_e",
        "Titer": "DM_igg_g",
        "NH3": "EX_nh4_e",
        "GLC": "EX_glc__D_e",
        "LAC": "EX_lac__L_e",
        "Growth_rate": "BIOMASS_cho_producing"
    },
    "iCHO3K_cho_prod_generic_unblocked": {
        "Ala": "EX_ala_L_e",
        "NH3": "EX_nh4_e",
        "Arg": "EX_arg_L_e",
        "Asn": "EX_asn_L_e",
        "Asp": "EX_asp_L_e",
        "GLC": "EX_glc_e",
        "Glu": "EX_glu_L_e",
        "Gln": "EX_gln_L_e",
        "Gly": "EX_gly_e",
        "His": "EX_his_L_e",
        "Ile": "EX_ile_L_e",
        "LAC": "EX_lac_L_e",
        "Leu": "EX_leu_L_e",
        "Lys": "EX_lys_L_e",
        "Met": "EX_met_L_e",
        "Phe": "EX_phe_L_e",
        "Pro": "EX_pro_L_e",
        "Ser": "EX_ser_L_e",
        "Thr": "EX_thr_L_e",
        "Trp": "EX_trp_L_e",
        "Tyr": "EX_tyr_L_e",
        "Val": "EX_val_L_e",
        "Titer": "DM_igg_g",
        "Growth_rate" : "biomass_cho_prod"
    },
    "iCHO2441": {
        "Ala": "EX_ala_L(e)",
        "Arg": "EX_arg_L(e)",
        "Asn": "EX_asn_L(e)",
        "Asp": "EX_asp_L(e)",
        #"Cys": "EX_cys_L(e)",
        "Glu": "EX_glu_L(e)",
        "Gln": "EX_gln_L(e)",
        "Gly": "EX_gly(e)",
        "His": "EX_his_L(e)",
        "Ile": "EX_ile_L(e)",
        "Leu": "EX_leu_L(e)",
        "Lys": "EX_lys_L(e)",
        "Met": "EX_met_L(e)",
        "Phe": "EX_phe_L(e)",
        "Pro": "EX_pro_L(e)",
        "Ser": "EX_ser_L(e)",
        "Thr": "EX_thr_L(e)",
        "Trp": "EX_trp_L(e)",
        "Tyr": "EX_tyr_L(e)",
        "Val": "EX_val_L(e)",
        "Titer": "DM_igg[g]",
        "NH3": "EX_nh4(e)",
        "GLC": "EX_glc(e)",
        "LAC": "EX_lac_L(e)",
        "Growth_rate": "biomass_cho_prod"
    }
}

filtered_rows = []

for model_name, model in models.items():
    print(f"\nChecking target exchanges for {model_name}...")

    uptake_names = model_uptake_names.get(model_name, {})

    for label, ex_id in uptake_names.items():
        # Try exact match first
        try:
            rxn = model.reactions.get_by_id(ex_id)
        except KeyError:
            rxn = None

        # Fallback: substring match
        if rxn is None:
            for r in model.reactions:
                if ex_id.strip("_").lower() in r.id.lower():
                    rxn = r
                    break

        # Append row
        if rxn is not None:
            met = list(rxn.metabolites.keys())[0]
            filtered_rows.append({
                "model": model_name,
                "label": label,
                "metabolite_id": met.id,
                "metabolite_name": met.name,
                "reaction_id": rxn.id,
                "lower_bound": rxn.lower_bound,
                "upper_bound": rxn.upper_bound
            })
        else:
            filtered_rows.append({
                "model": model_name,
                "label": label,
                "metabolite_id": None,
                "metabolite_name": None,
                "reaction_id": None,
                "lower_bound": None,
                "upper_bound": None
            })

exchange_igg_growth_rxn_df = pd.DataFrame(filtered_rows)
print(exchange_igg_growth_rxn_df )

# Filter out the Growth_rate row
exchange_igg_rxn_df = exchange_igg_growth_rxn_df [exchange_igg_growth_rxn_df ["label"] != "Growth_rate"]

# Save as CSV
csv_path = MODELS_DIR / "exchange_igg_rxn_id.csv"
exchange_igg_rxn_df.to_csv(csv_path, index=False)
print(f"Saved exchange reactions (with IGG) to CSV: {csv_path}")

# Save as Pickle (optional, faster reload in Python)
pkl_path = MODELS_DIR / "exchange_igg_rxn_id.pkl"
exchange_igg_rxn_df.to_pickle(pkl_path)
print(f"Saved exchange reactions (with IGG) to Pickle: {pkl_path}")


# Filter out the Titer row
exchange_growth_rxn_df = exchange_igg_growth_rxn_df [exchange_igg_growth_rxn_df ["label"] != "Titer"]

# Save as CSV
csv_path = MODELS_DIR / "exchange_growth_rxn_id.csv"
exchange_growth_rxn_df.to_csv(csv_path, index=False)
print(f"Saved exchange reactions (with Growth) to CSV: {csv_path}")

# Save as Pickle (optional, faster reload in Python)
pkl_path = MODELS_DIR / "exchange_growth_rxn_id.pkl"
exchange_growth_rxn_df.to_pickle(pkl_path)
print(f"Saved exchange reactions (with Growth) to Pickle: {pkl_path}")
