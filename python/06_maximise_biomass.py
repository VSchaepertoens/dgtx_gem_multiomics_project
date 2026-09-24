# %% Section 0: Imports
import pandas as pd
import cobra
from cobra.flux_analysis import pfba
import os
import re
from python.utils.model_utils import load_models, get_models_dir, get_common_model_name, parse_condition_window
from pathlib import Path

cobra.Configuration().solver = "glpk"

# %% Section 1: Load models
MODELS_DIR = get_models_dir(mode="exp_constrained")
models, summary_df = load_models(MODELS_DIR)

#%% # Project root is two levels up from python/
PROJECT_ROOT = Path.cwd()
SAVE_FOLDER = PROJECT_ROOT / "analysis" / "models" / "exp_constrained_ngam"
SAVE_FOLDER.mkdir(parents=True, exist_ok=True)
print(f"Analysis will be saved to: {SAVE_FOLDER}")

# %% Section 2.1: Load experimental data
rates = pd.read_csv("data/rates/rates_mM_gDCW_h_new.csv")
# Clean column names if needed
rates.columns = [c.strip() for c in rates.columns]

# Load exchange reaction mapping
ex_rxn_df = pd.read_csv("models/original/exchange_igg_rxn_id.csv")
ex_rxn_df.columns = [c.strip() for c in ex_rxn_df.columns]

# Build a lookup dictionary: model -> metabolite -> reaction ID
uptake_names = {}
for _, row in ex_rxn_df.iterrows():
    uptake_names.setdefault(row['model'], {})[row['label']] = row['reaction_id']

biomass_reactions = {
    "iCHO3K_cho_prod_generic_unblocked": "biomass_cho_prod",
    "iCHO2441": "biomass_cho_prod",
    "iCHOv1": "BIOMASS_cho_producing"
}

igg_reactions = {
    "iCHO3K_cho_prod_generic_unblocked": "DM_igg_g",
    "iCHO2441": "DM_igg[g]",
    "iCHOv1": "DM_igg_g"
}

ngam_reactions = {
    "iCHO3K_cho_prod_generic_unblocked": "ATPM",
    "iCHO2441": "DM_atp[c]",
    "iCHOv1": "ATPM"
}

# %% Define bounds for energy maintenance reaction
NGAM_MEAN = 4.3
NGAM_SD   = 1.7
# %%
# def parse_condition_window(model_name: str):
#     """
#     Extract Condition and Window from model name.
#
#     Example:
#     iCHOv1_Constant_4_constrained -> ("Constant", 4)
#     iCHOv1_Temp. shifted_2_constrained -> ("Temp. shifted", 2)
#     """
#     name = model_name.replace("_constrained", "")
#
#     if "_Constant_" in name:
#         parts = name.split("_Constant_")
#         return "Constant", int(parts[1])
#
#     elif "_Temp. shifted_" in name:
#         parts = name.split("_Temp. shifted_")
#         return "Temp. shifted", int(parts[1])
#
#     else:
#         print(f"⚠️ Could not parse condition/window from {model_name}")
#         return None, None


# %% Section 4.1: Loop through Models, Conditions, Windows. Set bounds, save models and get FBA/pFBA solutions
# maximize biomass production under observed titer + nutrient uptake

mus_records = []
reaction_records = []
models_solved = {}

# Initialize dictionaries for all models
mus_fba = {}
mus_pfba = {}
reaction_data_fba = {}
reaction_data_pfba = {}

for model_name, base_model in models.items():
    print(f"\nProcessing model {model_name}")

    # 🔹 Extract common model name
    common_name = get_common_model_name(model_name)

    # 🔹 Extract condition + window
    ex, w = parse_condition_window(model_name)

    if ex is None:
        continue

    print(f"  → Using Condition {ex}, Window {w}")

    # 🔹 Safety checks
    if common_name not in uptake_names:
        print(f"⚠️ No uptake mapping for {model_name}")
        continue

    if common_name not in biomass_reactions:
        print(f"⚠️ No biomass reaction for {model_name}")
        continue

    # Initialize dicts
    mus_fba[model_name] = {}
    mus_pfba[model_name] = {}
    reaction_data_fba[model_name] = {}
    reaction_data_pfba[model_name] = {}

    mus_fba[model_name][ex] = {}
    mus_pfba[model_name][ex] = {}
    reaction_data_fba[model_name][ex] = {}
    reaction_data_pfba[model_name][ex] = {}

    mus_fba[model_name][ex][w] = None
    mus_pfba[model_name][ex][w] = None
    reaction_data_fba[model_name][ex][w] = {}
    reaction_data_pfba[model_name][ex][w] = {}

    with base_model.copy() as model:

        # Reset bounds
        for reaction in model.reactions:
            reaction.upper_bound = 1000
            reaction.lower_bound = -1000 if reaction.reversibility else 0

        # Initialize storage
        for reaction in model.reactions:
            reaction_data_fba[model_name][ex][w][reaction.id] = {
                'LB': reaction.lower_bound,
                'UB': reaction.upper_bound,
                'flux': None
            }
            reaction_data_pfba[model_name][ex][w][reaction.id] = {
                'LB': reaction.lower_bound,
                'UB': reaction.upper_bound,
                'flux': None
            }

        # 🔹 Select ONLY matching rates
        one_set = rates[
            (rates.Condition == ex) &
            (rates.Window == w)
        ]

        # Apply constraints
        for _, row in one_set.iterrows():
            metabolite = row.Metabolite
            rate = row.rate_mM_gDCW_h
            sd = row.SD_mM_gDCW_h

            if metabolite not in uptake_names[common_name]:
                continue

            rxn_id = uptake_names[common_name][metabolite]

            if rxn_id not in model.reactions:
                print(f"⚠️ {rxn_id} not in {model_name}")
                continue

            r = model.reactions.get_by_id(rxn_id)

            picked = sorted([rate - sd, rate + sd])
            r.bounds = (picked[0], picked[1])

            reaction_data_fba[model_name][ex][w][rxn_id]['LB'] = picked[0]
            reaction_data_fba[model_name][ex][w][rxn_id]['UB'] = picked[1]
            reaction_data_pfba[model_name][ex][w][rxn_id]['LB'] = picked[0]
            reaction_data_pfba[model_name][ex][w][rxn_id]['UB'] = picked[1]

        # --- Apply NGAM constraint ---
        if common_name in ngam_reactions:
            ngam_id = ngam_reactions[common_name]
            if ngam_id in model.reactions:
                ngam_rxn = model.reactions.get_by_id(ngam_id)
                ngam_lb, ngam_ub = sorted([NGAM_MEAN - NGAM_SD, NGAM_MEAN + NGAM_SD])
                ngam_rxn.bounds = (ngam_lb, ngam_ub)
                print(f"    🔧 NGAM {ngam_id} bounded to [{ngam_lb}, {ngam_ub}]")

                reaction_data_fba[model_name][ex][w][ngam_id]['LB'] = ngam_lb
                reaction_data_fba[model_name][ex][w][ngam_id]['UB'] = ngam_ub
                reaction_data_pfba[model_name][ex][w][ngam_id]['LB'] = ngam_lb
                reaction_data_pfba[model_name][ex][w][ngam_id]['UB'] = ngam_ub
            else:
                print(f"⚠️ NGAM reaction {ngam_id} not in {model_name}")

        # Objective
        biomass_rxn = biomass_reactions[common_name]
        model.objective = biomass_rxn

        try:
            # FBA
            fba_solution = model.optimize()
            mus_fba[model_name][ex][w] = fba_solution.fluxes[biomass_rxn]

            # pFBA
            pfba_solution = pfba(model, fraction_of_optimum=1.0)
            mus_pfba[model_name][ex][w] = pfba_solution.fluxes[biomass_rxn]

            print(f"    ✅ FBA: {mus_fba[model_name][ex][w]:.4f}")
            print(f"    ✅ pFBA: {mus_pfba[model_name][ex][w]:.4f}")

        except cobra.exceptions.Infeasible:
            print(f"    ❌ Infeasible for {model_name}")
            continue

        # Store model for later saving
        models_solved[(model_name, ex, w)] = model.copy()

        # --- MUS (objective values) ---
        mus_records.append({
            'model_name': model_name,
            'common_name': common_name,
            'Condition': ex,
            'Window': w,
            'FBA_mu': fba_solution.fluxes[biomass_rxn],
            'pFBA_mu': pfba_solution.fluxes[biomass_rxn]
        })

        # --- Reaction data ---
        for r in model.reactions:
            reaction_records.append({
                'model_name': model_name,
                'common_name': common_name,
                'Condition': ex,
                'Window': w,
                'Reaction': r.id,
                'FBA_LB': reaction_data_fba[model_name][ex][w][r.id]['LB'],
                'FBA_UB': reaction_data_fba[model_name][ex][w][r.id]['UB'],
                'FBA_Flux': fba_solution.fluxes[r.id],
                'pFBA_LB': reaction_data_pfba[model_name][ex][w][r.id]['LB'],
                'pFBA_UB': reaction_data_pfba[model_name][ex][w][r.id]['UB'],
                'pFBA_Flux': pfba_solution.fluxes[r.id]
            })
# %% Save everything
mus_df = pd.DataFrame(mus_records)

mus_df.to_csv(
    SAVE_FOLDER / "mus_FBA_pFBA.csv",
    index=False
)

print("✅ Saved MUS results")

reaction_df = pd.DataFrame(reaction_records)

reaction_df.to_csv(
    SAVE_FOLDER / "reaction_fluxes_FBA_pFBA.csv",
    index=False
)

print("✅ Saved reaction flux data")

# MODEL_SAVE_FOLDER = SAVE_FOLDER / "models_json"
# MODEL_SAVE_FOLDER.mkdir(exist_ok=True)
#
# for model_name, model in models.items():
#
#     filename = f"{model_name}.json"
#
#     cobra.io.save_json_model(
#         model,
#         MODEL_SAVE_FOLDER / filename
#     )
#
#     print(f"✅ Saved {filename}")


