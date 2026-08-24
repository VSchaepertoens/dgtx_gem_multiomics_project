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
SAVE_FOLDER = PROJECT_ROOT / "analysis" / "models" / "exp_constrained"
SAVE_FOLDER.mkdir(parents=True, exist_ok=True)
print(f"Analysis will be saved to: {SAVE_FOLDER}")

# %% Section 2.1: Load experimental data
rates = pd.read_csv("data/rates/rates_mM_gDCW_h_new.csv")
# Clean column names if needed
rates.columns = [c.strip() for c in rates.columns]

rates = rates[rates.Metabolite != "Titer"] #line added for when igg is the objective

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

# %% Section 2.2: Load experimental data to constrain biomass
growth_df = pd.read_csv("data/rates/growth_rate.csv")
growth_df.columns = [c.strip() for c in growth_df.columns]

# %%
# def parse_condition_window(model_name: str):
#     """
#     Extract Condition and Window from model name.
#
#     Works for:
#     - iCHOv1_Constant_4_constrained
#     - iCHOv1_Temp. shifted_2_constrained
#     - iCHO3K_..._Constant_1_GIMME
#     - future variants with additional suffixes
#
#     Returns
#     -------
#     (condition: str, window: int)
#     """
#
#     # 🔹 Detect condition
#     if "Constant" in model_name:
#         condition = "Constant"
#     elif "Temp. shifted" in model_name:
#         condition = "Temp. shifted"
#     else:
#         print(f"⚠️ Could not detect condition in {model_name}")
#         return None, None
#
#     # 🔹 Extract window (last integer in string)
#     numbers = re.findall(r'\d+', model_name)
#
#     if not numbers:
#         print(f"⚠️ Could not extract window from {model_name}")
#         return condition, None
#
#     window = int(numbers[-1])
#
#     return condition, window

# %% Section 4.2: Loop through Models, Conditions, Windows. Set bounds, save models and get FBA/pFBA solutions
# maximize antibody production under observed growth + nutrient uptake

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

        #Apply experimental biomass constraints
        growth_row = growth_df[
            (growth_df.Condition == ex) &
            (growth_df.Window == w)
            ]

        if growth_row.empty:
            print("⚠️ No growth data found")
            continue

        mu = float(growth_row["avg_growth_µ"].values[0])
        if "sd_growth_µ" in growth_row.columns:
            sd = float(growth_row["sd_growth_µ"].values[0])
        else:
            sd = 0

        biomass_rxn = biomass_reactions[common_name]
        r_bio = model.reactions.get_by_id(biomass_rxn)

        lb, ub = sorted([mu - sd, mu + sd])

        r_bio.lower_bound = lb
        r_bio.upper_bound = ub

        print(f"  ✔ Biomass constrained: {lb:.4f} – {ub:.4f}")

        # Objective
        igg_rxn = igg_reactions[common_name]
        model.objective = igg_rxn

        try:
            # FBA (IGG maximization under growth constraint)
            fba_solution = model.optimize()
            mus_fba[model_name][ex][w] = fba_solution.fluxes[igg_rxn]

            # pFBA
            pfba_solution = pfba(model, fraction_of_optimum=1.0)
            mus_pfba[model_name][ex][w] = pfba_solution.fluxes[igg_rxn]

            print(f"    ✅ FBA: {mus_fba[model_name][ex][w]:.10f}")
            print(f"    ✅ pFBA: {mus_pfba[model_name][ex][w]:.10f}")

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
            'FBA_igg': float(fba_solution.fluxes[igg_rxn]),
            'pFBA_igg': float(pfba_solution.fluxes[igg_rxn])
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
    SAVE_FOLDER / "mus_igg_FBA_pFBA_exp_constrained.csv",
    index=False
)

print("✅ Saved MUS results")

# -----------------------------
# Reaction-level IGG results
# -----------------------------
reaction_df = pd.DataFrame(reaction_records)

# (optional but recommended) ensure consistent ordering
reaction_df = reaction_df.sort_values(
    ["model_name", "Condition", "Window", "Reaction"],
    na_position="last"
)

reaction_df.to_csv(
    SAVE_FOLDER / "reaction_igg_fluxes_FBA_pFBA_exp_constrained.csv",
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