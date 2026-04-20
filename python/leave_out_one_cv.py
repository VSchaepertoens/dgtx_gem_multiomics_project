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
MODELS_DIR = get_models_dir(mode="GIMME")
models, summary_df = load_models(MODELS_DIR)

#%% # Project root is two levels up from python/
PROJECT_ROOT = Path.cwd()
SAVE_FOLDER = PROJECT_ROOT / "analysis" / "models" / "GIMME"
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

# %% Section 5: Cross-validation (leave-one-exchange-out)
mus_records = []
reaction_records = []
models_solved = {}

mus_fba = {}
mus_pfba = {}
reaction_data_fba = {}
reaction_data_pfba = {}

for model_name, base_model in models.items():
    print(f"\nProcessing model {model_name}")

    common_name = get_common_model_name(model_name)
    ex, w = parse_condition_window(model_name)

    if ex is None:
        continue

    print(f"  → Using Condition {ex}, Window {w}")

    if common_name not in uptake_names:
        print(f"⚠️ No uptake mapping for {model_name}")
        continue

    if common_name not in biomass_reactions:
        print(f"⚠️ No biomass reaction for {model_name}")
        continue

    # initialize dicts (UNCHANGED)
    mus_fba[model_name] = {}
    mus_pfba[model_name] = {}
    reaction_data_fba[model_name] = {}
    reaction_data_pfba[model_name] = {}

    mus_fba[model_name][ex] = {}
    mus_pfba[model_name][ex] = {}
    reaction_data_fba[model_name][ex] = {}
    reaction_data_pfba[model_name][ex] = {}

    mus_fba[model_name][ex][w] = {}
    mus_pfba[model_name][ex][w] = {}
    reaction_data_fba[model_name][ex][w] = {}
    reaction_data_pfba[model_name][ex][w] = {}

    with base_model.copy() as model:

        # -----------------------------
        # CLEAN RESET
        # -----------------------------
        for r in model.reactions:
            r.upper_bound = 1000
            r.lower_bound = -1000 if r.reversibility else 0

        # -----------------------------
        # Apply ALL exchange constraints
        # -----------------------------
        one_set = rates[(rates.Condition == ex) & (rates.Window == w)]

        applied_rxns = []

        for _, row in one_set.iterrows():
            met = row.Metabolite
            rate = row.rate_mM_gDCW_h
            sd = row.SD_mM_gDCW_h

            if met not in uptake_names[common_name]:
                continue

            rxn_id = uptake_names[common_name][met]

            if rxn_id not in model.reactions:
                continue

            r = model.reactions.get_by_id(rxn_id)

            lb, ub = sorted([rate - sd, rate + sd])
            measured_flux = rate
            r.bounds = (lb, ub)

            applied_rxns.append((rxn_id, lb, ub))

        # objective
        biomass_rxn = biomass_reactions[common_name]
        model.objective = biomass_rxn

        # -----------------------------
        # BASELINE (all constraints ON)
        # -----------------------------
        try:
            sol_full = model.optimize()
            sol_pfba_full = pfba(model)

        except cobra.exceptions.Infeasible:
            print("❌ Infeasible full model")
            continue

        mu_full = sol_full.fluxes[biomass_rxn]
        mu_pfba_full = sol_pfba_full.fluxes[biomass_rxn]

        # -----------------------------
        # CROSS-VALIDATION LOOP
        # -----------------------------
        for rxn_id, lb, ub in applied_rxns:

            r = model.reactions.get_by_id(rxn_id)

            # save original bounds
            orig_lb, orig_ub = r.lower_bound, r.upper_bound

            # remove constraint (leave-one-out)
            r.upper_bound = 1000
            r.lower_bound = -1000 if r.reversibility else 0


            try:
                sol = model.optimize()
                sol_pfba = pfba(model)

                reaction_records.append({
                    "model_name": model_name,
                    "common_name": common_name,
                    "Condition": ex,
                    "Window": w,
                    "Reaction": rxn_id,

                    "measured_flux": rate,

                    "FBA_LB": lb,
                    "FBA_UB": ub,

                    "FBA_Flux": sol.fluxes[rxn_id],
                    "pFBA_Flux": sol_pfba.fluxes[rxn_id],

                    "mu_FBA_full": mu_full,
                    "mu_pFBA_full": mu_pfba_full,

                    "mu_FBA_LOO": sol.fluxes[biomass_rxn],
                    "mu_pFBA_LOO": sol_pfba.fluxes[biomass_rxn],
                })

            except cobra.exceptions.Infeasible:
                reaction_records.append({
                    "model_name": model_name,
                    "common_name": common_name,
                    "Condition": ex,
                    "Window": w,
                    "Reaction": rxn_id,
                    "error": "infeasible"
                })

            # restore constraint
            r.lower_bound, r.upper_bound = orig_lb, orig_ub

        # -----------------------------
        # store baseline MUS (UNCHANGED)
        # -----------------------------
        mus_records.append({
            "model_name": model_name,
            "common_name": common_name,
            "Condition": ex,
            "Window": w,
            "FBA_mu": mu_full,
            "pFBA_mu": mu_pfba_full
        })

        models_solved[(model_name, ex, w)] = model.copy()

        print(f"✔ Finished {model_name} | μ_FBA={mu_full:.4f} | μ_pFBA={mu_pfba_full:.4f}")

# %% Save everything

# # -----------------------------
# # MUS (baseline biomass only)
# # -----------------------------
# mus_df = pd.DataFrame(mus_records)
#
# mus_df.to_csv(
#     SAVE_FOLDER / "mus_FBA_pFBA_GIMME.csv",
#     index=False
# )
#
# print("✅ Saved MUS results")

# -----------------------------
# Reaction-level cross-validation results
# -----------------------------
reaction_df = pd.DataFrame(reaction_records)

reaction_df.to_csv(
    SAVE_FOLDER / "reaction_fluxes_CV_FBA_pFBA_GIMME.csv",
    index=False
)

print("✅ Saved reaction flux CV data")

# -----------------------------
# OPTIONAL: split baseline vs CV
# -----------------------------
cv_only_df = reaction_df[reaction_df.get("Reaction").notna()]
cv_only_df.to_csv(
    SAVE_FOLDER / "reaction_fluxes_CV_ONLY_GIMME.csv",
    index=False
)

print("✅ Saved CV-only reaction flux data")

# # -----------------------------
# # Save models (unchanged)
# # -----------------------------
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