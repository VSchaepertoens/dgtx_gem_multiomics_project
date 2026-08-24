# %% Section 0: Imports
import pandas as pd
import cobra
from cobra.flux_analysis import pfba
import os
from python.utils.model_utils import load_models, get_models_dir, get_common_model_name
from pathlib import Path

cobra.Configuration().solver = "glpk"

# %% Section 1: Load models
MODELS_DIR = get_models_dir(mode="original")
models, summary_df = load_models(MODELS_DIR)

#%% # Project root is two levels up from python/
PROJECT_ROOT = Path.cwd()
SAVE_FOLDER = PROJECT_ROOT / "analysis" / "models" / "exp_constrained"
SAVE_FOLDER.mkdir(parents=True, exist_ok=True)
print(f"Analysis will be saved to: {SAVE_FOLDER}")

# %% Section 2: Load experimental data
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

# %%
def parse_condition_window(model_name: str):
    """
    Extract Condition and Window from model name.

    Example:
    iCHOv1_Constant_4_constrained -> ("Constant", 4)
    iCHOv1_Temp. shifted_2_constrained -> ("Temp. shifted", 2)
    """
    name = model_name.replace("_constrained", "")

    if "_Constant_" in name:
        parts = name.split("_Constant_")
        return "Constant", int(parts[1])

    elif "_Temp. shifted_" in name:
        parts = name.split("_Temp. shifted_")
        return "Temp. shifted", int(parts[1])

    else:
        print(f"⚠️ Could not parse condition/window from {model_name}")
        return None, None

# %% Section 4: Loop through Models, Conditions, Windows. Set bounds, save models and get FBA/pFBA solutions
# List of models to process
#model_names = ["iCHO2441", "iCHOv1", "iCHO3K_cho_prod_generic_unblocked"]

# Initialize dictionaries for all models
mus_fba = {}
mus_pfba = {}
reaction_data_fba = {}
reaction_data_pfba = {}

# Get unique sets of Experiments and Windows
experiments = set(rates.Condition)
windows = set(rates.Window)

for model_name, base_model in models.items():
    print(f"\nProcessing model {model_name}")

    # 🔹 Extract common model name
    common_name = get_common_model_name(model_name)

    # 🔹 Safety checks
    if common_name not in uptake_names:
        print(f"⚠️ No uptake mapping for {model_name} (common: {common_name}) — skipping")
        continue

    if common_name not in biomass_reactions:
        print(f"⚠️ No biomass reaction for {model_name} (common: {common_name}) — skipping")
        continue

    mus_fba[model_name] = {}
    mus_pfba[model_name] = {}
    reaction_data_fba[model_name] = {}
    reaction_data_pfba[model_name] = {}

    for ex in experiments:
        mus_fba[model_name][ex] = {}
        mus_pfba[model_name][ex] = {}
        reaction_data_fba[model_name][ex] = {}
        reaction_data_pfba[model_name][ex] = {}

        for w in windows:
            # Initialize placeholders
            mus_fba[model_name][ex][w] = None
            mus_pfba[model_name][ex][w] = None
            reaction_data_fba[model_name][ex][w] = {}
            reaction_data_pfba[model_name][ex][w] = {}

            print(f"  Condition {ex}, Window {w}")

            # Fresh copy for each window
            with models[model_name].copy() as model:

                # Reset all reaction bounds
                for reaction in model.reactions:
                    reaction.upper_bound = 1000
                    reaction.lower_bound = -1000 if reaction.reversibility else 0

                # Initialize dictionary for all reactions
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

                # Select rates for this experiment and window
                one_set = rates[(rates.Condition == ex) & (rates.Window == w)]

                # Apply uptake/secretion rates using uptake_names
                for _, row in one_set.iterrows():
                    metabolite = row.Metabolite
                    rate = row.rate_mM_gDCW_h
                    sd = row.SD_mM_gDCW_h
                    if metabolite not in uptake_names[model_name]:
                        continue
                    rxn_id = uptake_names[model_name][metabolite]
                    r = model.reactions.get_by_id(rxn_id)

                    # Set bounds ± SD
                    picked = sorted([rate - sd, rate + sd])
                    r.bounds = (picked[0], picked[1])

                    # Store updated bounds
                    reaction_data_fba[model_name][ex][w][rxn_id]['LB'] = picked[0]
                    reaction_data_fba[model_name][ex][w][rxn_id]['UB'] = picked[1]
                    reaction_data_pfba[model_name][ex][w][rxn_id]['LB'] = picked[0]
                    reaction_data_pfba[model_name][ex][w][rxn_id]['UB'] = picked[1]

                # Set biomass objective dynamically
                biomass_rxn = biomass_reactions[model_name]
                model.objective = biomass_rxn

                try:
                    # FBA
                    fba_solution = model.optimize()
                    print(f"    FBA {biomass_rxn}:", fba_solution.fluxes[biomass_rxn])
                    mus_fba[model_name][ex][w] = fba_solution.fluxes[biomass_rxn]
                    for r in model.reactions:
                        reaction_data_fba[model_name][ex][w][r.id]['flux'] = fba_solution.fluxes[r.id]

                    # pFBA
                    pfba_solution = pfba(model, fraction_of_optimum=1.0)
                    print(f"    pFBA {biomass_rxn}:", pfba_solution.fluxes[biomass_rxn])
                    mus_pfba[model_name][ex][w] = pfba_solution.fluxes[biomass_rxn]
                    for r in model.reactions:
                        reaction_data_pfba[model_name][ex][w][r.id]['flux'] = pfba_solution.fluxes[r.id]

                    # Save rate-constrained model
                    filename_base = f"{model_name}_{ex}_{w}_constrained"
                    cobra.io.save_matlab_model(model, os.path.join(SAVE_FOLDER, filename_base + ".mat"))
                    cobra.io.write_sbml_model(model, os.path.join(SAVE_FOLDER, filename_base + ".xml"))

                    print(f"    ✅ Model saved for {model_name}, Condition {ex}, Window {w}")

                except cobra.exceptions.Infeasible:
                    print(f"    Infeasible solution for {ex} / {w}. Skipping.")
                    continue
