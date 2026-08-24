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
SAVE_FOLDER = PROJECT_ROOT / "models" / "exp_constrained"
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


# %% Section 1.5: Update IgG heavy and light chain reactions for all models

import re

# Define base equations with standard naming (using _L_c format)
HC_BASE_EQUATION = '20.0 ala_L_c + 11.0 arg_L_c + 19.0 asn_L_c + 21.0 asp_L_c + 460.0 atp_c + 11.0 cys_L_c + 16.0 gln_L_c + 20.0 glu_L_c + 25.0 gly_c + 916.0 gtp_c + 917.0 h2o_c + 10.0 his_L_c + 9.0 ile_L_c + 35.0 leu_L_c + 35.0 lys_L_c + 6.0 met_L_c + 15.0 phe_L_c + 37.0 pro_L_c + 49.0 ser_L_c + 41.0 thr_L_c + 10.0 trp_L_c + 16.0 tyr_L_c + 44.0 val_L_c --> adp_c + 459.0 amp_c + 916.0 gdp_c + 917.0 h_c + igg_hc_r + 917.0 pi_c + 459.0 ppi_c'

LC_BASE_EQUATION = '12.0 ala_L_c + 6.0 arg_L_c + 5.0 asn_L_c + 10.0 asp_L_c + 217.0 atp_c + 5.0 cys_L_c + 12.0 gln_L_c + 9.0 glu_L_c + 16.0 gly_c + 430.0 gtp_c + 431.0 h2o_c + 3.0 his_L_c + 6.0 ile_L_c + 14.0 leu_L_c + 14.0 lys_L_c + 2.0 met_L_c + 10.0 phe_L_c + 1.0 pro_L_c + 32.0 ser_L_c + 19.0 thr_L_c + 2.0 trp_L_c + 10.0 tyr_L_c + 15.0 val_L_c --> adp_c + 216.0 amp_c + 430.0 gdp_c + 431.0 h_c + igg_lc_r + 431.0 pi_c + 216.0 ppi_c'


def detect_metabolite_pattern(model):
    """
    Detect the metabolite naming pattern used in the model.
    Returns: 'bracket' for [c] notation, 'underscore' for _c notation, 'double_underscore' for __L_c
    """
    try:
        hc_rxn = model.reactions.get_by_id('igg_hc')
        reaction_str = hc_rxn.reaction

        # Check for [c] notation (iCHO2441 style)
        if re.search(r'[a-z]+_L\[c\]', reaction_str) or re.search(r'gly\[c\]', reaction_str):
            return 'bracket'
        # Check for double underscore (iCHOv1 style)
        elif re.search(r'[a-z]+__L_c', reaction_str):
            return 'double_underscore'
        # Check for single underscore (iCHO3K style)
        elif re.search(r'[a-z]+_L_c', reaction_str):
            return 'underscore'
        else:
            return 'underscore'  # default
    except:
        return 'underscore'  # default


def convert_metabolite_names(equation, pattern):
    """
    Convert metabolite names in the equation to the specified pattern.
    """
    # Define amino acids (with _L suffix)
    amino_acids = ['ala', 'arg', 'asn', 'asp', 'cys', 'gln', 'glu',
                   'his', 'ile', 'leu', 'lys', 'met', 'phe', 'pro',
                   'ser', 'thr', 'trp', 'tyr', 'val']

    # Define other metabolites (without _L suffix)
    other_mets = ['atp', 'gtp', 'adp', 'amp', 'gdp', 'h2o', 'h', 'pi', 'ppi']

    # Special case: glycine (gly) - no _L suffix
    special_mets = {'gly': 'gly'}  # gly is special

    converted = equation

    if pattern == 'bracket':
        # Convert to [c] notation (iCHO2441 style)
        # First handle amino acids with _L suffix
        for aa in amino_acids:
            converted = converted.replace(f'{aa}_L_c', f'{aa}_L[c]')

        # Handle glycine specially
        converted = converted.replace('gly_c', 'gly[c]')

        # Handle other metabolites
        for met in other_mets:
            converted = converted.replace(f'{met}_c', f'{met}[c]')

        # Handle product metabolites
        converted = converted.replace('igg_hc_r', 'igg_hc[r]')
        converted = converted.replace('igg_lc_r', 'igg_lc[r]')

    elif pattern == 'double_underscore':
        # Convert to double underscore (iCHOv1 style)
        for aa in amino_acids:
            converted = converted.replace(f'{aa}_L_c', f'{aa}__L_c')
        # glycine stays as gly_c (no _L)
        # Other metabolites stay as _c

    # For 'underscore' pattern (iCHO3K style), keep as is

    return converted


def verify_metabolites(model, equation, reaction_name):
    """
    Verify that all metabolites in the equation exist in the model.
    Returns list of missing metabolites.
    """
    # Extract metabolite IDs from the equation
    # Pattern matches: number.space.metabolite_id or just metabolite_id
    import re
    # Find all metabolite IDs (words that are not numbers or operators)
    # This is a simplified extraction - looks for patterns like 'ala_L_c' or 'atp[c]'
    met_pattern = r'([a-z][a-z0-9_]*[\[][a-z][\]])|([a-z][a-z0-9_]*_[a-z]_[a-z])|([a-z][a-z0-9_]*_[a-z])|([a-z][a-z0-9_]*(?=[\s\+>]))'

    # More specific: find all metabolite-like strings
    # Look for patterns like XXX_c, XXX[c], XXX__L_c, XXX_L_c
    met_pattern = r'([a-z]+_+[A-Z]?[a-z]?_?[a-z]?\[?[a-z]?\]?)'

    # Actually, let's just split by spaces and look for metabolite patterns
    import re
    # Find all tokens that look like metabolite IDs (contain letters and _ or [])
    tokens = re.findall(r'[a-z][a-z0-9_]*[\[c\]][a-z]?|[a-z][a-z0-9_]*_[a-z]_[a-z]|[a-z][a-z0-9_]*_[a-z]', equation)

    # Also find simple ones like gly_c, atp_c
    simple_mets = re.findall(r'[a-z][a-z0-9_]*_c', equation)
    tokens.extend(simple_mets)

    # Remove duplicates and clean
    tokens = list(set(tokens))

    missing = []
    for token in tokens:
        # Clean token (remove coefficients and spaces)
        token = token.strip()
        if token and not token.isdigit() and token not in ['-->']:
            try:
                model.metabolites.get_by_id(token)
            except KeyError:
                missing.append(token)

    return missing


# Loop through each model and update reactions
for model_name, model in models.items():
    print(f"\n{'=' * 50}")
    print(f"Processing model: {model_name}")
    print(f"{'=' * 50}")

    # Detect the metabolite naming pattern
    pattern = detect_metabolite_pattern(model)
    pattern_names = {
        'bracket': 'bracket notation [c]',
        'double_underscore': 'double underscore (__L_c)',
        'underscore': 'single underscore (_L_c)'
    }
    print(f"Detected metabolite pattern: {pattern_names.get(pattern, pattern)}")

    # Convert equations to the model's pattern
    hc_equation = convert_metabolite_names(HC_BASE_EQUATION, pattern)
    lc_equation = convert_metabolite_names(LC_BASE_EQUATION, pattern)

    # Verify metabolites before updating
    print("\nVerifying HC reaction metabolites...")
    missing_hc = verify_metabolites(model, hc_equation, 'HC')
    if missing_hc:
        print(f"  WARNING: Missing metabolites in HC equation: {missing_hc}")
    else:
        print("  All HC metabolites exist in model ✓")

    print("\nVerifying LC reaction metabolites...")
    missing_lc = verify_metabolites(model, lc_equation, 'LC')
    if missing_lc:
        print(f"  WARNING: Missing metabolites in LC equation: {missing_lc}")
    else:
        print("  All LC metabolites exist in model ✓")

    try:
        # Update heavy chain reaction
        print(f"\nOriginal HC reaction in {model_name}:")
        print(model.reactions.igg_hc.reaction)

        existing_reaction = model.reactions.get_by_id('igg_hc')
        existing_reaction.reaction = hc_equation

        print(f"\nUpdated HC reaction in {model_name}:")
        print(model.reactions.igg_hc.reaction)

    except KeyError as e:
        print(f"\nWarning: 'igg_hc' reaction not found in {model_name}: {e}")

    try:
        # Update light chain reaction
        print(f"\nOriginal LC reaction in {model_name}:")
        print(model.reactions.igg_lc.reaction)

        existing_reaction = model.reactions.get_by_id('igg_lc')
        existing_reaction.reaction = lc_equation

        print(f"\nUpdated LC reaction in {model_name}:")
        print(model.reactions.igg_lc.reaction)

    except KeyError as e:
        print(f"\nWarning: 'igg_lc' reaction not found in {model_name}: {e}")

    # Check for unknown metabolites after update
    print("\nChecking for newly created metabolites...")
    try:
        # Look for any metabolites with '_c' suffix that shouldn't exist in bracket models
        if pattern == 'bracket':
            hc_rxn = model.reactions.get_by_id('igg_hc')
            # Check if any _c metabolites were created
            import re

            c_metabolites = re.findall(r'[a-z][a-z0-9_]*_c', hc_rxn.reaction)
            if c_metabolites:
                print(f"  WARNING: Found _c metabolites in bracket model: {c_metabolites}")
            else:
                print("  No _c metabolites found (good!) ✓")
    except:
        pass

print(f"\n{'=' * 50}")
print("All models have been updated with the new IgG production reactions")
print(f"{'=' * 50}")

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
