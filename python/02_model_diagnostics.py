# %% Imports
import cobra
import pandas as pd
from python.utils.model_utils import load_models, get_models_dir

from cobra.flux_analysis import find_blocked_reactions
from cobra.flux_analysis import flux_variability_analysis
from cobra.medium import minimal_medium

import matplotlib.pyplot as plt
import seaborn as sns

cobra.Configuration().solver = "glpk"

# %% Hard-coded save directories
from pathlib import Path

# %% Section 1:  Load models
MODELS_DIR = get_models_dir(mode="exp_constrained")
models, summary_df = load_models(MODELS_DIR)

# %% Section 2: Initialize results containers
energy_results_minimal = []
energy_results_minimal_exchanges = []
energy_results_transporters = []
exchange_results = []
transporter_results = {}

# %% Section 3: Define key reactions for each model
biomass_reactions = {
    "iCHO3K_cho_prod_generic_unblocked": "biomass_cho_prod",
    "iCHO2441": "biomass_cho_prod",
    "iCHOv1": "BIOMASS_cho_producing"
}

igg_reactions = {
    "iCHO2441": "igg_formation",
    "iCHO3K_cho_prod_generic_unblocked": "igg_formation",
    "iCHOv1": "igg_formation"
}


oxygen_exchanges = {
    "iCHO3K_cho_prod_generic_unblocked": "EX_o2_e",
    "iCHO2441": "EX_o2(e)",
    "iCHOv1": "EX_o2_e"
}

glucose_exchanges = {
    "iCHO2441": "EX_glc(e)",
    "iCHO3K_cho_prod_generic_unblocked": "EX_glc_e",
    "iCHOv1": "EX_glc__D_e"
}

atp_reactions = {
    "iCHO2441": "DM_atp[c]",
    "iCHO3K_cho_prod_generic_unblocked": "ATPM",
    "iCHOv1": "ATPM"
}

# Define true glucose and oxygen transporters for each model
glucose_oxygen_transporters = {
    "iCHO2441": {
        "glucose": ["GLCt1r", "GLCt2r", "GLCter", "GLCNte"],
        "oxygen": ["O2t", "O2ter", "O2St"]
    },
    "iCHOv1": {
        "glucose": ["GLCter"],  # Only this one is actually present
        "oxygen": ["O2t", "O2ter", "O2St"]
    },
    "iCHO3K_cho_prod_generic_unblocked": {
        "glucose": ["GLCt1r", "GLCt2r", "GLCter", "GLCNte"],
        "oxygen": ["O2t", "O2ter", "O2St"]
    }
}
# %% Function to map contextualized filenames to common names
def get_common_model_name(model_name: str) -> str:
    """
    Strips suffixes like '_Constant_X_GIMME' to return the original/common model name
    """
    # Remove suffix starting with '_Constant'
    if "_Constant" in model_name:
        common_name = model_name.split("_Constant")[0]
    elif "_Temp. shifted" in model_name:
        common_name = model_name.split("_Temp. shifted")[0]
    elif "_GIMME" in model_name:
        common_name = model_name.split("_GIMME")[0]
    else:
        common_name = model_name
    return common_name


# # Example usage within diagnostics loop:
# for model_name, model in models.items():
#     common_name = get_common_model_name(model_name)
#
#     biomass_rxn = biomass_reactions.get(common_name)
#     o2_rxn = oxygen_exchanges.get(common_name)
#     glc_rxn = glucose_exchanges.get(common_name)
#     atp_rxn = atp_reactions.get(common_name)
#     transporters = glucose_oxygen_transporters.get(common_name, {})
#
#     if biomass_rxn is None:
#         print(f"Warning: No biomass reaction defined for model {model_name} (common: {common_name})")
#     if o2_rxn is None or glc_rxn is None:
#         print(f"Warning: Missing exchange reactions for {model_name} (common: {common_name})")
#
#     # Now you can use biomass_rxn, o2_rxn, glc_rxn, atp_rxn, transporters safely

# %% Function to map to common model name
# Check presence of key reactions
for model_name, model in models.items():
    common_name = get_common_model_name(model_name)

    o2_rxn = oxygen_exchanges.get(common_name)
    glc_rxn = glucose_exchanges.get(common_name)
    atp_rxn = atp_reactions.get(common_name)

    print(f"\n--- {model_name} ---")

    # Oxygen
    if o2_rxn:
        print(f"O2 exchange ({o2_rxn}):", o2_rxn in model.reactions)
    else:
        print("O2 exchange: not defined")

    # Glucose
    if glc_rxn:
        print(f"Glucose exchange ({glc_rxn}):", glc_rxn in model.reactions)
    else:
        print("Glucose exchange: not defined")

    # ATP
    if atp_rxn:
        print(f"ATP reaction ({atp_rxn}):", atp_rxn in model.reactions)
    else:
        print("ATP reaction: not defined")

# %% # Debugging the compartments issue
# for model_name, model in models.items():
#     print(f"\nModel: {model_name}")
#     print("Compartments:", model.compartments)
#     # Optional: count metabolites per compartment
#     compartment_counts = {}
#     for met in model.metabolites:
#         compartment_counts[met.compartment] = compartment_counts.get(met.compartment, 0) + 1
#     print("Metabolites per compartment:", compartment_counts)
#
# # Pick a model (adjust key if needed)
# model_name = "iCHO3K_cho_prod_generic_unblocked"
#
# # Original vs contextualized
# model_orig = models[model_name]
# model_ctx = models_contextualized[[k for k in models_contextualized if model_name in k][0]]
#
# print("\nContextualized metabolite IDs (first 20):")
# print([met.id for met in model_ctx.metabolites[:20]])
#
# print("\nContextualized metabolite compartments (unique):")
# print(set([met.compartment for met in model_ctx.metabolites]))

#%%
def reconstruct_compartments(model):
    """
    Reconstruct metabolite compartments from ID suffixes (e.g. atp_c, glc__D_m)
    and update model.compartments accordingly.
    """
    compartments = {}

    for met in model.metabolites:
        if "_" in met.id:
            comp = met.id.split("_")[-1]
        else:
            comp = "c"  # fallback if no suffix

        # Assign compartment to metabolite
        met.compartment = comp

        # Collect compartments
        compartments[comp] = ""

    # Update model compartments dictionary
    model.compartments = compartments

    return model

# List of model name substrings to apply reconstruction to
models_to_fix = ["iCHO3K_cho_prod_generic_unblocked", "iCHOv1"]

for model_name, model in models.items():
    if any(name in model_name for name in models_to_fix):
        reconstruct_compartments(model)
        print(f"Reconstructed compartments for {model_name}")

for model_name, model in models.items():
    print(f"\n--- {model_name} ---")

    print("Unique metabolite compartments:")
    print(set([met.compartment for met in model.metabolites]))

    print("Model compartments dictionary:")
    print(model.compartments)
    print(len(model.compartments))


# %% Section 4.1: Minimal ATP flux and exchange reactions
for model_name, model in models.items():
    print(f"\nRunning minimal ATP diagnostics for {model_name}...")
    common_name = get_common_model_name(model_name)

    # Basic statistics
    num_reactions = len(model.reactions)
    num_metabolites = len(model.metabolites)
    num_genes = len(model.genes)


    # Identify exchange reactions
    ex_reactions = [r.id for r in model.reactions if r.id.startswith("EX_")]
    exchange_results.append({
        "model_name": model_name,
        "num_exchanges": len(ex_reactions),
        "exchange_ids": ex_reactions
    })

    # Energy metabolism test: Maximize ATP with all exchanges closed
    with model:
        for rxn in model.exchanges:
            rxn.bounds = (0, 0)

            # # --- sanity check: print exchange reaction bounds ---
            # print(f"\nExchange reaction bounds for {model_name} after closure:")
            # for rxn in model.exchanges:
            #     print(f"{rxn.id}: lb={rxn.lower_bound}, ub={rxn.upper_bound}")

        #atp_rxn_id = atp_reactions.get(model_name)
        atp_rxn_id = atp_reactions.get(common_name)
        if atp_rxn_id and atp_rxn_id in model.reactions:
            model.objective = model.reactions.get_by_id(atp_rxn_id)
            sol = model.optimize()
            atp_flux = sol.objective_value
            print(f"{model_name}: Max ATP flux (minimal) = {atp_flux}")
        else:
            print(f"{model_name}: ATP maintenance reaction not found!")
            atp_flux = None

        energy_results_minimal.append({
            "model_name": model_name,
            "atp_flux_minimal": atp_flux
        })

# %% Section 4.2: ATP flux with O2 and glucose (including glucose uptake)
for model_name, model in models.items():
    print(f"\nRunning ATP diagnostics with O2 + glucose for {model_name}...")
    common_name = get_common_model_name(model_name)

    with model:
        # Close all exchange reactions
        for rxn in model.exchanges:
            rxn.bounds = (0, 0)

        # Open oxygen exchange
        o2_rxn_id = oxygen_exchanges.get(common_name)
        if o2_rxn_id in model.reactions:
            rxn = model.reactions.get_by_id(o2_rxn_id)
            rxn.lower_bound = -1000
            rxn.upper_bound = 1000
            print(f"{model_name}: Opened oxygen exchange -> {rxn.id}")

        # Open glucose exchange
        glc_rxn_id = glucose_exchanges.get(common_name)
        if glc_rxn_id in model.reactions:
            rxn = model.reactions.get_by_id(glc_rxn_id)
            rxn.lower_bound = -10
            rxn.upper_bound = 1000
            print(f"{model_name}: Opened glucose exchange -> {rxn.id}")

        # Maximize ATP
        atp_rxn_id = atp_reactions.get(common_name)
        if atp_rxn_id and atp_rxn_id in model.reactions:
            model.objective = model.reactions.get_by_id(atp_rxn_id)
            sol = model.optimize()
            atp_flux = sol.objective_value
            print(f"{model_name}: Max ATP flux (O2 + glucose) = {atp_flux}")
        else:
            atp_flux = None

        print("ATP reaction present:", atp_rxn_id in model.reactions)

        print("ATP reaction flux bounds:",
              model.reactions.get_by_id(atp_rxn_id).bounds)

        print("Number of reactions:", len(model.reactions))

        model.reactions.get_by_id(atp_rxn_id).summary()

        energy_results_transporters.append({
            "model_name": model_name,
            "atp_flux_glc_o2": atp_flux
        })

# %% Section 4.3: ATP flux with transporters (all models)
for model_name, model in models.items():
    common_name = get_common_model_name(model_name)

    glucose_tx = glucose_oxygen_transporters.get(common_name, {}).get("glucose", [])
    oxygen_tx = glucose_oxygen_transporters.get(common_name, {}).get("oxygen", [])

    if not glucose_tx and not oxygen_tx:
        print(f"{model_name}: No transporters found, skipping flux with transporters")
        continue

    with model:  # Use a context manager so changes are temporary
        # Close all exchange reactions first
        for rxn in model.exchanges:
            rxn.bounds = (0, 0)

        print(f"\n{model_name}: All exchange reactions closed.")

        # Open exchange reactions for glucose and oxygen
        glc_ex_id = glucose_exchanges.get(common_name)
        o2_ex_id = oxygen_exchanges.get(common_name)

        if glc_ex_id in model.reactions:
            rxn = model.reactions.get_by_id(glc_ex_id)
            rxn.lower_bound = -1000
            rxn.upper_bound = 1000
            print(f"{model_name}: Opened glucose exchange -> {rxn.id}")

        if o2_ex_id in model.reactions:
            rxn = model.reactions.get_by_id(o2_ex_id)
            rxn.lower_bound = -1000
            rxn.upper_bound = 1000
            print(f"{model_name}: Opened oxygen exchange -> {rxn.id}")

        # Open only the specified glucose and oxygen transporters
        for rxn_id in glucose_tx + oxygen_tx:
            if rxn_id in model.reactions:
                rxn = model.reactions.get_by_id(rxn_id)
                rxn.lower_bound = -1000
                rxn.upper_bound = 1000
                print(f"{model_name}: Opened transporter -> {rxn.id}")
            else:
                print(f"{model_name}: Transporter {rxn_id} NOT FOUND in model")

        # Set ATP reaction as objective
        atp_rxn_id = atp_reactions.get(common_name)
        if atp_rxn_id in model.reactions:
            model.objective = model.reactions.get_by_id(atp_rxn_id)
            sol = model.optimize()
            atp_flux = sol.objective_value
            print(f"{model_name}: Max ATP flux (with transporters) = {atp_flux}")
        else:
            atp_flux = None
            print(f"{model_name}: ATP reaction {atp_rxn_id} NOT FOUND!")

        # Store the result
        energy_results_transporters.append({
            "model_name": model_name,
            "atp_flux_with_transporters": atp_flux
        })

# %% Section 4.4: Minimal ATP flux and all exchange reactions open
for model_name, model in models.items():
    print(f"\nRunning minimal ATP diagnostics for {model_name}...")
    common_name = get_common_model_name(model_name)

    # Identify exchange reactions
    ex_reactions = [r.id for r in model.reactions if r.id.startswith("EX_")]
    exchange_results.append({
        "model_name": model_name,
        "num_exchanges": len(ex_reactions),
        "exchange_ids": ex_reactions
    })

    # Energy metabolism test: Maximize ATP with all exchanges open
    with model:
        for rxn in model.exchanges:
            rxn.bounds = (-1000, 1000)

        atp_rxn_id = atp_reactions.get(common_name)
        if atp_rxn_id and atp_rxn_id in model.reactions:
            model.objective = model.reactions.get_by_id(atp_rxn_id)
            sol = model.optimize()
            atp_flux = sol.objective_value
            print(f"{model_name}: Max ATP flux (all exchanges) = {atp_flux}")
        else:
            print(f"{model_name}: ATP maintenance reaction not found!")
            atp_flux = None

        energy_results_minimal_exchanges.append({
            "model_name": model_name,
            "atp_flux_all_exchanges": atp_flux
        })

# %% Section 5 : Blocked reaction analysis
cobra.Configuration().solver = "gurobi"
blocked_results = []

for model_name, model in models.items():
    print(f"\nAnalyzing blocked reactions in {model_name}...")

    blocked_rxns = find_blocked_reactions(model)

    print(f"{model_name}: {len(blocked_rxns)} blocked reactions")

    blocked_results.append({
        "model_name": model_name,
        "num_blocked_reactions": len(blocked_rxns),
        "blocked_reactions": blocked_rxns
    })

blocked_df = pd.DataFrame(blocked_results)

print("\nBlocked reaction summary:")
print(blocked_df[["model_name", "num_blocked_reactions"]])

# Is ATP maintenance reaction blocked?
for model_name, model in models.items():
    #atp_rxn = atp_reactions[common_name]
    atp_rxn = atp_reactions.get(common_name)

    if atp_rxn in blocked_df.loc[
        blocked_df["model_name"] == model_name, "blocked_reactions"
    ].values[0]:
        print(f"{model_name}: ATP reaction is BLOCKED")
    else:
        print(f"{model_name}: ATP reaction is NOT blocked")

# %% Section 6: Flux variability analysis for ATP (This section is now commented and moved after the section 8
# fva_results = []
#
# for model_name, model in models.items():
#     common_name = get_common_model_name(model_name)
#
#     atp_rxn_id = atp_reactions.get(common_name)
#
#     if atp_rxn_id not in model.reactions:
#         print(f"{model_name}: ATP reaction not found")
#         continue
#
#     print(f"\nRunning FVA for ATP in {model_name}...")
#
#     fva = flux_variability_analysis(
#         model,
#         reaction_list=[atp_rxn_id],
#         fraction_of_optimum=1
#     )
#
#     min_flux = fva.loc[atp_rxn_id, "minimum"]
#     max_flux = fva.loc[atp_rxn_id, "maximum"]
#
#     print(f"{model_name}: ATP flux range = [{min_flux}, {max_flux}]")
#
#     fva_results.append({
#         "model_name": model_name,
#         "ATP_min_flux": min_flux,
#         "ATP_max_flux": max_flux
#     })
#
# fva_df = pd.DataFrame(fva_results)
#
# print("\nATP flux variability:")
# print(fva_df)


 # %% Section 7: Minimal medium required for ATP production
#cobra.Configuration().solver = "glpk"
minimal_medium_results = []

for model_name, model in models.items():

    print(f"\nComputing minimal medium for ATP production in {model_name}...")
    common_name = get_common_model_name(model_name)

    with model:

        # Close all exchanges first
        for rxn in model.exchanges:
            rxn.lower_bound = -1000
            rxn.upper_bound = 1000

        # Set ATP reaction as objective
        atp_rxn_id = atp_reactions.get(common_name)

        if atp_rxn_id not in model.reactions:
            print(f"{model_name}: ATP reaction not found")
            continue

        model.objective = model.reactions.get_by_id(atp_rxn_id)

        # Calculate minimal medium
        medium = minimal_medium(model, 1)

        if medium is None:
            print(f"{model_name}: minimal medium not found")

            minimal_medium_results.append({
                "model_name": model_name,
                "required_exchanges": None
            })

        else:
            print(medium)

            minimal_medium_results.append({
                "model_name": model_name,
                "required_exchanges": list(medium.index),
                "fluxes": medium.to_dict()
            })

minimal_media_dict = {}
for entry in minimal_medium_results:
    if entry["required_exchanges"] is not None:
        # Use a pandas Series for easier iteration later
        minimal_media_dict[entry["model_name"]] = pd.Series(entry["fluxes"])
    else:
        minimal_media_dict[entry["model_name"]] = None

# %% Section 8: Max ATP flux with minimal medium open + glucose  + oxygen
minimal_medium_flux_results = []

for model_name, model in models.items():
    print(f"\nComputing ATP flux with minimal medium for {model_name}...")
    common_name = get_common_model_name(model_name)

    # Retrieve minimal medium for this model
    medium = minimal_media_dict.get(model_name)
    if medium is None or medium.empty:
        print(f"{model_name}: Minimal medium not available, skipping...")
        minimal_medium_flux_results.append({
            "model_name": model_name,
            "ATP_max_flux_minimal_medium": None
        })
        continue

    with model:
        # Close all exchange reactions first
        for rxn in model.exchanges:
            rxn.bounds = (0, 0)

        # Open only the reactions in the minimal medium
        for rxn_id in medium.index:
            if rxn_id in model.reactions:
                rxn = model.reactions.get_by_id(rxn_id)
                rxn.lower_bound = -1000
                rxn.upper_bound = 1000
                print(f"{model_name}: Opened minimal medium exchange -> {rxn.id}")
            else:
                print(f"{model_name}: Minimal medium reaction {rxn_id} NOT FOUND")

        # Open minimal medium + glucose + oxygen
        for rxn_id in medium.index.tolist() + [glucose_exchanges[common_name], oxygen_exchanges[common_name]]:
            if rxn_id in model.reactions:
                rxn = model.reactions.get_by_id(rxn_id)
                rxn.lower_bound = -1000
                rxn.upper_bound = 1000
                print(f"{model_name}: Opened minimal medium + glucose/O2 exchange -> {rxn.id}")
            else:
                print(f"{model_name}: Reaction {rxn_id} NOT FOUND")

        # Maximize ATP
        atp_rxn_id = atp_reactions.get(common_name)
        if atp_rxn_id in model.reactions:
            model.objective = model.reactions.get_by_id(atp_rxn_id)
            sol = model.optimize()
            atp_flux = sol.objective_value
            print(f"{model_name}: Max ATP flux (minimal medium) = {atp_flux}")
        else:
            atp_flux = None
            print(f"{model_name}: ATP reaction {atp_rxn_id} NOT FOUND")

        minimal_medium_flux_results.append({
            "model_name": model_name,
            "ATP_max_flux_minimal_medium": atp_flux
        })

# %% Section 6: ATP FVA with minimal medium open + glucose  + oxygen
fva_minimal_medium_results = []

for model_name, model in models.items():
    print(f"\nRunning FVA for ATP with minimal medium in {model_name}...")
    common_name = get_common_model_name(model_name)

    # Retrieve minimal medium for this model
    medium = minimal_media_dict.get(model_name)
    if medium is None or medium.empty:
        print(f"{model_name}: Minimal medium not available, skipping...")
        fva_minimal_medium_results.append({
            "model_name": model_name,
            "ATP_min_flux": None,
            "ATP_max_flux": None
        })
        continue

    with model:
        # Close all exchange reactions first
        for rxn in model.exchanges:
            rxn.bounds = (0, 0)

        # Open minimal medium exchanges
        for rxn_id in medium.index:
            if rxn_id in model.reactions:
                rxn = model.reactions.get_by_id(rxn_id)
                rxn.lower_bound = -1000
                rxn.upper_bound = 1000

        # Open minimal medium + glucose + oxygen
        for rxn_id in medium.index.tolist() + [glucose_exchanges[common_name], oxygen_exchanges[common_name]]:
            if rxn_id in model.reactions:
                rxn = model.reactions.get_by_id(rxn_id)
                rxn.lower_bound = -1000
                rxn.upper_bound = 1000
                print(f"{model_name}: Opened minimal medium + glucose/O2 exchange -> {rxn.id}")
            else:
                print(f"{model_name}: Reaction {rxn_id} NOT FOUND")

        # Set ATP objective
        atp_rxn_id = atp_reactions.get(common_name)
        if atp_rxn_id not in model.reactions:
            print(f"{model_name}: ATP reaction not found")
            fva_minimal_medium_results.append({
                "model_name": model_name,
                "ATP_min_flux": None,
                "ATP_max_flux": None
            })
            continue

        fva = flux_variability_analysis(
            model,
            reaction_list=[atp_rxn_id],
            fraction_of_optimum=1.0  # Use full ATP optimum
        )

        min_flux = fva.loc[atp_rxn_id, "minimum"]
        max_flux = fva.loc[atp_rxn_id, "maximum"]

        print(f"{model_name}: ATP flux range = [{min_flux}, {max_flux}]")

        fva_minimal_medium_results.append({
            "model_name": model_name,
            "ATP_min_flux": min_flux,
            "ATP_max_flux": max_flux
        })

# Convert to DataFrame
fva_minimal_medium_df = pd.DataFrame(fva_minimal_medium_results)
print("\nATP flux variability with minimal medium:")
print(fva_minimal_medium_df)

# %% Section 9: Biomass feasibility check
growth_results = []

for model_name, model in models.items():
    common_name = get_common_model_name(model_name)

    biomass_rxn_id = biomass_reactions[common_name]
    biomass_rxn = model.reactions.get_by_id(biomass_rxn_id)

    with model:
        model.objective = biomass_rxn

        solution = model.optimize()
        growth_rate = solution.objective_value

    growth_results.append({
        "model": model_name,
        "biomass_rxn": biomass_rxn_id,
        "growth_rate": growth_rate,
        "feasible": growth_rate > 1e-6
    })

growth_df = pd.DataFrame(growth_results)

print(growth_df)

# iCHO3K_cho_prod_generic_unblocked has glucose uptake unconstrained (lower_bound = 0, unconstrained = True).
#
# All other models have proper negative lower bounds for glucose (-10 or similar), which allows uptake.
#
# Oxygen is fine in all models.
#
# So the biomass infeasibility in iCHO3K is purely a medium/uptake issue, not a broken network.

# %% Section 10: Biomass precursor feasibility check
"""
Check whether all metabolites consumed by the biomass reaction can be produced by the model.

Conceptually:
- The biomass reaction is like a recipe for growth: it consumes precursors (amino acids, nucleotides, ATP, etc.).
- Growth fails if any single precursor cannot be produced.
- To diagnose this, we test each precursor individually:
    1. Add a temporary 'demand reaction' for the precursor (a sink allowing the metabolite to be exported).
    2. Maximize this demand using FBA to see if the network can produce the metabolite.
    3. Record whether the precursor is producible (max flux > 0).
- This identifies exactly which metabolites block biomass formation without altering the original model.

Outputs:
- A DataFrame listing, for each model and precursor: max flux and producibility (True/False).
"""
def check_biomass_precursors(models, biomass_rxns):
    """
    Check producibility of all metabolites consumed by the biomass reaction in each model.

    Args:
        models (dict): model_name -> cobra.Model
        biomass_rxns (dict): model_name -> biomass reaction ID

    Returns:
        pd.DataFrame: model, metabolite, max_flux, producible
    """
    precursor_results = []

    for model_name, model in models.items():
        common_name = get_common_model_name(model_name)
        biomass_rxn_id = biomass_rxns.get(common_name)
        biomass_rxn = model.reactions.get_by_id(biomass_rxn_id)

        # metabolites consumed in biomass (coeff < 0)
        biomass_precursors = [m for m, coeff in biomass_rxn.metabolites.items() if coeff < 0]
        print(f"{model_name}: Consumed by biomass = {biomass_precursors}")

        for met in biomass_precursors:
            with model:
                # add a temporary demand reaction
                demand_rxn = model.add_boundary(met, type="demand")
                model.objective = demand_rxn

                sol = model.optimize()

                precursor_results.append({
                    "model": model_name,
                    "metabolite": met.id,
                    "max_flux": sol.objective_value,
                    "producible": sol.objective_value > 1e-9
                })

    return pd.DataFrame(precursor_results)

# run the check
precursor_df = check_biomass_precursors(models, biomass_reactions)

# quickly inspect problematic precursors
problem_precursors = precursor_df[precursor_df["producible"] == False]
print("\nBlocked biomass precursors:")
print(problem_precursors)


# # Run the precursor check with temporary glucose/oxygen uptake
# icho3k = models["iCHO3K_cho_prod_generic_unblocked"]
#
# with icho3k:  # Changes inside this block are temporary
#     # Open glucose and oxygen uptake
#     icho3k.reactions.get_by_id("EX_glc_e").lower_bound = -1000
#     icho3k.reactions.get_by_id("EX_glc_e").upper_bound = 1000
#     icho3k.reactions.get_by_id("EX_o2_e").lower_bound = -1000
#     icho3k.reactions.get_by_id("EX_o2_e").upper_bound = 1000
#
#     # Run the precursor check
#     precursor_df = check_biomass_precursors(
#         {"iCHO3K_cho_prod_generic_unblocked": icho3k},
#         {"iCHO3K_cho_prod_generic_unblocked": biomass_reactions["iCHO3K_cho_prod_generic_unblocked"]}
#     )
#
# # Inspect problematic precursors
# problem_precursors = precursor_df[precursor_df["producible"] == False]
# print("\nBlocked biomass precursors:")
# print(problem_precursors)

# %% Section 11: Antibody (IgG) production feasibility check
igg_results = []

for model_name, model in models.items():
    common_name = get_common_model_name(model_name)
    igg_rxn_id = igg_reactions.get(common_name)

    if igg_rxn_id not in model.reactions:
        print(f"{model_name}: IgG reaction {igg_rxn_id} not found, skipping...")
        continue

    igg_rxn = model.reactions.get_by_id(igg_rxn_id)

    with model:
        model.objective = igg_rxn
        sol = model.optimize()
        igg_flux = sol.objective_value

    igg_results.append({
        "model_name": model_name,
        "igg_rxn": igg_rxn_id,
        "igg_flux": igg_flux,
        "feasible": igg_flux > 1e-6
    })

igg_df = pd.DataFrame(igg_results)
print("\nIgG production feasibility:")
print(igg_df)

# %% Section 12: Antibody (IgG) precursor feasibility check
def check_igg_precursors(models, igg_reactions):
    """
    Check producibility of all metabolites consumed by the IgG reaction in each model.

    Args:
        models (dict): model_name -> cobra.Model
        igg_reactions (dict): common_name -> IgG reaction ID

    Returns:
        pd.DataFrame: model, metabolite, max_flux, producible
    """
    precursor_results = []

    for model_name, model in models.items():
        common_name = get_common_model_name(model_name)
        igg_rxn_id = igg_reactions.get(common_name)

        if igg_rxn_id not in model.reactions:
            print(f"{model_name}: IgG reaction {igg_rxn_id} not found, skipping...")
            continue

        igg_rxn = model.reactions.get_by_id(igg_rxn_id)

        # metabolites consumed in IgG reaction (coeff < 0)
        igg_precursors = [m for m, coeff in igg_rxn.metabolites.items() if coeff < 0]
        print(f"{model_name}: Consumed by IgG = {[m.id for m in igg_precursors]}")

        for met in igg_precursors:
            with model:
                # add a temporary demand reaction
                demand_rxn = model.add_boundary(met, type="demand")
                model.objective = demand_rxn
                sol = model.optimize()

                precursor_results.append({
                    "model": model_name,
                    "metabolite": met.id,
                    "max_flux": sol.objective_value,
                    "producible": sol.objective_value > 1e-9
                })

    return pd.DataFrame(precursor_results)

# run the check
igg_precursor_df = check_igg_precursors(models, igg_reactions)

# # quickly inspect problematic precursors
# problem_igg_precursors = igg_precursor_df[igg_precursor_df["producible"] == False]
# print("\nBlocked IgG precursors:")
# print(problem_igg_precursors)
# %% Section 13: Compile final comparison table
comparison_results = []

# Create lookup dicts for easy merging
energy_min_dict = {r["model_name"]: r["atp_flux_minimal"] for r in energy_results_minimal}
energy_min_ex_dict = {
    r["model_name"]: r.get("atp_flux_all_exchanges", None)
    for r in energy_results_minimal_exchanges
}
# energy_min_ex_dict = {r["model_name"]: r["atp_flux_all_exchanges"] for r in energy_results_minimal_exchanges}
energy_glc_o2_dict = {r["model_name"]: r["atp_flux_glc_o2"] for r in energy_results_transporters if "atp_flux_glc_o2" in r}
energy_transporter_dict = {r["model_name"]: r["atp_flux_with_transporters"] for r in energy_results_transporters if "atp_flux_with_transporters" in r}
#blocked_dict = {r["model_name"]: r["num_blocked_reactions"] for r in blocked_results}
#blocked_rxn_ids = {r["model_name"]: r["blocked_reactions"] for r in blocked_results}
fva_min_dict = fva_minimal_medium_df.set_index("model_name")["ATP_min_flux"].to_dict()
fva_max_dict = fva_minimal_medium_df.set_index("model_name")["ATP_max_flux"].to_dict()
minimal_medium_flux_dict = {r["model_name"]: r["ATP_max_flux_minimal_medium"] for r in minimal_medium_flux_results}
# Minimal medium info
minimal_medium_series = minimal_media_dict.get(model_name)
if minimal_medium_series is not None:
    minimal_medium_reactions = minimal_medium_series.index.tolist()
    minimal_medium_fluxes = minimal_medium_series.to_dict()
    # Optional: formatted string for CSV readability
    minimal_medium_str = "; ".join(f"{rxn}:{flux:.6f}" for rxn, flux in minimal_medium_fluxes.items())
else:
    minimal_medium_reactions = []
    minimal_medium_fluxes = {}
    minimal_medium_str = None

# Number of blocked precursors per model
blocked_precursors = (
    precursor_df[precursor_df["producible"] == False]
    .groupby("model")
    .size()
    .to_dict()
)

growth_feasible_dict = {r["model"]: r["feasible"] for r in growth_results}
growth_rate_dict = {r["model"]: r["growth_rate"] for r in growth_results}
# Fraction of blocked precursors (optional)
total_precursors = precursor_df.groupby("model").size().to_dict()
blocked_fraction = {m: blocked_precursors.get(m, 0)/total_precursors[m] for m in models.keys()}
# IgG feasibility
igg_feasible_dict = {r["model_name"]: r["feasible"] for r in igg_results}
igg_rate_dict = {r["model_name"]: r["igg_flux"] for r in igg_results}
# IgG precursors (if checked)
igg_precursor_dict = {r["model"]: r["producible"] for r in igg_precursor_df.to_dict("records")} if 'igg_precursor_df' in globals() else {}



for model_name, model in models.items():
    # Basic model stats
    num_reactions = len(model.reactions)
    num_metabolites = len(model.metabolites)
    num_genes = len(model.genes)
    ex_entry = next((r for r in exchange_results if r["model_name"] == model_name), {})
    num_exchanges = ex_entry.get("num_exchanges", 0)

    comparison_results.append({
        "model_name": model_name,
        "num_reactions": num_reactions,
        "num_metabolites": num_metabolites,
        "num_genes": num_genes,
        "num_exchanges": num_exchanges,
        #"num_blocked_reactions": blocked_dict.get(model_name),
        #"blocked_reaction_ids": blocked_rxn_ids.get(model_name, []),
        "atp_flux_minimal": energy_min_dict.get(model_name),
        "atp_flux_glc_o2": energy_glc_o2_dict.get(model_name),
        "atp_flux_with_transporters": energy_transporter_dict.get(model_name),
        "atp_flux_all_exchanges": energy_min_ex_dict.get(model_name),
        "ATP_FVA_min": fva_min_dict.get(model_name),
        "ATP_FVA_max": fva_max_dict.get(model_name),
        "ATP_max_flux_minimal_medium": minimal_medium_flux_dict.get(model_name),
        # Add minimal medium columns
        "minimal_medium_reactions": minimal_medium_reactions,
        "minimal_medium_fluxes": minimal_medium_fluxes,
        "minimal_medium_str": minimal_medium_str,
        "biomass_feasible": growth_feasible_dict.get(model_name),
        "biomass_growth_rate": growth_rate_dict.get(model_name),
        "num_blocked_biomass_precursors": blocked_precursors.get(model_name, 0),
        "fraction_blocked_precursors": blocked_fraction.get(model_name, 0),
        "igg_feasible": igg_feasible_dict.get(model_name),
        "igg_formation_rate": igg_rate_dict.get(model_name),
        "igg_precursor_producible": igg_precursor_dict.get(model_name, True)  # True if not checked
    })


# Convert to DataFrame
comparison_df = pd.DataFrame(comparison_results)

print("\n✅ Final combined comparison table ready!")
print(comparison_df)

# %%  Save final comparison table
# Project root is two levels up from python/
PROJECT_ROOT = Path.cwd()
#ANALYSIS_MODELS_DIR = PROJECT_ROOT / "analysis" / "models" / "original"
ANALYSIS_MODELS_DIR = PROJECT_ROOT / "analysis" / "models" / "exp_constrained"
#ANALYSIS_MODELS_DIR = PROJECT_ROOT / "analysis" / "models" / "GIMME"

# Create directories if they don't exist
ANALYSIS_MODELS_DIR.mkdir(parents=True, exist_ok=True)

print(f"Results will be saved to: {ANALYSIS_MODELS_DIR}")

comparison_csv_path = ANALYSIS_MODELS_DIR / "model_comparison_final.csv"
comparison_df.to_csv(comparison_csv_path, index=False)
print(f"✅ Comparison table saved at {comparison_csv_path}")

# %% Section 12: Plotting model diagnostics
sns.set(style="whitegrid", context="talk")

# List of flux columns to plot
flux_columns = [
    "atp_flux_minimal",
    "atp_flux_all_exchanges",
    "atp_flux_glc_o2",
    "atp_flux_with_transporters",
    "ATP_max_flux_minimal_medium"  # if added to comparison_df
]

# Ensure minimal medium flux is in comparison_df
if "ATP_max_flux_minimal_medium" not in comparison_df.columns and 'minimal_medium_flux_results' in locals():
    minimal_flux_dict = {r["model_name"]: r["ATP_max_flux_minimal_medium"] for r in minimal_medium_flux_results}
    comparison_df["ATP_max_flux_minimal_medium"] = comparison_df["model_name"].map(minimal_flux_dict)

# Melt the DataFrame for easier seaborn plotting
df_melted = comparison_df.melt(
    id_vars="model_name",
    value_vars=flux_columns,
    var_name="Condition",
    value_name="ATP_flux"
)

# Replace underscores for nicer labels
df_melted["Condition"] = df_melted["Condition"].str.replace("_", " ").str.title()

 #%% Plot 1 ATP Flux Diagnostics Plot
# # -----------------------------
# # 1️⃣ ATP Flux Diagnostics Plot
# # -----------------------------
#
# # Create a mapping for shorter labels
# short_labels = {
#     "iCHO2441": "iCHO2441",
#     "iCHO3K_cho_prod_generic_unblocked": "iCHO3K",
#     "iCHOv1": "iCHOv1"
# }
#
# # Add a new column for plotting
# df_melted['model_short'] = df_melted['model_name'].map(short_labels)
#
# # Plot with shorter labels
# plt.figure(figsize=(10,6))
# sns.barplot(data=df_melted, x="model_short", y="ATP_flux", hue="Condition")
# plt.ylabel("ATP Flux")
# plt.xlabel("Model")
# plt.title("ATP Flux Diagnostics Across Models and Conditions")
# plt.xticks(rotation=0)  # rotation optional, can be 45 if needed
# plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
# plt.tight_layout()
#
# # Project root is two levels up from python/
# #FIGURES_MODELS_DIR = PROJECT_ROOT / "figures" / "models" / "original"
# FIGURES_MODELS_DIR = PROJECT_ROOT / "figures" / "models" / "GIMME"
#
# # Create directories if they don't exist
# FIGURES_MODELS_DIR.mkdir(parents=True, exist_ok=True)
#
# # Save plot
# flux_plot_path = FIGURES_MODELS_DIR / "ATP_flux_diagnostics.png"
# plt.savefig(flux_plot_path, dpi=300)
# plt.show()
# print(f"✅ ATP flux plot saved at {flux_plot_path}")

# %% Add base model + context columns

def extract_base_model(name):
    if "iCHO2441" in name:
        return "iCHO2441"
    elif "iCHO3K_cho_prod_generic_unblocked" in name:
        return "iCHO3K"
    elif "iCHOv1" in name:
        return "iCHOv1"
    else:
        return "Unknown"

def extract_context(name):
    # everything after base model
    return name.replace("iCHO2441_", "") \
               .replace("iCHO3K_cho_prod_generic_unblocked_", "") \
               .replace("iCHOv1_", "")

df_melted["base_model"] = df_melted["model_name"].apply(extract_base_model)
df_melted["context"] = df_melted["model_name"].apply(extract_context)

# %% Plot: ATP flux per model across contexts
g = sns.catplot(
    data=df_melted,
    x="context",
    y="ATP_flux",
    hue="Condition",
    col="base_model",
    kind="bar",
    col_wrap=1,
    height=5,
    aspect=2.5
)
g._legend.set_title("Condition")
g._legend.set_bbox_to_anchor((1, 1))  # push it slightly to the right

g._legend.set_bbox_to_anchor((0.95, 0.95))
# Titles and labels
g.set_titles("{col_name}")
g.set_axis_labels("Context", "ATP Flux")

# Rotate x labels
for ax in g.axes.flatten():
    ax.tick_params(axis='x', rotation=45)

# Adjust layout to make space for legend
plt.subplots_adjust(bottom=0.15)

plt.tight_layout()

# Project root is two levels up from python/
#FIGURES_MODELS_DIR = PROJECT_ROOT / "figures" / "models" / "original"
#FIGURES_MODELS_DIR = PROJECT_ROOT / "figures" / "models" / "GIMME"
FIGURES_MODELS_DIR = PROJECT_ROOT / "figures" / "models" / "exp_constrained"

# Create directories if they don't exist
FIGURES_MODELS_DIR.mkdir(parents=True, exist_ok=True)
#Save
flux_plot_path = FIGURES_MODELS_DIR / "ATP_flux_diagnostics_per_model.png"
plt.savefig(flux_plot_path, dpi=300)
plt.show()

print(f"✅ ATP flux per-model plot saved at {flux_plot_path}")

# %% Plot 2 Model Stats: Reactions, Metabolites, Genes
# # ---------------------------------
# # 2️⃣ Model Stats: Reactions, Metabolites, Genes
# # ---------------------------------
#
# # Add a new column for plotting
# comparison_df['model_short'] = comparison_df['model_name'].map(short_labels)
#
# # Number of reactions, metabolites, genes
# fig, axes = plt.subplots(1, 3, figsize=(16,5))
#
# sns.barplot(data=comparison_df, x="model_short", y="num_reactions", ax=axes[0])
# axes[0].set_title("Number of Reactions")
#
# sns.barplot(data=comparison_df, x="model_short", y="num_metabolites", ax=axes[1])
# axes[1].set_title("Number of Metabolites")
#
# sns.barplot(data=comparison_df, x="model_short", y="num_genes", ax=axes[2])
# axes[2].set_title("Number of Genes")
#
# for ax in axes:
#     ax.set_xlabel("Model")
#     ax.set_ylabel("Count")
#
# plt.tight_layout()
# plt.show()

# %% Plot 3:   Blocked Reactions Plot
# # ---------------------------------
# # 3️⃣  Blocked Reactions Plot
# # ---------------------------------
# plt.figure(figsize=(8,5))
# sns.barplot(data=blocked_df, x=blocked_df['model_name'].map(short_labels), y="num_blocked_reactions")
# plt.ylabel("Number of Blocked Reactions")
# plt.xlabel("Model")
# plt.title("Blocked Reactions per Model")
# plt.xticks(rotation=0)
# plt.tight_layout()
#
# blocked_plot_path = FIGURES_MODELS_DIR / "blocked_reactions.png"
# plt.savefig(blocked_plot_path, dpi=300)
# plt.show()
# print(f"✅ Saved blocked reactions plot at {blocked_plot_path}")

# %% Blocked reactions plot (per model across contexts)

# --- Extract base model + context (reuse if already done) ---
blocked_df["base_model"] = blocked_df["model_name"].apply(extract_base_model)
blocked_df["context"] = blocked_df["model_name"].apply(extract_context)

# --- Plot ---
g = sns.catplot(
    data=blocked_df,
    x="context",
    y="num_blocked_reactions",
    col="base_model",
    kind="bar",
    col_wrap=1,
    height=5,
    aspect=2.5
)

g.set_titles("{col_name}")
g.set_axis_labels("Context", "Number of Blocked Reactions")

# Rotate x labels
for ax in g.axes.flatten():
    ax.tick_params(axis='x', rotation=45)

plt.tight_layout()

# Save
blocked_plot_path = FIGURES_MODELS_DIR / "blocked_reactions_per_model.png"
plt.savefig(blocked_plot_path, dpi=300)

plt.show()
print(f"✅ Saved blocked reactions plot at {blocked_plot_path}")

# %% plot 4 :  FVA: ATP min/max fluxes per model
# # ---------------------------------
# # 3️⃣  FVA: ATP min/max fluxes per model
# # ---------------------------------
# plt.figure(figsize=(8,5))
#
# # Create error bars for min/max flux
# comparison_df['FVA_range'] = comparison_df['ATP_FVA_max'] - comparison_df['ATP_FVA_min']
#
# sns.barplot(
#     data=comparison_df,
#     x='model_short',
#     y='ATP_FVA_min',
#     color='skyblue',
#     label='Min ATP Flux'
# )
# sns.barplot(
#     data=comparison_df,
#     x='model_short',
#     y='FVA_range',
#     bottom=comparison_df['ATP_FVA_min'],
#     color='orange',
#     label='Max ATP Flux'
# )
#
# plt.ylabel("ATP Flux (FVA)")
# plt.xlabel("Model")
# plt.title("ATP Flux Variability Analysis (FVA) per Model")
# plt.legend()
# plt.tight_layout()
#
# # Save figure
# fva_plot_path = FIGURES_MODELS_DIR / "ATP_FVA.png"
# plt.savefig(fva_plot_path, dpi=300)
# plt.show()
# print(f"✅ Saved ATP FVA plot at {fva_plot_path}")

# %% FVA plot per model across contexts
comparison_df["base_model"] = comparison_df["model_name"].apply(extract_base_model)
comparison_df["context"] = comparison_df["model_name"].apply(extract_context)

# Compute FVA range
comparison_df["FVA_range"] = comparison_df["ATP_FVA_max"] - comparison_df["ATP_FVA_min"]

# --- Create FacetGrid ---
g = sns.FacetGrid(
    comparison_df,
    col="base_model",
    col_wrap=1,
    height=5,
    aspect=2.5,
    sharey=True
)

# --- Plot manually per facet ---
for ax, (model_name, df_sub) in zip(g.axes.flatten(), comparison_df.groupby("base_model")):
    x = range(len(df_sub))

    # Sort by context for consistency (optional but recommended)
    df_sub = df_sub.sort_values("context")

    ax.bar(
        x,
        df_sub["ATP_FVA_min"],
        label="Min ATP Flux"
    )

    ax.bar(
        x,
        df_sub["FVA_range"],
        bottom=df_sub["ATP_FVA_min"],
        label="Max ATP Flux"
    )

    ax.set_xticks(x)
    ax.set_xticklabels(df_sub["context"], rotation=45)
    ax.set_title(model_name)
    ax.set_ylabel("ATP Flux (FVA)")
    ax.set_xlabel("Context")

# --- Add legend (once) ---
handles, labels = ax.get_legend_handles_labels()
g.fig.legend(handles, labels, loc="upper right")

plt.tight_layout()

# Save
fva_plot_path = FIGURES_MODELS_DIR / "ATP_FVA_per_model.png"
plt.savefig(fva_plot_path, dpi=300)

plt.show()
print(f"✅ Saved ATP FVA plot at {fva_plot_path}")