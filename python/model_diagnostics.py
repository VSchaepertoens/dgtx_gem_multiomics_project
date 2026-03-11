# %% Imports
import cobra
import pandas as pd
from python.utils.model_utils import load_models, get_models_dir

from cobra.flux_analysis import find_blocked_reactions
from cobra.flux_analysis import flux_variability_analysis
from cobra.medium import minimal_medium

import matplotlib.pyplot as plt
import seaborn as sns

# %% Hard-coded save directories
from pathlib import Path

# %% Section 1:  Load models
models, summary_df = load_models()
MODELS_DIR = get_models_dir()

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

# %% Section 4.1: Minimal ATP flux and exchange reactions
for model_name, model in models.items():
    print(f"\nRunning minimal ATP diagnostics for {model_name}...")

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
            rxn.lower_bound = 0
            rxn.upper_bound = 0

        atp_rxn_id = atp_reactions.get(model_name)
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

    with model:
        # Close all exchange reactions
        for rxn in model.exchanges:
            rxn.lower_bound = 0
            rxn.upper_bound = 0

        # Open oxygen exchange
        o2_rxn_id = oxygen_exchanges.get(model_name)
        if o2_rxn_id in model.reactions:
            rxn = model.reactions.get_by_id(o2_rxn_id)
            rxn.lower_bound = -1000
            rxn.upper_bound = 1000
            print(f"{model_name}: Opened oxygen exchange -> {rxn.id}")

        # Open glucose exchange
        glc_rxn_id = glucose_exchanges.get(model_name)
        if glc_rxn_id in model.reactions:
            rxn = model.reactions.get_by_id(glc_rxn_id)
            rxn.lower_bound = -10
            rxn.upper_bound = 1000
            print(f"{model_name}: Opened glucose exchange -> {rxn.id}")

        # Maximize ATP
        atp_rxn_id = atp_reactions.get(model_name)
        if atp_rxn_id and atp_rxn_id in model.reactions:
            model.objective = model.reactions.get_by_id(atp_rxn_id)
            sol = model.optimize()
            atp_flux = sol.objective_value
            print(f"{model_name}: Max ATP flux (O2 + glucose) = {atp_flux}")
        else:
            atp_flux = None

        energy_results_transporters.append({
            "model_name": model_name,
            "atp_flux_glc_o2": atp_flux
        })

# %% Section 4.3: ATP flux with transporters (all models)
for model_name in ["iCHO2441", "iCHOv1", "iCHO3K_cho_prod_generic_unblocked"]:
    model = models[model_name]

    glucose_tx = glucose_oxygen_transporters[model_name]["glucose"]
    oxygen_tx = glucose_oxygen_transporters[model_name]["oxygen"]

    if not glucose_tx and not oxygen_tx:
        print(f"{model_name}: No transporters found, skipping flux with transporters")
        continue

    with model:  # Use a context manager so changes are temporary
        # Close all exchange reactions first
        for rxn in model.exchanges:
            rxn.lower_bound = 0
            rxn.upper_bound = 0

        print(f"\n{model_name}: All exchange reactions closed.")

        # Open exchange reactions for glucose and oxygen
        glc_ex_id = glucose_exchanges.get(model_name)
        o2_ex_id = oxygen_exchanges.get(model_name)

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
        atp_rxn_id = atp_reactions.get(model_name)
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
            rxn.lower_bound = -1000
            rxn.upper_bound = 1000

        atp_rxn_id = atp_reactions.get(model_name)
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
    atp_rxn = atp_reactions[model_name]

    if atp_rxn in blocked_df.loc[
        blocked_df["model_name"] == model_name, "blocked_reactions"
    ].values[0]:
        print(f"{model_name}: ATP reaction is BLOCKED")
    else:
        print(f"{model_name}: ATP reaction is NOT blocked")

# %% Section 6: Flux variability analysis for ATP
fva_results = []

for model_name, model in models.items():

    atp_rxn_id = atp_reactions.get(model_name)

    if atp_rxn_id not in model.reactions:
        print(f"{model_name}: ATP reaction not found")
        continue

    print(f"\nRunning FVA for ATP in {model_name}...")

    fva = flux_variability_analysis(
        model,
        reaction_list=[atp_rxn_id],
        fraction_of_optimum=0
    )

    min_flux = fva.loc[atp_rxn_id, "minimum"]
    max_flux = fva.loc[atp_rxn_id, "maximum"]

    print(f"{model_name}: ATP flux range = [{min_flux}, {max_flux}]")

    fva_results.append({
        "model_name": model_name,
        "ATP_min_flux": min_flux,
        "ATP_max_flux": max_flux
    })

fva_df = pd.DataFrame(fva_results)

print("\nATP flux variability:")
print(fva_df)


 # %% Section 7: Minimal medium required for ATP production
minimal_medium_results = []

for model_name, model in models.items():

    print(f"\nComputing minimal medium for ATP production in {model_name}...")

    with model:

        # Close all exchanges first
        for rxn in model.exchanges:
            rxn.lower_bound = -1000
            rxn.upper_bound = 1000

        # Set ATP reaction as objective
        atp_rxn_id = atp_reactions.get(model_name)

        if atp_rxn_id not in model.reactions:
            print(f"{model_name}: ATP reaction not found")
            continue

        model.objective = model.reactions.get_by_id(atp_rxn_id)

        # Calculate minimal medium
        medium = minimal_medium(model, 0.1)

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

# %% Section 8: Max ATP flux with minimal medium open
minimal_medium_flux_results = []

for model_name, model in models.items():
    print(f"\nComputing ATP flux with minimal medium for {model_name}...")

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
            rxn.lower_bound = 0
            rxn.upper_bound = 0

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
        for rxn_id in medium.index.tolist() + [glucose_exchanges[model_name], oxygen_exchanges[model_name]]:
            if rxn_id in model.reactions:
                rxn = model.reactions.get_by_id(rxn_id)
                rxn.lower_bound = -1000
                rxn.upper_bound = 1000
                print(f"{model_name}: Opened minimal medium + glucose/O2 exchange -> {rxn.id}")
            else:
                print(f"{model_name}: Reaction {rxn_id} NOT FOUND")

        # Maximize ATP
        atp_rxn_id = atp_reactions.get(model_name)
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


# %% Section 9: Biomass feasibility check
growth_results = []

for model_name, model in models.items():

    biomass_rxn_id = biomass_reactions[model_name]
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
        biomass_rxn_id = biomass_rxns[model_name]
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

# %% Section 11: Compile final comparison table
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
blocked_dict = {r["model_name"]: r["num_blocked_reactions"] for r in blocked_results}
blocked_rxn_ids = {r["model_name"]: r["blocked_reactions"] for r in blocked_results}
fva_min_dict = {r["model_name"]: r["ATP_min_flux"] for r in fva_results}
fva_max_dict = {r["model_name"]: r["ATP_max_flux"] for r in fva_results}
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

growth_feasible_dict = {r["model"]: r["feasible"] for r in growth_results}
growth_rate_dict = {r["model"]: r["growth_rate"] for r in growth_results}

# Number of blocked precursors per model
blocked_precursors = (
    precursor_df[precursor_df["producible"] == False]
    .groupby("model")
    .size()
    .to_dict()
)

# Fraction of blocked precursors (optional)
total_precursors = precursor_df.groupby("model").size().to_dict()
blocked_fraction = {m: blocked_precursors.get(m, 0)/total_precursors[m] for m in models.keys()}


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
        "num_blocked_reactions": blocked_dict.get(model_name),
        "blocked_reaction_ids": blocked_rxn_ids.get(model_name, []),
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
        "fraction_blocked_precursors": blocked_fraction.get(model_name, 0)
    })

# Convert to DataFrame
comparison_df = pd.DataFrame(comparison_results)

print("\n✅ Final combined comparison table ready!")
print(comparison_df)

# %%  Save final comparison table
# Project root is two levels up from python/
PROJECT_ROOT = Path.cwd()
ANALYSIS_MODELS_DIR = PROJECT_ROOT / "analysis" / "models" / "original"

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
# -----------------------------
# 1️⃣ ATP Flux Diagnostics Plot
# -----------------------------

# Create a mapping for shorter labels
short_labels = {
    "iCHO2441": "iCHO2441",
    "iCHO3K_cho_prod_generic_unblocked": "iCHO3K",
    "iCHOv1": "iCHOv1"
}

# Add a new column for plotting
df_melted['model_short'] = df_melted['model_name'].map(short_labels)

# Plot with shorter labels
plt.figure(figsize=(10,6))
sns.barplot(data=df_melted, x="model_short", y="ATP_flux", hue="Condition")
plt.ylabel("ATP Flux")
plt.xlabel("Model")
plt.title("ATP Flux Diagnostics Across Models and Conditions")
plt.xticks(rotation=0)  # rotation optional, can be 45 if needed
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()

# Project root is two levels up from python/
FIGURES_MODELS_DIR = PROJECT_ROOT / "figures" / "models" / "original"

# Create directories if they don't exist
FIGURES_MODELS_DIR.mkdir(parents=True, exist_ok=True)

# Save plot
flux_plot_path = FIGURES_MODELS_DIR / "ATP_flux_diagnostics.png"
plt.savefig(flux_plot_path, dpi=300)
plt.show()
print(f"✅ ATP flux plot saved at {flux_plot_path}")

# %% Plot 2 Model Stats: Reactions, Metabolites, Genes
# ---------------------------------
# 2️⃣ Model Stats: Reactions, Metabolites, Genes
# ---------------------------------

# Add a new column for plotting
comparison_df['model_short'] = comparison_df['model_name'].map(short_labels)

# Number of reactions, metabolites, genes
fig, axes = plt.subplots(1, 3, figsize=(16,5))

sns.barplot(data=comparison_df, x="model_short", y="num_reactions", ax=axes[0])
axes[0].set_title("Number of Reactions")

sns.barplot(data=comparison_df, x="model_short", y="num_metabolites", ax=axes[1])
axes[1].set_title("Number of Metabolites")

sns.barplot(data=comparison_df, x="model_short", y="num_genes", ax=axes[2])
axes[2].set_title("Number of Genes")

for ax in axes:
    ax.set_xlabel("Model")
    ax.set_ylabel("Count")

plt.tight_layout()
plt.show()

# %% Plot 3:   Blocked Reactions Plot
# ---------------------------------
# 3️⃣  Blocked Reactions Plot
# ---------------------------------
plt.figure(figsize=(8,5))
sns.barplot(data=blocked_df, x=blocked_df['model_name'].map(short_labels), y="num_blocked_reactions")
plt.ylabel("Number of Blocked Reactions")
plt.xlabel("Model")
plt.title("Blocked Reactions per Model")
plt.xticks(rotation=0)
plt.tight_layout()

blocked_plot_path = FIGURES_MODELS_DIR / "blocked_reactions.png"
plt.savefig(blocked_plot_path, dpi=300)
plt.show()
print(f"✅ Saved blocked reactions plot at {blocked_plot_path}")

# %% plot 4 :  FVA: ATP min/max fluxes per model
# ---------------------------------
# 3️⃣  FVA: ATP min/max fluxes per model
# ---------------------------------
plt.figure(figsize=(8,5))

# Create error bars for min/max flux
comparison_df['FVA_range'] = comparison_df['ATP_FVA_max'] - comparison_df['ATP_FVA_min']

sns.barplot(
    data=comparison_df,
    x='model_short',
    y='ATP_FVA_min',
    color='skyblue',
    label='Min ATP Flux'
)
sns.barplot(
    data=comparison_df,
    x='model_short',
    y='FVA_range',
    bottom=comparison_df['ATP_FVA_min'],
    color='orange',
    label='Max ATP Flux'
)

plt.ylabel("ATP Flux (FVA)")
plt.xlabel("Model")
plt.title("ATP Flux Variability Analysis (FVA) per Model")
plt.legend()
plt.tight_layout()

# Save figure
fva_plot_path = FIGURES_MODELS_DIR / "ATP_FVA.png"
plt.savefig(fva_plot_path, dpi=300)
plt.show()
print(f"✅ Saved ATP FVA plot at {fva_plot_path}")
