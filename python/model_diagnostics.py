# %% Imports
import cobra
import pandas as pd
from python.utils.model_utils import load_models, get_models_dir

from cobra.flux_analysis import find_blocked_reactions
from cobra.flux_analysis import flux_variability_analysis
from cobra.medium import minimal_medium

import matplotlib.pyplot as plt
import seaborn as sns


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

# # %% Section: Compile final comparison table
# comparison_results = []
#
# # Create lookup dicts for easy merging
# energy_min_dict = {r["model_name"]: r["atp_flux_minimal"] for r in energy_results_minimal}
# energy_min_ex_dict = {r["model_name"]: r["atp_flux_all_exchanges"] for r in energy_results_minimal_exchanges}
# energy_glc_o2_dict = {r["model_name"]: r["atp_flux_glc_o2"] for r in energy_results_transporters if "atp_flux_glc_o2" in r}
# energy_transporter_dict = {r["model_name"]: r["atp_flux_with_transporters"] for r in energy_results_transporters if "atp_flux_with_transporters" in r}
#
# for model_name, model in models.items():
#     num_reactions = len(model.reactions)
#     num_metabolites = len(model.metabolites)
#     num_genes = len(model.genes)
#     ex_entry = next((r for r in exchange_results if r["model_name"] == model_name), {})
#     num_exchanges = ex_entry.get("num_exchanges", 0)
#
#     # Transporters
#     transporters = transporter_results.get(model_name, {})
#     glucose_transporters = transporters.get("glucose", [])
#     oxygen_transporters = transporters.get("oxygen", [])
#
#     comparison_results.append({
#         "model_name": model_name,
#         "num_reactions": num_reactions,
#         "num_metabolites": num_metabolites,
#         "num_genes": num_genes,
#         "num_exchanges": num_exchanges,
#         "atp_flux_minimal": energy_min_dict.get(model_name),
#         "atp_flux_all_exchanges": energy_min_ex_dict.get(model_name),
#         "atp_flux_glc_o2": energy_glc_o2_dict.get(model_name),
#         "atp_flux_with_transporters": energy_transporter_dict.get(model_name)
#     })
#
# # Convert to DataFrame
# comparison_df = pd.DataFrame(comparison_results)
# #comparison_df.to_csv(f"{MODELS_DIR}/model_comparison_final.csv", index=False)
# print("\n✅ Final comparison table ready!")
# print(comparison_df)

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

# %% Section 9: Compile final comparison table
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
fva_min_dict = {r["model_name"]: r["ATP_min_flux"] for r in fva_results}
fva_max_dict = {r["model_name"]: r["ATP_max_flux"] for r in fva_results}
minimal_medium_flux_dict = {r["model_name"]: r["ATP_max_flux_minimal_medium"] for r in minimal_medium_flux_results}

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
        "atp_flux_minimal": energy_min_dict.get(model_name),
        "atp_flux_glc_o2": energy_glc_o2_dict.get(model_name),
        "atp_flux_with_transporters": energy_transporter_dict.get(model_name),
        "atp_flux_all_exchanges": energy_min_ex_dict.get(model_name),
        "ATP_FVA_min": fva_min_dict.get(model_name),
        "ATP_FVA_max": fva_max_dict.get(model_name),
        "ATP_max_flux_minimal_medium": minimal_medium_flux_dict.get(model_name)
    })

# Convert to DataFrame
comparison_df = pd.DataFrame(comparison_results)

# Save if needed
# comparison_df.to_csv(f"{MODELS_DIR}/model_comparison_final.csv", index=False)

print("\n✅ Final combined comparison table ready!")
print(comparison_df)

# %% Section 10: Plotting model diagnostics
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

# %% Plot: ATP fluxes across conditions

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
plt.show()

# %% Optional: Other diagnostics plots

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