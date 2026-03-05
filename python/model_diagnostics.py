# %% Imports
import cobra
import pandas as pd
from python.utils.model_utils import load_models, get_models_dir

# %% Load models
models, summary_df = load_models()
MODELS_DIR = get_models_dir()

# %% Initialize results containers
energy_results_minimal = []
energy_results_transporters = []
exchange_results = []
transporter_results = {}

# %% Define key reactions for each model
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

# %% Loop 1: Minimal ATP flux and exchange reactions
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

# %% Loop 2: ATP flux with O2 and glucose (including glucose uptake)
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

# %% Loop 3: ATP flux with transporters (all models)
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
# %% Section: Compile final comparison table
comparison_results = []

# Create lookup dicts for easy merging
energy_min_dict = {r["model_name"]: r["atp_flux_minimal"] for r in energy_results_minimal}
energy_glc_o2_dict = {r["model_name"]: r["atp_flux_glc_o2"] for r in energy_results_transporters if "atp_flux_glc_o2" in r}
energy_transporter_dict = {r["model_name"]: r["atp_flux_with_transporters"] for r in energy_results_transporters if "atp_flux_with_transporters" in r}

for model_name, model in models.items():
    num_reactions = len(model.reactions)
    num_metabolites = len(model.metabolites)
    num_genes = len(model.genes)
    ex_entry = next((r for r in exchange_results if r["model_name"] == model_name), {})
    num_exchanges = ex_entry.get("num_exchanges", 0)

    # Transporters
    transporters = transporter_results.get(model_name, {})
    glucose_transporters = transporters.get("glucose", [])
    oxygen_transporters = transporters.get("oxygen", [])

    comparison_results.append({
        "model_name": model_name,
        "num_reactions": num_reactions,
        "num_metabolites": num_metabolites,
        "num_genes": num_genes,
        "num_exchanges": num_exchanges,
        "atp_flux_minimal": energy_min_dict.get(model_name),
        "atp_flux_glc_o2": energy_glc_o2_dict.get(model_name),
        "atp_flux_with_transporters": energy_transporter_dict.get(model_name)
    })

# Convert to DataFrame
comparison_df = pd.DataFrame(comparison_results)
#comparison_df.to_csv(f"{MODELS_DIR}/model_comparison_final.csv", index=False)
print("\n✅ Final comparison table ready!")
print(comparison_df)


