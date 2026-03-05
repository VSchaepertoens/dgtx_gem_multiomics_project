# %% Imports
import cobra
import pandas as pd
from python.utils.model_utils import load_models, get_models_dir

# %% Load models
models, summary_df = load_models()
MODELS_DIR = get_models_dir()

# %% Initialize results containers
energy_results = []
exchange_results = []

# Define the oxygen exchange reactions for each model
oxygen_exchanges = {
    "iCHO3K_cho_prod_generic_unblocked": "EX_o2_e",
    "iCHO2441": "EX_o2(e)",
    "iCHOv1": "EX_o2_e"
}

# Define the ATP maintenance reactions for each model
atp_reactions = {
    "iCHO2441": "DM_atp[c]",
    "iCHO3K_cho_prod_generic_unblocked": "ATPM",
    "iCHOv1": "ATPM"
}

# %% Loop over models and run diagnostics
for model_name, model in models.items():
    print(f"\nRunning diagnostics for {model_name}...")

    # Basic statistics
    num_reactions = len(model.reactions)
    num_metabolites = len(model.metabolites)
    num_genes = len(model.genes)

    # Identify exchange reactions (EX_)
    ex_reactions = [r.id for r in model.reactions if r.id.startswith("EX_")]
    exchange_results.append({
        "model_name": model_name,
        "num_exchanges": len(ex_reactions),
        "exchange_ids": ex_reactions
    })

    # Energy metabolism test: Maximize ATP with only O2 and glucose available
    with model:

        # 1️⃣ Close all exchange reactions
        for rxn in model.exchanges:
            rxn.lower_bound = 0
            rxn.upper_bound = 0

        # 2️⃣ Open the correct oxygen exchange for this model
        o2_rxn_id = oxygen_exchanges.get(model_name)

        if o2_rxn_id and o2_rxn_id in model.reactions:
            o2_rxn = model.reactions.get_by_id(o2_rxn_id)
            o2_rxn.lower_bound = -1000
            print(f"{model_name}: Opened oxygen exchange -> {o2_rxn.id}")
            print(f"{model_name}: Reaction equation -> {o2_rxn.reaction}")
        else:
            print(f"{model_name}: Oxygen exchange reaction not found!")

        #  Set ATP maintenance reaction as objective
        atp_rxn_id = atp_reactions.get(model_name)

        if atp_rxn_id and atp_rxn_id in model.reactions:
            atp_rxn = model.reactions.get_by_id(atp_rxn_id)
            model.objective = atp_rxn

            print(f"{model_name}: Using ATP objective -> {atp_rxn.id}")
            print(f"{model_name}: Reaction equation -> {atp_rxn.reaction}")

            sol = model.optimize()
            atp_flux = sol.objective_value

            print(f"{model_name}: Max ATP flux = {atp_flux}")

        else:
            print(f"{model_name}: ATP maintenance reaction not found!")
            atp_flux = None

        energy_results.append({
            "model_name": model_name,
            "atp_flux_minimal": atp_flux
        })


# %%
# Define the glucose uptake reactions for each model
glucose_exchanges = {
    "iCHO2441": "EX_glc(e)",
    "iCHO3K_cho_prod_generic_unblocked": "EX_glc_e",
    "iCHOv1": "EX_glc__D_e"
}

 # Loop over models and run diagnostics
for model_name, model in models.items():
    print(f"\nRunning diagnostics for {model_name}...")

    # Basic statistics
    num_reactions = len(model.reactions)
    num_metabolites = len(model.metabolites)
    num_genes = len(model.genes)

    # Identify exchange reactions (EX_)
    ex_reactions = [r.id for r in model.reactions if r.id.startswith("EX_")]
    exchange_results.append({
        "model_name": model_name,
        "num_exchanges": len(ex_reactions),
        "exchange_ids": ex_reactions
    })

    # Energy metabolism test: Maximize ATP with only O2 and glucose available
    with model:

        # 1️⃣ Close all exchange reactions
        for rxn in model.exchanges:
            rxn.lower_bound = 0
            rxn.upper_bound = 0

        # 2️⃣ Open the correct oxygen exchange for this model
        o2_rxn_id = oxygen_exchanges.get(model_name)

        if o2_rxn_id and o2_rxn_id in model.reactions:
            o2_rxn = model.reactions.get_by_id(o2_rxn_id)
            o2_rxn.lower_bound = -1000
            print(f"{model_name}: Opened oxygen exchange -> {o2_rxn.id}")
            print(f"{model_name}: Reaction equation -> {o2_rxn.reaction}")
        else:
            print(f"{model_name}: Oxygen exchange reaction not found!")

        # 3 Add glucose uptake
        glc_rxn_id = glucose_exchanges.get(model_name)

        if glc_rxn_id in model.reactions:
            glc_rxn = model.reactions.get_by_id(glc_rxn_id)
            glc_rxn.lower_bound = -10
            print(f"{model_name}: Opened glucose exchange -> {glc_rxn.id}")
            print(f"{model_name}: Reaction equation -> {glc_rxn.reaction}")
        else:
            print(f"{model_name}: Glucose exchange not found!")

        #  Set ATP maintenance reaction as objective
        atp_rxn_id = atp_reactions.get(model_name)

        if atp_rxn_id and atp_rxn_id in model.reactions:
            atp_rxn = model.reactions.get_by_id(atp_rxn_id)
            model.objective = atp_rxn

            print(f"{model_name}: Using ATP objective -> {atp_rxn.id}")
            print(f"{model_name}: Reaction equation -> {atp_rxn.reaction}")

            sol = model.optimize()
            atp_flux = sol.objective_value

            print(f"{model_name}: Max ATP flux = {atp_flux}")

        else:
            print(f"{model_name}: ATP maintenance reaction not found!")
            atp_flux = None

        energy_results.append({
            "model_name": model_name,
            "atp_flux_minimal": atp_flux
        })

# %% Section: Identify glucose and oxygen transporters for iCHO2441 and iCHOv1
transporter_results = {}

for model_name in ["iCHO2441", "iCHOv1"]:
    model = models[model_name]
    print(f"\nSearching for transporters in {model_name}...")

    glucose_transporters = []
    oxygen_transporters = []

    # Loop over reactions to find likely glucose and oxygen transporters
    for rxn in model.reactions:
        # Check for glucose transporters: reactions that contain "glc" in the id and are reversible/transport
        if "glc" in rxn.id.lower() and "<=>" in str(rxn.reaction):
            glucose_transporters.append(rxn.id)

        # Check for oxygen transporters: reactions that contain "o2" in the id and are reversible/transport
        if "o2" in rxn.id.lower() and "<=>" in str(rxn.reaction):
            oxygen_transporters.append(rxn.id)

    transporter_results[model_name] = {
        "glucose": glucose_transporters,
        "oxygen": oxygen_transporters
    }

    print(f"{model_name}: Found {len(glucose_transporters)} glucose transporter(s): {glucose_transporters}")
    print(f"{model_name}: Found {len(oxygen_transporters)} oxygen transporter(s): {oxygen_transporters}")

# %% Section: New diagnostic loop for iCHO2441 and iCHOv1
energy_results_transporters = []

for model_name in ["iCHO2441", "iCHOv1"]:
    model = models[model_name]
    print(f"\nRunning ATP diagnostics with transporters for {model_name}...")

    glucose_tx = transporter_results[model_name]["glucose"]
    oxygen_tx = transporter_results[model_name]["oxygen"]

    with model:
        # 1️⃣ Close all exchange reactions initially
        for rxn in model.exchanges:
            rxn.lower_bound = 0
            rxn.upper_bound = 0

        # 2️⃣ Open glucose and oxygen transporters
        for rxn_id in glucose_tx + oxygen_tx:
            try:
                rxn = model.reactions.get_by_id(rxn_id)
                rxn.lower_bound = -1000
                rxn.upper_bound = 1000
                print(f"{model_name}: Opened transporter {rxn.id}")
            except KeyError:
                print(f"{model_name}: Transporter {rxn_id} not found in model!")

        # 3️⃣ Maximize ATP production
        # Use ATP reaction mapping you already defined
        atp_id = atp_reactions.get(model_name, "ATPM")
        try:
            atp_rxn = model.reactions.get_by_id(atp_id)
            model.objective = atp_rxn
            sol = model.optimize()
            atp_flux = sol.objective_value
        except KeyError:
            print(f"{model_name}: ATP reaction {atp_id} not found!")
            atp_flux = None

        energy_results_transporters.append({
            "model_name": model_name,
            "atp_flux_with_transporters": atp_flux
        })

        print(f"{model_name}: Max ATP flux = {atp_flux}")