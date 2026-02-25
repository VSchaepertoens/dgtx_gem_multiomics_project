# %% Imports
import os
from pathlib import Path
from datetime import datetime
import tempfile
import os
import subprocess
import sys

import pandas as pd
import cobra
import memote

# %% Define project root automatically (relative to script location)
BASE_DIR = Path().resolve()
MODELS_DIR = BASE_DIR / "models" / "original"

print(f"Project root: {BASE_DIR}")
print(f"Models directory: {MODELS_DIR}")

# %% Load all SBML models in models/original
models = {}
summary_data = []

for model_file in MODELS_DIR.glob("*.xml"):
    model_name = model_file.stem

    print(f"\nLoading {model_name}...")

    try:
        model = cobra.io.read_sbml_model(model_file)
        models[model_name] = model

        # Objective expression as string
        obj_expr = str(model.objective.expression)

        # Optional: try to extract first reaction ID from expression string
        # This assumes reaction IDs appear in the SBML IDs
        obj_rxn = obj_expr.split('*')[1].split()[0] if '*' in obj_expr else obj_expr

        summary_data.append({
            "model_name": model_name,
            "reactions": len(model.reactions),
            "metabolites": len(model.metabolites),
            "genes": len(model.genes),
            "compartments": len(model.compartments),
            "objective": obj_expr,
            "objective_rxn": obj_rxn,
        })

        print(f"Loaded {model_name} successfully!")

    except Exception as e:
        print(f"Error loading {model_name}: {e}")

# %% Create summary dataframe
summary_df = pd.DataFrame(summary_data)
# Print objective column only
print(summary_df[["model_name", "objective_rxn"]])

# %% Check biomass metabolites for each model
for model_name, model in models.items():
    # Get the objective reaction expression as string
    obj_expr = str(model.objective.expression)

    # Heuristic: get first reaction ID from expression
    if '*' in obj_expr:
        # Example: '1.0*biomass_cho - 1.0*biomass_cho_reverse'
        rxn_id = obj_expr.split('*')[1].split()[0]
    else:
        rxn_id = obj_expr  # fallback if expression is just a single ID

    # Safely get reaction and count metabolites
    try:
        biomass_rxn = model.reactions.get_by_id(rxn_id)
        n_mets = len(biomass_rxn.metabolites)
        print(f"{model_name}: Biomass '{rxn_id}' has {n_mets} metabolites")
    except KeyError:
        print(f"{model_name}: Reaction '{rxn_id}' not found in model")

# %% Should add the names of all exchange reactions across the models

# # Optional: Save summary
# output_path = BASE_DIR / "reports"
# output_path.mkdir(exist_ok=True)
#
# summary_df.to_csv(output_path / "phase1_model_structure_summary.csv", index=False)


#"C:\Users\b1095820\Anaconda3\envs\env\Scripts\memote.exe" report snapshot "C:\Users\b1095820\Documents\DGTX\dgtx_gem_multiomics_project\models\original\iCHOv1.xml" --filename "C:\Users\b1095820\Documents\DGTX\dgtx_gem_multiomics_project\models\memote_reports\iCHOv1_memote.html" --report-json "C:\Users\b1095820\Documents\DGTX\dgtx_gem_multiomics_project\models\memote_reports\iCHOv1_memote.json"