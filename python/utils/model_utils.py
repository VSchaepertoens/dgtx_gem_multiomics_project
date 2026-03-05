# model_utils.py
from pathlib import Path
import cobra
import pandas as pd

from pathlib import Path

def get_models_dir():
    """
    Return the path to the models/original folder relative to the project root.
    Assumes this utility is in python/utils/
    """
    # Project root is two levels up from utils/
    BASE_DIR = Path(__file__).resolve().parent.parent.parent
    MODELS_DIR = BASE_DIR / "models" / "original"
    return MODELS_DIR

def load_models():
    """
    Load all SBML models from models/original, create a summary dataframe,
    and return both the models dictionary, the summary, and the models folder path.

    Returns:
        models (dict): model_name -> cobra.Model
        summary_df (pd.DataFrame): summary information about each model
    """
    # Project root is two levels up from utils/
    BASE_DIR = Path(__file__).resolve().parent.parent.parent
    MODELS_DIR = BASE_DIR / "models" / "original"

    print(f"Project root: {BASE_DIR}")
    print(f"Models directory: {MODELS_DIR}")

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

            # Try to extract first reaction ID from expression string
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

    # Create summary dataframe
    summary_df = pd.DataFrame(summary_data)

    return models, summary_df