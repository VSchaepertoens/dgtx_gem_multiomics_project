# model_utils.py
import cobra
import pandas as pd
from pathlib import Path
from typing import Optional, Tuple, Dict
import re

def get_models_dir(mode: str = "original", custom_path: Optional[Path] = None) -> Path:
    """
    Return the path to the models folder.

    Parameters:
        mode (str): Either "original", "GIMME", or "mCADRE" (future). Default is "original".
        custom_path (Path, optional): If provided, this path is returned directly.

    Returns:
        Path: The folder containing the SBML models.
    """
    if custom_path is not None:
        return Path(custom_path).resolve()

    # Project root is two levels up from utils/
    BASE_DIR = Path(__file__).resolve().parent.parent.parent

    if mode == "original":
        MODELS_DIR = BASE_DIR / "models" / "original"
    elif mode == "GIMME":
        MODELS_DIR = BASE_DIR / "models" / "contextualized" / "GIMME"
    elif mode == "mCADRE":
        MODELS_DIR = BASE_DIR / "models" / "contextualized" / "mCADRE"
    elif mode == "exp_constrained":
        MODELS_DIR = BASE_DIR / "models" / "exp_constrained"
    elif mode == "GIMME25perc":
        MODELS_DIR = BASE_DIR / "models" / "contextualized" / "GIMME_25perc"
    elif mode == "GIMME_tpm_25perc":
        MODELS_DIR = BASE_DIR / "models" / "contextualized" / "GIMME_tpm_25perc"
    else:
        raise ValueError(f"Unknown mode '{mode}'. Choose 'original', 'GIMME', or 'mCADRE', or 'exp_constrained'.")

    return MODELS_DIR

def load_models(models_folder: Optional[Path] = None) -> Tuple[Dict[str, cobra.Model], pd.DataFrame]:
    """
    Load all SBML models from a given folder, create a summary dataframe,
    and return both the models dictionary and the summary.

    Parameters:
        models_folder (Path, optional): Path to folder containing SBML models.
                                        If None, defaults to get_models_dir(mode="original").

    Returns:
        models (dict): model_name -> cobra.Model
        summary_df (pd.DataFrame): Summary information about each model
    """
    if models_folder is None:
        models_folder = get_models_dir()

    models_folder = Path(models_folder).resolve()
    print(f"Loading models from folder: {models_folder}")

    models = {}
    summary_data = []

    for model_file in models_folder.glob("*.xml"):
        model_name = model_file.stem
        print(f"\nLoading {model_name}...")
        try:
            model = cobra.io.read_sbml_model(model_file)
            models[model_name] = model

            # Objective expression as string
            obj_expr = str(model.objective.expression)
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

    summary_df = pd.DataFrame(summary_data)
    return models, summary_df


# %% Section: Model naming utilities

def get_common_model_name(model_name: str) -> str:
    """
    Strips suffixes like '_Constant_X', '_Temp. shifted', '_GIMME'
    to return the original/common model name used in mappings.

    Examples
    --------
    iCHO2441_Constant_1      -> iCHO2441
    iCHO2441_GIMME            -> iCHO2441
    iCHO2441_Temp. shifted_2 -> iCHO2441
    """
    for suffix in ["_Constant", "_Temp. shifted", "_GIMME"]:
        if suffix in model_name:
            return model_name.split(suffix)[0]
    return model_name

# %% Section: Parse condition and window from model name
def parse_condition_window(model_name: str):
    """
    Extract Condition and Window from model name.

    Works for:
    - iCHOv1_Constant_4_constrained
    - iCHOv1_Temp. shifted_2_constrained
    - iCHO3K_..._Constant_1_GIMME
    - future variants with additional suffixes

    Returns
    -------
    (condition: str, window: int)
    """

    # 🔹 Detect condition
    if "Constant" in model_name:
        condition = "Constant"
    elif "Temp. shifted" in model_name:
        condition = "Temp. shifted"
    else:
        print(f"⚠️ Could not detect condition in {model_name}")
        return None, None

    # 🔹 Extract window (last integer in string)
    numbers = re.findall(r'\d+', model_name)

    if not numbers:
        print(f"⚠️ Could not extract window from {model_name}")
        return condition, None

    window = int(numbers[-1])

    return condition, window