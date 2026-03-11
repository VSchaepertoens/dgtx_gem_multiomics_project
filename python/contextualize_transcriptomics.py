# %% Imports
import cobra
import pandas as pd
from python.utils.model_utils import load_models, get_models_dir

import re

# %% Hard-coded save directories
from pathlib import Path

# %% Section 1:  Load models & set-up paths
models, summary_df = load_models()
MODELS_DIR = get_models_dir()

# %% Project root is two levels up from python/
PROJECT_ROOT = Path.cwd()
# ANALYSIS_MODELS_DIR = PROJECT_ROOT / "analysis" / "models" / "original"

# Create directories if they don't exist
# ANALYSIS_MODELS_DIR.mkdir(parents=True, exist_ok=True)

# %% Section 2: Inspect gene identifiers

for model_name, model in models.items():

    print(f"\n===== {model_name} =====")
    print(f"Number of genes: {len(model.genes)}")

    # show first few gene IDs
    print("\nFirst 10 gene IDs:")
    for g in list(model.genes)[:10]:
        print(g.id)

    # check annotations
    annotated_genes = [g for g in model.genes if g.annotation]
    print(f"\nGenes with annotations: {len(annotated_genes)}")

# %% Section 3. Check gene overlaps
    expression_csv_path = PROJECT_ROOT / "data/transcriptomics/charrun_FBA_RNA3_all_V02_20260311_VS.csv"
    expr = pd.read_csv(expression_csv_path)

    for model_name, model in models.items():
        model_genes = {g.id for g in model.genes}

        expr_genes = set(
            expr.loc[expr["model"] == model_name, "EntrezID"].astype(str)
        )

        overlap = model_genes & expr_genes

        print(f"\n{model_name}")
        print("Model genes:", len(model_genes))
        print("RNA genes:", len(expr_genes))
        print("Overlap:", len(overlap))

    if annotated_genes:
        print("\nExample gene annotations:")
        for g in annotated_genes[:5]:
            print(g.id, g.annotation)

# %% Section 4.  Contextualization - Threshold-based Boolean approach

# --- Parameters ---
THRESHOLD = 0.0  # gene expression threshold to keep reaction active

# --- Functions ---
def parse_gpr(rule_str):
    """
    Convert a GPR string into nested structure for Boolean evaluation.
    Returns: single gene string, ('and', [...]), ('or', [...])
    """
    import re
    rule_str = rule_str.strip()
    if not rule_str:
        return None
    # No parentheses: split by OR / AND
    if "(" not in rule_str:
        or_parts = re.split(r"\s+OR\s+", rule_str, flags=re.IGNORECASE)
        if len(or_parts) > 1:
            return ("or", [parse_gpr(p) for p in or_parts])
        and_parts = re.split(r"\s+AND\s+", rule_str, flags=re.IGNORECASE)
        if len(and_parts) > 1:
            return ("and", [parse_gpr(p) for p in and_parts])
        return rule_str.strip()
    # Handle parentheses with placeholders
    while "(" in rule_str:
        m = re.search(r"\(([^()]+)\)", rule_str)
        if not m:
            break
        inner = m.group(1)
        parsed_inner = parse_gpr(inner)
        placeholder = f"__PLACEHOLDER_{hash(inner)}__"
        rule_str = rule_str[:m.start()] + placeholder + rule_str[m.end():]
        parse_gpr.placeholders[placeholder] = parsed_inner
    # After placeholders, split by OR / AND
    or_parts = re.split(r"\s+OR\s+", rule_str, flags=re.IGNORECASE)
    if len(or_parts) > 1:
        return ("or", [parse_gpr.part_or_placeholder(p) for p in or_parts])
    and_parts = re.split(r"\s+AND\s+", rule_str, flags=re.IGNORECASE)
    if len(and_parts) > 1:
        return ("and", [parse_gpr.part_or_placeholder(p) for p in and_parts])
    return parse_gpr.part_or_placeholder(rule_str.strip())

parse_gpr.placeholders = {}
parse_gpr.part_or_placeholder = lambda token: parse_gpr.placeholders[token] if token.startswith("__PLACEHOLDER_") else token


def eval_gpr_boolean(rule, gene_expr_dict, context, threshold=THRESHOLD):
    """
    Recursively evaluate a parsed GPR rule as Boolean:
    - Single gene: True if expression > threshold
    - AND: True if all children True
    - OR: True if any child True
    """
    if rule is None:
        return False
    if isinstance(rule, str):
        return gene_expr_dict.get(rule, {}).get(context, 0.0) > threshold
    if isinstance(rule, tuple):
        op, args = rule[0], rule[1]
        if op.lower() == "and":
            return all(eval_gpr_boolean(a, gene_expr_dict, context, threshold) for a in args)
        elif op.lower() == "or":
            return any(eval_gpr_boolean(a, gene_expr_dict, context, threshold) for a in args)
    return False

# --- Main loop: build contextualized GEMs ---

summary_list = []

for model_name, model in models.items():
    print(f"\nProcessing model: {model_name}")

    # Filter RNA for genes present in the model
    model_genes = {g.id for g in model.genes}
    expr_model = expr[expr["model"] == model_name].copy()
    expr_model["EntrezID"] = expr_model["EntrezID"].astype(str)
    expr_model = expr_model[expr_model["EntrezID"].isin(model_genes)]
    contexts = [c for c in expr_model.columns if c not in ["model", "EntrezID"]]
    gene_expr_dict = expr_model.set_index("EntrezID")[contexts].to_dict(orient="index")

    missing_genes = model_genes - set(expr_model["EntrezID"])
    if missing_genes:
        print(f"WARNING: {len(missing_genes)} model genes have no RNA expression data")

    # Evaluate reactions per context-window
    for context in contexts:
        total_reactions = len(model.reactions)
        active_count = 0

        for rxn in model.reactions:
            if not rxn.genes:
                active_count += 1  # keep reactions with no genes by default
                continue

            gpr_rule = rxn.gene_reaction_rule
            try:
                parsed_rule = parse_gpr(gpr_rule)
            except Exception as e:
                parsed_rule = None

            is_active = eval_gpr_boolean(parsed_rule, gene_expr_dict, context, threshold=THRESHOLD)
            if is_active:
                active_count += 1

        removed_count = total_reactions - active_count

        summary_list.append({
            "model": model_name,
            "context": context,
            "total_reactions": total_reactions,
            "active_reactions": active_count,
            "removed_reactions": removed_count
        })

# Convert to DataFrame
context_summary_df = pd.DataFrame(summary_list)
print(context_summary_df)