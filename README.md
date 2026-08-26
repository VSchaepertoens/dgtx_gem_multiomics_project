# Multi-Model CHO GEM Contextualization & Validation using data from DGTX NISTCHO Characterization Runs
## 📄 Overview

This repository contains the computational workflow for **multi-model contextualization, simulation, and validation of CHO genome-scale metabolic models (GEMs)** using experimental data from DGTX NISTCHO characterization runs. Three CHO GEMs (iCHOv1/iCHO1766, iCHO2441, and iCHO3K) are systematically compared and evaluated, first in their original form and subsequently under experimentally constrained and transcriptomics-contextualized conditions.

The workflow integrates experimentally measured uptake, secretion, growth, and transcriptomics data across constant and temperature-shifted conditions and multiple time windows. It includes model diagnostics, experimental constraint integration, GIMME-based contextualization, FBA/pFBA flux simulations, leave-one-exchange-out cross-validation, and quantitative validation against experimental measurements. The final analyses characterize global and reaction-level metabolic changes and generate publication-ready visualizations of predicted flux distributions and pathway dynamics.

The repository is organized as a sequential analysis pipeline, with numbered Python, R, and MATLAB scripts corresponding to the major processing and analysis steps. Intermediate models, numerical results, and figures are stored separately in `models/`, `analysis/`, and `figures/`, respectively.

### 🧪 Experimental Design

The workflow uses experimental measurements from DGTX NISTCHO characterization runs covering two culture conditions:

- **Constant (CT)**
- **Temperature shifted (TS)**

Measurements are available across five experimental time windows and include uptake/secretion rates, growth, and transcriptomics data.

### 🔬 Analysis Workflow

The analysis follows four main stages:

1. **Original model comparison**  
   Load and diagnose the three CHO GEMs and establish model-specific reaction mappings.

2. **Experimental constraint & contextualization**  
   Apply experimentally measured exchange rates and contextualize the constrained models using transcriptomics with GIMME.

3. **Model simulation & validation**  
   Predict biomass and IgG production using FBA/pFBA and evaluate predictions against experimental measurements, including leave-one-exchange-out cross-validation.

4. **Flux analysis & visualization**  
   Characterize global and reaction-level flux patterns across models, conditions, and time windows.

## 🗂 Repository Structure
Scripts are numbered according to their recommended execution order. Utility scripts in `python/utils/` and `R/utils/` are not part of the sequential workflow.
  Original Models Comparison & Diagnostics
  - `python/00_create_SHA256.py` : Computes SHA-256 checksums of the original GEM model files to document model provenance and enable file-integrity verification.
  - `python/01_load_models.py` : Loads the three CHO genome-scale metabolic models using the repository utility functions and summarizes their model objectives. The script performs basic consistency checks on each model’s biomass reaction, identifies model-specific exchange reactions for experimentally measured metabolites, IgG production, and growth, and exports the corresponding reaction mappings as CSV and Pickle files for use in downstream analyses.
  - `python/02_model_diagnostics.py` : Performs cross-model diagnostics, assessing ATP production, blocked reactions, minimal media requirements, biomass and IgG feasibility, and precursor producibility across model contexts. Generates summary tables and diagnostic plots for model comparison.
  
  Experimental Constraints & Contextualization using GIMME
  - `python/03_set_bounds.py` : Applies experimental uptake and secretion rate constraints (±SD) to each GEM across experimental conditions and time windows, computes FBA/pFBA solutions, and saves the resulting constrained models for downstream contextualization and validation. Models saved as .xml and .mat in `models/exp_constrained`.
  - `R/04_preprocess_transcriptomics.R` : Preprocesses transcriptomics data by harmonizing model names, experimental conditions, and time windows, and combines the model-specific expression datasets into a unified GIMME-ready matrix with basic QC.
  - `Matlab/05_contextualize_GIMME.mat` : Contextualizes experimentally constrained CHO GEMs using GIMME and condition-specific transcriptomics. Applies a 25th-percentile expression threshold, generates context-specific models for each condition and time window, records contextualization diagnostics, and exports the resulting models and diagnostics.Generated models are saved in `models/contextualized/GIMME_tpm_25perc`.
  
  Models Solution & Validation
  - `python/06_maximise_biomass.py` : Maximizes biomass production in the selected CHO GEMs under experimentally measured condition- and window-specific uptake constraints. Runs both FBA and pFBA, records biomass growth rates and reaction fluxes, and exports the resulting solution and flux tables to analysis/models/.
  - `python/07_maximise_igg.py` : Maximizes igg production in the selected CHO GEMs under experimentally measured condition- and window-specific uptake constraints. Runs both FBA and pFBA, records igg production rates and reaction fluxes, and exports the resulting solution and flux tables to analysis/models/.
  - `python/08_leave_out_one_cv.py` : Performs leave-one-exchange-out cross-validation on the selected CHO GEMs by systematically removing each experimentally constrained exchange reaction, rerunning FBA and pFBA, and comparing the resulting biomass and reaction fluxes with the fully constrained baseline. Saves reaction-level cross-validation and leave-one-out results to analysis/models/.
  - `R/09_plot_exp_pred_rate.R` : Validates FBA and pFBA biomass growth predictions against experimentally measured growth rates for each CHO model, condition, and time window. Calculates R² and mean absolute error (MAE), generates publication- and poster-ready comparison and summary plots, and exports growth-validation statistics.
  - `R/10_validation_loocv.R` : Validates model-predicted exchange fluxes against experimentally measured uptake/secretion rates for each CHO GEM, condition, and time window. Maps measured metabolites to model-specific exchange reactions, calculates reaction-level R² and mean absolute error (MAE) for FBA and pFBA predictions, summarizes validation performance, and generates comparison and distribution plots.
  
  Metabolic Pathways & Final Plots
  - `R/11_global_flux_changes.R` : Analyzes predicted pFBA reaction fluxes across CHO GEMs, conditions, and time windows. Classifies reactions by activity and direction consistency, identifies consistently active and stable reactions, characterizes temporal flux trajectories using slope, range, cumulative change, and travel ratio, and performs PCA to assess global flux-state patterns. Generates reaction activity summaries, flux heatmaps, PCA plots, and trajectory/stability visualizations.
  - `R/12_reaction_level_dynamics.R` : Visualizes pFBA-predicted flux dynamics for selected central carbon metabolism, fermentation, transport, and TCA-cycle reactions in the iCHO3K model across constant and temperature-shifted conditions and five time windows. Generates individual reaction-level flux trajectory plots and a pathway-grouped heatmap with condition/time-window annotations.
  
  Helper scripts:
  - `python/utils/model_utils.py` : Provides reusable Python utilities for locating and loading CHO genome-scale models, summarizing model structure and objectives, standardizing model names, and parsing experimental condition and time-window information from model filenames. Supports original, experimentally constrained, and GIMME TPM contextualized model directories.
  - `R/utils/plotting_utils.R` : Defines global plotting parameters and reusable ggplot2 styling utilities for consistent figure formatting across the analysis. Provides standardized condition and event colour mappings, font sizes, line and point settings, significance annotation offsets, a clean theme_tidy() plotting theme, and a helper for converting point sizes to ggplot2 geometric units.


## 📦 Data Access 
  Original models were downloaded and SHA256 check-summed at the time of the download to ensure file integrity. 
  The original model files are located in /models/original. 
  
  Intermediate models and analysis results are organized as follows:
  | Directory | Contents |
|---|---|
| `data/rates/` | Experimentally measured exchange and growth rates |
| `data/transriptomics/` | Transcriptomics data already matched with genes in CHO GEMs |
| `models/original/` | Original CHO GEMs |
| `models/exp_constrained/` | Experimentally constrained models |
| `models/contextualized/` | Transcriptomics-contextualized models |
| `analysis/models/` | FBA/pFBA solutions, fluxes, and validation results |

## 🚀 Running the Workflow

The numbered scripts are intended to be executed sequentially from the repository root:

`00 → 01 → 02 → 03 → 04 → 05 → 06 → 07 → 08 → 09 → 10 → 11 → 12`

Some steps depend on external experimental datasets and MATLAB-based GIMME contextualization. Intermediate outputs generated by earlier steps are used as inputs for downstream analyses.

## 📝 Citation

## 🤝 Contact

#### Metabolic modelling & Computational analysis
**Veronika Schäpertöns**  
University of Salzburg
📧 **veronika.schaepertoens@plus.ac.at**  
🔗 **GitHub:** [@VSchaepertoens](https://github.com/VSchaepertoens) 
