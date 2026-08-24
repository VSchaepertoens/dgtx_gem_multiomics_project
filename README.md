# Multi-Model CHO GEM Contextualization & Validation using data from DGTX NISTCHO Characterization Runs
## 📄 Overview

## 🗂 Repository Structure
  Models Comparison & Diagnostics
  - `python/00_create_SHA256.py` :
  - `python/01_load_models.py` : 
  - `python/02_model_diagnostics.py` : 
  
  Experimental Constrains & Contextualization using GIMME
  - `python/03_set_bounds.py` : 
  - `R/04_preprocess_transcriptomics.R` : 
  - `Matlab/05_contextualize_GIMME.R` :
  
  Models Solution & Validation
  - `python/06_maximise_biomass.py` : 
  - `python/07_maximise_igg.py` : 
  - `python/08_leave_out_one_cv.py` : 
  - `R/09_plot_exp_pred_rate.R` : 
  - `R/10_validation_loocv.R` :
  
  Metabolic Pathways & Final Plots
  - `R/11_global_flux_changes.R` : 
  - `R/12_reaction_level_dynamics.R` : 
  
  Helper scripts:
  - `python/utils/model_utils.py` : 
  - `R/utils/plotting_utils.R` : 


## 📦 Data Access 
  Original models were downloaded and SHA256 check-summed at the time of the download to ensure file integrity. 
  The original model files are located in /models/original. 
  
  All intermediate files, including 

## 📝 Citation

## 🤝 Contact

#### Metabolic modelling & Computational analysis
**Veronika Schäpertöns**  
University of Salzburg
📧 **veronika.schaepertoens@plus.ac.at**  
🔗 **GitHub:** [@VSchaepertoens](https://github.com/VSchaepertoens) 
