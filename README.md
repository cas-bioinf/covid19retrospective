# covid19retrospective

Accompanying package to the paper "Detailed disease progression of 213 patients hospitalized with Covid-19 in the Czech Republic: An exploratory analysis". 

Please do get in touch if you want to reuse some of the code (we'll help you through the messy parts) or use some of the data in your project (there are limits due to privacy, but we can at least run your analysis locally and send summary results back). Either file an issue here or write to martin.modrak@biomed.cas.cz

Using R + knitr to build the complete paper. Requires `brms` version 2.14.2 and installed either `rstan` or `cmdstanr` (note that both require installation steps beyond just running `install.packages`, see the respective documentation)

Main contents of the package:

- `manuscript` folder contains the R Markdown files used to build the manuscript and supplementary material, including high-level code used to fit all the models used. Notable files:
  - `manuscript.Rmd` main manuscript
  - `analysis_brms.Rmd`, `competing_risks.Rmd`, `hmm.Rmd`, `jm.Rmd`, `single_outcome_coxph.Rmd` include code to fit individual model classes
  - `risk_scores.Rmd` code to evaluate the prognostic models for disease severity
  - `supplementary_figures.Rmd` 
  - `*_res.csv` are cached results of the model fits used in the manuscript and supplementary figures 
  - `all_auc.csv`, `all_auc_orig.csv` cached AUC data for risk scores (for our dataset/for the original publication)
  - `Covid19.json` - bibliography in JSON format
- `R` additional sources, notably:
  - Code to fit rate-based HMM models using `brms` (all files containing `hmm` in their name)
  - Various helper function for the individual model classes and analyses `*_tools.R`
  - Code to transform raw data in Excel and some special formats into .csv files (`raw_data_preprocessing.R`, `special_lab_data.R`, `translation.R`, `data_mapping.R`)
  - Code to load the anonymized .csv files (`data_preprocessing.R`)
  - Helpers to run Simulation-Based calibration (`sbc.R`, `sampling_multi.R`, `evaluation_tools.R`)
- `src` additional C++ sources for posterior prediction
- `dev` codes used during development, notably:
  - `raw_data_preprocessing.Rmd` runs the actual transformation from Excel to .csv files and performs checks
  - `hmm_model_devel.Rmd` code to test the HMM model on simulated data, including Simulation-based calibration
- `public_data` data that can be made public
  - `Data_collection.xlsx`, `Data_Collection_Cz.xlsx` - the Excel files used for data collection in English/Czech
  - `Translations.xlsx` mapping between the English and Czech versions of the data collection files
  - `Grein_et_al.xlsx` transcribed disease progression data from the early Remdesivir paper by Grein et al. in NEJM that motivated this project
- `data_description.md` description of the processed dataset
- `research_questions.md` is the initial list of questions we started the project with. Not all of them were addressed, because we ended up unable to have good answers even for the simpler ones.
