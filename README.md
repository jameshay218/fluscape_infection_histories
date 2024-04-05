# Reconstructed influenza A/H3N2 infection histories reveal variation in incidence and antibody dynamics over the life course
 ------------
Contains all of the code, scripts and data files to reconstruct the analyses and figures in [this paper](https://doi.org/10.1101/2024.03.18.24304371).

## Quick summary
The main `serosolver` fits can be run with `scripts/serosolver_fluscape_main.R`. This will take a very long time. You can also run the fits to the Ha Nam data (`scripts/serosolver_ha_nam_fits.R`) which will complete in much less time (still ~1 hour).

You can access the infection history posterior draws directly from the the RData files in `manuscript/r_data`. These will load in an object called `inf_chain`, which gives the posterior estimates of infection history states for each individual (`i` = individual; `j` = time period where 1 = Q1-1968; `x` = infection state; `sampno` and `chain_no` = posterior draw IDs). Note that missing combinations of `i`, `j`, `sampno` and `chain_no` imply a draw for `x`=0.


## Overview
The workflow is in three stages:

1) Cleaning the Fluscape data and antigenic map (internal use -- requires unpublished data)
2) Running the suite of simulation-recovery experiments
3) Fitting _serosolver_ to the Fluscape and Ha Nam data with various settings
4) Post-processing to generate summary statistics and plots (in the `manuscript` folder)

Although the scripts for part 1) are included, only the final, clean datasets are included in this repository. You may therefore skip to part 3).

All of part 4) can be run from the `manuscript/scripts/all_analyses.R` script, though you may wish to dive into the auxiliary scripts `manuscript/scripts/aux` to see how subsets of the analyses are performed.

INTERNAL USE: the `fluscape_serosolver` Git repo contains more complete versions of these scripts using the non-public data. Authors should use this version to reproduce the main manuscript results and figures. Only the final manuscript figures are included here.

## 1. Setup
It is recommended to clone this Git repository to your PC and set this as your working directory.

The main, custom R package used for these analyses is _serosolver_. This can be installed using `devtools::install_github("seroanalytics/serosolver",ref="fluscape")`. The most recent version of the code can be found [here](https://github.com/seroanalytics/serosolver), but you will need to use the `fluscape` branch. Note that you will need a C++ toolchain installed.

A number of other R packages are also used throughout:
```r
c("coda", "data.table", "plyr", "bayesplot", "ggplot2", "dplyr", 
"reshape2", "RColorBrewer", "gtable", "mgcv", "itsadug", "tidyverse", 
"patchwork", "viridi", "ggpubr", "NatParksPalettes", "lattice", 
"sp", "raster", "maptools", "ncf", "Hmisc", "corrplot", "gganimate", 
"cowplot")
```

NOTE: I originally wrote most of this code without `dplyr` and later converted to using the tidyverse. Therefore, there is a mix of tidyverse and pre-tidyverse syntax. This is largely a non-issue, but sometimes certain functions have conflicting namespaces (e.g., `select`, `summarize`). If you encounter an error like this, please just try to update the code with e.g., `select <- dplyr::select` at the start of the scripts before flagging an issue. It is likely a simple fix.

## 2. Data cleaning (internal use only)
### 2.1 HI titers
To clean the data from scratch, you will need the `fluscape` Git repository, which is currently only accessible by Fluscape collaborators. However, the outputs of `scripts/extract_data_initial.R` are available in the `data` folder (`fluscape_data_1_resolution.csv` and `fluscape_data_4_resolution.csv`, where the number refers to the time units) and all that is required to run the entirety of the analyses.

The Ha Nam cohort data can be found at `data/vietnam_data_primary.csv`. The extraction script is `scripts/extract_vietnam_data.R`, and uses the data as provided by Kucharski et al. PLOS Biology 2018 [here](https://github.com/adamkucharski/flu-model).

### IMPORTANT NOTE ON PUBLIC DATA
The age information in the Fluscape participant data have been slightly altered in the published version for privacy reasons -- for now, we have bucketed date of birth into 1-year bands and anonymised DOB for particularly old individuals (i.e., full DOB is removed). We have also modified the infection history posteriors post-hoc to prevent re-identification of date of birth. There will be some slight discrepancies in age-based estimates between the results in the paper and those obtained from re-running these scripts. You may also run into some errors with missing variables (e.g., running the `data_regressions.R` script). None of these errors will affect the ability to reproduce any analyses in the main results. We have also added random noise to the study location latitude and longitude coordinates. (Internal note -- this is the `remove_age_info.R` script).

### 2.2 Antigenic maps
The script `make_antigenic_maps.R` uses antigenic coordinates from our own data and Fonville et al. _Science_ 2014 to plot antigenic maps. We fit cubic smoothing splines to these coordinates and extrapolate to generate antigenic coordinates for each possible strain circulation time. The script can be modified to change the smoothing parameters, time resolution (ie. infection histories per year, quarter, month), and also to cluster strains by antigenic cluster or time period. Note that antigenic coordinates from Bedford et al. _eLife_ 2014 are also included here for comparison and the simulation-recovery experiments (see the figures in `data/antigenic_maps`), noting that the Fonville et al. coordinates are used the rest of the analyses.

## 3. Simulation-recovery experiments
All of the simulation-recovery experiments described in the long appendix can be run using the script `scripts/simulation_recovery/run_all_simulation.R`. This script runs each of the analyses from the numbered scripts in the same folder. The first script generates simulated serosurvey data which are saved in `data/simulated`. To run just the analyses presented in the supplementary materials, use the `supplement_1.R` and `supplement_2.R` scripts.

## 4. serosolver fits
All of the scripts to run the analyses described below are in the `scripts` folder. These scripts will generate MCMC chains saved to disk in the `chains` folder. The original chain files are very large files and thus are ignored by Git. Instead, smaller, pre-processed versions of the chains used for the main manuscript results are available in the `manuscript/r_data` folder.

### 4.0 Checking MCMC convergence
As many of these scripts take a long time to run, it can be helpful to check their progress. _serosolver_ periodically saves the current MCMC chain state to disk, so you can manually inspect the chain csv files in the `chains` folder. You can do this more formally using the `scripts/check_convergence.R` script, which reads in the chain csv files and plots density and MCMC trace plots.

### 4.1 Ha Nam data
Fits to the Ha Nam, Vietnam dataset can be generated using `scripts/serosolver_ha_nam_fits.R`. This only takes about 45 minutes to run locally (for the full 1.5 million iterations), so no HPC is needed. This uses the same antigenic map as the Fluscape fits later on, but note that a different map can be used. It is also possible to set the number of iterations lower in the `run_MCMC` call (set `iterations` and `adaptive_period`) to run a shorter chain.

There is also a simulation-recovery script designed to mimic the data dimensions of the Ha Nam data. This can be run in full using `scripts/serosolver_vietnam_sim_recovery.R`. This shows how the infection histories and antibody kinetics parameters are re-estimated accurately using a serosurvey with dimensions similar to the real Ha Nam data.

### 4.2 Base Fluscape data (full time resolution)
Fits of the base model (i.e., without considering strain-specific measurement offsets) to the Fluscape data can be generated using the `scripts/serosolver_fluscape_base_fits.R` script. This generates the results shown in the supplement of the manuscript, estimating infection histories in 3-month windows. This takes ~7 days to generate converged chains with high effective sample sizes.

### 4.3 Estimation of measurement offsets
As described in the manuscript, fitting the base _serosolver_ model to the Fluscape data gives clear over- and under-estimation of attack rates in some time periods. This appears to be a result of some A/H3N2 strains giving systematically higher or lower HI titer measurements, which biases the attack rates higher or lower respectively. To account for this, we fit a simpler version of the _serosolver_ model (using data at an annual resolution and using prior version 3) with strain-specific measurement offset terms estimated as free parameters (i.e., each observed titer against strain A is shifted up or down by X_A HI units). This analysis can be run with `scripts/serosolver_measurement_offset_estimates.R`, which takes ~1 day to run.

Once these fits are generated and converged chains have been saved into the `chains` folder, the script `scripts/extract_measurement_offset_estimates.R` is run to generate point estimates for each strain-specific offset term. These estimates are saved to `inputs/par_tab_offset_estimates.csv` and are used in section 4.4.

For particularly interested readers: we did attempt to estimate the measurement offset terms jointly with the main _serosolver_ fits, but were unable to generate converged chains due to identifiability issues. However, this approach did seem to work using simulated data. Scripts to fit this version of the model to the Fluscape data are `scripts/serosolver_joint.R` and `scripts/serosolver_joint_annual.R`, but are only included for completeness and are not used in the manuscript.

### 4.4 Full model using Fluscape data
We re-fit the full _serosolver_ model as described in section 4.2 but with the inclusion of *fixed* strain-specific measurement offsets. These are the main results presented in the manuscript.

The main script for this analysis is `scripts/serosolver_fluscape_main.R`, which takes ~7 days to run 50,000,000 iterations. An alternative version of this analysis at a coarser time resolution are also available: `scripts/serosolver_fuscape_main_annual.R` fits the same model but estimating infection histories for each 1-year window, and `scripts/serosolver_fuscape_main_mixedtime.R` estimates infection histories for each 1-year window prior to 2005, and per 3-month window thereafter.

## 5.0 Generating results and figures
All of the code and files required to reproduce the figures and statistics in the manuscript are in the `manuscript` folder. The main file to run is `manuscript/scripts/all_analyses.R`, which should run each of the analyses in turn by sourcing each respective auxiliary script in `manuscript/scripts/aux`. Some of these scripts take a while to run.

One crucial point is that you should not need to run most of the processing and data cleaning scripts. The `all_analyses.R` script will look for .RData objects in `manuscript/r_data` and only run certain scripts if these pre-processed objects do not exist -- these are `fluscape_dat.RData` and `measurement_quarter.RData`.

Although all of the required code and files are included, there will inevitably be some bugs and missing file dependencies (I have tried my best!). If you are trying to re-run everything but running into errors, please contact me or raise an issue here.
