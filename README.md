# GRP-GCB
Grinnell Resurvey Project Code (submission 1)

This respository contains data and code used for analysis in manuscript "A cloudy forecast for predicting species distributions..." submitted to Global Change Biology.  

## Contents 

**Inputs.RData** : This includes data required to fit nimble models. Specific contents include "Site_Data" (a data-frame populated with germane site-level "historic" [1900-1940] and "modern" [1980-2020] predictors), "y" (an array of site by species by era detection/non-detection data), and a species-list that corresponds to the order of indexing used for the observation data. 

**Mods.R** : This script includes nimble code used to fit models using MCMC. Note, not all (~40) models are included for space and organization. Rather, we present a set of 6 distinct models where both site-specific latent variables and species-specific coeffecients are time-varying (assuming that readers will be able to determine how to simplify the models to omit these terms or hold the coefficients time-invariant). 

