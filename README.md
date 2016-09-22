# LGR_SimulationStudy
Simulation study related to estimating escapement past Lower Granite Dam

## File descriptions
*  SimFunc.R: functions that simulate fish crossing over Lower Granite Dam, and generate "observed" datasets on a weekly time-scale
* SimScenarios: folder containing the R script for each scenario. Each script generates 500 similar simulations, and saves the resulting estimates from the ISEMP and SCOBI models.
*  SimulationFits: folder to store the results of simulations and estimates
* Run_AllSims.R: run all the scenario scripts, and save outputs
*  SummariseResults.Rmd: summarise the results of the simulations with tables and figures
* SCOBI: contains the R functions to run the SCOBI model (the TAC estimate)

### JAGS models
*  LGR_TotalEscape_JAGS.txt: JAGS model file to estimate wild escapement
*  LGR_TotalEscape_AllOrigins.txt: JAGS model file to estimate wild, hatchery and hatchery no-clip escapement
*  LGR_TotalEscape_NB2.txt: JAGS model file with extra over-dispersion parameter for the negative binomial, to help estimate the mean-variance relationship. Currently not tested or used.



### Outdated scripts
*  RunSims.R: run multiple simulations and save the resulting estimates and truth (based on the simulations). Focus is on wild fish only
*  RunSims_AllOrigins.R: same as above, but this model includes estimates of natural origin, hatchery origin and hatchery no-clip fish
* EvaluateSims.R: pre-cursor to SummarizeResults.Rmd
