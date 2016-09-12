# LGR_SimulationStudy
Simulation study related to estimating escapement past Lower Granite Dam

## File descriptions
*  SimFunc.R: functions that simulate fish crossing over Lower Granite Dam, and generate "observed" datasets on a weekly time-scale
*  RunSims.R: run multiple simulations and save the resulting estimates and truth (based on the simulations). Focus is on wild fish only
*  RunSims_AllOrigins.R: same as above, but this model includes estimates of natural origin, hatchery origin and hatchery no-clip fish
*  LGR_TotalEscape_JAGS.txt: JAGS model file to estimate wild escapement
*  LGR_TotalEscape_AllOrigins.txt: JAGS model file to estimate wild, hatchery and hatchery no-clip escapement
*  EvaluateSims.R: summarise the results of the simulations with tables and figures
*  SimulationFits: folder to store the results of simulations and estimates
