# Cycle through sims.

source.files <- c("RunSims_Baseline_WinErr_low.R","RunSims_NightReasc_WinErr_low.R","RunSims_NightReasc_WinErr_high.R")

for(i in 1:3){
  filename <- paste0("./SimScenarios/",source.files[i])
  source(filename)
  rm(list=setdiff(ls(), "source.files"))
}
