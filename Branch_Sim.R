#------------------------------------------------------------------------------
# Simulate 9 branches and 1 black box data set for testing branching model.
#------------------------------------------------------------------------------
# Ryan N. Kinzer
# Created: 9/23/2016
# Modified: 9/26/2016
#------------------------------------------------------------------------------

#source("./Script/lgd_run_timing.R") # generates the run timing parameters
# shown below.

# Branch run timing parameters
# A tibble: 3 × 6
#Group     n       mu    mu_sd    sd_mu    sd_sd
#<fctr> <int>    <dbl>    <dbl>    <dbl>    <dbl>
#1  Early    86 150.5457 19.91180 14.25784 6.954262
#2    Mid    42 170.9983 12.80495 22.21365 7.361752
#3   Late    21 182.9837 19.21272 22.16461 6.590507

#Black box run timing parameteres
# A tibble: 1 × 6 
#Species     n       mu    mu_sd    sd_mu    sd_sd
#<chr> <int>    <dbl>    <dbl>    <dbl>    <dbl>
#1 Chinook   149 160.8827 21.98749 17.52383 7.999656


SimulateBranchData <- function(trap.rate.df){

  source("SimFnc.R") # source LGR simulation function

  branch.size <- c(rep(250,3),rep(1000,3),rep(3000,3),12250)
  n.branch.pops <- c(rep(1,9),9)
  branch.det <- data.frame(lower = c(rep(c(.5,.5,.95),3),0),
                         upper = c(rep(c(.5,.95,.95),3),0))
  
  valid.list <- list()

  for(i in 1:length(branch.size)){
    tmp <- SimulateLGRdata(N.lgr = branch.size[i],
                           h.prob = 0,
                           hnc.prob = 0,
                           n.pops = n.branch.pops[i],
                           hatch.pop.prob = 0.5,
                           hnc.pop.prob = 0.15,                           
                           run.mu.mu = 171,
                           run.mu.sd = 13,
                           run.sd.mu = 22,
                           run.sd.sd = 7,
                           fallback.rate = 0.0,
                           reascension.rate = 1,
                           night.passage.rate = .06,
                           window.rate = 50/60,
                           marked.rate = 0.05,
                           ladder.det = 0.99,
                           start.date = ymd('20150101'),
                           trap.rate.df)

  branch.df <- tmp$sim %>%
    mutate(Lower.obs = rbinom(branch.size[i],1,branch.det[i,1]),
           Upper.obs = rbinom(branch.size[i],1,branch.det[i,2])) %>%
    filter(Trap == 1)

  valid.list[[i]] <- branch.df
} # end iloop

names(valid.list) <- c(paste0("Branch-",1:9),"Black-Box")

valid.df <- ldply(valid.list) %>%
          select(Branch = .id, everything())

return(valid.df)

} # end function

#--------------------------------------------------------------
# Test the function
#--------------------------------------------------------------

my_trap_rate = data.frame(Week = 1:52,
                          trap.rate = 0.08, # 0.16
                          trap.open = T)

# my_trap_rate %<>%
#   mutate(trap.rate = ifelse(Week %in% 30:32, 0, trap.rate),
#          trap.open = ifelse(Week %in% 30:32, F, trap.open))

valid.df <- SimulateBranchData(trap.rate.df = my_trap_rate)

valid.df %>%
  group_by(Branch) %>%
  summarise(n.tags = n(),
            lower = sum(Lower.obs)/n(), #obs. detection probs. at both sites
            upper = sum(Upper.obs)/n())

dim(valid.df)[1] # total valid tags
