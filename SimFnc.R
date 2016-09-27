# Author: Kevin See
# Purpose: Simulate data to estimate total escapement over Lower Granite Dam
# Created: 8/10/2016
# Last Modified: 9/12/2016
# Notes: Based on code provided by Ryan Kinzer. Currently works on a daily time-step, then summarises it to a weekly time-step.

#-----------------------------------------------------------------
# LGR Simulation Function
SimulateLGRdata = function(N.lgr = 100000,
                           h.prob = .70,
                           hnc.prob = .05,
                           n.pops = 25,
                           hatch.pop.prob = 0.5,
                           hnc.pop.prob = 0.15,                           
                           run.mu.mu = 171,
                           run.mu.sd = 13,
                           run.sd.mu = 22,
                           run.sd.sd = 7,
                           fallback.rate = .06,
                           reascension.rate = 1,
                           night.passage.rate = .06,
                           window.rate = 50/60,
                           marked.rate = 0.05,
                           ladder.det = 0.99,
                           start.date = ymd('20150101'),
                           trap.rate.df) {
  # N.lgr: Lower Granite total escapement
  # n.pops: number of spring pops and summer pops
  # run.mu.mu, run.mu.sd: Hyper-parameters for mean run-timing across spring and summer populations
  # run.sd.mu, run.sd.sd: Hyper-parameters for sd of run-timing across spring and summer populations
  # hatch.pop.prob: percentage of populations with hatchery populations
  # hnc.pop.prob: percentage of populations with hatchery no-clip populations
  # fallback.rate: probability of fish fall-back
  # reascension.rate: conditional probability of fish reascending given they fall-back
  # need to think about Tucannon R., reascension does not equal 1
  # possibly set at some weighted mean from previous studies
  # night.passage.rate: probability of fish passing during non window operation
  # window.rate: minutes per hour
  # marked.rate: proportion of fish previously tagged (i.e. marked)
  # ladder.rate: detection probability of PIT tag antenna in the ladder
  # start.date: date from which simulation begins
  # trap.rate.df: data.frame with Week and trap.rate columns
  
  # required packages to run this function
  library(MCMCpack)
  library(lubridate)
  library(plyr)
  library(dplyr)
  library(tidyr)
  library(magrittr)
  
  # quick clean-up of column names
  trap.rate.df %<>%
    rename(trap_open = trap.open)
  
  # R.V. of population mean run-timing given assigned hyper-parameters
  # R.V. of population mean run-timing given assigned hyper-parameters
  pop.mu <- rnorm(n.pops, run.mu.mu[1], run.mu.sd[1])
  
  # R.V. of population standard deviation of run-timing given assigned hyper-parameters
  pop.sd <- abs(rnorm(n.pops, run.sd.mu[1], run.sd.sd[1]))
  
  # LGR escapement by origin
  if(h.prob > 0){
    N.hat <- N.lgr * h.prob
    hat.pop.p <- rdirichlet(1, alpha = rbinom(n.pops,1,hatch.pop.prob))  
    N.hat.pop <- rmultinom(1, N.hat, hat.pop.p)  
  } else {
    N.hat <- 0
    N.hat.pop <- rep(0,n.pops)
  }
  
  if(hnc.prob >0){  
    N.hnc <- N.lgr * hnc.prob
    hnc.pop.p <- rdirichlet(1, alpha = rbinom(n.pops,1,hnc.pop.prob)) 
    N.hnc.pop <- rmultinom(1, N.hnc, hnc.pop.p)
  } else {
    N.hnc <- 0
    N.hnc.pop <- rep(0,n.pops)
  } 

  N.wild <- N.lgr - (N.hat + N.hnc)
  wild.pop.p <- rdirichlet(1, alpha = rep(1, n.pops))
  N.wild.pop <- rmultinom(1, N.wild, wild.pop.p)

  # LGR escapement by origin and population
  N.pop.origin <- cbind(N.wild.pop, N.hat.pop, N.hnc.pop)
  
  # True population escapement at Lower Granite Dam (includes all origins)  
  N.pop <- rowSums(N.pop.origin)
  
  Origin <- NULL
  # NEED TO NOW SPLIT COLUMN VALUES INTO INDIVIDUAL RECORDS.
  
  for(i in 1:n.pops){
    tmp <- c(rep("NOR", N.pop.origin[i,1]), 
             rep("HOR",N.pop.origin[i,2]), 
             rep("HNC",N.pop.origin[i,3]))
    Origin <- c(Origin, tmp)
  }
  
  # Empty list to collect passage dates of individuals by population
  obs.day <- list()
  
  # Generate random day of passage for each population
  for(i in 1:sum(n.pops)){
    obs.day[[i]] <- round(rnorm(N.pop[i], pop.mu[i], pop.sd[i]), 0)
  }
  
  names(obs.day) <- paste0("Pop-", 1:sum(n.pops))
  
  # create data frame with row for each fish
  LGR.df <- ldply(obs.day,data.frame, .id = 'Population') %>% tbl_df() %>%
    rename(Day = `X..i..`) %>%
    mutate(Week = week(ymd(paste0(year(start.date), '0101')) + days(Day)),
           Origin = factor(Origin,levels=c("NOR","HOR","HNC")),
           Run = ifelse(as.integer(gsub('Pop-', '', Population)) <= n.pops[1], 'Spring', 'Summer'),
           Run = factor(Run, levels = c('Spring', 'Summer')),
           Group = paste(Population, Origin, sep = '-'),
           Fallback = rbinom(N.lgr, 1, fallback.rate),
           Reascent = Fallback * rbinom(N.lgr, 1, reascension.rate),
           Night.passage = rbinom(N.lgr, 1, night.passage.rate),
           Day.passage = 1 - Night.passage,
           Window.passage = Day.passage * rbinom(N.lgr, 1, window.rate),
           Marked = rbinom(N.lgr, 1, marked.rate)) %>%
    left_join(trap.rate.df) %>%
    mutate(Trap = rbinom(N.lgr, 1, trap.rate),
           Ladder = Marked * rbinom(N.lgr, 1, ladder.det),
           id = 1:n()) %>%
    arrange(Day, Group, Origin) %>%
    select(id, everything())
  
  if(sum(LGR.df$Reascent) > 0) {
    LGR.reascents = LGR.df %>%
      filter(Reascent == 1) %>%
      mutate(Day = Day + 1 + rpois(nrow(.), 2),
             Week = week(ymd(paste0(year(start.date), '0101')) + days(Day)),
             Fallback = rbinom(nrow(.), 1, fallback.rate),
             Reascent = Fallback * rbinom(nrow(.), 1, reascension.rate),
             Night.passage = rbinom(nrow(.), 1, night.passage.rate),
             Day.passage = 1 - Night.passage,
             Window.passage = Day.passage * rbinom(nrow(.), 1, window.rate)) %>%
      select(-trap.rate, -trap_open) %>%
      left_join(trap.rate.df) %>%
      mutate(Trap = rbinom(nrow(.), 1, trap.rate),
             Ladder = Marked * rbinom(nrow(.), 1, ladder.det)) %>%
      arrange(Day, Group, Origin)
  } else {
    LGR.reascents = LGR.df %>%
      slice(1) %>%
      mutate_each(funs(ifelse(!is.na(.), NA, NA)))
  }
  
  if(sum(LGR.reascents$Reascent, na.rm = T) > 0) {
    LGR.reascents_2 = LGR.reascents %>%
      filter(Reascent == 1) %>%
      mutate(Day = Day + 1 + rpois(nrow(.), 2),
             Week = week(ymd(paste0(year(start.date), '0101')) + days(Day)),
             Fallback = rbinom(nrow(.), 1, fallback.rate),
             Reascent = Fallback * rbinom(nrow(.), 1, reascension.rate),
             Night.passage = rbinom(nrow(.), 1, night.passage.rate),
             Day.passage = 1 - Night.passage,
             Window.passage = Day.passage * rbinom(nrow(.), 1, window.rate)) %>%
      select(-trap.rate, -trap_open) %>%
      left_join(trap.rate.df) %>%
      mutate(Trap = rbinom(nrow(.), 1, trap.rate),
             Ladder = Marked * rbinom(nrow(.), 1, ladder.det)) %>%
      arrange(Day, Group, Origin)
  } else {
    LGR.reascents_2 = LGR.df %>%
      slice(1) %>%
      mutate_each(funs(ifelse(!is.na(.), NA, NA)))
  }
  
  if(sum(LGR.reascents_2$Reascent, na.rm = T) > 0) {
    LGR.reascents_3 = LGR.reascents_2 %>%
      filter(Reascent == 1) %>%
      mutate(Day = Day + 1 + rpois(nrow(.), 2),
             Week = week(ymd(paste0(year(start.date), '0101')) + days(Day)),
             Fallback = rbinom(nrow(.), 1, fallback.rate),
             Reascent = Fallback * rbinom(nrow(.), 1, reascension.rate),
             Night.passage = rbinom(nrow(.), 1, night.passage.rate),
             Day.passage = 1 - Night.passage,
             Window.passage = Day.passage * rbinom(nrow(.), 1, window.rate)) %>%
      select(-trap.rate, -trap_open) %>%
      left_join(trap.rate.df) %>%
      mutate(Trap = rbinom(nrow(.), 1, trap.rate),
             Ladder = Marked * rbinom(nrow(.), 1, ladder.det)) %>%
      arrange(Day, Group, Origin)
  } else {
    LGR.reascents_3 = LGR.df %>%
      slice(1) %>%
      mutate_each(funs(ifelse(!is.na(.), NA, NA)))
  }
  
  LGR.df %<>%
    bind_rows(LGR.reascents) %>%
    bind_rows(LGR.reascents_2) %>%
    bind_rows(LGR.reascents_3) %>%
    filter(!is.na(id)) %>%
    mutate(Date = start.date + days(Day)) %>%
    group_by(Week) %>%
    mutate(Start_Date = min(Date)) %>%
    ungroup() %>%
    select(-Date) %>%
    select(id, Day, Week, Start_Date, everything())
  
  
  return(list('sim' = LGR.df,
              'parameters' = pairlist(N.lgr = N.lgr,
                                      n.pops = n.pops,
                                      run.mu.mu = run.mu.mu,
                                      run.mu.sd = run.mu.sd,
                                      run.sd.mu = run.sd.mu,
                                      run.sd.sd = run.sd.sd,
                                      hatch.pop.prob = hatch.pop.prob,
                                      hnc.pop.prob = hnc.pop.prob,
                                      fallback.rate = fallback.rate,
                                      reascension.rate = reascension.rate,
                                      night.passage.rate = night.passage.rate,
                                      window.rate = window.rate,
                                      marked.rate = marked.rate,
                                      ladder.det = ladder.det,
                                      start.date = start.date,
                                      trap.rate.df = trap.rate.df)))
  
}

#-----------------------------------------------------------------
# Summarise on weekly basis, by origin and total
SummariseWeekly = function(sim_params,
                           lgr_truth,
                           samp.origin = 'NOR') {
  
  # these parameters are used elsewhere as well as the simulation function
  night.passage.rate = sim_params$night.passage.rate
  window.rate = sim_params$window.rate
  
  org_truth = lgr_truth %>%
    filter(Origin == samp.origin)
  
  org_week = org_truth %>%
    group_by(Week) %>%
    summarise(N_tot = length(id),
              N_uniq = n_distinct(id)) %>%
    left_join(org_truth %>%
                group_by(Week) %>%
                summarise_each(funs(sum), Night.passage, Day.passage, Window.passage, Fallback, Reascent, Trap, Ladder) %>%
                ungroup()) %>%
    left_join(org_truth %>%
                filter(Marked == 1,
                       Ladder == 1) %>%
                group_by(Week) %>%
                summarise_each(funs(sum), Reascent, Night.passage) %>%
                rename(mark_reasc = Reascent,
                       mark_night = Night.passage)) %>%
    mutate_each(funs(ifelse(is.na(.), 0, .)), -Week) %>%
    mutate(Origin = samp.origin) %>%
    select(Week, Origin, everything())
  
  return(org_week)
}



#-----------------------------------------------------------------
# Simulate observed values
SimulateLGRobs = function(sim_params,
                          lgr_truth,
                          # theta = 100,
                          error_rate = 0.05,
                          perfect.window = F) {
  
  window.rate = sim_params$window.rate
  
  wild_week = SummariseWeekly(sim_params, lgr_truth, samp.origin = 'NOR')
  hatch_week = SummariseWeekly(sim_params, lgr_truth, samp.origin = 'HOR')
  hnc_week = SummariseWeekly(sim_params, lgr_truth, samp.origin = 'HNC')
  
  lgr_week_org = wild_week %>%
    bind_rows(hatch_week) %>%
    bind_rows(hnc_week)
  
  lgr_week_obs = lgr_week_org %>%
    group_by(Week) %>%
    summarise_each(funs(sum), -Origin) %>%
    left_join(lgr_truth %>%
                select(Week, trap_open) %>%
                distinct())
  
  # simulate window counts. Add some observation error if !perfect.window
  if(perfect.window) {
    lgr_week_obs %<>%
      mutate(win_cnt = round(Window.passage / window.rate))
  }
  
  # add observation error to daily window counts, then sum by week
  if(!perfect.window) {
    win_cnt_week = lgr_truth %>%
      group_by(Day, Week) %>%
      summarise_each(funs(sum), Window.passage) %>%
      ungroup() %>%
      mutate(win_cnt = round(rnorm(nrow(.), mean = Window.passage / window.rate, sd = error_rate * (Window.passage / window.rate))),
             # prevent counts from being negative
             win_cnt = ifelse(win_cnt < 0, 0, win_cnt)) %>%
      group_by(Week) %>%
      summarise_each(funs(sum), win_cnt)
    lgr_week_obs %<>%
      left_join(win_cnt_week)
  }
  
  # if(!perfect.window) {
  #   lgr_week_obs %>%
  #     mutate(win_cnt = rnegbin(nrow(.), mu = Window.passage / window.rate, theta = theta))
  # }
  
  # mark-recapture estimate of trap rate
  trap_mr = lgr_truth %>%
    mutate(id = 1:n()) %>%
    filter(Marked == 1) %>%
    mutate(R = ifelse(Trap == 1 & Ladder == 1, 1, 0)) %>%
    filter(!(Trap == 0 & Ladder == 0)) %>%
    select(id, 
           Week,
           M = Trap,
           C = Ladder,
           R)
  
  # put together list of capture histories
  caphist_list = dlply(trap_mr, .(Week), capHistConvert, cols2use = c('M', 'C', 'R'), in.type = 'individual', out.type = 'frequency')
  
  # estimate trap rate based on mark-recapture
  trap_rate_mr = ldply(caphist_list, function(x) {
    mod = try(closedp.t(x, dfreq=T), silent = T)
    if(class(mod)[1] == 'try-error') {
      return(data.frame(N_hat = NA, N_se = NA, p = NA, p_se = NA))
    }
    p_mean = inv.logit(coef(mod$glm[['Mt']])[2:3])
    p_se = deltamethod(list(~ 1 / (1 + exp(-x1)), ~ 1 / (1 + exp(-x2))), mean = coef(mod$glm[['Mt']])[2:3], cov = vcov(mod$glm[['Mt']])[2:3,2:3])
    return(data.frame(N_hat = mod$results['Mt', 'abundance'],
                      N_se = mod$results['Mt', 'stderr'],
                      p = p_mean[1], 
                      p_se = p_se[1]))
  }, .id = 'Week') %>% tbl_df()
  
  
  # add info from trap
  lgr_week_obs %<>%
    left_join(lgr_week_org %>%
                group_by(Week, Origin) %>%
                summarise(n_fish = sum(Trap)) %>%
                ungroup() %>%
                mutate(Origin = factor(Origin, levels = c('HOR', 'HNC', 'NOR')),
                       Origin = revalue(Origin,
                                        c('NOR' = 'wild_fish',
                                          'HOR' = 'hatch_fish',
                                          'HNC' = 'HNC_fish'))) %>%
                spread(Origin, n_fish, fill = 0))
  if(!'wild_fish' %in% names(lgr_week_obs)) lgr_week_obs %<>%
    mutate(hatch_fish = 0)
  if(!'hatch_fish' %in% names(lgr_week_obs)) lgr_week_obs %<>%
    mutate(hatch_fish = 0)
  if(!'HNC_fish' %in% names(lgr_week_obs)) lgr_week_obs %<>%
    mutate(HNC_fish = 0)
  
  lgr_week_obs %<>%
    mutate_each(funs(ifelse(is.na(.), 0, .)), Ladder:HNC_fish) %>%
    left_join(trap_rate_mr %>%
                select(Week,
                       trap_rate = p,
                       trap_rate_se = p_se) %>%
                # set up parameters describing trap rate as a beta distribution
                mutate(trap_alpha = ((1 - trap_rate) / trap_rate_se^2 - 1 / trap_rate) * trap_rate^2,
                       trap_beta = trap_alpha * (1 / trap_rate - 1))) %>%
    mutate(trap_alpha = ifelse(is.na(trap_rate_se), 1, trap_alpha),
           trap_beta = ifelse(is.na(trap_rate_se), 1, trap_beta),
           trap_alpha = ifelse(trap_open, trap_alpha, 1e-12),
           trap_beta = ifelse(trap_open, trap_beta, 1)) %>%
    select(Week:Window.passage, win_cnt, everything())
  
  return(list('by_origin' = lgr_week_org,
              'obs' = lgr_week_obs))
}


# Re-formats weekly $obs object (output from SimulateLGRobs function) for use as SCOBI input
formatSCOBI_inputs <- function(weekly_obs, LGR_truth){
  # window_obs = output from SimulateLGRobs function; object$obs
  # LGR_truth = output from SimulateLGRdata function; obeject$sim
  
  win_cnt_tmp = weekly_obs %>%
    select(WeekNumber = Week,
           Counts = win_cnt) 
    
  fish_data <- LGR_truth %>%
    filter(Trap == 1) %>%
    select(WeekNumber = Week, Rear = Origin, GenStock = Population) %>%
    mutate(Rear = revalue(Rear, c("NOR" = "W",
                                  "HOR" = "H",
                                  "HNC" = "HNC"))) %>%
    arrange(WeekNumber)
    
  window_count <- fish_data %>%
    group_by(WeekNumber,Rear) %>%
    summarise(n = n()) %>%
    spread(Rear,n) %>%
    ungroup() %>%
    full_join(win_cnt_tmp) %>%
    mutate_each(funs(replace(.,which(is.na(.)),0))) %>%
    arrange(WeekNumber) %>%
    mutate(W_age = round(W/2,0),
           Cumsum = cumsum(W_age),
           Collaps = cumsum(W_age >= 100)+1L)  %>%  # Does not collapse following IDFG's rule set.
    select(Strata = WeekNumber, Counts, Collaps)
  
  return(list('window_count' = as.data.frame(window_count),
              'fish_data' = as.data.frame(fish_data)))
}

negbin_theta <- function(mu, error.rate) {
  
  # mu = window.passage / window.rate
  # error.rate = coefficient of variation proportion, ie. 0.05, 0.15
  
  theta = mu^2 / ((error.rate*mu)^2 - mu) 

  return(theta)
}
