# Author: Kevin See
# Purpose: Simulate data to estimate total escapement over Lower Granite Dam
# Created: 8/10/2016
# Last Modified: 9/12/2016
# Notes: Based on code provided by Ryan Kinzer. Currently works on a daily time-step, then summarises it to a weekly time-step.

#-----------------------------------------------------------------
# LGR Simulation Function
SimulateLGRdata = function(N.lgr = 100000,
                           n.pops = c(17,5),
                           run.mu.mu = c(140, 171),
                           run.mu.sd = c(16, 15),
                           run.sd.mu = c(16,25),
                           run.sd.sd = c(8,2),
                           hatch.pop.prob = 0.5,
                           hnc.pop.prob = 0.15,
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
  library(plyr)
  library(dplyr)
  library(tidyr)
  
  # R.V. of population mean run-timing given assigned hyper-parameters
  pop.mu <- c(rnorm(n.pops[1], run.mu.mu[1], run.mu.sd[1]),
              rnorm(n.pops[2], run.mu.mu[2], run.mu.sd[2]))
  
  # R.V. of population standard deviation of run-timing given assigned hyper-parameters
  pop.sd <- abs(c(rnorm(n.pops[1], run.sd.mu[1], run.sd.sd[1]),
                  rnorm(n.pops[2], run.sd.mu[2], run.sd.sd[2])))
  
  # Hyper-parameters for population proportions (includes all origins with pop.)
  pop.p.mu <- rdirichlet(1, alpha = rep(1, sum(n.pops))) 
  
  # True population escapement at Lower Granite Dam (includes all origins)
  N.pop <- rmultinom(1, N.lgr, pop.p.mu)
  
  # Assign Origin - Wild, Hatchery, Hatchery No Clip
  w.alpha  <- rep(1, sum(n.pops))
  h.alpha <- ifelse(rbinom(sum(n.pops), 1, hatch.pop.prob)==1,2,0)
  nc.alpha <- rbinom(sum(n.pops), 1, hnc.pop.prob)
  
  origin.alpha <- cbind(w.alpha, h.alpha, nc.alpha)
  
  origin.pop.p <- matrix(NA, nrow=sum(n.pops), ncol=3)
  N.pop.origin <- matrix(NA, nrow=sum(n.pops), ncol=3)
  Origin <- NULL
  
  # NEED TO NOW SPLIT COLUMN VALUES INTO INDIVIDUAL RECORDS.
  
  for(i in 1:sum(n.pops)){
    origin.pop.p[i,] <- rdirichlet(1, alpha = origin.alpha[i,])
    N.pop.origin[i,] <- t(rmultinom(1, N.pop[i], origin.pop.p[i,]))
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
    mutate(Week = ceiling((Day / 365)*52),
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
  
  LGR.reascents = LGR.df %>%
    filter(Reascent == 1) %>%
    mutate(Day = Day + 1 + rpois(nrow(.), 2),
           Week = ceiling((Day / 365)*52),
           Fallback = rbinom(nrow(.), 1, fallback.rate),
           Reascent = Fallback * rbinom(nrow(.), 1, reascension.rate),
           Night.passage = rbinom(nrow(.), 1, night.passage.rate),
           Day.passage = 1 - Night.passage,
           Window.passage = Day.passage * rbinom(nrow(.), 1, window.rate)) %>%
    left_join(trap.rate.df) %>%
    mutate(Trap = rbinom(nrow(.), 1, trap.rate),
           Ladder = Marked * rbinom(nrow(.), 1, ladder.det)) %>%
    arrange(Day, Group, Origin)
  
  LGR.reascents_2 = LGR.reascents %>%
    filter(Reascent == 1) %>%
    mutate(Day = Day + 1 + rpois(nrow(.), 2),
           Week = ceiling((Day / 365)*52),
           Fallback = rbinom(nrow(.), 1, fallback.rate),
           Reascent = Fallback * rbinom(nrow(.), 1, reascension.rate),
           Night.passage = rbinom(nrow(.), 1, night.passage.rate),
           Day.passage = 1 - Night.passage,
           Window.passage = Day.passage * rbinom(nrow(.), 1, window.rate)) %>%
    left_join(trap.rate.df) %>%
    mutate(Trap = rbinom(nrow(.), 1, trap.rate),
           Ladder = Marked * rbinom(nrow(.), 1, ladder.det)) %>%
    arrange(Day, Group, Origin)
  
  LGR.reascents_3 = LGR.reascents_2 %>%
    filter(Reascent == 1) %>%
    mutate(Day = Day + 1 + rpois(nrow(.), 2),
           Week = ceiling((Day / 365)*52),
           Fallback = rbinom(nrow(.), 1, fallback.rate),
           Reascent = Fallback * rbinom(nrow(.), 1, reascension.rate),
           Night.passage = rbinom(nrow(.), 1, night.passage.rate),
           Day.passage = 1 - Night.passage,
           Window.passage = Day.passage * rbinom(nrow(.), 1, window.rate)) %>%
    left_join(trap.rate.df) %>%
    mutate(Trap = rbinom(nrow(.), 1, trap.rate),
           Ladder = Marked * rbinom(nrow(.), 1, ladder.det)) %>%
    arrange(Day, Group, Origin)
  
  LGR.df %<>%
    bind_rows(LGR.reascents) %>%
    bind_rows(LGR.reascents_2) %>%
    bind_rows(LGR.reascents_3) %>%
    mutate(Date = ymd('20150101') + days(Day)) %>%
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
    # # deal with reascents. Assume all reascensions happen within the same week as initial passage
    # left_join(org_truth %>%
    #             filter(Reascent == 1) %>%
    #             group_by(Week) %>%
    #             summarise(N_uniq_reasc = n_distinct(id),
    #                       Day.reascents = rbinom(1, N_uniq_reasc, 1 - night.passage.rate),
    #                       Night.reascents = N_uniq_reasc - Day.reascents,
    #                       Window.reascents = rbinom(1, Day.reascents, window.rate))) %>%
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
                          theta = 100,
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
    summarise_each(funs(sum), -Origin)
  
  # simulate window counts. Add some observation error if !perfect.window
  if(perfect.window) {
    lgr_week_obs %<>%
      mutate(win_cnt = round(Window.passage / window.rate))
  }
  
  if(!perfect.window) {
    lgr_week_obs %<>%
      mutate(win_cnt = rnegbin(nrow(.), mu = Window.passage / window.rate, theta = theta))  
  }
  
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
                spread(Origin, n_fish, fill = 0)) %>%
    mutate_each(funs(ifelse(is.na(.), 0, .)), Ladder:wild_fish) %>%
    left_join(trap_rate_mr %>%
                select(Week,
                       trap_rate = p,
                       trap_rate_se = p_se) %>%
                # set up parameters describing trap rate as a beta distribution
                mutate(trap_alpha = ((1 - trap_rate) / trap_rate_se^2 - 1 / trap_rate) * trap_rate^2,
                       trap_beta = trap_alpha * (1 / trap_rate - 1))) %>%
    mutate(trap_alpha = ifelse(is.na(trap_rate_se), 1, trap_alpha),
           trap_beta = ifelse(is.na(trap_rate_se), 1, trap_beta)) %>%
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
