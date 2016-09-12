# LGR Simulation Script
library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(mvtnorm)
library(MCMCpack)

set.seed(21)
N.lgr <- 100000 # Lower Granite total escapement

n.pops <- c(17,5) # number of spring pops and summer pops

# Hyper-parameters for mean run-timing across spring and summer populations
run.mu.mu <- c(140, 171) # spring mu, summer mu
run.mu.sd <- c(16, 15) # spring sd, summer sd

# Hyper-parameters for sd of run-timing across spring and summer populations
run.sd.mu <- c(16,25)
run.sd.sd <- c(8,2)

# R.V. of population mean run-timing given assigned hyper-parameters
pop.mu <- c(rnorm(n.pops[1],run.mu.mu[1],run.mu.sd[1]),
            rnorm(n.pops[2],run.mu.mu[2],run.mu.sd[2]))

# R.V. of population standard deviation of run-timing given assigned hyper-parameters
pop.sd <- abs(c(rnorm(n.pops[1],run.sd.mu[1],run.sd.sd[1]),
            rnorm(n.pops[2],run.sd.mu[2],run.sd.sd[2])))

#------------------------------------------------------------------------------
# RK's fuzzy thoughts and scratch paper - SHOULD DELETE BEFORE FINAL
# Hyper-parameters for population proporitons
# p_i ~ Dir(alpha_1,....alpha_d) is generated from X_i ~ Gamma(alpha_i, 1)
# where p_i = X_i / sum(X_i).
# If alpha_1,..., alpha_d = 1, then X_i are iid Exp(1) since Gamma(1,1)=Exp(1)
# E[Gamma] = alpha*beta
# E[Exp] = theta
# if alpha = 1, then E[p_i] = 1/sum(n.pops)
# hist(rexp(1000,1))
# hist(rpois(1000,.5))
# x <- 4
# 1 / (sum(n.pops) + x)
# (1 + x)/ (sum(n.pops) + x)
# 
# hist(rgamma(1000,2,1))
#------------------------------------------------------------------------------

# Hyper-parameters for population proportions (includes all origins with pop.)
pop.p.mu <- rdirichlet(1,alpha=rep(1,sum(n.pops))) 


# True population escapement at Lower Granite Dam (includes all origins)
N.pop <- rmultinom(1,N.lgr,pop.p.mu)

# True population proportions crossing LGR
pop.p <- N.pop/N.lgr


# Assign Origin - Wild, Hatchery, Hatchery No Clip
w.alpha  <- rep(1,sum(n.pops))
h.alpha <- rbinom(sum(n.pops),1,.5)
nc.alpha <- rbinom(sum(n.pops),1,.2)

origin.alpha <- cbind(w.alpha,h.alpha,nc.alpha)

origin.pop.p <- matrix(NA,nrow=sum(n.pops),ncol=3)
N.pop.origin <- matrix(NA,nrow=sum(n.pops),ncol=3)
Origin <- NULL

# NEED TO NOW SPLIT COLUMN VALUES INTO INDIVIDUAL RECORDS.

for(i in 1:sum(n.pops)){
  origin.pop.p[i,] <- rdirichlet(1,alpha=origin.alpha[i,])
  N.pop.origin[i,] <- t(rmultinom(1,N.pop[i],origin.pop.p[i,]))
  tmp <- c(rep("NOR",N.pop.origin[i,1]),rep("HOR",N.pop.origin[i,2]),rep("HNC",N.pop.origin[i,3]))
  Origin <- c(Origin,tmp)
}

# Empty list to collect passage dates of individuals by population
obs.day <- list()

# Generate random day of passage for each population
for(i in 1:sum(n.pops)){
  obs.day[[i]] <- round(rnorm(N.pop[i],pop.mu[i],pop.sd[i]),0)
}

names(obs.day) <- paste0("Pop-",1:sum(n.pops))

# Create data frame to start appending information
LGR.df <- ldply(obs.day,data.frame, .id = 'Population') %>% tbl_df() %>%
  rename(Day = `X..i..`) %>%
  mutate(Week = ceiling((Day / 365)*52),
         Origin = factor(Origin,levels=c("NOR","HOR","HNC")),
         Run = ifelse(as.integer(gsub('Pop-', '', Population)) <= n.pops[1], 'Spring', 'Summer'),
         Run = factor(Run, levels = c('Spring', 'Summer')),
         Group = paste(Population, Origin, sep = '-'))


# Create data frame to start appending information
LGR.df <- ldply(obs.day,data.frame)

names(LGR.df) <- c("Population","Day")

LGR.df$Week <- ceiling((LGR.df$Day/365)*52)

LGR.df$Origin <- factor(Origin,levels=c("NOR","HOR","HNC"))

LGR.df$Run <- c(rep("Spring",sum(N.pop[1:n.pops[1]])),
                rep("Summer",sum(N.pop[(n.pops[1]+1):(n.pops[1]+n.pops[2])])))

LGR.df$Run <- factor(LGR.df$Run,levels=c("Spring","Summer"))
LGR.df$Group <- paste0(LGR.df$Population,"-",Origin)

fig_lgr_pop <- ggplot(data=LGR.df, aes(x = Day, group = Group, colour = Origin))+
                stat_ecdf(size=1)+
                scale_colour_manual(name="",values=c("green","blue","red"))+
                #facet_wrap(~Population)+
                theme_bw()
fig_lgr_pop


# Set other parameters.
fallback.rate <- .05 # probability of fish fall-back
reascension.rate <- 1 # conditional probability of fish reascending given they fall-back
                      # need to think about Tucannon R., reascension does not equal 1
                      # possibly set at some weighted mean from previous studies
night.passage.rate <- .05  # probability of fish passing during non window operation
window.rate <- 50/60 # minutes per hour

#window.obs <- .98 # probability observed at window; 1-window.obs is probability missed
# need to account for window observation error at a different point in simulation, account
# when tallying daily observation and place some error rate around the count which can 
# either be positively or negatively biased

marked.rate <- .05 # not sure what this should be; should we make the mark rate different by
                  # origin and population, or does one rate apply to all?
ladder.det <- .99

true.trap.rate <- c(0,.10,.15,.10,0,.10,0)
n.trap.rate <- c(3,10,10,10,2,10,7) #length should equal 52 weeks, and be same
                                    # dimensions of true.trap.rate
weekly.trap.rate <- c(rep(true.trap.rate[1],n.trap.rate[1]),
                      rep(true.trap.rate[2],n.trap.rate[2]),
                      rep(true.trap.rate[3],n.trap.rate[3]),
                      rep(true.trap.rate[4],n.trap.rate[4]),
                      rep(true.trap.rate[5],n.trap.rate[5]),
                      rep(true.trap.rate[6],n.trap.rate[6]),
                      rep(true.trap.rate[7],n.trap.rate[7]))

weekly.trap.rate = rep(0.07, 52)

trap.rate.df = data.frame(Week = 1:52,
                          trap.rate = 0.07) %>% tbl_df()


Fallback <- NA
Reascent <- NA
Night.passage <- NA
Window.passage <- NA
Marked <- NA
Ladder <- NA
Trap <- NA
#Window.observation <- NA


for(i in 1:N.lgr){
  
    Fallback[i] <- rbinom(1,1,fallback.rate)
  
    if(Fallback[i] == 1){Reascent[i] <- rbinom(1,1,reascension.rate)} else
                        {Reascent[i] <- 0}
    
    Night.passage[i] <- rbinom(1,1,night.passage.rate)
    
    if(Night.passage[i] == 0){Window.passage[i] <- rbinom(1,1,window.rate)} else
                             {Window.passage[i] <- 0}

    Marked[i] <- rbinom(1,1,marked.rate)

    Trap[i] <- rbinom(1,1,weekly.trap.rate[LGR.df$Week[i]])
    
    if(Marked[i] == 1 | Trap[i] == 1){Ladder[i] <- rbinom(1,1,ladder.det)} else
                                     {Ladder[i] <- 0}
    
    # Remove window.observation field.        
    #if(Window.passage[i] == 1){Window.observation[i] <- rbinom(1,1,window.obs)} else
    #                          {Window.observation[i] <- 0}  # account for obs error later.
}

LGR.df <- tbl_df(data.frame(LGR.df,cbind(Fallback,Reascent,Night.passage,Window.passage,
                                  Marked, Ladder, Trap))) #Window.observation))


LGR.df %<>% tbl_df() %>%
  mutate(Fallback = rbinom(N.lgr, 1, fallback.rate),
         Reascent = Fallback * rbinom(N.lgr, 1, reascension.rate),
         Night.passage = rbinom(N.lgr, 1, night.passage.rate),
         Window.passage = (1 - Night.passage) * rbinom(N.lgr, 1, window.rate),
         Marked = rbinom(N.lgr, 1, marked.rate)) %>%
  left_join(trap.rate.df) %>%
  mutate(Trap = rbinom(N.lgr, 1, trap.rate),
         Ladder = rbinom(N.lgr, 1, ladder.det))


# NEED TO REMOVE REASCENTS AND ADD TO BOTTOM OF LIST FOR WINDOW OBSERVATIONS

# Tweak the knobs at two population sizes: Low and High
# 1. All assumptions of both models are met.
# 2. Change trap rates across weekly sampling periods.
# 3. Trap closures.
# 4. mess with reascenion and night-time passage
# 5. Window count error.


