#------------------------------------------------------------------------------
# Estimate passage timing parameters for spring and summer Chinook runs.
# Run-timing parameters are then used to draw simulated population
# parameters from normal distributions.  The method produces theoretical
# run-timing curves for n.pop populations crossing Lower Granite Dam.
#------------------------------------------------------------------------------
# Ryan N. Kinzer
# Created: 7/8/2016
# Modified: 9/19/2016
#------------------------------------------------------------------------------
library(readxl)
library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)


#------------------------------------------------------------------------------
# Estimate run-timing for different stocks across LGR
#------------------------------------------------------------------------------
# Load Data - data was obtained from Rick Orme; includes all LGD
# pit detections from 2010-2015

dat <- read_excel("./Data/LGR chinook obs 2010-15.xlsx")
metadata <- read.csv("./Data/All MRR Sites Table.csv")

dat_names <- c("SbyC.fish","Tag.Code","unique.tag.year","Species","SRR.Code","Obs.Time.Value","Obs.Date",
               "Obs.Year","Spawn","Site.Name","Antenna.ID","Antenna.Location","Release.Date","Release.Site.Name",
               "Release.Site.Basin","Release.Site.Basin.Code","Release.Site.Subbasin.Code","Release.Site.Subbasin",
               "Release.Location","Release.Site.Code","Obs.Count")

names(dat) <- dat_names

metadata <- metadata %>%
              rename(Release.Site.Name = MRR.Site.Name)

# Create Week and Day (julian) fields from First.Obs.Date
#dat$Day = as.POSIXlt(dat$First.Obs.Date, format = "%m/%d/%Y")$yday
#datWeek = strftime(as.POSIXlt(dat$First.Obs.Date, format = "%m/%d/%Y"),format="%W")

#SRR <- c("Wild Spring Chinook", "Hat. Spring Chinook",
 #        "Wild Summer Chinook", "Hat. Summer Chinook",
  #       "Wild Fall Chinook",  "Hat. Fall Chinook",
   #      "Wild Summer Steelhead", "Hat. Summer Steelhead")

Basins <- c("Clearwater","Salmon","Lower Snake")

# Organize data
srr_df   <- left_join(dat,metadata) %>%
            filter(MRR.Site.Basin.Name %in% Basins & 
                   !is.null(Release.Site.Name)) %>%
            mutate(Origin = ifelse(grepl("H",SRR.Code),"Hatchery",
                                   ifelse(grepl("W",SRR.Code),"Natural","Unknown")),
                   Run = ifelse(grepl("11",SRR.Code),"Spring","Summer"),
                   Release.group = paste(Release.Site.Subbasin,Run,Origin,Obs.Year,sep="-"),
                   Day = yday(Obs.Date),
                   Date = paste(month(Obs.Date),mday(Obs.Date),sep="-"))

sp_sum_p <- srr_df %>%
              filter(!is.na(Species)) %>%
              group_by(Release.Site.Subbasin,Species,Run,Origin,Release.group,Obs.Year,Day) %>% # summarize by week
              summarise(n = n()) %>% # add up all observations at ISO
              mutate(cumsum_n = cumsum(n), # cumulative sum for ecdf
                     n_sample = sum(n), # yearly total
                     Notseen = n_sample - n, # create non seen for logit reg.
                     p = n/n_sample, 
                     ecdf = cumsum_n/n_sample)# %>%
              #filter(n_sample >= 5)

fig1 <- ggplot(data=sp_sum_p) +
        geom_line(aes(x=Day,y=ecdf,group=Release.group,colour= Run,linetype=Origin), size = 1)+
        scale_colour_brewer(palette = "Set1")+
        facet_wrap(~Obs.Year, scale="free_x")+
        theme_bw()
fig1

# Calculate the average passage day for spring and summer run chinook

srr_mean <-  srr_df %>%
              filter(!is.na(Species)) %>%
              group_by(Release.Site.Subbasin,Species,Run,Origin,Release.group,Obs.Year) %>% # summarize by week Species,Run,Origin,Release.code
              summarise(n_sample = n(),
                        mean = mean(Day),
                        stddev = sd(Day),
                        lower = mean - 1.96*stddev,
                        upper = mean + 1.96*stddev) %>%
              arrange(mean)

srr_mean$Release.Site.Subbasin = factor(srr_mean$Release.Site.Subbasin,levels=unique(srr_mean$Release.Site.Subbasin))

fig2 <- ggplot(data=srr_mean) +
          geom_point(aes(x = Release.Site.Subbasin, y = mean, colour=Run, shape = Origin),size = 3, position = position_dodge(width=1))+
          geom_errorbar(aes(x = Release.Site.Subbasin, ymin = lower, ymax = upper, colour = Run,linetype=Origin),size=1, position = position_dodge(width=1)) +
          scale_colour_manual(values=c("orange","blue"))+
          #facet_wrap(~Release.Site.Subbasin) +
          coord_flip() +
          theme_bw()
fig2

# Calculate means based on subbasin and run (exclude origin and different release codes)

srr_sub_mean <-  srr_df %>%
  filter(!is.na(Species)) %>%
  group_by(Release.Site.Subbasin,Species,Run,Obs.Year) %>% # summarize by week Species,Run,Origin,Release.code
  summarise(n_sample = n(),
            mean = mean(Day),
            stddev = sd(Day),
            lower = mean - 1.96*stddev,
            upper = mean + 1.96*stddev) %>%
  arrange(mean)

srr_sub_mean$Release.Site.Subbasin = factor(srr_sub_mean$Release.Site.Subbasin,
                                             levels=unique(srr_sub_mean$Release.Site.Subbasin))

fig3 <- ggplot(data=srr_sub_mean) +
  geom_point(aes(x = Release.Site.Subbasin, y = mean, colour=Run),size = 3, position = position_dodge(width=1))+
  geom_errorbar(aes(x = Release.Site.Subbasin, ymin = lower, ymax = upper, colour = Run),size=1, position = position_dodge(width=1)) +
  scale_colour_manual(values=c("orange","blue"))+
  coord_flip() +
  facet_wrap(~Obs.Year) +
  theme_bw()
fig3

# Calculate the average julian day of passage means for each population of spring and summer run chinook


run_mu <- srr_sub_mean %>%  # used to set hyper-parameters
        group_by(Release.Site.Subbasin) %>%
        summarise(n = n(),
                  mu = mean(mean),
                  mu_sd = sd(mean),
                  sd_mu = mean(stddev,na.rm=TRUE),
                  sd_sd = sd(stddev,na.rm=TRUE)) %>%
            arrange(mu)

# Prior Stuff
# Fit logistic model

library(lme4)
modfit <- glmer(formula = cbind(n,Notseen) ~ Day + I(Day^2) + (1 + Day|Release.code) ,
                family=binomial(link='logit'), data = srr_df)

summary(modfit)

# Create data frame for estimates
newdf <- expand.grid(Main.Branch.1 = unique(main_branch_grp$Main.Branch.1),
                    Day = unique(main_branch_grp$Day))

newdf$y.logit <- predict(modfit,newdata=newdf)
newdf$y.pred <- plogis(newdf$y.logit)

fig2 <- ggplot(data=newdf) +
  geom_line(aes(x=Day,y=y.pred,group=Main.Branch.1,colour=Main.Branch.1), size = 1)+
  #geom_point(aes(x=Week,y=p,colour=Main.Branch.1),size=2)+
  #facet_wrap(~Spawn.yr)+
  theme_bw()
fig2

# Create matrix of movement probabilities across LGR for each main branch.
LGR.df <- newdf %>%
            select(Week, Main.Branch.1, y.pred) %>%
            spread(Main.Branch.1,y.pred) %>%
            select(-Week)

weeks <- sort(unique(main_branch_grp$Week))

row.names(LGR.df) <- weeks

density.scale <- function(x){
  y <- x / sum(x)
  return(y)
}

LGR.mat <- apply(LGR.df,2,density.scale)
matplot(weeks,LGR.mat,type="l")






day <- 1:250
n.day <- length(day)
n.pops <- 20

beta1 <- -27.37
beta2 <- .3394
beta3 <- -.001135


beta <- c(beta1, beta2, beta3)
Sigma <- as.matrix(vcov(modfit))

beta.pops <- t(rmvnorm(n.pops,beta,Sigma))

day.p.mat <- matrix(data=NA,nrow=n.day,ncol=n.pops)

for(i in 1:n.day){
  for(j in 1:n.pops){
    day.p.mat[i,j] <- plogis(beta.pops[1,j] + beta.pops[2,j]*day[i] + beta.pops[3,j]*(day[i]^2))
  }
}

matplot(day.p.mat,xlim=c(75,225))
colSums(day.p.mat)

# matrix(c(2.77,-2.01,.00001,
#                   -.0201, .0001, -.0000001,
#                   .00001, -.0000001, .0000000005),nrow=3,byrow=TRUE)

# (Intercept)           Day      I(Day^2)
# (Intercept)  2.770388e+00 -2.013596e-02  1.261822e-05
# Day         -2.013596e-02  1.577248e-04 -1.649454e-07
# I(Day^2)     1.261822e-05 -1.649454e-07  5.347542e-10

pop.mu.int <- 0
pop.sd.int <- 1#5.42271
pop.beta <- 0
pop.sd.slope <- .001#.03564



beta1.pop <- rnorm(n.pops,beta1,pop.sd.int)
beta2.pop <- rnorm(n.pops,beta2,pop.sd.slope)


hist(beta1.pop)
hist(beta2.pop)

y <- beta1 + beta2*day + beta3*(day^2)
plot(day,plogis(y),type="l")



matplot(day.p.mat)


dmvnorm(x=c(0,0))
dmvnorm(x=c(0,0), mean=c(1,1))

sigma <- matrix(c(4,2,2,3), ncol=2)
x <- rmvnorm(n=500, mean=c(1,2), sigma=sigma)
colMeans(x)
var(x)

x <- rmvnorm(n=500, mean=c(1,2), sigma=sigma, method="chol")
colMeans(x)
var(x)

plot(x)




