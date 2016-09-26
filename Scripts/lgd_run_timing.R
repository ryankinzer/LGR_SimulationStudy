#------------------------------------------------------------------------------
# Estimate passage timing parameters for spring and summer Chinook runs.
# Run-timing parameters are then used to draw simulated population
# parameters from normal distributions.
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
                   Week = week(Obs.Date),
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

# Calculate means based on subbasin and run (exclude origin and different release codes)

srr_sub_mean <-  srr_df %>%
  filter(!is.na(Species)) %>%
  group_by(Release.Site.Subbasin,Species,Run,Obs.Year) %>% # summarize by week Species,Run,Origin,Release.code
  summarise(n_sample = n(),
            mean = mean(Day),
            stddev = sd(Day),
            lower = mean - 1.96*stddev,
            upper = mean + 1.96*stddev) %>%
  ungroup() %>%
  group_by(Obs.Year) %>%
  mutate(Group = factor(cut(mean,3,labels=FALSE),labels = c("Early","Mid","Late"))) %>%
  arrange(mean)

srr_sub_mean$Release.Site.Subbasin = factor(srr_sub_mean$Release.Site.Subbasin,
                                             levels=unique(srr_sub_mean$Release.Site.Subbasin))

fig2 <- ggplot(data=srr_sub_mean) +
  geom_point(aes(x = Release.Site.Subbasin, y = mean, colour=Group),size = 3, position = position_dodge(width=1))+
  geom_errorbar(aes(x = Release.Site.Subbasin, ymin = lower, ymax = upper, colour = Group),size=1, position = position_dodge(width=1)) +
  scale_colour_brewer(palette = "Dark2")+
  coord_flip() +
  facet_wrap(~Obs.Year) +
  theme_bw()
fig2

# Calculate the average julian day of passage means for each population of spring and summer run chinook

run_group <- srr_sub_mean %>%  # used to set hyper-parameters
        group_by(Group) %>%
        summarise(n = n(),
                  mu = mean(mean),
                  mu_sd = sd(mean),
                  sd_mu = mean(stddev,na.rm=TRUE),
                  sd_sd = sd(stddev,na.rm=TRUE)) %>%
            arrange(mu)

run_black <- srr_sub_mean %>%  # used to set hyper-parameters
  group_by(Species) %>%
  summarise(n = n(),
            mu = mean(mean),
            mu_sd = sd(mean),
            sd_mu = mean(stddev,na.rm=TRUE),
            sd_sd = sd(stddev,na.rm=TRUE)) %>%
  arrange(mu)
