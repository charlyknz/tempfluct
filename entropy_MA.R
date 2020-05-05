#### Script to calculate entropy in thesis data ####
## Charlotte Kunze

## orginial script at Datadryad: Burton et al. 2019 https://orcid.org/0000-0002-0215-0227 
##

##load required packages
library(Rmisc)
library(tidyverse)
library(AICcmodavg)
library(MuMIn)
library(lme4)
library(TSEntropies)


# import temperature data 
temp_data <- read.csv('~/Desktop/MA/MA_Rcode/project_data/all_temperatures.csv')
str(temp_data) #check import

treatments <- read.csv2('~/Desktop/MA/MA_Rcode/project_data/treatments_units.csv')
str(treatments)
treatments <- treatments %>%
  arrange(fluctuation) %>%
  mutate(Treatment =  fct_recode(as.character(fluctuation), "5" = "6",
                                            "1" = "12",
                                            "2" = "24",
                                            "3" = "36",
                                            "4" = "48"))
treatments$Treatment = as.numeric(treatments$Treatment)

temp <- left_join(temp_data, treatments, by = c('unit')) %>%
  select(unit, fluctuation, Treatment, timestamp, actual_tempmiddle) %>%
  mutate(datetime = timestamp) %>%
  separate(timestamp, into = c('date', 'time'), ' ') %>%
  mutate(HMS = time) %>%
  separate(time, into = c('H', 'M', 'S'), ':')%>%
  filter(as.character(datetime) > '2019-08-27 07:00:00') %>%
filter(as.character(datetime) < '2019-09-17 14:00:00') 

#### aggregate main dataframe by treatment
Summarytemp <- Rmisc::summarySE(temp, measurevar="actual_tempmiddle", groupvars=c("Treatment"), na.rm=TRUE)


#### calculate mean temperature etc per hr per treatment
mean <- temp %>%
  group_by(fluctuation, Treatment, date, H) %>%
  summarise(MeanT = mean(actual_tempmiddle, na.rm = T)) %>%
  arrange(fluctuation,date, H) 
str(mean)
#### entropy as a measure of unpredictability
#### lowest in predictable treatment and higher in UP treatments
entropy <- rep(NA,length(Summarytemp$Treatment))
Treatment <- entropy

#### for loop to calculate sample entropy for each of the 17 fluctuating treatments
j <- 1
for(t in 1:6){# for each of the 17 treatments
  temp1 <- subset(mean, mean$Treatment == t)
  j <- j + 1
  entropy[j] <- SampEn(temp1$MeanT, dim = 24, lag = 1, r = 0.2 * sd(temp1$MeanT))
  Treatment[j] <- mean(temp1$Treatment)}

entropy_data <- as.data.frame(cbind(entropy, Treatment))
entropy_data<- entropy_data[-1,]

entropy_MA <- left_join(entropy_data, treatments, by = c('Treatment'))
plot_entropy <- ggplot(entropy_MA, aes(x = fluctuation, y = entropy, col = as.factor(treatment) ))+
  geom_point()+
  scale_x_continuous(limits = c(0,48), breaks = c(0,6,12,24,36,48))+
  theme_classic()
plot_entropy
#ggsave(plot = plot_entropy, file = 'entropy_MA.png')
##################################################################
entropy_MA$entropy[is.na(entropy_MA$entropy)] <- 0

entropy_r <- left_join(entropy_MA, all_logistic_results, by = c('treatment'))
ggplot(entropy_r, aes(x = entropy, y = mumax, col = treatment))+
  geom_point()+
  theme_classic()
#ggsave(plot = last_plot(), file = 'no_predictability_effect.png')
