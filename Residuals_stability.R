## calculating stability as residulas ###
# By Charlotte Kunze 05.10.2020

#### Load packages & import data ####
# 1. load required packages
library(lme4)
library(nlme)
library(lmerTest)
library(vegan)
library(tidyverse)

#### Response data ####

#import  data
data <- read.csv2('~/Desktop/MA/MA_Rcode/project_data/cnp_data_treatments.csv')
names(data)
str(data)
#######################DATA WRANGLING############################################

data1 <- data %>%
  mutate(sampling = str_remove(sampling, 'S'),
         sampling = as.numeric(sampling),
         MC = str_remove(MC, 'P'),
         day = 2*sampling) %>%
  mutate(fluctuation = ifelse(treatment == 'control', 0, ifelse(treatment == 'Fluctuating_12', 12, ifelse(treatment == 'Fluctuating_24', 24, ifelse(treatment == 'Fluctuating_36', 36, ifelse(treatment == 'Fluctuating_48', 48,6 )))))) %>%
  select(-X) %>%
  mutate(dayname = as.factor(day)) %>%
  mutate(interval = 48/fluctuation)
data1$interval[!is.finite(data1$interval)] <- 0 
str(data1)

#first look at data
ggplot(data1, aes(x = sampling, y=c_umol_l, col = treatment, group = MC ))+
  geom_point()+
  geom_line()+
  scale_x_continuous(limits = c(0,18), breaks = seq(0,18,2))+
  scale_y_continuous(limits = c(100,300), breaks = seq(100,300,20))+
  theme_classic()

#look at data
hist(data1$c_umol_l)
qqnorm(data1$c_umol_l)
qqline(data1$c_umol_l)


# lmer for carbon data

# Response variable: Carbon mg (measured over time)
# Explanatory variables: fluctuation interval in 48 h + day (time)
# Random: MC number
#         dayname (categorical variable)

mod1 <- lmer((c_umol_l) ~ interval*day + (1|MC) + (1|dayname), data=data1)
summary(mod1)
anova(mod1)

library(plyr)
res_PC_all <- resid(mod1)
df_res <- as.data.frame(res_PC_all)
df_residuals <- cbind(data1, df_res)

PC_resVar_CT<-ddply(df_residuals, c('interval'), summarise,
                    N = length(res_PC_all),
                    mean = mean(res_PC_all),
                    sd = sd(res_PC_all),
                    se = sd/sqrt(N),
                    var=var(res_PC_all))


#Models for residual variance
PC_var = with(PC_resVar_CT,glm(var~interval,family=gaussian))
summary(PC_var)
summary(aov(PC_var))

####################################################################
#### Stoichiometry ####
# import data #
Mastertable_fluctron <- read_delim("~/Desktop/MA/MA_Rcode/project_data/Mastertable_fluctron.csv", 
                                   ";", escape_double = FALSE, locale = locale(decimal_mark = ","), 
                                   trim_ws = TRUE)
#View(Mastertable_fluctron)
str(Mastertable_fluctron)

# import treatment information


### subdataframe with nutrient informations only ####

nutrients_master <- dplyr::select(Mastertable_fluctron, fluctuation, sampling, planktotron, carbon_umol_l, nitrate_umol_l, 
                                  'diss_Nitrat+Nitrit_umol_l', diss_Silikat_umol_l, diss_Phosphat_umol_l,
                                  POP_micromol_l, SiP_micromol_l)
names(nutrients_master)



#### RUE ####
RUE <- nutrients_master %>%
  dplyr::select(fluctuation, planktotron, sampling, carbon_umol_l, nitrate_umol_l, POP_micromol_l, 'diss_Nitrat+Nitrit_umol_l', 'diss_Phosphat_umol_l') %>%
  dplyr::rename(diss_P = 'diss_Phosphat_umol_l',
                diss_N = 'diss_Nitrat+Nitrit_umol_l') %>%
  mutate(TP = diss_P + POP_micromol_l,
         TN = diss_N + nitrate_umol_l) %>%
  drop_na(carbon_umol_l) %>%
  mutate(P_RUE = carbon_umol_l/ TP,
         N_RUE = carbon_umol_l/ TN) %>%
  dplyr::rename( MC = planktotron) %>%
  mutate(day = 2* sampling, 
         dayname = as.factor(day),
         MC = as.character(MC),
         fluctuation = as.numeric(fluctuation)) %>%
  mutate(interval = 48/fluctuation) #%>%
#drop_na(SiP)
RUE$interval[!is.finite(RUE$interval)] <- 0

#### Model
mod_N <- lmer(N_RUE ~ interval*day + (1|MC) + (1|dayname), data=RUE)
summary(mod_N)
anova(mod_N) #interval:day 36.817  36.817     1 98.000 23.1380 5.461e-06 ***##

res_RUE_N <- resid(mod_N)
df_resN <- as.data.frame(res_RUE_N)
df_residualsN <- cbind(RUE, df_resN)

PC_resVar_N<-ddply(df_residualsN, c('interval'), summarise,
                   N = length(res_RUE_N),
                   mean = mean(res_RUE_N),
                   sd = sd(res_RUE_N),
                   se = sd/sqrt(N),
                   var=var(res_RUE_N))


#Models for residual variance
PC_var = with(PC_resVar_N,glm(var~interval,family=gaussian))
summary(PC_var)
summary(aov(PC_var))



#### rue p
mod_P <- lmer(P_RUE ~ interval*day + (1|MC) + (1|dayname), data=RUE)
summary(mod_P)
anova(mod_P) #interval:day 5598.4  5598.4     1 97.999  2.8345 0.09544 .


res_RUE_P <- resid(mod_P)
df_resP <- as.data.frame(res_RUE_P)
df_residualsP <- cbind(RUE, df_resP)

PC_resVar_P<-ddply(df_residualsP, c('interval'), summarise,
                    N = length(res_RUE_P),
                    mean = mean(res_RUE_P),
                    sd = sd(res_RUE_P),
                    se = sd/sqrt(N),
                    var=var(res_RUE_P))


#Models for residual variance
PC_var = with(PC_resVar_P,glm(var~interval,family=gaussian))
summary(PC_var)
summary(aov(PC_var))

