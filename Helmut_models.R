#### Models using lmer ####

# By Charlotte Kunze 06.08.2020

#### Load packages & import data ####
# 1. load required packages
library(lme4)
library(nlme)
library(lmerTest)
library(vegan)
library(tidyverse)

#### Response data ####

#import  data

Mastertable_fluctron <- read_delim("~/Desktop/MA/MA_Rcode/project_data/Mastertable_fluctron.csv", 
                                   ";", escape_double = FALSE, locale = locale(decimal_mark = ","), 
                                   trim_ws = TRUE)
#View(Mastertable_fluctron)
str(Mastertable_fluctron)

# use only carbon data
data <- Mastertable_fluctron %>%
  dplyr::select(fluctuation, sampling, planktotron, carbon_umol_l)
names(data)
#######################DATA WRANGLING############################################

data1 <- data %>%
  rename(MC = planktotron)%>%
  mutate(sampling = as.numeric(sampling),
         day = 2*sampling) %>%
  mutate(dayname = as.factor(day)) %>%
  mutate(interval = 48/fluctuation)
data1$interval[!is.finite(data1$interval)] <- 0 
str(data1)

#first look at data
ggplot(data1, aes(x = sampling, y=carbon_umol_l, col = fluctuation, group = MC ))+
  geom_point()+
  geom_line()+
  scale_x_continuous(limits = c(0,18), breaks = seq(0,18,2))+
  scale_y_continuous(limits = c(100,300), breaks = seq(100,300,20))+
  theme_classic()

#look at data
hist(data1$carbon_umol_l)
qqnorm(data1$carbon_umol_l)
qqline(data1$carbon_umol_l)


# lmer for carbon data

# Response variable: Carbon mg (measured over time)
# Explanatory variables: fluctuation interval in 48 h + day (time)
# Random: MC number
#         dayname (categorical variable)

mod1 <- lmer((carbon_umol_l) ~ interval*day + (1|MC) + (1|dayname), data=data1)
summary(mod1)
anova(mod1)

#plot residuals
par(mfrow=c(2,2),cex.axis=1.2, cex.lab=1.5)
plot(resid(mod1), ylab="residuales")
hist(resid(mod1), ylab="frecuencia",xlab="residuales", main="")
plot(fitted(mod1),resid(mod1),ylab="residuales")
qqnorm(resid(mod1), main=""); 
qqline(resid(mod1))


####################################################################
#######################Zooplankton###########################
Mastertable <- read_delim("~/Desktop/MA/MA_Rcode/project_data/Mastertable_fluctron.csv", 
                          ";", escape_double = FALSE, locale = locale(decimal_mark = ","), 
                          trim_ws = TRUE) %>%
  rename(MC = planktotron)
str(Mastertable)

## merge with treatment information

master <- Mastertable %>%
  dplyr::select(MC, sampling, fluctuation, C_Zoo_µmol_l ) %>%
  drop_na(C_Zoo_µmol_l) %>%
  mutate(day = 2* sampling, 
         dayname = as.factor(day),
         MC = as.character(MC)) %>%
  mutate(interval = 48/fluctuation)
master$interval[!is.finite(master$interval)] <- 0
  
# Response variable: Carbon mg (measured over time)
# Explanatory variables: fluctuation interval in 48 h + day (time)
# Random: MC number
#         dayname (categorical variable)
str(master)

zoo1 <- lmer(C_Zoo_µmol_l ~ interval*day + (1|MC) + (1|dayname), data=master)
summary(zoo1)
anova(zoo1)

#plot residuals
par(mfrow=c(2,2),cex.axis=1.2, cex.lab=1.5)
plot(resid(zoo1), ylab="residuales")
hist(resid(zoo1), ylab="frecuencia",xlab="residuales", main="")
plot(fitted(zoo1),resid(zoo1),ylab="residuales")
qqnorm(resid(zoo1), main=""); 
qqline(resid(zoo1))


########################Phytoplankton COMPOSITION LMM##########################################

#### composition data ####
CountsBV <-read_delim("~/Desktop/MA/MA_Rcode/project_data/CountsBV.csv", 
                                  ";", escape_double = FALSE, locale = locale(decimal_mark = ","), 
                                  trim_ws = TRUE)
names(CountsBV) 
str(CountsBV) 
  
# Data wrangling

data_comp <- CountsBV %>%
  dplyr::select(-X1) %>%
  drop_na(rel) %>%
  ungroup()%>%
  mutate(dummy = paste(treatment_id, MC, sampling, sep = ' ')) %>% #create dummy variable containing meta info
  dplyr::select(species, rel, dummy)%>% #remove non-numeric variables
  spread(key = species, value = rel) %>% #bring data in a wide format
  column_to_rownames('dummy') 

data_comp[is.na(data_comp)] <- 0 #substitute NAs with 0

#calculate bray distance using vegdist 
data.dist <- vegdist(data_comp, 'bray') %>% #distance calculation
  broom::tidy() %>% #change to df
  separate(item1, into = c('treatment', 'MC', 'sampling'), ' ') %>% #separate item1,2 into variable columns
  separate(item2, into = c('treatment2','MC2', 'sampling2'), ' ') %>%
  filter(sampling2 == 0)%>% #only start values to compare with
  filter(treatment == treatment2 & MC == MC2)%>% #extracts only matching treatment & MC observations
  separate(treatment, into = c('Fluctuating', 'Fluctuation'), '_') %>%
  dplyr::select(-Fluctuating) %>%
  mutate(sampling = as.numeric(sampling)) %>%
  select(-MC2,-sampling2)

data.dist$Fluctuation[is.na(data.dist$Fluctuation)] <-0

###################################################################################

#### LMM for bray distance ####
# from Linear Model to answer Hypothesis 3: 
# compositional turnover is more pronounced in quicker fluctuating treatments

#check data
names(data.dist)
str(data.dist)

# data wrangling to add interval and day information
dist_m <- data.dist %>%
  mutate(day = 2*sampling) %>%
  mutate(dayname = as.factor(day)) %>%
  mutate(Fluctuation = as.numeric(Fluctuation),
         interval = 48/Fluctuation)
dist_m$interval[!is.finite(dist_m$interval)] <- 0 
str(dist_m)

#model 
# Response variable: bray distance for each MC compared to T0
# Explanatory variables: interval (fluctuation interval in 48 h) + time
# Random: MC number
#         dayname (categorical variable)

bray1 <- lmer(log(distance) ~ interval*day + (1|MC) + (1|dayname), data=dist_m)
summary(bray1)
anova(bray1)


#plot residuals
par(mfrow=c(2,2),cex.axis=1.2, cex.lab=1.5)
plot(resid(bray1), ylab="residuales")
hist(resid(bray1), ylab="frecuencia",xlab="residuales", main="")
plot(fitted(bray1),resid(bray1),ylab="residuales")
qqnorm(resid(bray1), main=""); 
qqline(resid(bray1))


dist_m %>%
  group_by(Fluctuation, day) %>%
  summarise(mean = mean(distance, na.rm = T)) %>%
ggplot(., aes( x = day, y = mean, color = as.factor(Fluctuation), group = Fluctuation))+
  geom_point(size = 3)+
  geom_line()+
  theme_classic()
###################################################################################

####pigment data ####
pigment <- read.csv2("~/Desktop/MA/MA_Rcode/project_data/master_pigments.csv", sep = ";", dec = ',')
str(pigment)

#bring data in form
pigment_data <- pigment %>%
  gather(key = 'pigments',value = 'value',-X, -no, -date,-treatment, -sampling, -planktotron )%>%
  group_by(treatment,sampling, planktotron)%>%
  mutate(sum = sum(value)) %>%
  ungroup()%>%
  spread(pigments, value ) %>%
  mutate(rel_allo = Allo/sum)%>%
  mutate(rel_anth = Anth/sum)%>% ###berechnung der relativen werte in neuer spalte nach Schlueter et al.
  mutate(rel_bb.Car = bb.Car/sum)%>%
  mutate(rel_cNeo = c.Neo/sum) %>%
  mutate(rel_cantha = Cantha/sum)%>%
  mutate(rel_chl.a = Chl.a/sum)%>%
  mutate(rel_chl.b = Chl.b/sum)%>%
  mutate(rel_chl.c1 = Chl.c1/sum)%>%
  mutate(rel_chl.c2 = Chl.c2/sum)%>%
  mutate(rel_diad = Diadino/sum)%>%
  mutate(rel_diato = Dino/sum)%>%
  mutate(rel_dino = Diato/sum)%>%
  mutate(rel_echin = Echin/sum)%>%
  mutate(rel_fuco = Fuco/sum)%>%
  mutate(rel_lut = Lut/sum)%>%
  mutate(rel_myxo = Myxo/sum)%>%
  mutate(rel_peri = Peri/sum)%>%
  mutate(rel_phe.a = Phe.a/sum)%>%
  mutate(rel_viola = Viola/sum) %>%
  mutate(rel_zea = Zea/sum)%>%
  dplyr::select(-c("X","Allo", "Anth","bb.Car", "c.Neo", "Cantha","Chl.a","Chl.c1",
                   "Chl.c2","Cryp","Diadino", "Diato","Dino","Echin","Lut","Myxo","Peri",       
                   "Phe.a","Phe.b", "Viola", "Zea" , 'Fuco', 'Chl.b')) %>%
  gather(key = 'pigment', value = 'value', -planktotron, - sampling, -treatment,-sum, -no, -date) %>%
  mutate(fluctuation = paste(ifelse(treatment == 'control', 0, ifelse(treatment == 'Fluctuating_48', 48, ifelse(treatment == 'Fluctuating_36', 36,
                                                                                                                ifelse(treatment == 'Fluctuating_24', 24, ifelse(treatment == 'Fluctuating_12', 12, 6)))))))
pigment_data$fluctuation <- factor(as.factor(pigment_data$fluctuation),levels=c("0", "48", "36", '24', '12', '6'))
pigment_data$sampling = as.numeric(pigment_data$sampling)
pigment_data[is.na(pigment_data)] <-0 #make sure we don't have NAs

names(pigment_data)
# 1. decide which variables are active, remove character variables and turn them into row names 
pigment_data.active <- pigment_data %>% 
  ungroup()%>%
  dplyr::select(-treatment) %>% #remove columns which are not needed anymore
  mutate( id = paste(fluctuation, planktotron, sampling, sep = '_')) %>% #create dummy column
  dplyr::select(pigment, value, id) %>%
  spread(key = pigment, value = value) %>%
  column_to_rownames("id")%>% #explanatory variables as row names
  dplyr::select(-rel_bb.Car) #remove variable with only 0 entries

pigment_data.active[is.na(pigment_data.active)] <- 0 #substitute NAs with 0
  
  #calculate bray distance using vegdist 
pig.dist <- vegdist(pigment_data.active, "bray") %>% #distance calculation
    broom::tidy() %>% #change to df
    separate(item1, into = c('treatment', 'MC', 'sampling'), '_') %>% #separate item1,2 into variable columns
    separate(item2, into = c('treatment2','MC2', 'sampling2'), '_') %>%
    filter(sampling2 == 0)%>% #only start values to compare with
    filter(treatment == treatment2 & MC == MC2)%>% #extracts only matching treatment & MC observations
    mutate(sampling = as.numeric(sampling)) %>%
    select(-MC2,-sampling2) %>%
   mutate(day = 2*sampling,
          dayname = as.factor(day)) %>%
    mutate(treatment = as.numeric(treatment),
         interval = 48/treatment)
pig.dist$interval[!is.finite(pig.dist$interval)] <- 0 
str(pig.dist)
  
#check data
names(pig.dist)
str(pig.dist)

#model 
# Response variable: bray distance for each MC compared to T0
# Explanatory variables: interval (fluctuation interval in 48 h) + time
# Random: MC number
#         dayname (categorical variable)

bray1 <- lmer(distance ~ interval*day + (1|MC) + (1|dayname), data=pig.dist)
summary(bray1)
anova(bray1)

#plot residuals
par(mfrow=c(2,2),cex.axis=1.2, cex.lab=1.5)
plot(resid(bray1), ylab="residuales")
hist(resid(bray1), ylab="frecuencia",xlab="residuales", main="")
plot(fitted(bray1),resid(bray1),ylab="residuales")
qqnorm(resid(bray1), main=""); 
qqline(resid(bray1))
############################################################################################################
#### relative pigment diversity data ####

##
# R Skript to analyse pigment composition

##import master data
data <- read.csv2("~/Desktop/MA/MA_Rcode/project_data/master_pigments.csv", sep = ";", dec = ',')
names(data)

## data wrangling
pigData <- data %>%
  gather(key = 'pigments',value = 'value',-X, -no, -date,-treatment, -sampling, -planktotron )%>%
  group_by(treatment,sampling, planktotron)%>%
  mutate(sum = sum(value)) %>%
  ungroup()%>%
  spread(pigments, value ) %>%
  mutate(rel_allo = Allo/sum)%>%
  mutate(rel_anth = Anth/sum)%>% ###berechnung der relativen werte in neuer spalte nach Schlueter et al.
  mutate(rel_bb.Car = bb.Car/sum)%>%
  mutate(rel_cNeo = c.Neo/sum) %>%
  mutate(rel_cantha = Cantha/sum)%>%
  mutate(rel_chl.a = Chl.a/sum)%>%
  mutate(rel_chl.b = Chl.b/sum)%>%
  mutate(rel_chl.c1 = Chl.c1/sum)%>%
  mutate(rel_chl.c2 = Chl.c2/sum)%>%
  mutate(rel_diad = Diadino/sum)%>%
  mutate(rel_diato = Dino/sum)%>%
  mutate(rel_dino = Diato/sum)%>%
  mutate(rel_echin = Echin/sum)%>%
  mutate(rel_fuco = Fuco/sum)%>%
  mutate(rel_lut = Lut/sum)%>%
  mutate(rel_myxo = Myxo/sum)%>%
  mutate(rel_peri = Peri/sum)%>%
  mutate(rel_phe.a = Phe.a/sum)%>%
  mutate(rel_viola = Viola/sum) %>%
  mutate(rel_zea = Zea/sum)%>%
  dplyr::select(-c("X","Allo", "Anth","bb.Car", "c.Neo", "Cantha","Chl.a","Chl.c1",
                   "Chl.c2","Cryp","Diadino", "Diato","Dino","Echin","Lut","Myxo","Peri",       
                   "Phe.a","Phe.b", "Viola", "Zea" , 'Fuco', 'Chl.b')) %>%
  gather(key = 'pigment', value = 'value', -planktotron, - sampling, -treatment,-sum, -no, -date) %>%
  mutate(fluctuation = paste(ifelse(treatment == 'control', 0, ifelse(treatment == 'Fluctuating_48', 48, ifelse(treatment == 'Fluctuating_36', 36,
                                                                                                                ifelse(treatment == 'Fluctuating_24', 24, ifelse(treatment == 'Fluctuating_12', 12, 6)))))))
pigData$sampling = as.numeric(pigData$sampling)
pigData[is.na(pigData)] <-0 #make sure we don't have NAs

names(pigData)
# 1. decide which variables are active, remove character variables and turn them into row names 
pigment_data.active <- pigData %>% 
  ungroup()%>%
  dplyr::select(-treatment) %>% #remove columns which are not needed anymore
  spread(key = pigment, value = value) %>%
  dplyr::select(-rel_bb.Car, -date, -sum, -no) #remove variable with only 0 entries

pigment_data.active[is.na(pigment_data.active)] <- 0 #substitute NAs with 0


### calculate  diversity indices of pigment composition
diversity <- pigment_data.active 

#remove NAs and exchange with 0
diversity[is.na(diversity)] <- 0

##calculate shannon diversity index
diversity$shannon = diversity(diversity[, -c(1:4)], MARGIN = 1, index='shannon') #new column containing the calculated shannon index
diversity <- select(diversity, planktotron, sampling,fluctuation, shannon, everything())
diversity$simpson = diversity(diversity[, -c(1:5)], MARGIN = 1, index= 'invsimpson') #new column containing the calculated shannon index
diversity <- select(diversity, planktotron, sampling,fluctuation, shannon, simpson, everything())
## calculate species richness

ab_presence <- decostand(diversity[, -c(1:6)], method= 'pa', na.rm=T) #df giving absence/presence data using decostand function
diversity$rich = apply(ab_presence, MARGIN = 1, FUN = sum) #new column containing the sum of species present per side (by row = MARGIN = 1)

diversity$evenness = diversity$shannon/log(diversity$rich)

#### LMM on diversity ####
modell_dat <- diversity %>%
  dplyr::select(planktotron, sampling, fluctuation,evenness, simpson, rich) %>%
  mutate(day = sampling *2, #add day column
         dayname = as.factor(day)) %>%
  mutate(fluctuation = as.numeric(fluctuation),
         interval = 48/fluctuation) #calculate fluctuation interval in 48 h
modell_dat$interval[is.infinite(modell_dat$interval)] <-0
modell_dat$fluctuation[is.na(modell_dat$fluctuation)] <-0

div <- lmer(rich ~ interval*day + (1|planktotron) + (1|dayname), data=modell_dat)
summary(div)
anova(div)

############################################################################################################
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

############################################################################################################
#### calculate ratios ####

ratio <- nutrients_master %>%
  mutate(CN = carbon_umol_l/nitrate_umol_l,
         NP = nitrate_umol_l/POP_micromol_l,
         CP = carbon_umol_l/POP_micromol_l,
         CSi = carbon_umol_l/SiP_micromol_l,
         NSi = nitrate_umol_l/SiP_micromol_l,
         SiP = SiP_micromol_l/POP_micromol_l) %>%
  dplyr::select(nitrate_umol_l,fluctuation, sampling, planktotron,SiP_micromol_l,CN, NP, CP, CSi, NSi, NSi, SiP)%>%
  dplyr::rename( MC = planktotron) %>%
  mutate(day = 2* sampling, 
         dayname = as.factor(day),
         MC = as.character(MC)) %>%
  mutate(interval = 48/fluctuation) #%>%
  #drop_na(SiP)
ratio$interval[!is.finite(ratio$interval)] <- 0
  

#### Model
library(lmerTest)
CN1 <- lmer((CSi) ~ interval*day + (1|MC) + (1|dayname), data=ratio)
summary(CN1)
anova(CN1)

#output:

# CN interval:day 2.635e-05 ***
# CSi time **
# SiP time **
# NSi time **
# SiP_micromol_l time **
# P no sign
# NP no sign
# CP no sign
# NP no sign

#plot residuals
par(mfrow=c(2,2),cex.axis=1.2, cex.lab=1.5)
plot(resid(CN1), ylab="residuales")
hist(resid(CN1), ylab="frecuencia",xlab="residuales", main="")
plot(fitted(CN1),resid(CN1),ylab="residuales")
qqnorm(resid(CN1), main=""); 
qqline(resid(CN1))



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

mod_P <- lmer(P_RUE ~ interval*day + (1|MC) + (1|dayname), data=RUE)
summary(mod_P)
anova(mod_P) #interval:day 5598.4  5598.4     1 97.999  2.8345 0.09544 .


