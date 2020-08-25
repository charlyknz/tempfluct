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

mod1 <- lmer(c_umol_l ~ interval*day + (1|MC) + (1|dayname), data=data1)
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
####################################################################
Mastertable <- read_delim("~/Desktop/MA/MA_Rcode/project_data/Mastertable.csv", 
                          ";", escape_double = FALSE, locale = locale(decimal_mark = ","), 
                          trim_ws = TRUE) %>%
  rename(MC = planktotron)
str(Mastertable)

## merge with treatment information
treatments <- read.csv2('~/Desktop/MA/MA_Rcode/project_data/treatments_units.csv')
str(treatments)
names(treatments) = c('MC', 'fluctuation', 'treatment')
master_data <- left_join(Mastertable, treatments, by = c('MC')) 


master <- master_data %>%
  select(MC, sampling, treatment, fluctuation, C_Zoo_µmol_l ) %>%
  drop_na(C_Zoo_µmol_l) %>%
  mutate(day = 2* sampling, 
         dayname = as.factor(day)) %>%
  mutate(interval = 48/fluctuation)
master$interval[!is.finite(master$interval)] <- 0
  
# Response variable: Carbon mg (measured over time)
# Explanatory variables: fluctuation interval in 48 h + day (time)
# Random: MC number
#         dayname (categorical variable)

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


########################COMPOSITION LMM##########################################

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
data.dist <- vegdist(data_comp, "bray") %>% #distance calculation
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

bray1 <- lmer(distance ~ interval*day + (1|MC) + (1|dayname), data=dist_m)
summary(bray1)
anova(bray1)

#plot residuals
par(mfrow=c(2,2),cex.axis=1.2, cex.lab=1.5)
plot(resid(bray1), ylab="residuales")
hist(resid(bray1), ylab="frecuencia",xlab="residuales", main="")
plot(fitted(bray1),resid(bray1),ylab="residuales")
qqnorm(resid(bray1), main=""); 
qqline(resid(bray1))

###################################################################################

####PCA on pigment data only ####
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

