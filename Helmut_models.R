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
###################################################################################

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


###################################################################################
###################################################################################

#### composition data ####
counts <- read_delim("~/Desktop/MA/MA_Rcode/project_data/phyto_new.csv", 
                     ";", escape_double = FALSE, col_types = cols(date = col_character()), 
                     locale = locale(decimal_mark = ","), 
                     trim_ws = TRUE) %>%
  drop_na(MC)

str(counts) #check import
names(counts)

#metadata
df <- data.frame( MC = c('1', '4', '2', '10', '3', '5', '6', '7', '8', '9', '11', '12'),
                  treatment_id = c('control', 'control', 'Fluctuating_36','Fluctuating_36',
                                   'Fluctuating_6', 'Fluctuating_6', 'Fluctuating_48', 'Fluctuating_48',
                                   'Fluctuating_24', 'Fluctuating_24', 'Fluctuating_12','Fluctuating_12') )

#change variable format to numeric 
counts$cells_ml <- as.numeric(gsub(",", ".", counts$cells_ml))
counts$MC = as.character(counts$MC)
counts$cells_ml[is.na(counts$cells_ml)] <-0


###################################################################################
#### species contribution to BioVolume ####

#### Biovolume ####
algal_measurement <- read_delim("~/Desktop/MA/MA_Rcode/project_data/algal_measurement.csv", 
                                ";", escape_double = FALSE, locale = locale(decimal_mark = ","), 
                                trim_ws = TRUE)
str(algal_measurement)

#calculate mean biovolume (mean of measurements)
BioV <- algal_measurement %>%
  group_by(genus, species, id) %>%
  summarise(mean_BioV = mean(BioV, na.rm = T),
            sd = sd(BioV, na.rm = T)) %>%
  ungroup() %>%
  mutate(species = paste(genus, species, sep = ' ')) %>%
  dplyr::select(-genus)


#### merge all_data ####
all_data <- left_join(counts, df, by = c('MC')) %>%
  left_join(., BioV, by = c('species') ) %>%
  dplyr::select(treatment_id, MC, sampling, species, phylum, cells_ml,mean_BioV) %>%
  mutate(volume = cells_ml*mean_BioV)  #calculate biovolume

###################################################################################

#copy sampling 0 data to replicate MC
data0 <- filter(all_data,sampling == 0)
data0$MC[data0$MC == 4] <-1
data0$MC[data0$MC == 11] <-12
data0$MC[data0$MC == 8] <-9
data0$MC[data0$MC == 3] <-5
data0$MC[data0$MC == 6] <-7
data0$MC[data0$MC == 10] <-2

#new df to calculate mean biovolume and data wrangling
pca_B <- all_data%>%
  bind_rows(., data0) %>% #join data 
  group_by(sampling, treatment_id, MC) %>%
  mutate(sum = sum(volume, na.rm  =T ), #calculate relative contribution
         rel = volume/sum)%>%
  group_by(treatment_id, MC, sampling, species) %>%
  summarise(mean_V = mean(rel, na.rm = T), #calculate mean rel contr
            sd_V = sd(rel, na.rm = T),
            se_V = sd_V/sqrt(n())) %>% #calculate mean cells with biovolume
  dplyr::select(-sd_V) %>%
  drop_na(mean_V) %>%
  ungroup()%>%
  mutate(dummy = paste(treatment_id, MC, sampling, sep = ' ')) %>% #create dummy variable containing meta info
  dplyr::select(-treatment_id, -MC, -sampling, -se_V)%>% #remove non-numeric variables
  spread(key = species, value = mean_V) %>% #bring data in a wide format
  column_to_rownames('dummy') 

pca_B[is.na(pca_B)] <- 0 #substitute NAs with 0

#calculate bray distance using vegdist 
data.dist <- vegdist(pca_B, "bray") %>% #distance calculation
  broom::tidy() %>% #change to df
  separate(item1, into = c('treatment', 'MC', 'sampling'), ' ') %>% #separate item1,2 into variable columns
  separate(item2, into = c('treatment2','MC2', 'sampling2'), ' ') %>%
  filter(sampling2 == 0)%>% #only start values
  filter(treatment == treatment2 & MC == MC2)%>% #extracts only treatment & MC observations
  separate(treatment, into = c('Fluctuating', 'Fluctuation'), '_') %>%
  dplyr::select(-Fluctuating) %>%
  mutate(sampling = as.numeric(sampling)) %>%
  select(-MC2,-sampling2)

data.dist$Fluctuation[is.na(data.dist$Fluctuation)] <-0

###################################################################################

#### LMM for bray distance ####
# form Linear Model to answer Hypothesis 3: 
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

#model using bray distance for each MC 
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


