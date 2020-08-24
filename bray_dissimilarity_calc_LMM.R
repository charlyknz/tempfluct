# scatter plot with bray-curits diss and LRR of function (carbon)
# By Charlotte Kunze 06.08.2020

#### Load packages & import data ####
library(vegan)

## import data ###
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


#### species contribution to BioVolume####

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

#copy sampling 0 data to replicates
data0 <- filter(all_data,sampling == 0)
data0$MC[data0$MC == 4] <-1
data0$MC[data0$MC == 11] <-12
data0$MC[data0$MC == 8] <-9
data0$MC[data0$MC == 3] <-5
data0$MC[data0$MC == 6] <-7
data0$MC[data0$MC == 10] <-2

#new df to calculate mean biovolume and data wrangling
pca_B <- all_data%>%
  bind_rows(., data0) %>%
  group_by(sampling, treatment_id, MC) %>%
  mutate(sum = sum(volume, na.rm  =T ),
         rel = volume/sum)%>%
  group_by(treatment_id, MC, sampling, species) %>%
  summarise(mean_V = mean(rel, na.rm = T),
            sd_V = sd(rel, na.rm = T),
            se_V = sd_V/sqrt(n())) %>% #calculate mean cells with biovolume
  dplyr::select(-sd_V) %>%
  drop_na(mean_V) %>%
  ungroup()%>%
  mutate(dummy = paste(treatment_id, MC, sampling, sep = ' ')) %>%
  dplyr::select(-treatment_id, -MC, -sampling, -se_V)%>%
  spread(key = species, value = mean_V) %>% #bring data in a wide format
  column_to_rownames('dummy')

pca_B[is.na(pca_B)] <- 0

#calculate bray distance using vegdist 
data.dist <- vegdist(pca_B, "bray") %>%
  broom::tidy() %>% #change to df
  separate(item1, into = c('treatment', 'MC', 'sampling'), ' ') %>%
  separate(item2, into = c('treatment2','MC2', 'sampling2'), ' ') %>%
  filter(sampling2 == 0)%>%
  filter(treatment == treatment2 & MC == MC2)%>% #extracts only treatment observations
  separate(treatment, into = c('Fluctuating', 'Fluctuation'), '_') %>%
  dplyr::select(-Fluctuating) %>%
  mutate(sampling = as.numeric(sampling)) %>%
  select(-MC2,-sampling2)

data.dist$Fluctuation[is.na(data.dist$Fluctuation)] <-0


#### STAS LMM for bray distance ####
# form Linear Model to answer Hypothesis 3: compositional turnover is more pronounec in quicker fluctuating

# Helmuts Model
# 1. load required packages
library(lme4)
library(nlme)
library(tidyverse)
library(lmerTest)

# data wrangling to add interval and day information
names(data.dist)
str(data.dist)
dist_m <- data.dist %>%
  mutate(day = 2*sampling) %>%
  mutate(dayname = as.factor(day)) %>%
  mutate(Fluctuation = as.numeric(Fluctuation),
         interval = 48/Fluctuation)
dist_m$interval[!is.finite(dist_m$interval)] <- 0 
str(dist_m)

bray1 <- lmer(distance ~ interval*day + (1|MC) + (1|dayname), data=dist_m)
summary(bray1)
anova(bray1)

par(mfrow=c(2,2),cex.axis=1.2, cex.lab=1.5)
plot(resid(bray1), ylab="residuales")
hist(resid(bray1), ylab="frecuencia",xlab="residuales", main="")
plot(fitted(bray1),resid(bray1),ylab="residuales")
qqnorm(resid(bray1), main=""); 
qqline(resid(bray1))




################################################################################
#### import LRRs for function data ####
stability <- read_delim("~/Desktop/MA/MA_Rcode/project_data/stability_MA.csv", 
                        ";", escape_double = FALSE, locale = locale(decimal_mark = ","), 
                        trim_ws = TRUE) %>%
  mutate(Fluctuation = str_remove(caseID, 'F')) %>%
  mutate(sampling = DAY/2) %>%
  dplyr::select(Fluctuation, sampling, LRR, var.lrr)

#calculate mean distance for each treatment
data.dist <- data.dist%>%
  group_by(Fluctuation, sampling, treatment2) %>%
  summarise(mean = mean(distance, na.rm = T))

dissimi_data <- left_join(stability, data.dist, by = c('Fluctuation', 'sampling')) %>%
  drop_na(mean)
dissimi_data$Fluctuation <- factor(as.factor(dissimi_data$Fluctuation),levels=c("48", "36", '24', '12', '6'))

ggplot(dissimi_data, aes(x = mean, y = LRR, fill = Fluctuation))+
  geom_point(aes(alpha = as.factor(sampling)), size = 3, pch = 21, col = 'black')+
  geom_hline(yintercept = 0)+
  scale_y_continuous(limits = c(-0.5, 0.5), breaks = c(-0.4,-0.2,0,0.2,0.4))+
  #scale_x_continuous(limits = c(0.1, 0.5), breaks = c(0,0.2,0.4))+
  scale_fill_manual(values = c('#0868ac','#41b6c4','#31a354','#addd8e','#fed976'))+
  labs( x = 'compositional turnover (bray-curtis)', y = 'LRR of POC', title = 'correlation funct and comp changes relative to T0', alpha = 'Sampling')+
  theme( panel.background = element_rect(fill = NA), #loescht den Hintergrund meines Plots/ fuellt ihn mit nichts
         #panel.grid.major.y = element_line(color='grey', linetype = 'dashed', size=0.2),
         panel.border= element_rect(colour = "black", fill=NA, size=0.5),
         strip.background = element_rect(color ='black', fill = 'white'),
         strip.text = element_text(face = 'italic'),
         legend.background = element_blank(),
         legend.key = element_blank(),
         text = element_text(size=12))
#ggsave(plot = last_plot(), file = 'dissimilarity plot_txt0.png')

