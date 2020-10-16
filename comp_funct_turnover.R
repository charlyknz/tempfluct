### Functional and compositional stability ###
library(vegan)
library(tidyverse)
library(ggpubr)
library(cowplot)

### 1. biomass and  compositional turnover ###

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
         fluctuation = treatment,
         MC = as.numeric(MC)) %>%
  dplyr::select(MC, fluctuation, day, sampling, distance)
str(pig.dist)
pig.dist <- as.data.frame(pig.dist)


####1. Biomass data ####
Mastertable_fluctron <- read_delim("~/Desktop/MA/MA_Rcode/project_data/Mastertable_fluctron.csv", 
                                   ";", escape_double = FALSE, locale = locale(decimal_mark = ","), 
                                   trim_ws = TRUE)
#View(Mastertable_fluctron)
str(Mastertable_fluctron)

# use only carbon data
cnp_data <- Mastertable_fluctron %>%
  dplyr::select(fluctuation, sampling, planktotron, carbon_umol_l)
names(cnp_data)

#carbon accumulated
co_data <- cnp_data %>%
  ungroup()%>%
  rename(MC= planktotron) %>%
  mutate(sampling = as.numeric(sampling),
         day = sampling *2) %>%
  ungroup() %>%
  mutate(fluctuation = as.numeric(fluctuation)) %>%
  dplyr::select(MC, fluctuation, day, sampling, carbon_umol_l)

all <- right_join(co_data, pig.dist, by = c('sampling', 'fluctuation', 'MC', 'day')) %>%
  dplyr::select(sampling, fluctuation, MC,distance, carbon_umol_l ) %>%
  drop_na(carbon_umol_l)

formula <- y ~ x

all$fluctuation <- factor(as.factor(all$fluctuation),levels=c("0", "48", "36", '24', '12', '6'))
BioComp<-ggplot(all, aes(x = distance, y = carbon_umol_l,group = fluctuation, fill = as.factor(fluctuation)))+
  geom_point(pch = 21, size = 3)+
  #geom_errorbar(aes(ymin = lower.ci.mpg, ymax = upper.ci.mpg), color = '#404040', width = .2, size = .5)+
  #stat_smooth( method = "lm", se = T, size = 0.5,formula = formula, color = 'black') +
 # stat_regline_equation(formula = formula, 
#                        aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~"))) +     
  scale_fill_manual(values = c( '#000000','#0868ac','#41b6c4','#31a354','#addd8e','#fed976'))+
  labs(x = 'Compositional turnover', y = '', fill = 'Fluctuation  \nfrequency [h]')+
  theme( panel.background = element_rect(fill = NA), #loescht den Hintergrund meines Plots/ fuellt ihn mit nichts
         #panel.grid.major.y = element_line(color='grey', linetype = 'dashed', size=0.2),
         panel.border= element_rect(colour = "black", fill=NA, size=0.5),
         strip.background = element_rect(color ='black', fill = 'white'),
         strip.text = element_text(face = 'italic'),
         legend.background = element_blank(),
         legend.position  ='none',
         legend.key = element_blank(),
         text = element_text(size=17))
BioComp
#ggsave(plot = BioComp, file = 'comp_funct_turnover_biom.png', width = 5, height = 4)

ggscatter(all, x = "distance", y = "carbon_umol_l", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman", xlab = 'compositional turnover',
          ylab = 'Phytoplankton biomass')


#####################################################################################

####2. diversity of pigments and carbon ####

##import master data
data <- read.csv2("~/Desktop/MA/MA_Rcode/project_data/master_pigments.csv", sep = ";", dec = ',')
names(data)

#### diversity for relative data ####
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

names(diversity)

#### Eveness of species phytoplankton biomass relationship####
data_comp  <- diversity %>%
  select(sampling, fluctuation,planktotron,shannon, evenness, rich, simpson) %>%
  gather(key = 'index', value = 'value', -fluctuation, -sampling, -planktotron) %>%
  mutate(day = sampling *2,
         fluctuation= as.numeric(fluctuation)) %>%
  rename(MC = planktotron)
str(data_comp)  
#######

####2. Biomass data ####
Mastertable_fluctron <- read_delim("~/Desktop/MA/MA_Rcode/project_data/Mastertable_fluctron.csv", 
                                   ";", escape_double = FALSE, locale = locale(decimal_mark = ","), 
                                   trim_ws = TRUE)
#View(Mastertable_fluctron)
str(Mastertable_fluctron)

# use only carbon data
cnp_data <- Mastertable_fluctron %>%
  dplyr::select(fluctuation, sampling, planktotron, carbon_umol_l)
names(cnp_data)

#carbon accumulated
co_data <- cnp_data %>%
  ungroup()%>%
  rename(MC= planktotron) %>%
  mutate(sampling = as.numeric(sampling),
         day = sampling *2) %>%
  ungroup() %>%
  mutate(fluctuation = as.numeric(fluctuation)) %>%
  dplyr::select(MC, fluctuation, day, sampling, carbon_umol_l)
str(co_data)


## melt biomass and compositionap data together 
all_rich <- right_join(co_data, data_comp, by = c('sampling', 'fluctuation', 'MC', 'day')) %>%
  dplyr::select(day, sampling,fluctuation, MC,value,index, carbon_umol_l ) %>%
  drop_na(carbon_umol_l)

formula <- y ~ x

all_rich$fluctuation <- factor(as.factor(all_rich$fluctuation),levels=c("0", "48", "36", '24', '12', '6'))
BioRich <- ggplot(subset(all_rich, index == 'rich'), aes(x = value, y = carbon_umol_l,group = fluctuation, fill = as.factor(fluctuation)))+
  geom_point(pch = 21, size = 3)+
  #geom_errorbar(aes(ymin = lower.ci.mpg, ymax = upper.ci.mpg), color = '#404040', width = .2, size = .5)+
  #stat_smooth( method = "lm", se = T, size = 0.5,formula = formula, color = 'black') +
  # stat_regline_equation(formula = formula, 
  #                        aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~"))) +     
  scale_fill_manual(values = c( '#000000','#0868ac','#41b6c4','#31a354','#addd8e','#fed976'))+
  scale_x_continuous(limits = c(5.5, 18), breaks = seq(6,18,3))+
  labs(x = 'Pigment richness', y = '', fill = 'Fluctuation  \nfrequency [h]')+
  theme( panel.background = element_rect(fill = NA), #loescht den Hintergrund meines Plots/ fuellt ihn mit nichts
         #panel.grid.major.y = element_line(color='grey', linetype = 'dashed', size=0.2),
         panel.border= element_rect(colour = "black", fill=NA, size=0.5),
         strip.background = element_rect(color ='black', fill = 'white'),
         strip.text = element_text(face = 'italic'),
         legend.background = element_blank(),
         legend.position  ='none',
         legend.key = element_blank(),
         text = element_text(size=17))
BioRich

BioEven <- ggplot(subset(all_rich, index == 'evenness'), aes(x = value, y = carbon_umol_l,group = fluctuation, fill = as.factor(fluctuation)))+
  geom_point(pch = 21, size = 3)+
  #geom_errorbar(aes(ymin = lower.ci.mpg, ymax = upper.ci.mpg), color = '#404040', width = .2, size = .5)+
  #stat_smooth( method = "lm", se = T, size = 0.5,formula = formula, color = 'black') +
  # stat_regline_equation(formula = formula, 
  #                        aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~"))) +     
  scale_fill_manual(values = c( '#000000','#0868ac','#41b6c4','#31a354','#addd8e','#fed976'))+
  labs(x = 'Evenness', y = expression(Phytoplankton~C~'['~mu*mol*~L^-1~']'), fill = 'Fluctuation  \nfrequency [h]')+
  theme( panel.background = element_rect(fill = NA), #loescht den Hintergrund meines Plots/ fuellt ihn mit nichts
         #panel.grid.major.y = element_line(color='grey', linetype = 'dashed', size=0.2),
         panel.border= element_rect(colour = "black", fill=NA, size=0.5),
         strip.background = element_rect(color ='black', fill = 'white'),
         strip.text = element_text(face = 'italic'),
         legend.background = element_blank(),
         legend.position  ='none',
         legend.key = element_blank(),
         text = element_text(size=17))
BioEven

## make plots

#scatter plots to extract spearman corr values
ggscatter(subset(all_rich, index == 'evenness' ), x = "value", y = "carbon_umol_l", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman", xlab = 'pigment richness',
          ylab = 'Phytoplankton biomass')

# plot for results part: 
plot_grid( BioEven, BioComp, labels=c("(a)","(b)"),ncol = 2, label_size = 18, hjust = 0, vjust = 0.95)
#ggsave(plot = last_plot(), file = 'biomass_comp_relationship.tiff', width = 9, height = 4)

#######################################################################################
#######################################################################################

####CV and pigment compositinal turnover ####
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
         fluctuation = treatment) %>%
  group_by(fluctuation, sampling) %>%
  summarise(mean_pig = mean(distance, na.rm = T))
str(pig.dist)


#### Filter data ####
Mastertable_fluctron <- read_delim("~/Desktop/MA/MA_Rcode/project_data/Mastertable_fluctron.csv", 
                                   ";", escape_double = FALSE, locale = locale(decimal_mark = ","), 
                                   trim_ws = TRUE)
#View(Mastertable_fluctron)
str(Mastertable_fluctron)

# use only carbon data
cnp_data <- Mastertable_fluctron %>%
  dplyr::select(fluctuation, sampling, planktotron, carbon_umol_l)
names(cnp_data)

#carbon accumulated
co_data <- cnp_data %>%
  ungroup()%>%
  mutate(sampling = as.numeric(sampling)) %>%
  group_by(fluctuation, sampling) %>%
  summarise(mean = mean(carbon_umol_l, na.rm = T),
            sd = sd(carbon_umol_l, na.rm = T),
            se = sd/sqrt(n()),
            CV = sd/mean,
            se_CV = CV/sqrt(2*n()),
            n = n()) %>%
  mutate(lower.ci.mpg = CV - 1.96*se_CV/sqrt(n),
         upper.ci.mpg = CV + 1.96*se_CV/sqrt(n))%>%
  ungroup() %>%
  mutate(fluctuation = as.numeric(fluctuation))

all <- left_join(co_data, pig.dist, by = c('sampling', 'fluctuation')) %>%
  dplyr::select(sampling, CV, fluctuation, mean_pig)

formula <- y ~ x

all$fluctuation <- factor(as.factor(all$fluctuation),levels=c("0", "48", "36", '24', '12', '6'))
ggplot(all, aes(x = mean_pig, y = CV,group = fluctuation, fill = as.factor(fluctuation)))+
  geom_point(pch = 21, size = 3)+
  #geom_errorbar(aes(ymin = lower.ci.mpg, ymax = upper.ci.mpg), color = '#404040', width = .2, size = .5)+
  stat_smooth( method = "lm", se = T, size = 0.5,formula = formula, color = 'black') +
  stat_regline_equation(formula = formula, 
                        aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~"))) +     
  scale_fill_manual(values = c( '#000000','#0868ac','#41b6c4','#31a354','#addd8e','#fed976'))+
  labs(x = 'Compositional turnover', y = expression(coefficient~of~variance), fill = 'Fluctuation  \nfrequency [h]')+
  theme( panel.background = element_rect(fill = NA), #loescht den Hintergrund meines Plots/ fuellt ihn mit nichts
         #panel.grid.major.y = element_line(color='grey', linetype = 'dashed', size=0.2),
         panel.border= element_rect(colour = "black", fill=NA, size=0.5),
         strip.background = element_rect(color ='black', fill = 'white'),
         strip.text = element_text(face = 'italic'),
         legend.background = element_blank(),
         legend.position  ='bottom',
         legend.key = element_blank(),
         text = element_text(size=18))
#ggsave(plot = last_plot(), file = 'comp_funct_turnover.png', width = 5, height = 4)

ggscatter(all, x = "mean_pig", y = "CV", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman", xlab = 'compositional turnover',
          ylab = 'coefficient of variation')

#####
