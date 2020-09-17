## script to analyse nutrients

library(tidyverse)
library(viridis)

# load dis nut data
Charly_disnut <- read_delim("~/Desktop/MA/MA_Rcode/project_data/Charly_disnut.csv", 
                            ";", escape_double = FALSE, trim_ws = TRUE)
str(Charly_disnut)
#change names
data <- Charly_disnut[1:8] #remove empty columns, keep only 1:8
names(data) <- c('sample', 'nitrate_nitrit', 'nitrate_nitrit_mol', 'silicate', 'silicate_mol',
                 'phosphate', 'phosphate_mol', 'day')

#subsitute comma with dot and change format to numeric for nutrients
data$nitrate_nitrit = as.numeric(gsub("\\,", ".",  data$nitrate_nitrit))
data$nitrate_nitrit_mol = as.numeric(gsub("\\,", ".", data$nitrate_nitrit_mol))
data$silicate = as.numeric(gsub("\\,", ".",  data$silicate))
data$silicate_mol = as.numeric(gsub("\\,", ".",  data$silicate_mol))
data$phosphate = as.numeric(gsub("\\,", ".",  data$phosphate))
data$phosphate_mol = as.numeric(gsub("\\,", ".", data$phosphate_mol))

str(data)# check new structure

df <- data.frame( sample = c('P1', 'P4', 'P2', 'P10', 'P3', 'P5', 'P6', 'P7', 'P8', 'P9', 'P11', 'P12'),
                  treatment = c('control', 'control', 'Fluctuating_36','Fluctuating_36',
                                'Fluctuating_6', 'Fluctuating_6', 'Fluctuating_48', 'Fluctuating_48',
                                'Fluctuating_24', 'Fluctuating_24', 'Fluctuating_12','Fluctuating_12') )

nut <- left_join(data, df, by = 'sample') %>%
  drop_na(treatment) %>%
  gather(key = 'nutrient', value = 'value', - sample,-day,-treatment) %>%
  filter(nutrient %in% c('nitrate_nitrit_mol', 'phosphate_mol', 'silicate_mol')) %>%
  group_by(treatment, day, nutrient) %>%
  summarise(mean = mean(value, na.rm = T),
            sd = sd(value, na.rm =T),
            se = sd/sqrt(n())) %>%
  mutate(days = day *2,
         lower.ci = mean - 1.96*se/sqrt(n()),
         upper.ci = mean + 1.96*se/sqrt(n())) %>%
  filter(day != 8) #remove outlier


nut$treatment = factor(as.factor(nut$treatment), levels = c('control', 'Fluctuating_48', 'Fluctuating_36', 'Fluctuating_24', 'Fluctuating_12', 'Fluctuating_6'))

label_n <- c(nitrate_nitrit_mol= 'Nitrate Nitrite' ,phosphate_mol= 'Phosphate', silicate_mol='Silicate')
ggplot(nut, aes(x = days, y = mean, group = treatment))+
  geom_errorbar(aes(ymin = lower.ci, ymax = upper.ci, color = treatment), width = .5, position = position_dodge2(width = .5))+
  geom_line(linetype = 'dashed')+
  geom_point(position = position_dodge(width = 0.5), aes(fill = treatment), pch = 21, size = 3)+
  #scale_x_continuous(limits = c(0,18), breaks = seq(0,18,2))+
  scale_fill_manual(values = c( '#000000','#0868ac','#41b6c4','#31a354','#addd8e','#fed976'))+
  scale_color_manual(values = c( '#000000','#0868ac','#41b6c4','#31a354','#addd8e','#fed976'))+
  facet_wrap(~nutrient,scales = 'free_y', labeller = labeller(nutrient = label_n))+
  labs(x = 'Time [days]', y = expression(Dissolved~nutrients~'['~mu*mol*~L^-1~']'))+
  theme_classic()+
  theme( panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         strip.background = element_blank(),
         panel.border= element_rect(colour = "black", fill=NA, size=0.5),
         strip.text = element_text(face = 'bold'),
         legend.background = element_blank(),
         legend.position  ='none',
         legend.key = element_blank(),
         text = element_text(size=18))
#ggsave(plot = last_plot(), file = 'nutrients.png', width = 9, height = 5) 


###########################################################################################
#### Filter data ####
cnp_data <- read.csv2('~/Desktop/MA/MA_Rcode/project_data/cnp_data_treatments.csv')
str(cnp_data)

#carbon accumulated
co_data <- cnp_data %>%
  ungroup()%>%
  mutate(sampling = str_remove(sample, 'S'),
         sampling = as.numeric(sampling)) %>%
  group_by(treatment) %>%
  summarise(mean = mean(c_umol_l, na.rm = T),
            sd = sd(c_umol_l, na.rm = T),
            se = sd/sqrt(n()),
            CV = sd/mean,
            se_CV = CV/sqrt(2*n()),
            n = n()) %>%
  mutate(lower.ci.mpg = CV - 1.96*se_CV/sqrt(n),
         upper.ci.mpg = CV + 1.96*se_CV/sqrt(n))%>%
  mutate(fluctuation = paste(ifelse(treatment == 'control', 0, ifelse(treatment == 'Fluctuating_48', 48, ifelse(treatment == 'Fluctuating_36', 36,
  ifelse(treatment == 'Fluctuating_24', 24, ifelse(treatment == 'Fluctuating_12', 12, 6)))))))%>%
  ungroup() 

co_data$fluctuation <- factor(as.factor(co_data$fluctuation),levels=c("0", "48", "36", '24', '12', '6'))
ggplot(co_data, aes(x = fluctuation, y = CV,group = treatment, fill = as.factor(fluctuation)))+
  geom_col()+
  geom_errorbar(aes(ymin = lower.ci.mpg, ymax = upper.ci.mpg), color = '#404040', width = .2, size = .5)+
  scale_fill_manual(values = c( '#000000','#0868ac','#41b6c4','#31a354','#addd8e','#fed976'))+
  labs(x = 'Fluctuation frequency [h]', y = expression(coefficient~of~variance), fill = 'Fluctuation  \nfrequency [h]')+
  theme( panel.background = element_rect(fill = NA), #loescht den Hintergrund meines Plots/ fuellt ihn mit nichts
         #panel.grid.major.y = element_line(color='grey', linetype = 'dashed', size=0.2),
         panel.border= element_rect(colour = "black", fill=NA, size=0.5),
         strip.background = element_rect(color ='black', fill = 'white'),
         strip.text = element_text(face = 'italic'),
         legend.background = element_blank(),
         legend.position  ='bottom',
         legend.key = element_blank(),
         text = element_text(size=14))
#ggsave(plot = last_plot(), file = 'CV_carbon_accumulated.png', width = 9, height = 6)

##########################################################################
#### Si and P ####
part_si_data <- read_delim("~/Desktop/MA/MA_Rcode/project_data/part_si_data.csv", 
                           ";", escape_double = FALSE, locale = locale(decimal_mark = ","), 
                           trim_ws = TRUE) 
part_si_data$Planktotron = as.character(part_si_data$Planktotron)
POP_data <- read_delim("~/Desktop/MA/MA_Rcode/project_data/POP_data.csv", 
                       ";", escape_double = FALSE, locale = locale(decimal_mark = ","), 
                       trim_ws = TRUE)
POP_data$Planktotron = as.character(POP_data$Planktotron)

df <- data.frame( Planktotron = c('1', '4', '2', '10', '3', '5', '6', '7', '8', '9', '11', '12'),
                  treatment_id = c('control', 'control', 'Fluctuating_36','Fluctuating_36',
                                   'Fluctuating_6', 'Fluctuating_6', 'Fluctuating_48', 'Fluctuating_48',
                                   'Fluctuating_24', 'Fluctuating_24', 'Fluctuating_12','Fluctuating_12') )
names(part_si_data)
Psidata <- part_si_data %>%
  left_join(., df, by = c('Planktotron')) %>%
  group_by(Sampling, treatment_id) %>%
  summarise(mean = mean(conc_mol_l, na.rm = T),
            sd = sd(conc_mol_l, na.rm = T),
            se = sd/sqrt(n())) %>%
  drop_na(mean) %>%
  mutate(nutrient = paste('PSi'))

POPdata <-POP_data %>%
  left_join(., df, by = c('Planktotron')) %>%
  group_by(Sampling, treatment_id) %>%
  summarise(mean = mean(POP_conc_micromolP_L, na.rm = T),
            sd = sd(POP_conc_micromolP_L, na.rm = T),
            se = sd/sqrt(n())) %>%
  drop_na(mean)%>%
  mutate(nutrient = paste('POP'))

Pdata <- bind_rows(POPdata, Psidata) %>%
  mutate(day = Sampling *2) %>%
  mutate(lower.ci = mean - 1.96*se/sqrt(n()),
         upper.ci = mean + 1.96*se/sqrt(n())) %>%
  ungroup()%>%
  separate(treatment_id, into = c('treatment', 'fluctuation'), sep = '_')%>%
dplyr::select(day, lower.ci, upper.ci, mean, nutrient, fluctuation) 

Pdata$fluctuation[is.na(Pdata$fluctuation)] <-0

#### add particulate N #####
n_data <- cnp_data %>%
  ungroup()%>%
  mutate(sampling = str_remove(sample, 'S'),
         sampling = as.numeric(sampling)) %>%
  group_by(treatment, sampling) %>%
  summarise(mean = mean(c_umol_l, na.rm = T),
            sd = sd(n_umol_l, na.rm = T),
            se = sd/sqrt(n())) %>%
  mutate(fluctuation = paste(ifelse(treatment == 'control', 0, ifelse(treatment == 'Fluctuating_48', 48, ifelse(treatment == 'Fluctuating_36', 36,
                                                                                                                ifelse(treatment == 'Fluctuating_24', 24, ifelse(treatment == 'Fluctuating_12', 12, 6))))))) %>%
  mutate(day = sampling *2) %>%
  mutate(lower.ci = mean - 1.96*se/sqrt(n()),
         upper.ci = mean + 1.96*se/sqrt(n()),
         nutrient = paste('PN')) %>%
  ungroup() %>%
  dplyr:: select(day, lower.ci, upper.ci, mean, nutrient, fluctuation)


all <- bind_rows(Pdata, n_data)
all$fluctuation <- factor(as.factor(all$fluctuation),levels=c("0", "48", "36", '24', '12', '6'))
                             
ggplot(all, aes(x= day, y = mean, group = fluctuation))+
  geom_line(linetype = 'dashed', aes(color = fluctuation))+
  geom_errorbar(aes(ymin = lower.ci, ymax = upper.ci, color = fluctuation), width = .5, position = position_dodge(width = .5))+
  geom_point(aes(fill = fluctuation), pch = 21, col = 'black', size = 3, position = position_dodge(width = .5))+
  scale_fill_manual(values = c( '#000000','#0868ac','#41b6c4','#31a354','#addd8e','#fed976'))+
  scale_color_manual(values = c( '#000000','#0868ac','#41b6c4','#31a354','#addd8e','#fed976'))+
  facet_wrap(~nutrient,scales = 'free_y')+
  labs(x = 'Time [days]', y = expression(Particulate~nutrients~'['~mu*mol*~L^-1~']'))+
  theme_classic()+
  theme( panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         strip.background = element_blank(),
         panel.border= element_rect(colour = "black", fill=NA, size=0.5),
         strip.text = element_text(face = 'bold'),
         legend.background = element_blank(),
         legend.position  ='none',
         legend.key = element_blank(),
         text = element_text(size=18))
#ggsave(plot = last_plot(), file = 'P_nutrient.png', width = 9, height = 5)
  
