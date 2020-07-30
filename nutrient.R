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
  filter(nutrient %in% c('nitrate_nitrit', 'phosphate', 'silicate'))

nut %>%
  filter(nutrient %in% c('phosphate', 'nitrate_nitrit'))%>%
ggplot(., aes(x = day, y = value, col = nutrient, group = nutrient))+
  geom_point()+
  geom_smooth()+
  facet_grid(~treatment, scales = 'free')+
  theme_bw()
#ggsave(plot = last_plot(), file = 'nitrate_nitrite_phosphate.png')  

nut %>%
  filter(nutrient == 'silicate')%>%
  ggplot(., aes(x = day, y = value, group = treatment))+
  geom_point()+
  geom_smooth()+
  facet_grid(~treatment, scales = 'free')+
  theme_bw()
#ggsave(plot = last_plot(), file = 'silicate.png') 


###########################################################################################

cnp_data <- read.csv2('~/Desktop/MA/MA_Rcode/project_data/cnp_data_treatments.csv')
str(cnp_data)
cnp_data <- cnp_data %>%
  ungroup()%>%
  mutate(sampling = str_remove(sample, 'S'),
         sampling = as.numeric(sampling)) %>%
  group_by(treatment, sampling) %>%
  summarise(mean = mean(c_umol_l, na.rm = T),
         sd = sd(c_umol_l, na.rm = T),
         se = sd/sqrt(n()))

ggplot(cnp_data, aes(x = sampling, y = mean,group = treatment))+
  geom_line(linetype = 'dashed', aes(color = treatment))+
  geom_point(aes(fill=treatment), colour="black",pch=21, size=3, position = position_dodge(width = 1))+
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = .5, size = .4, position = position_dodge(width = 1))+
  scale_fill_manual(values = c( '#000000', '#addd8e','#31a354','#41b6c4','#0868ac','#fed976'))+
  scale_color_manual(values = c( '#000000', '#addd8e','#31a354','#41b6c4','#0868ac','#fed976'))+
  labs(x = 'sampling', y = expression(carbon~mu*mol*~L^-1))+
  scale_x_continuous(limits = c(-1,19), breaks = seq(0,18,2))+
  theme( panel.background = element_rect(fill = NA), #loescht den Hintergrund meines Plots/ fuellt ihn mit nichts
         #panel.grid.major.y = element_line(color='grey', linetype = 'dashed', size=0.2),
         panel.border= element_rect(colour = "black", fill=NA, size=0.5),
         strip.background = element_rect(color ='black', fill = 'white'),
         strip.text = element_text(face = 'italic'),
         legend.background = element_blank(),
         legend.position  ='bottom',
         legend.key = element_blank(),
         text = element_text(size=14))
#ggsave(plot = last_plot(), file = 'carbon_over_time.png')





