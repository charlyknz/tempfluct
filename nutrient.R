## script to analyse nutrients

library(tidyverse)

# load dis nut data
Charly_disnut <- read_delim("~/Desktop/MA/MA_Rcode/tempfluct/Charly_disnut.csv", 
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
  ggplot(., aes(x = day, y = value, group = sample))+
  geom_point()+
  geom_smooth()+
  facet_grid(~treatment, scales = 'free')+
  theme_bw()
#ggsave(plot = last_plot(), file = 'silicate.png') 


###########################################################################################

cnp_data <- read.csv2('~/Desktop/MA/MA_Rcode/tempfluct/CNP_filter_data.csv')
str(cnp_data)

df <- data.frame( MC = c('P1', 'P4', 'P2', 'P10', 'P3', 'P5', 'P6', 'P7', 'P8', 'P9', 'P11', 'P12'),
                  treatment = c('control', 'control', 'Fluctuating_36','Fluctuating_36',
                                'Fluctuating_6', 'Fluctuating_6', 'Fluctuating_48', 'Fluctuating_48',
                                'Fluctuating_24', 'Fluctuating_24', 'Fluctuating_12','Fluctuating_12') )

cnp_data <- cnp_data[1:10] %>%
  separate(Probe, c('MC', 'sampling', 'no'), ' ')
cnp <- cnp_data %>%
  drop_na(no) %>%
  select(-Dateiname, -Inj..Datum, -Uhrzeit, -Typ)

all_data <- left_join(cnp, df, by = 'MC') %>% 
  group_by(sampling, treatment) %>%
  mutate(mean = mean(C.µg, na.rm = T))
all_data$sample = factor(all_data$sampling, levels = c('S0', 'S2', 'S4', 'S6', 'S8', 'S10', 'S12', 'S14', 'S16', 'S18'))


ggplot(all_data, aes(x = sample, y = C.µg, col = treatment, group = treatment))+
  geom_point(size = 3)+
 # geom_line()+
  theme_bw()






