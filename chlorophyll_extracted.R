
library(tidyverse)

data <- read.csv2("~/Desktop/MA/MA_data/composition/Pigments/Kalibrierung_extracted_chl.csv", sep = ';', dec = ',')

df <- data.frame( sample = c('P1', 'P4', 'P2', 'P10', 'P3', 'P5', 'P6', 'P7', 'P8', 'P9', 'P11', 'P12'),
                  treatment = c('control', 'control', 'Fluctuating_36','Fluctuating_36',
                                   'Fluctuating_6', 'Fluctuating_6', 'Fluctuating_48', 'Fluctuating_48',
                                   'Fluctuating_24', 'Fluctuating_24', 'Fluctuating_12','Fluctuating_12') )

chl_data <- left_join(data, df, by = c('sample')) %>%
  filter(treatment != 'Blank') %>%
  group_by(treatment, sampling) %>%
  mutate(mean = mean(conc_liter, na.rm = T),
         log = log(mean),
         sd = sd(conc_liter, na.rm = T),
         se = sd/sqrt(n()))
#write.csv2(x = chl_data, file = 'extracted_chl.csv')
ggplot(chl_data, aes( x = sampling, y = log, col = treatment, group = treatment))+
  geom_point()+
  geom_line()+
 # geom_errorbar(aes(ymax= mean+se, ymin =mean-se))+
  labs(x = 'sampling', y = 'mean Chlorophyll a')+
  theme_bw()
#ggsave(plot = last_plot(), file = 'extracted_chl.png')
