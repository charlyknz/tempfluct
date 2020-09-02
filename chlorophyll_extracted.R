
library(tidyverse)

data <- read.csv2("~/Desktop/MA/MA_Rcode/project_data/Kalibrierung_extracted_chl.csv", sep = ';', dec = ',')

df <- data.frame( sample = c('P1', 'P4', 'P2', 'P10', 'P3', 'P5', 'P6', 'P7', 'P8', 'P9', 'P11', 'P12'),
                  treatment = c('control', 'control', 'Fluctuating_36','Fluctuating_36',
                                   'Fluctuating_6', 'Fluctuating_6', 'Fluctuating_48', 'Fluctuating_48',
                                   'Fluctuating_24', 'Fluctuating_24', 'Fluctuating_12','Fluctuating_12') )

chl_data <- left_join(data, df, by = c('sample')) %>%
  filter(treatment != 'Blank') %>%
  group_by(treatment, sampling, sample) %>%
  summarise(mean = mean(conc_liter, na.rm = T),
         sd = sd(conc_liter, na.rm = T),
         se = sd/sqrt(n()),
         se_log = log(se))%>%
arrange(sampling, sample)
#write.csv2(x = chl_data, file = 'extracted_chl.csv')
ggplot(chl_data, aes( x = sampling, y = mean, group = treatment))+
  geom_line(linetype = 'dashed', aes(color = treatment))+
  geom_point(aes(fill=treatment), colour="black",pch=21, size=3, position = position_dodge(width = 0.8))+
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = .5, size = .4, position = position_dodge(width = 0.8))+
  scale_fill_manual(values = c( '#000000', '#addd8e','#31a354','#41b6c4','#0868ac','#fed976'))+
  scale_color_manual(values = c( '#000000', '#addd8e','#31a354','#41b6c4','#0868ac','#fed976'))+
  labs(x = 'sampling', y = expression(Chlorophyll~mu*g*~L^-1))+
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
#ggsave(plot = last_plot(), file = 'extracted_chl.png')
