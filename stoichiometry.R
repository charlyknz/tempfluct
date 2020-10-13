#### new script to analyse stoichiometry ####
# by Charlotte Kunze 26.8.2020

library(tidyverse)
library(cowplot)


#### import data ####
Mastertable_fluctron <- read_delim("~/Desktop/MA/MA_Rcode/project_data/Mastertable_fluctron.csv", 
                                   ";", escape_double = FALSE, locale = locale(decimal_mark = ","), 
                                   trim_ws = TRUE)
#View(Mastertable_fluctron)
str(Mastertable_fluctron)

# import treatment information


### subdataframe with nutrient informations only ####

nutrients_master <- dplyr::select(Mastertable_fluctron, fluctuation, sampling, planktotron, carbon_umol_l, Chl.a_microg_l, nitrate_umol_l, n_ug_l,
                                  'diss_Nitrat+Nitrit_umol_l', 'diss_Nitrat+Nitrit_ug_l', diss_Silikat_umol_l, diss_Phosphat_umol_l,diss_Phosphat_ug_l,
                                  POP_micromol_l, POP_microg_l, SiP_micromol_l, srp_micromol_l)
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
  dplyr::select(fluctuation, sampling, planktotron,CN, NP, CP, CSi, NSi, NSi, SiP)%>%
  gather(key = 'ratio', value = 'value', -fluctuation, -sampling, -planktotron) %>%
 # mutate(log = log10(value))%>%
  group_by(ratio, sampling, fluctuation) %>%
  summarise(mean = mean(value, na.rm = T),
         sd = sd(value, na.rm = T),
         se = sd/sqrt(n())) %>%
  drop_na(mean) %>%
  mutate(day = sampling *2) %>%
  mutate(lower.ci = mean - 1.96*se/sqrt(n()),
         upper.ci = mean + 1.96*se/sqrt(n()),
         trans = 10^mean, 
         ci_l = 10^lower.ci,
         ci_u = 10^upper.ci) #create new column named trans_pred with transformed predictions.
  
ratio$fluctuation <- factor(as.factor(ratio$fluctuation),levels=c("0", "48", "36", '24', '12', '6'))

ggplot(ratio, aes( x = day, y = mean))+
  geom_point(aes(fill = fluctuation), pch = 21, size = 3.5, col = 'black')+
  geom_line(linetype = 'dashed' , aes(col = fluctuation,  group = fluctuation))+
  geom_errorbar(aes(ymin = mean -se, ymax = mean + se ), width = .5)+
  facet_wrap(~ratio, scales = 'free_y')+
  scale_fill_manual(values = c( '#000000', '#addd8e','#31a354','#41b6c4','#0868ac','#fed976'))+
  scale_color_manual(values = c( '#000000', '#addd8e','#31a354','#41b6c4','#0868ac','#fed976'))+
  labs(x = 'Time (in days)', y = 'mean Molar ratio')+
  theme_classic()
#ggsave(plot = last_plot(), file = 'molar_ratios.png')

#### calculate RUE ####
# RUE = Biomass (micromol) / TP bzw. TN 

# 1. calculate Total Phosphoros/ Nitrate
# diss Nut + Part Nut
# RUE = biomass / TP/N
# mean of RUE per treatment

## use molar values
RUE12 <- nutrients_master %>%
  dplyr::select(fluctuation, planktotron, sampling, srp_micromol_l, carbon_umol_l, nitrate_umol_l, POP_micromol_l, 'diss_Nitrat+Nitrit_umol_l', 'diss_Phosphat_umol_l') %>%
  dplyr::rename(diss_P = 'diss_Phosphat_umol_l',
         diss_N = 'diss_Nitrat+Nitrit_umol_l') %>%
  mutate(TP = diss_P + POP_micromol_l,
         TN = diss_N + nitrate_umol_l) %>%
  drop_na(carbon_umol_l) %>%
  mutate(P_RUE = carbon_umol_l/ TP,
         N_RUE = carbon_umol_l/ TN) %>%
  gather(key = 'nutrient', value = 'value', -fluctuation, -sampling, -planktotron) %>%
  dplyr::group_by(fluctuation, sampling, nutrient) %>%
  dplyr::summarise(mean_N = mean(value, na.rm = T),
            sd_N = sd(value, na.rm = T),
            n = dplyr::n(),
            se_N = sd_N/sqrt(n)) %>%
  mutate(lower.ci = mean_N - 1.96*se_N/sqrt(n),
         upper.ci = mean_N + 1.96*se_N/sqrt(n)) %>%
  filter(mean_N > 0)%>%
  mutate(day = sampling *2)

RUE12$fluctuation <- factor(as.factor(RUE12$fluctuation),levels=c("0", "48", "36", '24', '12', '6'))


##########################################################
RUE_N <- ggplot(subset(RUE12, nutrient == 'N_RUE'), aes(x = day, y = mean_N))+
  geom_line(linetype = 'dashed', aes(col = fluctuation))+
 # geom_errorbar(aes(ymin = mean_N -se_N, ymax =  mean_N +se_N), width = .5)+
  geom_errorbar(aes(ymin = lower.ci, ymax = upper.ci, color = fluctuation), width = .5,position = position_dodge2(width = .5))+
  geom_point(aes(fill = fluctuation), pch = 21, col = 'black', size = 3,position = position_dodge2(width = .5))+
  scale_fill_manual(values = c( '#000000','#0868ac','#41b6c4','#31a354','#addd8e','#fed976'))+
  scale_color_manual(values = c( '#000000','#0868ac','#41b6c4','#31a354','#addd8e','#fed976'))+
  labs(x = 'Time [days]', y = 'RUE (Total N)')+
  #facet_wrap(~nutrient, scales = 'free_y')+
  theme( panel.background = element_rect(fill = NA), #loescht den Hintergrund meines Plots/ fuellt ihn mit nichts
         #panel.grid.major.y = element_line(color='grey', linetype = 'dashed', size=0.2),
         panel.border= element_rect(colour = "black", fill=NA, size=1),
         strip.background = element_rect(color ='black', fill = 'white'),
         strip.text = element_text(face = 'bold'),
         legend.background = element_blank(),
         legend.position  ='none',
         legend.key = element_blank(),
         text = element_text(size=18))
RUE_N


RUE_P <- ggplot(subset(RUE12, nutrient == 'P_RUE'), aes(x = day, y = mean_N))+
  geom_line(linetype = 'dashed', aes(col = fluctuation))+
  #geom_errorbar(aes(ymin = mean_N -se_N, ymax =  mean_N +se_N), width = .5,position = position_dodge2(width = .5))+
  geom_errorbar(aes(ymin = lower.ci, ymax = upper.ci, color = fluctuation), width = .5,position = position_dodge2(width = .5))+
  geom_point(aes(fill = fluctuation), pch = 21, col = 'black', size = 3,position = position_dodge2(width = .5))+
  scale_fill_manual(values = c( '#000000','#0868ac','#41b6c4','#31a354','#addd8e','#fed976'))+
  scale_color_manual(values = c( '#000000','#0868ac','#41b6c4','#31a354','#addd8e','#fed976'))+
  labs(x = 'Time [days]', y = 'RUE (Total P)')+
  #facet_wrap(~nutrient, scales = 'free_y')+
  theme( panel.background = element_rect(fill = NA), #loescht den Hintergrund meines Plots/ fuellt ihn mit nichts
         #panel.grid.major.y = element_line(color='grey', linetype = 'dashed', size=0.2),
         panel.border= element_rect(colour = "black", fill=NA, size=1),
         strip.background = element_rect(color ='black', fill = 'white'),
         strip.text = element_text(face = 'bold'),
         legend.background = element_blank(),
         legend.position  ='none',
         legend.key = element_blank(),
         text = element_text(size=18))
RUE_P


CN <- ggplot(subset(ratio, ratio == 'CN') , aes( x = day, y = mean))+
  geom_line(linetype = 'dashed' , aes(col = fluctuation,  group = fluctuation))+
  geom_errorbar(aes(ymin = lower.ci, ymax =upper.ci , color = fluctuation), width = .5,position = position_dodge2(width = .5))+
  geom_point(aes(fill = fluctuation), pch = 21, size = 3, col = 'black',position = position_dodge2(width = .5))+
  scale_fill_manual(values = c( '#000000','#0868ac','#41b6c4','#31a354','#addd8e','#fed976'))+
  scale_color_manual(values = c( '#000000','#0868ac','#41b6c4','#31a354','#addd8e','#fed976'))+
  labs(x = ' ', y = 'C:N ratio')+
  theme( panel.background = element_rect(fill = NA), #loescht den Hintergrund meines Plots/ fuellt ihn mit nichts
         #panel.grid.major.y = element_line(color='grey', linetype = 'dashed', size=0.2),
         panel.border= element_rect(colour = "black", fill=NA, size=1),
         strip.background = element_rect(color ='black', fill = 'white'),
         strip.text = element_text(face = 'bold'),
         legend.background = element_blank(),
         legend.position  ='none',
         legend.key = element_blank(),
         text = element_text(size=18))
CN

CP <- ggplot(subset(ratio, ratio == 'CP') , aes( x = day, y = mean))+
  geom_line(linetype = 'dashed' , aes(col = fluctuation,  group = fluctuation))+
  geom_errorbar(aes(ymin = lower.ci, ymax = upper.ci, col = fluctuation), width = .5,position = position_dodge2(width = .5))+
  geom_point(aes(fill = fluctuation), pch = 21, size = 3, col = 'black',position = position_dodge2(width = .5))+
  scale_fill_manual(values = c( '#000000','#0868ac','#41b6c4','#31a354','#addd8e','#fed976'))+
  scale_color_manual(values = c( '#000000','#0868ac','#41b6c4','#31a354','#addd8e','#fed976'))+
  labs(x = '', y = 'C:P ratio')+
  theme( panel.background = element_rect(fill = NA), #loescht den Hintergrund meines Plots/ fuellt ihn mit nichts
         #panel.grid.major.y = element_line(color='grey', linetype = 'dashed', size=0.2),
         panel.border= element_rect(colour = "black", fill=NA, size=1),
         strip.background = element_rect(color ='black', fill = 'white'),
         strip.text = element_text(face = 'bold'),
         legend.background = element_blank(),
         legend.position  ='none',
         legend.key = element_blank(),
         text = element_text(size=18))
CP

CSi <- ggplot(subset(ratio, ratio == 'CSi') , aes( x = day, y = mean))+
  geom_line(linetype = 'dashed' , aes(col = fluctuation,  group = fluctuation))+
  geom_errorbar(aes(ymin = lower.ci, ymax = upper.ci, col = fluctuation), width = .5,position = position_dodge2(width = .5))+
  geom_point(aes(fill = fluctuation), pch = 21, size = 3, col = 'black',position = position_dodge2(width = .5))+
  scale_fill_manual(values = c( '#000000','#0868ac','#41b6c4','#31a354','#addd8e','#fed976'))+
  scale_color_manual(values = c( '#000000','#0868ac','#41b6c4','#31a354','#addd8e','#fed976'))+
  labs(x = 'Time [days]', y = 'C:Si ratio', color = 'Fluctuation frequency (in h)',fill = 'Fluctuation frequency (in h)')+
  theme( panel.background = element_rect(fill = NA), #loescht den Hintergrund meines Plots/ fuellt ihn mit nichts
         #panel.grid.major.y = element_line(color='grey', linetype = 'dashed', size=0.2),
         panel.border= element_rect(colour = "black", fill=NA, size=1),
         strip.background = element_rect(color ='black', fill = 'white'),
         strip.text = element_text(face = 'bold'),
         legend.background = element_blank(),
         legend.position  ='none',
         legend.key = element_blank(),
         text = element_text(size=18))
CSi

NP <- ggplot(subset(ratio, ratio == 'NP') , aes( x = day, y = mean))+
  geom_line(linetype = 'dashed' , aes(col = fluctuation,  group = fluctuation))+
  geom_errorbar(aes(ymin = lower.ci, ymax = upper.ci, col = fluctuation), width = .5,position = position_dodge2(width = .5))+
  geom_point(aes(fill = fluctuation), pch = 21, size = 3, col = 'black',position = position_dodge2(width = .5))+
  scale_fill_manual(values = c( '#000000','#0868ac','#41b6c4','#31a354','#addd8e','#fed976'))+
  scale_color_manual(values = c( '#000000','#0868ac','#41b6c4','#31a354','#addd8e','#fed976'))+
  labs(x = 'Time [days]', y = 'N:P ratio', color = 'Fluctuation frequency (in h)',fill = 'Fluctuation frequency (in h)')+
  theme( panel.background = element_rect(fill = NA), #loescht den Hintergrund meines Plots/ fuellt ihn mit nichts
         #panel.grid.major.y = element_line(color='grey', linetype = 'dashed', size=0.2),
         panel.border= element_rect(colour = "black", fill=NA, size=1),
         strip.background = element_rect(color ='black', fill = 'white'),
         strip.text = element_text(face = 'bold'),
         legend.background = element_blank(),
         legend.position  ='none',
         legend.key = element_blank(),
         text = element_text(size=18))
NP


plot_grid(CN, CP,CSi,NP, labels=c("(a)","(b)", '(c)', '(d)'),ncol = 2, label_size = 18, hjust = 0, vjust = 1)

#ggsave(plot = last_plot(), file = 'MolarRatio.tiff', width = 9, height = 7)

