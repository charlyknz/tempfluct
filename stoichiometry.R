#### new script to analyse stoichiometry ####
# by Charlotte Kunze 26.8.2020

library(tidyverse)

#### import data ####
Mastertable_fluctron <- read_delim("~/Desktop/MA/MA_Rcode/project_data/Mastertable_fluctron.csv", 
                                   ";", escape_double = FALSE, locale = locale(decimal_mark = ","), 
                                   trim_ws = TRUE)
#View(Mastertable_fluctron)
str(Mastertable_fluctron)

# import treatment information


### subdataframe with nutrient informations only ####

nutrients_master <- dplyr::select(Mastertable_fluctron, fluctuation, sampling, planktotron, carbon_umol_l, nitrate_umol_l, 
                                  'diss_Nitrat+Nitrit_umol_l', diss_Silikat_umol_l, diss_Phosphat_umol_l,
                                  POP_micromol_l, SiP_micromol_l)
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
  group_by(ratio, sampling, fluctuation) %>%
  summarise(mean = mean(value, na.rm = T),
         sd = sd(value, na.rm = T),
         se = sd/sqrt(n())) %>%
  drop_na(mean) %>%
  mutate(day = sampling *2)
ratio$fluctuation <- factor(as.factor(ratio$fluctuation),levels=c("0", "48", "36", '24', '12', '6'))

ggplot(ratio, aes( x = day, y = mean))+
  geom_point(aes(fill = fluctuation), pch = 21, size = 3.5, col = 'black')+
  geom_line(linetype = 'dashed' , aes(col = fluctuation,  group = fluctuation))+
  geom_errorbar(aes(ymin = mean -se, ymax = mean + se ), width = .5)+
  facet_wrap(~ratio, scales = 'free_y')+
  scale_fill_manual(values = c( '#000000', '#addd8e','#31a354','#41b6c4','#0868ac','#fed976'))+
  scale_color_manual(values = c( '#000000', '#addd8e','#31a354','#41b6c4','#0868ac','#fed976'))+
  labs(x = 'time (in days)', y = 'mean Molar ratio')+
  theme_classic()
ggsave(plot = last_plot(), file = 'molar_ratios.png')

