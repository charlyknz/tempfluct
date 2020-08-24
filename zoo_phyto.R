## zoo and phyto data in one plot
#import all data

Mastertable <- read_delim("~/Desktop/MA/MA_Rcode/project_data/Mastertable.csv", 
                          ";", escape_double = FALSE, locale = locale(decimal_mark = ","), 
                          trim_ws = TRUE) %>%
  rename(MC = planktotron)
str(Mastertable)

## merge with treatment information
treatments <- read.csv2('~/Desktop/MA/MA_Rcode/project_data/treatments_units.csv')
str(treatments)
names(treatments) = c('MC', 'fluctuation', 'treatment')
master_data <- left_join(Mastertable, treatments, by = c('MC')) 

  
master <- master_data %>%
  select(MC, sampling, treatment, fluctuation, carbon_umol_l, C_Zoo_µmol_l ) %>%
  gather(key = 'variable', value = 'carbon', -MC, -sampling, -treatment, -fluctuation, ) %>%
  mutate(dummy = paste(ifelse(variable == 'carbon_umol_l', 'phytoplankton', 'zooplankton'))) %>%
  group_by(sampling, treatment, fluctuation, dummy) %>%
  summarise(mean = mean(carbon, na.rm = T),
            sd = sd(carbon, na.rm = T),
            se = sd/sqrt(n())) %>%
  mutate(treatment_dummy = paste(treatment, dummy, sep = '_')) %>%
  drop_na(mean)

master$treatment = factor(as.factor(master$treatment), levels = c('con', 'F48', 'F36', 'F24', 'F12', 'F6'))
ggplot(master, aes(x = sampling, y = mean,group = treatment_dummy))+
  geom_point(aes(fill = treatment),size = 3, pch = 21, color = 'black')+
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = .5)+
  geom_line(linetype = 'dashed')+
  labs(y = 'mean carbon umol/l')+
  scale_fill_manual(values = c( '#000000','#0868ac','#41b6c4','#31a354','#addd8e','#fed976'))+
  facet_wrap(~dummy, scales = 'free_y')+
  theme_classic()+
  theme(legend.position = 'bottom')
ggsave(plot = last_plot(), file = 'zoo_phyto.png', width = 8, height = 5)


#############################################################
nutri <- master_data %>%
  select(MC, sampling, treatment, fluctuation, mean_POP_mg_l, mean_part_si_ug_l, mean_srp_mg_l,hand_si_mg_l) %>%
  gather(key = 'variable', value = 'value', -MC, -sampling, -treatment, -fluctuation, ) %>%
  group_by(sampling, treatment, fluctuation, variable) %>%
  summarise(mean = log10(mean(value, na.rm = T)),
            sd = sd(value, na.rm = T),
            se = log10(sd/sqrt(n()))) %>%
  drop_na(mean)

nutri$treatment = factor(as.factor(nutri$treatment), levels = c('con', 'F48', 'F36', 'F24', 'F12', 'F6'))
ggplot(subset(nutri, !variable %in% c('hand_si_mg_l', 'mean_srp_mg_l')), aes(x = sampling, y = mean,  group = treatment, col = treatment))+
  geom_point()+
  geom_line()+
  facet_wrap(~variable, scales = 'free_y')+
  labs(y = 'mean values (in ug/mg/l)')+
  scale_color_manual(values = c( '#000000','#0868ac','#41b6c4','#31a354','#addd8e','#fed976'))+
  theme_classic()
#ggsave(plot = last_plot(), file = 'part_si_POP.png')
