##
# R Skript to analyse pigment composition
library(tidyverse)
## load data
all_pigments <- read.csv2("~/Desktop/MA/MA_data/composition/Pigments/Pigmente_R_script/all_pigments.csv", 
                           sep = ";", dec = ',')

#View(all_pigments)
df <- data.frame( planktotron = c('1', '4', '2', '10', '3', '5', '6', '7', '8', '9', '11', '12'),
                  treatment = c('control', 'control', 'Fluctuating_36','Fluctuating_36',
                                'Fluctuating_6', 'Fluctuating_6', 'Fluctuating_48', 'Fluctuating_48',
                                'Fluctuating_24', 'Fluctuating_24', 'Fluctuating_12','Fluctuating_12') )

pigi <- all_pigments %>%
  filter(planktotron != 'Blank') %>%
  select(-X) %>%
  spread(key = 'pigments', value = 'value') %>%
    mutate(rel_anth = Anth/Chl.a)%>% ###berechnung der relativen werte in neuer spalte nach Schlueter et al.
    mutate(rel_cantha = Cantha/Chl.a)%>%
    mutate(rel_chl.a = Chl.a/Chl.a)%>%
    mutate(rel_chl.b = Chl.b/Chl.a)%>%
    mutate(rel_chl.c2 = Chl.c2/Chl.a)%>%
    mutate(rel_diato = Diato/Chl.a)%>%
    mutate(rel_echin = Echin/Chl.a)%>%
    mutate(rel_fuco = Fuco/Chl.a)%>%
    mutate(rel_myxo = Myxo/Chl.a)%>%
    mutate(rel_peri = Peri/Chl.a)%>%
    mutate(rel_phe.a = Phe.a/Chl.a)%>%
    mutate(rel_viola = Viola/Chl.a) %>%
  select(-c("Allo", "Anth","bb.Car", "c.Neo", "Cantha","Chl.a","Chl.b","Chl.c1",
            "Chl.c2","Cryp","Diadino", "Diato","Dino","Echin","Fuco", "Lut","Myxo","Peri",       
            "Phe.a","Phe.b", "Viola", "Zea" )) %>%
  gather(key = 'pigment', value = 'value', -planktotron, - sampling)

data <- left_join(pigi, df, by = c('planktotron'))

data%>%
  filter(sampling <14)%>%
ggplot(aes(x = sampling, y = value, fill = pigment))+
  geom_bar(stat = 'identity')+
  labs(x = 'sampling', y = 'relative_contribution')+
  facet_wrap(~treatment)+
  theme_bw()

#ggsave(plot = last_plot(), file = 'pigments.png')
#