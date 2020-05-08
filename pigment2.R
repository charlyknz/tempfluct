##
# R Skript to analyse pigment composition
library(tidyverse)
## load data
all_pigments <- read.csv2("~/Desktop/MA/MA_Rcode/project_data/all_pigments.csv", 
                          sep = ";", dec = ',')

#View(all_pigments)
df <- data.frame( planktotron = c('1', '4', '2', '10', '3', '5', '6', '7', '8', '9', '11', '12'),
                  treatment = c('control', 'control', 'Fluctuating_36','Fluctuating_36',
                                'Fluctuating_6', 'Fluctuating_6', 'Fluctuating_48', 'Fluctuating_48',
                                'Fluctuating_24', 'Fluctuating_24', 'Fluctuating_12','Fluctuating_12') )

pigi <- all_pigments %>%
  filter(planktotron != 'Blank') %>%
  select(-X) 
data <- left_join(pigi, df, by = c('planktotron'))

master_data <- data %>%
  group_by(treatment,sampling) %>%
  #  mutate(sum = sum(value)) %>%
  spread(key = 'pigments', value = 'value') %>%
  arrange(sampling, planktotron) %>%
  mutate(planktotron = as.numeric(planktotron))%>%
  left_join(jacuzzi_no, by = c('planktotron', 'sampling')) %>%
  select(no, planktotron, sampling, date, everything())
#write.csv2(x = master_data, file = 'master_pigments.csv')

pig_data <- data %>%
  group_by(treatment,sampling) %>%
  mutate(sum = sum(value)) %>%
  spread(key = 'pigments', value = 'value') %>%
  mutate(rel_anth = Anth/sum)%>% ###berechnung der relativen werte in neuer spalte nach Schlueter et al.
  mutate(rel_cantha = Cantha/sum)%>%
  mutate(rel_chl.a = Chl.a/sum)%>%
  mutate(rel_chl.b = Chl.b/sum)%>%
  mutate(rel_chl.c2 = Chl.c2/sum)%>%
  mutate(rel_diato = Diato/sum)%>%
  mutate(rel_echin = Echin/sum)%>%
  mutate(rel_fuco = Fuco/sum)%>%
  mutate(rel_myxo = Myxo/sum)%>%
  mutate(rel_peri = Peri/sum)%>%
  mutate(rel_phe.a = Phe.a/sum)%>%
  mutate(rel_viola = Viola/sum) %>%
  select(-c("Allo", "Anth","bb.Car", "c.Neo", "Cantha","Chl.a","Chl.b","Chl.c1",
           "Chl.c2","Cryp","Diadino", "Diato","Dino","Echin","Fuco", "Lut","Myxo","Peri",       
           "Phe.a","Phe.b", "Viola", "Zea" )) %>%
  gather(key = 'pigment', value = 'value', -planktotron, - sampling, -treatment,-sum) %>%
  group_by(treatment, sampling, pigment) %>% 
  mutate(mean = mean(value, na.rm = T))
  

data %>%
  filter(pigments %in% c('Peri', 'Fuco', 'Chl.b', 'Chl.a', 'Diadino', 'Diato')) %>%
ggplot(., aes(x = sampling, y = value, col = planktotron, group = planktotron))+
  geom_point()+
  geom_line()+
  labs(x = 'sampling', y = 'pigment concentration in ug/L')+
  facet_wrap(~pigments*treatment, scales = 'free')+
  theme_bw()

data%>%
  filter(sampling <14)%>%
  ggplot(aes(x = sampling, y = value, fill = pigment))+
  geom_bar(stat = 'identity')+
  labs(x = 'sampling', y = 'relative_contribution')+
  facet_wrap(~treatment)+
  theme_bw()

data%>%
  filter(sampling %in% c(0,6,10,14,18))%>%
  ggplot(aes(x = treatment, y = value, fill = pigment))+
  geom_bar(stat = 'identity')+
  labs(x = 'sampling', y = 'relative_contribution')+
  facet_wrap(~sampling)+
  theme_bw()
#ggsave(plot = last_plot(), file = 'pigments.png')
#
