##
# R Skript to analyse pigment composition
library(tidyverse)
library(vegan)
library(cowplot)
## load data


#### Data exploration ####
##import master data
data <- read.csv2("~/Desktop/MA/MA_Rcode/project_data/master_pigments.csv", sep = ";", dec = ',')
names(data)
pig_data <- data %>%
  gather(key = 'pigments',value = 'value',-X, -no, -date,-treatment, -sampling, -planktotron )%>%
  group_by(treatment,sampling, planktotron)%>%
  mutate(sum = sum(value)) %>%
  ungroup()%>%
  spread(pigments, value ) %>%
  mutate(rel_allo = Allo/sum)%>%
  mutate(rel_anth = Anth/sum)%>% ###berechnung der relativen werte in neuer spalte nach Schlueter et al.
  mutate(rel_bb.Car = bb.Car/sum)%>%
  mutate(rel_cNeo = c.Neo/sum) %>%
  mutate(rel_cantha = Cantha/sum)%>%
  mutate(rel_chl.a = Chl.a/sum)%>%
  mutate(rel_chl.b = Chl.b/sum)%>%
  mutate(rel_chl.c1 = Chl.c1/sum)%>%
  mutate(rel_chl.c2 = Chl.c2/sum)%>%
  mutate(rel_diad = Diadino/sum)%>%
  mutate(rel_diato = Dino/sum)%>%
  mutate(rel_dino = Diato/sum)%>%
  mutate(rel_echin = Echin/sum)%>%
  mutate(rel_fuco = Fuco/sum)%>%
  mutate(rel_lut = Lut/sum)%>%
  mutate(rel_myxo = Myxo/sum)%>%
  mutate(rel_peri = Peri/sum)%>%
  mutate(rel_phe.a = Phe.a/sum)%>%
  mutate(rel_viola = Viola/sum) %>%
  mutate(rel_zea = Zea/sum)%>%
  dplyr::select(-c("X","Allo", "Anth","bb.Car", "c.Neo", "Cantha","Chl.a","Chl.c1",
           "Chl.c2","Cryp","Diadino", "Diato","Dino","Echin","Lut","Myxo","Peri",       
           "Phe.a","Phe.b", "Viola", "Zea" , 'Fuco', 'Chl.b')) %>%
  gather(key = 'pigment', value = 'value', -planktotron, - sampling, -treatment,-sum, -no, -date) %>%
  group_by(treatment, sampling, pigment) %>% 
  summarise(mean = mean(value, na.rm = T))%>%
  mutate(fluctuation = paste(ifelse(treatment == 'control', 0, ifelse(treatment == 'Fluctuating_48', 48, ifelse(treatment == 'Fluctuating_36', 36,
                                                                                                                ifelse(treatment == 'Fluctuating_24', 24, ifelse(treatment == 'Fluctuating_12', 12, 6)))))))
pig_data$fluctuation <- factor(as.factor(pig_data$fluctuation),levels=c("0", "48", "36", '24', '12', '6'))

  

pig_data %>%
  filter(pigment %in% c('rel_fuco', 'rel_chl.b')) %>%
ggplot(., aes(x = sampling, y = mean, col = treatment, group = treatment))+
  geom_point()+
  geom_line()+
  labs(x = 'sampling', y = 'pigment concentration in ug/L')+
  facet_wrap(~treatment*pigment)+
  theme_bw()
#ggsave(plot = last_plot(), file = 'mean_pigments.png')
#
pig_data%>%
  filter(sampling %in% c(0,6,10,14,18))%>%
  ggplot(aes(x = fluctuation, y = mean, fill = pigment))+
  geom_col(col = 'black')+
  labs(x = 'sampling', y = 'relative_contribution of pigment abundance')+
  facet_wrap(~sampling)+
  theme_bw()+
  theme(legend.position = 'bottom')
#ggsave(plot = last_plot(), file = 'rel_pigments.png')

##################################################################################################################################
##################################################################################################################################
#### diversity for relative data ####
pigData <- data %>%
  gather(key = 'pigments',value = 'value',-X, -no, -date,-treatment, -sampling, -planktotron )%>%
  group_by(treatment,sampling, planktotron)%>%
  mutate(sum = sum(value)) %>%
  ungroup()%>%
  spread(pigments, value ) %>%
  mutate(rel_allo = Allo/sum)%>%
  mutate(rel_anth = Anth/sum)%>% ###berechnung der relativen werte in neuer spalte nach Schlueter et al.
  mutate(rel_bb.Car = bb.Car/sum)%>%
  mutate(rel_cNeo = c.Neo/sum) %>%
  mutate(rel_cantha = Cantha/sum)%>%
  mutate(rel_chl.a = Chl.a/sum)%>%
  mutate(rel_chl.b = Chl.b/sum)%>%
  mutate(rel_chl.c1 = Chl.c1/sum)%>%
  mutate(rel_chl.c2 = Chl.c2/sum)%>%
  mutate(rel_diad = Diadino/sum)%>%
  mutate(rel_diato = Dino/sum)%>%
  mutate(rel_dino = Diato/sum)%>%
  mutate(rel_echin = Echin/sum)%>%
  mutate(rel_fuco = Fuco/sum)%>%
  mutate(rel_lut = Lut/sum)%>%
  mutate(rel_myxo = Myxo/sum)%>%
  mutate(rel_peri = Peri/sum)%>%
  mutate(rel_phe.a = Phe.a/sum)%>%
  mutate(rel_viola = Viola/sum) %>%
  mutate(rel_zea = Zea/sum)%>%
  dplyr::select(-c("X","Allo", "Anth","bb.Car", "c.Neo", "Cantha","Chl.a","Chl.c1",
                   "Chl.c2","Cryp","Diadino", "Diato","Dino","Echin","Lut","Myxo","Peri",       
                   "Phe.a","Phe.b", "Viola", "Zea" , 'Fuco', 'Chl.b')) %>%
  gather(key = 'pigment', value = 'value', -planktotron, - sampling, -treatment,-sum, -no, -date) %>%
  mutate(fluctuation = paste(ifelse(treatment == 'control', 0, ifelse(treatment == 'Fluctuating_48', 48, ifelse(treatment == 'Fluctuating_36', 36,
                                                                                                                ifelse(treatment == 'Fluctuating_24', 24, ifelse(treatment == 'Fluctuating_12', 12, 6)))))))
pigData$sampling = as.numeric(pigData$sampling)
pigData[is.na(pigData)] <-0 #make sure we don't have NAs

names(pigData)
# 1. decide which variables are active, remove character variables and turn them into row names 
pigment_data.active <- pigData %>% 
  ungroup()%>%
  dplyr::select(-treatment) %>% #remove columns which are not needed anymore
  spread(key = pigment, value = value) %>%
  dplyr::select(-rel_bb.Car, -date, -sum, -no) #remove variable with only 0 entries

pigment_data.active[is.na(pigment_data.active)] <- 0 #substitute NAs with 0


### calculate  diversity indices of pigment composition
diversity <- pigment_data.active 

#remove NAs and exchange with 0
diversity[is.na(diversity)] <- 0

##calculate shannon diversity index
diversity$shannon = diversity(diversity[, -c(1:4)], MARGIN = 1, index='shannon') #new column containing the calculated shannon index
diversity <- select(diversity, planktotron, sampling,fluctuation, shannon, everything())
diversity$simpson = diversity(diversity[, -c(1:5)], MARGIN = 1, index= 'invsimpson') #new column containing the calculated shannon index
diversity <- select(diversity, planktotron, sampling,fluctuation, shannon, simpson, everything())
## calculate species richness

ab_presence <- decostand(diversity[, -c(1:6)], method= 'pa', na.rm=T) #df giving absence/presence data using decostand function
diversity$rich = apply(ab_presence, MARGIN = 1, FUN = sum) #new column containing the sum of species present per side (by row = MARGIN = 1)

diversity$evenness = diversity$shannon/log(diversity$rich)

names(diversity)
#data wrangling
data_plot  <- diversity %>%
  select(sampling, fluctuation,planktotron,shannon, evenness, rich, simpson) %>%
  gather(key = 'index', value = 'value', -fluctuation, -sampling, -planktotron) %>%
  dplyr::group_by(sampling, fluctuation, index) %>%
  dplyr::summarise(mean_index = mean(value, na.rm = T),
            sd.mpg = sd(value, na.rm = T),
            n.mpg = dplyr::n()) %>%
  mutate(se_index = sd.mpg/sqrt(n.mpg),
         lower.ci.mpg = mean_index - 1.96*se_index/sqrt(n.mpg),
         upper.ci.mpg = mean_index + 1.96*se_index/sqrt(n.mpg)) %>%
  mutate(day = sampling *2) %>%
  group_by(mean_index) %>%
  mutate(width = 0.5 * n())

#plot
data_plot$fluctuation <- factor(as.factor(data_plot$fluctuation),levels=c("0", "48", "36", '24', '12', '6'))

simpson <- ggplot(subset(data_plot, index %in% c('simpson')), aes(x =day, y = mean_index))+
  geom_line(linetype = 'dashed', size = 0.5, aes(color = fluctuation))+
  #geom_errorbar(aes(ymin = mean_index - se_index, ymax = mean_index+se_index, color = fluctuation), width = .5)+
  geom_errorbar(aes(ymin = lower.ci.mpg, ymax = upper.ci.mpg, color = fluctuation),position = position_dodge(width = .9),width = .6, size=.5)+
  geom_point(aes(fill = fluctuation), pch=21, size=3, col = '#000000', position = position_dodge(width = .9))+
  scale_fill_manual(values = c( '#000000','#0868ac','#41b6c4','#31a354','#addd8e','#fed976'))+
  scale_color_manual(values = c( '#000000','#0868ac','#41b6c4','#31a354','#addd8e','#fed976'))+
  labs(x = 'Time [Days]', y= 'Simpson diversity index', fill = 'Fluctuation  \nfrequency [h]', col = 'Fluctuation  \nfrequency [h]')+
  theme_classic()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),legend.position = 'none', 
        panel.background = element_rect(fill = NA), #loescht den Hintergrund meines Plots/ fuellt ihn mit nichts
        #panel.grid.major.y = element_line(color='grey', linetype = 'dashed', size=0.2),
        panel.border= element_rect(colour = "black", fill=NA, size=1),
        strip.text = element_text(face = 'bold'),
        legend.background = element_blank(),
        legend.key = element_blank(),
        text = element_text(size=17))#
simpson 
#ggsave(plot = last_plot(), file = 'pigment_simpson.png', width = 7, height = 6)

even <- ggplot(subset(data_plot, index %in% c('evenness') & fluctuation %in% c(0, 48,36,24,12,6)), aes(x =day, y = mean_index))+
  geom_line(linetype = 'dashed', size = 0.5, aes(color = fluctuation))+
 # geom_errorbar(aes(ymin = mean_index - se_index, ymax = mean_index+se_index), width = .5)+
  geom_errorbar(aes(ymin = lower.ci.mpg, ymax = upper.ci.mpg, color = fluctuation), width = .6,position = position_dodge2(width = .6))+
  geom_point(aes(fill = fluctuation), pch=21, size=3, col = '#000000',position = position_dodge2(width = .6))+
  scale_fill_manual(values = c( '#000000','#0868ac','#41b6c4','#31a354','#addd8e','#fed976'))+
  scale_color_manual(values = c( '#000000','#0868ac','#41b6c4','#31a354','#addd8e','#fed976'))+
  scale_y_continuous(limits = c(0.4, 1), breaks = seq(0.5,1,0.25))+
  labs(x = '', y= 'Evenness', fill = 'Fluctuation  \nfrequency [h]', col = 'Fluctuation  \nfrequency [h]')+
  theme_classic()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),legend.position = 'none', 
        panel.background = element_rect(fill = NA), #loescht den Hintergrund meines Plots/ fuellt ihn mit nichts
        #panel.grid.major.y = element_line(color='grey', linetype = 'dashed', size=0.2),
        panel.border= element_rect(colour = "black", fill=NA, size=1),
        strip.text = element_text(face = 'bold'),
        legend.background = element_blank(),
        legend.key = element_blank(),
        text = element_text(size=17))#
even 
#ggsave(plot = even, file = 'evenness_0483624126.png', width = 5, height = 4)



rich <- ggplot(subset(data_plot, index %in% c('rich')), aes(x =day, y = mean_index))+
  geom_line(linetype = 'dashed', size = 0.5, aes(color = fluctuation))+
  # geom_errorbar(aes(ymin = mean_index - se_index, ymax = mean_index+se_index), width = .5)+
  geom_errorbar(aes(ymin = lower.ci.mpg, ymax = upper.ci.mpg, color = fluctuation), width = .6,position = position_dodge2(width = .6))+
  geom_point(aes(fill = fluctuation), pch=21, size=3, col = '#000000',position = position_dodge2(width = .6))+
  scale_fill_manual(values = c( '#000000','#0868ac','#41b6c4','#31a354','#addd8e','#fed976'))+
  scale_color_manual(values = c( '#000000','#0868ac','#41b6c4','#31a354','#addd8e','#fed976'))+
  labs(x = 'Time [days]', y= 'Richness', fill = 'Fluctuation  \nfrequency [h]', col = 'Fluctuation  \nfrequency [h]')+
  theme_classic()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),legend.position = 'none', 
        panel.background = element_rect(fill = NA), #loescht den Hintergrund meines Plots/ fuellt ihn mit nichts
        #panel.grid.major.y = element_line(color='grey', linetype = 'dashed', size=0.2),
        panel.border= element_rect(colour = "black", fill=NA, size=1),
        strip.text = element_text(face = 'bold'),
        legend.background = element_blank(),
        legend.key = element_blank(),
        text = element_text(size=17))#
rich 

#ggsave(plot = rich, file = 'richness of pigment composition.tiff', width = 5, height = 4)
plot_grid(even, rich,  labels=c("(a)","(b)"),ncol = 1, label_size = 18, hjust = 0, vjust = 0.95)
#ggsave(plot = last_plot(), file = 'pigments.tiff', width = 9, height = 9)

#######



