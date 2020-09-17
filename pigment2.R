##
# R Skript to analyse pigment composition
library(tidyverse)
library(vegan)
library(cowplot)
## load data

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
#approach to calculate an overall diversity index 
data_try <- data %>%
  select(-X, -no) 
#remove NAs and exchange with 0
data_try[is.na(data_try)] <- 0
data_try$shannon = diversity(data_try[, -c(1:4)], MARGIN = 1, index='shannon') #new column containing the calculated shannon index
## calculate species richness
data_try <- select(data_try, treatment, sampling, date, planktotron, shannon, everything())
ab_presence <- decostand(data_try[, -c(1:5)], method= 'pa', na.rm=T) #df giving absence/presence data using decostand function
data_try$richness = apply(ab_presence, MARGIN = 1, FUN = sum) #new column containing the sum of species present per side (by row = MARGIN = 1)

data_try$evenness = data_try$shannon/log(data_try$richness)

data_try1 <- data_try %>%
  gather(key = 'pigment', value = 'value', -date,-planktotron,-sampling,-treatment, -shannon, -evenness, -richness) %>%
  gather(key = 'index', value = 'values', -date,-planktotron, -sampling,-treatment, -pigment,-treatment,-value) %>%
  group_by(treatment, index) %>%
  summarise(mean = mean(values, na.rm = T),
            sd = sd(values, na.rm = T),
            se = sd/sqrt(n()))%>%
  mutate(fluctuation = paste(ifelse(treatment == 'control', 0, ifelse(treatment == 'Fluctuating_48', 48, ifelse(treatment == 'Fluctuating_36', 36,
                                                                                                                                            ifelse(treatment == 'Fluctuating_24', 24, ifelse(treatment == 'Fluctuating_12', 12, 6)))))))
data_try1$fluctuation <- factor(as.factor(data_try1$fluctuation),levels=c("0", "48", "36", '24', '12', '6'))
ggplot(data_try1, aes( x = fluctuation, y = mean, fill = as.factor(fluctuation))) +
  geom_col()+
  geom_errorbar(aes(ymin = mean -se, ymax = mean +se), width = .2)+
  scale_fill_manual(values = c( '#000000','#0868ac','#41b6c4', '#31a354','#addd8e','#fed976'))+
  labs( x = 'Fluctuation frequency (in h)', y = 'shannon H index (pigment)', fill = 'Fluctuation')+
  facet_wrap(~index, scales = 'free_y')+
  theme_bw()
#ggsave(plot = last_plot(), file = 'accumulated_indices_pigments.png')
##################################################################################################################################

### calculate  diversity indices of pigment composition
data_diversity <- data 

#remove NAs and exchange with 0
data_diversity[is.na(data_diversity)] <- 0

##calculate shannon diversity index
data_diversity$shannon = diversity(data_diversity[, -c(1:6)], MARGIN = 1, index='shannon') #new column containing the calculated shannon index
data_diversity <- select(data_diversity, X, no, planktotron, sampling,date, treatment, shannon, everything())
data_diversity$simpson = diversity(data_diversity[, -c(1:7)], MARGIN = 1, index= 'invsimpson') #new column containing the calculated shannon index
data_diversity <- select(data_diversity, X, no, planktotron, sampling,date, treatment, shannon, simpson, everything())
## calculate species richness

ab_presence <- decostand(data_diversity[, -c(1:8)], method= 'pa', na.rm=T) #df giving absence/presence data using decostand function
data_diversity$rich = apply(ab_presence, MARGIN = 1, FUN = sum) #new column containing the sum of species present per side (by row = MARGIN = 1)

data_diversity$evenness = data_diversity$shannon/log(data_diversity$rich)

#### LMM on diversity ####
modell_dat <- data_diversity %>%
  dplyr::select(planktotron, sampling, treatment,evenness, simpson) %>%
  separate(treatment, into = c('mist', 'fluctuation'), '_') %>%
  mutate(day = sampling *2,
         dayname = as.factor(day)) %>%
  mutate(fluctuation = as.numeric(fluctuation),
         interval = 48/fluctuation)

modell_dat$interval[is.na(modell_dat$interval)] <-0
modell_dat$fluctuation[is.na(modell_dat$fluctuation)] <-0
div <- lmer(evenness ~ interval*day + (1|planktotron) + (1|dayname), data=modell_dat)
summary(div)
anova(div)

#data wrangling
data_plot  <- data_diversity %>%
  select(sampling, treatment,planktotron,shannon, evenness, rich, simpson) %>%
  gather(key = 'index', value = 'value', -treatment, -sampling, -planktotron) %>%
  dplyr::group_by(sampling, treatment, index) %>%
  dplyr::summarise(mean_index = mean(value, na.rm = T),
            sd.mpg = sd(value, na.rm = T),
            n.mpg = dplyr::n()) %>%
  mutate(se_index = sd.mpg/sqrt(n.mpg),
         lower.ci.mpg = mean_index - 1.96*se_index/sqrt(n.mpg),
         upper.ci.mpg = mean_index + 1.96*se_index/sqrt(n.mpg)) %>%
  separate(treatment, into = c('mist', 'fluctuation'), '_') %>%
  mutate(day = sampling *2) %>%
  group_by(mean_index) %>%
  mutate(width = 0.5 * n())
  data_plot$fluctuation[is.na(data_plot$fluctuation)] <- 0

#plot
data_plot$fluctuation <- factor(as.factor(data_plot$fluctuation),levels=c("0", "48", "36", '24', '12', '6'))

simpson <- ggplot(subset(data_plot, index %in% c('simpson')), aes(x =day, y = mean_index))+
  geom_line(linetype = 'dashed', size = 0.5, aes(color = fluctuation))+
  #geom_errorbar(aes(ymin = mean_index - se_index, ymax = mean_index+se_index, color = fluctuation), width = .5)+
  geom_errorbar(aes(ymin = lower.ci.mpg, ymax = upper.ci.mpg, color = fluctuation),position = position_dodge(width = .9),width = .6, size=.5)+
  geom_point(aes(fill = fluctuation), pch=21, size=3, col = '#000000', position = position_dodge(width = .9))+
  scale_fill_manual(values = c( '#000000','#0868ac','#41b6c4','#31a354','#addd8e','#fed976'))+
  scale_color_manual(values = c( '#000000','#0868ac','#41b6c4','#31a354','#addd8e','#fed976'))+
  labs(x = 'Time [days]', y= 'simpson index of pigment diversity', fill = 'Fluctuation  \nfrequency [h]', col = 'Fluctuation  \nfrequency [h]')+
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

even <- ggplot(subset(data_plot, index %in% c('evenness')), aes(x =day, y = mean_index))+
  geom_line(linetype = 'dashed', size = 0.5, aes(color = fluctuation))+
 # geom_errorbar(aes(ymin = mean_index - se_index, ymax = mean_index+se_index), width = .5)+
  geom_errorbar(aes(ymin = lower.ci.mpg, ymax = upper.ci.mpg, color = fluctuation), width = .6,position = position_dodge2(width = .6))+
  geom_point(aes(fill = fluctuation), pch=21, size=3, col = '#000000',position = position_dodge2(width = .6))+
  scale_fill_manual(values = c( '#000000','#0868ac','#41b6c4','#31a354','#addd8e','#fed976'))+
  scale_color_manual(values = c( '#000000','#0868ac','#41b6c4','#31a354','#addd8e','#fed976'))+
  labs(x = ' ', y= 'evenness of pigment diversity', fill = 'Fluctuation  \nfrequency [h]', col = 'Fluctuation  \nfrequency [h]')+
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

plot_grid(even, simpson,  labels=c("(a)","(b)"),ncol = 1, label_size = 18, hjust = 0, vjust = 0.95)
#ggsave(plot = last_plot(), file = 'pigments.png', width = 9, height = 9)

#######

ci <- data_plot %>%
  dplyr::select(day, fluctuation, index, mean_index, lower.ci.mpg, upper.ci.mpg) %>%
  arrange(index, day)
#write.csv2(x = ci, file = 'ci_pigments.csv')
