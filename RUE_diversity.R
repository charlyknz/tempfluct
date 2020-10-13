#### new script to analyse RUE diversity relationship ####
# by Charlotte Kunze 21.9.2020

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
  dplyr::select(fluctuation, sampling, planktotron, P_RUE, N_RUE) %>%
  gather(key = 'RUE', value = 'value', -fluctuation, -sampling, -planktotron)

RUE12$fluctuation <- factor(as.factor(RUE12$fluctuation),levels=c("0", "48", "36", '24', '12', '6'))
names(RUE12)

#### pigment diversity ####

data <- read.csv2("~/Desktop/MA/MA_Rcode/project_data/master_pigments.csv", sep = ";", dec = ',')
names(data)

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
  gather(key = 'index', value = 'value_div', -fluctuation, -sampling, -planktotron) 
names(data_plot)


all_RUE_div <- left_join(RUE12, data_plot, by = c('sampling', 'fluctuation', 'planktotron'))

all_RUE_div$fluctuation <- factor(as.factor(all_RUE_div$fluctuation),levels=c("0", "48", "36", '24', '12', '6'))

levels(as.factor(all_RUE_div$fluctuation))

##########################################################
library(ggpubr)

formula <- y ~ x
RUE_div <- ggplot(subset(all_RUE_div, RUE == 'P_RUE' & index == 'evenness'), aes(x = value_div, y = value))+
  geom_point(aes(color = fluctuation), size = 3)+
  scale_color_manual(values = c( '#000000','#0868ac','#41b6c4','#31a354','#addd8e','#fed976'))+
  labs(x = 'evenness', y = expression(RUE['P']))+
 # stat_smooth(aes(fill = fluctuation, color = fluctuation), method = "lm", se = F, size = 0.5,formula = formula) +
  stat_smooth( method = "lm", se = T, size = 0.5,formula = formula, color = 'black') +
  stat_regline_equation(formula = formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~"))) +     
  theme( panel.background = element_rect(fill = NA), #loescht den Hintergrund meines Plots/ fuellt ihn mit nichts
         #panel.grid.major.y = element_line(color='grey', linetype = 'dashed', size=0.2),
         panel.border= element_rect(colour = "black", fill=NA, size=1),
         strip.background = element_rect(color ='black', fill = 'white'),
         strip.text = element_text(face = 'bold'),
         legend.background = element_blank(),
         legend.position  ='none',
         legend.key = element_blank(),
         text = element_text(size=18))
RUE_div
#ggsave(plot = RUE_div, file = 'RUE_P_evenness_1lm.tiff', width = 5, height = 4)
RUE_divN <- ggplot(subset(all_RUE_div, RUE == 'N_RUE' & index == 'evenness'), aes(x = value_div, y = value))+
  geom_point(aes(fill = fluctuation), pch = 21, color = 'black',size = 3)+
  scale_fill_manual(values = c( '#000000','#0868ac','#41b6c4','#31a354','#addd8e','#fed976'))+
  labs(x = 'pigment evenness', y =expression(RUE['N']))+
  #stat_smooth(aes(fill = fluctuation, color = fluctuation), method = "lm", se = F, size = 0.5,formula = formula) +
  #stat_regline_equation(formula = formula, 
   #                  aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")))+     
  #stat_smooth( method = "lm", se = T, size = 0.5,formula = formula, color = 'black') +
  theme( panel.background = element_rect(fill = NA), #loescht den Hintergrund meines Plots/ fuellt ihn mit nichts
         #panel.grid.major.y = element_line(color='grey', linetype = 'dashed', size=0.2),
         panel.border= element_rect(colour = "black", fill=NA, size=1),
         strip.background = element_rect(color ='black', fill = 'white'),
         strip.text = element_text(face = 'bold'),
         legend.background = element_blank(),
         legend.position  ='none',
         legend.key = element_blank(),
         text = element_text(size=18))
RUE_divN
ggsave(RUE_divN, file = 'RUE_N_diversity.tiff', width = 5, height = 4)
plot_grid(RUE_div, RUE_divN, labels = c('(a)', '(b)', label_size = 18))
#ggsave(plot = last_plot(), file = 'even_RUE_overall_lm.png', width = 9, height = 5)


ggscatter(subset(all_RUE_div, RUE == 'N_RUE' & index == 'evenness'), x = "value_div", y = "value", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman", xlab = 'pigment evenness',
          ylab = 'RUE (Total N)')
ggsave(plot = last_plot(), 'scatter_RUE_Even_plot.png')
