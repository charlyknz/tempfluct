### Phytoplankton counts analysis ###
# Charlotte Kunze
# June 29, 2020
library(tidyverse)
library(vegan)
########################################################################################

# import dataset
data <- read_delim("~/Desktop/MA/MA_Rcode/project_data/phyto_new.csv", 
                   ";", escape_double = FALSE, col_types = cols(date = col_character()), 
                   locale = locale(decimal_mark = ","), 
                   trim_ws = TRUE) %>%
  drop_na(MC)
View(data)
str(data)


#change variable format to numeric 
data$cells_ml <- as.numeric(gsub(",", ".", data$cells_ml))
#data1 <- filter(data, !comment %in% 'recount')

#first look at the data
ggplot(data, aes(x = sampling, y = cells_ml, col = species))+
  geom_point()+
  facet_wrap(~MC, scales = 'free_y')+
  theme_bw()


## merge with treatment information
treatments <- read.csv2('~/Desktop/MA/MA_Rcode/project_data/treatments_units.csv')
str(treatments)
names(treatments) = c('MC', 'fluctuation', 'treatment')
rich <- left_join(data, treatments, by = c('MC')) %>%
  separate(species, into = c('genus', 'species'), ' ')

########################################################################################
# calculate species abundance over time 
spec <- rich %>%
  mutate(species_id = paste(genus, species, sep = ' '))  %>%
  filter(sampling != 0)%>%
  group_by(treatment, sampling, species_id) %>%
  summarise(mean_cells = mean(cells_ml, na.rm = T), 
            sd = sd(cells_ml, na.rm = T),
            se = sd/sqrt(n())) %>%
  select(-sd)

spec0 <- rich %>% 
  filter(sampling == 0) %>%
  mutate(species_id = paste(genus, species, sep = ' '))  %>%
  group_by( sampling, species_id) %>%
  summarise(con = mean(cells_ml, na.rm = T),
            F6 = mean(cells_ml, na.rm = T),
            F12 = mean(cells_ml, na.rm = T),
            F24 = mean(cells_ml, na.rm = T),
            F36 = mean(cells_ml, na.rm = T),
            F48 = mean(cells_ml, na.rm = T),
            sd_con = sd(cells_ml, na.rm = T),
            sd_F6 = sd(cells_ml, na.rm = T),
            sd_F12 = sd(cells_ml, na.rm = T),
            sd_F24 = sd(cells_ml, na.rm = T),
            sd_F36 = sd(cells_ml, na.rm = T),
            sd_F48 = sd(cells_ml, na.rm = T),
            se_con= sd_con/sqrt(n()),
            se_F6= sd_F6/sqrt(n()),
            se_F12= sd_F12/sqrt(n()),
            se_F24= sd_F24/sqrt(n()),
            se_F36= sd_F36/sqrt(n()),
            se_F48= sd_F48/sqrt(n()))%>%
  select(-sd_con,-sd_F6, -sd_F12, -sd_F24, -sd_F36, -sd_F48)%>%
  gather(key = 'treatment', value = 'mean_cells',  -sampling, -species_id, -se_con,-se_F6, -se_F12, -se_F24, -se_F36, -se_F48) %>%
  gather(key = 'variable', value = 'se', -sampling, -species_id, -treatment,-mean_cells) %>%
  select(-variable) %>%
  distinct(sampling, species_id,treatment, mean_cells, se) %>%
  na.omit()
  
all_spec <- rbind(spec0, spec) %>%
  filter(species_id != 'Rhizosolemia indet.')%>%
  mutate(fluctuation = as.numeric(paste(ifelse(treatment == 'con', 0, 
         ifelse(treatment == 'F48', 48, ifelse(treatment == 'F36', 36, ifelse(treatment == 'F24', 24,
         ifelse(treatment == 'F12',12,6))))))))

all_spec$fluctuation <- factor(as.factor(all_spec$fluctuation),levels=c("0", "48", "36", '24', '12', '6'))

ggplot(subset(all_spec, sampling < 11), aes(x = sampling, y = mean_cells,  group = fluctuation))+
  geom_line(linetype = 2)+
  geom_point(aes(fill = fluctuation), size = 3,col = '#000000', pch = 21)+
  geom_errorbar(aes(ymin = mean_cells - se,  ymax = mean_cells + se), width = .5)+
  scale_fill_manual(values = c( '#000000','#0868ac','#41b6c4','#31a354','#addd8e','#fed976'))+
  labs(x = 'sampling', y= expression(cells~mL^{-1}))+
  facet_wrap(~species_id, scales = 'free_y')+
  theme_bw()+
  theme(legend.position = 'bottom')
#ggsave(plot = last_plot(), file = 'species_abundance_overtime10.png')
########################################################################################
#### indices abundance based ####
#calculate richness, evenness, shannon and simpson based on abundance data

shannon <- rich %>%
  select(treatment, fluctuation,  genus, species, sampling,cells_ml, MC)%>%
  mutate(species_id = paste(genus, species, sep = '_'))%>%
 # group_by(treatment, sampling, species_id) %>%
  #summarise(mean_cells = mean(cells_ml, na.rm = T)) %>%
 select(-genus,-species)%>%
  spread(key = species_id, value = cells_ml) %>%
  #group_by(treatment, sampling) %>%
  mutate(fluctuation = as.numeric(paste(ifelse(treatment == 'con', 0, 
                                               ifelse(treatment == 'F48', 48, ifelse(treatment == 'F36', 36, ifelse(treatment == 'F24', 24,
                                                ifelse(treatment == 'F12',12,6)))))))) %>%
  select(treatment, sampling, fluctuation, everything())
 
#remove NAs and exchange with 0
shannon[is.na(shannon)] <- 0

##calculate shannon diversity index
shannon$shan = diversity(shannon[, -c(1:4)], MARGIN = 1, index='shannon') #new column containing the calculated shannon index

## calculate species richness
shannon <- select(shannon, treatment, sampling, fluctuation, MC, shan, everything())
absence_presence <- decostand(shannon[, -c(1:5)], method= 'pa', na.rm=T) #df giving absence/presence data using decostand function
shannon$no = apply(absence_presence, MARGIN = 1, FUN = sum) #new column containing the sum of species present per side (by row = MARGIN = 1)

#evenness
shannon$evenness = shannon$shan/log(shannon$no)

#simpson
shannon <- select(shannon, treatment, sampling, fluctuation,MC, shan, no, evenness, everything())
shannon$simpson = diversity(shannon[, -c(1:6)], MARGIN = 1, index='simpson') #new column containing the calculated shannon index


#data for export
data_index <- shannon %>%
  select(treatment, sampling,MC, fluctuation, shan, no, evenness,simpson) %>%
  gather(key = 'index', value = 'value', -treatment, -fluctuation,-sampling, -MC) %>%
  group_by(treatment, fluctuation, sampling, index )%>%
  summarise(mean = mean(value, na.rm = T),
            sd = sd(value, na.rm =T),
            se = sd/sqrt(n()))
#write.csv2(x = data_index, file = 'richness_s10.csv')
shannon$fluctuation <- factor(as.factor(shannon$fluctuation),levels=c("0", "48", "36", '24', '12', '6'))
data_index$fluctuation <- factor(as.factor(data_index$fluctuation),levels=c("0", "48", "36", '24', '12', '6'))

#plot for species richness or shannon
ggplot(subset(data_index, sampling == 10), aes(x = fluctuation, y = mean)) +
  geom_point(aes(fill = fluctuation), pch =21, size=3,col = '#000000')+
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = .5)+
  scale_fill_manual(values = c( '#000000','#0868ac','#41b6c4','#31a354','#addd8e','#fed976'))+
  labs(x = 'Fluctuation frequency (in h)', y = 'species richness', fill = 'Fluctuation')+
 # geom_hline(yintercept = 15, linetype = 'dashed', size = 0.5)+
  facet_wrap(~index, scales = 'free_y')+
  theme( panel.background = element_rect(fill = NA), #loescht den Hintergrund meines Plots/ fuellt ihn mit nichts
         #panel.grid.major.y = element_line(color='grey', linetype = 'dashed', size=0.2),
         panel.border= element_rect(colour = "black", fill=NA, size=0.5),
         strip.background = element_rect(color ='black', fill = 'white'),
         strip.text = element_text(face = 'italic'),
         legend.background = element_blank(),
         legend.position  ='bottom',
         legend.key = element_blank(),
         text = element_text(size=12))
#ggsave(plot = last_plot(), file = 'indices_samp10_counting.png')

########################################################################################
########################################################################################
#### Biovolume ####
algal_measurement <- read_delim("~/Desktop/MA/MA_Rcode/project_data/algal_measurement.csv", 
                                ";", escape_double = FALSE, locale = locale(decimal_mark = ","), 
                                trim_ws = TRUE)
str(algal_measurement)

#calculate mean biovolume (mean of measurements)
BioV <- algal_measurement %>%
  group_by(genus, species, id) %>%
  summarise(mean_BioV = mean(BioV, na.rm = T),
            sd = sd(BioV, na.rm = T)) %>%
  ungroup()


#### merge all_data ####
all_data <- left_join(rich, BioV, by = c('genus', 'species') ) %>%
 select(treatment, fluctuation, MC, sampling, genus, species, id, cells_ml,mean_BioV) %>%
 mutate(volume = cells_ml*mean_BioV) #calculate biovolume

#new df to calculate mean biovolume
rel_BV <- all_data%>%
  #filter(sampling ==6)%>%
  mutate(species_id = paste(genus, species, sep = '_'))%>%
  group_by(treatment, fluctuation, sampling, MC)%>%
  mutate(sum = sum(volume, na.rm = T),
         rel_V = volume/sum) %>%
  group_by(treatment, sampling, fluctuation, species_id) %>%
  summarise(mean_contr = mean(rel_V, na.rm = T),
            sd_contr = sd(rel_V, na.rm =T),
            se_contr = sd_contr/ sqrt(n()),
            rel_contr = mean_contr*100,
            se_rel = se_contr*100) %>% #calculate mean cells with biovolume
    filter(mean_contr >0.01)

rel_BV$fluctuation <- factor(as.factor(rel_BV$fluctuation),levels=c("0", "48", "36", '24', '12', '6'))


#plot
ggplot(rel_BV, aes( x = sampling, y = rel_contr))+
  geom_col(aes(fill = species_id), col = 'black')+
  #scale_y_continuous(limits = c(0,100), breaks = seq(0,100, 20))+
  facet_wrap(~treatment)+
  labs(x = 'sampling', y = 'mean contribution of species to BioV', fill = 'treatment')+
  theme( panel.background = element_rect(fill = NA), #loescht den Hintergrund meines Plots/ fuellt ihn mit nichts
         #panel.grid.major.y = element_line(color='grey', linetype = 'dashed', size=0.2),
         panel.border= element_rect(colour = "black", fill=NA, size=0.5),
         strip.background = element_rect(color ='black', fill = 'white'),
         strip.text = element_text(face = 'bold'),
         legend.background = element_blank(),
         legend.position  ='bottom',
         legend.key = element_blank(),
         text = element_text(size=10))
#ggsave(plot=last_plot(), file = 'rel_V_perspecies.png', width = 11, height = 8)
 

########################################################################################
#### indices based on BioV ####

shannon_BV <- all_data %>%
  mutate(species_id = paste(genus, species, sep = '_'))%>%
  select(-genus, -species,-id,-mean_BioV, -cells_ml) %>%
  #group_by(treatment, sampling, species_id) %>%
 # summarise(mean_V = mean(volume, na.rm = T)) %>%
#  ungroup()%>%
  spread(key = species_id, value = volume) %>%
  group_by(treatment, sampling) %>%
  mutate(fluctuation = as.numeric(paste(ifelse(treatment == 'con', 0, 
                                               ifelse(treatment == 'F48', 48, ifelse(treatment == 'F36', 36, ifelse(treatment == 'F24', 24,
                                                                                                                    ifelse(treatment == 'F12',12,6)))))))) %>%
  select(treatment, sampling, fluctuation, everything()) 

#remove NAs and exchange with 0
shannon_BV[is.na(shannon_BV)] <- 0

##calculate shannon diversity index
shannon_BV$shan = diversity(shannon_BV[, -c(1:4)], MARGIN = 1, index='shannon') #new column containing the calculated shannon index
shannon_BV <- select(shannon_BV, treatment, sampling, fluctuation, MC, shan,everything() )
shannon_BV$simpson = diversity(shannon_BV[, -c(1:5)], MARGIN = 1, index='invsimpson') #new column containing the calculated shannon index
shannon_BV <- select(shannon_BV, treatment, sampling, fluctuation, MC, shan, simpson,everything() )
## calculate species richness
absence_presence <- decostand(shannon_BV[, -c(1:6)], method= 'pa', na.rm=T) #df giving absence/presence data using decostand function
shannon_BV$no = apply(absence_presence, MARGIN = 1, FUN = sum) #new column containing the sum of species present per side (by row = MARGIN = 1)

shannon_BV$evenness = shannon_BV$shan/log(shannon_BV$no)

diversity_BV <- shannon_BV %>%
  ungroup() %>%
  select(MC, fluctuation, sampling, evenness, no, shan,simpson) %>%
  gather(key = 'index', value = 'value', -sampling, -fluctuation, -MC) %>%
  group_by(fluctuation, sampling, index) %>%
  summarise(mean_index = mean(value, na.rm = T),
         sd = sd(value, na.rm =T),
         se = sd/ sqrt(n()),
         n = n()) %>%
  mutate(lower.ci.mpg = mean_index - 1.96*se/sqrt(n),
         upper.ci.mpg = mean_index + 1.96*se/sqrt(n),
         day = sampling *2)

diversity_BV$fluctuation <- factor(as.factor(diversity_BV$fluctuation),levels=c("0", "48", "36", '24', '12', '6'))
even_p <- ggplot(subset(diversity_BV, index == 'evenness'), aes(x = day, y = mean_index, group = fluctuation)) +
  geom_line(linetype = 'dashed', size = 0.5, aes(col = fluctuation))+
  geom_errorbar(aes(ymin = lower.ci.mpg , ymax = upper.ci.mpg, col = fluctuation), width = .5, position = position_dodge(width = .9))+
  geom_point(aes(fill = fluctuation), pch = 21,size= 3,col = '#000000', position = position_dodge(width = .9))+
 # scale_x_continuous(limits = c(-1, 20), breaks = seq(0,18,2))+
  scale_fill_manual(values = c( '#000000','#0868ac','#41b6c4','#31a354','#addd8e','#fed976'))+
  scale_color_manual(values = c( '#000000','#0868ac','#41b6c4','#31a354','#addd8e','#fed976'))+
  labs(x = '', y= 'Evenness species composition', fill = 'Fluctuation  \nfrequency [h]', col = 'Fluctuation  \nfrequency [h]')+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),legend.position = 'none', 
        panel.background = element_rect(fill = NA), #loescht den Hintergrund meines Plots/ fuellt ihn mit nichts
        #panel.grid.major.y = element_line(color='grey', linetype = 'dashed', size=0.2),
        panel.border= element_rect(colour = "black", fill=NA, size=1),
        strip.text = element_text(face = 'bold'),
        legend.background = element_blank(),
        legend.key = element_blank(),
        text = element_text(size=17))
even_p

even_no <- ggplot(subset(diversity_BV, index == 'no'), aes(x = day, y = mean_index, group = fluctuation)) +
  geom_line(linetype = 'dashed', size = 0.5, aes(col = fluctuation))+
  geom_errorbar(aes(ymin = lower.ci.mpg , ymax = upper.ci.mpg, col = fluctuation), width = .5, position = position_dodge(width = .9))+
  geom_point(aes(fill = fluctuation), pch = 21,size= 3,col = '#000000', position = position_dodge(width = .9))+
  # scale_x_continuous(limits = c(-1, 20), breaks = seq(0,18,2))+
  scale_fill_manual(values = c( '#000000','#0868ac','#41b6c4','#31a354','#addd8e','#fed976'))+
  scale_color_manual(values = c( '#000000','#0868ac','#41b6c4','#31a354','#addd8e','#fed976'))+
  labs(x = 'Time [days]', y= 'Species richness', fill = 'Fluctuation  \nfrequency [h]', col = 'Fluctuation  \nfrequency [h]')+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),legend.position = 'none', 
        panel.background = element_rect(fill = NA), #loescht den Hintergrund meines Plots/ fuellt ihn mit nichts
        #panel.grid.major.y = element_line(color='grey', linetype = 'dashed', size=0.2),
        panel.border= element_rect(colour = "black", fill=NA, size=1),
        strip.text = element_text(face = 'bold'),
        legend.background = element_blank(),
        legend.key = element_blank(),
        text = element_text(size=17))
even_no

library(cowplot)
plot_grid(even_p, even_no,  labels=c("(a)","(b)"),ncol = 2, label_size = 18, hjust = 0, vjust = 0.95)
#ggsave(plot = last_plot(), file = 'even_no_specBV.png', width = 9, height = 5)



#ggsave(plot = last_plot(), file = 'shannon_div_BioV_time.png')  
#####################################################################################

#### STATS ####

diversity <- shannon_BV %>%
  ungroup() %>%
  select(MC, fluctuation, sampling, evenness, no, shan,simpson) %>%
  mutate(day = 2*sampling,
         dayname = as.factor(day)) %>%
  mutate(fluctuation = as.numeric(fluctuation),
         interval = 48/fluctuation)
diversity$interval[!is.finite(diversity$interval)] <- 0 
even1 <- lmer(evenness ~ interval*day + (1|MC) + (1|dayname), data=diversity)
summary(even1)
anova(even1)  
  
  
  