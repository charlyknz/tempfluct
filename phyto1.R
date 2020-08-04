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
ggplot(data, aes(x = sampling, y = cells_ml, fill = species))+
  geom_col()+
  facet_wrap(~MC)+
  theme_bw()

########################################################################################

## merge with treatment information
treatments <- read.csv2('~/Desktop/MA/MA_Rcode/project_data/treatments_units.csv')
str(treatments)
names(treatments) = c('MC', 'fluctuation', 'treatment')
rich <- left_join(data, treatments, by = c('MC')) %>%
  separate(species, into = c('genus', 'species'), ' ')

spec <- rich %>%
  mutate(species_id = paste(genus, species, sep = ' '))  %>%
  group_by(treatment, sampling, species_id) %>%
  summarise(mean_cells = mean(cells_ml, na.rm = T), 
            sd = sd(cells_ml, na.rm = T),
            se = sd/sqrt(n())) %>%
  ggplot(., aes(x = treatment, y = mean_cells, fill = species_id))+
  geom_col()+
  facet_wrap(~sampling)+
  theme_bw()
spec
  
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
            rel_contr = mean_contr*100 ) %>% #calculate mean cells with biovolume
    filter(mean_contr >0.01)

rel_BV$fluctuation <- factor(as.factor(rel_BV$fluctuation),levels=c("0", "48", "36", '24', '12', '6'))


#plot
ggplot(rel_BV, aes( x = fluctuation, y = rel_contr))+
  geom_col(aes(fill = species_id), col = 'black')+
  #scale_y_continuous(limits = c(0,100), breaks = seq(0,100, 20))+
  facet_wrap(~sampling)+
  labs(x = 'Fluctuation frequency (in H)', y = 'mean contribution of species to BioV', fill = 'species')+
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
 

#####################################################################
#shannon based on BioV

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
shannon_BV$simpson = diversity(shannon_BV[, -c(1:5)], MARGIN = 1, index='simpson') #new column containing the calculated shannon index
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
  mutate(mean_index = mean(value, na.rm = T),
         sd = sd(value, na.rm =T),
         se = sd/ sqrt(n()))

diversity_BV$fluctuation <- factor(as.factor(diversity_BV$fluctuation),levels=c("0", "48", "36", '24', '12', '6'))
ggplot(subset(diversity_BV, sampling == 10), aes(x = fluctuation, y = mean_index)) +
  geom_line(linetype = 'dashed', size = 0.5)+
  geom_point(aes(fill = fluctuation), pch = 21,size= 3,col = '#000000')+
  geom_errorbar(aes(ymin = mean_index - se, ymax = mean_index + se), width = .5)+
 # scale_x_continuous(limits = c(-1, 20), breaks = seq(0,18,2))+
  scale_fill_manual(values = c( '#000000','#0868ac','#41b6c4','#31a354','#addd8e','#fed976'))+
  labs(x = 'Fluctuation frequency (in h)', y = 'shannon diversity', fill = 'Fluctuation')+
  theme_bw()+
  facet_wrap(~index, scales = 'free_y')
  
#ggsave(plot = last_plot(), file = 'shannon_div_BioV_time.png')  
#####################################################################################
