### Phytoplankton counts analysis ###
# Charlotte Kunze
# June 29, 2020
library(tidyverse)
library(vegan)


## ARF: This response factor is defined as the biovolume fraction of one
#taxonomic group within the community at t1 and t2, respectively,
#divided by the group’s initial biovolume fraction at t0.


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
  select(-species)  %>%
  group_by(treatment, sampling, genus) %>%
  summarise(mean_cells = mean(cells_ml, na.rm = T)) %>%
  ungroup()%>%
  spread(key = genus, value = mean_cells) %>%
  group_by(treatment, sampling) %>%
  mutate(fluctuation = as.numeric(paste(ifelse(treatment == 'con', 0, 
                                               ifelse(treatment == 'F48', 48, ifelse(treatment == 'F36', 36, ifelse(treatment == 'F24', 24,
                                                ifelse(treatment == 'F12',12,6)))))))) %>%
  select(treatment, sampling, fluctuation, everything())
 
#remove NAs and exchange with 0
shannon[is.na(shannon)] <- 0

##calculate shannon diversity index
shannon$shan = diversity(shannon[, -c(1:3)], MARGIN = 1, index='shannon') #new column containing the calculated shannon index

## calculate species richness
absence_presence <- decostand(shannon[, -c(1:3)], method= 'pa', na.rm=T) #df giving absence/presence data using decostand function
shannon$no = apply(absence_presence, MARGIN = 1, FUN = sum) #new column containing the sum of species present per side (by row = MARGIN = 1)

#write.csv2(x = shannon, file = 'diversity_indices.csv')

data_index <- shannon %>%
  select(-c(4:23)) %>%
  filter(sampling == 10)
#write.csv2(x = data_index, file = 'richness_s10.csv')
shannon$fluctuation <- factor(as.factor(shannon$fluctuation),levels=c("0", "48", "36", '24', '12', '6'))

#plot for species richness or shannon
ggplot(subset(shannon, sampling == 10), aes(x = fluctuation, y = no)) +
  geom_col(aes(fill = fluctuation), col = '#000000')+
  scale_y_continuous(limits = c(-1, 20), breaks = seq(0,18,3))+
  scale_fill_manual(values = c( '#000000','#0868ac','#41b6c4','#31a354','#addd8e','#fed976'))+
  labs(x = 'Fluctuation frequency (in h)', y = 'species richness', fill = 'Fluctuation')+
  geom_hline(yintercept = 16, linetype = 'dashed', size = 0.5)+
  theme( panel.background = element_rect(fill = NA), #loescht den Hintergrund meines Plots/ fuellt ihn mit nichts
         #panel.grid.major.y = element_line(color='grey', linetype = 'dashed', size=0.2),
         panel.border= element_rect(colour = "black", fill=NA, size=0.5),
         strip.background = element_rect(color ='black', fill = 'white'),
         strip.text = element_text(face = 'italic'),
         legend.background = element_blank(),
         legend.position  ='bottom',
         legend.key = element_blank(),
         text = element_text(size=12))
#ggsave(plot = last_plot(), file = 'richness.png')
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

ggplot(rel_BV, aes( x = fluctuation, y = mean_contr))+
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
  