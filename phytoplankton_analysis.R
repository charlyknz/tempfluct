### Phytoplankton counts analysis ###
# Charlotte Kunze
# June 29, 2020
library(tidyverse)
library(vegan)


## ARF: This response factor is defined as the biovolume fraction of one
#taxonomic group within the community at t1 and t2, respectively,
#divided by the group’s initial biovolume fraction at t0.


# import dataset
data <- read_delim("~/Desktop/MA/MA_Rcode/project_data/phytoplankton_counts_may2020.csv", 
                   ";", escape_double = FALSE, col_types = cols(date = col_character()), 
                   locale = locale(decimal_mark = ","), 
                   trim_ws = TRUE) 
View(data)
str(data)


#change variable format to numeric 
data$cells_ml <- as.numeric(gsub(",", ".", data$cells_ml))
data1 <- filter(data, !comment %in% 'recount')

#first look at the data
ggplot(data, aes(x = sampling, y = cells_ml, fill = species))+
  geom_col()+
  facet_wrap(~MC)+
  theme_bw()

## merge with treatment information
treatments <- read.csv2('~/Desktop/MA/MA_Rcode/project_data/treatments_units.csv')
str(treatments)
names(treatments) = c('MC', 'fluctuation', 'treatment')
rich <- left_join(data1, treatments, by = c('MC')) 
rich_dat <- rich%>%
  drop_na(treatment)%>%
  group_by(fluctuation, treatment, species, sampling) %>%
  summarise(mean = mean(cells_ml, na.rm = T)) 

ggplot(rich_dat, aes(x = sampling, y = mean, fill = species))+
  geom_col()+
  facet_wrap(~treatment)+
  scale_x_continuous(breaks = c(0,6,10,14,18), limits = c(-3,20))+
  theme_bw()

#only treatment specific
darich <- rich %>%
  drop_na(treatment)%>%
  select(sampling, MC, treatment, fluctuation, species, cells_ml) %>%
  group_by(fluctuation, treatment, species, sampling) %>%
  summarise(mean = mean(cells_ml, na.rm = T)) 
  
#calculate species richness
data_rich <- rich %>%
  drop_na(treatment)%>%
  select(sampling, MC, treatment, fluctuation, species, cells_ml) %>%
  #group_by(fluctuation, treatment, species, sampling) %>%
  #summarise(mean = mean(cells_ml, na.rm = T)) %>%
  #group_by(treatment, sampling, species) %>%
 # mutate(sum_cells_ml = sum(cells_ml, na.rm = T))%>%
  #mutate('Odontella species' = replace_na('Odontella species', 0))%>%
  #select(-cells_ml, -volume, -magnification, -grid_length, -counts) %>%
  ungroup()%>%  
  spread(key = species, value = cells_ml, convert = T) #%>%
  group_by(treatment, sampling) %>%
  mutate(Gymnodinium_all = sum(Gymnodinium, Gymnodinium_Lepodinium)) %>%
  select(-Gymnodinium, -Gymnodinium_Lepodinium) 

#replace nas with 0
data_rich[is.na(data_rich)] <- 0

data_rich$shannon = diversity(data_rich[, -c(1:3)], MARGIN = 1, index='shannon') #new column containing the calculated shannon index

absence_presence <- decostand(data_rich[, -c(1:3)], method= 'pa', na.rm=T) #df giving absence/presence data using decostand function
data_rich$no = apply(absence_presence, MARGIN = 1, FUN = sum) #new column containing the sum of species present per side (by row = MARGIN = 1)

#species richness
richness <- group_by(data_rich, fluctuation, treatment, sampling) %>%
  summarise(richness = mean(no)) 
richness$treatment = factor(as.factor(richness$treatment), levels = c('con', 'F48', 'F36','F24', 'F12','F6'))
#plot for species richness
ggplot(subset(richness, sampling == 10), aes(x = treatment, y = richness, fill = as.factor(treatment)))+
  geom_col()+
  labs(x = 'Treatment', y = 'species richness')+
  theme_bw()
#ggsave(file = 'richenss_plot.png', plot = last_plot())

#stacked barplot with species contributing >5 % to the community composition
stack <- data_rich %>%
  gather(key = 'species', value = 'abundance_ml', -treatment, -sampling, -date, -MC, -fluctuation,-'stripes/grids', -no, -shannon) %>%
  mutate(sum = sum(abundance_ml, na.rm = T),
         rel_ab = (abundance_ml/sum)*100) 
  
ggplot(stack, aes(x = sampling, y = abundance_ml, fill = species))+
  geom_col()+
  #scale_x_continuous(breaks = seq(1,12,1), limits = c(0,13))+
  facet_wrap(~treatment)+
  theme_bw()
#ggsave(plot = last_plot(), file = 'stacked_species.png')
