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


shannon <- rich %>%
  select(-species)  %>%
  group_by(treatment, sampling, genus) %>%
  summarise(mean_cells = mean(cells_ml, na.rm = T)) %>%
  ungroup()%>%
  spread(key = genus, value = mean_cells) %>%
  group_by(treatment, sampling) 
 
#remove NAs and exchange with 0
shannon[is.na(shannon)] <- 0

##calculate shannon diversity index
shannon$shan = diversity(shannon[, -c(1:2)], MARGIN = 1, index='shannon') #new column containing the calculated shannon index

## calculate species richness
absence_presence <- decostand(shannon[, -c(1:2)], method= 'pa', na.rm=T) #df giving absence/presence data using decostand function
shannon$no = apply(absence_presence, MARGIN = 1, FUN = sum) #new column containing the sum of species present per side (by row = MARGIN = 1)

#plot for species richness or shannon
ggplot(shannon, aes(x = sampling, y = shan, col = treatment)) +
  #geom_point()+
  scale_x_continuous(limits = c(-1, 19), breaks = c(0,6,10,14,18))+
  geom_jitter(width = 0.25,size = 2)+
  theme_bw()


#### Biovolume ####
algal_measurement <- read_delim("~/Desktop/MA/MA_Rcode/project_data/algal_measurement.csv", 
                                ";", escape_double = FALSE, locale = locale(decimal_mark = ","), 
                                trim_ws = TRUE)
str(algal_measurement)

#calculate mean biovolume
BioV <- algal_measurement %>%
  group_by(genus, species, id) %>%
  summarise(mean_BioV = mean(BioV, na.rm = T),
            sd = sd(BioV, na.rm = T)) %>%
  ungroup()

## all_data
all_data <- left_join(rich, BioV, by = c('genus', 'species') ) %>%
 select(treatment, fluctuation, MC, sampling, genus, species, id, cells_ml,mean_BioV) %>%
  mutate(volume = cells_ml*mean_BioV) 

#new df to calculate mean biovolume
calc_BV <- all_data%>%
  filter(sampling != 0) %>%
  group_by(treatment, fluctuation, sampling, genus, species, id) %>%
  summarise(mean_V = mean(volume, na.rm = T)) %>% #calculate mean cells with biovolume
  mutate(sp_id = paste(genus,species, sep = '_'))

#calculate mean Volume for sampling 0, assuming all treatments are the same at the beginning
data0 <- all_data %>% 
  filter(sampling == 0) %>%
  mutate(sp_id = paste(genus,species, sep = '_')) %>%
  group_by(genus, species, sp_id, sampling, id) %>%
  summarise(con = mean(volume, na.rm = T),
            F6 = mean(volume, na.rm = T),
            F12 = mean(volume, na.rm = T),
            F24 = mean(volume, na.rm = T),
            F36 = mean(volume, na.rm = T),
            F48 = mean(volume, na.rm = T)) %>%
  gather(key = 'treatment', value = 'mean_V', -genus, -species,-sp_id, -sampling, -id)

##join zero-values and rest
all_BV <- dplyr::bind_rows(calc_BV, data0) %>%
  drop_na(mean_V) %>%
  group_by(id, sampling, treatment) %>%
  summarise(sum = sum(mean_V))

ggplot(all_BV, aes(x = sampling, y = sum, col = treatment))+
  geom_jitter(width = 0.25,size = 2)+
  scale_x_continuous(limits = c(-1, 19), breaks = c(0,6,10,14,18))+
  facet_wrap(~id, scales = 'free_y')+
  theme_bw()
