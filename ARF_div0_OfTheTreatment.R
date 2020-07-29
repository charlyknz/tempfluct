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
  mutate(volume = cells_ml*mean_BioV)  #calculate biovolume

#Schritt 1: new df to calculate mean biovolume
calc_BV <- all_data%>%
  filter(sampling == 10) %>% #remove sampling 0 
  group_by(treatment, fluctuation, sampling, id) %>%
  summarise(mean_V = mean(volume, na.rm = T),
            sd_V = sd(volume, na.rm = T),
            se_V = sd_V/sqrt(n())) %>% #calculate mean cells with biovolume
  ungroup()%>%
  select(-sd_V, -sampling)%>%
  drop_na(mean_V) #mean_V is the mean of two replicate MC

# Schritt 2: subset of sampling 0
ar0 <- all_data %>%
  filter(sampling %in% c(0)) %>% #remove sampling 0 
  group_by(treatment, fluctuation, sampling, id) %>%
  summarise(mean_V0 = mean(volume, na.rm = T),
            sd_V0 = sd(volume, na.rm = T),
            se_V0 = sd_V0/sqrt(n())) %>% #calculate mean cells with biovolume
  ungroup()%>%
  select(-sd_V0, -sampling)%>%
  drop_na(mean_V0)

# Schritt 3. merge 0 and normal values again
BV_each <- left_join(calc_BV, ar0, by = c('id', 'treatment', 'fluctuation')) %>% #mean_V is the mean of two replicate MC
          group_by(id, treatment, fluctuation) %>%
          summarise(BV_010 = mean_V/mean_V0,
                 se_010 = se_V/se_V0)
#Schritt 4: divide values by t0 values of the treatment 

###subset of control data to later merge them as a new column with the other df
data_co <- BV_each %>%
  filter(treatment == 'con') %>%
  rename(con_div = BV_010,
         se_div = se_010) %>%
  ungroup()%>%
  select(-treatment, -fluctuation)

## Schritt 4. Merge all data together, we have now two columns of values (both divdided by t0):
# Treatment and control (divdided by t0 respectively), those will be divided now for ARF
AR <- left_join(BV_each, data_co, by = c('id')) %>%
  filter(treatment != 'con') %>%
  mutate(AR = BV_010/con_div,
         se_AR = se_010/se_div)

#Visulasitation
ggplot(AR, aes(x = fluctuation, y = AR, col = id))+
  geom_point(size = 1.8, position=position_dodge(width=2))+
  geom_errorbar(aes(ymin = AR - se_AR, ymax = AR + se_AR), width = .8, position=position_dodge(width=2))+
  scale_x_continuous(limits = c(4, 52), breaks = c(6,12,24,36,48))+
  geom_hline(yintercept = 0)+  
  facet_wrap(~id, scales = 'free_y')+
  labs(y = 'ARF divided by treatment 0 values')+
  theme_bw()
#ggsave(plot = last_plot(), file = 'ARF_divided0_treatment.png' )
