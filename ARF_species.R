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

#new df to calculate mean biovolume
calc_BV <- all_data%>%
  filter(sampling == 10) %>% #remove sampling 0 
  group_by(treatment, fluctuation, sampling, id, genus, species) %>%
  summarise(mean_V = mean(volume, na.rm = T),
            sd_V = sd(volume, na.rm = T),
            se_V = sd_V/sqrt(n())) %>% #calculate mean cells with biovolume
  select(-sd_V)
#calculate mean Volume for sampling 0, assuming all treatments are the same at the beginning
data0 <- all_data %>% 
  filter(sampling == 0) %>%
  group_by( sampling, id, genus, species) %>%
  summarise(con = mean(volume, na.rm = T),
            F6 = mean(volume, na.rm = T),
            F12 = mean(volume, na.rm = T),
            F24 = mean(volume, na.rm = T),
            F36 = mean(volume, na.rm = T),
            F48 = mean(volume, na.rm = T),
            sd_con = sd(volume, na.rm = T),
            sd_F6 = sd(volume, na.rm = T),
            sd_F12 = sd(volume, na.rm = T),
            sd_F24 = sd(volume, na.rm = T),
            sd_F36 = sd(volume, na.rm = T),
            sd_F48 = sd(volume, na.rm = T),
            se_con= sd_con/sqrt(n()),
            se_F6= sd_F6/sqrt(n()),
            se_F12= sd_F12/sqrt(n()),
            se_F24= sd_F24/sqrt(n()),
            se_F36= sd_F36/sqrt(n()),
            se_F48= sd_F48/sqrt(n()))%>%
  select(-sd_con,-sd_F6, -sd_F12, -sd_F24, -sd_F36, -sd_F48)%>%
  gather(key = 'treatment', value = 'mean_V',  -sampling, -id,-genus,-species, -se_con,-se_F6, -se_F12, -se_F24, -se_F36, -se_F48) %>%
  gather(key = 'variable', value = 'se_V', -sampling, -id, -genus,-species,-treatment,-mean_V) %>%
  select(-variable) %>%
  distinct(sampling, id, genus, species, treatment, mean_V, se_V) %>%
  na.omit()

##merge zero-values and data-rest
all_BV <- dplyr::bind_rows(calc_BV, data0) %>%
  drop_na(mean_V) #mean_V is the mean of two replicate MC



#### calculate Algal Response Factor ###

#Step 1: subset of sampling 0 values, renaming to merge with other data 
# aiming to divide t10 / t0 values 

ARF0 <- all_BV %>%
  filter(sampling == 0) %>%
  ungroup()%>%
  rename(t0_volume = mean_V,
         t0_se = se_V) %>%
  select(-sampling,-fluctuation ) 

# subset of sampling 10, remove columns and calculate the mean biovolume per phylum
# Schritt 2. merge then with 0 data to calculate value /divided by 0 for ARF 
data_ARF <- all_BV %>%
  filter(sampling ==10) %>%
  ungroup()%>%
  select(-sampling) %>%
  left_join(., ARF0, by = c('treatment', 'id', 'genus', 'species')) %>%
  group_by(id, genus, species, treatment, fluctuation) %>%
  summarise( mean_t10_0 = mean_V/t0_volume,
             se_t10_0 = se_V/t0_se)  %>%
  mutate(species_id = paste(genus, species, sep = ' '))

#visualise mean divided by t0 values:
ggplot(data_ARF, aes(x = fluctuation, y = mean_t10_0, col = treatment))+
  geom_point()+
  geom_errorbar(aes(ymin = mean_t10_0 - se_t10_0, ymax = mean_t10_0 + se_t10_0), width = .8)+
  scale_x_continuous(limits = c(-1, 52), breaks = c(0,6,12,24,36,48))+
  facet_wrap(~species_id, scales = 'free_y')+
  theme_bw()

# Schritt 3: subset of control data to later merge them as a new column with the other df
data_con <- data_ARF %>%
  filter(treatment == 'con') %>%
  rename(con_div0 = mean_t10_0,
         se_div0 = se_t10_0) %>%
  ungroup()%>%
  select(-treatment, -fluctuation)

## Schritt 4. Merge all data together, we have now two columns of values (both divdided by t0):
# Treatment and control (divdided by t0 respectively), those will be divided now for ARF
ARF <- left_join(data_ARF, data_con, by = c('id', 'species_id','genus', 'species')) %>%
  filter(treatment != 'con') %>%
  mutate(ARF = mean_t10_0/con_div0,
         se_ARF = se_t10_0/se_div0) %>%
  filter(!genus %in% c('Alexandrium', 'Pinnularia', 'Peridiniella', 'Odontella') )

#Visulasitation
library(viridis)
ggplot(subset(ARF, id == 'Ciliophora'), aes(x = fluctuation, y = ARF, col = fluctuation))+
  geom_point(size = 1.8)+
  geom_errorbar(aes(ymin = ARF - se_ARF, ymax = ARF + se_ARF), width = .8)+
  scale_x_continuous(limits = c(4, 52), breaks = c(6,12,24,36,48))+
  geom_hline(yintercept = 0)+ 
  facet_wrap(~species_id, scales = 'free_y')+
  scale_color_viridis(option = "D", direction = -1)+
theme_bw()
#ggsave(plot = last_plot(), file = 'ARF_plot_Ciliophora.png' )

####  together ####
ggplot(subset(ARF, id != 'Ciliophora'), aes(x = fluctuation, y = ARF))+
  geom_point(aes(fill = as.factor(fluctuation)), pch=21, size=4, col = 'black',position = position_dodge(width = 3))+
  geom_errorbar(aes(ymin = ARF - se_ARF, ymax = ARF + se_ARF), width = .8,position = position_dodge(width = 3))+
  scale_x_continuous(limits = c(4, 52), breaks = c(6,12,24,36,48))+
  geom_hline(yintercept = 0)+ 
  #scale_shape_manual(values = c(1, 15,6,17))+
  facet_wrap(~species_id, scales = 'free_y', ncol = 5)+
  #scale_color_viridis(option = "D", direction = -1, discrete = F)+
  scale_fill_manual(values = c( '#fed976','#addd8e','#31a354','#41b6c4','#0868ac'))+
  labs( x = 'Fluctuation frequency (in h)', col = 'treatment', title = 'species specific ARF at sampling 10', fill = 'Fluctuation frequency')+
  theme( panel.background = element_rect(fill = NA), #loescht den Hintergrund meines Plots/ fuellt ihn mit nichts
         #panel.grid.major.y = element_line(color='grey', linetype = 'dashed', size=0.2),
         panel.border= element_rect(colour = "black", fill=NA, size=0.5),
         strip.background = element_rect(color ='black', fill = 'white'),
         strip.text = element_text(face = 'italic'),
         legend.background = element_blank(),
         legend.position  ='bottom',
         legend.key = element_blank(),
         text = element_text(size=14))
#ggsave(plot = last_plot(), file = 'ARF_allspec.png', width = 11, height = 7 )
