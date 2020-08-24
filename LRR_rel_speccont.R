#### alternative ARF for relative species contribution ####
# calculate LRR

library(tidyverse)
library(vegan)


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
  group_by(sampling, treatment, MC) %>%
  mutate(sum = sum(volume, na.rm  =T ),
         rel = volume/sum)%>%
  ungroup()%>%
  group_by(treatment, fluctuation, sampling, id, genus, species) %>%
  summarise(mean_V = mean(rel, na.rm = T),
            sd_V = sd(rel, na.rm = T),
            se_V = sd_V/sqrt(n())) %>% #calculate mean cells with biovolume
  select(-sd_V)

data_ARF <- calc_BV %>%
  ungroup()%>%
  mutate(species_id = paste(genus, species, sep = ' '))

# Schritt 3: subset of control data to later merge them as a new column with the other df
data_con <- data_ARF %>%
  filter(treatment == 'con') %>%
  rename(con_div = mean_V,
         se_div = se_V) %>%
  ungroup()%>%
  select(-treatment, -fluctuation)

LRR <- left_join(data_ARF, data_con, by = c('id', 'species_id','genus', 'species', 'sampling')) %>%
  filter(treatment != 'con') %>%
  mutate(LRR = log10(mean_V/con_div),
         se_LRR = log10(se_V/se_div)) %>%
  filter(!genus %in% c('Alexandrium', 'Pinnularia', 'Peridiniella', 'Odontella') ) %>%
  filter(!LRR %in% c(Inf, -Inf) )

ggplot(subset(LRR, id != 'Ciliophora'  & sampling == 10), aes(x = fluctuation, y = LRR, group = fluctuation))+
  geom_point(aes(fill = as.factor(fluctuation)), pch=21, size=3, col = 'black', position = position_dodge(width = 2))+
  geom_errorbar(aes(ymin = LRR - se_LRR, ymax = LRR + se_LRR), width = .8,position = position_dodge(width = 2))+
  #geom_line(linetype = 'dashed')+
  # scale_x_continuous(limits = c(4, 52), breaks = c(6,12,24,36,48))+
  geom_hline(yintercept = 0)+ 
  #scale_shape_manual(values = c(1, 15,6,17))+
  facet_wrap(~species_id,ncol = 5)+
  #scale_color_viridis(option = "D", direction = -1, discrete = F)+
  scale_fill_manual(values = c( '#fed976','#addd8e','#31a354','#41b6c4','#0868ac'))+
  labs( x = 'Fluctuation frequency (in h)', col = 'treatment', title = 'species specific LRR at sampling 10', fill = 'Fluctuation frequency')+
  theme( panel.background = element_rect(fill = NA), #loescht den Hintergrund meines Plots/ fuellt ihn mit nichts
         #panel.grid.major.y = element_line(color='grey', linetype = 'dashed', size=0.2),
         panel.border= element_rect(colour = "black", fill=NA, size=0.5),
         strip.background = element_rect(color ='black', fill = 'white'),
         strip.text = element_text(face = 'italic'),
         legend.background = element_blank(),
         legend.position  ='bottom',
         legend.key = element_blank(),
         text = element_text(size=14))#
ggsave(plot = last_plot(), file = 'rel_LRR_allspec_s10.png', width = 11, height = 7 )
