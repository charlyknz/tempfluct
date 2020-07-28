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
rich <- left_join(data1, treatments, by = c('MC'))  %>%
  select(-comment)


shannon <- rich %>%
  group_by(treatment, sampling, species) %>%
  summarise(mean_cells = mean(cells_ml, na.rm = T)) %>%
  ungroup()%>%
  spread(key = species, value = mean_cells) %>%
  group_by(treatment, sampling) %>%
  mutate(Gymnodinium_sp = sum(Gymnodinium, Gymnodinium_Lepodinium)) %>%
  select(-Gymnodinium_Lepodinium, -Gymnodinium)
 
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


  