# scatter plot with bray-curits diss and LRR of function (carbon)
# By Charlotte Kunze 06.08.2020

#### Load packages & import data ####
library(vegan)

## import data ###
counts <- read_delim("~/Desktop/MA/MA_Rcode/project_data/phyto_new.csv", 
                     ";", escape_double = FALSE, col_types = cols(date = col_character()), 
                     locale = locale(decimal_mark = ","), 
                     trim_ws = TRUE) %>%
  drop_na(MC)
str(counts) #check import
names(counts)
#metadata
df <- data.frame( MC = c('1', '4', '2', '10', '3', '5', '6', '7', '8', '9', '11', '12'),
                  treatment_id = c('control', 'control', 'Fluctuating_36','Fluctuating_36',
                                   'Fluctuating_6', 'Fluctuating_6', 'Fluctuating_48', 'Fluctuating_48',
                                   'Fluctuating_24', 'Fluctuating_24', 'Fluctuating_12','Fluctuating_12') )

#change variable format to numeric 
counts$cells_ml <- as.numeric(gsub(",", ".", counts$cells_ml))
counts$MC = as.character(counts$MC)
counts$cells_ml[is.na(counts$cells_ml)] <-0


#### PCA on species contribution to BioVolume####

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
  ungroup() %>%
  mutate(species = paste(genus, species, sep = ' ')) %>%
  dplyr::select(-genus)


#### merge all_data ####
all_data <- left_join(counts, df, by = c('MC')) %>%
  left_join(., BioV, by = c('species') ) %>%
  dplyr::select(treatment_id, MC, sampling, species, phylum, cells_ml,mean_BioV) %>%
  mutate(volume = cells_ml*mean_BioV)  #calculate biovolume


multi_data<- all_data%>%
  group_by(sampling, treatment_id, MC) %>%
  mutate(sum = sum(volume, na.rm  =T ),
         rel = volume/sum) %>%
  drop_na(rel)
bray <-multivariate_change(multi_data, time.var='sampling', species.var='species', abundance.var= 'rel',
                    replicate.var = 'MC', treatment.var = 'treatment_id', reference.time = NULL)

#new df to calculate mean biovolume
pca_BV <- all_data%>%
  group_by(sampling, treatment_id, MC) %>%
  mutate(sum = sum(volume, na.rm  =T ),
         rel = volume/sum)%>%
  group_by(treatment_id, sampling, species) %>%
  summarise(mean_V = mean(rel, na.rm = T),
            sd_V = sd(rel, na.rm = T),
            se_V = sd_V/sqrt(n())) %>% #calculate mean cells with biovolume
  dplyr::select(-sd_V) %>%
  drop_na(mean_V) %>%
  ungroup()%>%
  mutate(dummy = paste(treatment_id, sampling, sep = ' ')) %>%
  dplyr::select(-treatment_id, -sampling, -se_V)%>%
  spread(key = species, value = mean_V) %>%
  column_to_rownames('dummy')

pca_BV[is.na(pca_BV)] <- 0

data.dist <- vegdist(pca_BV, "bray") %>%
  broom::tidy() %>%
  separate(item1, into = c('treatment', 'sampling'), ' ') %>%
  separate(item2, into = c('treatment2', 'sampling2'), ' ') %>%
  group_by(treatment, sampling) %>%
  summarise(average = mean(distance, na.rm = T),
            sd = sd(distance, na.rm = T),
            se = sd/sqrt(n())) %>%
  filter(treatment!= 'control') %>%
  separate(treatment, into = c('Fluctuating', 'Fluctuation'), '_') %>%
  dplyr::select(-Fluctuating) %>%
  mutate(sampling = as.numeric(sampling))


#### import LRRs for function data ####
stability <- read_delim("~/Desktop/MA/MA_Rcode/project_data/stability_MA.csv", 
                           ";", escape_double = FALSE, locale = locale(decimal_mark = ","), 
                           trim_ws = TRUE) %>%
  mutate(Fluctuation = str_remove(caseID, 'F')) %>%
  mutate(sampling = DAY/2) %>%
  dplyr::select(Fluctuation, sampling, LRR, var.lrr)


dissimi_data <- left_join(stability, data.dist, by = c('Fluctuation', 'sampling')) %>%
  drop_na(average)
dissimi_data$Fluctuation <- factor(as.factor(dissimi_data$Fluctuation),levels=c( "48", "36", '24', '12', '6'))

ggplot(dissimi_data, aes(x = average, y = LRR, fill = Fluctuation))+
  geom_point(aes(alpha = as.factor(sampling)), size = 3, pch = 21, col = 'black')+
  geom_hline(yintercept = 0)+
  scale_y_continuous(limits = c(-0.5, 0.5), breaks = c(-0.4,-0.2,0,0.2,0.4))+
  scale_fill_manual(values = c( '#0868ac','#41b6c4','#31a354','#addd8e','#fed976'))+
  labs( x = 'compositional turnover (bray-curtis)', y = 'LRR of POC', title = 'correlation funct and comp changes', alpha = 'Sampling')+
  theme( panel.background = element_rect(fill = NA), #loescht den Hintergrund meines Plots/ fuellt ihn mit nichts
         #panel.grid.major.y = element_line(color='grey', linetype = 'dashed', size=0.2),
         panel.border= element_rect(colour = "black", fill=NA, size=0.5),
         strip.background = element_rect(color ='black', fill = 'white'),
         strip.text = element_text(face = 'italic'),
         legend.background = element_blank(),
         legend.key = element_blank(),
         text = element_text(size=14))
#ggsave(plot = last_plot(), file = 'dissimilarity plot.png')

