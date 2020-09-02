#Clear all
rm(list=ls(all=TRUE))
dev.off()
#https://jkzorz.github.io/2019/06/06/NMDS.html
##### Script for ndms plots 
# by Laura Hennigs & Charlotte Kunze
##-----------------------------------------------------------------------------------------------##
#check working directory
getwd()
##-----------------------------------------------------------------------------------------------##
# Load required packages
library(tidyverse)
library(lubridate)
library(scales)
library(vegan)

####import data####
#using the manual import function (readr) and after that copying the console output
phyto_new <- read_delim("~/Desktop/MA/MA_Rcode/project_data/phyto_new.csv", 
                        ";", escape_double = FALSE, locale = locale(decimal_mark = ","), 
                        trim_ws = TRUE)
phyto_new$cells_ml <- as.numeric(gsub(",", ".", phyto_new$cells_ml))

####a data frame#### 
#with treatment information, which has the same column planktotron_no as in data.scores
df=data.frame(treat=c('0','36','6','0','6','48','48',
                      '24','24','36','12','12'),  
              MC = c(1:12) )

#combine the two df

# new data frame including data wrangling
data_nd <- phyto_new %>%
  #rename(MC = planktotron) %>%
  select(MC, species, sampling, cells_ml) %>% 
  na.omit() %>% #select columns we want to keep
  filter(MC != 'Blank') %>%
  full_join(df, df, by = c('MC')) %>%
  group_by(sampling, treat, species) %>%
  summarise(mean_cells = mean(cells_ml, na.rm = T)) %>%
  spread(key = 'species', value = 'mean_cells',  fill = NA, convert=F)       #convert the column 'species' to many columns, for each species one. Fill these columns with the cells per ml for that species

#fill all N.A.s with 0
data_nd[3:ncol(data_nd)][is.na(data_nd[3:ncol(data_nd)])] <-0          
str(data_nd)


com = data_nd[, 3:ncol(data_nd)]   #create a data frame from column 4 till 18 (all the species)
m_com = as.matrix(com)             #the data frame we just created should be seen as matrix

#nmds
nmds = metaMDS(m_com, k=2)         #k are the dimensions for the plots
plot(nmds)                         #not a nice plot, so we continue in order to plot it with ggplot
#within the plot red crosses = species, open circles = communities

#extract nmds scores
data.scores = as.data.frame(scores(nmds)) #obtain the coordinates for the axes and put them in a new data frame

#add columns with further information
data.scores$sampling = data_nd$sampling
#data.scores$date = data_nd$date
data.scores$treat = data_nd$treat


head(data.scores)  #see our results (but only the first rows)

#plot with ggplot
ggplot(data.scores, aes(x = NMDS1, y = NMDS2, colour = treat, shape = as.factor(sampling)))+           #treatments in same colour
  geom_point(size = 4, aes(color = treat))+
  scale_color_manual(values = c( '#000000','#0868ac','#41b6c4','#31a354','#addd8e','#fed976'))+
  #scale_shape_manual(values = c(15,10,17))+
  labs(shape = 'sampling', color = 'Fluctuatin frequency')+  
  annotate('text', x = 0.35, y = 0.38, label = paste('Stress:', round(nmds$stress, 3)))+   #3 gives number of digits
  theme_bw()        
#ggsave(plot = last_plot(), file = 'NMDS_counting.png')  
## stress
# stress < 0.05 excellent representation
# stress < 0.1 great
# stress < 0.2 okay
# stress < 0.3 poor represented
# stress > 0.3 choose another method


################################################################################
##### Ordiellipse Graph####
NMDS = data.frame(NMDS1 = data$NMDS1, NMDS2=data$NMDS2, group=as.factor(data$treat))#sets up data frame 

g1<-ggplot(data = NMDS, aes(NMDS1, NMDS2))+   
  geom_point(aes(color = group), size=3)+  
  stat_ellipse(geom = 'polygon', type = "norm", linetype = 2, aes(fill=group, col = group), alpha = 0.2, level = 0.99) +
  theme_bw() 
g1

####################################################################################

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


#### NMDS on species contribution to BioVolume####

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

#new df to calculate mean biovolume
nmds_BV <- all_data%>%
  group_by(sampling, treatment_id, MC) %>%
  mutate(sum = sum(volume, na.rm  =T ),
         rel = volume/sum)%>%
  group_by(treatment_id, sampling, species) %>%
  summarise(mean_V = mean(rel, na.rm = T)) %>% #calculate mean cells with biovolume
  drop_na(mean_V) %>%
  separate(treatment_id, into = c('treat', 'fluctuation'), '_') %>%
  dplyr::select(-treat) %>%
  spread(key = species, value = mean_V) 
nmds_BV$fluctuation[is.na(nmds_BV$fluctuation)] <-0
#fill all N.A.s with 0
nmds_BV[3:ncol(nmds_BV)][is.na(nmds_BV[3:ncol(nmds_BV)])] <-0          
str(nmds_BV)


com1 = nmds_BV[, 3:ncol(nmds_BV)]   #create a data frame from column 4 till 18 (all the species)
m_com1 = as.matrix(com1)             #the data frame we just created should be seen as matrix

#nmds
nmds1 = metaMDS(m_com1, k=2)         #k are the dimensions for the plots
plot(nmds1)                         #not a nice plot, so we continue in order to plot it with ggplot
#within the plot red crosses = species, open circles = communities

#extract nmds scores
data.scores1 = as.data.frame(scores(nmds1)) #obtain the coordinates for the axes and put them in a new data frame

#add columns with further information
data.scores1$sampling = nmds_BV$sampling
data.scores1$day = nmds_BV$sampling*2
#data.scores$date = data_nd$date
data.scores1$treatment = nmds_BV$fluctuation


head(data.scores1)  #see our results (but only the first rows)

#plot with ggplot
data.scores1$treatment <- factor(data.scores1$treatment, levels= c('0', '48', '36', '24', '12', '6'))
ggplot(data.scores1, aes(x = NMDS1, y = NMDS2, colour = treatment, shape = as.factor(day)))+           #treatments in same colour
  geom_point(size = 4, aes(color = treatment))+
  scale_color_manual(values = c( '#000000','#0868ac','#41b6c4','#31a354','#addd8e','#fed976'))+
  #scale_shape_manual(values = c(15,10,17))+
  labs(shape = 'Time [days]', color = 'Fluctuation \nfrequency [h]')+  
  annotate('text', x = 0.2, y = 0.38, label = paste('Stress:', round(nmds$stress, 3)))+   #3 gives number of digits
  theme( panel.background = element_rect(fill = NA), #loescht den Hintergrund meines Plots/ fuellt ihn mit nichts
         #panel.grid.major.y = element_line(color='grey', linetype = 'dashed', size=0.2),
         panel.border= element_rect(colour = "black", fill=NA, size=1),
         strip.background = element_rect(color ='black', fill = 'white'),
         strip.text = element_text(face = 'bold'),
         legend.background = element_blank(),
         legend.position  ='right',
         legend.key = element_blank(),
         text = element_text(size=18))        
ggsave(plot = last_plot(), file = 'NMDS_BV.png', width = 8, height = 6)  
## stress
# stress < 0.05 excellent representation
# stress < 0.1 great
# stress < 0.2 okay
# stress < 0.3 poor represented
# stress > 0.3 choose another method

