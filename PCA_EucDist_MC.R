#### Skript to perform a Principal Component Analysis (PCA)####
# for a detailed tutorial see: http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/118-principal-component-analysis-in-r-prcomp-vs-princomp/

# By Charlotte Kunze 06.08.2020

#### Load packages & import data ####
library(tidyverse)


#### helpful links: 
#https://uc-r.github.io/pca
#https://cmdlinetips.com/2019/04/introduction-to-pca-with-r-using-prcomp/
#https://www.analyticsvidhya.com/blog/2016/03/pca-practical-guide-principal-component-analysis-python/

###################################################################################################################  
#### PCA 1. Phytoplankton counts ####
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

data0 <- filter(all_data,sampling == 0)
data0$MC[data0$MC == 4] <-1
data0$MC[data0$MC == 11] <-12
data0$MC[data0$MC == 8] <-9
data0$MC[data0$MC == 3] <-5
data0$MC[data0$MC == 6] <-7
data0$MC[data0$MC == 10] <-2

#new df to calculate mean biovolume
pca_BV <- all_data%>%
  bind_rows(., data0)%>%
  group_by(sampling, treatment_id, MC) %>%
  mutate(sum = sum(volume, na.rm  =T ),
         rel = volume/sum)%>%
  group_by(treatment_id, MC, sampling, species) %>%
  summarise(mean_V = mean(rel, na.rm = T),
            sd_V = sd(rel, na.rm = T),
            se_V = sd_V/sqrt(n())) %>% #calculate mean cells with biovolume
  dplyr::select(-sd_V) %>%
  drop_na(mean_V) %>%
  ungroup()%>%
  mutate(dummy = paste(treatment_id, MC, sampling, sep = ' ')) %>%
  dplyr::select(-treatment_id, -sampling, -se_V, -MC)%>%
  spread(key = species, value = mean_V) %>%
  column_to_rownames('dummy')

pca_BV[is.na(pca_BV)] <- 0

#### perform PCA #### 
prin_comp <- prcomp(pca_BV, scale. = T)
names(prin_comp)

# see how many dimensions we have
dim(prin_comp$x)
biplot(prin_comp, scale = 0)

#2. The prcomp() function also provides the facility to compute standard deviation 
#of each principal component. sdev refers to the standard deviation of principal components.

#compute standard deviation of each principal component
std_dev <- prin_comp$sdev

#compute variance
pr_var <- std_dev^2

#check variance of first 10 components
pr_var[1:10]

#3. To compute the proportion of variance explained by each component, we simply divide the variance by sum of total variance. 
#This results in:

#proportion of variance explained
prop_varex <- pr_var/sum(pr_var)
prop_varex[1:6]
plot(prop_varex, xlab = "Principal Component",
     ylab = "Proportion of Variance Explained",
     type = "b")
#cumulative scree plot
plot(cumsum(prop_varex), xlab = "Principal Component",
     ylab = "Cumulative Proportion of Variance Explained",
     type = "b")

PCA_BV <- prin_comp$x %>% 
  as.data.frame %>%
  rownames_to_column("fluctuation_sampling") %>%
  separate(fluctuation_sampling,into = c("fluctuation", 'MC','sampling'),' ') %>%
  dplyr::select(PC1, PC2, PC3, PC4, PC5, PC6,fluctuation,MC, sampling) %>%
  mutate(sampling = as.numeric(sampling))%>%
  arrange(sampling)

par(mfrow=c(1,1),cex.axis=1.2, cex.lab=1.5)

PCA_BV$fluctuation <- factor(as.factor(PCA_BV$fluctuation),levels=c("control", "Fluctuating_48", "Fluctuating_36", 'Fluctuating_24', 'Fluctuating_12', 'Fluctuating_6'))
ggplot(PCA_BV,aes(x=PC1,y=PC2)) + 
  geom_hline(yintercept = 0, linetype = 'dashed', size = 0.5)+
  geom_vline(xintercept = 0, linetype = 'dashed', size = 0.5)+
  geom_point(size=1.5,aes( group = fluctuation , col = fluctuation)) +
  #geom_line(arrow = arrow( ends = "both", type = "open"), aes(linetype = fluctuation))+
  geom_path(aes(x = PC1, y = PC2, group = fluctuation, col = fluctuation, linetype = fluctuation), 
            arrow = arrow(ends = 'last',type = 'closed',length = unit(0.35, "cm")))+
  scale_color_manual(values = c( '#000000', '#0868ac','#41b6c4','#31a354','#addd8e','#fed976'))+
  #scale_y_continuous(limits = c(-6.2, 6.2), breaks = c(-6,-4,-2,0,2,4,6))+
  #scale_x_continuous(limits = c(-6.2, 6.2), breaks = c(-6,-4,-2,0,2,4,6))+
  labs(x=paste0("PC1: ",round(prop_varex[1]*100,1),"%"),
       y=paste0("PC2: ",round(prop_varex[2]*100,1),"%")) +
  facet_wrap(~MC)+
  theme( panel.background = element_rect(fill = NA), #loescht den Hintergrund meines Plots/ fuellt ihn mit nichts
         #panel.grid.major.y = element_line(color='grey', linetype = 'dashed', size=0.2),
         panel.border= element_rect(colour = "black", fill=NA, size=0.5),
         strip.background = element_rect(color ='black', fill = 'white'),
         strip.text = element_text(face = 'bold'),
         legend.background = element_blank(),
         legend.position  ='bottom',
         legend.key = element_blank(),
         text = element_text(size=14))
#ggsave(plot = last_plot(), file = 'PCA_rel_specBV_MC.png', width = 8, height = 7)


#### Euclidian distance ####
dist <- PCA_BV %>%
  mutate(id = paste(fluctuation, MC, sampling, sep =' ')) %>%
  column_to_rownames('id')%>%
  dplyr::select(-fluctuation, -sampling)%>%
  dist() %>%
  broom::tidy() %>% 
  separate(item1, c('fluctuation', 'MC','sampling'), ' ')%>%
  separate(item2, c('fluctuation2', 'MC2','sampling2'), ' ') %>%
  group_by(sampling, fluctuation)%>%
  filter(fluctuation2 == 'control' & fluctuation != 'control')%>%
  filter(sampling == sampling2) %>%
  summarise(mean = mean(distance, na.rm =T),
            sd = sd(distance, na.rm = T),
            se = sd/sqrt(n()))
dist$sampling =as.numeric(dist$sampling)
dist <- arrange(dist, sampling)
dist$fluctuation <- factor(as.factor(dist$fluctuation),levels=c("control", "Fluctuating_48", "Fluctuating_36", 'Fluctuating_24', 'Fluctuating_12', 'Fluctuating_6'))

ggplot(subset(dist, fluctuation != 'control'), aes(x = sampling, y = mean, group = fluctuation))+
  geom_point(aes(fill = fluctuation), pch = 21, size = 3)+
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width = .8)+
  geom_line(linetype = 'dashed', aes(col = fluctuation))+
  scale_y_continuous(limits = c(0, 15), breaks = seq(0,15,5))+
  scale_fill_manual(values = c('#0868ac','#41b6c4','#31a354','#addd8e','#fed976'))+
  scale_color_manual(values = c('#0868ac','#41b6c4','#31a354','#addd8e','#fed976'))+
  labs(x = 'sampling', y = 'Euclidian distance to control (based on PC 1:6)', title = 'dissimilarity of species contribution to BV')+
  theme( panel.background = element_rect(fill = NA), #loescht den Hintergrund meines Plots/ fuellt ihn mit nichts
         #panel.grid.major.y = element_line(color='grey', linetype = 'dashed', size=0.2),
         panel.border= element_rect(colour = "black", fill=NA, size=0.5),
         strip.background = element_rect(color ='black', fill = 'white'),
         strip.text = element_text(face = 'italic'),
         legend.background = element_blank(),
         legend.key = element_blank(),
         text = element_text(size=12))
#ggsave(plot = last_plot(), file = 'distance_treatment_to_control.png', width = 10, height = 6)

