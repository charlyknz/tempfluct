## PCA on species data 
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
# 1. Phytoplankton counts
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


#### PCA on counts only ####
count_dat <- counts %>%
  ungroup()%>%
  dplyr::select(MC, sampling, species, phylum, cells_ml) %>%
  # filter(sampling != 0 & species != 'Rhizosolemia indet.')%>%
  #bind_rows(., spec0) %>%
  left_join(., df, by = c('MC')) %>%
  mutate(dummy = paste(treatment_id, sampling, sep = ' ')) %>%
  group_by(treatment_id, sampling, species, dummy) %>%
  summarise(cell_per_ml = mean(cells_ml, na.rm = T)) %>%
  filter(species!= 'Rhizosolemia indet.')%>%
  spread(key = species, value = cell_per_ml) %>%
  ungroup()%>%
  column_to_rownames("dummy") %>% #explanatory variables as row names
  dplyr::select(-sampling,-treatment_id)
count_dat[is.na(count_dat)] <-0 #make sure we don't have NAs
str(count_dat)



#### perform PCA #### 
prin_comp <- prcomp(count_dat, scale. = T)
names(prin_comp)

# see how many dimensions we have
dim(prin_comp$x)
biplot(prin_comp, scale = 0)
#[1] 30 30

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
prop_varex[1:5]
plot(prop_varex, xlab = "Principal Component",
     ylab = "Proportion of Variance Explained",
     type = "b")

#cumulative scree plot
plot(cumsum(prop_varex), xlab = "Principal Component",
     ylab = "Cumulative Proportion of Variance Explained",
     type = "b")


#### plot PC1 & 2 ####
#PCA dataframe including our PC1:30 axes

PCA <- prin_comp$x %>% 
  as.data.frame %>%
  rownames_to_column("dummy")%>%
  separate(dummy, into = c("fluctuation", 'sampling'), ' ', convert = TRUE) %>%
  dplyr::select(PC1, PC2, PC3, PC4, fluctuation, sampling)

##change format to bring variables in the right direction
PCA$fluctuation <- factor(as.factor(PCA$fluctuation),levels=c("control", "Fluctuating_48", "Fluctuating_36", 'Fluctuating_24', 'Fluctuating_12', 'Fluctuating_6'))
PCA$sampling <- gsub("6", "06", PCA$sampling)
PCA$sampling =as.numeric(PCA$sampling)
PCA <- arrange(PCA, sampling)

ggplot(PCA,aes(x=PC1,y=PC2)) + 
  geom_point(size=1.5,aes(color = fluctuation,  group = fluctuation)) +
  #geom_line(arrow = arrow( ends = "both", type = "open"), aes(linetype = fluctuation))+
  geom_path(aes(x = PC1, y = PC2, group = fluctuation, col = fluctuation, linetype = fluctuation), 
            arrow = arrow(ends = 'last',type = 'closed',length = unit(0.35, "cm")))+
  scale_color_manual(values = c( '#000000', '#0868ac','#31a354','#41b6c4','#addd8e','#fed976'))+
  labs(x=paste0("PC1: ",round(prop_varex[1]*100,1),"%"),
       y=paste0("PC2: ",round(prop_varex[2]*100,1),"%")) +
  facet_wrap(~fluctuation)+
  theme( panel.background = element_rect(fill = NA), #loescht den Hintergrund meines Plots/ fuellt ihn mit nichts
         #panel.grid.major.y = element_line(color='grey', linetype = 'dashed', size=0.2),
         panel.border= element_rect(colour = "black", fill=NA, size=0.5),
         strip.background = element_rect(color ='black', fill = 'white'),
         strip.text = element_text(face = 'bold'),
         legend.background = element_blank(),
         legend.position  ='bottom',
         legend.key = element_blank(),
         text = element_text(size=14))
#ggsave(plot = last_plot(), file = 'PCA_species_abs.png')


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
prop_varex[1:5]
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
  separate(fluctuation_sampling,into = c("fluctuation", 'sampling'),' ') %>%
  dplyr::select(PC1, PC2, PC3, PC4, fluctuation, sampling) %>%
  mutate(sampling = as.numeric(sampling))%>%
  arrange(sampling)

PCA_BV$fluctuation <- factor(as.factor(PCA_BV$fluctuation),levels=c("control", "Fluctuating_48", "Fluctuating_36", 'Fluctuating_24', 'Fluctuating_12', 'Fluctuating_6'))
ggplot(PCA_BV,aes(x=PC1,y=PC2)) + 
  geom_hline(yintercept = 0, linetype = 'dashed', size = 0.5)+
  geom_vline(xintercept = 0, linetype = 'dashed', size = 0.5)+
  geom_point(size=1.5,aes( group = fluctuation , col = fluctuation)) +
  #geom_line(arrow = arrow( ends = "both", type = "open"), aes(linetype = fluctuation))+
  geom_path(aes(x = PC1, y = PC2, group = fluctuation, col = fluctuation, linetype = fluctuation), 
            arrow = arrow(ends = 'last',type = 'closed',length = unit(0.35, "cm")))+
  scale_color_manual(values = c( '#000000', '#0868ac','#41b6c4','#31a354','#addd8e','#fed976'))+
  scale_y_continuous(limits = c(-6.2, 6.2), breaks = c(-6,-4,-2,0,2,4,6))+
  scale_x_continuous(limits = c(-6.2, 6.2), breaks = c(-6,-4,-2,0,2,4,6))+
  labs(x=paste0("PC1: ",round(prop_varex[1]*100,1),"%"),
       y=paste0("PC2: ",round(prop_varex[2]*100,1),"%")) +
  facet_wrap(~fluctuation)+
   theme( panel.background = element_rect(fill = NA), #loescht den Hintergrund meines Plots/ fuellt ihn mit nichts
         #panel.grid.major.y = element_line(color='grey', linetype = 'dashed', size=0.2),
         panel.border= element_rect(colour = "black", fill=NA, size=0.5),
         strip.background = element_rect(color ='black', fill = 'white'),
         strip.text = element_text(face = 'bold'),
         legend.background = element_blank(),
         legend.position  ='bottom',
         legend.key = element_blank(),
         text = element_text(size=14))
#ggsave(plot = last_plot(), file = 'PCA_rel_specBV.png', width = 8, height = 7)


#### Euclidian distance ####
dist <- PCA %>%
  mutate(id = paste(fluctuation, sampling, sep =' ')) %>%
  column_to_rownames('id')%>%
  dplyr::select(-fluctuation, -sampling)%>%
  dist() %>%
  broom::tidy() %>% 
  separate(item1, c('fluctuation', 'sampling'), ' ')%>%
  separate(item2, c('fluctuation2', 'sampling2'), ' ') %>%
  group_by(sampling)%>%
  summarise(mean = mean(distance, na.rm =T),
            sd = sd(distance, na.rm = T),
            se = sd/sqrt(n()))
dist$sampling =as.numeric(dist$sampling)
dist <- arrange(dist, sampling)

ggplot(dist, aes(x = sampling, y = mean))+
  geom_point()+
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width = .5)+
  geom_line()+
  scale_y_continuous(limits = c(0, 6), breaks = seq(0,6,2))+
  labs(x = 'sampling', y = 'distance among treatments', title = 'species contribution to BV')+
  theme_bw()
#ggsave(plot = last_plot(), file = 'distance_spec.png')
