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


#the same initial cell abundance for every treatment/MC
spec0 <- counts %>% 
  filter(sampling == 0) %>%
  group_by( sampling, species, phylum) %>%
  summarise('1' = mean(cells_ml, na.rm = T),
            '2' = mean(cells_ml, na.rm = T),
            '3'= mean(cells_ml, na.rm = T),
            '4' = mean(cells_ml, na.rm = T),
            '5'= mean(cells_ml, na.rm = T),
            '6' = mean(cells_ml, na.rm = T),
            '7' =mean(cells_ml, na.rm = T),
            '8' =mean(cells_ml, na.rm = T),
            '9'=mean(cells_ml, na.rm = T),
            '10'=mean(cells_ml, na.rm = T),
            '11'=mean(cells_ml, na.rm = T),
            '12'=mean(cells_ml, na.rm = T)) %>%
  gather(key = 'MC', value = 'cells_ml',  -sampling, -species, -phylum)  %>%
  distinct(MC, cells_ml, sampling, species, phylum)

#subdata set 
count_data <- counts %>%
  ungroup()%>%
  dplyr::select(MC, sampling, species, phylum, cells_ml) %>%
 # filter(sampling != 0 & species != 'Rhizosolemia indet.')%>%
  #bind_rows(., spec0) %>%
  left_join(., df, by = c('MC')) %>%
  group_by(MC, sampling, phylum) %>%
  summarise(cell_per_ml = mean(cells_ml, na.rm = T)) %>%
  spread(key = phylum, value = cell_per_ml)


#2. Pigment data 
pigment <- read.csv2("~/Desktop/MA/MA_Rcode/project_data/master_pigments.csv", sep = ";", dec = ',')
str(pigment)

#bring data in form
pigment_data <- pigment %>%
  dplyr::select(-X, -no, -date)%>%
  rename(MC = planktotron)%>%
  mutate(MC = as.character(MC))

#3. join datasets by MC and sampling column
comp_data <- left_join(pigment_data, count_data, by = c('MC', 'sampling')) %>%
  separate(treatment, into = c('mist', 'fluctuation'), '_') 
comp_data$fluctuation[is.na(comp_data$fluctuation)] <-0
comp_data$fluctuation=as.numeric(comp_data$fluctuation)
comp_data$sampling = as.numeric(comp_data$sampling)
comp_data[is.na(comp_data)] <-0 #make sure we don't have NAs


############################## PCA ###############################################
#### Data preparation #### 
# 1. decide which variables are active, remove character variables and turn them into row names 
comp_data.active <- comp_data %>% 
  dplyr::select(-MC, -mist) %>% #remove columns which are not needed anymore
  filter(sampling %in% c(0,6,10,14,18)) %>% #use only data, which we have pigments and coutns
  mutate( id = paste(fluctuation, sampling, sep = '_')) %>% #create dummy column
  dplyr::select(-fluctuation, -sampling) %>%
  gather(key = 'var', value = 'value', -id) %>% 
  group_by(id, var) %>%
  summarise(val = mean(value, na.rm = T)) %>% #calculate mean cells/concentration for each treatment and sampling
  spread(key = var, value = val) %>%
  column_to_rownames("id") %>% #explanatory variables as row names
  dplyr::select(-bb.Car) #remove variable with only 0 entries
  
#### perform PCA #### 
prin_comp <- prcomp(comp_data.active, scale. = T)
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
  rownames_to_column("fluctuation_sampling") %>%
  separate(fluctuation_sampling,into = c("fluctuation", 'sampling'),'_') %>%
  dplyr::select(PC1, PC2, fluctuation, sampling)

##change format to bring variables in the right direction
PCA$fluctuation <- factor(as.factor(PCA$fluctuation),levels=c("0", "48", "36", '24', '12', '6'))
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
#ggsave(plot = last_plot(), file = 'PCA_phylum.png')

############################################################################################################################
#calculate euclidian distance

pca_dat <- pigment_data.active %>%
  dist()%>% ## rawdata
  prcomp(scale=T)%>% 
  broom::tidy() %>% 
  mutate(PC=sprintf("PC%d", PC)) %>% 
  rename(SampleID=row) %>%
  separate(SampleID,into = c("fluctuation", 'sampling'),'_') %>%
  spread(key = PC, value = value) %>%
  dplyr::select(PC1, PC2, fluctuation, sampling)

pca_dat$fluctuation <- factor(as.factor(pca_dat$fluctuation),levels=c("0", "48", "36", '24', '12', '6'))
pca_dat$sampling <- gsub("6", "06", pca_dat$sampling)
pca_dat$sampling =as.numeric(pca_dat$sampling)
pca_dat <- arrange(pca_dat, sampling)

ggplot(pca_dat,aes(x=PC1,y=PC2)) + 
  geom_point(size=1.5,aes(color = fluctuation,  group = fluctuation, shape = as.factor(sampling))) +
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
ggsave(plot = last_plot(), file = 'PCA_euclid_phylum.png')
####################################################################################################


####PCA on pigment data only ####
pigment <- read.csv2("~/Desktop/MA/MA_Rcode/project_data/master_pigments.csv", sep = ";", dec = ',')
str(pigment)

#bring data in form
pigment_data <- pigment %>%
  gather(key = 'pigments',value = 'value',-X, -no, -date,-treatment, -sampling, -planktotron )%>%
  group_by(treatment,sampling, planktotron)%>%
  mutate(sum = sum(value)) %>%
  ungroup()%>%
  spread(pigments, value ) %>%
  mutate(rel_allo = Allo/sum)%>%
  mutate(rel_anth = Anth/sum)%>% ###berechnung der relativen werte in neuer spalte nach Schlueter et al.
  mutate(rel_bb.Car = bb.Car/sum)%>%
  mutate(rel_cNeo = c.Neo/sum) %>%
  mutate(rel_cantha = Cantha/sum)%>%
  mutate(rel_chl.a = Chl.a/sum)%>%
  mutate(rel_chl.b = Chl.b/sum)%>%
  mutate(rel_chl.c1 = Chl.c1/sum)%>%
  mutate(rel_chl.c2 = Chl.c2/sum)%>%
  mutate(rel_diad = Diadino/sum)%>%
  mutate(rel_diato = Dino/sum)%>%
  mutate(rel_dino = Diato/sum)%>%
  mutate(rel_echin = Echin/sum)%>%
  mutate(rel_fuco = Fuco/sum)%>%
  mutate(rel_lut = Lut/sum)%>%
  mutate(rel_myxo = Myxo/sum)%>%
  mutate(rel_peri = Peri/sum)%>%
  mutate(rel_phe.a = Phe.a/sum)%>%
  mutate(rel_viola = Viola/sum) %>%
  mutate(rel_zea = Zea/sum)%>%
  dplyr::select(-c("X","Allo", "Anth","bb.Car", "c.Neo", "Cantha","Chl.a","Chl.c1",
                   "Chl.c2","Cryp","Diadino", "Diato","Dino","Echin","Lut","Myxo","Peri",       
                   "Phe.a","Phe.b", "Viola", "Zea" , 'Fuco', 'Chl.b')) %>%
  gather(key = 'pigment', value = 'value', -planktotron, - sampling, -treatment,-sum, -no, -date) %>%
  group_by(treatment, sampling, pigment) %>% 
  summarise(mean = mean(value, na.rm = T))%>%
  mutate(fluctuation = paste(ifelse(treatment == 'control', 0, ifelse(treatment == 'Fluctuating_48', 48, ifelse(treatment == 'Fluctuating_36', 36,
                                                                                                                ifelse(treatment == 'Fluctuating_24', 24, ifelse(treatment == 'Fluctuating_12', 12, 6)))))))
pigment_data$fluctuation <- factor(as.factor(pigment_data$fluctuation),levels=c("0", "48", "36", '24', '12', '6'))
pigment_data$sampling = as.numeric(pigment_data$sampling)
pigment_data[is.na(pigment_data)] <-0 #make sure we don't have NAs

names(pigment_data)
# 1. decide which variables are active, remove character variables and turn them into row names 
pigment_data.active <- pigment_data %>% 
  ungroup()%>%
  dplyr::select(-treatment) %>% #remove columns which are not needed anymore
  #filter(sampling %in% c(0,6,10,14,18)) %>% #use only data, which we have pigments and coutns
  mutate( id = paste(fluctuation, sampling, sep = '_')) %>% #create dummy column
  dplyr::select(-fluctuation, -sampling) %>%
  spread(key = pigment, value = mean) %>%
  column_to_rownames("id")%>% #explanatory variables as row names
  dplyr::select(-rel_bb.Car) #remove variable with only 0 entries


#### perform PCA #### 
prin_comp <- prcomp(pigment_data.active, scale. = T)
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
  rownames_to_column("fluctuation_sampling") %>%
  separate(fluctuation_sampling,into = c("fluctuation", 'sampling'),'_') %>%
  dplyr::select(PC1, PC2, PC3, PC4, fluctuation, sampling) 
 
 

##change format to bring variables in the right direction
PCA$fluctuation <- factor(as.factor(PCA$fluctuation),levels=c("0", "48", "36", '24', '12', '6'))

PCA$sampling =as.numeric(PCA$sampling)
PCA <- arrange(PCA, sampling)

ggplot(PCA,aes(x=PC1,y=PC2)) + 
  geom_hline(yintercept = 0, linetype = 'dashed', size = 0.5)+
  geom_vline(xintercept = 0, linetype = 'dashed', size = 0.5)+
  geom_point(size=1.5,aes( group = fluctuation , col = as.factor(fluctuation))) +
  #geom_line(arrow = arrow( ends = "both", type = "open"), aes(linetype = fluctuation))+
  geom_path(aes(x = PC1, y = PC2, group = fluctuation, col = as.factor(fluctuation), linetype = as.factor(fluctuation)), 
            arrow = arrow(ends = 'last',type = 'closed',length = unit(0.35, "cm")))+
  scale_color_manual(values = c( '#000000', '#0868ac','#41b6c4','#31a354','#addd8e','#fed976'))+
  scale_y_continuous(limits = c(-6.2, 6.2), breaks = c(-6,-4,-2,0,2,4,6))+
  scale_x_continuous(limits = c(-6.2, 6.2), breaks = c(-6,-4,-2,0,2,4,6))+
  labs(x=paste0("PC1: ",round(prop_varex[1]*100,1),"%"),
       y=paste0("PC2: ",round(prop_varex[2]*100,1),"%")) +
  #facet_wrap(~fluctuation)+
  theme( panel.background = element_rect(fill = NA), #loescht den Hintergrund meines Plots/ fuellt ihn mit nichts
         #panel.grid.major.y = element_line(color='grey', linetype = 'dashed', size=0.2),
         panel.border= element_rect(colour = "black", fill=NA, size=0.5),
         strip.background = element_rect(color ='black', fill = 'white'),
         strip.text = element_text(face = 'bold'),
         legend.background = element_blank(),
         legend.position  ='bottom',
         legend.key = element_blank(),
         text = element_text(size=14))
ggsave(plot = last_plot(), file = 'PCA_all_rel_pigments.png', width = 6, height = 7)

#### calculate Euclidian distance of PC1-4 (Fukami et al. 2005)
dist <- PCA %>%
  mutate(id = paste(fluctuation, sampling, sep ='_')) %>%
  column_to_rownames('id')%>%
  dplyr::select(-fluctuation, -sampling)%>%
  dist() %>%
  broom::tidy() %>% 
  separate(item1, c('fluctuation', 'sampling'), '_')%>%
  separate(item2, c('fluctuation2', 'sampling2'), '_') %>%
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
  labs(x = 'sampling', y = 'distance among treatments')+
  theme_bw()
#ggsave(plot = last_plot(), file = 'distance_rel_pigments.png')
