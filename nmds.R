#Clear all
rm(list=ls(all=TRUE))

#https://jkzorz.github.io/2019/06/06/NMDS.html
##### Script for ndms plots 
# by Laura Hennigs & Charlotte Kunze
##-----------------------------------------------------------------------------------------------##

#set working directory
setwd("~/Desktop/MA/ndms_script")

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
 # group_by(sampling, treat, species) %>%
  #summarise(mean_cells = mean(cells_ml, na.rm = T)) %>%
  spread(key = 'species', value = 'cells_ml',  fill = NA, convert=F)       #convert the column 'species' to many columns, for each species one. Fill these columns with the cells per ml for that species

#fill all N.A.s with 0
data_nd[4:ncol(data_nd)][is.na(data_nd[4:ncol(data_nd)])] <-0          
str(data_nd)


com = data_nd[, 4:ncol(data_nd)]   #create a data frame from column 4 till 18 (all the species)
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
data <- subset(data.scores, !sampling %in% c( 1:4, 14, 18) )
data <- data.scores
#plot with ggplot
ggplot(data, aes(x = NMDS1, y = NMDS2, colour = treat, shape = as.factor(sampling)))+           #treatments in same colour
  geom_point(size = 4, aes(fill = treat))+
  scale_color_manual(values = c( '#000000','#0868ac','#41b6c4','#31a354','#addd8e','#fed976'))+
  scale_shape_manual(values = c(15,10,17))+
  labs(shape = 'sampling', col = 'Fluctuatin frequency')+  
  annotate('text', x = 0.1, y = 0.25, label = paste('Stress:', round(nmds$stress, 3)))+   #3 gives number of digits
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

