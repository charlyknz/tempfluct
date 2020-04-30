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

#import data using the manual import function (readr) and after that copying the console output
tidy_counts <- read_delim("~/Desktop/MA/Data/cell_counts_sampling3.csv", 
                          ";", escape_double = FALSE, col_types = cols(date = col_date(format = "%d.%m.%Y")), 
                           trim_ws = TRUE)

# new data frame including data wrangling
data_nd <- tidy_counts %>%
  select(planktotron_no, species, sampling, date, cells_ml) %>%            #select columns we want to keep
  distinct(planktotron_no, species, sampling, date, cells_ml) %>%          #keep only distinct combinations of given columns
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
data.scores$date = data_nd$date
data.scores$planktotron_no = data_nd$planktotron_no

#write a data frame with treatment information, which has the same column planktotron_no as in data.scores
df=data.frame(treat=c('constant','fluct 36','fluct 6','constant','fluct 6','fluct 48','fluct 48',
                      'fluct 24','fluct 24','fluct 36','fluct 12','fluct 12'),  
              planktotron_no = c(1:12) )

#combine the two df
data.scores <- full_join(data.scores, df, by = c('planktotron_no'))

head(data.scores)  #see our results (but only the first rows)


#plot with ggplot
ggplot(data.scores, aes(x = NMDS1, y = NMDS2, colour = treat))+           #treatments in same colour
  geom_point(size = 4)+
  scale_y_continuous(breaks = seq(-0.5, 1.0, 0.25), limits=c(-0.51, 1.0))+
  annotate('text', x = 0.5, y = 0.6, label = paste('Stress:', round(nmds$stress, 3)))+   #3 gives number of digits
  theme_bw()                                                          
  
## stress
# stress < 0.05 excellent representation
# stress < 0.1 great
# stress < 0.2 okay
# stress < 0.3 poor represented
# stress > 0.3 choose another method