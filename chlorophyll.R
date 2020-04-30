## Script to analyse chlorophyll data
# Jacuzzitron 2019

# by Charlotte Kunze 
# 04.09.2019 
setwd("~/Desktop/MA/MA_data/function/Chlorophyll_a")
#load required packages
library(tidyverse)

# import data file
chlorophylla_data <- read_csv("chlorophylla.csv")
write.csv2(x = chlorophylla_data, file = 'chlorophyll.csv')
#change column structure
chl_data <- chlorophylla_data %>% #new df
  select(-X1) %>% #remove empty row
  select(planktotronNo, treatment_id, everything()) %>%#change column order
  group_by(treatment_id, sampling)%>% #group our data to tell R which variables it should take into account when calculating mean, sd
  mutate(mean_chl = mean(chl, na.rm = T),#new column calculating mean and sd
         sd = sd(chl, na.rm = T),
         se = sd/sqrt(n())) 
chl_control <- chl_data %>%
  filter(treatment == 'control')

chl <- chl_data %>%
  filter(treatment != 'control')
str(chl_data) #check variable format
chl_data$sampling = as.numeric(chl_data$sampling)
#plot 
ggplot(chl_data, aes(x = sampling, y = mean_chl, col = chl_data$treatment_id))+ #aes gives data that should be displayed
  geom_point()+ #we want a scatter ploint so our geom is point
  geom_errorbar(aes(ymin = mean_chl - se, ymax = mean_chl + se), width=.2)+
  geom_line(linetype = 'dashed', size = 0.5)+#lines to connect our points, linetype changes the lines
  #geom_point(data = chl_control, colour = 'black')+
  scale_x_continuous(breaks = seq(0,18,2), limits = c(0,18))+
  labs(x = 'Sampling', y = expression(Chlorophyll~'in'~mu*gL^{-1}), col = 'Treatment')+
  #facet_wrap(~treatment_id)+
  theme_bw()

#ggsave(plot = last_plot(), file = 'chlorophyll_facet_se.png', width = 9, height = 7)


#anova

# only sampling 4 is relevant 

chl_stats <- chl_data %>%
  filter(sampling == 4)

# 2 Bedingungen

#1. Normalverteilung
hist(log10(chl_stats$chl))
qqnorm(chl_stats$chl)
qqline(chl_stats$chl)

# 2. homogenitaet der Varianzen
bartlett.test(chl ~ treatment_id, chl_stats) #p>0.05

aov_chl <- aov(chl ~ treatment_id, chl_stats)
summary(aov_chl)

TukeyHSD(aov_chl)
plot(aov_chl)

## ideas for plots

#ndms with species presence/ absence
# chl over time
# pigment composition in ndms
# stacked barplots with species composition