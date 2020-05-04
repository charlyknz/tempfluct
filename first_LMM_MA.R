## Model formulation for thesis analysis

# 1. load required packages
library(lme4)
library(nlme)
library(tidyverse)

#import your data
data <- read.csv2('~/Desktop/MA/MA_Rcode/project_data/cnp_data_treatments.csv')
names(data)
data1 <- data %>%
  mutate(sampling = str_remove(sampling, 'S'),
         sampling = as.numeric(sampling),
         MC = str_remove(MC, 'P')) %>%
  mutate(fluctuation = ifelse(treatment == 'control', 0, ifelse(treatment == 'Fluctuating_12', 12, ifelse(treatment == 'Fluctuating_24', 24, ifelse(treatment == 'Fluctuating_36', 36, ifelse(treatment == 'Fluctuating_48', 48,6 )))))) %>%
  select(-X)

#plot your data
ggplot(data1, aes(x = sampling, y= C.mg, col = treatment ))+
  geom_point()+
  scale_x_continuous(limits = c(0,18), breaks = seq(0,18,2))

hist(data1$C.mg)
qqnorm(data1$C.mg)
qqline(data1$C.mg)

##### fixed and random effects ######
# Response variable: Carbon mg (measured over time)
# Explanatory variables: treatments (fluctuation hours) + time
# Random: MC number

#Step 1. Test possible cases:Fit models using REML (use ML in simplifications)
C_m1 = lme(C.mg ~ fluctuation*sampling, random = ~1|MC, method = 'REML', data = data1)
C_m2 = lme(C.mg ~ fluctuation*sampling, random = ~0+sampling|MC, method = 'REML', data = data1)
C_m3 = lme(C.mg ~ fluctuation*sampling, random = ~sampling|MC, method = 'REML',control= lmeControl(niterEM =5000, msMaxIter =5000, msMaxEval =5000), data = data1)
anova(C_m1, C_m2, C_m3)

