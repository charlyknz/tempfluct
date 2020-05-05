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
ggplot(data1, aes(x = sampling, y=c_umol_l, col = treatment, group = MC ))+
  geom_point()+
  geom_line()+
  scale_x_continuous(limits = c(0,18), breaks = seq(0,18,2))+
  scale_y_continuous(limits = c(100,300), breaks = seq(100,300,20))+
  theme_classic()

hist(data1$c_umol_l)
qqnorm(data1$c_umol_l)
qqline(data1$c_umol_l)

## remove data after sampling 10 to exclude zoopl effect ##
data10 <- data1 %>%
  filter(sampling < 11)

##### fixed and random effects ######
# Response variable: Carbon mg (measured over time)
# Explanatory variables: treatments (fluctuation hours) + time
# Random: MC number

#Step 1. Test possible cases:Fit models using REML (use ML in simplifications)
C_m1 = lme(c_umol_l ~ fluctuation*sampling, random = ~1|MC, method = 'REML', data = data10)
C_m2 = lme(c_umol_l ~ fluctuation*sampling, random = ~0+sampling|MC, method = 'REML', data = data10)
C_m3 = lme(c_umol_l ~ fluctuation*sampling, random = ~sampling|MC, method = 'REML',control= lmeControl(niterEM =5000, msMaxIter =5000, msMaxEval =5000), data = data10)
#compare models
anova(C_m1,C_m2, C_m3)

#save fitted values
data10$fit_InterceptOnly2 <- predict(C_m2)


## fit model output ##
ggplot(data10, aes(x = sampling, y=c_umol_l, col = treatment, group = MC ))+
  geom_point()+
  geom_line(aes(y = fit_InterceptOnly2), size = 1) +
  scale_x_continuous(limits = c(0,10), breaks = seq(0,10,2))+
  theme_classic()


# mixed model vs. model without random component: gls 
C_m0=gls(c_umol_l~fluctuation*sampling, method="REML",data =data10, na.action=na.omit)
anova(C_m2, C_m0) #model with random effect is better


#Autocorrelation test (data are not independent) - only if do not have nas
plot(ACF(C_m2), alpha=0.05)

#Residuals
par(mfrow=c(2,2),cex.axis=1.2, cex.lab=1.5)
plot(resid(C_m2, type = "normalized"), ylab="residuales")
hist(resid(C_m2, type = "normalized"), ylab="frecuencia",xlab="residuales", main="")
plot(fitted(C_m2),resid(C_m2, type = "normalized"),ylab="residuales")
qqnorm(resid(C_m2, type = "normalized"), main=""); qqline(resid(C_m2, type = "normalized"))


### LMM with autocorrelation 

#Fit final model with REML
C_m4=lme(c_umol_l ~ fluctuation*sampling, random = ~0+sampling|MC, method = 'REML',
         corr=corAR1(0.1,form=~sampling|MC),data = data10)
  
par(mfrow=c(2,2),cex.axis=1.2, cex.lab=1.5)
plot(resid(C_m4, type = "normalized"), ylab="residuales")
hist(resid(C_m4, type = "normalized"), ylab="frecuencia",xlab="residuales", main="")
plot(fitted(C_m4),resid(C_m4, type = "normalized"),ylab="residuales")
qqnorm(resid(C_m4, type = "normalized"), main=""); qqline(resid(C_m4, type = "normalized"))


qqnorm(C_m4, ~resid(., type
                      ="p")|MC, abline = c(0, 1),
       cex.axis=1.2, cex.lab=1.5,
       pch=20, col="blue", cex=2)

#Corr vs non-corr
anova(C_m2,C_m4) #lower AIC C_m2

summary(C_m2) #only time is signifficant
anova(C_m2)

#Residuals for random effects
qqnorm(C_m2,~ranef(.),abline = c(0, 1),
       cex.axis=1.2, cex.lab=1.5,
       pch=20, col="blue", cex=2)

# Size effects
library(MuMIn)
r.squaredGLMM(C_m4) #

#Simplification:
#Fit model with ML
C_m5=lme(c_umol_l ~ fluctuation*sampling, random=~0+sampling|MC, method="ML",
         corr=corAR1(0.1,form=~sampling|MC),na.action=na.omit,data=data10)

summary(C_m5)

#Delete two-way interaction
C_m6=update(C_m5,~.-sampling:fluctuation)
anova(C_m6,C_m5) #No-significant

C_m7=lme(c_umol_l ~ fluctuation+sampling, random=~0+sampling|MC, method="REML",
         corr=corAR1(0.1,form=~sampling|MC),na.action=na.omit,data=data10)

summary(C_m7)


#save fitted values
data10$fit_InterceptOnly7 <- predict(C_m7)


## fit model output ##
ggplot(data10, aes(x = sampling, y=c_umol_l, col = treatment, group = MC ))+
  geom_point()+
  geom_line(aes(y = fit_InterceptOnly7), size = 1) +
  scale_x_continuous(limits = c(0,10), breaks = seq(0,10,2))+
  #facet_wrap(~treatment)+
  theme_classic()
#ggsave(plot = last_plot(), file = '1_LMM_MA.png')
