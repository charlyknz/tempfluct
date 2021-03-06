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
         MC = str_remove(MC, 'P'),
         day = 2*sampling) %>%
  mutate(fluctuation = ifelse(treatment == 'control', 0, ifelse(treatment == 'Fluctuating_12', 12, ifelse(treatment == 'Fluctuating_24', 24, ifelse(treatment == 'Fluctuating_36', 36, ifelse(treatment == 'Fluctuating_48', 48,6 )))))) %>%
  select(-X) %>%
  mutate(dayname = as.factor(day)) %>%
  mutate(interval = 48/fluctuation)
data1$interval[!is.finite(data1$interval)] <- 0 
str(data1)

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

##### fixed and random effects ######
# Response variable: Carbon mg (measured over time)
# Explanatory variables: treatments (fluctuation hours) + time
# Random: MC number

#Step 1. Test possible cases:Fit models using REML (use ML in simplifications)
C_m1 = lme(c_umol_l ~ interval*day, random = ~1|MC, method = 'REML', data = data1)
C_m2 = lme(c_umol_l ~ interval*day, random = ~0+day|MC, method = 'REML', data = data1)
C_m3 = lme(c_umol_l ~ interval*day, random = ~day|MC, method = 'REML',control= lmeControl(niterEM =5000, msMaxIter =5000, msMaxEval =5000), data = data1)


##model using dayname and MC as random effects

#compare models
anova(C_m1,C_m2, C_m3)

#save fitted values
data1$fit_InterceptOnly2 <- predict(C_m2)


## fit model output ##
data1$fluctuation <- factor(data1$fluctuation, levels = c('0', '48','36','24','12','6'))
ggplot(data1, aes(x = fluctuation, y=c_umol_l, col = treatment ))+
  geom_point(aes(alpha = sampling))+
  #geom_line(aes(y = fit_InterceptOnly2), size = 1) +
  #scale_x_continuous(limits = c(0,10), breaks = seq(0,10,2))+
  theme_classic()


# mixed model vs. model without random component: gls 
C_m0=gls(c_umol_l~interval*day, method="REML",data =data1, na.action=na.omit)
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
C_m4=lme(c_umol_l ~ interval*day, random = ~0+day|MC, method = 'REML',
         corr=corAR1(0.1,form=~day|MC),data = data1)
  
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

summary(C_m4) #only time is signifficant
anova(C_m4)

#Residuals for random effects
qqnorm(C_m2,~ranef(.),abline = c(0, 1),
       cex.axis=1.2, cex.lab=1.5,
       pch=20, col="blue", cex=2)

# Size effects
library(MuMIn)
r.squaredGLMM(C_m4) #

#Simplification:
#Fit model with ML
C_m5=lme(c_umol_l ~ interval*day, random=~0+day|MC, method="ML",
         corr=corAR1(0.1,form=~day|MC),na.action=na.omit,data=data1)

summary(C_m5)

#Delete two-way interaction
C_m6=update(C_m5,~.-day:interval)
anova(C_m6,C_m5) #No-significant

C_m7=lme(c_umol_l ~ interval+day, random=~0+day|MC, method="REML",
         corr=corAR1(0.1,form=~day|MC),na.action=na.omit,data=data1)

summary(C_m7)


#save fitted values
data1$fit_InterceptOnly7 <- predict(C_m7)


## fit model output ##
ggplot(data1, aes(x = day, y=c_umol_l, col = treatment, group = MC))+
  geom_point()+
  geom_line(aes(y = fit_InterceptOnly7), size = 1) +
  scale_x_continuous(limits = c(0,10), breaks = seq(0,10,2))+
  #facet_wrap(~treatment)+
  theme_classic()
#ggsave(plot = last_plot(), file = '1_LMM_MA.png')

###########################################################################
# Helmuts Model
library(lmerTest)

mod1 <- lmer(c_umol_l ~ interval*day + (1|MC) + (1|dayname), data=data1)
summary(mod1)
anova(mod1)

par(mfrow=c(2,2),cex.axis=1.2, cex.lab=1.5)
plot(resid(mod1), ylab="residuales")
hist(resid(mod1), ylab="frecuencia",xlab="residuales", main="")
plot(fitted(mod1),resid(mod1),ylab="residuales")
qqnorm(resid(mod1), main=""); 
qqline(resid(mod1))


###########################################################################
##### Explore Zooplankton Phytoplankton realtionship #####
library(ggpubr)

#import data
zoo_data <- read.csv2('~/Desktop/MA/MA_Rcode/project_data/Zoo_all.csv')
str(zoo_data)


zoo <- zoo_data %>%
  select(n_sampling, planktotron, fluct, chl_a,C_water_µmol_l, overall_abundance_l, dryweight_l ) %>%
  rename(sampling = n_sampling, MC = planktotron, treatment = fluct, c_umol_l = C_water_µmol_l) %>%
  mutate(fluctuation = ifelse(treatment == 'Con', 0, ifelse(treatment == 'F12', 12, ifelse(treatment == 'F24', 24,
                                                       ifelse(treatment == 'F36', 36, ifelse(treatment == 'F48', 48,6 )))))) %>%
  drop_na(overall_abundance_l) 

## correlation scatter
ggscatter(zoo, x = "dryweight_l", y = "c_umol_l", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson", xlab = 'zooplankton biomass',
          ylab = 'phytoplankton carbon content')
  
## per treatment 
ggscatter(zoo, x = "overall_abundance_l", y = "c_umol_l", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson", xlab = 'zooplankton abundance',
          color = "treatment", palette = "jco",           # Color by groups "cyl"
          shape = "treatment", facet.by = 'treatment',
          ylab = 'phytoplankton carbon content')+
  stat_cor(aes(color = treatment), label.x = 3) 

ggsave(plot = last_plot(), file = 'MA_correlation_plot_1.png')
res <- cor.test(zoo$dryweight_l,zoo$c_umol_l, 
                method = "pearson")
res
head(zoo)


## see also http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/78-perfect-scatter-plots-with-correlation-and-marginal-histograms/

