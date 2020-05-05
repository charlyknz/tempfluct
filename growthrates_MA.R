
#####growth models ####
library(growthrates)
#https://tpetzoldt.github.io/growthrates/doc/Introduction.html

# import temperature data 
data <- read.csv2('~/Desktop/MA/MA_Rcode/project_data/cnp_data_treatments.csv')
names(data)
data1 <- data %>%
  mutate(sampling = str_remove(sampling, 'S'),
         sampling = as.numeric(sampling),
         MC = str_remove(MC, 'P')) %>%
  mutate(fluctuation = ifelse(treatment == 'control', 0, ifelse(treatment == 'Fluctuating_12', 12, ifelse(treatment == 'Fluctuating_24', 24, ifelse(treatment == 'Fluctuating_36', 36, ifelse(treatment == 'Fluctuating_48', 48,6 )))))) %>%
  select(-X)%>%
  filter(sampling <12) 

p     <- c(y0 = 5, mumax = 0.5, K = 200)
lower <- c(y0 = -20, mumax = 0,   K = 10)
upper <- c(y0 = 200, mumax = 5,   K = 300)

all_logistic_fits <- all_growthmodels(
  c_umol_l ~ grow_logistic(sampling, parms) | treatment ,
  data = data1,
  p = p, lower = lower, upper = upper,
  log = "y", ncores=1)
#sets up the framing for the plot, and plots every individual growth curve with logistic fits

par(mfrow = c(4, 3))
par(mar=c(1,1,1,1))
plot(all_logistic_fits)
#dev.off()

#takes the parameter estimates from the previous bit and puts them into a table:
all_logistic_results <- results(all_logistic_fits) %>% 
  mutate(treatment = fct_recode(as.factor(treatment), "con" = "control",
                                          "F6" = "Fluctuating_6",
                                          "F12" = "Fluctuating_12",
                                          "F24" = "Fluctuating_24",
                                          "F36" = "Fluctuating_36",
                                          "F48" = "Fluctuating_48"))
#write.csv2(x = all_logistic_results, file = 'growthrates.csv')
