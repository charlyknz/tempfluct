## package
library(tidyverse)

##import data

part_si_data <- read_delim("~/Desktop/MA/MA_Rcode/project_data/part_si_data.csv", 
                           ";", escape_double = FALSE, locale = locale(decimal_mark = ","), 
                           trim_ws = TRUE)
srp_data <- read_delim("~/Desktop/MA/MA_Rcode/project_data/srp_data.csv", 
                     ";", escape_double = FALSE, locale = locale(decimal_mark = ","), 
                      trim_ws = TRUE)
silicate_data <- read_delim("~/Desktop/MA/MA_Rcode/project_data/silicate_data.csv", 
                      ";", escape_double = FALSE, locale = locale(decimal_mark = ","), 
                      trim_ws = TRUE)

POP_data <- read_delim("~/Desktop/MA/MA_Rcode/project_data/POP_data.csv", 
                      ";", escape_double = FALSE, locale = locale(decimal_mark = ","), 
                      trim_ws = TRUE)

chlorophyll <- read_delim("~/Desktop/MA/MA_Rcode/project_data/chlorophyll.csv", 
                                   ";", escape_double = FALSE, col_types = cols(X1 = col_skip(), 
                                    X1_1 = col_skip()), locale = locale(decimal_mark = ","), 
                                   trim_ws = TRUE) %>%
  rename('Planktotron' = 'planktotronNo', 'Sampling' ='sampling') %>%
  mutate(Sampling = as.numeric(Sampling)) 

chlA <- chlorophyll %>%
  group_by(Sampling, Planktotron) %>%
  summarise(chl_a = mean(chl, na.rm = T))%>%
  left_join(., srp, by = c('Planktotron', 'Sampling')) %>%
  select(no, Sampling, Planktotron, date, chl_a)
#write.csv2(x = chlA, file = 'invivo_Chlorophyll_no.csv')

part_si <- part_si_data %>%
  group_by(no, Planktotron,Sampling, date) %>%
  summarise(mean_part_si = mean(part_si_conc_ug_l,na.rm = T)) 

srp <- srp_data %>%
  group_by(no, Planktotron,Sampling, date) %>%
  summarise(mean_srp = mean(srp_conc_mg_l,na.rm = T))

silicate <- silicate_data %>%
  group_by(no, Planktotron,Sampling, date) %>%
  summarise(mean_si = mean(Si_conc_mg_l, na.rm = T)) %>%
  select(-date)

pop <- POP_data %>%
  group_by(no, Planktotron,Sampling, date) %>%
  summarise(mean_POP = mean(POP_conc_mg_l, na.rm = T))

all_nutrients <- full_join(silicate, srp, by = c('no', 'Planktotron', 'Sampling' ))

all_nutrients <- full_join(all_nutrients, pop, by = c('no','Planktotron', 'Sampling', 'date'))
all_nutrients <- full_join(all_nutrients, part_si, by = c('no','Planktotron', 'Sampling', 'date'))

#write.csv2(x = all_nutrients, file = 'nutrients_without_pseudoreplicates.csv')
  