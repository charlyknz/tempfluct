## package
library(tidyverse)

##import data

part_si_data <- read_delim("Desktop/MA/MA_data/abiotics/part_si/part_si_data.csv", 
                     ";", escape_double = FALSE, locale = locale(decimal_mark = ","), 
                     trim_ws = TRUE)
srp_data <- read_delim("Desktop/MA/MA_data/abiotics/SRP /srp_data.csv", 
                     ";", escape_double = FALSE, locale = locale(decimal_mark = ","), 
                      trim_ws = TRUE)
silicate_data <- read_delim("Desktop/MA/MA_data/abiotics/silicate/silicate_data.csv", 
                      ";", escape_double = FALSE, locale = locale(decimal_mark = ","), 
                      trim_ws = TRUE)

POP_data <- read_delim("Desktop/MA/MA_data/abiotics/POP/POP_data.csv", 
                      ";", escape_double = FALSE, locale = locale(decimal_mark = ","), 
                      trim_ws = TRUE)


