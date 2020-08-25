#### ANOSIM AND SIMPER TEST ####
#https://jkzorz.github.io/2019/06/11/ANOSIM-test.html

# install packages
library(tidyverse)
library(vegan)

#### How to ####
# first import your data and bring them in the correct format
# to implement them in your anosim, we need a matrix with our values (species/pigment names) as headers
# second we need a dataframe, which still contains all our explanatory variables such as sampling, treatment usw


#### Step 1: import datasets####
data <- read_delim("~/Desktop/MA/MA_Rcode/project_data/phyto_new.csv", 
                   ";", escape_double = FALSE, col_types = cols(date = col_character()), 
                   locale = locale(decimal_mark = ","), 
                   trim_ws = TRUE) %>%
  drop_na(MC)

#check your importt:
View(data)
str(data)

#change variable format to numeric 
data$cells_ml <- as.numeric(gsub(",", ".", data$cells_ml))
## merge with treatment information
treatments <- read.csv2('~/Desktop/MA/MA_Rcode/project_data/treatments_units.csv')
str(treatments)
names(treatments) = c('MC', 'fluctuation', 'treatment')
rich <- left_join(data, treatments, by = c('MC')) %>%
  separate(species, into = c('genus', 'species'), ' ')


#import measurements 
algal_measurement <- read_delim("~/Desktop/MA/MA_Rcode/project_data/algal_measurement.csv", 
";", escape_double = FALSE, locale = locale(decimal_mark = ","), 
trim_ws = TRUE)
str(algal_measurement)

#calculate mean biovolume (mean of all measurements)
BioV <- algal_measurement %>%
  group_by(genus, species, id) %>%
  summarise(mean_BioV = mean(BioV, na.rm = T),
            sd = sd(BioV, na.rm = T)) %>%
  ungroup()


#### merge all_data ####
all_data <- left_join(rich, BioV, by = c('genus', 'species') ) %>%
  select(treatment, fluctuation, MC, sampling, genus, species, id, cells_ml,mean_BioV) %>%
  mutate(volume = cells_ml*mean_BioV) #calculate biovolume


##### Step 2: get data in the right format ####
# create matrix and dataframe with necessary informations

rel_BV <- all_data%>%
  #filter(sampling ==6)%>%
  mutate(species_id = paste(genus, species, sep = '_'))%>% #create a species id column
  group_by(treatment, fluctuation, sampling, MC)%>% #groups after MC, treatment, sampling etc
  mutate(sum = sum(volume, na.rm = T), #calculate sum and relative contribution of each species/pigment
         rel_V = volume/sum) %>%
  ungroup()%>%
  dplyr::select(sampling, fluctuation, MC, species_id,rel_V) %>% #select only important columns
  drop_na(rel_V) %>% #remove NAs from the dataset
  mutate(interval = 48/fluctuation, #interval and data columns as explanatory variables
         day = sampling *2) %>% 
  spread(key = species_id, value = rel_V) #wide format

rel_BV[is.na(rel_BV)] <-0 #exchange NA with 0
rel_BV$interval[is.infinite(rel_BV$interval)] <-0 #infinite values shall be 0 for interval

#create my data Matrix
Data <- dplyr::select(rel_BV, -sampling,-day, -fluctuation, -MC, -interval)#remove grouping variables
Data <- as.matrix(Data) #create matrix

####Step 3: ANOSIM####
anosim(Data, rel_BV$interval, permutations = 999, distance = "bray", strata = NULL,
       parallel = getOption("mc.cores"))

### interpretation: The divisor is chosen so that R will be in the interval -1 … +1, value 0 indicating completely random grouping.
#An R value close to “1.0” suggests dissimilarity between groups while an R value close to “0” suggests an even distribution of high and low ranks within and between groups”
## R signif greater than 0.05, means that there is no statistical difference between the microbial communities in your groups.

# time as grouping
anosim(Data, rel_BV$day, permutations = 999, distance = "bray", strata = NULL,
       parallel = getOption("mc.cores"))
#An R value close to  “0” suggests an even distribution of high and low ranks within and between groups”
#My significance value is much lower than 0.05
#Therefore, there is a statistically significant difference in my microbial communities based on the grouping “Time”.

simper(Data, rel_BV$day, permutations = 0, trace = FALSE,
       parallel = getOption("mc.cores"))
# greater than 70% is required to say that groups are different from each other



#### Diversity indices ####
# reorder by using select() after every index because the calculation covers column after column

##calculate shannon diversity index
rel_BV$shan = diversity(rel_BV[, -c(1:4)], MARGIN = 1, index='shannon') #new column containing the calculated shannon index

rel_BV <- select(rel_BV, interval, day, sampling, fluctuation, MC, shan,everything() )

#calculate simpson index
rel_BV$simpson = diversity(rel_BV[, -c(1:6)], MARGIN = 1, index='simpson') #new column containing the calculated shannon index

rel_BV <- select(rel_BV, interval, day, sampling, fluctuation, MC, shan, simpson,everything() )

## calculate species richness
absence_presence <- decostand(rel_BV[, -c(1:7)], method= 'pa', na.rm=T) #df giving absence/presence data using decostand function
rel_BV$no = apply(absence_presence, MARGIN = 1, FUN = sum) #new column containing the sum of species present per side (by row = MARGIN = 1)

rel_BV$evenness = rel_BV$shan/log(rel_BV$no)


