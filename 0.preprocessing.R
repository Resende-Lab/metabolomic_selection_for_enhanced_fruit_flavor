# 
# METABOLOMIC SELECTION FOR ENHANCED FRUIT FLAVOR
# Colantonio and Ferrao, et al., 2022
# https://doi.org/10.1073/pnas.2115865119
#

###
# Initial Preprocessing
###

library(tidyverse)
library(readxl)

###
# Tomato Data
###

# Extract SD1_Original_Data sheet from SD1_dataset_tomato.xlsx file
tom_og <- read_excel("data/supplemental_datasets/SD1_dataset_tomato.xlsx", sheet = "SD1_Original_Data")

# Remove columns not needed for most computation
tom_og <- tom_og %>% select(-species, -`panel number`)

# Scale all features and responses to N(0,1)
tom_og[,2:length(tom_og)]  <- scale(tom_og[,2:length(tom_og)])

# Write new scaled data to .csv file for clean input to other scripts
write.csv(tom_og, "./data/input/tom_imputed_scaled.csv", row.names = F)


###
# Verify data scaled here is same as in SD1_Imputed_Scaled
tom_imputed_scaled = read.csv("./data/input/tom_imputed_scaled.csv", header = T, check.names = F)
sd1_tom_imputed_scaled = read_excel("data/supplemental_datasets/SD1_dataset_tomato.xlsx", sheet = "SD1_Imputed_Scaled")

sd1_tom_imputed_scaled <- sd1_tom_imputed_scaled %>% select(-species, -`panel number`)

# Find the differences between the two datasets
differences <- tom_imputed_scaled[,2:ncol(tom_imputed_scaled)] - 
                sd1_tom_imputed_scaled[,c(2:ncol(sd1_tom_imputed_scaled))]

# Plot out the differences to make sure they're mostly 0 or negligibly small 
## (negligible values (e.g. 10^-16) could be due to random differences between computer processors)
differences %>% pivot_longer(cols = colnames(differences), 
                             names_to = "metabolite", 
                             values_to = "difference") %>% 
  ggplot(., aes(x = difference)) + geom_histogram(bins = 100)


# Although the differences are small, to try and be as reproducible as possible, 
#   we will use the scaled data from the SD1_dataset_tomato.xlsx file
write.csv(sd1_tom_imputed_scaled, "./data/input/tom_imputed_scaled.csv", row.names = F)

#~~~





####
# Blueberry Data
####

# Extract SD2_Original_Data sheet from SD2_dataset_blueberry.xlsx file
bb_og <- read_excel("data/supplemental_datasets/SD2_dataset_blueberry.xlsx", sheet = "SD2_Original_Data", na = "NA")

# Remove columns not needed for most computation
bb_og <- bb_og %>% select(-species)

# Make id column name the same as in the tomato dataset
colnames(bb_og)[1] <- "id"

# For the following metabolites, not every sample in the population has been quantified (e.g. missing data)
mets_with_missing_data <- c("firmness",
                            "citric",
                            "fructose",
                            "glucose",
                            "sucrose",
                            "Limonene",
                            "Nerylacetone")

# Thus, we will impute any missing data with the mean of the non-missing samples in the population

# Calculate means of non-missing samples for each metabolite containing missing data
means <- bb_og %>% 
  select(mets_with_missing_data) %>% 
  colMeans(na.rm = T)

# For each metabolite with missing data we find which samples are NAs 
#   and replace the NAs with the mean of the non-NA samples

for(i in 1:length(mets_with_missing_data)){
  bb_og[is.na(bb_og[, mets_with_missing_data[i]]), mets_with_missing_data[i]] <- means[i]
}

# Write new imputed data to .csv file
write.csv(bb_og, "./data/input/bb_imputed.csv", row.names = F)

# We then scale all features and responses to N(0,1)
bb_og[,2:ncol(bb_og)]  <- scale(bb_og[,2:ncol(bb_og)])

# Write new imputed scaled data to .csv file for clean input to other scripts
write.csv(bb_og, "./data/input/bb_imputed_scaled.csv", row.names = F)

###

###
# Verify data scaled here is same as in SD2_Imputed_Scaled
bb_imputed_scaled = read.csv("./data/input/bb_imputed_scaled.csv", header = T, check.names = F)
sd2_bb_imputed_scaled = read_excel("data/supplemental_datasets/SD2_dataset_blueberry.xlsx", sheet = "SD2_Imputed_Scaled")

sd2_bb_imputed_scaled <- sd2_bb_imputed_scaled %>% select(-species)
colnames(sd2_bb_imputed_scaled)[1] <- "id"

# Find the differences between the two datasets
differences <- bb_imputed_scaled[, 2:ncol(bb_imputed_scaled)] - 
  sd2_bb_imputed_scaled[,c(2:ncol(sd2_bb_imputed_scaled))]

# Plot out the differences to make sure they're small
differences %>% pivot_longer(cols = colnames(differences), 
                             names_to = "metabolite", 
                             values_to = "difference") %>% 
  ggplot(., aes(x = difference)) + geom_histogram(bins = 100)


# Although the differences are small, to try and be as reproducible as possible, 
#   we will use the scaled data from the SD1_dataset_tomato.xlsx file
write.csv(sd2_bb_imputed_scaled, "./data/input/bb_imputed_scaled.csv", row.names = F)

#~~~

#~~~








































