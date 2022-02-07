# 
# METABOLOMIC SELECTION FOR ENHANCED FRUIT FLAVOR
# Colantonio and Ferrao et al., 2022
# https://doi.org/10.1073/pnas.2115865119
#

library(tidyverse)

####
# Histograms
####

# Set results directory and seed
resdir <- "results/fig1/"
set.seed(100)

# Load data
#key = read_excel("data/input/Metabolites_cluster.xlsx", sheet = 1)
key = read.csv("data/input/tom_metabolites_clusters_key.csv")
tom = tibble(read.csv("data/input/tom_imputed.csv", check.names = F))


# Make unique genotype ids
tom$id <- paste(1:nrow(tom), tom$id)

# Gathering names of traits and metabolites
sensory = c("liking","sweetness", "sour", "umami", "intensity")
mets = colnames(tom)[!colnames(tom) %in% c(sensory, "id")]

# Check for the ones that don't overlap
key$Metabolite[!key$Metabolite %in% mets]
mets[!mets %in% key$Metabolite]

# Creating a metabolite key by merging
idx.tom = data.frame(Metabolite = mets) %>%
  merge(., key, by="Metabolite")

# Selecting only metabolite concentrations
tmp <- tom %>% select(-all_of(sensory))

# Creating the optimal data frame for plotting
tmp2 <- data.frame(t(tmp[, -1]))
colnames(tmp2) <- tmp$id
tmp2$id <- rownames(data.frame(tmp2))

tmp3 <- merge(tmp2, idx.tom[,c("Metabolite", "fig1b_histogram")], by.x = "id", by.y = "Metabolite", all.x = T)

tmp3$fig1b_histogram <- as.factor(tmp3$fig1b_histogram)
colnames(tmp3)[1] <- "metabolite"

tmp4 <- tmp3 %>% #select(-fig1b_histogram) %>% 
  pivot_longer(cols = c(-metabolite, -fig1b_histogram), 
               names_to = "genotype", 
               values_to = "concentration") %>% 
  drop_na() %>% 
  filter(!(fig1b_histogram %in% c("Acid/Sugar", "unknown"))) %>% 
  filter(concentration != 0)


# Plotting
p <- tmp4 %>% 
  ggplot(., aes(x = fig1b_histogram, 
                y = log((concentration), 10), 
                fill = fig1b_histogram)) + 
  geom_violin(color = "black", size = 0.5) +
  ylab(expression(paste(log[10], bgroup("[", ng/gfw/hr, "]")))) +
  scale_y_continuous(limits = c(-3.25, 3.25), 
                     breaks = seq(-3, 3, by = 1)) +
  scale_fill_manual(values = c("#9999ffff", 
                               "#74a9ccff", 
                               "#ffcc99ff", 
                               "#e19582ff", 
                               "#ccccccff"),
                    guide = "none") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 14),
        axis.ticks.x = element_blank()) 


## Saving plot as SVG file

svg(paste0(resdir, "fig1_panelB.svg"), width = 6, height = 2.5)
  plot(p)
dev.off()














































