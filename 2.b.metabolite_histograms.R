# 
# METABOLOMIC SELECTION FOR ENHANCED FRUIT FLAVOR
# Colantonio and Ferrao et al., 2022
# https://doi.org/10.1073/pnas.2115865119
#

library(tidyverse)

####
# Blueberry Histograms 
####


# Set results directory and seed
resdir <- "results/fig2/"
set.seed(100)


# Load data
key = read.csv("data/input/bb_metabolites_clusters_key.csv")
bb = tibble(read.csv("data/input/bb_imputed.csv", check.names = F))

# Make unique genotype ids
bb$id <- paste0("bb", 1:nrow(bb), "_", bb$id)

# Gathering names of traits and metabolites
sensory = c("liking","sweetness", "sour", "intensity")
mets = colnames(bb)[!colnames(bb) %in% c(sensory, "id")]

# Check for the ones that don't overlap
key$Metabolite[!key$Metabolite %in% mets]
mets[!mets %in% key$Metabolite]

# Creating a metabolite key by merging
idx.bb = data.frame(Metabolite = mets) %>%
  merge(., key, by="Metabolite")

# Selecting only metabolite concentrations
tmp <- bb %>% select(-all_of(sensory))

# Creating the optimal data frame for plotting
tmp2 <- data.frame(t(tmp[, -1]))
colnames(tmp2) <- tmp$id
tmp2$id <- rownames(data.frame(tmp2))

tmp3 <- merge(tmp2, idx.bb[,c("Metabolite", "fig2b_histogram")], by.x = "id", by.y = "Metabolite", all.x = T)

tmp3$fig2b_histogram <- as.factor(tmp3$fig2b_histogram)
colnames(tmp3)[1] <- "metabolite"

tmp4 <- tmp3 %>% 
  pivot_longer(cols = c(-metabolite, -fig2b_histogram), 
               names_to = "genotype", 
               values_to = "concentration") %>% 
  drop_na() %>% 
  filter(!(fig2b_histogram %in% c("Acid/Sugar", "unknown"))) %>% 
  filter(concentration != 0)


# Plotting
p <- tmp4 %>% 
  ggplot(., aes(x = fig2b_histogram, 
                y = log((concentration), 10), 
                fill = fig2b_histogram)) + 
  geom_violin(color = "black", size = 0.5) +
  ylab(expression(paste(log[10], bgroup("[", ng/gfw/hr, "]")))) +
  scale_y_continuous(limits = c(-2.25, 4.25), 
                   breaks = seq(-2, 4, by = 1)) +
  scale_fill_manual(values = c("#9999ffff",
                               "#74a9ccff", 
                               "#9fccc2ff",
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

svg(paste0(resdir, "fig2_panelB.svg"), width = 6, height = 2.5)
plot(p)
dev.off()

































