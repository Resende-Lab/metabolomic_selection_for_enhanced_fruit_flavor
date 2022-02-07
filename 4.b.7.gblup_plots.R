# 
# METABOLOMIC SELECTION FOR ENHANCED FRUIT FLAVOR
# Colantonio and Ferrao et al., 2022
# https://doi.org/10.1073/pnas.2115865119
#

###
# Figure 4b Metabolomic Selection Models
###

library(ggplot2)
library(dplyr)
library(tidyr)

###############
V <- theme_bw() + 
  theme(#legend.position = c(0.75, 0.75),  
    panel.grid.major.x  = element_blank(),   
    panel.grid.minor.x  = element_blank(),
    panel.grid.major.y  = element_blank(),   
    panel.grid.minor.y  = element_blank(),       
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    text = element_text(size = 18),
    #   text = element_text(size = 18, face = "bold"),
    #   text = element_text(color = cols[12], size = 18, family="Harding Text Web Regular"),
    #   panel.border = element_rect(color = cols[12]),
    #   axis.ticks = element_line(color = cols[12]),
    #   legend.key = element_rect(fill = cols[9], color = cols[9])
  ) 



gen <- read.table("rfGenomic.txt", header = T)
meta <- read.table("rfMetabolomicRun1.txt", header = T)

gen$type <- "genomic"
meta$type <- "metabolomic"

df <- rbind(gen, meta)
#str(df)

df <- df[df$trait != "umami",]
df <- df[df$model != "randomForest",]

df$experiment <- as.numeric(df$experiment)
df$cv <- as.numeric(df$cv)
df$accuracy <- as.numeric(df$accuracy)

df$trait <- factor(df$trait, levels = c("intensity", "sour", "sweetness", "liking"))
#df$trait <- as.factor(df$trait)
df$model <- as.factor(df$model)
df$type <- as.factor(df$type)

write.csv(df, "ST6_genomic_tom.csv", row.names = F, quote = F)


df <- read.csv("data/supplemental_datasets/SD6_gBLUP_tomato.csv", header = T)

levels(df$trait) <-  c("Flavor\nIntensity", "Sourness\nIntensity", "Sweetness\nIntensity", "Overall\nLiking")
#levels(df$trait) <- levels(df$trait)]

df %>%
  group_by(model, trait, experiment) %>%
  summarize(accuracy = mean(accuracy)) %>%
  ggplot(aes(x = experiment, y = accuracy, color = model, group = model)) + geom_col() + facet_wrap(~trait+ model)

p <- df %>%   group_by(model, trait, experiment, type) %>%
  summarize(accuracy = mean(accuracy)) %>%
  ggplot(aes(x = type, y = accuracy, fill = interaction(model, type))) + 
  geom_violin(position = "dodge") + 
  theme_bw() + 
  ylim(c(-0.5, 1.2)) +
  scale_y_continuous(breaks = seq(-0.2, 1.0, by = 0.2)) +
  scale_fill_manual(values = c("#8cb5ec83", "#dc9074b7"), labels = c("Genomic", "Metabolomic"), name = "") +
  facet_wrap(~trait, ncol = 4) + V + theme(legend.position = "none") + ylab("Accuracy") + xlab("") + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

pid

gblup <- p

svg("gblup.svg", width = 8, height = 4)
p
dev.off()

svg("fig4bc.svg", width = 8, height = 8)
plot_grid(gblup, subsamp, align = "hv", ncol = 1)
dev.off()





#pdf(width = 10, height = 4)
#p
#dev.off()

res <- df %>%
  group_by(model, trait, experiment, type) %>%
  summarize(accuracy = mean(accuracy)) 

df2 <- data.frame(res)
df2[df2$model == "gBLUP",]
df2[df2$experiment == 1 & df2$type == "metabolomic",]

df2 %>% select(experiment = 1)



res2 <- res %>% 
  pivot_wider(names_from = trait, values_from = accuracy) %>% 
  arrange(type)  

write.csv(res2, "ST6_gBLUP_tomato.csv", row.names = F)

res %>%     
  group_by(model, trait, type) %>% 
  summarize(accuracy2 = mean(accuracy)) %>% 
  pivot_wider(names_from = trait, values_from = accuracy) %>% 
  arrange(type) 










