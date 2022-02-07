# 
# METABOLOMIC SELECTION FOR ENHANCED FRUIT FLAVOR
# Colantonio and Ferrao et al., 2022
# https://doi.org/10.1073/pnas.2115865119
#

###
# Figure 4c Metabolomic Selection Models
###


library(ggplot2)
library(dplyr)
library(RColorBrewer)


pal=pnw_palette("Shuksan2",100)
pal.1=c(pnw_palette("Bay",6),"gray")
pal.1 = pal.1[-3] 

cols <- c("#50fa7b", "#ffb86c", "#bd93f9", "#ff79c6", 
          "#ff5555", "#f1fa8c", "#6272a4", "#8be9fd", 
          "#282a36", "#44475a", "#44475a", "#f8f8f2")

drac <- c("#50fa7b", "#ffb86c", "#bd93f9", "#ff79c6", 
          "#ff5555", "#f1fa8c", "#6272a4", "#8be9fd", 
          "#282a36", "#44475a", "#44475a", "#f8f8f2")

V <- theme_bw() + 
  theme(#legend.position = c(0.75, 0.75),  
    panel.grid.major.x  = element_blank(),   
    panel.grid.minor.x  = element_blank(),
    panel.grid.major.y  = element_blank(),   
    panel.grid.minor.y  = element_blank(),       
    #      plot.title = element_text(hjust = 0.5),
    #        plot.background = element_rect(fill = cols[9], color = cols[9]),
    #        panel.background = element_rect(fill = cols[9]),
    #        legend.background = element_rect(fill = cols[9]),
    #        axis.text = element_text(color = cols[12]),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    text = element_text(size = 18),
    #   text = element_text(size = 18, face = "bold"),
    #   text = element_text(color = cols[12], size = 18, family="Harding Text Web Regular"),
    #   panel.border = element_rect(color = cols[12]),
    #   axis.ticks = element_line(color = cols[12]),
    #   legend.key = element_rect(fill = cols[9], color = cols[9])
  ) 


df <- read.table("resultsSub.txt", header = T)
df <- df[df$experiment != "experiment",]

#str(df)

df$experiment <- as.numeric(df$experiment)
df$subSize <- as.numeric(df$subSize)
df$accuracy <- as.numeric(df$accuracy)

#df <- df[df$trait != "umami",]
df <- df[df$trait == "umami",]

#df$trait <- factor(df$trait, levels = c("intensity", "sour", "sweetness", "liking"))
df$trait <- factor(df$trait, levels = c("umami"))
df$model <- as.factor(df$model)

levels(df$trait) <-  c("Umami")

levels(df$trait) <-  c("Flavor\nIntensity", "Sourness\nIntensity", "Sweetness\nIntensity", "Overall\nLiking")

p <- df[df$model %in% c( "GradientBoostMachine", "gBLUP", "BayesA"),] %>%
  group_by(model, trait, subSize) %>%
  summarize(accuracy = mean(accuracy)) %>%
  ggplot(aes(x = subSize, y = accuracy, color = model, group = model)) + 
  geom_line(size = 1) + theme_bw() +
  scale_color_manual(values = brewer.pal(12, "Set3")[c(10,6,7)], name = "") +
  scale_x_continuous(breaks = seq(50, 170, by = 40)) +
  ylim(c(0.2, 0.7)) +
  facet_wrap(~trait, ncol = 4) + 
  V + xlab("# of Samples") + ylab("Accuracy") + theme(legend.position = "none")

p

subsamp <- p

svg("umamiSub.svg", width = 4, height = 4)
p
dev.off()


svg("subSamp.svg", width = 8, height = 4)
p
dev.off()




res <- df %>%
  group_by(model, trait, experiment, type) %>%
  summarize(accuracy = mean(accuracy)) 


res2 <- df %>% 
  pivot_wider(names_from = model, values_from = accuracy) %>% 
  arrange(experiment, trait, subSize)  

write.csv(res2, "ST7_subSampling_tomato.csv", row.names = F)




p <- df %>%
  group_by(model, trait, subSize) %>%
  summarize(accuracy = mean(accuracy)) %>%
  ggplot(aes(x = subSize, y = accuracy, color = model, group = model)) + geom_line() + facet_wrap(~trait)

pdf("subSampFull.pdf", width = 10, height = 8)
p
dev.off()


df[df$subSize == 170,] %>%
  group_by(model) %>%
  summarize(accuracy = mean(accuracy))  %>%
  arrange(accuracy)




































