# 
# METABOLOMIC SELECTION FOR ENHANCED FRUIT FLAVOR
# Colantonio and Ferrao et al., 2022
# https://doi.org/10.1073/pnas.2115865119
#


library(ggrepel)
library(ggplot2)
library(RColorBrewer)

weights <- read.csv("tomatoFinalWeights.csv")

weights <- aggregate(weights[,4:5], by = list(weights$metabolite, weights$trait), FUN = mean)
colnames(weights)[1:2] <- c("metabolite", "trait")

nodecols <- read.csv("./tomNetKey.csv")
nodecols <- nodecols[order(match(nodecols[,1], weights$metabolite)),]
nodecols <- nodecols[-c(69,70),]

labs <- list()

#bitter
labs[[1]] <- c("citric", "glutamic acid", "malic acid", "benzyl cyanide", "salicylaldehyde", "glucose", "fructose", "cis-4-decenal", "guaiacol", "b-ionone")

#intensity
labs[[2]] <- c("citric", "trans-2-pentenal", "Soluble solids", "glucose", "fructose", "trans-2-heptenal", "trans-2-hexenal", "glutamic acid", "benzaldehyde", "guaiacol")

#liking
labs[[3]] <- c("glucose", "fructose", "Soluble solids", "malic acid", "4-carene", "benzyl alcohol", "trans-3-hexen-1-ol", "3-pentenone", "citric","trans-2-pentenal")

#salty
labs[[4]] <- c("citric", "trans-2-pentenal", "Soluble solids", "isovaleric acid", "3-pentenone", "glutamic acid", "3-pentanone")

#sour
labs[[5]] <- c("citric", "malic acid", "2-butylacetate", "propyl acetate", "cis-4-decenal", "fructose", "glutamic acid", "Soluble solids","benzothiazole")

#sweetness
labs[[6]] <- c("glucose", "fructose", "trans-2-pentenal", "4-carene", "trans-2-hexenal", "prenyl acetate", "Soluble solids", "citric", "2-phenyl ethanol")

#umami
labs[[7]] <- c("trans-2-pentenal", "glutamic acid", "1-nitro-2-phenylethane", "benzyl cyanide", "2-phenyl ethanol", "citric", "geranial", "2-methylbutyraldehyde", "trans-3-hexen-1-ol", "trans-2-hexenal")


firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

weights$trait <- firstup(as.character(weights$trait))
weights$label <- ""

for(i in 1:7){
  weights$label[weights$metabolite %in% labs[[i]] & weights$trait == unique(weights$trait)[i]] <- as.character(weights$metabolite[weights$metabolite %in% labs[[i]] & weights$trait == unique(weights$trait)[i]])
}


#For true classifications
#weights$color <- rep(nodecols$type, 7)

#For classification by network assigned modules
weights$color <- rep(nodecols$type2, 7)

weights$size <- ((abs(weights$BayesA) * weights$GradientBoostMachine * 0.1) + 1)
weights$size2 <- ((abs(weights$BayesA) * weights$GradientBoostMachine * 0.2) + 4.5)

p <- ggplot(weights, aes(x = GradientBoostMachine, y = BayesA, color = color, label = (label)))+
  geom_point(size = weights$size)+
  theme_bw()+  
  
  geom_text_repel(size = weights$size2, 
                  min.segment.length = 0,
                  point.padding = 0.25, 
                  segment.color = "grey", 
                  show.legend = F) +
  
  facet_wrap(~trait, ncol = 4) +
  
  theme(strip.background =element_rect(fill=brewer.pal(n = 4, name = "Set3")[4])) +
  scale_color_manual(name = "Groups",values = brewer.pal(n = 12, name = "Set3")[c(8,1,11,6,3:4,7,9:10,5)]) +
  
  theme(strip.text = element_text(size = 24),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 18),
        legend.text = element_text(size = 18), 
        legend.spacing.y = unit(0.5, 'cm'),
        legend.key = element_rect(color = NA, fill = NA),
        legend.key.size = unit(1, "cm"),
        legend.title = element_text(size = 24),
        legend.position = c(0.85, 0.2))

ggsave(file="TomatoWeights.svg", plot=p, width=22, height=10)
