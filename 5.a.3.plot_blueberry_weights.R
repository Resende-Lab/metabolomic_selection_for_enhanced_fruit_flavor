# 
# METABOLOMIC SELECTION FOR ENHANCED FRUIT FLAVOR
# Colantonio and Ferrao et al., 2022
# https://doi.org/10.1073/pnas.2115865119
#



library(ggrepel)
library(ggplot2)
library(RColorBrewer)

weights <- read.csv("blueberryFinalWeights.csv")

weights <- aggregate(weights[,4:5], by = list(weights$metabolite, weights$trait), FUN = mean)
colnames(weights)[1:2] <- c("metabolite", "trait")

nodecols <- read.csv("./bbNetKey.csv")
nodecols <- nodecols[order(match(nodecols[,1], weights$metabolite)),]
nodecols <- nodecols[-c(69,70),]

labs <- list()

#intensity
labs[[1]] <- c("firmness", "2-Undecanone", "Soluble solids", "Ethyl propionate", "Methyl isovalerate", "Caryophyllene", "Valeraldehyde", "Methyl salicylate")

#liking
labs[[2]] <- c("2-Undecanone", "fructose", "Soluble solids", "Eucalyptol", "Hexanal", "Ethyl propionate", "firmness", "Valeraldehyde", "Phenylacetaldehyde","Methyl isovalerate", "Methyl salicylate", "TA", "Linalool")

#sour
labs[[3]] <- c("TA", "Eucalyptol", "2-Hexenyl-acetate", "sucrose", "Linalool", "Neral", "3-Hexenyl-acetate", "2-Hexen-1-al", "citric")

#sweetness
labs[[4]] <- c("TA", "fructose", "Valeraldehyde", "Eucalyptol", "2-Nonanone", "firmness", "2-Undecanone", "Phenylacetaldehyde", "Ethyl propionate", "Soluble solids")


firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

weights$trait <- firstup(as.character(weights$trait))
weights$label <- ""

for(i in 1:4){
  weights$label[weights$metabolite %in% labs[[i]] & weights$trait == unique(weights$trait)[i]] <- as.character(weights$metabolite[weights$metabolite %in% labs[[i]] & weights$trait == unique(weights$trait)[i]])
}


#For true classifications
#weights$color <- rep(nodecols$type, 7)

#For classification by network assigned modules
weights$color <- rep(nodecols$type2, 4)

weights$size <- ((abs(weights$BayesA) * weights$GradientBoostMachine * 0.1) + 1)
weights$size2 <- ((abs(weights$BayesA) * weights$GradientBoostMachine * 0.15) + 8)


p <- ggplot(weights, aes(x = GradientBoostMachine, y = BayesA, color = color, label = (label)))+
  geom_point(size = weights$size)+
  theme_bw()+  
  
  geom_text_repel(size = weights$size2, 
                  min.segment.length = 0,
                  point.padding = 0.25, 
                  segment.color = "grey", 
                  show.legend = F) +
  
  facet_wrap(~trait, ncol = 4) +
  
  theme(strip.background =element_rect(fill=brewer.pal(n = 4, name = "Set3")[8])) +
  scale_color_manual(name = "Groups",values = brewer.pal(n = 12, name = "Set3")[c(8,1,11,6,3:4,7,9:10,5)]) +
  
  theme(strip.text = element_text(size = 24),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 18),
        legend.text = element_text(size = 18), 
        legend.spacing.y = unit(0.5, 'cm'),
        legend.key = element_rect(color = NA, fill = NA),
        legend.key.size = unit(1, "cm"),
        legend.title = element_text(size = 24),
        legend.position = "right")

ggsave(file="bbWeights.svg", plot=p, width=22, height=10)
