library(tidyverse)
library(sp)
library(gridExtra)

human = read.table("~/Documents/sapiens_files/sorted_tab_3_cut_gencodes_exons.bed",
                   fill = T,
                   sep = "\t")

human.s1 = human[human$V4 == "TRAIP;" ,]
human.s2 = human[human$V4 == "FOXP1;" ,]

hgraph = ggplot()+
  geom_density(data = human, aes(x=V2))+
  geom_density(data = human.s1, aes(x=V2, color = 'TRAIP'))+
  geom_density(data = human.s2, aes(x=V2, color = 'FOXP1'))+
  labs(colour = "Genes", x='Position along CHR3', y='Density')+
  theme(legend.position = c(0,1), legend.justification = c(0,1))+
  scale_color_manual(values = c('blue', 'red'))+
  theme_minimal()+
  ggtitle("Human")


nea = read.table("~/Documents/gene_bed_files_chag/genes_chag_3_tab.txt",
                 fill = T,
                 sep = "\t")
nea.s1 = nea[nea$V4 == "TRAIP;" ,]
nea.s2 = nea[nea$V4 == "FOXP1;" ,]

ngraph = ggplot()+
  geom_density(data = nea, aes(x=V2))+
  geom_density(data = nea.s1, aes(x=V2, color = "TRAIP"))+
  geom_density(data = nea.s2, aes(x=V2, color = "FOXP1"))+
  labs(x= ' ', y= ' ')+
  theme_minimal()+
  theme(legend.position = 'none')+
  scale_color_manual(values = c('blue', 'red'))+
  ggtitle("Neanderthal")

grid.arrange(hgraph, ngraph, nrow = 2)
  
