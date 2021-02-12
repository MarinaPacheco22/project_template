library(BiocManager)
library(STRINGdb)
library(igraph)
library(ggplot2)
library(linkcomm)

setwd("C:/Users/Usuario/OneDrive/Universidad/Curso 4/1º Cuatrimestre/Biología de Sistemas/Proyecto Sars-CoV2")
setwd("D:/UMA/4 curso/Bioinformatica/Biologia de sistemas/Trabajo final/project_template-master")
string_db <- STRINGdb$new(
  version="11",
  species=9606,
  score_threshold=400,
  input_directory="")

#Proteínas humanas que pueden ser relevantes para el SARS-CoV2
proteins.table <- read.delim(file = "media-6.csv", sep = ';', )
proteins.name <- data.frame(proteins.table[3][-1,])
names(proteins.name) <- "gene"

#Mapeo de las proteínas anteriores
proteins_mapped <- string_db$map(proteins.name, "gene", removeUnmappedRows = TRUE)

#Representación de nuestra red
png("proteins_network.png")
proteins_network <- string_db$plot_network(proteins_mapped$STRING_id)
dev.off()

bet = betweenness(graph = proteins_network, normalized = T)
clo = closeness(graph = proteins_network, mode="all", normalized=T)

medias <- unlist(lapply(1:length(bet), FUN = function(x)(bet[x]+clo[x])/2))

order(medias, decreasing = TRUE)
