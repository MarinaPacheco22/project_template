if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("STRINGdb", quietly = TRUE))
  install.packages("STRINGdb")

if (!requireNamespace("igraph", quietly = TRUE))
  install.packages("igraph")

if (!requireNamespace("ggplot2", quietly = TRUE))
  install.packages("ggplot2")

if (!requireNamespace("dplyr", quietly = TRUE))
  install.packages("dplyr")

library(xlsx)
library(BiocManager)
library(STRINGdb)
library(igraph)
library(ggplot2)
library(dplyr)

#setwd("C:/Users/Usuario/OneDrive/Universidad/Curso 4/1º Cuatrimestre/Biología de Sistemas/Proyecto Sars-CoV2")
#setwd("D:/UMA/4 curso/Bioinformatica/Biologia de sistemas/Trabajo final/project_template-master")
#setwd(C:/Users/anton/Desktop/trabajoBio)
#setwd(C:/Users/bogda/OneDrive - Universidad de Málaga/4º Ing Salud/Biologia de Sistemas/ProyectoBioSist

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
proteins_network <- string_db$get_subnetwork(proteins_mapped$STRING_id)

png("proteins_network.png")
string_db$plot_network(proteins_mapped$STRING_id)
dev.off()

#Simplificamos la red eliminando los nodos no conectados

proteins_network_components <- components(proteins_network)
nodes_to_remove <- names(proteins_network_components$membership[proteins_network_components$membership!=1])
proteins_network_simple <- delete_vertices(proteins_network, nodes_to_remove)

bet = betweenness(graph = proteins_network_simple, normalized = T)
clo = closeness(graph = proteins_network_simple, mode="all", normalized=T)

medias <- unlist(lapply(1:length(bet), FUN = function(x)(bet[x]+clo[x])/2))

medias_ordenadas <- order(medias, decreasing = TRUE)

#Seleccionamos las 50 proteínas más relevantes de la red.
final_proteins <- c()
final_proteins <- names(unlist(lapply(1:50, FUN = function(x) final_proteins[x] <- medias[medias_ordenadas[x]])))
final_proteins_dataframe <- data.frame(final_proteins)

proteins_mapped_50 <- string_db$map(final_proteins_dataframe, "final_proteins", removeUnmappedRows = TRUE)
proteins_network_50 <- string_db$get_subnetwork(proteins_mapped_50$STRING_id)

png("string_network_50.png")
string_db$plot_network(proteins_mapped_50$STRING_id)
dev.off()

V(proteins_network_50)$label <- NA
V(proteins_network_50)$name <- NA
V(proteins_network_50)$color <- "pink"

png("igraph_network_50.png")
plot(proteins_network_50)
dev.off()


#Top 50 proteínas + Uniprot Function in Disease
tabla_proteins_mapped <- data.frame(proteins_mapped)
tabla_proteins_top50 <- tabla_proteins_mapped[which(tabla_proteins_mapped$STRING_id %in% final_proteins),]
final_table_complete <- proteins.table[which(proteins.table$X.1 %in% tabla_names_id_top50$gene),]
final_table_complete2 <- final_table_complete[final_table_50$X.8 != "",]
final_table <- select(final_table_complete2, X.1, X.8)
#--------------------------

#Proteinas, enfermedades y genes asociados
diseases.prot.associated <- read.delim(file = "selected_prot_disease.csv", sep = ';', )

#Enanismo primordial osteodisplásico microcefálico, tipo II
genes_enanismo <- diseases.prot.associated$Associated.genes[diseases.prot.associated$Disease == 'Microcephalic Osteodysplastic Primordial Dwarfism, Type Ii']
genes_enanismo_vector <- strsplit(genes_enanismo, ", ")[[1]]
genes_enanismo_tabla <- data.frame(genes_enanismo_vector)

string_db_enanismo <- STRINGdb$new(
  version="11",
  species=9606,
  score_threshold=400,
  input_directory="")

proteins_mapped_enanismo <- string_db_enanismo$map(genes_enanismo_tabla, "genes_enanismo_vector", removeUnmappedRows = TRUE)
proteins_network_enanismo <- string_db_enanismo$get_subnetwork(proteins_mapped_enanismo$STRING_id)

png("red_string_enanismo.png")
string_db_enanismo$plot_network(proteins_mapped_enanismo$STRING_id)
dev.off()

V(proteins_network_enanismo)$label <- NA
V(proteins_network_enanismo)$name <- NA
V(proteins_network_enanismo)$color <- "tomato"

png("red_enanismo.png")
plot(proteins_network_enanismo)
dev.off()

#Esquizofrenia 1
genes_esquizo <- diseases.prot.associated$Associated.genes[diseases.prot.associated$Disease == 'Schizophrenia 1']
genes_esquizo_vector <- strsplit(genes_esquizo, ", ")[[1]]
genes_esquizo_tabla <- data.frame(genes_esquizo_vector)

string_db_esquizo <- STRINGdb$new(
  version="11",
  species=9606,
  score_threshold=400,
  input_directory="")

proteins_mapped_esquizo <- string_db_esquizo$map(genes_esquizo_tabla, "genes_esquizo_vector", removeUnmappedRows = TRUE)
proteins_network_esquizo <- string_db_esquizo$get_subnetwork(proteins_mapped_esquizo$STRING_id)

png("red_string_esquizo.png")
string_db_esquizo$plot_network(proteins_mapped_esquizo$STRING_id)
dev.off()

V(proteins_network_esquizo)$label <- NA
V(proteins_network_esquizo)$name <- NA
V(proteins_network_esquizo)$color <- "blue"

png("red_esquizo.png")
plot(proteins_network_esquizo)
dev.off()

#Microcefalia 3, primaria, autosómica recesiva
genes_microcef <- diseases.prot.associated$Associated.genes[diseases.prot.associated$Disease == 'Microcephaly 3, Primary, Autosomal Recessive']
genes_microcef_vector <- strsplit(genes_microcef, ", ")[[1]]
genes_microcef_tabla <- data.frame(genes_microcef_vector)

string_db_microcef <- STRINGdb$new(
  version="11",
  species=9606,
  score_threshold=400,
  input_directory="")

proteins_mapped_microcef <- string_db_microcef$map(genes_microcef_tabla, "genes_microcef_vector", removeUnmappedRows = TRUE)
proteins_network_microcef <- string_db_microcef$get_subnetwork(proteins_mapped_microcef$STRING_id)

png("red_string_microcef.png")
string_db_microcef$plot_network(proteins_mapped_microcef$STRING_id)
dev.off()

V(proteins_network_microcef)$label <- NA
V(proteins_network_microcef)$name <- NA
V(proteins_network_microcef)$color <- "green"

png("red_microcef.png")
plot(proteins_network_microcef)
dev.off

# Las 3 enfermedades juntas
genes_enfermedades <- append(genes_enanismo_tabla[,1], genes_esquizo_tabla[,1])
genes_enfermedades_vector <- append(genes_enfermedades, genes_microcef_tabla[,1])
genes_enfermedades_tabla <- data.frame(genes_enfermedades_vector)

string_db_enfermedades <- STRINGdb$new(
  version="11",
  species=9606,
  score_threshold=400,
  input_directory="")

proteins_mapped_enfermedades <- string_db_enfermedades$map(genes_enfermedades_tabla, "genes_enfermedades_vector", removeUnmappedRows = TRUE)
proteins_network_enfermedades <- string_db_enfermedades$get_subnetwork(proteins_mapped_enfermedades$STRING_id)

string_db_enfermedades$plot_network(proteins_mapped_enfermedades$STRING_id)

nodos <- names(V(proteins_network_enfermedades))
n_enanismo <- which(nodos %in% proteins_mapped_enanismo$STRING_id)
n_esquizo <- which(nodos %in% proteins_mapped_esquizo$STRING_id)
n_microcef <- which(nodos %in% proteins_mapped_microcef$STRING_id)

lista_nodos <- list(n_enanismo, n_esquizo, n_microcef)

V(proteins_network_enfermedades)$label <- NA
V(proteins_network_enfermedades)$name <- NA
V(proteins_network_enfermedades)$color <- "tomato"
V(proteins_network_enfermedades)$size <- 5

png("genes_enfermedades_conjuntas.png")
plot(proteins_network_enfermedades, mark.groups = lista_nodos)
legend(x = "topright", legend = c("Enanismo", "Microcefalia", "Esquizofrenia"), fill = c("red", "blue", "green"), 
       title = "Enfermedades")
dev.off()

#Mapeo de los genes humanos que están asociados a los genes del SARS-CoV2 (final_table_complete)
final_network_vector <- append(genes_enfermedades_vector, final_table_complete$X.1)
final_network_table <- data.frame(final_network_vector)

string_db_final <- STRINGdb$new(
  version="11",
  species=9606,
  score_threshold=400,
  input_directory="")

final_network_proteins_mapped <- string_db_final$map(final_network_table, "final_network_vector", removeUnmappedRows = TRUE)
final_network <- string_db_final$get_subnetwork(final_network_proteins_mapped$STRING_id)
string_db_final$plot_network(final_network_proteins_mapped$STRING_id)

nodos_final <- names(V(final_network))
n_enanismo <- which(nodos_final %in% proteins_mapped_enanismo$STRING_id)
n_esquizo <- which(nodos_final %in% proteins_mapped_esquizo$STRING_id)
n_microcef <- which(nodos_final %in% proteins_mapped_microcef$STRING_id)
n_human <- which(nodos_final %in% tabla_names_id_top50$STRING_id)

lista_nodos_final <- list(n_enanismo, n_esquizo, n_microcef, n_human)

V(final_network)$label <- NA
V(final_network)$name <- NA
V(final_network)$color <- "tomato"
V(final_network)$size <- 5

png("red final.png")
plot(final_network, mark.groups = lista_nodos_final)
legend(x = "topright", legend = c("Enanismo", "Microcefalia", "Esquizofrenia", "Sars-Human"), fill = c("red", "cadetblue2", "green", "purple"), 
       title = "Enfermedades", xjust = 0, yjust = 0)
dev.off()
