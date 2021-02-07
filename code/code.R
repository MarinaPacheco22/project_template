library(BiocManager)
library(STRINGdb)
library(igraph)
library(ggplot2)
library(EBImage)

string_db <- STRINGdb$new(
  version="11",
  species=9606,
  score_threshold=400,
  input_directory="")

#Sars-Cov2 Network
img1 = readImage("string_sars_cov2_image.png")
display(img1, method = "raster")

#Proteinas del Sars-Cov2
proteins.sars_cov2.table <- read.delim(file = "C:/Users/Usuario/OneDrive/Universidad/Curso 4/1º Cuatrimestre/Biología de Sistemas/Proyecto Sars-CoV2/string_protein_annotations_sars_cov2.tsv", sep = '\t')
proteins.sars_cov2.name <- data.frame(proteins.sars_cov2.table[1])

#Sars-Human Network
img2 = readImage("string_sars_human_image.png")
display(img2, method = "raster")

#Proteinas del Sars-Human
proteins.sars_human.table <- read.delim(file = "C:/Users/Usuario/OneDrive/Universidad/Curso 4/1º Cuatrimestre/Biología de Sistemas/Proyecto Sars-CoV2/string_protein_annotations_sars_human.tsv", sep = '\t')
proteins.sars_human.name <- data.frame(proteins.sars_human.table[1])