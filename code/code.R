library(BiocManager)
library(STRINGdb)
library(igraph)
library(ggplot2)

string_db <- STRINGdb$new(
  version="11",
  species=9606,
  score_threshold=995,
  input_directory="")

g = string_db$get_graph()

save(g, file = "string_graph.RData")

string_db$get_png(g$gene, file="string_normal_image.png")
