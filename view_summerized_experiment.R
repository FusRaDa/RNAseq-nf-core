# Run this code to get libraries
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("SummarizedExperiment")
# install.packages("tidyverse")
# install.packages("ggplot2")
# install.packages("ggpubr")

library(SummarizedExperiment)
library(ggplot2)
library(ggpubr)
library(this.path)

# The objective of this script is to generate gene expression of all six samples
# and produce 6 tables ranking expresion levels of these genes
# as well as to annotate these genes for further reference

# We will be only using TPM values
# TPM values should only be used to compare gene expression in the SAME sample.


# Load the RDS object
setwd(this.dir())
rds_file <- readRDS("results/Homo_sapiens.GRCh38.dna.chromosome.22.fa/star_salmon/null.merged.gene.SummarizedExperiment.rds")

# View available assays
assayNames(rds_file)
colData(rds_file)

# save in variable
rds_gene_counts = assay(rds_file, "salmon.merged.gene_counts")
rds_gene_lengths = assay(rds_file, "salmon.merged.gene_lengths")
rds_gene_counts_length_scaled = assay(rds_file, "salmon.merged.gene_counts_length_scaled")
rds_gene_tpm = assay(rds_file, "salmon.merged.gene_tpm")
rds_gene_counts_scaled = assay(rds_file, "salmon.merged.gene_counts_scaled")

# create function to create graphs
compare_expression <- function(rds_matrix, count_type, filter_count) {
  # Filter sum of counts per row that is greater than 1000
  row_sums <- rowSums(rds_matrix)
  filtered <- rds_matrix[row_sums > filter_count, ]
  rds_df <- data.frame(filtered)
  
  # Create a vertical bar plot first
  b1 <- ggplot(rds_df, aes(x = reorder(row.names(rds_df), ENCSR000COQ1), y = ENCSR000COQ1, fill = row.names(rds_df))) +
    geom_bar(stat = "identity") +
    labs(
      title = "ENCSR000COQ1",
      y = NULL,
      x = NULL
    ) +
    theme(plot.title = element_text(hjust = 0.5), legend.title = element_blank()) + 
    coord_flip()
  
  b2 <- ggplot(rds_df, aes(x = reorder(row.names(rds_df), ENCSR000COQ2), y = ENCSR000COQ2, fill = row.names(rds_df))) +
    geom_bar(stat = "identity") +
    labs(
      title = "ENCSR000COQ2",
      y = NULL,
      x = NULL
    ) +
    theme(plot.title = element_text(hjust = 0.5), legend.title = element_blank()) + 
    coord_flip()
  
  b3 <- ggplot(rds_df, aes(x = reorder(row.names(rds_df), ENCSR000COR1), y = ENCSR000COR1, fill = row.names(rds_df))) +
    geom_bar(stat = "identity") +
    labs(
      title = "ENCSR000COR1",
      y = NULL,
      x = NULL
    ) +
    theme(plot.title = element_text(hjust = 0.5), legend.title = element_blank()) + 
    coord_flip()
  
  b4 <- ggplot(rds_df, aes(x = reorder(row.names(rds_df), ENCSR000COR2), y = ENCSR000COR2, fill = row.names(rds_df))) +
    geom_bar(stat = "identity") +
    labs(
      title = "ENCSR000COR2",
      y = NULL,
      x = NULL
    ) +
    theme(plot.title = element_text(hjust = 0.5), legend.title = element_blank()) + 
    coord_flip()
  
  b5 <- ggplot(rds_df, aes(x = reorder(row.names(rds_df), ENCSR000CPO1), y = ENCSR000CPO1, fill = row.names(rds_df))) +
    geom_bar(stat = "identity") +
    labs(
      title = "ENCSR000CPO1",
      y = NULL,
      x = NULL
    ) +
    theme(plot.title = element_text(hjust = 0.5), legend.title = element_blank()) + 
    coord_flip()
  
  b6 <- ggplot(rds_df, aes(x = reorder(row.names(rds_df), ENCSR000CPO2), y = ENCSR000CPO2, fill = row.names(rds_df))) +
    geom_bar(stat = "identity") +
    labs(
      title = "ENCSR000CPO2",
      y = NULL,
      x = NULL
    ) +
    theme(plot.title = element_text(hjust = 0.5), legend.title = element_blank()) + 
    coord_flip()
  
  figure <- ggarrange(b1, b2, b3, b4, b5, b6, ncol = 2, nrow = 3, common.legend = TRUE, legend = "bottom")
  
  graph_title <- paste("Sample Top Gene Expression (", count_type, ")")
  
  annotated_figure <- annotate_figure(
    figure,
    top = text_grob(graph_title, color = "black", face = "bold", size = 16)
  )
  return(annotated_figure)
}

# plot RDS matrices
rds_gene_counts_graph <- compare_expression(rds_gene_counts, "Raw Counts", 1000)
rds_gene_lengths_graph <- compare_expression(rds_gene_lengths, "Raw Lengths", 100000)
rds_gene_counts_length_scaled_graph <- compare_expression(rds_gene_counts_length_scaled, "Scaled Lengths", 1000)
rds_gene_counts_scaled_graph <- compare_expression(rds_gene_counts_scaled, "Scaled Counts", 1000)
rds_gene_tpm_graph <- compare_expression(rds_gene_tpm, "TPM", 200000)

ggsave("rds_gene_counts_graph.pdf", plot = rds_gene_counts_graph)
ggsave("rds_gene_lengths_graph.pdf", plot = rds_gene_lengths_graph)
ggsave("rds_gene_counts_length_scaled_graph.pdf", plot = rds_gene_counts_length_scaled_graph)
ggsave("rds_gene_counts_scaled_graph.pdf", plot = rds_gene_counts_scaled_graph)
ggsave("rds_gene_tpm_graph.pdf", plot = rds_gene_tpm_graph)
