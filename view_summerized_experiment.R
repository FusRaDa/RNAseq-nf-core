library(SummarizedExperiment)
library(ggplot2)
library(ggpubr)

# The objective of this script is to generate gene expression of all six samples
# and produce 6 tables ranking expresion levels of these genes
# as well as to annotate these genes for further reference

# We will be only using TPM values
# TPM values should only be used to compare gene expression in the SAME sample.


# Load the RDS object
rds_file <- readRDS("~/Projects/RNA-seq-nf-core/results/Homo_sapiens.GRCh38.dna.chromosome.22.fa/star_salmon/null.merged.gene.SummarizedExperiment.rds")

# View available assays
assayNames(rds_file)
colData(rds_file)

# save in variable
tpm = assay(rds_file, "salmon.merged.gene_counts")

# Filter sum of counts per row that is greater than 1000
row_sums <- rowSums(tpm)
tpm_filtered <- tpm[row_sums > 1000, ]
tpm_df <- data.frame(tpm_filtered)

# Create a vertical bar plot first
b1 <- ggplot(tpm_df, aes(x = reorder(row.names(tpm_df), ENCSR000COQ1), y = ENCSR000COQ1, fill = row.names(tpm_df))) +
  geom_bar(stat = "identity") +
  labs(
    title = "ENCSR000COQ1",
    y = NULL,
    x = NULL
  ) +
  theme(plot.title = element_text(hjust = 0.5), legend.title = element_blank()) + 
  coord_flip()

b2 <- ggplot(tpm_df, aes(x = reorder(row.names(tpm_df), ENCSR000COQ2), y = ENCSR000COQ2, fill = row.names(tpm_df))) +
  geom_bar(stat = "identity") +
  labs(
    title = "ENCSR000COQ2",
    y = NULL,
    x = NULL
  ) +
  theme(plot.title = element_text(hjust = 0.5), legend.title = element_blank()) + 
  coord_flip()

b3 <- ggplot(tpm_df, aes(x = reorder(row.names(tpm_df), ENCSR000COR1), y = ENCSR000COR1, fill = row.names(tpm_df))) +
  geom_bar(stat = "identity") +
  labs(
    title = "ENCSR000COR1",
    y = NULL,
    x = NULL
  ) +
  theme(plot.title = element_text(hjust = 0.5), legend.title = element_blank()) + 
  coord_flip()

b4 <- ggplot(tpm_df, aes(x = reorder(row.names(tpm_df), ENCSR000COR2), y = ENCSR000COR2, fill = row.names(tpm_df))) +
  geom_bar(stat = "identity") +
  labs(
    title = "ENCSR000COR2",
    y = NULL,
    x = NULL
  ) +
  theme(plot.title = element_text(hjust = 0.5), legend.title = element_blank()) + 
  coord_flip()

b5 <- ggplot(tpm_df, aes(x = reorder(row.names(tpm_df), ENCSR000CPO1), y = ENCSR000CPO1, fill = row.names(tpm_df))) +
  geom_bar(stat = "identity") +
  labs(
    title = "ENCSR000CPO1",
    y = NULL,
    x = NULL
  ) +
  theme(plot.title = element_text(hjust = 0.5), legend.title = element_blank()) + 
  coord_flip()

b6 <- ggplot(tpm_df, aes(x = reorder(row.names(tpm_df), ENCSR000CPO2), y = ENCSR000CPO2, fill = row.names(tpm_df))) +
  geom_bar(stat = "identity") +
  labs(
    title = "ENCSR000CPO2",
    y = NULL,
    x = NULL
  ) +
  theme(plot.title = element_text(hjust = 0.5), legend.title = element_blank()) + 
  coord_flip()

figure <- ggarrange(b1, b2, b3, b4, b5, b6, ncol = 2, nrow = 3, common.legend = TRUE, legend = "bottom")

annotated_figure <- annotate_figure(
  figure,
  top = text_grob("Sample Top Gene Expression (Raw Counts)", color = "black", face = "bold", size = 16)
)
annotated_figure

