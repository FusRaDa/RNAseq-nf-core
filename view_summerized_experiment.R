library(SummarizedExperiment)
library(ggplot2)

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
tpm = assay(rds_file, "salmon.merged.gene_tpm")

# Filter sum of TPM per row to greater than 0 and turn to a dataframe
row_sums <- rowSums(tpm)
tpm_filtered <- tpm[row_sums != 0, ]
tpm_df <- data.frame(tpm_filtered)

# Create dataframes for all 6 samples
sample_1 <- head(tpm_df[order(tpm_df$ENCSR000COQ1, decreasing = TRUE), ][1])
sample_2 <- head(tpm_df[order(tpm_df$ENCSR000COQ2, decreasing = TRUE), ][2])
sample_3 <- head(tpm_df[order(tpm_df$ENCSR000COR1, decreasing = TRUE), ][3])
sample_4 <- head(tpm_df[order(tpm_df$ENCSR000COR2, decreasing = TRUE), ][4])
sample_5 <- head(tpm_df[order(tpm_df$ENCSR000CPO1, decreasing = TRUE), ][5])
sample_6 <- head(tpm_df[order(tpm_df$ENCSR000CPO2, decreasing = TRUE), ][6])


# Create a vertical bar plot first
ggplot(sample_1, aes(x = reorder(row.names(sample_1), ENCSR000COQ1), y = ENCSR000COQ1)) +
  geom_bar(stat = "identity") +
  
  # Flip the coordinates to make it horizontal
  coord_flip()


# Create a vertical bar plot first
ggplot(sample_2, aes(x = reorder(row.names(sample_2), ENCSR000COQ2), y = ENCSR000COQ2)) +
  geom_bar(stat = "identity") +
  
  # Flip the coordinates to make it horizontal
  coord_flip()



