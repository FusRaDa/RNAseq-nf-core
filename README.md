## NF-CORE/RNA-SEQ pipeline

This notebook was written in Google Collab and intended to be ran in that environment. View the markdown of the notebook for more details into the specific hardware requirement.

This pipeline will produce a total of 22 directories consisting of 6 samples each aligned to all 22 of the human chromasomes from the notebook: `RNA-seq-nf.ipynb`

The `view_summerized_experiment.R` then takes each .RDS file produced by the pipeline and ranks each gene of every sample by their count/expression level.


### Citation
`The nf-core framework for community-curated bioinformatics pipelines.
Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
Nat Biotechnol. 2020 Feb 13. doi: 10.1038/s41587-020-0439-x.`

### Understanding each output from pipeline with command:
'''

!nextflow run nf-core/rnaseq \
    --input /content/drive/MyDrive/Colab_Notebooks/RNA-seq-nf/single-end.csv \
    --outdir nf_temp/$folder_name \
    --gtf /content/drive/MyDrive/Colab_Notebooks/RNA-seq-nf/Homo_sapiens.GRCh38.115.gtf.gz \
    --fasta $fasta_file_path \
    -profile conda

'''

### Basic processes are:
Overview: https://nf-co.re/rnaseq/3.22.2/docs/output/

Visit this link to see all the possible outputs of this pipeline. For now I am interested in:
1. Alignment and quantification
    - STAR & Salmon
2. Alignment post-processing
    - SAMtools
    - picard MarkDuplicates
3. Quality Control 
    - RSeQC
    - QualiMap
    - DupRadar
    - featureCounts
    - DESeq2
    - MultiQC

There are many more features in this pipeline, but the primary objective of this project is to:
1. perform RNAseq analysis on reads from the: https://www.encodeproject.org/
2. validate process with quality control parameters
3. rank gene expression for each sample in each aligned chromasome