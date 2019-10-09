# scRNAseq_MMB
Single-Cell Capture, RNA-seq, and Transcriptome Analysis from the Neural Retina

This GitHub page provides scRNA-seq data and R script for the demonstration of scRNA-seq data analysis from human retinal organoids , which was published in Methods Molecular Biology, Vol. 2092, 2020.

Here, we describe a workflow of scRNA-seq data analysis using open-source tools. we provide an example procedure for analyzing scRNA-seq data starting from read count or the unique molecular identifier(UMI) profile matrix to decipher the biological significance, using the Seurat package (version 2.3.4; https://satijalab.org/seurat/), one of the most commonly used tools. 

In brief, the data processing step includes quality control (QC), classifying cell population, identifying cell type-specific genes, and pseudo-time analysis. Different workflow and tools can be used depending on the study design and biological questions. Briefly, after mapping the read sequences to the reference genome using alignment software, such STAR, HISAT, and Bowtie, the read count or UMI matrix with genes in row and cells in column can be generated [1–3]. Data QC is performed to filter low-quality cells and genes by the proportion of unexpressed genes, followed by normalization and scaling of the filtered matrix data. Next, the high-dimensional data are visualized using PCA, canonical correlation analysis (CCA), or t-Distributed Stochastic Neighbor Embedding (t-SNE, https://lvdmaaten.github.io/tsne/) to investigate cell subpopulation structures. The identity of each cluster is determined based on the expression pattern of highly variable genes or marker genes. Significantly expressed genes unique to each cell cluster can be identified based on comparing with the rest of the cluster. In addition, Gene Ontology (GO) or pathway analysis of the significant genes can be performed for differentially regulated signaling pathways in each cluster. 


# Data and R script
- scRNAseq_analysis_demo.R : R script for processing the demo read count data( scRNAseq_readCount_organoid8M.zip)
- scRNAseq_readCount_organoid8M.zip: Count matrix data
- utilities_sc.R: Functions for data processing, used by "source(utilities_sc.R)"





# References

1. Dobin A, Davis CA, Schlesinger F, Drenkow J,Zaleski C, Jha S, Batut P, Chaisson M, Gingeras TR (2013) STAR: ultrafast universal RNA-seq aligner. Bioinformatics 29(1):15–21. https://doi.org/10.1093/bioinformatics/bts635 
2. Kim D, Langmead B, Salzberg SL (2015) HISAT: a fast spliced aligner with low memory requirements. Nat Methods 12(4):357–360.
https://doi.org/10.1038/nmeth.3317
3. Langmead B, Trapnell C, Pop M, Salzberg SL (2009) Ultrafast and memory-efficient align937 ment of short DNA sequences to the human genome. Genome Biol 10(3):R25. https://doi.org/10.1186/gb-2009-10-3-r25
