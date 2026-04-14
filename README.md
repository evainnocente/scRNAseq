# scRNAseq

## Introduction

Single-cell RNA sequencing is a useful method to investigate the transcriptome of each cell. Millions of cells can be studied at one time, allowing cell types to be accurately identified and advancing many areas of research including health and disease (Jovic et al., 2022). 

A typical scRNAseq workflow consists of isolating RNA from cells of interest, reverse transcribing it into cDNA, amplifying it, and sequencing it (Jovic et al., 2022). Once sequenced, computational methods are applied to evaluate differences in gene expression. Typical steps of this process include preprocessing and quality checks, gene filtering, normalisation, differential expression analysis, feature selection, dimensionality reduction, batch correction, cluster analysis, cluster annotation, and more (He, Lin, & Chen, 2022; Billato et al., 2025).  

Several tools and packages exist for these purposes (He, Lin, & Chen, 2022). A popular and useful R package is Seurat (He, Lin, & Chen, 2022; Butler et al., 2018), which notably exhibits good clustering performance (Germain, Sonrel, & Robinson; 2020). While other available tools like OSCA (Amezquita et al., 2020) and scrapper (Lun & Kancherla, 2023) have been shown to cluster data more accurately than Seurat (Billato et al., 2025), I will use Seurat for my analysis. I chose Seurat given its ease of use (comprehensive vignettes and tutorials) and ability to perform all necessary aspects of the analysis while only using one package. 
Kazer et al. (2024) studied the differences in transcription between different tissue types in the nasal mucosa of mice before and during primary and secondary infection by influenza A virus, in order to capture changes in the nasal mucosa over the course of infection. The three tissue types included the respiratory mucosa (RM), the olfactory mucosa (OM), and the lateral nasal gland (LNG). I will use data from this study (in Seurat object format) to partially recreate the work of Kazer et al. (2024) and evaluate differences in gene expression between the three tissue types. 


## Methods

I loaded the data as a Seurat object into R Studio v4.5.2. For reproducibility, I used as many of the same filtering parameters as Kazer et al. (2024) as possible. I retained only cells with unique feature counts over 500, percentage of mitochondrial DNA less than 15%, and cells with RNA counts between 750 and 100,000, cells with less than 10,000 hashes and removed genes not expressed in at least ten cells. However, I also chose to filter out cells with unique feature counts over 2500, as per the “Seurat Guided Clustering Tutorial” (https://satijalab.org/seurat/articles/pbmc3k_tutorial.html#cluster-the-cells). I normalised the data using the Log Normalise method, and the variance stabilising transformation with 5000 features. I then scaled the data using the scaleData function from Seurat. I used Seurat v5.4.0.

I ran a PCA to reduce dimensionality and used the elbow method to determine the optimal number of principal components. I used the optimal number of PCs to find clusters in the data using a resolution value of 0.5, as specified by the “Seurat Guided Clustering Tutorial” (https://satijalab.org/seurat/articles/pbmc3k_tutorial.html#cluster-the-cells). I annotated the clusters with the SingleR package v2.12.0 (Aran et al., 2019), using the fine-scale labels for more detail as per Lun (2025). I focused on only one cell type cluster (neurons) for further analyses.

I pseudobulked the samples in order to assess differential expression of genes (as per Murphy & Skene, 2022) in the neurons cluster, which I did with DESeq2 v1.50.2 (Love, Hubers, & Anders, 2014). I looked at the expression of genes between time points: comparing gene expression at two, five, eight, and fourteen days post influenza infection compared to naïve state (0 days post infection; Kazer et al., 2025). I then performed Gene Set Enrichment Analysis with the clusterProfiler package v.4.18.4 (Yu et al., 2012), focusing on comparing naïve vs fourteen days post infection time points. 

## Results

The optimal number of PCs was 17 according to the elbow plot (Fig. 1). Clustering analysis found that there are 34 cell-type clusters in the data (Fig. 2).

<img width="2130" height="1467" alt="elbow" src="https://github.com/user-attachments/assets/85d1f0eb-3693-430c-8d26-6ece6881b174" />

Figure 1. Elbow plot used to determine optimal number of principal components.

<img width="3000" height="2400" alt="clusters_timepoints" src="https://github.com/user-attachments/assets/e1d47b25-ca10-4053-a20e-82e40f4ba73d" />


<img width="3000" height="2400" alt="clusters_UMAP" src="https://github.com/user-attachments/assets/17cb6229-ac4d-4e8f-85ab-c19939f7179c" />


<img width="3000" height="2400" alt="DE_d14_vs_naive_volcanoplot" src="https://github.com/user-attachments/assets/c3865670-ac9e-498b-b388-add70802b0e0" />

<img width="3000" height="2400" alt="gsea_d14_vs_naive_dotplot" src="https://github.com/user-attachments/assets/744b66d6-5927-435a-9f9e-b22475c1bda1" />

Each time point of the experiment had a different composition of cell type clusters (Fig. 3). A UMAP of the labelled clusters is shown in Fig. 4. One of the largest clusters was the cell type neurons, which is what I chose to focus on for further analysis.


While I calculated the differential expression between Naïve and all other time points, I will only report the results of the comparison between Naïve and 14 days post-infection for the sake of clarity. The top 20 most differentially expressed significant genes between naïve and day 14 of infection are found in Table 1. Gene Set Enrichment Analysis detected several enriched GO terms, seen in Fig. 5. 


-	Find top marker genes in neurons

  
## Discussion

I detected 34 clusters, which after annotation, represented only 14 distinct clusters of cell types. This could represent a limitation of the analyses, and in the future I would annotate the clusters manually in order to capture as much distinction between cell types as possible. Given that Kazer et al. (2025) detected 42 cell type clusters, it seems that my analysis is unable to distinguish some cell types from others.

However, similar to Kazer et al. (2025), neural cells and fibroblasts were some of the largest clusters (Fig. **). I also found that cell types were represented unevenly across tissue types (OM, RM, LNG) and time points (0, 2, 5, 8, 14 days post infection). (Fig. **). Neurons were identified as the top cluster (cluster 0). Within the neurons cluster, the top 20 significant genes with the largest positive log-fold change can be seen in Table **. These genes are upregulated in neurons compared to all other cell types in the data. Some of the top genes include Pth2 (parathyroid hormone 2) which is expressed in the brain and involved in pituitary hormone release, anxiety, and nociception (RefSeq: NM_053256.2), as well as Cnga4 (cyclic nucleotide gated channel alpha 4) which is involved in smell (RefSeq: NM_001033317.3). That these genes are included in the neurons cluster corroborates that the cell types in the data have been clustered correctly. In fact, Kazer et al. (2025) reported that olfactory sensory neurons, and neurons in general, represent a large fraction of their dataset, representative of their importance in the mouse nasal mucosa. 

I also assessed differential expression of genes and GO term enrichment in the neuron cluster across time points of the experiment, and in particular between 0 and 14 days post infection, as by this time, cellular composition of the tissue has changed from the naïve state (Kazer et al., 2025). The significantly differentially expressed gene with the largest absolute log-fold change between 0 and 14 days post infection was mt-Atp8 (LFC = 2.129425, adjusted p-value = 1.39433673989524 x 10-16). This is a mitochondrial gene that is involved in ATP synthase activity (source). A study of peripheral blood mononuclear cells found that mitochondrial genes are upregulated in response to Sars-COV2 infection (Shao et al., 2006), and mitochondria may play a role in antiviral immunity more broadly (Koshiba, Bashiruddin, & Kawabata, 2011). That mt-Atp8 is upregulated in neurons after influenza-A infection could hint at its involvement in viral immunity, and presents an interesting avenue for future research. Further support for the involvement of mitochondria in antiviral immunity is the significant enrichment of the GO categories “mitochondrial ATP synthesis”, “mitochondrial gene expression”, and “mitochondrial translation” (Fig. **). 

Overal...

