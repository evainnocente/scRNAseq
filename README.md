# scRNAseq

## Introduction

Single-cell RNA sequencing is a useful method to investigate the transcriptome of each cell. Millions of cells can be studied at one time, allowing cell types to be accurately identified and advancing many areas of research including health and disease (Jovic et al., 2022). 

A typical scRNAseq workflow consists of isolating RNA from cells of interest, reverse transcribing it into cDNA, amplifying it, and sequencing it (Jovic et al., 2022). Once sequenced, computational methods are applied to evaluate differences in gene expression. Typical steps of this process include preprocessing and quality checks, gene filtering, normalisation, differential expression analysis, feature selection, dimensionality reduction, batch correction, cluster analysis, cluster annotation, and more (He, Lin, & Chen, 2022; Billato et al., 2025).  

Several tools and packages exist for these purposes (He, Lin, & Chen, 2022). A popular and useful R package is Seurat (He, Lin, & Chen, 2022; Butler et al., 2018), which notably exhibits good clustering performance (Germain, Sonrel, & Robinson; 2020). While other available tools like OSCA (Amezquita et al., 2020) and scrapper (Lun & Kancherla, 2023) have been shown to cluster data more accurately than Seurat (Billato et al., 2025), I will use Seurat for my analysis. I chose Seurat given its ease of use (comprehensive vignettes and tutorials) and ability to perform all necessary aspects of the analysis while only using one package. 
Kazer et al. (2024) studied the differences in transcription between different tissue types in the nasal mucosa of mice before and during primary and secondary infection by influenza A virus, in order to capture changes in the nasal mucosa over the course of infection. The three tissue types included the respiratory mucosa (RM), the olfactory mucosa (OM), and the lateral nasal gland (LNG). I will use data from this study (in Seurat object format) to partially recreate the work of Kazer et al. (2024) and evaluate differences in gene expression between the three tissue types. 


## Methods

I loaded the data as a Seurat object into R Studio v4.5.2. For reproducibility, I used as many of the same filtering parameters as Kazer et al. (2024) as possible. I retained only cells with unique feature counts over 500, percentage of mitochondrial DNA less than 15%, and cells with RNA counts between 750 and 100,000, cells with less than 10,000 hashes and removed genes not expressed in at least ten cells. However, I also chose to filter out cells with unique feature counts over 2500, as per the “Seurat Guided Clustering Tutorial” (https://satijalab.org/seurat/articles/pbmc3k_tutorial.html#cluster-the-cells). I normalised the data using the Log Normalise method, and the variance stabilising transformation with 5000 features. I then scaled the data using the scaleData function from Seurat. I used Seurat v5.4.0.

I ran a PCA to reduce dimensionality and used the elbow method to determine the optimal number of principal components. I used the optimal number of PCs to find clusters in the data using a resolution value of 0.5, as specified by the “Seurat Guided Clustering Tutorial” (https://satijalab.org/seurat/articles/pbmc3k_tutorial.html#cluster-the-cells). I annotated the clusters with the SingleR package v2.12.0 (Aran et al., 2019), using the fine-scale labels for more detail as per Lun (2025). I used FindAllMarkers() from Seurat to identify genes that are significantly differentially expressed in the top cluster (cluster 0, annoatted as Neurons) compared to all other clusters.

I pseudobulked the samples in order to assess differential expression of genes (as per Murphy & Skene, 2022) in the neurons cluster, which I did with DESeq2 v1.50.2 (Love, Hubers, & Anders, 2014). I looked at the expression of genes between time points: comparing gene expression at two, five, eight, and fourteen days post influenza infection compared to naïve state (0 days post infection; Kazer et al., 2025). I then performed Gene Set Enrichment Analysis with the clusterProfiler package v.4.18.4 (Yu et al., 2012), focusing on comparing naïve vs fourteen days post infection time points. 

## Results

The optimal number of PCs was 17 according to the elbow plot (Fig. 1). Clustering analysis found that there are 34 cell-type clusters in the data (Fig. 2).

<img width="2130" height="1467" alt="elbow" src="https://github.com/user-attachments/assets/85d1f0eb-3693-430c-8d26-6ece6881b174" />

Figure 1. Elbow plot used to determine optimal number of principal components.


Figure 2. UMAP of cell type clusters in the data.

Each time point of the experiment had a different composition of cell type clusters (Fig. 3A), as did each tissue type (Fig. 3B). A UMAP of the labelled clusters is shown in Fig. 4. The top cluster (cluster 0) was annotated as Neurons, so I focused on this cluster only for further analyses. The top 20 most significantly upregulated markers in this cluster compared to all other clusters are found in Table 1.

<img width="3000" height="2400" alt="clusters_timepoints" src="https://github.com/user-attachments/assets/e1d47b25-ca10-4053-a20e-82e40f4ba73d" />

Figure 3A. Composition of cell type clusters in the data at different experimental timepoints and B) different tissue types.

<img width="3000" height="2400" alt="clusters_UMAP" src="https://github.com/user-attachments/assets/17cb6229-ac4d-4e8f-85ab-c19939f7179c" />

Figure 4. UMAP of cell type clusters in the data, annotated by SingleR. After annotation, the data only separated into 14 distinct clusters instead of 34. 

Table 1. Top 20 significant upregulated markers in Neurons cluster, compared to all other clusters in the data.

| Gene          | p_val | avg_log2FC  | pct.1 | pct.2 | p_val_adj | cluster |
|---------------|-------|-------------|-------|-------|-----------|---------|
| Cnga4         | 0     | 7.501545559 | 0.729 | 0.007 | 0         | Neurons |
| Gm11992       | 0     | 7.323761053 | 0.682 | 0.007 | 0         | Neurons |
| Pth2          | 0     | 7.172649403 | 0.504 | 0.005 | 0         | Neurons |
| Ano2          | 0     | 7.121006938 | 0.656 | 0.007 | 0         | Neurons |
| Kcnc4         | 0     | 6.947227081 | 0.775 | 0.01  | 0         | Neurons |
| Umodl1        | 0     | 6.889571972 | 0.925 | 0.015 | 0         | Neurons |
| Nrn1l         | 0     | 6.846576326 | 0.885 | 0.013 | 0         | Neurons |
| Adcy3         | 0     | 6.815156114 | 0.94  | 0.022 | 0         | Neurons |
| Slc1a2        | 0     | 6.743434866 | 0.656 | 0.009 | 0         | Neurons |
| Faim2         | 0     | 6.742757849 | 0.863 | 0.013 | 0         | Neurons |
| Cnga2         | 0     | 6.709677889 | 0.879 | 0.014 | 0         | Neurons |
| Gm6878        | 0     | 6.704290364 | 0.273 | 0.003 | 0         | Neurons |
| Cngb1         | 0     | 6.649576135 | 0.703 | 0.01  | 0         | Neurons |
| A730046J19Rik | 0     | 6.533680827 | 0.574 | 0.009 | 0         | Neurons |
| Gldc          | 0     | 6.501306585 | 0.363 | 0.005 | 0         | Neurons |
| Ric8b         | 0     | 6.398302096 | 0.927 | 0.025 | 0         | Neurons |
| Dlg2          | 0     | 6.379007892 | 0.419 | 0.015 | 0         | Neurons |
| Atp8b3        | 0     | 6.35937516  | 0.259 | 0.003 | 0         | Neurons |
| Fam81b        | 0     | 6.341692151 | 0.771 | 0.012 | 0         | Neurons |
| S100a5        | 0     | 6.326618528 | 0.951 | 0.363 | 0         | Neurons |


Within the Neurons cluster, I calculated the differential expression between Naïve cells and all other time points, but I will only report the results of the comparison between Naïve and 14 days post-infection for the sake of clarity. The top 20 most differentially expressed significant genes between naïve and day 14 of infection are found in Table 2 and Fig. 5. Gene Set Enrichment Analysis detected several enriched GO terms, seen in Fig. 6. 

<img width="3000" height="2400" alt="DE_d14_vs_naive_volcanoplot" src="https://github.com/user-attachments/assets/c3865670-ac9e-498b-b388-add70802b0e0" />

Figure 5. Volcano plot of the differentially epxressed genes in the neuron cell type cluster, between 0 and 14 days post influenza infection.

Table 2. Top 20 significantly differentially expressed genes at 14 days post infetcion compared to 0 days post infection, in the neurons cell type cluster.  

| gene     | baseMean    | log2FoldChange | lfcSE       | pvalue               | padj                 |
|----------|-------------|----------------|-------------|----------------------|----------------------|
| mt-Atp8  | 105.6692536 | 2.129425375    | 0.243848985 | 6.56776608523429e-20 | 1.39433673989524e-16 |
| Fosb     | 7.942069444 | 1.739571416    | 0.4754735   | 7.6071511956996e-06  | 0.000259832          |
| Nme7     | 106.3530869 | 1.632024251    | 0.250927881 | 4.36121903077514e-12 | 1.54314466705594e-09 |
| Gm42418  | 3675.193617 | 1.611080655    | 0.483364177 | 3.369222103422e-05   | 0.000870642          |
| Brd1     | 35.04680657 | -1.586043381   | 0.116556887 | 6.28652536385363e-43 | 5.3385173389845e-39  |
| Plac8    | 5.643845691 | -1.53069498    | 0.584934853 | 0.000188875          | 0.003405368          |
| Pbxip1   | 4.639589648 | 1.485346741    | 0.397263034 | 1.21610310685091e-05 | 0.000376903          |
| Dlg2     | 106.8577606 | -1.48443532    | 0.240949333 | 4.19538992559544e-11 | 1.07961367418656e-08 |
| Eps8l1   | 4.207638555 | 1.408787253    | 0.264521206 | 5.35365813099171e-09 | 7.10363513255962e-07 |
| Ptprz1   | 9.191059344 | 1.349243471    | 0.426288308 | 4.69787342439799e-05 | 0.001108176          |
| Dio2     | 3.767114851 | -1.338203737   | 0.400056998 | 4.58857587110066e-05 | 0.00109149           |
| Trib3    | 9.030108082 | -1.336478      | 0.309371767 | 6.32209126827373e-07 | 3.70256545173659e-05 |
| Cast     | 4.588221237 | 1.333023699    | 0.426718397 | 7.44169078716634e-05 | 0.001569342          |
| Mbd2     | 47.04601168 | -1.287415524   | 0.128899988 | 1.36433633936844e-24 | 5.79297209695838e-21 |
| Csnk2a2  | 16.34294147 | -1.240864426   | 0.155759421 | 5.05243218160248e-17 | 5.36315676077103e-14 |
| Tsix     | 14.18492881 | 1.230461353    | 0.297177217 | 1.73299943919601e-06 | 8.13073549041575e-05 |
| Hist1h4d | 30.81437522 | 1.151968995    | 0.19060641  | 8.08472500988106e-11 | 1.96158527954028e-08 |
| Nmt1     | 17.1838069  | 1.141979929    | 0.318071277 | 1.69010592393376e-05 | 0.000501831          |
| Kirrel2  | 126.6786381 | -1.092980902   | 0.191710263 | 5.79115447425691e-10 | 1.11769281353158e-07 |
| Akap12   | 4.099115697 | 1.086789144    | 0.394813007 | 0.000310509          | 0.004928677          |


<img width="3000" height="2400" alt="gsea_d14_vs_naive_dotplot" src="https://github.com/user-attachments/assets/744b66d6-5927-435a-9f9e-b22475c1bda1" />

Figure 6. Dotplot of significantly enriched GO categories between 0 and 14 days post infection in the neurons cell type cluster. 

  
## Discussion

I detected 34 clusters, which after annotation, represented only 14 distinct clusters of cell types. This could represent a limitation of the analyses, and in the future I would annotate the clusters manually in order to capture as much distinction between cell types as possible. Given that Kazer et al. (2025) detected 42 cell type clusters, it seems that my analysis is unable to distinguish some cell types from others.

However, similar to Kazer et al. (2025), neural cells and fibroblasts were some of the largest clusters (Fig. 4). I also found that cell types were represented unevenly across tissue types (OM, RM, LNG) and time points (0, 2, 5, 8, 14 days post infection) (Fig. 3A & 3B). Neurons were identified as the top cluster (cluster 0). Within the neurons cluster, the top 20 significant genes with the largest positive log-fold change can be seen in Table 1. These genes are upregulated in neurons compared to all other cell types in the data. Some of the top genes include Pth2 (parathyroid hormone 2) which is expressed in the brain and involved in pituitary hormone release, anxiety, and nociception (RefSeq: NM_053256.2), as well as Cnga4 (cyclic nucleotide gated channel alpha 4) which is involved in smell (RefSeq: NM_001033317.3). That these genes are included in the neurons cluster supports that cell types in the data have been clustered correctly. In fact, Kazer et al. (2025) reported that olfactory sensory neurons, and neurons in general, represent a large fraction of their dataset, representative of their importance in the mouse nasal mucosa. 

I also assessed differential expression of genes and GO term enrichment in the neuron cluster across time points of the experiment, and in particular between 0 and 14 days post infection, as by this time, cellular composition of the tissue had changed from the naïve state (Kazer et al., 2025). The significantly differentially expressed gene with the largest absolute log-fold change between 0 and 14 days post infection was mt-Atp8 (LFC = 2.129425, adjusted p-value = 1.39433673989524 x 10-16) (Table 2). This is a mitochondrial gene that is involved in ATP synthase activity (Accession: GCF_000001635.27). A study of peripheral blood mononuclear cells found that mitochondrial genes are upregulated in response to Sars-COV2 infection (Shao et al., 2006), and mitochondria may play a role in antiviral immunity more broadly (Koshiba, Bashiruddin, & Kawabata, 2011). That mt-Atp8 is upregulated in neurons after influenza-A infection could hint at its involvement in viral immunity, and presents an interesting avenue for future research. Further support for the involvement of mitochondria in antiviral immunity is the significant enrichment of the GO categories “mitochondrial ATP synthesis”, “mitochondrial gene expression”, and “mitochondrial translation” (Fig. 6). 

Overall, this analysis provides an insight into the cell types that make up each tissue type of the mouse nasal mucosa, and provide a small window of insight into gene expression changes over the course of influenza infection. Further work could characterise could expand cell type cluster identification and resolve gene expressed profiles of different cell types over infection time in more detail, to provide further insight into mouse viral immunity.

## References

Amezquita, R. A., Lun, A. T. L., Becht, E., Carey, V. J., Carpp, L. N., Geistlinger, L., Marini, F., Rue-Albrecht, K., Risso, D., Soneson, C., Waldron, L., Pagès, H., Smith, M. L., Huber, W., Morgan, M., Gottardo, R., & Hicks, S. C. (2020). Orchestrating single-cell analysis with Bioconductor. Nature Methods, 17(2), 137–145. https://doi.org/10.1038/s41592-019-0654-x

Aran, D., Looney, A. P., Liu, L., Wu, E., Fong, V., Hsu, A., Chak, S., Naikawadi, R. P., Wolters, P. J., Abate, A. R., Butte, A. J., & Bhattacharya, M. (2019). Reference-based analysis of lung single-cell sequencing reveals a transitional profibrotic macrophage. Nature Immunology, 20(2), 163–172. https://doi.org/10.1038/s41590-018-0276-y

ATP synthase F0 subunit 8. (n.d.). NCBI. Retrieved April 14, 2026, from https://www.ncbi.nlm.nih.gov/datasets/gene/17706/

Billato, I., Pages, H., Carey, V., Waldron, L., Sales, G., Romualdi, C., & Risso, D. (2025). Benchmarking large-scale single-cell RNA-seq analysis. https://doi.org/10.1101/2025.10.28.681564

Butler, A., Hoffman, P., Smibert, P., Papalexi, E., & Satija, R. (2018). Integrating single-cell transcriptomic data across different conditions, technologies, and species. Nature Biotechnology, 36(5), 411–420. https://doi.org/10.1038/nbt.4096

Cyclic nucleotide gated channel alpha 4. (n.d.). NCBI. Retrieved April 14, 2026, from https://www.ncbi.nlm.nih.gov/datasets/gene/233649/

Germain, P.-L., Sonrel, A., & Robinson, M. D. (2020). Pipecomp, a general framework for the evaluation of computational pipelines, reveals performant single cell rna-seq preprocessing tools. Genome Biology, 21(1), 227. https://doi.org/10.1186/s13059-020-02136-7

He, J., Lin, L., & Chen, J. (2022). Practical bioinformatics pipelines for single-cell RNA-seq data analysis. Biophysics Reports, 8(3), 158–169. https://doi.org/10.52601/bpr.2022.210041

Jovic, D., Liang, X., Zeng, H., Lin, L., Xu, F., & Luo, Y. (2022). Single‐cell RNA sequencing technologies and applications: A brief overview. Clinical and Translational Medicine, 12(3), e694. https://doi.org/10.1002/ctm2.694

Kazer, S. W., Match, C. M., Langan, E. M., Messou, M.-A., LaSalle, T. J., O’Leary, E., Marbourg, J., Naughton, K., von Andrian, U. H., & Ordovas-Montanes, J. (2024). Primary nasal influenza infection rewires tissue-scale memory response dynamics. Immunity, 57(8), 1955-1974.e8. https://doi.org/10.1016/j.immuni.2024.06.005

Koshiba, T., Bashiruddin, N., & Kawabata, S. (2011). Mitochondria and antiviral innate immunity. International Journal of Biochemistry and Molecular Biology, 2(3), 257–262. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3193288/

Love, M. I., Huber, W., & Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology, 15(12), 550. https://doi.org/10.1186/s13059-014-0550-8

Lun, A. (2025). Assigning cell types with SingleR. https://bioconductor.org/books/release/SingleRBook/

Lun, A. T. L., & Kancherla, J. (2023). Powering single-cell analyses in the browser with WebAssembly. Journal of Open Source Software, 8(89), 5603. https://doi.org/10.21105/joss.05603

Murphy, A. E., & Skene, N. G. (2022). A balanced measure shows superior performance of pseudobulk methods in single-cell RNA-sequencing analysis. Nature Communications, 13(1), 7851. https://doi.org/10.1038/s41467-022-35519-4

Parathyroid hormone 2. (n.d.). NCBI. Retrieved April 14, 2026, from https://www.ncbi.nlm.nih.gov/datasets/gene/114640/

Shao, H., Lan, D., Duan, Z., Liu, Z., Min, J., Zhang, L., Huang, J., Su, J., Chen, S., & Xu, A. (2006). Upregulation of mitochondrial gene expression in pbmc from convalescent sars patients. Journal of Clinical Immunology, 26(6), 546–554. https://doi.org/10.1007/s10875-006-9046-y

Yu, G., Wang, L.-G., Han, Y., & He, Q.-Y. (2012). Clusterprofiler: An r package for comparing biological themes among gene clusters. OMICS : A Journal of Integrative Biology, 16(5), 284–287. https://doi.org/10.1089/omi.2011.0118


