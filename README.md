# RM_PBMC_B_cell_scRNA_VDJseq
This repository contains codes and datasets for the publication: The B Cell Subtypes and VDJ Repertoire of Young Adult Rhesus Macaques (Journal of Immunology, 2025, in review).

In brief, we created a bioinformatic pipeline to describe both the gene expression and VDJ repertoire of single PBMC B cells in macaques (~13,000 B cells from six RMs), resulting in the annotation of 12 RM B cell clusters of distinct transcriptional profiles paired with their VDJ properties, such as IGH/K/L-V(D)JC gene family and subtype usage, CSR status, SHM rate and levels, CDRH3 AA length and hydrophobicity scores. Both gene expression and the VDJ repertoire of the global RM B cell population were also described. 

Here are all the original (non-downsampled) datasets for the generation of the main figures. 
1. scReoertoire_B_cell files were generated from filtered_contig_annotations.csv files from 6RMs via scRepertoire combineBCR function.
2. SHM_with_scPrepertoire_B_cell_object.RDS represented data combining from both the scRepertoire object and SHM information (derived from IMGT_StatClonotype and IGMT_HighVQuest. 
3. vgm_baseline_VDJ.RDS represented the VDJ information from the Platypus pipeline.
4. all_clone_pid files were created from combining both the platypus VDJ object, with the scRepertoire object and the SHM information. 
5. find_all_markers_latest_baseline_B_0.4.RDS contained all the differentially expressed genes between IGHD, IGHM, IGHG1 and IGHA, B cells for Figure 8, generated via FindMarkers function of Seurat

All the Seurat objects were downsampled to 5000 cells for easy reproduction, and they are accessible here [here](https://drive.google.com/drive/folders/13vwu4l_OMG-6dXU3cWPqLD0h93fYBzVe?usp=sharing)
Raw sequencing reads and other processed Seurat files (features.tsv, barcodes.tsv, matrix.mtx, and cell_metadata.txt.gz) are availble in NCBI GEO (accession number GSE302169).
