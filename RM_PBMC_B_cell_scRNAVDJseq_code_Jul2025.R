# RM_PBMC_B_cell_scRNA/VDJseq_code
# Author: Andy Chan 
# Date: 10 July 2025
# Description: This R document contained all the processed datasets and the codes employed in the production of all the main figures of the manuscript. All Seurat objects were downsampled to 5000 cells while all other datasets are in their original size. 

# Library -------------
library("seqinr")
library("writexl")
library("readxl")
library("dplyr")
library("Seurat") 
library("SeuratData")
library("patchwork")
library("purrr")
library("tidyr")
library("ggpubr")
library("spatstat.utils")
library("tibble")
library("cowplot")
library("ggplot2")
library("fields")
library("ROCR")
library("KernSmooth")
library("Matrix")
library("parallel")
library("clustree")
library("DoubletFinder")
library('magrittr')
library('pheatmap')
library("EnhancedVolcano")
library("metap")
library("stringr")
library("stringi")
library("SeuratWrappers")
library("gridExtra")
library("grid")
library("presto")
library("BPCells")
library('SingleCellExperiment')
library('escape')
library('dittoSeq')
library('DESeq2')
library("Azimuth")
library("scRepertoire")
library("gghighlight")
library("RColorBrewer")
library("seqinr")
library("gghighlight")
library("RColorBrewer")
library("ggpubr")
library("rstatix")
library("patchwork")
library("AnnotationDbi")
library("clusterProfiler")
library("AnnotationHub")
library("org.Mmu.eg.db")
library("org.Hs.eg.db")
library("alakazam")



# Load datasets ----------
Integrated_6RMs_PBMC_Seurat_downsample <- readRDS("Integrated_6RMs_PBMC_Seurat_downsample.RDS")
Integrated_6RMs_B_cells_Seurat_downsample <- readRDS("Integrated_6RMs_B_cells_Seurat_downsample.RDS")
vgm_baseline_VDJ <- readRDS("vgm_baseline_VDJ.RDS") #Nr_of_VDJ_chains <= 1 & Nr_of_VJ_chains <= 1
vgm_baseline_GEX_downsample <- readRDS("vgm_baseline_GEX_downsample.RDS")
combined_baseline_B <- readRDS("scRepertoire_B_cell_object.RDS")
combined_baseline_B_renamed <- readRDS("scRepertoire_B_cell_object_renamed.RDS")
SHM_results_combined_baseline_B <- readRDS("SHM_with_scRepertoire_B_cell_object.RDS")
all_clone_pid_baseline_IGH <- readRDS("all_clone_pid_baseline_IGH.RDS")
all_clone_pid_baseline_IGH_mutational_state <- readRDS("all_clone_pid_baseline_IGH_mutational_state.RDS")
find_all_markers_latest_baseline_B_0.4 <- readRDS ("find_all_markers_latest_baseline_B_0.4.RDS")

# Fig1 ------------
# Fig1A and B
Fig1A <- DimPlot(integrated_MeV_baseline, reduction = "umap.rpca", group.by = "predicted.celltype.l1", label = FALSE) +
  ggtitle ("Predicted PBMC \nSeurat clusters")+ 
  theme(axis.text  = element_text(size = 6))

Fig1B <- DimPlot(Integrated_6RMs_B_cells_Seurat_downsample, reduction = "umap", group.by = c("RNA_snn_res_B_0.4"), label = TRUE, label.size = 3) +
  scale_color_discrete(labels = c(
    "0 - M-MBC 1",
    "1 - CD11c+ MBC 1",
    "2 - Naive B cells 1",
    "3 - C-MBC 1",
    "4 - M-MZL BC",
    "5 - CD11c+ MBC 2",
    "6 - M-MBC 2",
    "7 - Activated B cells",
    "8 - C-MBC 2",
    "9 - Contamination",
    "10 - Naive B cells 2",
    "11 - Naive B cells 3",
    "12 - Antibody-secreting cells"
  ))+ 
  theme(
    legend.text  = element_text(size = 7), 
  )+
  ggtitle ("B cell \nSeurat clusters") + 
  tag_theme +
  labs(tag = "B")+
  theme(axis.text  = element_text(size = 6))

# Fig1C
vln_list <- VlnPlot(
  Integrated_6RMs_B_cells_Seurat_downsample,
  features = c("percent.mt", "percent.rb", "nFeature_RNA", "nCount_RNA"),
  pt.size  = 0,
  combine  = FALSE                    
)

titles <- c("Mitochondrial RNA (%)",
            "Ribosomal RNA (%)",
            "Genes Detected (nFeature_RNA)",
            "UMI Count (nCount_RNA)")

for (i in seq_along(vln_list)) {
  vln_list[[i]] <- vln_list[[i]] +
    ggtitle(titles[i])
}

vln_list <- lapply(vln_list, function(p)
  
  p + theme(
    text        = element_text(size = 7), 
    axis.text   = element_text(size = 7),
    axis.title.x = element_blank(),
    legend.position = "none",
    plot.title  = element_text(size = 7, face = "bold")
  )
)

Fig1C <- grid.arrange(grobs = vln_list, ncol = 4, bottom = textGrob("Clusters (celltypes)", gp = gpar(fontsize = 7)))

# Fig1D and E
cell_count_per_txgrp_timepoint_1 <- table(Integrated_6RMs_B_cells_Seurat_downsample$RNA_snn_res_B_0.4, Integrated_6RMs_B_cells_Seurat_downsample$RM_ID)
cell_count_per_txgrp_timepoint_1 <- as.data.frame(cell_count_per_txgrp_timepoint_1)
cell_count_per_txgrp_timepoint_1$Var1 <- as.character(cell_count_per_txgrp_timepoint_1$Var1)
cell_count_per_txgrp_timepoint_1$Var2 <- as.character(cell_count_per_txgrp_timepoint_1$Var2)

sample_id_name <- unique(cell_count_per_txgrp_timepoint_1$Var2)

datatable <- cell_count_per_txgrp_timepoint_1 %>%
  group_by(Var2)%>%
  mutate(proportions = Freq/sum(Freq)*100) 

datatable_summary <- datatable %>%
  group_by(Var1)%>%
  summarise(total_count = sum(Freq), average_proportions = mean(proportions)) 

cluster_levels <- 0:12

Fig1D <- ggplot(datatable, aes(x = factor(Var1, level= cluster_levels), y= Freq, fill= factor(Var1, level= cluster_levels))) +
  geom_boxplot()+
  geom_jitter(alpha = 0.3)+
  xlab("Clusters (celltypes)") +
  ylab("Cell count") +
  labs(fill = "B cell subtypes")+
  theme_bw(base_size = 8)+
  stat_summary(fun = "mean", geom = "point", shape = 21, size = 2, color = "red")+ 
  theme(legend.position = "none")

Fig1D_2 <-  ggplot(datatable, aes(x = factor(Var1, level= cluster_levels), y= proportions, fill= factor(Var1, level= cluster_levels))) +
  geom_boxplot()+
  geom_jitter(alpha = 0.3)+
  xlab("Clusters (celltypes)") +
  ylab("Proportion") +
  labs(fill = "B cell subtypes")+
  theme_bw(base_size = 8)+
  stat_summary(fun = "mean", geom = "point", shape = 21, size = 2, color = "red")

Fig1D+ Fig1D_2


ggsave("Fig1A.tiff",
       plot = Fig1A,
       dpi = 600, 
       height = 3,
       width  = 3, 
       units = "in")   

ggsave("Fig1B.tiff",
       plot = Fig1B,
       dpi = 600, 
       height = 3,
       width  = 3.9, 
       units = "in")  

ggsave("Fig1C.tiff",
       plot = Fig1C,
       dpi = 600, 
       height = 1.5,
       width  = 6.9, 
       units = "in")  

ggsave("Fig1D.tiff",
       plot = Fig1D+ Fig1D_2,
       dpi = 600, 
       height = 2,
       width  = 6.9, 
       units = "in")  

# Fig2 -----------

# Fig2A 
Fig2A <- FeaturePlot (Integrated_6RMs_B_cells_Seurat_downsample, features = c("CD19", "MS4A1","IGHM", "TCL1A","ITGAX", "JCHAIN"), reduction = "umap", ncol = 6)&
  theme(plot.title = element_text(size = 10))  &
  theme(text  = element_text(size = 8), 
        axis.text  = element_text(size = 6))

Fig2A_2 <- VlnPlot(Integrated_6RMs_B_cells_Seurat_downsample, features = c("CD19", "MS4A1","IGHM", "TCL1A","ITGAX", "JCHAIN"),group.by = "RNA_snn_res_B_0.4", assay = "RNA", pt.size = 0, ncol = 6) & 
  theme(plot.title = element_text(size=10))&
  theme(text  = element_text(size = 8), 
        axis.text  = element_text(size = 6))

# Fig2B
genes_to_plot <- c(
  
  "FCER2", "TCL1A", "IL4R", "SELL","CXCR4", "CD40", "PLAAT4", "MPP7",  "SNX29", "ZEB1", "LYST","BACH2","FOXO1", "EBF1",
  
  "NFATC1","NFKB1", "CCL4L1","CCL5", "CCR7", "IRF4", "MYC", "CD69", "PANX1", "BCL2", "EGR1", "EGR2", "EGR3", "NFKBIE", "NFKBID","BCL2A1","MARCKS", "DUSP2", "SEMA7A", "FCRL5", "RILPL2", "TAF4B",  "NR4A1",
  
  "CD72", "LILRB1", "NOTCH2", "CEMIP2", "TNS3", "TYROBP", "DLGAP1", "HDAC9",
  
  "FCRL2", "ENSMMUG00000056792", "CR1","CMTM7","JCHAIN", "MZB1", "A4GALT","FUT8",
  
  "FOS", "TOX","TESC", "CD86", "ZBTB32", "FOSB" , "VPREB3","CENPM", "SSBP3", "SSPN",
  
  "ITGAX", "ITGB1","LITAF","FGR", "TOX2", "ZEB2",
  
  "IGHE", "PRDM1", "ZBP1", "DERL3", "XBP1", "MKI67", "TNFRSF17")

alldata <- ScaleData(Integrated_6RMs_B_cells_Seurat_downsample, features = genes_to_plot, assay = "RNA")

alldata$new_identity <- "Unassigned"

alldata$new_identity[alldata$RNA_snn_res_B_0.4 == "2"] <- "Naive B cells 1 (Cluster 2)"
alldata$new_identity[alldata$RNA_snn_res_B_0.4 == "10"] <- "Naive B cells 2 (Cluster 10)"
alldata$new_identity[alldata$RNA_snn_res_B_0.4 == "11"] <- "Naive B cells 3 (Cluster 11)"
alldata$new_identity[alldata$RNA_snn_res_B_0.4 == "7"] <- "Activated B cells (Cluster 7)"
alldata$new_identity[alldata$RNA_snn_res_B_0.4 == "0"] <- "M-MBC 1 (Cluster 0)"
alldata$new_identity[alldata$RNA_snn_res_B_0.4 == "6"] <- "M-MBC 2 (Cluster 6)"
alldata$new_identity[alldata$RNA_snn_res_B_0.4 == "4"] <- "M-MZL BC (Cluster 4)"
alldata$new_identity[alldata$RNA_snn_res_B_0.4 == "3"] <- "C-MBC 1 (Cluster 3)"
alldata$new_identity[alldata$RNA_snn_res_B_0.4 == "8"] <- "C-MBC 2 (Cluster 8)"
alldata$new_identity[alldata$RNA_snn_res_B_0.4 == "1"] <- "CD11c+ MBC 1 (Cluster 1)"
alldata$new_identity[alldata$RNA_snn_res_B_0.4 == "5"] <- "CD11c+ MBC 2 (Cluster 5)"
alldata$new_identity[alldata$RNA_snn_res_B_0.4 == "12"] <- "Antibody-secreting cells (Cluster 12)"
alldata$new_identity[alldata$RNA_snn_res_B_0.4 == "9"] <- "Contamination (Cluster 9)"

alldata$new_identity <- factor (alldata$new_identity, levels = rev(c("Naive B cells 1 (Cluster 2)","Naive B cells 2 (Cluster 10)","Naive B cells 3 (Cluster 11)", "Activated B cells (Cluster 7)", "M-MBC 1 (Cluster 0)", "M-MBC 2 (Cluster 6)", "M-MZL BC (Cluster 4)", "C-MBC 1 (Cluster 3)", "C-MBC 2 (Cluster 8)", "CD11c+ MBC 1 (Cluster 1)", "CD11c+ MBC 2 (Cluster 5)", "Antibody-secreting cells (Cluster 12)", "Contamination (Cluster 9)")))

Idents(alldata) <- "new_identity"

Fig2B <- DotPlot(subset(alldata, downsample = 500), features = genes_to_plot, group.by = "new_identity", assay = "RNA", dot.scale = 4) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))+
  scale_colour_gradient2(low = "blue", mid = "white", high = "red")& 
  theme(text  = element_text(size = 8), 
        axis.text  = element_text(size = 8), 
        axis.text.x  = element_text(size = 7), 
        legend.position = "right")

# Fig2C
vgm_gene_subtype_function_VDJ_cgene_1 <- function (IG_chain_gene_segment) {
  dataset <- vgm_baseline_VDJ %>%
    filter(!!sym(IG_chain_gene_segment) != '') %>%
    filter(seurat_clusters != '') %>%
    mutate (class_switch_status = VDJ_cgene)%>%
    mutate(across(class_switch_status, ~str_replace (., '^IGHM$', 'unswitched'))) %>%
    mutate(across(class_switch_status, ~str_replace (., '^IGHD$', 'unswitched'))) %>%
    mutate(across(class_switch_status, ~str_replace (., '^IGHG1$', 'switched'))) %>%
    mutate(across(class_switch_status, ~str_replace (., '^IGHG2$', 'switched'))) %>%
    mutate(across(class_switch_status, ~str_replace (., '^IGHG3$', 'switched'))) %>%
    mutate(across(class_switch_status, ~str_replace (., '^IGHG4$', 'switched'))) %>%
    mutate(across(class_switch_status, ~str_replace (., '^IGHA$', 'switched'))) %>%
    mutate(across(class_switch_status, ~str_replace (., '^IGHE$', 'switched'))) %>%
    mutate(VDJ_cgene = factor (VDJ_cgene, levels = c("IGHM", "IGHD", "IGHG1", "IGHG2", "IGHG3","IGHG4", "IGHA", "IGHE")))%>%
    group_by(orig.ident, sample_id,RM_ID, seurat_clusters, !!sym(IG_chain_gene_segment), class_switch_status) %>%
    summarize(n = n())  %>%
    group_by(seurat_clusters, sample_id) %>%
    mutate(ttl = sum(n), Freq = n/ttl*100)

  return(dataset)
}

vgm_VDJ_cgene <- vgm_gene_subtype_function_VDJ_cgene_1("VDJ_cgene")

vgm_VDJ_cgene_plot_1 <- ggplot(vgm_VDJ_cgene, aes(x=seurat_clusters , y= n))+
  geom_bar(stat = "identity", aes(fill = VDJ_cgene))+
  theme_bw(base_size=10)+
  xlab("B cell subtypes")+
  ylab("Cell Counts")+
  labs(fill="IGHC")+
  theme(legend.position = "bottom")+
  theme(legend.text = element_text(size = 6))

vgm_VDJ_cgene_plot_2 <- ggplot(vgm_VDJ_cgene, aes(x=seurat_clusters , y= n))+
  geom_col(position = "fill", width = 0.5, aes(fill = VDJ_cgene))+
  theme_bw(base_size=10)+
  xlab("B cell subtypes")+
  ylab("Proportions")+
  labs(fill="IGHC") +
  theme(legend.position = "bottom")+
  theme(legend.text = element_text(size = 6))

vgm_VDJ_cgene_plot_3 <- ggplot(vgm_VDJ_cgene, aes(x=seurat_clusters , y= n))+
  geom_col(position = "fill", width = 0.5, aes(fill = class_switch_status))+
  theme_bw(base_size=10)+
  xlab("B cell subtypes")+
  ylab("Proportions")+
  labs(fill="IGHC") +
  theme(axis.text.x=element_text(angle=0))+
  theme(legend.position = "bottom")+
  theme(legend.text = element_text(size = 6))

# Fig2D
mutational_state_levels <- c("highly mutated", "moderately mutated", "lowly mutated", "unmutated" )

SHM_per_cluster_1 <- ggboxplot(all_clone_pid_baseline_IGH, x = "seurat_clusters", y = "shm", color = "seurat_clusters")+
  geom_jitter(position=position_jitter(0.15), cex=1.5, alpha=0)+
  ylab ("SHM rate")+
  xlab ("B cell subtypes")  +
  theme_bw(base_size=10)+
  theme(legend.position="none") 

SHM_per_cluster_2 <- ggplot(all_clone_pid_baseline_IGH_mutational_state, aes(x = seurat_clusters,  y = n, fill= factor(mutational_state, levels = mutational_state_levels)))+
  geom_col()+
  labs(fill="SHM levels")+
  scale_fill_manual(values = c( "#FC4E07", "#E7B800", "#99CC66", "#00AFBB"))+
  theme_bw(base_size=10)+
  ylab ("Count")+
  xlab ("B cell subtypes")+
  theme(legend.position="none") 

SHM_per_cluster_3 <- ggplot(all_clone_pid_baseline_IGH_mutational_state, aes(x = seurat_clusters,  y = Freq ))+
  geom_col(position = "fill", width = 0.5, aes(fill = factor(mutational_state, levels = mutational_state_levels)))+
  labs(fill="SHM levels")+
  scale_fill_manual(values = c( "#FC4E07", "#E7B800", "#99CC66", "#00AFBB"))+
  theme_bw(base_size=10)+
  ylab ("Proportions")+
  xlab ("B cell subtypes")+
  theme(legend.position="right")


ggsave("Fig2A.tiff",
       plot = Fig2A,
       dpi = 600, 
       height = 1.5,
       width  = 10, 
       units = "in") 

ggsave("Fig2A_2.tiff",
       plot = Fig2A_2,
       dpi = 600, 
       height = 1.5,
       width  = 10, 
       units = "in") 

ggsave("Fig2B.tiff",
       plot = Fig2B,
       dpi = 600, 
       height = 4.5,
       width  = 10, 
       units = "in") 

ggsave("Fig2C.tiff",
       plot = vgm_VDJ_cgene_plot_1 + vgm_VDJ_cgene_plot_2 + vgm_VDJ_cgene_plot_3,
       dpi = 600, 
       height = 3,
       width  = 9, 
       units = "in") 

ggsave("Fig2D.tiff",
       plot = SHM_per_cluster_1 + SHM_per_cluster_2 +SHM_per_cluster_3,
       dpi = 600, 
       height = 2,
       width  = 9, 
       units = "in") 

# Fig3 --------------------
# Fig3A, B, D, E, F
DefaultAssay(Integrated_6RMs_B_cells_Seurat_downsample) <- "RNA"
Idents(Integrated_6RMs_B_cells_Seurat_downsample) <- "RNA_snn_res_B_0.4"

FindMarkers_function <- function (cell_type_1, cell_type_2, avglog2fc_value= a) {
  
  baseline_B_deg <- FindMarkers(Integrated_6RMs_B_cells_Seurat_downsample, ident.1 = cell_type_2, ident.2 = cell_type_1, test.use="wilcox", min.pct = 0.25, logfc.threshold = 0)
  
  my_volcano_plot <- function(clus) {
    clus$levels <- "NOT SIGNIFICANT" # Default state is no change
    clus$levels[clus$avg_log2FC > avglog2fc_value & clus$p_val_adj < 1e-10] <- "UP" # Upregulated
    clus$levels[clus$avg_log2FC < -avglog2fc_value & clus$p_val_adj < 1e-10] <- "DOWN" # Downregulated
    # clus$levels[clus$p_val_adj > 1e-10] <- "NO" # Not significant
    
    clus$delabel <- NA
    clus$gene_symbol <- rownames(clus)
    
    # Label only significant differential expressed genes
    clus$delabel[clus$levels == "UP"| clus$levels == "DOWN"] <- clus$gene_symbol[clus$levels == "UP" | clus$levels == "DOWN"]
    
    # Count the number of UP and DOWN levels
    num_up <- sum(clus$levels == "UP")
    num_down <- sum(clus$levels == "DOWN")
    
    clus <- clus %>%
      mutate(p_val_adj = ifelse(p_val_adj < 1e-300, 1e-300, p_val_adj))
    
    # Create the plot
    p <- ggplot(data=clus, aes(x=avg_log2FC, y=-log10(p_val_adj), color=levels, label=delabel)) +
      geom_point(alpha=0.5) + 
      theme_minimal() +
      theme(axis.line = element_line(colour = "black"))+
      theme(plot.title = element_text(size = 9, face = "bold"), 
            axis.text=element_text(size=6),
            axis.title=element_text(size=9,face="bold"), 
            legend.title= element_text(size = 9, face = "bold"), 
            legend.text = element_text(size = 8))+
      theme(legend.position="top")+
      geom_text_repel(max.overlaps = Inf, size = 2.2) +
      scale_color_manual(values=c( "DOWN" = "blue", "UP"= "red","NOT SIGNIFICANT" = "grey"),
                         breaks = c("DOWN", "UP", "NOT SIGNIFICANT"),
                         labels = c(paste0("cluster ", cell_type_1), paste0("cluster ", cell_type_2),"other \ngenes"), 
                         name = "Highly expressed \nin") +
      geom_vline(xintercept=c(-avglog2fc_value, avglog2fc_value), col="black", linetype = "dashed") +
      geom_hline(yintercept=-log10(10e-10), col="black", linetype = "dashed") 
      # annotate("text", x=Inf, y=Inf, label=paste("# of genes in blue:", num_up, "\n# of genes in red:", num_down), 
      #          hjust=1.1, vjust=2, size=3, color="black")
    
    return(p)
  }
  
  plot <- list()
  
  plot[[1]] <- my_volcano_plot (baseline_B_deg)+ggtitle(paste0("B cell clusters ", cell_type_1, " vs ", cell_type_2))
  
  do.call(gridExtra::grid.arrange, c(plot, ncol=1))
  
}
GO_plot_function <- function (cluster_number) {
  DefaultAssay(Integrated_6RMs_B_cells_Seurat_downsample) <- "RNA"
  
  Idents(Integrated_6RMs_B_cells_Seurat_downsample) <- "RNA_snn_res_B_0.4"
  
  top_genes_lfc_output_subtype <- find_all_markers_latest_baseline_B_0.4 %>%
    group_by(cluster) %>%
    filter (avg_log2FC > 0.40) %>%
    filter (pct.1 > 0.25) %>%
    arrange(cluster, desc(avg_log2FC), p_val_adj) %>%
    dplyr :: slice(1:50) %>%
    mutate(category = row_number()) %>% 
    ungroup()%>% 
    filter (cluster == cluster_number)
  
  genes_to_test <- top_genes_lfc_output_subtype$gene
  
  GO_results <- enrichGO(gene= genes_to_test, OrgDb = "org.Mmu.eg.db", ont = "BP", keyType = "SYMBOL") #BP, MP or CC
  as.data.frame(GO_results)
  
  graph <- plot(barplot(GO_results, showCategory = 5)) +
    ggtitle(paste0("Cluster ", cluster_number, " GO plot")) +
    theme(axis.text.x=element_text(size = 8), 
          axis.text.y=element_text(size = 8), 
          text=element_text(size = 8), 
          axis.title.x = element_text(size = 8), 
          legend.position = "bottom")
  return (graph)
}

Fig3A <- GO_plot_function("7")

Fig3B <- FindMarkers_function ("2", "7", 3) # threshold 3
Fig3D <- GO_plot_function("0")
Fig3E <- FindMarkers_function ("0", "6", 1.5) # threshold 1.5
Fig3F <- GO_plot_function("4")

# Fig3C
Fig3C_func <- function (qc_variables) {
  VlnPlot(subset(Integrated_6RMs_B_cells_Seurat_downsample, RNA_snn_res_B_0.4 %in% c("0", "6")), features = qc_variables ,pt.size = 0) +
    theme(plot.title = element_text(size=8))+
    stat_compare_means(method = "wilcox.test", label="p.signif", hide.ns = FALSE, label.x = 1.5)+
    stat_summary(fun = "mean", geom = "point", shape = 21, size = 2, color = "red")+
    theme(text  = element_text(size = 7), 
          axis.text  = element_text(size = 7))
}

titles <- c("Mitochondrial RNA (%)",
            "Ribosomal RNA (%)",
            "Genes Detected",
            "UMI Count")

Fig3C_plot_list <- list(Fig3C_func("percent.mt"), Fig3C_func("percent.rb"), Fig3C_func("nFeature_RNA"), Fig3C_func("nCount_RNA"))

for (i in seq_along(Fig3C_plot_list)) {
  Fig3C_plot_list[[i]] <- Fig3C_plot_list[[i]] +
    ggtitle(titles[i])
}

Fig3C <- grid.arrange(grobs = Fig3C_plot_list, ncol =4) 

ggsave("Fig3A.tiff",
       plot = Fig3A,
       dpi = 600, 
       height = 2.5,
       width  = 3.2, 
       units = "in") 

ggsave("Fig3B.tiff",
       plot = Fig3B,
       dpi = 600, 
       height = 3.7,
       width  = 3.7, 
       units = "in") 

ggsave("Fig3C.tiff",
       plot = Fig3C,
       dpi = 600, 
       height = 1.5,
       width  = 6.9, 
       units = "in") 

ggsave("Fig3D.tiff",
       plot = Fig3D,
       dpi = 600, 
       height = 2,
       width  = 3.2, 
       units = "in") 

ggsave("Fig3E.tiff",
       plot = Fig3E,
       dpi = 600, 
       height = 4,
       width  = 3.7, 
       units = "in")

ggsave("Fig3F.tiff",
       plot = Fig3F,
       dpi = 600, 
       height = 2.2,
       width  = 3.2, 
       units = "in") 
# Fig4 -------------
# Fig4A

Fig4A <-quantContig(combined_baseline_B_renamed, cloneCall="gene+nt", scale = F)+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(), text = element_text(size = 10), legend.position = "none")  +
  
  quantContig(combined_baseline_B_renamed, cloneCall="gene+nt", scale = T)+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(), text = element_text(size = 10))

Fig4B <-abundanceContig(combined_baseline_B_renamed, cloneCall = "gene+nt", scale = F) + theme (text = element_text(size = 10))

Fig4C <-compareClonotypes(combined_baseline_B_renamed, 
                  numbers = 5,
                  cloneCall="aa", 
                  graph = "alluvial")+
  theme (legend.position = "none")+
  ylab ("Clonal proportion")+
  theme (text = element_text(size = 10))

Fig4D <- clonalDiversity(combined_baseline_B_renamed, cloneCall = "gene+nt", group.by = "sample", n.boots = 100) +
  theme (text = element_text(size = 10))
  

ggsave("Fig4A.tiff",
       plot = Fig4A,
       dpi = 600, 
       height = 2,
       width  = 5, 
       units = "in") 

ggsave("Fig4B.tiff",
       plot = Fig4B,
       dpi = 600, 
       height = 2,
       width  = 3.45, 
       units = "in") 

ggsave("Fig4C.tiff",
       plot = Fig4C,
       dpi = 600, 
       height = 2,
       width  = 3.45, 
       units = "in")

ggsave("Fig4D.tiff",
       plot = Fig4D,
       dpi = 600, 
       height = 2,
       width  = 6.9, 
       units = "in") 

# Fig5 -------------
vizGenes_function_baseline_B_1 <- function (IG_chain, gene_segment) {
  dataset <-  vizGenes(combined_baseline_B, gene = gene_segment, chain = IG_chain, plot = "bar", order = "variance", scale = TRUE, exportTable = TRUE) %>%
    separate(element.names, into = c("position", "RM_ID", "DPI"), sep = "-", remove=FALSE)%>%
    mutate(sample_id = element.names)%>%
    mutate(txgrp = ifelse (RM_ID %in% c("13F", "47G", "69G"), "LAMV", "WTMeV"))%>%
    mutate(sex = ifelse (RM_ID %in% c("13F", "43F"), "M", "F"))%>%
    mutate(RM_ID = factor(RM_ID, levels = c("13F", "47G", "69G", "43F", "83H", "84H")))%>%    
    select(-c("sd", "varcount", "element.names"))
  
  dataset_2 <- dataset %>%
    group_by (Var1)%>%
    summarize(meanvalue = mean (n))%>%
    arrange(desc(meanvalue))
  
  dataset$Var1 <- factor(dataset$Var1, levels = dataset_2$Var1)
  
  return(dataset)
  
}
vizGenes_function_baseline_B_2 <- function (dataset){
  dataset_2 <- dataset %>%
    group_by (Var1)%>%
    filter (any(n > 0.025)) %>%
    ungroup () %>%
    mutate(RM_ID = factor(RM_ID, levels= c("13F", "47G", "69G", "43F", "83H", "84H")))
  
  return (dataset_2)
  
}
vizGenes_function_baseline_B_3 <- function (dataset){
  dataset_3 <-  dataset%>%
    mutate (substring_gene = substr(Var1,1,6),
            substring_gene = gsub ("S", "", substring_gene), 
            substring_gene = gsub("-.*", "", substring_gene))%>%
    group_by(sample_id,RM_ID,substring_gene, DPI) %>%
    summarise(Freq = sum(n))
  
}

vizGenes_IGHV_baseline_B_original <- vizGenes_function_baseline_B_1("IGH", "V")
vizGenes_IGHV_baseline_B_selected <- vizGenes_function_baseline_B_2(vizGenes_IGHV_baseline_B_original)
vizGenes_IGHV_baseline_B_substring_gene <- vizGenes_function_baseline_B_3(vizGenes_IGHV_baseline_B_original) 

vizGenes_IGHD_baseline_B_original <- vizGenes_function_baseline_B_1("IGH", "D")
vizGenes_IGHD_baseline_B_selected <- vizGenes_function_baseline_B_2(vizGenes_IGHD_baseline_B_original)
vizGenes_IGHD_baseline_B_substring_gene <- vizGenes_function_baseline_B_3(vizGenes_IGHD_baseline_B_original) 

vizGenes_IGHJ_baseline_B_original <- vizGenes_function_baseline_B_1("IGH", "J")
vizGenes_IGHJ_baseline_B_selected <- vizGenes_function_baseline_B_2(vizGenes_IGHJ_baseline_B_original)
vizGenes_IGHJ_baseline_B_substring_gene <- vizGenes_function_baseline_B_3(vizGenes_IGHJ_baseline_B_original) 

vizGenes_IGHC_baseline_B_original <- vizGenes_function_baseline_B_1("IGH", "C")
vizGenes_IGHC_baseline_B_selected <- vizGenes_function_baseline_B_2(vizGenes_IGHC_baseline_B_original)
vizGenes_IGHC_baseline_B_substring_gene <- vizGenes_function_baseline_B_3(vizGenes_IGHC_baseline_B_original) 

vizGenes_IGKV_baseline_B_original <- vizGenes_function_baseline_B_1("IGL", "V")%>%
  filter(substr(Var1, 1,3 ) == "IGK")
vizGenes_IGKV_baseline_B_selected <- vizGenes_function_baseline_B_2(vizGenes_IGKV_baseline_B_original)
vizGenes_IGKV_baseline_B_substring_gene <- vizGenes_function_baseline_B_3(vizGenes_IGKV_baseline_B_original) 

vizGenes_IGKJ_baseline_B_original <- vizGenes_function_baseline_B_1("IGL", "J")%>%
  filter(substr(Var1, 1,3 ) == "IGK")
vizGenes_IGKJ_baseline_B_selected <- vizGenes_function_baseline_B_2(vizGenes_IGKJ_baseline_B_original)
vizGenes_IGKJ_baseline_B_substring_gene <- vizGenes_function_baseline_B_3(vizGenes_IGKJ_baseline_B_original) 

vizGenes_IGKC_baseline_B_original <- vizGenes_function_baseline_B_1("IGL", "C")%>%
  filter(substr(Var1, 1,3 ) == "IGK")
vizGenes_IGKC_baseline_B_selected <- vizGenes_function_baseline_B_2(vizGenes_IGKC_baseline_B_original)
vizGenes_IGKC_baseline_B_substring_gene <- vizGenes_function_baseline_B_3(vizGenes_IGKC_baseline_B_original) 

vizGenes_IGLV_baseline_B_original <- vizGenes_function_baseline_B_1("IGL", "V")%>%
  filter(substr(Var1, 1,3 ) == "IGL")
vizGenes_IGLV_baseline_B_selected <- vizGenes_function_baseline_B_2(vizGenes_IGLV_baseline_B_original)
vizGenes_IGLV_baseline_B_substring_gene <- vizGenes_function_baseline_B_3(vizGenes_IGLV_baseline_B_original) %>%
  mutate (substring_gene = factor (substring_gene, levels = c("IGLV1", "IGLV2", "IGLV3", "IGLV4", "IGLV5", "IGLV6", "IGLV7", "IGLV8", "IGLV10", "IGLV11")))

vizGenes_IGLJ_baseline_B_original <- vizGenes_function_baseline_B_1("IGL", "J")%>%
  filter(substr(Var1, 1,3 ) == "IGL")
vizGenes_IGLJ_baseline_B_selected <- vizGenes_function_baseline_B_2(vizGenes_IGLJ_baseline_B_original)
vizGenes_IGLJ_baseline_B_substring_gene <- vizGenes_function_baseline_B_3(vizGenes_IGLJ_baseline_B_original) 

vizGenes_IGLC_baseline_B_original <- vizGenes_function_baseline_B_1("IGL", "C")%>%
  filter(substr(Var1, 1,3 ) == "IGL")
vizGenes_IGLC_baseline_B_selected <- vizGenes_function_baseline_B_2(vizGenes_IGLC_baseline_B_original)
vizGenes_IGLC_baseline_B_substring_gene <- vizGenes_function_baseline_B_3(vizGenes_IGLC_baseline_B_original) 

plot_1_gene_prop_per_gene_function_B <- function (dataset, IG_chain_gene_segment){
  
  plot <-  ggplot(dataset) + aes(x= Var1, y= n*100)+
    geom_boxplot()+
    geom_point (aes(group = RM_ID), size = 1)+
    theme_bw(base_size=10)+
    xlab(IG_chain_gene_segment)+
    ylab("Proportion (%)") +
    theme(axis.text.x=element_text(angle=90)) +
    stat_summary(fun = "mean", geom = "point", shape = 20, size = 2, color = "red")
  
  
  return (plot)
  
}
plot_2_substring_gene_function_B <- function (dataset, IG_chain_gene_segment){
  
  plot <- ggplot(dataset, aes (x = substring_gene, y = Freq*100))+
    geom_boxplot()+
    geom_point (aes(group = RM_ID), size = 1)+
    xlab(IG_chain_gene_segment)+
    ylab("Proportion (%)")+
    labs(fill="RM_ID")+
    theme_bw(base_size=10)+    
    theme(axis.text.x=element_text(angle=90))+
    stat_summary(fun = "mean", geom = "point", shape = 20, size = 3, color = "red")
  
  return (plot)
}

Fig5A_1 <- plot_2_substring_gene_function_B(vizGenes_IGHV_baseline_B_substring_gene, "IGHV")+
  plot_2_substring_gene_function_B(vizGenes_IGHD_baseline_B_substring_gene, "IGHD")

Fig5A_2 <- plot_1_gene_prop_per_gene_function_B(vizGenes_IGHJ_baseline_B_original, "IGHJ")+ 
  plot_1_gene_prop_per_gene_function_B(vizGenes_IGHC_baseline_B_original, "IGHC")

Fig5A_3 <- plot_2_substring_gene_function_B(vizGenes_IGKV_baseline_B_substring_gene, "IGKV")+
  plot_1_gene_prop_per_gene_function_B(vizGenes_IGKJ_baseline_B_original, "IGKJ")+
  plot_1_gene_prop_per_gene_function_B(vizGenes_IGKC_baseline_B_original, "IGKC")

Fig5A_4 <-plot_2_substring_gene_function_B(vizGenes_IGLV_baseline_B_substring_gene, "IGLV")+
  plot_1_gene_prop_per_gene_function_B(vizGenes_IGLJ_baseline_B_original, "IGLJ")+
  plot_1_gene_prop_per_gene_function_B(vizGenes_IGLC_baseline_B_original, "IGLC")

Fig5B_1 <- plot_1_gene_prop_per_gene_function_B(vizGenes_IGHV_baseline_B_selected, "IGHV")+
  plot_1_gene_prop_per_gene_function_B(vizGenes_IGHD_baseline_B_selected, "IGHD")

Fig5B_2 <- plot_1_gene_prop_per_gene_function_B(vizGenes_IGKV_baseline_B_selected, "IGKV")+
  plot_1_gene_prop_per_gene_function_B(vizGenes_IGLV_baseline_B_selected, "IGLV")

ggsave("Fig5A_1.tiff",
       plot = Fig5A_1,
       dpi = 600, 
       height = 2,
       width  = 3.45, 
       units = "in") 

ggsave("Fig5A_2.tiff",
       plot = Fig5A_2,
       dpi = 600, 
       height = 2,
       width  = 3.45, 
       units = "in") 

ggsave("Fig5A_3.tiff",
       plot = Fig5A_3,
       dpi = 600, 
       height = 2,
       width  = 6.9, 
       units = "in")

ggsave("Fig5A_4.tiff",
       plot = Fig5A_4,
       dpi = 600, 
       height = 2,
       width  = 6.9, 
       units = "in") 



ggsave("Fig5B_1.tiff",
       plot = Fig5B_1,
       dpi = 600, 
       height = 2,
       width  = 6.9, 
       units = "in") 

ggsave("Fig5B_2.tiff",
       plot = Fig5B_2,
       dpi = 600, 
       height = 2,
       width  = 6.9, 
       units = "in") 


# Fig6 --------------
vizGenes_function_baseline_B_4a <- function (IG_chain, gene_segment_1, gene_segment_2) {
  dataset <-  vizGenes(combined_baseline_B, gene = gene_segment_1, y.axis = gene_segment_2, chain = IG_chain, plot = "bar", order = "variance", scale = TRUE, exportTable = TRUE) %>%
    separate(element.names, into = c("position", "RM_ID", "DPI"), sep = "-", remove=FALSE)%>%
    mutate(sample_id = element.names)%>%
    mutate(txgrp = ifelse (RM_ID %in% c("13F", "47G", "69G"), "LAMV", "WTMeV"))%>%
    mutate(sex = ifelse (RM_ID %in% c("13F", "43F"), "M", "F"))%>%
    mutate(RM_ID = factor(RM_ID, levels = c("13F", "47G", "69G", "43F", "83H", "84H")))%>%    
    select(-c("sd", "varcount", "element.names"))%>%
    mutate (IGHV_IGHJ = paste0(Var1,"_", Var2)) %>%
    group_by(sample_id, RM_ID, IGHV_IGHJ) %>%
    summarise(Freq = sum(n))%>%
    group_by (IGHV_IGHJ)%>%
    filter (any(Freq > 0.01))
  
  dataset_2 <- dataset %>%
    group_by (IGHV_IGHJ)%>%
    summarize(meanvalue = mean (Freq))%>%
    arrange(desc(meanvalue))
  
  dataset$IGHV_IGHJ <- factor(dataset$IGHV_IGHJ, levels = dataset_2$IGHV_IGHJ)
  
  return(dataset)
  
}
vizGenes_function_baseline_B_4b <- function (IG_chain, gene_segment_1, gene_segment_2) {
  dataset <-  vizGenes(combined_baseline_B, gene = gene_segment_1, y.axis = gene_segment_2, chain = IG_chain, plot = "bar", order = "variance", scale = TRUE, exportTable = TRUE) %>%
    separate(element.names, into = c("position", "RM_ID", "DPI"), sep = "-", remove=FALSE)%>%
    mutate(sample_id = element.names)%>%
    mutate(txgrp = ifelse (RM_ID %in% c("13F", "47G", "69G"), "LAMV", "WTMeV"))%>%
    mutate(sex = ifelse (RM_ID %in% c("13F", "43F"), "M", "F"))%>%
    mutate(RM_ID = factor(RM_ID, levels = c("13F", "47G", "69G", "43F", "83H", "84H")))%>%    
    select(-c("sd", "varcount", "element.names"))%>%
    mutate (substring_gene = substr(Var1,1,6),
            substring_gene = gsub ("S", "", substring_gene), 
            substring_gene = gsub("-.*", "", substring_gene))%>%
    group_by(sample_id, RM_ID, Var2, substring_gene) %>%
    summarise(Freq = sum(n)) %>%
    ungroup()%>%
    group_by (RM_ID, Var2)%>%
    mutate (total = sum(Freq), Freq_per_IGHC = Freq/total)
  
  return(dataset)
  
}

plot_1b_gene_prop_per_gene_function_B <- function (dataset, IG_chain_gene_segment){
  
  plot <-  ggplot(dataset) + aes(x= IGHV_IGHJ, y= Freq*100)+
    geom_boxplot()+
    geom_point (aes(group = RM_ID))+
    theme_bw(base_size=10)+
    xlab("IGHV-IGHJ gene combination")+
    ylab("Proportion (%)") +
    theme(axis.text.x=element_text(angle=90)) +
    stat_summary(fun = "mean", geom = "point", shape = 20, size = 3, color = "red")
  
  
  return (plot)
  
}
plot_2b_substring_gene_function_B <- function (dataset, IG_chain_gene_segment){
  
  plot <- ggplot(subset(dataset, Var2 %in% c("IGHD", "IGHM", "IGHG1", "IGHA" )), aes (x = factor(Var2, levels = c("IGHD", "IGHM", "IGHG1", "IGHA" )), y = Freq_per_IGHC*100))+
    geom_boxplot(aes(fill = factor(Var2, levels = c("IGHD", "IGHM", "IGHG1", "IGHA" ))))+
    xlab ("IGHC")+
    ylab("Proportion (%)")+
    labs(fill="Selected IGHC")+
    theme_bw(base_size=10)+    
    theme(axis.text.x=element_text(angle=90), 
          legend.position = "bottom")+
    stat_compare_means(method = "wilcox.test", label="p.signif", hide.ns = TRUE, comparisons = list(c('IGHD', 'IGHM'),c('IGHD', 'IGHG1'), c('IGHD', 'IGHA'), c('IGHM', 'IGHG1'), c('IGHM', 'IGHA'), c('IGHG1', 'IGHA')), color="black")+
    facet_grid (~substring_gene)
  
  return (plot)
}

vizGenes_IGHV_IGHJ_baseline_B_original <- vizGenes_function_baseline_B_4a("IGH", "V", "J")
vizGenes_IGHV_IGHC_baseline_B_original <- vizGenes_function_baseline_B_4b("IGH", "V", "C")

Fig6A <- plot_1b_gene_prop_per_gene_function_B(vizGenes_IGHV_IGHJ_baseline_B_original, "IGHV_IGHJ")
Fig6B <- plot_2b_substring_gene_function_B (vizGenes_IGHV_IGHC_baseline_B_original, "IGHC")


ggsave("Fig6A.tiff",
       plot = Fig6A,
       dpi = 600, 
       height = 3.5,
       width  = 6.9, 
       units = "in") 

ggsave("Fig6B.tiff",
       plot = Fig6B,
       dpi = 600, 
       height = 3,
       width  = 6.9, 
       units = "in") 

# Fig7 ----------

# CDR3 AA length and hydrophobicity score 
all_clone_pid_baseline_IGH$hydrophobicity_score <- gravy(all_clone_pid_baseline_IGH$VDJ_cdr3s_aa)

lengthContig_output <- lengthContig(combined_baseline_B, cloneCall="aa", chain = "IGH", exportTable = TRUE)
lengthContig_output_IGH_function <- function (dataset) {
  
  lengthContig_output_IGH <- data.frame()
  
  for (i in 1:length(dataset)) {
    aaa <- combined_baseline_B[[i]] %>%
      select (c("IGH", "cdr3_aa1", "sample", "ID"))%>%
      mutate (sample_id = sample)
    lengthContig_output_IGH <- rbind(lengthContig_output_IGH, aaa)
  }
  
  lengthContig_output_IGH_2 <- lengthContig_output_IGH %>%
    mutate(length = nchar(lengthContig_output_IGH$cdr3_aa1)) %>%
    filter (complete.cases(.[,6 ])) %>%
    separate(IGH, into = c("IGHV", "IGHD", "IGHJ", "IGHC"), sep = "\\.") %>%
    mutate (IGHV_new = substring(IGHV, 1,5)) %>%
    separate(sample, into = c("position", "RM_ID", "DPI_orig"), sep = "-", remove=FALSE)%>%
    filter(IGHC != "NA")
  
  lengthContig_output_IGH_2_order_IGHV <- lengthContig_output_IGH_2 %>%
    group_by (IGHV)%>%
    summarize(meanvalue_IGHV = mean (length))%>%
    arrange(desc(meanvalue_IGHV))
  
  lengthContig_output_IGH_2_order_IGHV_family <- lengthContig_output_IGH_2 %>%
    group_by (IGHV_new)%>%
    summarize(meanvalue_IGHV_family = mean (length))%>%
    arrange(desc(meanvalue_IGHV_family))
  
  lengthContig_output_IGH_2_order_IGHC <- lengthContig_output_IGH_2 %>%
    group_by (IGHC)%>%
    summarize(meanvalue_IGHC = mean (length))%>%
    arrange(desc(meanvalue_IGHC))
  
  lengthContig_output_IGH_2 <- lengthContig_output_IGH_2 %>%
    mutate (IGHV = factor (IGHV, levels = lengthContig_output_IGH_2_order_IGHV$IGHV)) %>%
    mutate (IGHV_new = factor (IGHV_new, levels = lengthContig_output_IGH_2_order_IGHV_family$IGHV_new)) %>%
    mutate (IGHC = factor (IGHC, levels = lengthContig_output_IGH_2_order_IGHC$IGHC))
  
  return (lengthContig_output_IGH_2)
}
lengthContig_output_IGH_2 <- lengthContig_output_IGH_function(combined_baseline_B)
lengthContig_output_IGL <- data.frame()
for (i in 1:length(combined_baseline_B)) {
  aaa <- combined_baseline_B[[i]] %>%
    select (c("IGLC", "cdr3_aa2", "sample", "ID")) %>%
    mutate (sample_id = sample)
  lengthContig_output_IGL <- rbind(lengthContig_output_IGL, aaa)
}
lengthContig_output_IGL<- lengthContig_output_IGL %>%
  mutate (IGlight = substring(lengthContig_output_IGL$IGLC, 1,3)) %>%
  separate(sample, into = c("position", "RM_ID", "DPI_orig"), sep = "-", remove=FALSE)
lengthContig_output_IGL$length <- nchar(lengthContig_output_IGL$cdr3_aa2)
lengthContig_output_IgLambda <- lengthContig_output_IGL %>%
  filter (IGlight == "IGL")
lengthContig_output_IgKappa<- lengthContig_output_IGL %>%
  filter (IGlight == "IGK")

lengthContig_output_IGH_3_function <- function (dataset){
  
  dataset$hydrophobicity_score <- gravy(dataset$cdr3_aa1)
  
  lengthContig_output_IGHV_order <- dataset %>%
    group_by (IGHV)%>%
    summarize(meanvalue = mean (hydrophobicity_score))%>%
    arrange(desc(meanvalue))
  
  lengthContig_output_IGHV_substring_order <- dataset %>%
    group_by (IGHV_new)%>%
    summarize(meanvalue = mean (hydrophobicity_score))%>%
    arrange(desc(meanvalue))
  
  lengthContig_output_IGHC_order <- dataset %>%
    group_by (IGHC)%>%
    summarize(meanvalue = mean (hydrophobicity_score))%>%
    arrange(desc(meanvalue))
  
  dataset$IGHV <- factor(dataset$IGHV, levels = lengthContig_output_IGHV_order$IGHV)
  
  dataset$IGHV_new <- factor(dataset$IGHV_new, levels = lengthContig_output_IGHV_substring_order$IGHV_new)
  
  dataset$IGHC <- factor(dataset$IGHC, levels = lengthContig_output_IGHC_order$IGHC)
  
  return (dataset)
}

lengthContig_output_IGH_3 <- lengthContig_output_IGH_3_function(lengthContig_output_IGH_2)





# Fig7A-C
Fig7A <- ggboxplot(SHM_results_combined_baseline_B_1, x = "chain", y = "shm", color = "chain", size= 0.5)+
  ggtitle("SHM per chain")+
  ylab ("SHM rate (%)")+
  theme_bw(base_size=10)+
  theme (legend.position = "none")

plot_1_CDR3_length_of_all_samples <- function (dataset, IG_chain){
  
  plot  <- ggboxplot(dataset, x="RM_ID", y="length", fill = "RM_ID", size=0.5, add = c("boxplot", "mean_sd")) +
    geom_jitter(position=position_jitter(0.15), cex=1.5, alpha=0.0001) +
    ggtitle(paste0("CDR", IG_chain, "3 aa")) +
    theme_bw(base_size=10)+
    labs(x ="", y = "Length")+
    scale_x_discrete (labels = c("13F" = "RM_1", 
                                 "47G" = "RM_2", 
                                 "69G" = "RM_3", 
                                 "43F" = "RM_4", 
                                 "83H"= "RM_5", 
                                 "84H" = "RM_6"))+
    theme (legend.position = "none")+
    theme(axis.text.x=element_text(angle=90))
  
  
  return (plot) 
  
}

Fig7B_1 <- plot_1_CDR3_length_of_all_samples(lengthContig_output_IGH_2, "H")
Fig7B_2 <- plot_1_CDR3_length_of_all_samples(lengthContig_output_IgKappa, "K")+ ylim (0,18)
Fig7B_3 <- plot_1_CDR3_length_of_all_samples(lengthContig_output_IgLambda, "L")

Fig7C <- ggboxplot(lengthContig_output_IGH_3, x="RM_ID", y="hydrophobicity_score", size=0.5, add = c("boxplot", "mean_sd"), fill = "RM_ID") +
  ggtitle("CDRH3 aa") +
  theme_bw(base_size=10)+
  labs(y = "Hydrophobicity score")+
  scale_x_discrete (labels = c("13F" = "RM_1", 
                               "47G" = "RM_2", 
                               "69G" = "RM_3", 
                               "43F" = "RM_4", 
                               "83H"= "RM_5", 
                               "84H" = "RM_6"))+
  theme (legend.position = "none")+
  theme(axis.text.x=element_text(angle=90), 
        axis.title.x = element_blank())


ggsave("Fig7A.tiff",
       plot = Fig7A,
       dpi = 600, 
       height = 2,
       width  = 2, 
       units = "in") 

ggsave("Fig7B_1.tiff",
       plot = Fig7B_1,
       dpi = 600, 
       height = 2,
       width  = 1.6, 
       units = "in") 

ggsave("Fig7B_2.tiff",
       plot = Fig7B_2,
       dpi = 600, 
       height = 2,
       width  = 1.6, 
       units = "in") 

ggsave("Fig7B_3.tiff",
       plot = Fig7B_3,
       dpi = 600, 
       height = 2,
       width  = 1.6, 
       units = "in") 

ggsave("Fig7C.tiff",
       plot = Fig7C,
       dpi = 600, 
       height = 2,
       width  = 2, 
       units = "in") 

# Fig7D-G
Fig7D <- ggboxplot(SHM_results_combined_baseline_B_IGH, x = "substring_gene", y = "shm", color = "substring_gene")+
  theme_bw(base_size=10)+
  xlab("IGHV")+
  ylab ("SHM rate (%)")+
  theme(axis.text.x=element_text(angle=90))+
  theme(legend.position="none")

Fig7E <- ggplot(SHM_results_combined_baseline_B_IGH_mutational_state_v_gene_family, aes (x = substring_gene, y = Freq*100, fill= factor(mutational_state, levels = mutational_state_levels)))+
  geom_col()+
  theme_bw(base_size = 10)+
  labs(fill="SHM levels")+
  xlab("IGHV")+
  ylab ("Proportion (%)")+
  scale_fill_manual(values = c( "#FC4E07", "#E7B800", "#99CC66", "#00AFBB"))+
  theme(legend.position="left")+
  theme(axis.text.x=element_text(angle=90))

Fig7F <- ggboxplot(SHM_results_combined_baseline_B_IGH, x = "c_gene", y = "shm", color = "c_gene")+
  theme_bw(base_size = 10)+
  theme(axis.text.x=element_text(angle=90))+
  xlab("IGHC")+
  ylab ("SHM rate (%)")+
  theme(legend.position="none")


Fig7G <- ggplot(SHM_results_combined_baseline_B_IGH_mutational_state_c_gene, aes (x = c_gene, y = Freq*100, fill= factor(mutational_state, levels = mutational_state_levels)))+
  geom_col()+
  labs(fill="SHM levels")+
  theme_bw(base_size = 10)+
  theme(axis.text.x=element_text(angle=90))+
  xlab("IGHC")+
  ylab ("Proportion (%)")+
  theme(legend.position="none")+
  scale_fill_manual(values = c( "#FC4E07", "#E7B800", "#99CC66", "#00AFBB"))

ggsave("Fig7D.tiff",
       plot = Fig7E,
       dpi = 600, 
       height = 2,
       width  = 2.45, 
       units = "in") 

ggsave("Fig7E.tiff",
       plot = Fig7F,
       dpi = 600, 
       height = 2,
       width  = 4.45,
       units = "in") 

ggsave("Fig7F.tiff",
       plot = Fig7G,
       dpi = 600, 
       height = 2,
       width  = 2.45, 
       units = "in") 

ggsave("Fig7G.tiff",
       plot = Fig7H,
       dpi = 600, 
       height = 2,
       width  = 2.45, 
       units = "in") 


# Fig7H-M
Fig7H <- ggboxplot(lengthContig_output_IGH_2, x="IGHV_new", y="length", color = "IGHV_new", size=0.5, add = c("boxplot", "mean_sd")) +
  theme_bw(base_size = 10)+
  theme(legend.position="none")+
  labs(x ="IGHV", y = "CDRH3 length")+
  theme(axis.text.x=element_text(angle=90))

Fig7I<- ggboxplot(lengthContig_output_IGH_2, x="IGHC", y="length", color = "IGHC", size=0.5, add = c("boxplot", "mean_sd")) +
  theme_bw(base_size = 10)+
  theme(legend.position="none")+
  labs(x ="IGHC", y = "CDRH3 length")+
  theme(axis.text.x=element_text(angle=90))

Fig7J <- ggboxplot(subset(all_clone_pid_baseline_IGH, !is.na(mutational_state)), x = "mutational_state",  y = "CDRH3_length", color = "mutational_state", size=0.5, add = c("boxplot", "mean_sd"))+
  theme_bw(base_size = 10)+
  theme(legend.position="none", 
        axis.text.x = element_text(angle=90))+
  xlab("SHM levels")+
  ylab("CDRH3 length")+
  scale_color_manual(values = c( "#FC4E07", "#E7B800", "#99CC66", "#00AFBB"))+
  scale_x_discrete (labels = c("highly mutated" = "High", 
                               "moderately mutated" = "Moderate", 
                               "lowly mutated" = "Low", 
                               "unmutated" = "Unmutated"))

Fig7K <- ggboxplot(lengthContig_output_IGH_3, x="IGHV_new", y="hydrophobicity_score", color = "IGHV_new", size=0.5, add = c("boxplot", "mean_sd")) +
  theme_bw(base_size = 10)+
  theme(legend.position="none")+
  labs(x ="IGHV", y = "Hydrophobicity score")+
  theme(axis.text.x=element_text(angle=90), 
        axis.title.y = element_text(size = 8))
    
Fig7L <- ggboxplot(lengthContig_output_IGH_3, x="IGHC", y="hydrophobicity_score", color = "IGHC", size=0.5, add = c("boxplot", "mean_sd")) +
  theme_bw(base_size = 10)+
  theme(legend.position="none")+
  labs(x ="IGHC", y = "Hydrophobicity score")+ 
  theme(axis.text.x=element_text(angle=90), 
        axis.title.y = element_text(size = 8))

Fig7M <- ggboxplot(subset(all_clone_pid_baseline_IGH, !is.na(mutational_state)), x = "mutational_state",  y = "hydrophobicity_score", color = "mutational_state", size=0.5, add = c("boxplot", "mean_sd"))+
  theme_bw(base_size=10)+
  theme(legend.position="none")  +
  xlab("SHM levels")+
  ylab("Hydrophobicity \nscore")+
  scale_color_manual(values = c( "#FC4E07", "#E7B800", "#99CC66", "#00AFBB"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        axis.title.y = element_text(size = 8))+
  scale_x_discrete (labels = c("highly mutated" = "High", 
                               "moderately mutated" = "Moderate", 
                               "lowly mutated" = "Low", 
                               "unmutated" = "Unmutated"))

ggsave("Fig7H.tiff",
       plot = Fig7I,
       dpi = 600, 
       height = 2,
       width  = 2.3, 
       units = "in") 

ggsave("Fig7I.tiff",
       plot = Fig7J,
       dpi = 600, 
       height = 2,
       width  = 2.3, 
       units = "in") 

ggsave("Fig7J.tiff",
       plot = Fig7K,
       dpi = 600, 
       height = 2,
       width  = 2.3,
       units = "in") 

ggsave("Fig7K.tiff",
       plot = Fig7L,
       dpi = 600, 
       height = 2,
       width  = 2.3, 
       units = "in") 

ggsave("Fig7L.tiff",
       plot = Fig7M,
       dpi = 600, 
       height = 2,
       width  = 2.3, 
       units = "in") 

ggsave("Fig7M.tiff",
       plot = Fig7N,
       dpi = 600, 
       height = 2,
       width  = 2.3, 
       units = "in") 

# Fig8 -------------
Fig8A <- Seurat::DimPlot(vgm_baseline_GEX_downsample,reduction = "umap", group.by = "VDJ_available", shuffle = T)+
  guides(color = "none", fill = "none")+
  theme_bw(base_size = 8) +
  ggtitle ("VDJ chain availability")

pt <- table(vgm_baseline_GEX_downsample$VDJ_available, vgm_baseline_GEX_downsample$RM_ID, vgm_baseline_GEX_downsample$seurat_clusters)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)
pt$Var2 <- factor (pt$Var2, levels = c("13F", "47G", "69G", "43F", "83H", "84H"))


Fig8B <- ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 8) +
  geom_col(position = "fill", width = 0.5) +
  xlab("RM_ID") +
  ylab("Proportion") +
  theme(legend.title = element_blank(), legend.position = "none")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_fill_discrete (labels = c("Unavailable", "Available"))+    
  scale_x_discrete (labels = c("13F" = "RM_1", 
                               "47G" = "RM_2",
                               "69G" = "RM_3",
                               "43F" = "RM_4",
                               "83H" = "RM_5", 
                               "84H" = "RM_6"))+
  
  ggplot(pt, aes(x = Var3, y = Freq*100, fill = Var1)) +
  theme_bw(base_size = 8) +
  geom_col(position = "fill", width = 0.5) +
  xlab("B cell cluster") +
  ylab("Proportion") +
  theme(legend.title = element_blank())+
    scale_fill_discrete (labels = c("Unavailable", "Available"))


Fig8C <- ggplot(pt, aes(x = Var3, y = Freq, fill = Var1)) +
  theme_bw(base_size = 8) +
  geom_col(position = "fill", width = 0.5) +
  xlab("B cell cluster") +
  ylab("Proportion of VDJ availability") +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 5))+
  facet_grid(~Var2, labeller = labeller (Var2 =
                                           c("13F" = "RM_1",
                                             "47G" = "RM_2",
                                             "69G" = "RM_3",
                                             "43F" = "RM_4",
                                             "83H" = "RM_5",
                                             "84H" = "RM_6")))+
  scale_fill_discrete (labels = c("Unavailable", "Available"))

ggsave("Fig8A.tiff",
       plot = Fig8A,
       dpi = 600, 
       height = 2,
       width  = 2, 
       units = "in") 

ggsave("Fig8B.tiff",
       plot = Fig8B,
       dpi = 600, 
       height = 2,
       width  = 4.9, 
       units = "in") 

ggsave("Fig8C.tiff",
       plot = Fig8C,
       dpi = 600, 
       height = 1.7,
       width  = 6.9,
       units = "in") 


# Fig8D
vgm_baseline_GEX_subset <- subset(vgm_baseline_GEX_downsample, Nr_of_VDJ_chains == 1 & Nr_of_VJ_chains == 1)
table(vgm_baseline_GEX_subset$Nr_of_VDJ_chains)

vgm_baseline_GEX_remove_VDJ_cgene_missing_value <- subset(vgm_baseline_GEX_subset, subset = VDJ_cgene == "IGHM" | VDJ_cgene == "IGHD"| VDJ_cgene == "IGHG1"| VDJ_cgene == "IGHG2"| VDJ_cgene == "IGHG3"| VDJ_cgene == "IGHG4"| VDJ_cgene == "IGHA"| VDJ_cgene == "IGHE")

table(vgm_baseline_GEX_remove_VDJ_cgene_missing_value$Nr_of_VDJ_chains)

Idents(vgm_baseline_GEX_remove_VDJ_cgene_missing_value) <- "RNA_snn_res_B_0.4" 

Fig8D <- DimPlot(vgm_baseline_GEX_remove_VDJ_cgene_missing_value,reduction = "umap", group.by = "VDJ_cgene", shuffle = T) +
  theme_bw(base_size = 8) +
  ggtitle ("IGHC genes on \nB cell Seurat object")

# Fig8E
vgm_baseline_GEX_four_heavy_chains <-  subset(vgm_baseline_GEX_remove_VDJ_cgene_missing_value, subset = VDJ_cgene == "IGHM" | VDJ_cgene == "IGHD"| VDJ_cgene == "IGHG1"| VDJ_cgene == "IGHA")

vgm_baseline_GEX_four_heavy_chains$VDJ_cgene <- factor (vgm_baseline_GEX_four_heavy_chains$VDJ_cgene, levels = c("IGHD", "IGHM", "IGHG1", "IGHA"))

Idents(vgm_baseline_GEX_four_heavy_chains) <- "VDJ_cgene" 
find_all_markers_VDJ_cgene <- FindAllMarkers(object = vgm_baseline_GEX_four_heavy_chains,
                                             only.pos = FALSE,
                                             logfc.threshold = 0.25)

top_genes_lfc_output_function <- function (find_all_markers_dataset, slice_val) {
  
  top_genes_lfc_output_dataset <- find_all_markers_dataset %>%
    group_by(cluster) %>%
    filter (avg_log2FC > 0.40) %>%
    filter (pct.1 > 0.27) %>%
    arrange(cluster, desc(avg_log2FC), p_val_adj) %>%
    dplyr :: slice(1:slice_val) %>%
    mutate(category = row_number()) %>% 
    ungroup()
  
  return (top_genes_lfc_output_dataset)
  
}

top_genes_lfc_output <- top_genes_lfc_output_function(find_all_markers_VDJ_cgene, 15)

Idents(vgm_baseline_GEX_four_heavy_chains) <- "VDJ_cgene"

alldata <- ScaleData(vgm_baseline_GEX_four_heavy_chains, features = as.character(unique(top_genes_lfc_output$gene)), assay = "RNA")

alldata$VDJ_cgene <- factor (alldata$VDJ_cgene, levels = c("IGHA", "IGHG1", "IGHM", "IGHD"))

Fig8E <-DotPlot(subset(alldata, downsample = 1000), features = (as.character(unique(top_genes_lfc_output$gene))), group.by = "VDJ_cgene", assay = "RNA")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))+
  scale_colour_gradient2(low = "blue", mid = "white", high = "red")& 
  theme(text  = element_text(size = 8), 
        axis.text  = element_text(size = 8), 
        axis.text.x  = element_text(size = 7)) +
  theme (legend.position = "bottom")

ggsave("Fig8D.tiff",
       plot = Fig8D,
       dpi = 600, 
       height = 2,
       width  = 2.5, 
       units = "in")

ggsave("Fig8E.tiff",
       plot = Fig8E,
       dpi = 600, 
       height = 3,
       width  = 6.9, 
       units = "in")


# Fig9 --------------

vgm_gene_subtype_function <- function (IG_chain_gene_segment) {
  dataset <- vgm_baseline_VDJ %>%
    filter(!!sym(IG_chain_gene_segment) != '') %>%
    filter(seurat_clusters != '') %>%
    group_by(orig.ident, sample_id, RM_ID, seurat_clusters, !!sym(IG_chain_gene_segment)) %>%
    summarize(n = n())  %>%
    group_by(seurat_clusters, sample_id) %>%
    mutate(ttl = sum(n), Freq = n/ttl*100)
  
  dataset_2 <- dataset %>%
    group_by (!!sym(IG_chain_gene_segment))%>%
    summarize(meanvalue = mean (n))%>%
    arrange(desc(meanvalue))
  
  dataset[[IG_chain_gene_segment]] <- factor(dataset[[IG_chain_gene_segment]], levels = dataset_2[[IG_chain_gene_segment]])
  
  return(dataset)
  
}

vgm_gene_subtype_function_VDJ_cgene_1 <- function (IG_chain_gene_segment) {
  dataset <- vgm_baseline_VDJ %>%
    filter(!!sym(IG_chain_gene_segment) != '') %>%
    filter(seurat_clusters != '') %>%
    mutate (class_switch_status = VDJ_cgene)%>%
    mutate(across(class_switch_status, ~str_replace (., '^IGHM$', 'unswitched'))) %>%
    mutate(across(class_switch_status, ~str_replace (., '^IGHD$', 'unswitched'))) %>%
    mutate(across(class_switch_status, ~str_replace (., '^IGHG1$', 'switched'))) %>%
    mutate(across(class_switch_status, ~str_replace (., '^IGHG2$', 'switched'))) %>%
    mutate(across(class_switch_status, ~str_replace (., '^IGHG3$', 'switched'))) %>%
    mutate(across(class_switch_status, ~str_replace (., '^IGHG4$', 'switched'))) %>%
    mutate(across(class_switch_status, ~str_replace (., '^IGHA$', 'switched'))) %>%
    mutate(across(class_switch_status, ~str_replace (., '^IGHE$', 'switched'))) %>%
    mutate(VDJ_cgene = factor (VDJ_cgene, levels = c("IGHM", "IGHD", "IGHG1", "IGHG2", "IGHG3","IGHG4", "IGHA", "IGHE")))%>%
    group_by(orig.ident, sample_id,RM_ID, seurat_clusters, !!sym(IG_chain_gene_segment), class_switch_status) %>%
    summarize(n = n())  %>%
    group_by(seurat_clusters, sample_id) %>%
    mutate(ttl = sum(n), Freq = n/ttl*100)
  
  # dataset_2 <- dataset %>%
  #   group_by (!!sym(IG_chain_gene_segment))%>%
  #   summarize(meanvalue = mean (n))%>%
  #   arrange(desc(meanvalue))
  # 
  # dataset[[IG_chain_gene_segment]] <- factor(dataset[[IG_chain_gene_segment]], levels = dataset_2[[IG_chain_gene_segment]])
  
  return(dataset)
  
}

vgm_gene_family_function <- function (dataset, IG_chain_gene_segment){
  vgm_all_total <- dataset %>%
    mutate (substring_gene = substr(!!sym(IG_chain_gene_segment),1,6),
            substring_gene = gsub ("S", "", substring_gene), 
            substring_gene = gsub("-.*", "", substring_gene), 
            substring_gene = gsub(";*", "", substring_gene))%>%
    group_by(seurat_clusters, substring_gene) %>%
    summarize(n = sum(n))%>%
    group_by(seurat_clusters) %>%
    mutate(total_2 = sum(n), Freq = n/total_2*100)
  
  return (vgm_all_total)
}

vgm_VDJ_vgene <- vgm_gene_subtype_function("VDJ_vgene")
vgm_VDJ_dgene <- vgm_gene_subtype_function("VDJ_dgene")
vgm_VDJ_jgene <- vgm_gene_subtype_function("VDJ_jgene")
vgm_VDJ_cgene <- vgm_gene_subtype_function_VDJ_cgene_1("VDJ_cgene")

vgm_VJ_kvgene <- vgm_gene_subtype_function("VJ_vgene")%>%
  filter(substr(VJ_vgene, 1,3 ) == "IGK")
vgm_VJ_kjgene <- vgm_gene_subtype_function("VJ_jgene")%>%
  filter(substr(VJ_jgene, 1,3 ) == "IGK")
vgm_VJ_kcgene <- vgm_gene_subtype_function("VJ_cgene")%>%
  filter(substr(VJ_cgene, 1,3 ) == "IGK")

vgm_VJ_lvgene <- vgm_gene_subtype_function("VJ_vgene")%>%
  filter(substr(VJ_vgene, 1,3 ) == "IGL")
vgm_VJ_ljgene <- vgm_gene_subtype_function("VJ_jgene")%>%
  filter(substr(VJ_jgene, 1,3 ) == "IGL")
vgm_VJ_lcgene <- vgm_gene_subtype_function("VJ_cgene")%>%
  filter(substr(VJ_cgene, 1,3 ) == "IGL")

vgm_gene_family_VDJ_vgene <- vgm_gene_family_function(vgm_VDJ_vgene, "VDJ_vgene")
vgm_gene_family_VDJ_dgene <- vgm_gene_family_function(vgm_VDJ_dgene, "VDJ_dgene")
vgm_gene_family_VDJ_jgene <- vgm_gene_family_function(vgm_VDJ_jgene, "VDJ_jgene")
vgm_gene_family_VDJ_cgene <- vgm_gene_family_function(vgm_VDJ_cgene, "VDJ_cgene")

vgm_gene_family_VJ_kvgene <- vgm_gene_family_function(vgm_VJ_kvgene, "VJ_vgene")
vgm_gene_family_VJ_kjgene <- vgm_gene_family_function(vgm_VJ_kjgene, "VJ_jgene")
vgm_gene_family_VJ_kcgene <- vgm_gene_family_function(vgm_VJ_kcgene, "VJ_cgene")

vgm_gene_family_VJ_lvgene <- vgm_gene_family_function(vgm_VJ_lvgene, "VJ_vgene") %>%
  mutate(substring_gene = factor (substring_gene, levels = c("IGLV1", "IGLV2", "IGLV3", "IGLV4", "IGLV5", "IGLV6", "IGLV7", "IGLV8", "IGLV9", "IGLV10", "IGLV11")))
vgm_gene_family_VJ_ljgene <- vgm_gene_family_function(vgm_VJ_ljgene, "VJ_jgene")
vgm_gene_family_VJ_lcgene <- vgm_gene_family_function(vgm_VJ_lcgene, "VJ_cgene")


plot_1_gene_family_prop_per_cluster_function <- function (dataset, IG_chain_gene_segment) {
  plot_d_prop <- ggplot(subset(dataset, seurat_clusters %in% c(0:8, 12)), aes(x=seurat_clusters , y= n))+
    geom_col(position = "fill", width = 0.5, aes(fill = substring_gene))+
    theme_bw(base_size=8)+
    xlab("B cell cluster")+
    ylab("Proportions")+
    labs(fill=IG_chain_gene_segment) 
  return (plot_d_prop)
  
}

Fig9A <- plot_1_gene_family_prop_per_cluster_function (vgm_gene_family_VDJ_vgene, "IGHV")+
  plot_1_gene_family_prop_per_cluster_function (vgm_gene_family_VDJ_dgene, "IGHD")+
  plot_1_gene_family_prop_per_cluster_function (vgm_gene_family_VDJ_jgene, "IGHJ")

Fig9B <- plot_1_gene_family_prop_per_cluster_function (vgm_gene_family_VJ_kvgene, "IGKV")+
  plot_1_gene_family_prop_per_cluster_function (vgm_gene_family_VJ_kjgene, "IGKJ")+
  plot_1_gene_family_prop_per_cluster_function (vgm_gene_family_VJ_kcgene, "IGKC")

Fig9C <- plot_1_gene_family_prop_per_cluster_function (vgm_gene_family_VJ_lvgene, "IGLV")+
  plot_1_gene_family_prop_per_cluster_function (vgm_gene_family_VJ_ljgene, "IGLJ")+
  plot_1_gene_family_prop_per_cluster_function (vgm_gene_family_VJ_lcgene, "IGLC")


Fig9D <- ggboxplot(all_clone_pid_baseline_IGH, x="seurat_clusters", y="CDRH3_length", color = "seurat_clusters", size=0.5, add = c("boxplot", "mean_sd")) +
  labs(x ="B cell cluster", y = "CDRH3 length")+ theme(legend.position="right")+
  stat_summary(fun = "mean", geom = "point", shape = 21, size = 2, color = "red")+
  theme_bw(base_size=10)+
  theme(legend.position="none")

Fig9E <- ggboxplot(all_clone_pid_baseline_IGH, x = "seurat_clusters", y = "hydrophobicity_score", color = "seurat_clusters")+
  stat_summary(fun = "mean", geom = "point", shape = 21, size = 2, color = "red")+
  theme_bw(base_size=10)+
  theme(legend.position="none")+
  labs(x ="B cell cluster", y = "Hydrophobicity \nscore")

ggsave("Fig9A.tiff",
       plot = Fig9A,
       dpi = 600, 
       height = 2.2,
       width  = 6.9, 
       units = "in") 

ggsave("Fig9B.tiff",
       plot = Fig9B,
       dpi = 600, 
       height = 2.2,
       width  = 6.9, 
       units = "in") 

ggsave("Fig9C.tiff",
       plot = Fig9C,
       dpi = 600, 
       height = 3,
       width  = 6.9,
       units = "in") 

ggsave("Fig9D.tiff",
       plot = Fig9D,
       dpi = 600, 
       height = 2.2,
       width  = 3.45, 
       units = "in") 

ggsave("Fig9E.tiff",
       plot = Fig9E,
       dpi = 600, 
       height = 2.2,
       width  = 3.45,
       units = "in") 
