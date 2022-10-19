#this is tested using Cicero conda environment (with R 4 and Seurat 4) on a high-performance computing node with 16 cores and at least 160 gb or ram. 
library("Seurat") #load Seurat 4
library("matrixStats")
library('tidyverse')
library('magrittr')
library(ggridges)
library(viridis)
library(scico)
library('dplyr')
library('pheatmap')
library('gplots')

ObjList_filt <- readRDS("SeuratObjects/ObjList_filt.preprocessed.ADTclustered.wMeta.Rds")

dir.create("DownsampleLn6_plots")

## the B cells sorted and put in 10x lane 6 are gated on S+ cells
ObjList_filt$MergeBcells$Spositve <- ifelse(ObjList_filt$MergeBcells$Lane == "06", TRUE, FALSE)

## pre-downsampling counts per cell
ggsave("DownsampleLn6_plots/nFeature_RNA.NonDSlane6.pdf", width=5, height=4,
      plot =VlnPlot(ObjList_filt$MergeBcells, "nFeature_RNA", group.by = "Tissue", split.by = "Spositve")
)
ggsave("DownsampleLn6_plots/nCount_RNA.NonDSlane6.pdf", width=5, height=4,
		plot = VlnPlot(ObjList_filt$MergeBcells, "nCount_RNA", group.by = "Tissue", split.by = "Spositve")
)
ggsave("DownsampleLn6_plots/nCount_CITE.NonDSlane6.pdf", width=5, height=4,
		plot = VlnPlot(ObjList_filt$MergeBcells, "nCount_CITE", group.by = "Tissue", split.by = "Spositve")
)

# import the cellranger counts of mRNA and features for the Downsampled lane 6
SPosBcells_data = list()
rawdir = "H5DataFromGEO"
for(i in 1:1){
  SPosBcells_data[[i]] = Read10X_h5(paste(rawdir,"Pam1_CITE_multi5P06_DS/outs/multi/count/raw_feature_bc_matrix.h5", sep="")
  )
}
SPosBcells_SeuratObj <- CreateSeuratObject(counts = SPosBcells_data[[1]]$'Gene Expression', assay = "RNA", min.feature = 20)
SPosBcells_SeuratObj <- RenameCells(SPosBcells_SeuratObj, new.names = paste(substr(colnames(SPosBcells_SeuratObj), start = 1, stop = 17),"06", sep = ""))

Lane6cells <- colnames(ObjList_filt$MergeBcells)[which(ObjList_filt$MergeBcells$Lane == "06")]

##adding the downsampled RNA data for the lane 6 (S+ sorted B cells)
ObjList_filt$MergeBcells@assays$RNA@counts[,Lane6cells] <- SPosBcells_SeuratObj@assays$RNA@counts[,Lane6cells]
ObjList_filt$MergeBcells <- NormalizeData(ObjList_filt$MergeBcells)
ObjList_filt$MergeBcells <- ScaleData(ObjList_filt$MergeBcells)
ObjList_filt$MergeBcells$nCount_RNA = colSums(x = ObjList_filt$MergeBcells, slot = "counts")  # nCount_RNA
ObjList_filt$MergeBcells$nFeature_RNA = colSums(x = GetAssayData(object = ObjList_filt$MergeBcells, slot = "counts") > 0)  # nFeatureRNA

ggsave("DownsampleLn6_plots/nFeature_RNA.DSlane6.pdf", width=5, height=4,
		plot = VlnPlot(ObjList_filt$MergeBcells, "nFeature_RNA", group.by = "Tissue", split.by = "Spositve")
)
ggsave("DownsampleLn6_plots/nCount_RNA.DSlane6.pdf", width=5, height=4,
		plot = VlnPlot(ObjList_filt$MergeBcells, "nCount_RNA", group.by = "Tissue", split.by = "Spositve")
)
ggsave("DownsampleLn6_plots/nCount_CITE.DSlane6.pdf", width=5, height=4,
		plot = VlnPlot(ObjList_filt$MergeBcells, "nCount_CITE", group.by = "Tissue", split.by = "Spositve")
)


COVIDsampleBcells <- subset(ObjList_filt$MergeBcells, Class == "COVIDpos")

Idents(COVIDsampleBcells) <- "Spositve"

dir.create("DownsampleLn6_plots/CellDE")

PBMCSposvsNegDEgenes = FindMarkers(subset(COVIDsampleBcells, Tissue == "PBMC"), ident.1 = "TRUE", ident.2 = "FALSE", test.use = "MAST", latent.vars = "Donor")
AdenoidSposvsNegDEgenes = FindMarkers(subset(COVIDsampleBcells, Tissue == "Adenoid"), ident.1 = "TRUE", ident.2 = "FALSE", test.use = "MAST", latent.vars = "Donor")
TonsilSposvsNegDEgenes = FindMarkers(subset(COVIDsampleBcells, Tissue == "Tonsil"), ident.1 = "TRUE", ident.2 = "FALSE", test.use = "MAST", latent.vars = "Donor")

write.csv(PBMCSposvsNegDEgenes, "DownsampleLn6_plots/CellDE/PBMCSposvsNegDEgenes_DownSample_MASTDonorLatentVar.csv")
write.csv(AdenoidSposvsNegDEgenes, "DownsampleLn6_plots/CellDE/AdenoidSposvsNegDEgenes_DownSample_MASTDonorLatentVar.csv")
write.csv(TonsilSposvsNegDEgenes, "DownsampleLn6_plots/CellDE/TonsilSposvsNegDEgenes_DownSample_MASTDonorLatentVar.csv")

COVIDsampleBcells$Sample_Tissue_Class_Spos <- paste(COVIDsampleBcells$Sample_Tissue_Class, COVIDsampleBcells$Spositve)

Idents(COVIDsampleBcells) <- "adt_snn_res.1"
clustRNAmarkers <- readRDS("clustRNAmarkers.Rds")
celltypeheatgenes = clustRNAmarkers[[6]] %>% group_by(cluster) %>% top_n(10, avg_log2FC)

BmemGenes <- c("BAIAP3", "COL4A3", "GDPD5", "SAMSN1", "CD80", "RORA", "COL4A4", "TNFRSF13B", "RASSF6",
				"MUC16", "TEX9", "C10orf128", "TOX", "AIM2", "TRERF1", "TGM2", "TRPV3",
				"SPRY1", "GABARAPL1", "APBB2", "ZBTB16", "PCDH9", "IKZF2", "TCL1A")

COVIDsampleBcells <- ScaleData(COVIDsampleBcells, features = BmemGenes)
ggsave("DownsampleLn6_plots/DotPlot.BmemGenes.pdf", width=8, height=6,
		DotPlot(COVIDsampleBcells, features = BmemGenes, group.by = "Sample_Tissue_Class_Spos") + RotatedAxis()
)

shared_s1_clone_cells <- read.csv("shared_s1_clone_cells.csv")
shared_s1_clone_cells$Cluster <- COVIDsampleBcells$seurat_clusters[shared_s1_clone_cells$cell_id]

write.csv(shared_s1_clone_cells, "shared_s1_clone_cells_withCITEclust.csv")


COVIDsampleBcells$SharedClone <- ifelse(colnames(COVIDsampleBcells) %in% shared_s1_clone_cells$cell_id, TRUE, FALSE)
table(COVIDsampleBcells$seurat_clusters, COVIDsampleBcells$SharedClone)
COVIDsampleBcells <- NormalizeData(COVIDsampleBcells)
TACovidBcells <- subset(COVIDsampleBcells, Tissue == "PBMC", invert = TRUE)
TACovidBcells$SPos_Shared <- paste(TACovidBcells$Spositve, TACovidBcells$SharedClone,  sep = "_")
TACovidBcells$SPos_Shared_Cluster <- paste(TACovidBcells$Spositve, TACovidBcells$SharedClone, TACovidBcells$seurat_clusters, sep = "_")
TACovidBcells$SPos_Cluster <- paste(TACovidBcells$Spositve, TACovidBcells$seurat_clusters, sep = "_")

Idents(TACovidBcells) <- "SPos_Shared_Cluster"

Clust2_SnegvsSposNonsharedDE <- FindMarkers(TACovidBcells, ident.1 = "TRUE_FALSE_2", ident.2 = "FALSE_FALSE_2", test.use = "MAST", latent.vars = "Donor")
Clust2_SnegvsSposSharedDE <- FindMarkers(TACovidBcells, ident.1 = "TRUE_TRUE_2", ident.2 = "FALSE_FALSE_2", test.use = "MAST", latent.vars = "Donor")
Clust2_SposNonsharedvsSposSharedDE <- FindMarkers(TACovidBcells, ident.1 = "TRUE_TRUE_2", ident.2 = "TRUE_FALSE_2", test.use = "MAST", latent.vars = "Donor")

Idents(TACovidBcells) <- "SPos_Cluster"

Clust2_SnegvsSposSharedandNonSharedDE <- FindMarkers(TACovidBcells, ident.1 = "TRUE_2", ident.2 = "FALSE_2", test.use = "MAST", latent.vars = "Donor")
Clust9_SnegvsSposSharedandNonSharedDE <- FindMarkers(TACovidBcells, ident.1 = "TRUE_9", ident.2 = "FALSE_9", test.use = "MAST", latent.vars = "Donor")

VlnPlot(subset(TACovidBcells, seurat_clusters == "2"), "nCount_RNA", group.by = "Tissue", split.by = "Spositve")
VlnPlot(subset(TACovidBcells, seurat_clusters == "2"), "percent.mito", group.by = "Tissue", split.by = "Spositve")

write.csv(Clust2_SnegvsSposNonsharedDE, "DownsampleLn6_plots/CellDE/Clust2_SnegvsSposNonsharedDE_DownSample_MASTDonorLatentVar-2.csv")
write.csv(Clust2_SnegvsSposSharedDE, "DownsampleLn6_plots/CellDE/Clust2_SnegvsSposSharedDE_DownSample_MASTDonorLatentVar-2.csv")
write.csv(Clust2_SposNonsharedvsSposSharedDE, "DownsampleLn6_plots/CellDE/Clust2_SposNonsharedvsSposSharedDE_DownSample_MASTDonorLatentVar-2.csv")
write.csv(Clust2_SnegvsSposSharedandNonSharedDE, "DownsampleLn6_plots/CellDE/Clust2_SnegvsSposSharedandNonSharedDE_DownSample_MASTDonorLatentVar-2.csv")
write.csv(Clust9_SnegvsSposSharedandNonSharedDE, "DownsampleLn6_plots/CellDE/Clust9_SnegvsSposSharedandNonSharedDE_DownSample_MASTDonorLatentVar-2.csv")

Clust2_SnegvsSposNonsharedDE <- read.csv("DownsampleLn6_plots/CellDE/Clust2_SnegvsSposNonsharedDE_DownSample_MASTDonorLatentVar-2.csv", header=TRUE)
Clust2_SnegvsSposSharedDE <- read.csv("DownsampleLn6_plots/CellDE/Clust2_SnegvsSposSharedDE_DownSample_MASTDonorLatentVar-2.csv", header=TRUE)
Clust2_SposNonsharedvsSposSharedDE <- read.csv("DownsampleLn6_plots/CellDE/Clust2_SposNonsharedvsSposSharedDE_DownSample_MASTDonorLatentVar-2.csv", header=TRUE)
Clust2_SnegvsSposSharedandNonSharedDE <- read.csv("DownsampleLn6_plots/CellDE/Clust2_SnegvsSposSharedandNonSharedDE_DownSample_MASTDonorLatentVar-2.csv", header=TRUE)
Clust9_SnegvsSposSharedandNonSharedDE <- read.csv("CellDE/Clust9_SnegvsSposSharedandNonSharedDE_DownSample_MASTDonorLatentVar-2.csv", header=TRUE)

##VolcanoPlotTop = function(contrast.result.list, contrast.name, save.path,
 ##                                                  ngenestolabel = 8, fig.height, fig.width) {
    library(ggrepel)
  
      ## control how many genes to label here in the top10 column
    contrast.result.list =
        Clust2_SnegvsSposSharedandNonSharedDE %>%
        arrange(p_val_adj)
    contrast.result.list = contrast.result.list %>% mutate(top10 = if_else(X %in% contrast.result.list[1:25, ]$X, "1", "0"))

      ## Define n genes plot by subsetting contrast result list here ^^^^^

      p = ggplot(data = contrast.result.list,
                 aes(x = avg_log2FC, y= -log(p_val))) +
        geom_point(aes(col=top10, alpha = 0.9), show.legend = FALSE) +
        scale_color_manual(values=c("darkgrey", "firebrick4")) +
        geom_text_repel(data=contrast.result.list %>% filter(top10 == "1"),
                        aes(label=X), size = 4.5, fontface = "bold", nudge_x = -0.1,
                        segment.size = 0.1,point.padding = unit(0.2, "lines")) +
        ggtitle("Clust2_SnegvsSposSharedandNonSharedDE") +
        labs(x = "avg_log2FC", y = "-log(p_val)") +
        theme(legend.position = "none") +
        theme_classic()
        ggsave(plot = p, filename = "DownsampleLn6_plots/CellDE/Clust2_SnegvsSposSharedandNonSharedDE_Volcano_DownSample_MAST_donorLatentVar-2.pdf",
                width = 6, height = 6)
      #print(p)
  
Cluster2CovidBcells <- subset(TACovidBcells, seurat_clusters == "2")
Cluster2CovidBcells <- ScaleData(Cluster2CovidBcells, features = contrast.result.list[1:50, ]$X)
contrast.result.list.top40 <- contrast.result.list %>% slice_min(p_val_adj, n=40)
contrast.result.list.top40 <- contrast.result.list.top40 %>% arrange(desc(avg_log2FC))

ggsave(filename = "DownsampleLn6_plots/CellDE/Clust2_SnegvsSposSharedandNonSharedDE_Top40genes_Heatmap_DownSample_MAST_donorLatentVar-2.pdf",
                width = 6, height = 6,
	plot = DoHeatmap(Cluster2CovidBcells, features = contrast.result.list.top40$X, group.by="Spositve")
)
Cluster2CovidBcells$S_Donor_Tissue <- paste(Cluster2CovidBcells$Spositve, Cluster2CovidBcells$Donor, Cluster2CovidBcells$Tissue, sep="_")
ggsave(filename = "DownsampleLn6_plots/CellDE/Clust2_SnegvsSposSharedandNonSharedDE_Top40genes_Dotplot_DownSample_MAST_donorLatentVar-2.pdf",
                width = 12, height = 6,
	plot = DotPlot(Cluster2CovidBcells, features = contrast.result.list.top40$X, group.by="S_Donor_Tissue") +RotatedAxis()
)

MZBcellGenes_GSM538210_500_OverlapwithDE <- c("TRAF5", "RPL39",	"HLA-DOA", "HLA-DOB", "STRBP", "TMEM154", "BCL11A", "CD83", "FCER2",
                                               "CIITA", "IFI30", "BCAR3", "FGD2", "BMP2K", "DENND3", "KMO", "CTSS", "ARHGAP24",
                                               "SPIB", "SWAP70", "GAPT", "STAP1", "BLK", "NEDD9", "CXCR5", "MALT1",	"FCRL1",
											   "RAB30", "BTK", "TMEM163", "KYNU", "PLD4", "ICOSLG", "KTN1", "RALGPS2", "EBF1",
 												"CYB561A3", "MS4A1", "CD38", "TLR1", "CD40", "TNFRSF13C", "MARCKS", "RPL13", "SP140")
Cluster2CovidBcells <- ScaleData(Cluster2CovidBcells, features = MZBcellGenes_GSM538210_500_OverlapwithDE)

ggsave(filename = "DownsampleLn6_plots/CellDE/Cluster2CovidBcellsMZBcellGenes_GSM538210_500_OverlapwithDE.DotPlot.pdf",
		DotPlot(Cluster2CovidBcells, features = MZBcellGenes_GSM538210_500_OverlapwithDE) + RotatedAxis()
)
ggsave(filename = "DownsampleLn6_plots/CellDE/Cluster2CovidBcellsMZBcellGenes_GSM538210_500_OverlapwithDE.Heatmap.pdf",
		DoHeatmap(Cluster2CovidBcells, features = MZBcellGenes_GSM538210_500_OverlapwithDE)
)

ArtS_GCBgenes <- list(Bmem = read.table("BmemExclusiveGenes_ArtSchaffer.txt", sep="\t", header=FALSE),
					LZGC = read.table("LZGCBExclusiveGenes_ArtSchaffer.txt", sep="\t", header=FALSE),
					DZGC = read.table("DZGCBExclusiveGenes_ArtSchaffer.txt", sep="\t", header=FALSE))
					
Cluster2CovidBcells <- ScaleData(Cluster2CovidBcells, features = as.character(unlist(ArtS_GCBgenes)))

Cluster2CovidBcells <- AddModuleScore(Cluster2CovidBcells, features = list(intersect(ArtS_GCBgenes[[1]]$V1, rownames(Cluster2CovidBcells))), name = "Bmem", nbin=12)					
Cluster2CovidBcells <- AddModuleScore(Cluster2CovidBcells, features = list(intersect(ArtS_GCBgenes[[2]]$V1, rownames(Cluster2CovidBcells))), name = "LZGC", nbin=12)					
Cluster2CovidBcells <- AddModuleScore(Cluster2CovidBcells, features = list(intersect(ArtS_GCBgenes[[3]]$V1, rownames(Cluster2CovidBcells))), name = "DZGC", nbin=12)					

ggsave(filename = "DownsampleLn6_plots/CellDE/Cluster2CovidBcellsSchafferBcellGenesScore.pdf",
		VlnPlot(Cluster2CovidBcells, features = c("Bmem1","LZGC1","DZGC1"))		
)
### Create heatmap of selected pathway genes
library(ComplexHeatmap)
library(viridis)
viridis(9)

				
Milpeid2020FrontImmunolgenesets <- list(GC = c("TCL1A","ELL3","STMN1","SERPINA9","LPP","HMGB2","PTTG1","KIAA010","AICDA","NEIL1"),
										Mem = c("CD52","BANK1","FCMR","CAPG","SELL","CCR7","IFITM2","PARP15","TRBC2","ARHGAP24"),
										PB.PC = c("MZB1","SSR4","DERL3","XBP1","HSP90B1","PRDX4","FKBP11","ITM2C","DNAJB9","SLAMF7"))	
										
CovidBcells.forheatmap <- subset(COVIDsampleBcells, Spositve=="TRUE" & Tissue %in% c("Adenoid","Tonsil","PBMC"))
CovidBcells.forheatmap <- ScaleData(CovidBcells.forheatmap, features = as.character(unlist(Milpeid2020FrontImmunolgenesets)))
label <- which(rownames(GetAssayData(CovidBcells.forheatmap, assay = "RNA", slot = "scale.data")) %in% as.character(unlist(Milpeid2020FrontImmunolgenesets)))
label_genes <- rownames(GetAssayData(CovidBcells.forheatmap, assay = "RNA", slot = "scale.data"))[label]

library("circlize")
col_fun = colorRamp2(c(-1.5, 0, 1.5), c("#FF00FF", "#000000", "#FFFF00"))

DiscretePalette(2, palette = "alphabet")

ha.COVIDBcells = HeatmapAnnotation(
  SharedClone = CovidBcells.forheatmap$SharedClone,
  Tissue = CovidBcells.forheatmap$Tissue,
  col = list(SharedClone = c("FALSE" = "cornsilk", "TRUE" = "black"),
               Tissue = c("Adenoid" = "yellow", "Tonsil" = "purple", "PBMC"="red"))
  )

GC = rownames(GetAssayData(CovidBcells.forheatmap, assay = "RNA", slot = "scale.data")) %in% Milpeid2020FrontImmunolgenesets$GC
Mem = rownames(GetAssayData(CovidBcells.forheatmap, assay = "RNA", slot = "scale.data")) %in% Milpeid2020FrontImmunolgenesets$Mem
PB.PC = rownames(GetAssayData(CovidBcells.forheatmap, assay = "RNA", slot = "scale.data")) %in% Milpeid2020FrontImmunolgenesets$PB.PC

pdf("DownsampleLn6_plots/SpositiveCovidBcells.Milpeid2020FrontImmunolgenesets_donwsample.pdf", width = 12, height = 8)
Heatmap(as.matrix(GetAssayData(CovidBcells.forheatmap, assay = "RNA", slot = "scale.data")), name = "COVIDBcells", show_column_names = FALSE,
        top_annotation = ha.COVIDBcells, cluster_columns = FALSE, cluster_rows= TRUE, column_title_rot = 90,
        column_split = CovidBcells.forheatmap$adt_snn_res.1,
        row_km = 3,
        col = col_fun,
        row_names_gp = gpar(fontsize = 10)
) +
  Heatmap(GC + 0, name = "GC", col = c("0" = "white", "1" = "purple"),
          show_heatmap_legend = FALSE, width = unit(2, "mm"), column_names_gp = gpar(fontsize = 6)) +
  Heatmap(Mem + 0, name = "Mem", col = c("0" = "white", "1" = "red"),
          show_heatmap_legend = FALSE, width = unit(2, "mm"), column_names_gp = gpar(fontsize = 6)) +
  Heatmap(PB.PC + 0, name = "PB/PC", col = c("0" = "white", "1" = "blue"),
          show_heatmap_legend = FALSE, width = unit(2, "mm"), column_names_gp = gpar(fontsize = 6)) +
  rowAnnotation(foo = anno_mark(at = label,
                              labels = label_genes,
                              labels_gp = gpar(fontsize = 14)))
dev.off()

CovidBcells.forheatmap <- ScaleData(CovidBcells.forheatmap, assay = "CITE", features = rownames(CovidBcells.forheatmap@assays$CITE))

col_fun2 = colorRamp2(c( 0, 10), c("#EBEDEF", "#1B2631"))

pdf("DownsampleLn6_plots/SpositiveCovidBcells.CITEMarkers_DownsampleRNA.pdf", width = 12, height = 3)
Heatmap(as.matrix(GetAssayData(CovidBcells.forheatmap, assay = "CITE", slot = "data"))[ c("IgD", "CD27", "CD38"), ], name = "dsb CITEseq", show_column_names = FALSE,
        cluster_columns = FALSE, cluster_rows = FALSE,
        column_split = CovidBcells.forheatmap$adt_snn_res.1,
        col = col_fun2,
        row_names_gp = gpar(fontsize = 10)
)
dev.off()


pdf("DownsampleLn6_plots/SpositiveCovidBcells.Milpeid2020FrontImmunolgenesets.clusterColumns_donwsample.pdf", width = 12, height = 8)
Heatmap(as.matrix(GetAssayData(CovidBcells.forheatmap, assay = "RNA", slot = "scale.data")), name = "COVIDBcells", show_column_names = FALSE,
        top_annotation = ha.COVIDBcells, cluster_columns = TRUE, column_title_rot = 90,
        column_split = CovidBcells.forheatmap$adt_snn_res.1,
        col = col_fun,
        row_km = 3, 
        row_names_gp = gpar(fontsize = 5)
) +
  Heatmap(GC + 0, name = "GC", col = c("0" = "white", "1" = "purple"),
          show_heatmap_legend = FALSE, width = unit(2, "mm"), column_names_gp = gpar(fontsize = 10)) +
  Heatmap(Mem + 0, name = "Mem", col = c("0" = "white", "1" = "red"),
          show_heatmap_legend = FALSE, width = unit(2, "mm"), column_names_gp = gpar(fontsize = 10)) +
  Heatmap(PB.PC + 0, name = "PB/PC", col = c("0" = "white", "1" = "blue"),
          show_heatmap_legend = FALSE, width = unit(2, "mm"), column_names_gp = gpar(fontsize = 10)) +
  rowAnnotation(foo = anno_mark(at = label,
                              labels = label_genes,
                              labels_gp = gpar(fontsize = 10)))
dev.off()

##heatmap for all the COVID sample B cells

CovidBcells.forheatmap <- subset(COVIDsampleBcells)
CovidBcells.forheatmap <- ScaleData(CovidBcells.forheatmap, features = as.character(unlist(Milpeid2020FrontImmunolgenesets)))
label <- which(rownames(GetAssayData(CovidBcells.forheatmap, assay = "RNA", slot = "scale.data")) %in% as.character(unlist(Milpeid2020FrontImmunolgenesets)))
label_genes <- rownames(GetAssayData(CovidBcells.forheatmap, assay = "RNA", slot = "scale.data"))[label]

library("circlize")
col_fun = colorRamp2(c(-1.5, 0, 1.5), c("#FF00FF", "#000000", "#FFFF00"))
ha.COVIDBcells = HeatmapAnnotation(
  SharedClone = CovidBcells.forheatmap$SharedClone,
  Tissue = CovidBcells.forheatmap$Tissue
  )

GC = rownames(GetAssayData(CovidBcells.forheatmap, assay = "RNA", slot = "scale.data")) %in% Milpeid2020FrontImmunolgenesets$GC
Mem = rownames(GetAssayData(CovidBcells.forheatmap, assay = "RNA", slot = "scale.data")) %in% Milpeid2020FrontImmunolgenesets$Mem
PB.PC = rownames(GetAssayData(CovidBcells.forheatmap, assay = "RNA", slot = "scale.data")) %in% Milpeid2020FrontImmunolgenesets$PB.PC

pdf("DownsampleLn6_plots/AllCovidBcells.Milpeid2020FrontImmunolgenesets_downsampleLn6.pdf", width = 12, height = 8)
Heatmap(as.matrix(GetAssayData(CovidBcells.forheatmap, assay = "RNA", slot = "scale.data")), name = "COVIDBcells", show_column_names = FALSE,
        cluster_columns = FALSE, column_title_rot = 90,
        column_split = CovidBcells.forheatmap$adt_snn_res.1,
        col = col_fun,
        row_km = 3,
        row_names_gp = gpar(fontsize = 10)
) +
  Heatmap(GC + 0, name = "GC", col = c("0" = "white", "1" = "purple"),
          show_heatmap_legend = FALSE, width = unit(2, "mm"), column_names_gp = gpar(fontsize = 6)) +
  Heatmap(Mem + 0, name = "Mem", col = c("0" = "white", "1" = "red"),
          show_heatmap_legend = FALSE, width = unit(2, "mm"), column_names_gp = gpar(fontsize = 6)) +
  Heatmap(PB.PC + 0, name = "PB/PC", col = c("0" = "white", "1" = "blue"),
          show_heatmap_legend = FALSE, width = unit(2, "mm"), column_names_gp = gpar(fontsize = 6)) +
  rowAnnotation(foo = anno_mark(at = label,
                              labels = label_genes,
                              labels_gp = gpar(fontsize = 14)))
dev.off()

pdf("DownsampleLn6_plots/AllCovidBcells.Milpeid2020FrontImmunolgenesets.clusterColumns_downsampleLn6.pdf", width = 12, height = 8)
Heatmap(as.matrix(GetAssayData(CovidBcells.forheatmap, assay = "RNA", slot = "scale.data")), name = "COVIDBcells", show_column_names = FALSE,
        cluster_columns = TRUE, column_title_rot = 90,
        column_split = CovidBcells.forheatmap$adt_snn_res.1,
        col = col_fun,
        row_km = 3, 
        row_names_gp = gpar(fontsize = 5)
) +
  Heatmap(GC + 0, name = "GC", col = c("0" = "white", "1" = "purple"),
          show_heatmap_legend = FALSE, width = unit(2, "mm"), column_names_gp = gpar(fontsize = 10)) +
  Heatmap(Mem + 0, name = "Mem", col = c("0" = "white", "1" = "red"),
          show_heatmap_legend = FALSE, width = unit(2, "mm"), column_names_gp = gpar(fontsize = 10)) +
  Heatmap(PB.PC + 0, name = "PB/PC", col = c("0" = "white", "1" = "blue"),
          show_heatmap_legend = FALSE, width = unit(2, "mm"), column_names_gp = gpar(fontsize = 10)) +
  rowAnnotation(foo = anno_mark(at = label,
                              labels = label_genes,
                              labels_gp = gpar(fontsize = 10)))
dev.off()
	
###run RNA clustering within the B cells
DefaultAssay(COVIDsampleBcells) <- "RNA"
COVIDsampleBcells <- FindVariableFeatures(COVIDsampleBcells)
COVIDsampleBcells = ScaleData(COVIDsampleBcells, vars.to.regress=c("Donor","Tissue"))
COVIDsampleBcells <- RunPCA(COVIDsampleBcells, assay="RNA", slot = "scale.data", 
                     reduction.name = "pcarna")
COVIDsampleBcells <- FindNeighbors(COVIDsampleBcells, dims = 1:15, reduction = "pcarna")
COVIDsampleBcells <- FindClusters(COVIDsampleBcells, resolution = c(1.15), algorithm = 1)

COVIDsampleBcellsRNAmarkers <- FindAllMarkers(COVIDsampleBcells, max.cells.per.ident = 1000)

write.csv(COVIDsampleBcellsRNAmarkers, "DownsampleLn6_plots/COVIDsampleBcellsRNAClusteringmarkers_downsampleLn6_2.csv")
COVIDsampleBcellsRNAmarkers <- read.csv("DownsampleLn6_plots/COVIDsampleBcellsRNAClusteringmarkers_downsampleLn6_2.csv", header=TRUE)

RNAclusterheatgenes = COVIDsampleBcellsRNAmarkers %>% group_by(cluster) %>% top_n(10, avg_log2FC)

##heatmap for all the COVID sample B cells showing RNA clusters

CovidBcells.forheatmap <- subset(COVIDsampleBcells)
CovidBcells.forheatmap <- ScaleData(CovidBcells.forheatmap, features = RNAclusterheatgenes$gene)

library("circlize")
col_fun = colorRamp2(c(-1.5, 0, 1.5), c("#FF00FF", "#000000", "#FFFF00"))
CITEseqClusterColors <- DiscretePalette(15, palette = "alphabet2")
names(CITEseqClusterColors) <- unique(CovidBcells.forheatmap$adt_snn_res.1)
ha.COVIDBcells = HeatmapAnnotation(
  CITEseqCluster = CovidBcells.forheatmap$adt_snn_res.1,
  Spositive = CovidBcells.forheatmap$Spositve,
  Tissue = CovidBcells.forheatmap$Tissue,
  col = list(CITEseqCluster = CITEseqClusterColors,
  			Spositive = c("FALSE" = "darkseagreen", "TRUE" = "coral4"),
            Tissue = c("Adenoid" = "darkslategray3", "Tonsil" = "darkolivegreen2", "PBMC"="deeppink1")
  			) 
  )

pdf("DownsampleLn6_plots/AllCovidBcells.RNAclusterTop10GenesPerCluster_downsampleLn6.pdf", width = 12, height = 15)
Heatmap(as.matrix(GetAssayData(CovidBcells.forheatmap, assay = "RNA", slot = "scale.data")), name = "COVIDBcells", show_column_names = FALSE,
        top_annotation = ha.COVIDBcells, cluster_columns = FALSE, column_title_rot = 90,
        column_split = CovidBcells.forheatmap$seurat_clusters,
        col = col_fun,
        row_km = 18,
        row_names_gp = gpar(fontsize = 6)
) 
dev.off()


### heatmap of the top RNA for the CITEseq clusters
CovidBcells.forheatmap <- ScaleData(CovidBcells.forheatmap, features = celltypeheatgenes$gene)
col_fun = colorRamp2(c(-1.5, 0, 1.5), c("#FF00FF", "#000000", "#FFFF00"))
RNAClusterColors <- DiscretePalette(19, palette = "alphabet2")
names(RNAClusterColors) <- unique(CovidBcells.forheatmap$RNA_snn_res.1.15)

ha.COVIDBcells = HeatmapAnnotation(
  RNA_Cluster = CovidBcells.forheatmap$RNA_snn_res.1.15,
  Spositive = CovidBcells.forheatmap$Spositve,
  Tissue = CovidBcells.forheatmap$Tissue,
  col = list(RNA_Cluster = RNAClusterColors,
  			Spositive = c("FALSE" = "darkseagreen", "TRUE" = "coral4"),
            Tissue = c("Adenoid" = "darkslategray3", "Tonsil" = "darkolivegreen2", "PBMC"="deeppink1")
 			)
)

pdf("DownsampleLn6_plots//AllCovidBcells.CITEseqclusterTop10GenesPerCluster_downsampleLn6.pdf", width = 12, height = 15)
Heatmap(as.matrix(GetAssayData(CovidBcells.forheatmap, assay = "RNA", slot = "scale.data")), name = "COVIDBcells", show_column_names = FALSE,
        top_annotation = ha.COVIDBcells, cluster_columns = FALSE, column_title_rot = 90,
        column_split = CovidBcells.forheatmap$adt_snn_res.1,
        col = col_fun,
        row_km = 18,
        row_names_gp = gpar(fontsize = 6)
) 
dev.off()

CovidBcells.forheatmap <- ScaleData(CovidBcells.forheatmap, features = as.character(unlist(Milpeid2020FrontImmunolgenesets)))
GC = rownames(GetAssayData(CovidBcells.forheatmap, assay = "RNA", slot = "scale.data")) %in% Milpeid2020FrontImmunolgenesets$GC
Mem = rownames(GetAssayData(CovidBcells.forheatmap, assay = "RNA", slot = "scale.data")) %in% Milpeid2020FrontImmunolgenesets$Mem
PB.PC = rownames(GetAssayData(CovidBcells.forheatmap, assay = "RNA", slot = "scale.data")) %in% Milpeid2020FrontImmunolgenesets$PB.PC
label <- which(rownames(GetAssayData(CovidBcells.forheatmap, assay = "RNA", slot = "scale.data")) %in% as.character(unlist(Milpeid2020FrontImmunolgenesets)))
label_genes <- rownames(GetAssayData(CovidBcells.forheatmap, assay = "RNA", slot = "scale.data"))[label]
	
pdf("DownsampleLn6_plots/AllCovidBcells.Milpeid2020FrontImmunolgenesets.RNAclusters_downsampleLn6.pdf", width = 12, height = 8)
Heatmap(as.matrix(GetAssayData(CovidBcells.forheatmap, assay = "RNA", slot = "scale.data")), name = "COVIDBcells", show_column_names = FALSE,
        top_annotation = ha.COVIDBcells, cluster_columns = FALSE, column_title_rot = 90,
        column_split = CovidBcells.forheatmap$RNA_snn_res.1.15,
        col = col_fun,
        row_km = 3, 
        row_names_gp = gpar(fontsize = 5)
) +
  Heatmap(GC + 0, name = "GC", col = c("0" = "white", "1" = "purple"),
          show_heatmap_legend = FALSE, width = unit(2, "mm"), column_names_gp = gpar(fontsize = 10)) +
  Heatmap(Mem + 0, name = "Mem", col = c("0" = "white", "1" = "red"),
          show_heatmap_legend = FALSE, width = unit(2, "mm"), column_names_gp = gpar(fontsize = 10)) +
  Heatmap(PB.PC + 0, name = "PB/PC", col = c("0" = "white", "1" = "blue"),
          show_heatmap_legend = FALSE, width = unit(2, "mm"), column_names_gp = gpar(fontsize = 10)) +
  rowAnnotation(foo = anno_mark(at = label,
                              labels = label_genes,
                              labels_gp = gpar(fontsize = 10)))
dev.off()
