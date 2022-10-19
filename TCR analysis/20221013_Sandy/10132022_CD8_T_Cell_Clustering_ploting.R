library(Seurat) #load Seurat 4
library(matrixStats)
library(tidyverse)
library(viridis)
library(scico)
library(pheatmap)
library(reshape)
library(stringi)
library(readxl)
library(plyr)


##read CD8.rds
##############################################################
##Clustering with CITEseq surface antibodies of CD8 T cells
adt.dist.T.CD8<- dist(t(GetAssayData(CD8, assay = "CITE2", slot = "data")))
CD8[["adt_CD8_snn"]] <- FindNeighbors(adt.dist.T.CD8, nn.eps = 1)$snn
CD8 <- FindClusters(CD8, resolution = c(1), graph.name = "adt_CD8_snn", algorithm = 1)

##Run UMAP of CD8 T cells
CD8 <- RunUMAP(CD8, assay = "CITE2", features = rownames(CD8[["CITE2"]]), n_neighbors=50L, min_dist=1)
##Get UMAP plot of CD8 T cells
active.assay <- "CITE2"
UMAPciteCD8 <- DimPlot(CD8, pt.size = 0.5, reduction = "umap", group.by = "adt_CD8_snn_res.1", label=TRUE,cols = c25)
UMAPciteCD8
ggsave("adt_CD8_snn.pdf", width = 5, height = 5)

UMAPciteCD8_tissue <- DimPlot(CD8, pt.size = 0.5, reduction = "umap", group.by = "Tissue", order=c("Tonsil","Adenoid","PBMC"),label=TRUE)
UMAPciteCD8_tissue
ggsave("CD8_UMAP_group_by_tissue.pdf", width = 5, height = 5)
# Get metadata
mdCD8 = CD8@meta.data %>% select(adt_CD8_snn_res.1)

# Get data
adtCD8 = GetAssayData(CD8[["CITE2"]]) %>%
  t %>%
  as.data.frame %>%
  rownames_to_column("cell")
mdCD8 <- mdCD8 %>% rownames_to_column("cell")

adtCD8 <- adtCD8 %>%
  select(-c(cell)) %>%
  mutate(cell_type = mdCD8$adt_CD8_snn_res.1) %>%
  select(cell_type, everything())

##source("AverageExpression_MeanOnly.r")
source(file.choose())
environment(AverageExpression_MeanOnly) <- asNamespace('Seurat') # need this to allow the function to call hidden seurat package functions

Idents(CD8) <- "adt_CD8_snn_res.1"
aver_CD8 = AverageExpression_MeanOnly(CD8, return.seurat=T)

quantile_breaks <- function(xs, n = 100) {
  breaks <- quantile(xs, probs = seq(0.5, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

mat_breaks_CD8 <- quantile_breaks(as.matrix(GetAssayData(aver_CD8[["CITE2"]])), n = 101)


plotCD8 <- pheatmap(GetAssayData(aver_CD8[["CITE2"]]), 
                    border_color=NA,
                    inferno(length(mat_breaks_CD8) - 1),
                    kmeans_k = NA,
                    breaks = mat_breaks_CD8,
                    scale = "none", ### CHANGE t() ABOVE TO CONTROL BETTER ###
                    angle_col = 90,
                    cluster_rows = TRUE,
                    cellwidth = 10,
                    cellheight = 10,
                    treeheight_row = 10,
                    treeheight_col = 10,
                    silent = FALSE,
                    main = "Columns: clusters")
ggsave("CD8_heatmap_06072022.pdf", width = 4, height = 3,dpi = 300)
###

meta <- CD8@meta.data
TRBmatch <-if_else(meta$IC_hit_TRB_1 != "Yes" & meta $IC_hit_TRB_2 != "Yes" & meta$VDJDB_hit_TRB_1 != "Yes"& meta$VDJDB_hit_TRB_2 != "Yes", "neg", "pos") 
CD8$TRBmatch <- TRBmatch

cloneType_Type_Sample <- ifelse(meta$TRB_Frequency_Type_Sample>0.01,"Large (0.01 < X <= 0.1)","Small_Medium (0 < X <= 0.01)")
CD8$cloneType_Type_Sample <- cloneType_Type_Sample

Idents(CD8) <- "cloneType_Type_Sample"
##Highlight CD8 cells which are clonal expanded
Large_expanded_CD8 <- WhichCells(CD8, idents = c("Large (0.01 < X <= 0.1)"))
Non_NA_CD8<- WhichCells(CD8, idents = c( "Large (0.01 < X <= 0.1)","Small_Medium (0 < X <= 0.01)"))
CD8_remove_na <- subset(CD8, cells = Non_NA_CD8) 
p1 <- DimPlot(CD8_remove_na, pt.size = 0.5, reduction = "umap", 
              cells.highlight= list(Large_expanded_CD8),
              cols.highlight = "#6A3D9A",
              sizes.highlight = 0.5,
              label.size = 7,
              na.value = element_blank(),
              group.by = "cloneType_Type_Sample",
              label=FALSE,
              cols = "lightgrey")+NoLegend()
p1+theme(plot.title = element_blank())

ggsave("large clonal expanded CD8 T cells among all CD8 T cells.pdf",
       width = 5, height = 5,dpi =300 )

## Define colors
c25 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)

##The distribution of cells in each cluster among expanded vs. unexpanded clonotypes of CD8 T cells
Idents(CD8_remove_na) <- "adt_CD8_snn_res.1"

pt2 <- table(Idents(CD8_remove_na), CD8_remove_na$cloneType_Type_Sample)%>%as.data.frame()
pt2[,2] <- factor(pt2[,2],                 # Relevel group factor
                  levels = c("Small_Medium (0 < X <= 0.01)","Large (0.01 < X <= 0.1)"))
p2 <- ggplot(pt2, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 12) +
  geom_col(position = "fill", width = 0.5) +
  xlab("ColoneType") +
  ylab("Proportion") +
  scale_fill_manual(values = c25) +
  theme(legend.title = element_blank(),
        axis.text.x = element_text(angle = 5, vjust = 0.5, hjust=0.5))

p2
ggsave("Cluster distribution of expanded and unexpanded CD8 T cells.pdf", 
       width = 3,height =4, dpi = 300)

CD8_remove_na_large <- subset(CD8_remove_na, cells = Large_expanded_CD8)


### UMAP of TRB_CDR3aa matched CD8 T cells among large expanded CD8
Idents(CD8_remove_na_large) <- "TRBmatch"
TCRb_matched_CD8_in_large_expanded <- WhichCells(CD8_remove_na_large, idents = c("pos"))
p3 <- DimPlot(CD8_remove_na_large, pt.size = 1, reduction = "umap", 
              cells.highlight= list(TCRb_matched_CD8_in_large_expanded),
              cols.highlight = c("dodgerblue2"),
              sizes.highlight = 1,
              label.size = 7,
              na.value = element_blank(),
              group.by = "cloneType_Type_Sample",
              label=FALSE,
              cols = "lightgrey")+NoLegend()
p3+theme(plot.title = element_blank())
ggsave("TRB_CDR3aa matched CD8 T cells among large clonal expanded CD8 T cells.pdf",
       width = 5, height = 5,dpi =300 )

####Frequency of large expanded CD8 T cells among total CD8
pt4 <- prop.table(table(CD8_remove_na$cloneType_Type_Sample,CD8_remove_na$Type_Sample),2)%>%as.data.frame()
pt4[,1] <- factor(pt4[,1],                 # Relevel group factor
                  levels = c("Small_Medium (0 < X <= 0.01)","Large (0.01 < X <= 0.1)"))
pt4[,2] <- factor(pt4[,2],                 # Relevel group factor
                  levels = c("CD8_CNMC_71_PBMC", 
                             "CD8_CNMC_71_Adenoid",
                             "CD8_CNMC_71_Tonsil",
                             "CD8_CNMC_89_PBMC",
                             "CD8_CNMC_89_Adenoid",
                             "CD8_CNMC_89_Tonsil",
                             "CD8_CNMC_99_PBMC",
                             "CD8_CNMC_99_Adenoid",
                             "CD8_CNMC_99_Tonsil"))

p4 <- ggplot(pt4,aes(x=factor(Var2),y=Freq)) + 
  geom_col(aes(fill=Var1))
p4+ theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  labs(x = "Sample_type",y = "Frequency of the sample")+
  scale_fill_manual(values = c( "grey","#6A3D9A"))

ggsave("Frequency of large expanded CD8 T cells among total CD8.pdf",
       width = 5,
       height =4,
       dpi=300)


###Frequency of TRB_CDR3aa matched CD8 T cells among large expanded CD8
pt5 <- prop.table(table(CD8_remove_na_large$TRBmatch,CD8_remove_na_large$Type_Sample))%>%as.data.frame()
pt5[,2] <- factor(pt5[,2],                 # Relevel group factor
                  levels = c("CD8_CNMC_71_PBMC", 
                             "CD8_CNMC_71_Adenoid",
                             "CD8_CNMC_71_Tonsil",
                             "CD8_CNMC_89_PBMC",
                             "CD8_CNMC_89_Adenoid",
                             "CD8_CNMC_89_Tonsil",
                             "CD8_CNMC_99_PBMC",
                             "CD8_CNMC_99_Adenoid",
                             "CD8_CNMC_99_Tonsil"))
pt5[,1] <- factor(pt5[,1],                 # Relevel group factor
                  levels = c("pos","neg"))
p5 <- ggplot(pt5,aes(x=factor(Var2),y=Freq)) + 
  geom_col(aes(fill=Var1))
p5+ theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  labs(x = "Sample_type",y = "Frequency of the sample")+
  scale_fill_manual(values = c( "dodgerblue2","grey"))

ggsave("Frequency of TRB_CDR3aa matched cells among total large expaned CD8.pdf",
       width =4,
       height =4)


#