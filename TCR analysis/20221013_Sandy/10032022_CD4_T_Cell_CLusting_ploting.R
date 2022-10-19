library(Seurat) #load Seurat 4
library(matrixStats)
library(tidyverse)
library(viridis)
library(scico)
library(pheatmap)
library(reshape)
library(stringi)
library(readxl)
library(plyr)###Remove the package of plyr when using dplyr.It causes problems of using dplyr.

#######Read the T cell Seurat Object 
##Tcell_path <-"/Users/xuq6/Documents/Tonsil and adenoid Project/Cihan_TCR_analysis/20220628_Cihan/Human_CD8_CD4_integrated_cDNA_VDJ_jun28_2022.rds"
Tcell_path <- "/Users/xuq6/Documents/Cousera/R_code_representation/Code_All/Cihan/20221016_Cihan/TCR_analysis_code_and_input_files_july12_2022/Human_CD8_CD4_integrated_cDNA_VDJ.rds"
Tcell <- readRDS(Tcell_path)
#######Read the new metadata
Meta_new_path <- "/Users/xuq6/Documents/Cousera/R_code_representation/Code_All/Cihan/20221016_Cihan/TCR_scov2_database_integration/TCR_VDJDB_TRA_TRB_metadata_abundance_pediatric_samples_oct16_2022.xlsx"
Meta_new <- read_excel(path = Meta_new_path,sheet =1)

#######Remove unused columns
Tcell_meta <- Tcell@meta.data
Tcell_meta$cell_barcode <- rownames(Tcell_meta)

##Check if these two tables are with the same order of the cell_barcode
all(Meta_new$cell_barcode == Tcell_meta$cell_barcode)
##Returned FAlSE

###Join new metadata to the Seurat Object
meta_update <- plyr::join(Tcell_meta,Meta_new,by="cell_barcode",type = "left")
names(meta_update)

##Split TRA into TRA1 and TRA2
TRA1_2<- str_split_fixed(meta_update$TRA, ';', n=2)%>%as.data.frame()
colnames(TRA1_2) <- c("TRA1","TRA2")
meta_update$TRA1 <- TRA1_2$TRA1
meta_update$TRA2 <- TRA1_2$TRA2
## Add tissue 
tissue <- str_split_fixed(Tcell@meta.data$Sample, '_', n=3)%>%as.data.frame()
meta_update$Tissue<- tissue$V3 


##Cell name should be included before AddMetaData
rownames(meta_update) <- meta_update$cell_barcode
names(Tcell@meta.data)
all(meta_update$cell_barcode == Tcell_meta$cell_barcode)
Tcell <- AddMetaData(Tcell, metadata=meta_update)

head(Tcell@meta.data$VDJDB_hit_TRB_1)
head(Meta_new$VDJDB_hit_TRB_1)

####The staining of CCR6(CD196) and CCR7(CD197) is not good. Remove them from CITE list ##"CD197","CD196"
####Remove the B cell markers
CITEmarkers <- c("CD38","CD27","CD62L" ,"CD69","HLA-DR","CD279","CD57","CD103","CD183","CD185")
SubsetCITE <- Tcell@assays$CITE[CITEmarkers,colnames(Tcell)]
Tcell[["CITE2"]] <- CreateAssayObject(data = as.matrix(SubsetCITE))

#######################################
##Subset CD4 and CD8
Idents(Tcell) <- "Type"
CD4 <- subset(Tcell,idents = "CD4Tcells")
CD8 <- subset(Tcell,idents = "CD8Tcells")

##Clustering with new CITE2 surface antibodies within CD4
adt.dist.T.CD4<- dist(t(GetAssayData(CD4, assay = "CITE2", slot = "data")))
CD4[["adt_CD4_snn"]] <- FindNeighbors(adt.dist.T.CD4, nn.eps = 1)$snn
CD4 <- FindClusters(CD4, resolution = c(1), graph.name = "adt_CD4_snn", algorithm = 1)

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
##Run UMAP for CD4 T cells
CD4 <- RunUMAP(CD4, assay = "CITE2", features = rownames(CD4[["CITE2"]]), n_neighbors=50L, min_dist=1)
active.assay <- "CITE2"
UMAPciteCD4 <- DimPlot(CD4, pt.size = 0.5, reduction = "umap", group.by = "adt_CD4_snn_res.1", label=TRUE, cols=c25)
UMAPciteCD4
ggsave("adt_CD4_snn.pdf", width = 5, height = 5)

UMAPciteCD4_tissue <- DimPlot(CD4, pt.size = 0.5, reduction = "umap", group.by = "Tissue", order=c("Tonsil","Adenoid","PBMC"),label=TRUE)
UMAPciteCD4_tissue
ggsave("CD4_UMAP_group_by_tissue.pdf", width = 5, height = 5)

##Get metadata
mdCD4 = CD4@meta.data %>% select(adt_CD4_snn_res.1)

##Get data
adtCD4 = GetAssayData(CD4[["CITE2"]]) %>%
  t %>%
  as.data.frame %>%
  rownames_to_column("cell")
mdCD4 <- mdCD4 %>% rownames_to_column("cell")

adtCD4 <- adtCD4 %>%
  select(-c(cell)) %>%
  mutate(cell_type = mdCD4$adt_CD4_snn_res.1) %>%
  select(cell_type, everything())

source(file.choose())
environment(AverageExpression_MeanOnly) <- asNamespace('Seurat') # need this to allow the function to call hidden Seurat package functions

Idents(CD4) <- "adt_CD4_snn_res.1"
aver_CD4 = AverageExpression_MeanOnly(CD4, return.seurat=T)

quantile_breaks <- function(xs, n = 100) {
  breaks <- quantile(xs, probs = seq(0.1, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

mat_breaks_CD4 <- quantile_breaks(as.matrix(GetAssayData(aver_CD4[["CITE2"]])), n = 101)

##Get heatmap of surface protein antibodies of the CD4 T cells
plotCD4 <- pheatmap(GetAssayData(aver_CD4[["CITE2"]]), 
                     border_color=NA,
                     inferno(length(mat_breaks_CD4) - 1),
                     kmeans_k = NA,
                     breaks = mat_breaks_CD4,
                     scale = "none", ### CHANGE t() ABOVE TO CONTROL BETTER ###
                     angle_col = 90,
                     cluster_rows = TRUE,
                     cellwidth = 10,
                     cellheight = 10,
                     treeheight_row = 10,
                     treeheight_col = 10,
                     silent = FALSE,
                     main = "Columns: clusters")
ggsave("CD4_heatmap_08242022.pdf",
       width = 4,
       height = 3,
       dpi = 300)


#######Remove TRB with na
Idents(CD4) <- "TRB_Abundance_Type_Sample"
table(is.na(CD4$TRB_Abundance_Type_Sample))
Non_NA_CD4_TRB <- filter(CD4@meta.data,TRB_Abundance_Type_Sample!="NA")
CD4_remove_na_TRB <- subset(CD4, cells = rownames(Non_NA_CD4_TRB))
CD4_TRBmatch <-ifelse(CD4_remove_na_TRB$VDJDB_hit_TRB_1 == "Yes" |CD4_remove_na_TRB$VDJDB_hit_TRB_2== "Yes"|CD4_remove_na_TRB$IC_hit_TRB_1 == "Yes"|CD4_remove_na_TRB$IC_hit_TRB_2 == "Yes", "pos", "neg") 
CD4_remove_na_TRB$CD4_TRBmatch  <- CD4_TRBmatch 

###TRB analysis: Expanded CD4 with medium and larger expansion, and number of cells within one clone more than 3 from an individual sample.
Expanded_Type_Sample_CD4_TRB <- ifelse(CD4_remove_na_TRB@meta.data[["TRB_Frequency_Type_Sample"]]>0.001&CD4_remove_na_TRB@meta.data[["TRB_Abundance_Type_Sample"]]>2,"Medium_large_expansion_TRB","Non_expanded_TRB")
CD4_remove_na_TRB$Expanded_Type_Sample_CD4_TRB <- Expanded_Type_Sample_CD4_TRB
Medium_large_TRB_CD4 <- filter(CD4_remove_na_TRB@meta.data,Expanded_Type_Sample_CD4_TRB == "Medium_large_expansion_TRB")

#Idents(CD4_remove_na_TRB) <- "Expanded_Type_Sample_CD4_TRB"
Medium_large_TRB_CD4_cells <- subset(CD4_remove_na_TRB, cells = rownames(Medium_large_TRB_CD4))

###Medium_Large expanded CD4 among total non NA CD4
pt1 <- prop.table(table(CD4_remove_na_TRB$Expanded_Type_Sample_CD4_TRB,CD4_remove_na_TRB$Type_Sample),2)%>%as.data.frame()
pt1[,2] <- factor(pt1[,2],                 # Relevel group factor
                   levels = c("CD4_CNMC_71_PBMC", 
                              "CD4_CNMC_71_Adenoid",
                              "CD4_CNMC_71_Tonsil",
                              "CD4_CNMC_89_PBMC",
                              "CD4_CNMC_89_Adenoid",
                              "CD4_CNMC_89_Tonsil",
                              "CD4_CNMC_99_PBMC",
                              "CD4_CNMC_99_Adenoid",
                              "CD4_CNMC_99_Tonsil"))
pt1[,1] <- factor(pt1[,1],                 # Relevel group factor
                   levels = c("Non_expanded_TRB","Medium_large_expansion_TRB"))

p1 <- ggplot(pt1,aes(x=factor(Var2),y=Freq)) + 
  geom_col(aes(fill=Var1))
p1+ theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  labs(x = "Sample_type",y = "Frequency of the sample")+
  scale_fill_manual(values = c("grey","#6A3D9A"))

ggsave("Frequency of Medium large expanded CD4 among non NA CD4 TRB.pdf",
       width = 5,
       height =4,
       dpi=300)



###Read table matching to the ref papers
Ref_path <- "/Users/xuq6/Documents/Cousera/R_code_representation/Code_All/Sandy/Update_10132022/CD4_new_meta_ref.csv"
Ref<- read.csv(file = Ref_path)
###TRAmatch for CD4
table(is.na(CD4$TRA_Frequency_Type_Sample))
CD4_meta <- CD4@meta.data

Non_NA_CD4_TRA <- filter(CD4_meta, TRA_Abundance_Type_Sample!="NA")

###Subset CD4 
Idents(CD4) <- "TRA_Abundance_Type_Sample"
CD4_Non_na_TRA <- subset(CD4, cells =rownames(Non_NA_CD4_TRA))

CD4_TRA_meta<- CD4_Non_na_TRA@meta.data
write.csv(CD4_TRA_meta,"CD4_remove_na_TRA_meta.csv")
rownames(Ref) <- Ref$cell_barcode
Ref <- Ref[,c(35,71:74)]

all(Ref$cell_barcode == CD4_TRA_meta$cell_barcode)

CD4_meta_update <- plyr::join(CD4_TRA_meta,Ref,by="cell_barcode",type = "left")
rownames(CD4_meta_update) <- CD4_meta_update$cell_barcode


##Add Ref to CD4
CD4_Non_na_TRA <- AddMetaData(CD4_Non_na_TRA,metadata = CD4_meta_update)
CD4_meta_na_TRA <- CD4_Non_na_TRA@meta.data
CD4_TRAmatch <-ifelse(CD4_meta_na_TRA$VDJDB_hit_TRA == "Yes" |CD4_meta_na_TRA$ID_all_Ref_matched_TRAaa == "pos", "pos", "neg") 
CD4_Non_na_TRA$CD4_TRAmatch  <- CD4_TRAmatch 

###TRA analysis: Expanded CD4 with medium and larger expansion, and number of cells within one clone more than 3 from an individual sample.
#Idents(CD4_Non_na_TRA) <- "TRA_Frequency_Type_Sample"
Expanded_Type_Sample_CD4_TRA <- ifelse(CD4_Non_na_TRA@meta.data[["TRA_Frequency_Type_Sample"]]>0.001&CD4_Non_na_TRA@meta.data[["TRA_Abundance_Type_Sample"]]>2,"Medium_large_expansion_TRA","Non_expanded_TRA")
CD4_Non_na_TRA$Expanded_Type_Sample_CD4_TRA <- Expanded_Type_Sample_CD4_TRA

Medium_large_TRA_CD4 <- filter(CD4_Non_na_TRA@meta.data,Expanded_Type_Sample_CD4_TRA=="Medium_large_expansion_TRA")
Medium_large_TRA_CD4_cells <- subset(CD4_Non_na_TRA, cells = rownames(Medium_large_TRA_CD4))

Idents(CD4_Non_na_TRA) <- "Expanded_Type_Sample_CD4_TRA"
CD4_medium_large_expanded_TRA <- WhichCells(CD4_Non_na_TRA, idents = c("Medium_large_expansion_TRA"))


###TRA Medium_Large expanded CD4 among total non NA CD4
pt2 <- prop.table(table(CD4_Non_na_TRA$Expanded_Type_Sample_CD4_TRA,CD4_Non_na_TRA$Type_Sample),2)%>%as.data.frame()
pt2[,2] <- factor(pt2[,2],                 # Relevel group factor
                  levels = c("CD4_CNMC_71_PBMC", 
                             "CD4_CNMC_71_Adenoid",
                             "CD4_CNMC_71_Tonsil",
                             "CD4_CNMC_89_PBMC",
                             "CD4_CNMC_89_Adenoid",
                             "CD4_CNMC_89_Tonsil",
                             "CD4_CNMC_99_PBMC",
                             "CD4_CNMC_99_Adenoid",
                             "CD4_CNMC_99_Tonsil"))
pt2[,1] <- factor(pt2[,1],                 # Relevel group factor
                  levels = c("Non_expanded_TRA","Medium_large_expansion_TRA"))

p2 <- ggplot(pt2,aes(x=factor(Var2),y=Freq)) + 
  geom_col(aes(fill=Var1))
p2+ theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  labs(x = "Sample_type",y = "Frequency of the sample")+
  scale_fill_manual(values = c("grey","#6A3D9A"))

ggsave("Frequency of Medium large expanded CD4 among non NA CD4 TRA.pdf",
       width = 5,
       height =4,
       dpi=300)




TRA_pos <- filter(Medium_large_TRA_CD4_cells@meta.data,CD4_TRAmatch=="pos")
table(TRA_pos$Type_Sample,TRA_pos$TRA_Abundance_Type_Sample)


##TRA: The distribution of cells in each cluster among TRB matched vs. unmatched of expanded CD4 T cells
Idents(Medium_large_TRA_CD4_cells) <- "adt_CD4_snn_res.1"

## Define colors, no expanded cells in cluster 7
c13 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", 
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)



table <- table(CD4_Non_na_TRA$TRA1)%>%as.data.frame()
table <- filter(table,Freq>2)
TRA_among_all_table <- filter(Medium_large_TRA_CD4_cells@meta.data,TRA %in% table$Var1)
write.csv(TRA_among_all_table, "TRA_among_all_table.csv")

###Merge TRA and TRB medium and expanded together
Idents(CD4) <- "adt_CD4_snn_res.1"

Shared_list <-filter(Medium_large_TRA_CD4_cells@meta.data,cell_barcode%in%Medium_large_TRB_CD4_cells@meta.data$cell_barcode)
CD4_medium_large_expanded_TRA_minus_shared <- Medium_large_TRA_CD4_cells@meta.data[!rownames(Medium_large_TRA_CD4_cells@meta.data) %in% rownames(Shared_list),]
CD4_medium_large_expanded_TRB_minus_shared <-Medium_large_TRB_CD4_cells@meta.data[!rownames(Medium_large_TRB_CD4_cells@meta.data) %in% rownames(Shared_list),]


p3 <- DimPlot(CD4, pt.size = 0.5, reduction = "umap", 
              cells.highlight= list(rownames(CD4_medium_large_expanded_TRA_minus_shared),
                                    rownames(CD4_medium_large_expanded_TRB_minus_shared),
                                    rownames(Shared_list)),
              cols.highlight = c("Turquoise4","#6A3D9A","red2"),
              sizes.highlight = 1,
              label.size = 7,
              na.value = element_blank(),
              label=FALSE,
              cols = "lightgrey")
p3+theme(plot.title = element_blank())
ggsave("TRA and TRB expanded CD4 T cells.pdf", 
       width = 5,height =5, dpi = 300)


######Distribution of expanded TRA and TRB among CD4
####read table
pt4 <- read_excel(path = "/Users/xuq6/Documents/Tonsil and adenoid Project/Cihan_TCR_analysis/20220819_updated_analysis/TRA_TRB_expansion.xlsx",sheet =1)
pt4 <- as.data.frame(pt4)
pt4[,2] <- factor(pt4[,2],                 # Relevel group factor
                  levels = c("Total_CD4","Medium_large_expansion_TRA","Medium_large_expansion_TRB"))
pt4[,1] <- factor(pt4[,1])
p4 <- ggplot(pt4, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 12) +
  geom_col(position = "fill", width = 0.5) +
  xlab("ColoneType") +
  ylab("Proportion") +
  scale_fill_manual(values = c25) +
  theme(legend.title = element_blank(),
        axis.text.x = element_text(angle = 5, vjust = 0.5, hjust=0.5))


p4
ggsave("CD4 expanded TRA and TRB cluster distribution.pdf", 
       width = 4,height =4, dpi = 300)



saveRDS(CD4, "CD4.rds")
saveRDS(CD8, "CD8.rds")