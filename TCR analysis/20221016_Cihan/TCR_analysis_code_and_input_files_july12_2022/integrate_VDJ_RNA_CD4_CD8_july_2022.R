setwd("/Users/oguzc/Downloads/CITEseq_TCR/CITEseq_TCR_code_for_github")

library(scRepertoire)
library(Seurat)

data("contig_list")
cols_example<-as.matrix(colnames(contig_list[[1]]))

immdata_10x_VDJ<-list()
colclass_example<-sapply(contig_list[[1]], class)

immdata_10x_VDJ[[1]]<-read.csv("filtered_contig_annotations_09.csv",colClasses=colclass_example)

ind_match<-which((colnames(immdata_10x_VDJ[[1]]) %in% colnames(contig_list[[1]]))==TRUE)
immdata_10x_VDJ[[1]]<-immdata_10x_VDJ[[1]][,c(ind_match)]


immdata_10x_VDJ[[2]]<-read.csv("filtered_contig_annotations_10.csv",colClasses=colclass_example)
ind_match<-which((colnames(immdata_10x_VDJ[[2]]) %in% colnames(contig_list[[2]]))==TRUE)
immdata_10x_VDJ[[2]]<-immdata_10x_VDJ[[2]][,c(ind_match)]

immdata_10x_VDJ[[3]]<-read.csv("filtered_contig_annotations_11.csv",colClasses=colclass_example)
ind_match<-which((colnames(immdata_10x_VDJ[[3]]) %in% colnames(contig_list[[3]]))==TRUE)
immdata_10x_VDJ[[3]]<-immdata_10x_VDJ[[3]][,c(ind_match)]


immdata_10x_VDJ[[4]]<-read.csv("filtered_contig_annotations_12.csv",colClasses=colclass_example)
ind_match<-which((colnames(immdata_10x_VDJ[[4]]) %in% colnames(contig_list[[4]]))==TRUE)
immdata_10x_VDJ[[4]]<-immdata_10x_VDJ[[4]][,c(ind_match)]

#09-10 for CD4 and 11-12 for CD8

combinedVDJ <- combineTCR(immdata_10x_VDJ, samples = c("CD4Tcells","CD4Tcells","CD8Tcells","CD8Tcells"), ID = c("09","10","11","12"), cells = "T-AB")

##########
##########

### "CD4Tcells.rds" and "CD8Tcells.rds" are subseted from the ObjList "ObjList_filt.preprocessed.ADTclustered.wMeta.NoUS.Rds".
CD4_sc<-readRDS("CD4Tcells.rds")

CD8_sc<-readRDS("CD8Tcells.rds")


ind_CD4_Singlet<-which(CD4_sc@meta.data$HTO_classification.global=="Singlet")


#CNMC 71 and 89 are COVID-convalescent and CNMC 99 is a control
#There are tonsil cells, adenoid cells, and PBMCs from each donor

HTO_sample_mapping<-read.table("HTO_sample_mapping_tonsils.txt",header = T)
ind_CD4_Singlet<-which(CD4_sc@meta.data$HTO_classification.global=="Singlet")


Idents(CD4_sc) <- "HTO_classification.global"
# Extract the singlets
Human_CD4_sc <- subset(CD4_sc, idents = "Singlet")
Human_CD4_sc[["percent.mt"]] <- PercentageFeatureSet(Human_CD4_sc, pattern = "^MT-")
Idents(CD8_sc) <- "HTO_classification.global"
# Extract the singlets
Human_CD8_sc <- subset(CD8_sc, idents = "Singlet")
Human_CD8_sc[["percent.mt"]] <- PercentageFeatureSet(Human_CD8_sc, pattern = "^MT-")

Human_CD4_sc@meta.data[["orig.ident"]]<-"Human_CD4_sc"
Human_CD8_sc@meta.data[["orig.ident"]]<-"Human_CD8_sc"


#######
#######
Human_CD4_sc@meta.data[["orig.ident"]]<-paste0("CD4_",Human_CD4_sc@meta.data[["HTO_maxID"]])

Human_CD8_sc@meta.data[["orig.ident"]]<-paste0("CD8_",Human_CD8_sc@meta.data[["HTO_maxID"]])

Human_CD4_sc@meta.data$Sample<-"Sample"

for (rowno in 1:nrow(HTO_sample_mapping))
{
#######
#######
ind_selec<-NULL
#######
#######
ind_selec<-which(Human_CD4_sc@meta.data$HTO_maxID==HTO_sample_mapping$HTO_number[rowno])
#######
#######
Human_CD4_sc@meta.data$Sample[ind_selec]<-HTO_sample_mapping$Sample[rowno]
#######
#######
}


Human_CD8_sc@meta.data$Sample<-"Sample"

for (rowno in 1:nrow(HTO_sample_mapping))
{
#######
#######
ind_selec<-NULL
#######
#######
ind_selec<-which(Human_CD8_sc@meta.data$HTO_maxID==HTO_sample_mapping$HTO_number[rowno])
#######
#######
Human_CD8_sc@meta.data$Sample[ind_selec]<-HTO_sample_mapping$Sample[rowno]
#######
#######
}


features <- SelectIntegrationFeatures(object.list = list(Human_CD4_sc,Human_CD8_sc))

Human_CD4_sc <- ScaleData(Human_CD4_sc, features = features, verbose = TRUE)
Human_CD8_sc <- ScaleData(Human_CD8_sc, features = features, verbose = TRUE)

Human_CD4_sc <- RunPCA(Human_CD4_sc, features = features, verbose = TRUE)
Human_CD8_sc <- RunPCA(Human_CD8_sc, features = features, verbose = TRUE)


Human_CD8_sc_CD4_sc.anchors <- FindIntegrationAnchors(object.list = list(Human_CD4_sc,Human_CD8_sc), dims = 1:20,verbose=TRUE, anchor.features = features, reduction = "rpca",k.anchor = 40, reference = c(1, 2))

Human_CD8_CD4_integrated <- IntegrateData(anchorset = Human_CD8_sc_CD4_sc.anchors, dims = 1:20)

DefaultAssay(Human_CD8_CD4_integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
Human_CD8_CD4_integrated <- ScaleData(Human_CD8_CD4_integrated, verbose = TRUE)

Human_CD8_CD4_integrated <- RunPCA(Human_CD8_CD4_integrated, npcs = 20, verbose = FALSE)
# UMAP and Clustering
Human_CD8_CD4_integrated <- RunUMAP(Human_CD8_CD4_integrated, reduction = "pca", dims = 1:20)
Human_CD8_CD4_integrated <- FindNeighbors(Human_CD8_CD4_integrated, reduction = "pca", dims = 1:20)
Human_CD8_CD4_integrated <- FindClusters(Human_CD8_CD4_integrated, resolution = 0.3)

####################

#Human_CD8_CD4_integrated@meta.data$Sample<-gsub("Adneoid","Adenoid",Human_CD8_CD4_integrated@meta.data$Sample)


#####################
ind_CD8<-grep("CD8",Human_CD8_CD4_integrated@meta.data[["Type"]])
ind_CD4<-grep("CD4",Human_CD8_CD4_integrated@meta.data[["Type"]])

Human_CD8_CD4_integrated$Sample_type<-Human_CD8_CD4_integrated$orig.ident

Human_CD8_CD4_integrated$orig.ident[ind_CD8]<-"CD8_Human"
Human_CD8_CD4_integrated$orig.ident[ind_CD4]<-"CD4_Human"


cell.names <- rownames(Human_CD8_CD4_integrated@meta.data)


#(09-10 for CD4 and 11-12 for CD8)

ind_group1<-grep('CD4', Human_CD8_CD4_integrated@meta.data[["orig.ident"]])
ind_group2<-grep('CD8', Human_CD8_CD4_integrated@meta.data[["orig.ident"]])

cell.names[ind_group1] <- paste0("CD4_Human_", cell.names[ind_group1])
cell.names[ind_group2] <- paste0("CD8_Human_", cell.names[ind_group2])


cell.names.matrix<-as.matrix(cell.names)

cell.names.matrix[grep("09",cell.names.matrix)]<-gsub("CD4_Human_","CD4Tcells_09_",cell.names.matrix[grep("09",cell.names.matrix)])
cell.names.matrix[grep("10",cell.names.matrix)]<-gsub("CD4_Human_","CD4Tcells_10_",cell.names.matrix[grep("10",cell.names.matrix)])
cell.names.matrix[grep("11",cell.names.matrix)]<-gsub("CD8_Human_","CD8Tcells_11_",cell.names.matrix[grep("11",cell.names.matrix)])
cell.names.matrix[grep("12",cell.names.matrix)]<-gsub("CD8_Human_","CD8Tcells_12_",cell.names.matrix[grep("12",cell.names.matrix)])


cell.names.matrix[grep("-09",cell.names.matrix)]<-gsub("-09","",cell.names.matrix[grep("-09",cell.names.matrix)])
cell.names.matrix[grep("-10",cell.names.matrix)]<-gsub("-10","",cell.names.matrix[grep("-10",cell.names.matrix)])
cell.names.matrix[grep("-11",cell.names.matrix)]<-gsub("-11","",cell.names.matrix[grep("-11",cell.names.matrix)])
cell.names.matrix[grep("-12",cell.names.matrix)]<-gsub("-12","",cell.names.matrix[grep("-12",cell.names.matrix)])


Human_CD8_CD4_integrated <- RenameCells(Human_CD8_CD4_integrated, new.names = cell.names.matrix)

#########
#########
#########
#########

combinedVDJ[[1]][["barcode"]]<-gsub("-1","",combinedVDJ[[1]][["barcode"]])
combinedVDJ[[2]][["barcode"]]<-gsub("-1","",combinedVDJ[[2]][["barcode"]])
combinedVDJ[[3]][["barcode"]]<-gsub("-1","",combinedVDJ[[3]][["barcode"]])
combinedVDJ[[4]][["barcode"]]<-gsub("-1","",combinedVDJ[[4]][["barcode"]])


#########
#########
#########

Human_CD8_CD4_integrated <- combineExpression(combinedVDJ, Human_CD8_CD4_integrated,cloneCall="aa", group.by = "none")

#########
#########

Human_CD8_CD4_integrated@meta.data$Type_Sample<-paste0(Human_CD8_CD4_integrated@meta.data$Type,"_",Human_CD8_CD4_integrated@meta.data$Sample)

Human_CD8_CD4_integrated@meta.data$Type_Sample<-gsub("Tcells","",Human_CD8_CD4_integrated@meta.data$Type_Sample)

#########
#########
#CNMC 71 and 89 are COVID-convalescent and CNMC 99 is a control; there are tonsil cells, adenoid cells, and PBMCs from each donor


###Do case-control based splitting here
Human_CD8_CD4_integrated@meta.data$Type_Group<-Human_CD8_CD4_integrated@meta.data$Type_Sample

Human_CD8_CD4_integrated@meta.data$Type_Group<-gsub("CNMC_99","Control",Human_CD8_CD4_integrated@meta.data$Type_Group)

Human_CD8_CD4_integrated@meta.data$Type_Group<-gsub("CNMC_89","Convalescent",Human_CD8_CD4_integrated@meta.data$Type_Group)

Human_CD8_CD4_integrated@meta.data$Type_Group<-gsub("CNMC_71","Convalescent",Human_CD8_CD4_integrated@meta.data$Type_Group)


saveRDS(Human_CD8_CD4_integrated, "Human_CD8_CD4_integrated_cDNA_VDJ.rds")
