library(scRepertoire)
library(Seurat)


Human_CD8_CD4_integrated<-readRDS("/Users/xuq6/Documents/Tonsil and adenoid Project/Cihan_TCR_analysis/20221016_Cihan/TCR_scov2_database_integration/Human_CD8_CD4_integrated_cDNA_VDJ.rds")
Idents(Human_CD8_CD4_integrated) <- "orig.ident"
# Extract the CD4 and CD8 data
Human_CD4_integrated <- subset(Human_CD8_CD4_integrated, idents = "CD4_Human")
Human_CD8_integrated <- subset(Human_CD8_CD4_integrated, idents = "CD8_Human")

#########
##########Clonal overlap among CD8 T cells

Human_CD8_integrated@meta.data$cloneType<-Human_CD8_integrated@meta.data$cloneType_Type_Group

Human_CD8_integrated_sc <- as.SingleCellExperiment(Human_CD8_integrated)
newList_Type_Group_CD8 <-expression2List(Human_CD8_integrated_sc, group = "Type_Group")

p1<-clonalOverlap(newList_Type_Group_CD8, cloneCall = "aa", method = "morisita")
p1<-p1 + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+labs(title="Clonotype overlap (Morisita index): TRA-TRB amino acid
")

pdf(file = "oropharyngeal_CD8_plot_Type_Group_clonalOverlap_Morisita_aa_jun19_2022.pdf",width=9, height=5)
print(p1)
dev.off()

###########
###########
###########
###########Clonal overlap among CD4 T cells

Human_CD4_integrated@meta.data$cloneType<-Human_CD4_integrated@meta.data$cloneType_Type_Group

newList_Type_Group_CD4 <-expression2List(Seurat::as.SingleCellExperiment(Human_CD4_integrated), group = "Type_Group")

p2<-clonalOverlap(newList_Type_Group_CD4, cloneCall = "aa", method = "morisita")


p2<-p2 + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+labs(title="Clonotype overlap (Morisita index): TRA-TRB amino acid
")

pdf(file = "oropharyngeal_CD4_plot_Type_Group_clonalOverlap_Morisita_aa_jun19_2022.pdf",width=9, height=5)
print(p1)
dev.off()


