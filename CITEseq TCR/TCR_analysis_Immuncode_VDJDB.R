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
library(ggplot2)
library(reshape)

###
CD8 <- readRDS(file.choose())
meta <- CD8@meta.data
TRBmatch <-if_else(meta$IC_hit_TRB_1 != "Yes" & meta $IC_hit_TRB_2 != "Yes" & meta$VDJDB_hit_TRB_1 != "Yes"& meta$VDJDB_hit_TRB_2 != "Yes", "neg", "pos") 
CD8$TRBmatch <- TRBmatch

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


##Cluster distribution of IC and VDJDB matched TCRb_CDR3aa CD8 T cells among large expanded CD8 cells
###Mapping to the database of ImmuneCODE and VDJDB
CD8_remove_na_large <- subset(CD8_remove_na, cells = Large_expanded_CD8) 
Idents(CD8_remove_na_large) <- "adt_CD8_snn_res.1"
pt3 <- table(Idents(CD8_remove_na_large), CD8_remove_na_large$TRBmatch)%>%as.data.frame()
p3 <- ggplot(pt3, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 12) +
  geom_col(position = "fill", width = 0.5) +
  xlab("TCRb_match") +
  ylab("Proportion") +
  scale_fill_manual(values = c25) +
  theme(legend.title = element_blank())
p3
ggsave("Cluster distrubution of IC and VDJDB matched TCRb_CDR3aa from large expanded CD8 cells.pdf", 
       width = 3,height =4, dpi = 300)

### UMAP of TRB_CDR3aa matched CD8 T cells among large expanded CD8
Idents(CD8_remove_na_large) <- "TRBmatch"
TCRb_matched_CD8_in_large_expanded <- WhichCells(CD8_remove_na_large, idents = c("pos"))
p4 <- DimPlot(CD8_remove_na_large, pt.size = 0.5, reduction = "umap", 
              cells.highlight= list(TCRb_matched_CD8_in_large_expanded),
              cols.highlight = c("dodgerblue2"),
              sizes.highlight = 1,
              label.size = 7,
              na.value = element_blank(),
              group.by = "cloneType_Type_Sample",
              label=FALSE,
              cols = "lightgrey")+NoLegend()
p4+theme(plot.title = element_blank())
ggsave("TRB_CDR3aa matched CD8 T cells among large clonal expanded CD8 T cells.pdf",
       width = 5, height = 5,dpi =300 )

####Frequency of large expanded CD8 T cells among total CD8
pt5 <- prop.table(table(CD8_remove_na$cloneType_Type_Sample,CD8_remove_na$Type_Sample),2)%>%as.data.frame()
pt5[,1] <- factor(pt5[,1],                 # Relevel group factor
                  levels = c("Small_Medium (0 < X <= 0.01)","Large (0.01 < X <= 0.1)"))
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

p5 <- ggplot(pt5,aes(x=factor(Var2),y=Freq)) + 
  geom_col(aes(fill=Var1))
p5+ theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  labs(x = "Sample_type",y = "Frequency of the sample")+
  scale_fill_manual(values = c( "grey","#6A3D9A"))

ggsave("Frequency of large expanded CD8 T cells among total CD8.pdf",
       width = 5,
       height =4,
       dpi=300)


###Frequency of TRB_CDR3aa matched CD8 T cells among large expanded CD8
pt6 <- prop.table(table(CD8_remove_na_large$TRBmatch,CD8_remove_na_large$Type_Sample))%>%as.data.frame()
pt6[,2] <- factor(pt6[,2],                 # Relevel group factor
                  levels = c("CD8_CNMC_71_PBMC", 
                             "CD8_CNMC_71_Adenoid",
                             "CD8_CNMC_71_Tonsil",
                             "CD8_CNMC_89_PBMC",
                             "CD8_CNMC_89_Adenoid",
                             "CD8_CNMC_89_Tonsil",
                             "CD8_CNMC_99_PBMC",
                             "CD8_CNMC_99_Adenoid",
                             "CD8_CNMC_99_Tonsil"))
pt6[,1] <- factor(pt6[,1],                 # Relevel group factor
                  levels = c("pos","neg"))
p6 <- ggplot(pt6,aes(x=factor(Var2),y=Freq)) + 
  geom_col(aes(fill=Var1))
p6+ theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  labs(x = "Sample_type",y = "Frequency of the sample")+
  scale_fill_manual(values = c( "dodgerblue2","grey"))

ggsave("Frequency of TRB_CDR3aa matched cells among total large expaned CD8.pdf",
       width =4,
       height =4)


###visualization in CD4
CD4 <- readRDS(file.choose())

Idents(CD4) <- "cloneType_Type_Sample"
Non_NA_CD4<- WhichCells(CD4, idents = c( "Large (0.01 < X <= 0.1)","Small_Medium (0 < X <= 0.01)"))
CD4_remove_na <- subset(CD4, cells = Non_NA_CD4)  
###Large expanded CD4 among total CD4
pt7 <- prop.table(table(CD4_remove_na$cloneType_Type_Sample,CD4_remove_na$Type_Sample),2)%>%as.data.frame()
pt7[,2] <- factor(pt7[,2],                 # Relevel group factor
                   levels = c("CD4_CNMC_71_PBMC", 
                              "CD4_CNMC_71_Adenoid",
                              "CD4_CNMC_71_Tonsil",
                              "CD4_CNMC_89_PBMC",
                              "CD4_CNMC_89_Adenoid",
                              "CD4_CNMC_89_Tonsil",
                              "CD4_CNMC_99_PBMC",
                              "CD4_CNMC_99_Adenoid",
                              "CD4_CNMC_99_Tonsil"))
pt7[,1] <- factor(pt7[,1],                 # Relevel group factor
                   levels = c(levels = c("Small_Medium (0 < X <= 0.01)"),"Large (0.01 < X <= 0.1)"))

p7 <- ggplot(pt7,aes(x=factor(Var2),y=Freq)) + 
  geom_col(aes(fill=Var1))
p7+ theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  labs(x = "Sample_type",y = "Frequency of the sample")+
  scale_fill_manual(values = c("grey","#6A3D9A"))

ggsave("Frequency of large expanded CD4 among total CD4.pdf",
       width = 5,
       height =4,
       dpi=300)

saveRDS(CD4,"CD4.rds")
saveRDS(CD8,"CD8.rds")
