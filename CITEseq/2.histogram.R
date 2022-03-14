#this is tested using Cicero conda environment (with R 4 and Seurat 4) on a high-performance computing node with 16 cores and at least 160 gb or ram. 
library("Seurat") #load Seurat 4
library("matrixStats")
library('tidyverse')
library('magrittr')
library('DescTools')
library(ggridges)
library(viridis)
library(scico)
library('dplyr')
library('pheatmap')
library('gplots')

ObjList_filt <- readRDS("SeuratObjects/ObjList_filt.preprocessed.ADTclustered.wMeta.Rds")

i = 1
	# get metadata
	md = ObjList_filt[[i]]@meta.data %>% select(adt_snn_res.1)

	# get data
	adt = GetAssayData(ObjList_filt[[i]][["CITE"]]) %>%
  	t %>%
  	as.data.frame %>%
  	rownames_to_column("cell")
	md <- md %>% rownames_to_column("cell")
	
	adt <- adt %>%
  		select(-c(cell)) %>%
  		mutate(cell_type = md$adt_snn_res.1) %>%
  		select(cell_type, everything())

	source("AverageExpression_MeanOnly.r")
	environment(AverageExpression_MeanOnly) <- asNamespace('Seurat') # need this to allow the function to call hidden seurat package functions
	
	Idents(ObjList_filt[[i]]) <- "adt_snn_res.1"
	aver = AverageExpression_MeanOnly(ObjList_filt[[i]], return.seurat=T)

	quantile_breaks <- function(xs, n = 100) {
  		breaks <- quantile(xs, probs = seq(0.5, 1, length.out = n))
  		breaks[!duplicated(breaks)]
		}

	mat_breaks <- quantile_breaks(as.matrix(GetAssayData(aver[["CITE"]])), n = 101)
	library("viridis")
	ggsave(filename = paste0("ClusterVis/", names(ObjList_filt[i]), ".bycluster_infernoColor_quantileBreaks.pdf"), width = 8, height = 14,
        plot = pheatmap(GetAssayData(aver[["CITE"]]), scale = "none", border_color=NA, inferno(length(mat_breaks) - 1), breaks = mat_breaks))

	# get hclust order of proteins and cell types
    x = pheatmap::pheatmap(GetAssayData(aver[["CITE"]]),
             cluster_rows = T, cluster_cols = T,
             fontsize_col = 10, fontsize_row = 8, border_color = NA)

    prot_order = rownames(GetAssayData(aver[["CITE"]])[x$tree_row$order, ])
    celltype_order = colnames(GetAssayData(aver[["CITE"]])[,x$tree_col$order ])

	library('genefilter')
	#filter out proteins that do not have at least 1 cluster with an average over 2
    	f1 <- kOverA(1, 2)
    	ffun <- filterfun(f1)
    	filtaver.2 <- genefilter(GetAssayData(aver[["CITE"]])[x$tree_row$order, ], ffun)

   		prot_plot = prot_order[which(filtaver.2)]

	adt$cell_type <- as.character(adt$cell_type)
	adt.l = adt %>%
  		gather(key = prot, value = normalized_count, CD86:CD197) %>%
  		filter(prot %in% prot_plot)

  	adt.l$prot = factor(adt.l$prot)
  	adt.l$prot =  reorder.factor(adt.l$prot, new.order = prot_plot)
  	adt.l$cell_type = factor(adt.l$cell_type)
  	adt.l$cell_type =  reorder.factor(adt.l$cell_type, new.order = celltype_order)


	p =  ggplot(adt.l, aes(x=normalized_count, y=prot, fill = ..x..)) +
   			geom_density_ridges_gradient(alpha = 0.7, scale = 3, rel_min_height = 0.02) +
     		scale_fill_viridis(name = "normalized count", option = "B") +
    		scale_x_continuous(limits = c(-1,25)) +
    		theme_ridges(font_size = 8, line_size = 0.01) +
  		facet_grid(~cell_type)

	
	ggsave(paste0("ClusterVis/", names(ObjList_filt[i]), ".histogram.bycluster.pdf"), width = 8.5, height = 11,
    plot = p)
 
 

rm(i)
for(i in 2:length(ObjList_filt)){
	# get metadata
	md = ObjList_filt[[i]]@meta.data %>% select(adt_snn_res.1)

	# get data
	adt = GetAssayData(ObjList_filt[[i]][["CITE"]]) %>%
  	t %>%
  	as.data.frame %>%
  	rownames_to_column("cell")
	md <- md %>% rownames_to_column("cell")
	
	adt <- adt %>%
  		select(-c(cell)) %>%
  		mutate(cell_type = md$adt_snn_res.1) %>%
  		select(cell_type, everything())

	source("AverageExpression_MeanOnly.r")
	environment(AverageExpression_MeanOnly) <- asNamespace('Seurat') # need this to allow the function to call hidden seurat package functions
	
	Idents(ObjList_filt[[i]]) <- "adt_snn_res.1"
	aver = AverageExpression_MeanOnly(ObjList_filt[[i]], return.seurat=T)

	quantile_breaks <- function(xs, n = 100) {
  		breaks <- quantile(xs, probs = seq(0.5, 1, length.out = n))
  		breaks[!duplicated(breaks)]
		}

	mat_breaks <- quantile_breaks(as.matrix(GetAssayData(aver[["CITE"]])), n = 101)
	library("viridis")
	ggsave(filename = paste0("ClusterVis/", names(ObjList_filt[i]), ".bycluster_infernoColor_quantileBreaks.pdf"), width = 10, height = 11,
        plot = pheatmap(GetAssayData(aver[["CITE"]]), scale = "none", border_color=NA, inferno(length(mat_breaks) - 1), breaks = mat_breaks))

	# get hclust order of proteins and cell types
    x = pheatmap::pheatmap(GetAssayData(aver[["CITE"]]),
             cluster_rows = T, cluster_cols = T,
             fontsize_col = 10, fontsize_row = 8, border_color = NA)

    prot_order = rownames(GetAssayData(aver[["CITE"]])[x$tree_row$order, ])
    celltype_order = colnames(GetAssayData(aver[["CITE"]])[,x$tree_col$order ])

   	prot_plot = prot_order

	adt$cell_type <- as.character(adt$cell_type)
	adt.l = adt %>%
  		gather(key = prot, value = normalized_count, IgD:CD185) %>%
  		filter(prot %in% prot_plot)

  	adt.l$prot = factor(adt.l$prot)
  	adt.l$prot =  reorder.factor(adt.l$prot, new.order = prot_plot)
  	adt.l$cell_type = factor(adt.l$cell_type)
  	adt.l$cell_type =  reorder.factor(adt.l$cell_type, new.order = celltype_order)


	p =  ggplot(adt.l, aes(x=normalized_count, y=prot, fill = ..x..)) +
   			geom_density_ridges_gradient(alpha = 0.7, scale = 3, rel_min_height = 0.02) +
     		scale_fill_viridis(name = "normalized count", option = "B") +
    		scale_x_continuous(limits = c(-1,25)) +
    		theme_ridges(font_size = 8, line_size = 0.01) +
  		facet_grid(~cell_type)

	
	ggsave(paste0("ClusterVis/", names(ObjList_filt[i]), ".histogram.bycluster.pdf"), width = 8.5, height = 11,
    plot = p)
 
 
 }


US.CD4 = subset(ObjList_filt[[1]], adt_snn_res.1 %in% c("2","8","17"))
US.B = subset(ObjList_filt[[1]], adt_snn_res.1 %in% c("9","6","1","7","10","16","0"))
ObjList_select = list(US.CD4 = US.CD4, US.B = US.B)

for(i in 1:length(ObjList_select)){
	# get metadata
	md = ObjList_select[[i]]@meta.data %>% select(Tissue)

	# get data
	adt = GetAssayData(ObjList_select[[i]][["CITE"]]) %>%
  	t %>%
  	as.data.frame %>%
  	rownames_to_column("cell")
	md <- md %>% rownames_to_column("cell")
	
	adt <- adt %>%
  		select(-c(cell)) %>%
  		mutate(cell_type = md$Tissue) %>%
  		select(cell_type, everything())

	source("AverageExpression_MeanOnly.r")
	environment(AverageExpression_MeanOnly) <- asNamespace('Seurat') # need this to allow the function to call hidden seurat package functions
	
	Idents(ObjList_select[[i]]) <- "Tissue"
	aver = AverageExpression_MeanOnly(ObjList_select[[i]], return.seurat=T)

	quantile_breaks <- function(xs, n = 100) {
  		breaks <- quantile(xs, probs = seq(0.5, 1, length.out = n))
  		breaks[!duplicated(breaks)]
		}

	mat_breaks <- quantile_breaks(as.matrix(GetAssayData(aver[["CITE"]])), n = 101)
	library("viridis")
	ggsave(filename = paste0("ClusterVis/", names(ObjList_select[i]), ".byTissue_infernoColor_quantileBreaks.pdf"), width = 5, height = 11,
        plot = pheatmap(GetAssayData(aver[["CITE"]]), scale = "none", border_color=NA, inferno(length(mat_breaks) - 1), breaks = mat_breaks))
	
	# get hclust order of proteins and cell types
    x = pheatmap::pheatmap(GetAssayData(aver[["CITE"]]),
             cluster_rows = T, cluster_cols = T,
             fontsize_col = 10, fontsize_row = 8, border_color = NA)
	dev.off()
    prot_order = rownames(GetAssayData(aver[["CITE"]])[x$tree_row$order, ])
    celltype_order = colnames(GetAssayData(aver[["CITE"]])[,x$tree_col$order ])

   	prot_plot = prot_order

	adt$cell_type <- as.character(adt$cell_type)
	adt.l = adt %>%
  		gather(key = prot, value = normalized_count, CD86:CD197) %>%
  		filter(prot %in% prot_plot)

  	adt.l$prot = factor(adt.l$prot)
  	adt.l$prot =  reorder.factor(adt.l$prot, new.order = prot_plot)
  	adt.l$cell_type = factor(adt.l$cell_type)
  	adt.l$cell_type =  reorder.factor(adt.l$cell_type, new.order = celltype_order)


	p =  ggplot(adt.l, aes(x=normalized_count, y=prot, fill = ..x..)) +
   			geom_density_ridges_gradient(alpha = 0.7, scale = 3, rel_min_height = 0.02) +
     		scale_fill_viridis(name = "normalized count", option = "B") +
    		scale_x_continuous(limits = c(-1,25)) +
    		theme_ridges(font_size = 8, line_size = 0.01) +
  		facet_grid(~cell_type)

	
	ggsave(paste0("ClusterVis/", names(ObjList_select[i]), ".histogram.by.tissue.pdf"), width = 5, height = 16,
    plot = p)
 
 
 }

FeatureScatter(US.B, "cite_CD183", "cite_CD19")
FeatureScatter(US.CD4, "cite_CD45RA", "cite_CD197")
FeatureScatter(US.CD4, "cite_CD183", "cite_CD197")
FeatureScatter(US.B, "cite_CD183", "cite_CD197")
FeatureScatter(US.CD4, "cite_CD183", "cite_CD197")
FeatureScatter(US.B, "cite_CD183", "cite_CD27")
FeatureScatter(US.B, "cite_CD183", "cite_IgD")
