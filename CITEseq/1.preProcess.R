#this is tested using Cicero conda environment (with R 4 and Seurat 4) on a high-performance computing node with 16 cores and at least 160 gb or ram. 
library("Seurat") #load Seurat 4
library("matrixStats")
library('tidyverse')

#set up how the lanes of the 10x ("01", "02", etc) will be assigned to different Seurat objects depending on the different unsorted or sorted cell groups in the experiment. 

LanesUnsorted = c("01","02","03","04","05")
LanesSPosBcells = c("06")
LanesSNegBcells = c("07","08")
LanesCD4Tcells = c("09","10")
LanesCD8Tcells = c("11","12")

US_data = list()
SPosBcells_data = list()
SNegBcells_data = list()
CD4Tcells_data = list()
CD8Tcells_data = list()

US_SeuratObj = list()
SPosBcells_SeuratObj = list()
SNegBcells_SeuratObj = list()
CD4Tcells_SeuratObj = list()
CD8Tcells_SeuratObj = list()

rawdir = "/hpcdata/sg/sg_data/illumina_NCI_runs/PamKalpanaTonsilCITEseq/"

### import data and Seurat object creation
## Unsorted cells
# import the cellranger counts of mRNA and features
for(i in 1:length(LanesUnsorted)){
  US_data[[i]] = Read10X_h5(paste(rawdir,"Pam1_CITE_multi5P",LanesUnsorted[[i]],"/outs/multi/count/raw_feature_bc_matrix.h5", sep="")
  )
}

# since the some of the antibody features have the same names as the mRNA, cellranger added a ".1" to them; removing it. Careful with this if you have antibody names that alread have ".1" in them
for(i in 1:length(LanesUnsorted)){
  rownames(US_data[[i]]$`Antibody Capture`) = gsub("\\.1$","", rownames(US_data[[i]]$`Antibody Capture`)) # this regex matches all ".1" at end of string, replaces with nothing
}

# create Seurat objects for each data lane
for(i in 1:length(LanesUnsorted)){
  US_SeuratObj[[i]] <- CreateSeuratObject(counts = US_data[[i]]$'Gene Expression', assay = "RNA", min.feature = 20)
  US_SeuratObj[[i]][["CITE"]] <- CreateAssayObject(counts = US_data[[i]]$`Antibody Capture`[11:148,colnames(US_SeuratObj[[i]])])
  US_SeuratObj[[i]][["HTO"]] <- CreateAssayObject(counts = US_data[[i]]$`Antibody Capture`[c(1,3:10),colnames(US_SeuratObj[[i]])])
  US_SeuratObj[[i]] <- RenameCells(US_SeuratObj[[i]], new.names = paste(substr(colnames(US_SeuratObj[[i]]), start = 1, stop = 17),LanesUnsorted[[i]], sep = ""))
  US_SeuratObj[[i]]$Batch  <- rep("B1UnSort", length(colnames(US_SeuratObj[[i]])))
  US_SeuratObj[[i]]$Type  <- rep("Unsorted", length(colnames(US_SeuratObj[[i]])))
}

US_merge <- merge(US_SeuratObj[[1]], US_SeuratObj[2:length(US_SeuratObj)])

## SPosBcells cells
# import the cellranger counts of mRNA and features
for(i in 1:length(LanesSPosBcells)){
  SPosBcells_data[[i]] = Read10X_h5(paste(rawdir,"Pam1_CITE_multi5P",LanesSPosBcells[[i]],"/outs/multi/count/raw_feature_bc_matrix.h5", sep="")
  )
}

# since the some of the antibody features have the same names as the mRNA, cellranger added a ".1" to them; removing it. Careful with this if you have antibody names that alread have ".1" in them
for(i in 1:length(LanesSPosBcells)){
  rownames(SPosBcells_data[[i]]$`Antibody Capture`) = gsub("\\.1$","", rownames(SPosBcells_data[[i]]$`Antibody Capture`)) # this regex matches all ".1" at end of string, replaces with nothing
}

# create Seurat objects for each data lane !!! only one lane in this one so not using the for loop

  SPosBcells_SeuratObj[[1]] <- CreateSeuratObject(counts = SPosBcells_data[[1]]$'Gene Expression', assay = "RNA", min.feature = 20)
  SPosBcells_SeuratObj[[1]][["CITE"]] <- CreateAssayObject(counts = SPosBcells_data[[1]]$`Antibody Capture`[11:32,colnames(SPosBcells_SeuratObj[[1]])])
  SPosBcells_SeuratObj[[1]][["HTO"]] <- CreateAssayObject(counts = SPosBcells_data[[1]]$`Antibody Capture`[c(1,3:10),colnames(SPosBcells_SeuratObj[[1]])])
  SPosBcells_SeuratObj[[1]] <- RenameCells(SPosBcells_SeuratObj[[1]], new.names = paste(substr(colnames(SPosBcells_SeuratObj[[1]]), start = 1, stop = 17), LanesSPosBcells[[1]], sep = ""))
  SPosBcells_SeuratObj[[1]]$Batch  <- rep("B1SPosBcells", length(colnames(SPosBcells_SeuratObj[[1]])))
  SPosBcells_SeuratObj[[1]]$Type  <- rep("SpositiveB", length(colnames(SPosBcells_SeuratObj[[1]])))


SPosBcells_merge <- SPosBcells_SeuratObj[[1]]

## LanesSNegBcells cells
# import the cellranger counts of mRNA and features
for(i in 1:length(LanesSNegBcells)){
  SNegBcells_data[[i]] = Read10X_h5(paste(rawdir,"Pam1_CITE_multi5P",LanesSNegBcells[[i]],"/outs/multi/count/raw_feature_bc_matrix.h5", sep="")
  )
}

# since the some of the antibody features have the same names as the mRNA, cellranger added a ".1" to them; removing it. Careful with this if you have antibody names that alread have ".1" in them
for(i in 1:length(LanesSNegBcells)){
  rownames(SNegBcells_data[[i]]$`Antibody Capture`) = gsub("\\.1$","", rownames(SNegBcells_data[[i]]$`Antibody Capture`)) # this regex matches all ".1" at end of string, replaces with nothing
}

# create Seurat objects for each data lane
# note that these cells have less antibodies in panel
for(i in 1:length(LanesSNegBcells)){
  SNegBcells_SeuratObj[[i]] <- CreateSeuratObject(counts = SNegBcells_data[[i]]$'Gene Expression', assay = "RNA", min.feature = 20)
  SNegBcells_SeuratObj[[i]][["CITE"]] <- CreateAssayObject(counts = SNegBcells_data[[i]]$`Antibody Capture`[11:32,colnames(SNegBcells_SeuratObj[[i]])])
  SNegBcells_SeuratObj[[i]][["HTO"]] <- CreateAssayObject(counts = SNegBcells_data[[i]]$`Antibody Capture`[c(1,3:10),colnames(SNegBcells_SeuratObj[[i]])])
  SNegBcells_SeuratObj[[i]] <- RenameCells(SNegBcells_SeuratObj[[i]], new.names = paste(substr(colnames(SNegBcells_SeuratObj[[i]]), start = 1, stop = 17),LanesSNegBcells[[i]], sep = ""))
  SNegBcells_SeuratObj[[i]]$Batch  <- rep("B1SNegBcells", length(colnames(SNegBcells_SeuratObj[[i]])))
  SNegBcells_SeuratObj[[i]]$Type  <- rep("SnegativeB", length(colnames(SNegBcells_SeuratObj[[i]])))

}

SNegBcells_merge <- merge(SNegBcells_SeuratObj[[1]], SNegBcells_SeuratObj[2])

# LanesCD4Tcells cells
# import the cellranger counts of mRNA and features
for(i in 1:length(LanesCD4Tcells)){
  CD4Tcells_data[[i]] = Read10X_h5(paste(rawdir,"Pam1_CITE_multi5P",LanesCD4Tcells[[i]],"/outs/multi/count/raw_feature_bc_matrix.h5", sep="")
  )
}

# since the some of the antibody features have the same names as the mRNA, cellranger added a ".1" to them; removing it. Careful with this if you have antibody names that alread have ".1" in them
for(i in 1:length(LanesCD4Tcells)){
  rownames(CD4Tcells_data[[i]]$`Antibody Capture`) = gsub("\\.1$","", rownames(CD4Tcells_data[[i]]$`Antibody Capture`)) # this regex matches all ".1" at end of string, replaces with nothing
}

# create Seurat objects for each data lane
# note that these cells have less antibodies in panel
for(i in 1:length(LanesCD4Tcells)){
  CD4Tcells_SeuratObj[[i]] <- CreateSeuratObject(counts = CD4Tcells_data[[i]]$'Gene Expression', assay = "RNA", min.feature = 20)
  CD4Tcells_SeuratObj[[i]][["CITE"]] <- CreateAssayObject(counts = CD4Tcells_data[[i]]$`Antibody Capture`[11:32,colnames(CD4Tcells_SeuratObj[[i]])])
  CD4Tcells_SeuratObj[[i]][["HTO"]] <- CreateAssayObject(counts = CD4Tcells_data[[i]]$`Antibody Capture`[c(1,3:10),colnames(CD4Tcells_SeuratObj[[i]])])
  CD4Tcells_SeuratObj[[i]] <- RenameCells(CD4Tcells_SeuratObj[[i]], new.names = paste(substr(colnames(CD4Tcells_SeuratObj[[i]]), start = 1, stop = 17),LanesCD4Tcells[[i]], sep = ""))
  CD4Tcells_SeuratObj[[i]]$Batch  <- rep("B1CD4Tcells", length(colnames(CD4Tcells_SeuratObj[[i]])))
  CD4Tcells_SeuratObj[[i]]$Type  <- rep("CD4Tcells", length(colnames(CD4Tcells_SeuratObj[[i]])))

}

CD4Tcells_merge <- merge(CD4Tcells_SeuratObj[[1]], CD4Tcells_SeuratObj[2])

# LanesCD8Tcells cells
# import the cellranger counts of mRNA and features
for(i in 1:length(LanesCD8Tcells)){
  CD8Tcells_data[[i]] = Read10X_h5(paste(rawdir,"Pam1_CITE_multi5P",LanesCD8Tcells[[i]],"/outs/multi/count/raw_feature_bc_matrix.h5", sep="")
  )
}

# since the some of the antibody features have the same names as the mRNA, cellranger added a ".1" to them; removing it. Careful with this if you have antibody names that alread have ".1" in them
for(i in 1:length(LanesCD8Tcells)){
  rownames(CD8Tcells_data[[i]]$`Antibody Capture`) = gsub("\\.1$","", rownames(CD8Tcells_data[[i]]$`Antibody Capture`)) # this regex matches all ".1" at end of string, replaces with nothing
}

# create Seurat objects for each data lane
# note that these cells have less antibodies in panel
for(i in 1:length(LanesCD8Tcells)){
  CD8Tcells_SeuratObj[[i]] <- CreateSeuratObject(counts = CD8Tcells_data[[i]]$'Gene Expression', assay = "RNA", min.feature = 20)
  CD8Tcells_SeuratObj[[i]][["CITE"]] <- CreateAssayObject(counts = CD8Tcells_data[[i]]$`Antibody Capture`[11:32,colnames(CD8Tcells_SeuratObj[[i]])])
  CD8Tcells_SeuratObj[[i]][["HTO"]] <- CreateAssayObject(counts = CD8Tcells_data[[i]]$`Antibody Capture`[c(1,3:10),colnames(CD8Tcells_SeuratObj[[i]])])
  CD8Tcells_SeuratObj[[i]] <- RenameCells(CD8Tcells_SeuratObj[[i]], new.names = paste(substr(colnames(CD8Tcells_SeuratObj[[i]]), start = 1, stop = 17),LanesCD8Tcells[[i]], sep = ""))
  CD8Tcells_SeuratObj[[i]]$Batch  <- rep("B1CD8Tcells", length(colnames(CD8Tcells_SeuratObj[[i]])))
  CD8Tcells_SeuratObj[[i]]$Type  <- rep("CD8Tcells", length(colnames(CD8Tcells_SeuratObj[[i]])))

}

CD8Tcells_merge <- merge(CD8Tcells_SeuratObj[[1]], CD8Tcells_SeuratObj[2])

### Filter cells for antibody clumps before doing HTOdemux. Checking data distributions to set cutoff.
quantile(US_merge$nCount_CITE, probs = c(0, .1, .25, .5, .75, .95, .99, .999, .9999, 1))
quantile(SPosBcells_merge$nCount_CITE, probs = c(0, .1, .25, .5, .75, .95, .99, .999, .9999, 1))
quantile(SNegBcells_merge$nCount_CITE, probs = c(0, .1, .25, .5, .75, .95, .99, .999, .9999, 1))
quantile(CD4Tcells_merge$nCount_CITE, probs = c(0, .1, .25, .5, .75, .95, .99, .999, .9999, 1))
quantile(CD8Tcells_merge$nCount_CITE, probs = c(0, .1, .25, .5, .75, .95, .99, .999, .9999, 1))

US_merge <- subset(US_merge, subset = nCount_CITE < quantile(US_merge$nCount_CITE, probs = c(0.995)))
SPosBcells_merge <- subset(SPosBcells_merge, subset = nCount_CITE < quantile(SPosBcells_merge$nCount_CITE, probs = c(0.995)))
SNegBcells_merge <- subset(SNegBcells_merge, subset = nCount_CITE < quantile(SNegBcells_merge$nCount_CITE, probs = c(0.995)))
CD4Tcells_merge <- subset(CD4Tcells_merge, subset = nCount_CITE <  quantile(CD4Tcells_merge$nCount_CITE, probs = c(0.995)))
CD8Tcells_merge <- subset(CD8Tcells_merge, subset = nCount_CITE <  quantile(CD8Tcells_merge$nCount_CITE, probs = c(0.995)))

quantile(US_merge$nCount_HTO, probs = c(0, .1, .25, .5, .75, .95, .99, .999, .9999, 1))
quantile(SPosBcells_merge$nCount_HTO, probs = c(0, .1, .25, .5, .75, .95, .99, .999, .9999, 1))
quantile(SNegBcells_merge$nCount_HTO, probs = c(0, .1, .25, .5, .75, .95, .99, .999, .9999, 1))
quantile(CD4Tcells_merge$nCount_HTO, probs = c(0, .1, .25, .5, .75, .95, .99, .999, .9999, 1))
quantile(CD8Tcells_merge$nCount_HTO, probs = c(0, .1, .25, .5, .75, .95, .99, .999, .9999, 1))

US_merge <- subset(US_merge, subset = nCount_HTO <  quantile(US_merge$nCount_HTO, probs = c(0.995)))
SPosBcells_merge <- subset(SPosBcells_merge, subset = nCount_HTO <  quantile(SPosBcells_merge$nCount_HTO, probs = c(0.995)))
SNegBcells_merge <- subset(SNegBcells_merge, subset = nCount_HTO <  quantile(SNegBcells_merge$nCount_HTO, probs = c(0.995))) 
CD4Tcells_merge <- subset(CD4Tcells_merge, subset = nCount_HTO <  quantile(CD4Tcells_merge$nCount_HTO, probs = c(0.995)))
CD8Tcells_merge <- subset(CD8Tcells_merge, subset = nCount_HTO <  quantile(CD8Tcells_merge$nCount_HTO, probs = c(0.995)))

# Do HTODemux to get the negatives list for DSB normalization

ObjList <- list(US = US_merge,
				SPosBcells = SPosBcells_merge,
				SNegBcells = SNegBcells_merge,
				CD4Tcells = CD4Tcells_merge,
				CD8Tcells = CD8Tcells_merge)
				 
ObjList["MergeBcells"] <- merge(ObjList[[2]], ObjList[3])
 
for(i in 1:length(ObjList)){
	ObjList[[i]] <- NormalizeData(ObjList[[i]], assay= "HTO", normalization.method = "CLR", margin=2)
	ObjList[[i]] <- ScaleData(ObjList[[i]], assay = "HTO", model.use = "negbinom")
	ObjList[[i]] = HTODemux(ObjList[[i]], positive.quantile = 0.9999)
	Idents(ObjList[[i]]) <- "hash.ID"
}

#add lane number as metadata (this is already in the cell names, taking it from there). "\\D" matches all non-digits, then gsub replaces that with nothing
for(i in 1:length(ObjList)){
	ObjList[[i]]$Lane <- gsub("\\D","",colnames(ObjList[[i]]))
}	

#visually inspect HTO calls
dir.create("HTOvisualization")
for(i in 1:length(ObjList)){
	ggsave(filename = paste0("HTOvisualization/", names(ObjList[i]),".CLRdemux.HTO4vs10.pdf"),
		   plot = FeatureScatter(subset(ObjList[[i]], idents = c("HTO4","HTO10", "HTO1")), feature1="hto_HTO4", feature2="hto_HTO10"),
		   width = 4, height = 4)
	ggsave(filename = paste0("HTOvisualization/", names(ObjList[i]),".CLRdemux.HTO1vs10.pdf"),
		   plot = FeatureScatter(subset(ObjList[[i]], idents = c("HTO4","HTO10", "HTO1")), feature1="hto_HTO1", feature2="hto_HTO10"),
		   width = 4, height = 4)
	ggsave(filename = paste0("HTOvisualization/", names(ObjList[i]),".CLRdemux.HTO3vs10.pdf"),
		   plot = FeatureScatter(subset(ObjList[[i]], idents = c("HTO4","HTO10", "HTO3")), feature1="hto_HTO3", feature2="hto_HTO10"),
		   width = 4, height = 4)
	ggsave(filename = paste0("HTOvisualization/", names(ObjList[i]),".CLRdemux.HTO5vs10.pdf"),
           plot = FeatureScatter(subset(ObjList[[i]], idents = c("HTO4","HTO10", "HTO5")), feature1="hto_HTO5", feature2="hto_HTO10"),
           width = 4, height = 4)
    ggsave(filename = paste0("HTOvisualization/", names(ObjList[i]),".CLRdemux.HTO6vs10.pdf"),
		   plot = FeatureScatter(subset(ObjList[[i]], idents = c("HTO4","HTO10", "HTO6")), feature1="hto_HTO6", feature2="hto_HTO10"),
		   width = 4, height = 4)
	ggsave(filename = paste0("HTOvisualization/", names(ObjList[i]),".CLRdemux.HTO7vs10.pdf"),
		   plot = FeatureScatter(subset(ObjList[[i]], idents = c("HTO4","HTO10", "HTO7")), feature1="hto_HTO7", feature2="hto_HTO10"),
		   width = 4, height = 4)
    ggsave(filename = paste0("HTOvisualization/", names(ObjList[i]),".CLRdemux.HTO8vs10.pdf"),
		   plot = FeatureScatter(subset(ObjList[[i]], idents = c("HTO4","HTO10", "HTO8")), feature1="hto_HTO8", feature2="hto_HTO10"),
		   width = 4, height = 4)
    ggsave(filename = paste0("HTOvisualization/", names(ObjList[i]),".CLRdemux.HTO9vs10.pdf"),
           plot = FeatureScatter(subset(ObjList[[i]], idents = c("HTO4","HTO10", "HTO9")), feature1="hto_HTO9", feature2="hto_HTO10"),
		   width = 4, height = 4)
    ggsave(filename = paste0("HTOvisualization/", names(ObjList[i]),".CLRdemux.HTOheatmap.pdf"),
 			plot = HTOHeatmap(ObjList[[i]]),
 			width = 4, height = 4)
}

# DSB Normalization
library("dsb")


#move negatives to seperate object
quantile(ObjList[[1]]$nFeature_RNA, probs = c(0, .2, .4, .6, .8, .9, 1)) ## most of the droplets per lane probably empty (expected 12k per lane, forced cellranger to output 30k)

ObjList_neg = list()
for(i in 1:length(ObjList)){
	ObjList_neg[[i]] <- subset(ObjList[[i]], subset = hash.ID == "Negative" & nFeature_RNA < 100)
	ObjList_neg[[i]] <- subset(ObjList_neg[[i]], cells = sample(Cells(ObjList_neg[[i]]),1000))
	ObjList[[i]] <- subset(ObjList[[i]], nFeature_RNA > 100)

	
}

DSB.HTO.list = list()
for(i in 1:length(ObjList)){
	DSB.HTO.list[[i]] = DSBNormalizeProtein(cell_protein_matrix = as.matrix(GetAssayData(ObjList[[i]][["HTO"]], slot = "counts")),
                                     empty_drop_matrix = as.matrix(GetAssayData(ObjList_neg[[i]][["HTO"]], slot = "counts")), denoise.counts = FALSE)
	ObjList[[i]] <- SetAssayData(ObjList[[i]], assay = "HTO", new.data = DSB.HTO.list[[i]], slot = "data")
}


#visualize the DSB norm HTOdemux results
for(i in 1:length(ObjList)){
	ggsave(filename = paste0("HTOvisualization/", names(ObjList[i]),".DSBdemux.HTO4vs10.pdf"),
		   plot = FeatureScatter(subset(ObjList[[i]], idents = c("HTO4","HTO10", "HTO1")), feature1="hto_HTO4", feature2="hto_HTO10"),
		   width = 4, height = 4)
	ggsave(filename = paste0("HTOvisualization/", names(ObjList[i]),".DSBdemux.HTO1vs10.pdf"),
		   plot = FeatureScatter(subset(ObjList[[i]], idents = c("HTO4","HTO10", "HTO1")), feature1="hto_HTO1", feature2="hto_HTO10"),
		   width = 4, height = 4)
	ggsave(filename = paste0("HTOvisualization/", names(ObjList[i]),".DSBdemux.HTO3vs10.pdf"),
		   plot = FeatureScatter(subset(ObjList[[i]], idents = c("HTO4","HTO10", "HTO3")), feature1="hto_HTO3", feature2="hto_HTO10"),
		   width = 4, height = 4)
	ggsave(filename = paste0("HTOvisualization/", names(ObjList[i]),".DSBdemux.HTO5vs10.pdf"),
           plot = FeatureScatter(subset(ObjList[[i]], idents = c("HTO4","HTO10", "HTO5")), feature1="hto_HTO5", feature2="hto_HTO10"),
           width = 4, height = 4)
    ggsave(filename = paste0("HTOvisualization/", names(ObjList[i]),".DSBdemux.HTO6vs10.pdf"),
		   plot = FeatureScatter(subset(ObjList[[i]], idents = c("HTO4","HTO10", "HTO6")), feature1="hto_HTO6", feature2="hto_HTO10"),
		   width = 4, height = 4)
	ggsave(filename = paste0("HTOvisualization/", names(ObjList[i]),".DSBdemux.HTO7vs10.pdf"),
		   plot = FeatureScatter(subset(ObjList[[i]], idents = c("HTO4","HTO10", "HTO7")), feature1="hto_HTO7", feature2="hto_HTO10"),
		   width = 4, height = 4)
    ggsave(filename = paste0("HTOvisualization/", names(ObjList[i]),".DSBdemux.HTO8vs10.pdf"),
		   plot = FeatureScatter(subset(ObjList[[i]], idents = c("HTO4","HTO10", "HTO8")), feature1="hto_HTO8", feature2="hto_HTO10"),
		   width = 4, height = 4)
    ggsave(filename = paste0("HTOvisualization/", names(ObjList[i]),".DSBdemux.HTO9vs10.pdf"),
           plot = FeatureScatter(subset(ObjList[[i]], idents = c("HTO4","HTO10", "HTO9")), feature1="hto_HTO9", feature2="hto_HTO10"),
		   width = 4, height = 4)
    ggsave(filename = paste0("HTOvisualization/", names(ObjList[i]),".DSBdemux.HTOheatmap.pdf"),
 			plot = HTOHeatmap(ObjList[[i]]),
 			width = 4, height = 4)
}

# Some demux results don't look good. Since the DSB norm defines the difference of the signal above the background of negative drops, we can try and use a cutoff based approach

HTONames = rownames(ObjList[[1]][["HTO"]])

for(i in 1:length(ObjList)){
	pdf(paste0("HTOvisualization/", names(ObjList[i]),".HTOcheck.pdf"))
	for (j in 1:length(HTONames)) {
  		hist(as.matrix(GetAssayData(ObjList[[i]], assay = "HTO", slot = "data"))[j,], breaks=100, main = HTONames[j])
	}
	dev.off()
}

# manual threshold HTO calls
hithres = list( US = c(15, 10, 20, 18, 4, 12, 5, 15, 5),
				SPosBcells = c(40, 20, 20, 15, 13, 6, 10, 50, 30),
				SNegBcells = c(6, 10, 10, 8, 9, 7, 6, 7.5, 7.5),
				CD4Tcells = c(12, 20, 10, 15, 8, 12, 5, 12.5, 12),
				CD8Tcells = c(10, 15, 10, 15, 17, 16, 10, 17, 17),
				MergeBcells = c(7, 8, 10, 8, 8.5, 5, 6, 8, 7)
			)
lothres = list( US = c(10, 7, 10, 10, 3.5, 10, 4.5, 12, 4.5),
				SPosBcells = c(25, 18, 15, 10, 12.5, 5, 8, 40, 25),
				SNegBcells = c(5, 7, 8, 6, 8, 6, 5, 6, 7),
				CD4Tcells = c(10, 15, 8, 13, 7, 10, 4, 10, 10),
				CD8Tcells = c(7, 12.5, 8, 12.5, 15, 15, 7, 15, 15),
				MergeBcells = c(5, 6, 8, 6, 8, 4, 5, 7, 6)
			)

automaticHTOsHiThres = list()
rm(i)
rm(j)
for(i in 1:length(ObjList)){
	automaticHTOsHiThres[[i]] <- character()
	for(j in 1:length(colnames(ObjList[[i]]))){
  		automaticHTOsHiThres[[i]][j] = paste(HTONames[GetAssayData(ObjList[[i]][["HTO"]], slot = "data")[,j] > hithres[[i]]], collapse="+")
	}
}

# manual lower threshold HTO calls, for doublet removal
automaticHTOsloThres = list()
rm(i)
rm(j)
for(i in 1:length(ObjList)){
	automaticHTOsloThres[[i]] <- character()
	for(j in 1:length(colnames(ObjList[[i]]))){
  		automaticHTOsloThres[[i]][j] = paste(HTONames[GetAssayData(ObjList[[i]][["HTO"]], slot = "data")[,j] > lothres[[i]]], collapse="+")
	}
}

singlets <- list()
autoHashcalls <- list()
for(i in 1:length(ObjList)){
	singlets[[i]] = automaticHTOsHiThres[[i]] %in% HTONames & automaticHTOsloThres[[i]] %in% HTONames
	autoHashcalls[[i]] = ifelse(singlets[[i]], automaticHTOsHiThres[[i]], "NonSinglet")
	ObjList[[i]] <- AddMetaData(object = ObjList[[i]], metadata = autoHashcalls[[i]], col.name = "autoHashcalls")
}

table(ObjList[[1]]$hash.ID, ObjList[[1]]$autoHashcalls)
table(ObjList[[2]]$hash.ID, ObjList[[2]]$autoHashcalls)
table(ObjList[[3]]$hash.ID, ObjList[[3]]$autoHashcalls)
table(ObjList[[4]]$hash.ID, ObjList[[4]]$autoHashcalls)
table(ObjList[[5]]$hash.ID, ObjList[[5]]$autoHashcalls)
table(ObjList[[6]]$hash.ID, ObjList[[6]]$autoHashcalls)

Idents(ImmDevSeurat) <- "autoHashcalls"

#visualize the autoHashcalls results
for(i in 1:length(ObjList)){
	Idents(ObjList[[i]]) <- "autoHashcalls"
	ggsave(filename = paste0("HTOvisualization/", names(ObjList[i]),".AutohashCalls.HTO4vs10.pdf"),
		   plot = FeatureScatter(subset(ObjList[[i]], idents = c("HTO4","HTO10", "HTO1")), feature1="hto_HTO4", feature2="hto_HTO10"),
		   width = 4, height = 4)
	ggsave(filename = paste0("HTOvisualization/", names(ObjList[i]),".AutohashCalls.HTO1vs10.pdf"),
		   plot = FeatureScatter(subset(ObjList[[i]], idents = c("HTO4","HTO10", "HTO1")), feature1="hto_HTO1", feature2="hto_HTO10"),
		   width = 4, height = 4)
	ggsave(filename = paste0("HTOvisualization/", names(ObjList[i]),".AutohashCalls.HTO3vs10.pdf"),
		   plot = FeatureScatter(subset(ObjList[[i]], idents = c("HTO4","HTO10", "HTO3")), feature1="hto_HTO3", feature2="hto_HTO10"),
		   width = 4, height = 4)
	ggsave(filename = paste0("HTOvisualization/", names(ObjList[i]),".AutohashCalls.HTO5vs10.pdf"),
           plot = FeatureScatter(subset(ObjList[[i]], idents = c("HTO4","HTO10", "HTO5")), feature1="hto_HTO5", feature2="hto_HTO10"),
           width = 4, height = 4)
    ggsave(filename = paste0("HTOvisualization/", names(ObjList[i]),".AutohashCalls.HTO6vs10.pdf"),
		   plot = FeatureScatter(subset(ObjList[[i]], idents = c("HTO4","HTO10", "HTO6")), feature1="hto_HTO6", feature2="hto_HTO10"),
		   width = 4, height = 4)
	ggsave(filename = paste0("HTOvisualization/", names(ObjList[i]),".AutohashCalls.HTO7vs10.pdf"),
		   plot = FeatureScatter(subset(ObjList[[i]], idents = c("HTO4","HTO10", "HTO7")), feature1="hto_HTO7", feature2="hto_HTO10"),
		   width = 4, height = 4)
    ggsave(filename = paste0("HTOvisualization/", names(ObjList[i]),".AutohashCalls.HTO8vs10.pdf"),
		   plot = FeatureScatter(subset(ObjList[[i]], idents = c("HTO4","HTO10", "HTO8")), feature1="hto_HTO8", feature2="hto_HTO10"),
		   width = 4, height = 4)
    ggsave(filename = paste0("HTOvisualization/", names(ObjList[i]),".AutohashCalls.HTO9vs10.pdf"),
           plot = FeatureScatter(subset(ObjList[[i]], idents = c("HTO4","HTO10", "HTO9")), feature1="hto_HTO9", feature2="hto_HTO10"),
		   width = 4, height = 4)
    ggsave(filename = paste0("HTOvisualization/", names(ObjList[i]),".AutohashCalls.HTOheatmap.pdf"),
 			plot = HTOHeatmap(ObjList[[i]]),
 			width = 4, height = 4)
}


#remove Non-Singlets from the data into new object
 
ObjList_filt = list()
for(i in 1:length(ObjList)){
	ObjList_filt[[i]] <- subset(ObjList[[i]], subset = autoHashcalls == "NonSinglet", invert=TRUE)
	
}

names(ObjList_filt) <- names(ObjList)
dir.create("CountsPerSample")

write.csv(table(ObjList[[1]]$hash.ID, ObjList[[1]]$autoHashcalls), "CountsPerSample/Unsorted.CellsPerSample.HashCalls.Overlap.csv")
write.csv(table(ObjList[[2]]$hash.ID, ObjList[[2]]$autoHashcalls), "CountsPerSample/SPosBcells.CellsPerSample.HashCalls.Overlap.csv")
write.csv(table(ObjList[[3]]$hash.ID, ObjList[[3]]$autoHashcalls), "CountsPerSample/SNegBcells.CellsPerSample.HashCalls.Overlap.csv")
write.csv(table(ObjList[[4]]$hash.ID, ObjList[[4]]$autoHashcalls), "CountsPerSample/CD4Tcells.CellsPerSample.HashCalls.Overlap.csv")
write.csv(table(ObjList[[5]]$hash.ID, ObjList[[5]]$autoHashcalls), "CountsPerSample/CD8Tcells.CellsPerSample.HashCalls.Overlap.csv")
write.csv(table(ObjList[[6]]$hash.ID, ObjList[[6]]$autoHashcalls), "CountsPerSample/MergeBcells.CellsPerSample.HashCalls.Overlap.csv")

write.csv(table(ObjList_filt[[1]]$autoHashcalls), "CountsPerSample/Unsorted.CellsPerSample.gatingHashCalls.csv")
write.csv(table(ObjList_filt[[2]]$autoHashcalls), "CountsPerSample/SPosBcells.CellsPerSample.gatingHashCalls.csv")
write.csv(table(ObjList_filt[[3]]$autoHashcalls), "CountsPerSample/SNegBcells.CellsPerSample.gatingHashCalls.csv")
write.csv(table(ObjList_filt[[4]]$autoHashcalls), "CountsPerSample/CD4Tcells.CellsPerSample.gatingHashCalls.csv")
write.csv(table(ObjList_filt[[5]]$autoHashcalls), "CD8Tcells.CellsPerSample.gatingHashCalls.csv")
write.csv(table(ObjList_filt[[6]]$autoHashcalls), "CountsPerSample/MergeBcells.CellsPerSample.gatingHashCalls.csv")


#DSB on surface proteins
DSB.ADT.list = list()
#doing unsorted cells seperately, since it has isotype controls
isotype.control.name.vec = c("IgG2bKiso", "IgG1Kiso", "IgG2aKiso", "ratIgG2bKiso", "ratIgG1Kiso", "ratIgG2aKiso", "ArmHamsterIgGiso")
DSB.ADT.list[[1]] = DSBNormalizeProtein(cell_protein_matrix = as.matrix(GetAssayData(ObjList_filt[[1]][["CITE"]], slot = "counts")),
                                     empty_drop_matrix = as.matrix(GetAssayData(ObjList_neg[[1]][["CITE"]], slot = "counts")), denoise.counts = TRUE, isotype.control.name.vec = isotype.control.name.vec)
ObjList_filt[[1]] <- SetAssayData(ObjList_filt[[1]], assay = "CITE", new.data = DSB.ADT.list[[1]], slot = "data")

for(i in 2:length(ObjList_filt)){
	DSB.ADT.list[[i]] = DSBNormalizeProtein(cell_protein_matrix = as.matrix(GetAssayData(ObjList_filt[[i]][["CITE"]], slot = "counts")),
                                     empty_drop_matrix = as.matrix(GetAssayData(ObjList_neg[[i]][["CITE"]], slot = "counts")), denoise.counts = FALSE)
	ObjList_filt[[i]] <- SetAssayData(ObjList_filt[[i]], assay = "CITE", new.data = DSB.ADT.list[[i]], slot = "data")
}


##QC plots and further filtering
mito.genes = grep(pattern = "^MT-", x = rownames(ObjList_filt[[1]]), value = TRUE)

for(i in 1:length(ObjList_filt)){
	ggsave(filename = paste0("HTOvisualization/", names(ObjList[i]),".VlnPlot_nFeature_RNA.pdf"),
		   plot = VlnPlot(ObjList_filt[[i]], features = "nFeature_RNA", group.by = "autoHashcalls", pt.size = 0),
		   width = 8, height = 4)
	ggsave(filename = paste0("HTOvisualization/", names(ObjList[i]),".VlnPlot_nCount_RNA.pdf"),
		   plot = VlnPlot(ObjList_filt[[i]], features = "nCount_RNA", group.by = "autoHashcalls", pt.size = 0),
		   width = 8, height = 4)
	ggsave(filename = paste0("HTOvisualization/", names(ObjList[i]),".VlnPlot_nCount_HTO.pdf"),
		   plot = VlnPlot(ObjList_filt[[i]], features = "nCount_HTO", group.by = "autoHashcalls", pt.size = 0),
		   width = 8, height = 4)
	ggsave(filename = paste0("HTOvisualization/", names(ObjList[i]),".VlnPlot_nCount_CITE.pdf"),
		   plot = VlnPlot(ObjList_filt[[i]], features = "nCount_CITE", group.by = "autoHashcalls", pt.size = 0),
		   width = 8, height = 4)
	ObjList_filt[[i]] <- AddMetaData(object = ObjList_filt[[i]], metadata = Matrix::colSums(ObjList_filt[[i]][mito.genes,])/Matrix::colSums(ObjList_filt[[i]]), col.name = "percent.mito")
	ggsave(filename = paste0("HTOvisualization/", names(ObjList[i]),".VlnPlot_percentMito.pdf"),
		   plot = VlnPlot(ObjList_filt[[i]], features = "percent.mito", group.by = "autoHashcalls", pt.size = 0),
		   width = 8, height = 4)
}

for(i in 1:length(ObjList_filt)){
	ObjList_filt[[i]] = subset(ObjList_filt[[i]], subset = nCount_RNA < 25000 & percent.mito < 0.30)
}

## Clustering
dir.create("ClusterVis")
for(i in 1:length(ObjList_filt)){
  adt.dist.temp <- dist(t(GetAssayData(ObjList_filt[[i]], assay = "CITE", slot = "data")))
  ObjList_filt[[i]][["adt_snn"]] <- FindNeighbors(adt.dist.temp, nn.eps = 1)$snn
  ObjList_filt[[i]] <- FindClusters(ObjList_filt[[i]], resolution = c(1), graph.name = "adt_snn", algorithm = 1)
  ObjList_filt[[i]] <- DietSeurat(ObjList_filt[[i]], scale.data = TRUE)
}

  
for(i in 1:length(ObjList_filt)){
  ObjList_filt[[i]] <- RunUMAP(ObjList_filt[[i]], assay = "CITE", features = rownames(ObjList_filt[[i]][["CITE"]]), n_neighbors=50L, min_dist=1)
  ggsave(paste0("ClusterVis/", names(ObjList_filt[i]), ".umap.png"), width = 5, height = 5,
    plot = AugmentPlot(DimPlot(ObjList_filt[[i]], pt.size = 0.5, reduction = "umap", group.by = "adt_snn_res.1", label=FALSE)) + NoLegend()
 )
}

library('pheatmap')
for(i in 1:length(ObjList_filt)){
  Idents(ObjList_filt[[i]]) <- "adt_snn_res.1"
  ObjList_filt[[i]] <- SetAssayData(ObjList_filt[[i]], assay = "CITE", slot = "scale.data", new.data = GetAssayData(ObjList_filt[[i]], assay = "CITE", slot = "data"))
  aver.temp = AverageExpression(ObjList_filt[[i]], assays = "CITE", slot = "scale.data", return.seurat=F)
  protstoplot = rownames(aver.temp$CITE)[which(rowMaxs(as.matrix(aver.temp$CITE)) > 0.5)]
  
  ggsave(filename = paste0("ClusterVis/", names(ObjList_filt[i]),".averExp_res.1_CITE.pdf"), width = 5, height = 12,
         plot = pheatmap(aver.temp$CITE[protstoplot,], scale = "none", border_color=NA))
  
}

for(i in 1:length(ObjList_filt)){
  pheatmap(prop.table(x = table(ObjList_filt[[i]]$adt_snn_res.1, ObjList_filt[[i]]$autoHashcalls), margin =  2)  %>% '*'(100) %>% round(2),
           filename = paste0("ClusterVis/", names(ObjList_filt[i]),".ClustPercentageperSample.pdf"), width = 3, height = 6,
           cluster_cols = FALSE, cluster_rows = FALSE, display_numbers = TRUE

  )
}


dir.create("SeuratObjects")
saveRDS(ObjList_filt, "SeuratObjects/ObjList_filt.preprocessed.ADTclustered.Rds")
saveRDS(ObjList, "SeuratObjects/ObjList.preprocessed.Rds")


#DE genes per cluster
clustRNAmarkers = list()
aver.celltype.CITEclust = list()
for(i in 1:length(ObjList_filt)){
  clustRNAmarkers[[i]] = FindAllMarkers(ObjList_filt[[i]], assay = "RNA")
  
}

for(i in 1:length(ObjList_filt)){
celltypeheatgenes = clustRNAmarkers[[i]] %>% group_by(cluster) %>% top_n(10, avg_log2FC)
aver.celltype.CITEclust[[i]] = AverageExpression(ObjList_filt[[i]], return.seurat=T)
ggsave(filename = paste0("ClusterVis/", names(ObjList_filt[i]),".averExp_CITEclust_DE.RNA.pdf"), width = 8, height = 18,
       plot = pheatmap(aver.celltype.CITEclust[[i]][["RNA"]][unique(celltypeheatgenes$gene),], scale = "row", border_color=NA))
}

## add metadata with sample names etc.
for(i in 1:length(ObjList_filt)){
	Idents(ObjList_filt[[i]]) <- "autoHashcalls"
    ObjList_filt[[i]] <- RenameIdents(ObjList_filt[[i]],
    								  'HTO1' = "CNMC89_PBMC_COVIDpos",
    								  'HTO3' = "CNMC99_PBMC_COVIDneg",
    								  'HTO4' = "CNMC71_PBMC_COVIDpos",
    								  'HTO5' = "CNMC71_Adenoid_COVIDpos",
    								  'HTO6' = "CNMC71_Tonsil_COVIDpos",
    								  'HTO7' = "CNMC89_Adenoid_COVIDpos",
    								  'HTO8' = "CNMC89_Tonsil_COVIDpos",
    								  'HTO9' = "CNMC99_Adenoid_COVIDneg",
    								  'HTO10' = "CNMC99_Tonsil_COVIDneg")
    
}
for(i in 1:length(ObjList_filt)){
	ObjList_filt[[i]]$Sample_Tissue_Class <- Idents(ObjList_filt[[i]])	
	ObjList_filt[[i]]$Donor <- sapply(strsplit(as.character(ObjList_filt[[i]]$Sample_Tissue_Class), split = "_"),'[',1)	
	ObjList_filt[[i]]$Tissue <- sapply(strsplit(as.character(ObjList_filt[[i]]$Sample_Tissue_Class), split = "_"),'[',2)							  
	ObjList_filt[[i]]$Class <- sapply(strsplit(as.character(ObjList_filt[[i]]$Sample_Tissue_Class), split = "_"),'[',3)							  
						  
}


write.csv(table(ObjList_filt[[1]]$Sample_Tissue_Class), "CountsPerSample/Unsorted.CellsPerSampleTissueClass.gatingHashCalls.csv")
write.csv(table(ObjList_filt[[2]]$Sample_Tissue_Class), "CountsPerSample/SPosBcells.CellsPerSampleTissueClass.gatingHashCalls.csv")
write.csv(table(ObjList_filt[[3]]$Sample_Tissue_Class), "CountsPerSample/SNegBcells.CellsPerSampleTissueClass.gatingHashCalls.csv")
write.csv(table(ObjList_filt[[4]]$Sample_Tissue_Class), "CountsPerSample/CD4Tcells.CellsPerSampleTissueClass.gatingHashCalls.csv")
write.csv(table(ObjList_filt[[5]]$Sample_Tissue_Class), "CountsPerSample/CD8Tcells.CellsPerSampleTissueClass.gatingHashCalls.csv")
write.csv(table(ObjList_filt[[6]]$Sample_Tissue_Class), "CountsPerSample/MergeBcells.CellsPerSampleTissueClass.gatingHashCalls.csv")
write.csv(table(ObjList_filt[[6]]$Sample_Tissue_Class,  ObjList_filt[[6]]$Type), "CountsPerSample/MergeBcells.CellsPerSampleTissueClassPersortedtype.gatingHashCalls.csv")

saveRDS(ObjList_filt, "SeuratObjects/ObjList_filt.preprocessed.ADTclustered.wMeta.Rds")
saveRDS(ObjList_filt[[6]], "SeuratObjects/SortedBcells.preprocessed.ADTclustered.wMeta.Rds")

