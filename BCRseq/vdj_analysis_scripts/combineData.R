# Kenneth B. Hoehn
# 2/1/2021
# Combine BCR data from different sources into a single file
# annotated with metadata

library(alakazam)
library(dplyr)
library(stringr)
library(stringdist)
library(Seurat)

seurat = readRDS("raw/SortedBcells.preprocessed.ADTclustered.wMeta.Rds")

metadata = seurat@meta.data

lanes = paste0("0",6:8)

data = tibble()
for(lane in lanes){
	dir = paste0("Pam1_CITE_multi5P",lane)
	print(dir)
	h = readChangeoDb(paste0(
		dir,
		"/filtered_contig_heavy_productive-T.tsv"))
	l = readChangeoDb(paste0(
		dir,
		"/filtered_contig_light_productive-T.tsv"))

	comb = bind_rows(h,l)

	meta = filter(metadata,Lane == lane)

	# bind metadata columns to data
	print(paste(dir,"matched metadata table"))
	comb$orig.barcode = comb$cell_id
	comb$cell_id = paste0(
		unlist(lapply(strsplit(comb$cell_id,split="-"),function(x)x[1])),
		"-",lane)

	print(mean(comb$cell_id %in% rownames(meta)))
	print(mean(rownames(meta) %in% comb$cell_id))

	comb = filter(comb, cell_id %in% rownames(meta))
	meta$mcell_id = rownames(meta)

	m = match(comb$cell_id, meta$mcell_id)
	mean(comb$cell_id == meta[m,]$mcell_id)

	comb = bind_cols(comb, meta[m,])

	data = bind_rows(data,comb)
}
data = filter(data,productive == TRUE)

data$sequence_id = paste0(data$sequence_id, "-", 1:nrow(data))
writeChangeoDb(data,file="processed/alldata.tsv")
