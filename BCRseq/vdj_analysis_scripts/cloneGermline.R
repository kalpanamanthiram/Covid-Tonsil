# Kenneth B. Hoehn
# 2/1/2021
# Match cells to barcodes and annotations
# Cluster sequences into clonal clusters
# Reconstruct germline sequences

library(alakazam)
library(shazam)
library(dplyr)
library(tidyr)
library(scoper)
library(dowser)

data = readChangeoDb(file="processed/alldata.tsv")

data = filter(data, !hash.ID %in% c("Doublet", "Negative"))

# check clonal thresholds
dist_cross = distToNearest(filter(data,locus=="IGH"),
        sequenceColumn="junction", 
        vCallColumn="v_call", jCallColumn="j_call",
        model="ham", normalize="len", nproc=1,
        cross="Donor")

pdf("intermediates/crossDistance.pdf",height=20,width=8)
ggplot(subset(dist_cross, !is.na(cross_dist_nearest)), 
             aes(x=cross_dist_nearest)) + 
    theme_bw() + 
    xlab("Cross-sample_id Hamming distance") + 
    ylab("Count") +
    geom_histogram(color="white", binwidth=0.02) +
    geom_vline(xintercept=0.12, color="firebrick", linetype=2) +
    facet_grid(Donor ~ ., scales="free_y")
dev.off()

plots = list()
pdf("intermediates/dist_to_nearest.pdf",width=6,height=6)
for(Donor in unique(dist_cross$Donor)){
   print(Donor)
   temp = filter(dist_cross,!!Donor==Donor)
   dist_ham <- distToNearest(filter(temp,locus=="IGH"), sequenceColumn="junction", 
                       vCallColumn="v_call", jCallColumn="j_call",
                       model="ham", normalize="len", nproc=2)
   output <- findThreshold(dist_ham$dist_nearest, method="density")
   threshold <- output@threshold
   g = ggplot(subset(dist_ham, !is.na(dist_nearest)),aes(x=dist_nearest,
   ,y = ..density..)) + 
     theme_bw() + 
     xlab("Hamming distance") + 
     ylab("Count") +
     scale_x_continuous(breaks=seq(0, 1, 0.1)) +
     geom_histogram(color="white", binwidth=0.02) +
     ggtitle(paste(Donor,threshold))+
     geom_histogram(
       aes(x=cross_dist_nearest,y = -..density..),
       color="white", binwidth=0.02,fill="black")+
     xlim(0,max(filter(dist_cross,
       !is.na(cross_dist_nearest))$cross_dist_nearest))+
     geom_vline(xintercept=0.1,color="grey")
   if(!is.na(threshold)){
       g = g + geom_vline(xintercept=threshold, color="firebrick", linetype=2)
   }
   plots[[Donor]] = g
   print(g)
}
dev.off()

# Scoper dies if multiple heavy chains
multi_heavy = table(filter(data, locus=="IGH")$cell_id)
multi_heavy_cells = names(multi_heavy)[multi_heavy > 1]
data = filter(data, !cell_id %in% multi_heavy_cells)

# Identify clones
clones = data %>%
    group_by(Donor) %>%
    do(as.data.frame(
    hierarchicalClones(., 0.1, cell_id = "cell_id", locus = "locus",
    only_heavy = FALSE, split_light = TRUE, cdr3 = TRUE, nproc = 2,
    verbose = FALSE, log = NULL, summarize_clones = TRUE)))

clones = ungroup(clones)

clones$clone_id = paste0(clones$Donor,"-",clones$clone_id)

comb = clones

# Create Germlines
references = readIMGT(dir = "~/share/germlines/imgt/human/vdj")
h = createGermlines(filter(comb,locus=="IGH"),references)
k = createGermlines(filter(comb,locus=="IGK"),references,locus="IGK")
l = createGermlines(filter(comb,locus=="IGL"),references,locus="IGL")

comb_germline = bind_rows(h, l, k)

# calculate SHM frequency in the V gene
comb_germline <- observedMutations(comb_germline, sequenceColumn="sequence_alignment",
       germlineColumn="germline_alignment_d_mask",
       regionDefinition=IMGT_V,
       frequency=TRUE,
       combine=TRUE, 
       nproc=3)

writeChangeoDb(comb_germline,"processed/all_cloned_data.tsv")

