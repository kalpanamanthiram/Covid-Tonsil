# Kenneth B. Hoehn
# 5/29/2022
# Perform BCR analysis 

library(alakazam)
library(dplyr)
library(tidyr)
library(ggpubr)
library(dowser)
library(ggpubr)
library(ggtree)
library(gridExtra)
library(scoper)
library(RColorBrewer)

# set up palettes
subisotype_levels = c("IGHM", "IGHD", "IGHG3", "IGHG1", "IGHA1", "IGHG2", "IGHG4", "IGHE", "IGHA2")
subisotype_palette = brewer.pal(length(subisotype_levels), "Paired")
names(subisotype_palette) = subisotype_levels

isotype_palette = c()
isotype_palette["IGHM"] = subisotype_palette["IGHM"]
isotype_palette["IGHD"] = subisotype_palette["IGHD"]
isotype_palette["IGHG"] = subisotype_palette["IGHG1"]
isotype_palette["IGHA"] = subisotype_palette["IGHA1"]
isotype_palette["IGHE"] = subisotype_palette["IGHE"]

tissue_palette = c(
"PBMC"="#ff0000",
"Adenoid"="#ffff00",
"Tonsil"="#a020f0"
  )

cloned = readChangeoDb("processed/all_cloned_data.tsv")
cloned = filter(cloned, locus == "IGH")
cloned = filter(cloned, !is.na(c_call))
cloned$major_c = substr(cloned$c_call,1,4)

cloned = filter(cloned, !is.na(germline_alignment_d_mask))

cloned$c_call = factor(cloned$c_call,
	levels = c("IGHM","IGHD","IGHG3","IGHG1","IGHA1",
		"IGHG2","IGHG4","IGHE","IGHA2"))

cloned$major_c = factor(cloned$major_c,
  levels = c("IGHM","IGHD","IGHG","IGHA","IGHE"))

cloned$SType = cloned$Type
cloned$SType[cloned$Type == "SpositiveB"] = "S1+"
cloned$SType[cloned$Type == "SnegativeB"] = "S1-"

# Cluster calibration plots
pdf("results/cluster_plots.pdf",height = 9.5, width =8.5, useDingbats=FALSE)
g1 = ggplot(cloned,aes(x=factor(seurat_clusters),y=mu_freq))+
    geom_boxplot(outlier.shape=NA)+
    geom_jitter(width=0.1,height=0,size=0.3,aes(color=Type))+theme_bw()+
    xlab("Cluster")+ylab("SHM frequency")+
    theme(axis.text.x = element_text(angle = 0, vjust = 0.5, 
      color = "black"), axis.text.y = element_text(color = "black"))+
    theme(legend.key.size = unit(0.5,"line"))
  
cluster_isotypes = cloned %>%
      filter(!is.na(major_c)) %>%
    group_by(seurat_clusters,major_c, .drop=FALSE) %>%
    summarize(n = n()) %>%
    mutate(Frequency = n/sum(n))
  
g2 = ggplot(filter(cluster_isotypes),
    aes(x=factor(seurat_clusters),y=Frequency,fill=major_c))+
    geom_bar(stat='identity',color="black",size=0.3,width=1)+
    scale_fill_brewer(palette="Set1",name="Isotype")+theme_bw()+
    ylab("Proportion of cells")+
    xlab("Cluster")+
    theme(axis.text.x = element_text(angle = 0, vjust = 0.5, 
      color = "black"), axis.text.y = element_text(color = "black"))+
    theme(legend.key.size = unit(0.5,"line"))

cluster_tissue = cloned %>%
    group_by(seurat_clusters,Tissue) %>%
    summarize(n = n()) %>%
    mutate(Frequency = n/sum(n))  

g3 = ggplot(filter(cluster_tissue),
    aes(x=factor(seurat_clusters),y=Frequency,fill=Tissue))+
    geom_bar(stat='identity',color="black",size=0.3,width=1)+theme_bw()+
    ylab("Proportion of cells")+
    xlab("Cluster")+
    theme(axis.text.x = element_text(angle = 0, vjust = 0.5, 
      color = "black"), axis.text.y = element_text(color = "black"))+
    theme(legend.key.size = unit(0.5,"line"))  

cluster_donor = cloned %>%
    group_by(seurat_clusters,Donor) %>%
    summarize(n = n()) %>%
    mutate(Frequency = n/sum(n))  

g4 = ggplot(filter(cluster_donor),
    aes(x=factor(seurat_clusters),y=Frequency,fill=Donor))+
    geom_bar(stat='identity',color="black",size=0.3,width=1)+theme_bw()+
    ylab("Proportion of cells")+
    xlab("Cluster")+
    theme(axis.text.x = element_text(angle = 0, vjust = 0.5, 
      color = "black"), axis.text.y = element_text(color = "black"))+
    theme(legend.key.size = unit(0.5,"line"))


cluster_type = cloned %>%
    group_by(seurat_clusters,Type) %>%
    summarize(n = n()) %>%
    mutate(Frequency = n/sum(n))  

g5 = ggplot(filter(cluster_type),
    aes(x=factor(seurat_clusters),y=Frequency,fill=Type))+
    geom_bar(stat='identity',color="black",size=0.3,width=1)+theme_bw()+
    ylab("Proportion of cells")+
    xlab("Cluster")+
    theme(axis.text.x = element_text(angle = 0, vjust = 0.5, 
      color = "black"), axis.text.y = element_text(color = "black"))+
    theme(legend.key.size = unit(0.5,"line"))

grid.arrange(g1,g2,g3,g4,g5,ncol=2)

cluster_type_2 = cloned %>%
    group_by(Type, seurat_clusters) %>%
    summarize(n = n()) %>%
    mutate(Frequency = n/sum(n))  

g6 = ggplot(filter(cluster_type_2),
    aes(y=Frequency,fill=factor(seurat_clusters),x=Type))+
    geom_bar(stat='identity',color="black",size=0.3,width=1)+theme_bw()+
    geom_text(aes(label=seurat_clusters),position = position_stack(vjust = 0.5))+
    ylab("Proportion of cells")+
    xlab("Cluster")+
    theme(axis.text.x = element_text(angle = 0, vjust = 0.5, 
      color = "black"), axis.text.y = element_text(color = "black"))+
    theme(legend.key.size = unit(0.5,"line"))+
    labs(fill="Cluster")
  ggsave(g6, file="results/cluster_frequency.pdf",width=4,height=8)

dev.off()

# match sequences against known COVID database
cloned$cdr3_aa = translateDNA(cloned$cdr3)
cloned$cdr3_aal = nchar(cloned$cdr3_aa)
covid = read.csv("CoV-AbDab_200422.csv")
covid = filter(covid, grepl("uman",Heavy.V.Gene))
names(covid)[1] = "Name"
matches = tibble()
for(i in 1:nrow(covid)){
  temp = covid[i,]
  v = getGene(temp$Heavy.V.Gene)
  j = getGene(temp$Heavy.J.Gene)
  cdrh3 = gsub("\\s+","",temp$CDRH3)
  pmatch = filter(cloned, locus == "IGH" &
    grepl(paste0(v,"\\*"), v_call) & 
    grepl(paste0(j,"\\*"), j_call) &
    cdr3_aal == nchar(cdrh3))
  if(nrow(pmatch) > 0){
    pmatch$aa_dist = unlist(lapply(pmatch$cdr3_aa, function(x)
      seqDist(cdrh3, x, getAAMatrix())))/nchar(cdrh3)
    pmatch$match = temp$Name
    matches = bind_rows(matches, pmatch)
  }
}

epitope = covid$Protein...Epitope
names(epitope) = covid$Name

neut = covid$Neutralising.Vs
names(neut) = covid$Name

clone_size = table(filter(cloned, locus=="IGH")$clone_id)

matches_long = matches %>%
  filter(aa_dist <= 0.2) %>%
  select(sequence_id, Donor, Tissue, SType, c_call, seurat_clusters, clone_id, match, aa_dist) %>%
  mutate(clone_size=as.numeric(clone_size[clone_id]), 
    epitope = epitope[match], neutralizing=neut[match]) %>%
  arrange(desc(clone_size), clone_id)

matches_seq = matches_long %>%
  group_by(sequence_id, Donor, Tissue, SType, c_call, seurat_clusters, clone_id, clone_size) %>%
  summarize(
    epitopes = paste(unique(epitope),collapse=","),
    matches = paste(unique(match),collapse=","),
    neutralizing = paste(unique(neutralizing),collapse=",")) %>%
  arrange(desc(clone_size), clone_id)

write.csv(matches_long, file="results/matches_long.csv",row.names=FALSE)
write.csv(matches_seq,  file="results/matches.csv",row.names=FALSE)
matches_long = read.csv("results/matches_long.csv")
matches_seq =  read.csv("results/matches.csv")

matches_seq %>%
  filter(SType == "S1+") %>%
  group_by(c_call) %>%
  summarize(n=n())

matches_seq %>%
  filter(SType == "S1+") %>%
  group_by(seurat_clusters) %>%
  summarize(n=n())

matches_seq %>%
  filter(SType == "S1+") %>%
  group_by(epitopes) %>%
  summarize(n=n())

n_distinct(filter(matches_long, clone_id == "CNMC71-2781")$match)

matches_long %>%
  filter(clone_id == "CNMC71-2781") %>%
  group_by(epitope) %>%
  summarize(n=n())

contingency = cloned %>%
  filter(locus=="IGH") %>%
  mutate(match=sequence_id %in% filter(matches, aa_dist <= 0.2)$sequence_id) %>%
  group_by(SType, match) %>%
  summarize(Cells=n_distinct(cell_id))

write.csv(contingency, file="results/S1vsMatch.csv",row.names=FALSE)

spcon = contingency %>% spread(SType, match)
chisq.test(as.matrix(contingency[,2:3]))

sc = filter(cloned, seurat_clusters==9)

counts = cloned %>%
  group_by(seurat_clusters, SType, Donor) %>%
  summarize(n=n()) %>%
  group_by(seurat_clusters) %>%
  spread(Donor, n)

write.csv(counts, file="results/cluster_counts.csv")

boxsize=0.15
bracketsize=0.15
axissize=0.15
panelsize=0.25
sigsize=2
pointsize=0.1
fontsize=9

# Compare mutation frequency across donors, tissues, and cell types
my_comparisons <- list( unique(cloned$SType))
max_shm = max(cloned$mu_freq)
cloned$SType = factor(cloned$SType, levels=c("S1-", "S1+"))

p1 = ggboxplot(filter(cloned, !is.na(mu_freq)),
    x = "SType", y = "mu_freq", outlier.shape=NA,size=boxsize) +
  geom_jitter(width = 0.1, size= pointsize, height=0, aes(color=major_c))+
  stat_compare_means(comparisons = my_comparisons,
   method="wilcox",method.args=list(correct=FALSE),
   size=sigsize,bracket.size=bracketsize)+
  ylab("SHM frequency")+
  xlab("Sorted B cell population")+
  facet_grid(Tissue~Donor)+
  theme_bw()+
  theme(text=element_text(size=fontsize))+
  theme(strip.background =element_rect(fill="white"))+
  ylim(0,max_shm+0.04)+
  labs(color="Isotype")+
  guides(colour = guide_legend(override.aes = list(size=2)))+
  scale_color_manual(values=isotype_palette)

p2 = ggboxplot(filter(cloned, major_c=="IGHG"),
    x = "SType", y = "mu_freq", outlier.shape=NA,size=boxsize) +
  geom_jitter(width = 0.1, size= pointsize, height=0, aes(color=major_c))+
  stat_compare_means(comparisons = my_comparisons,
   method="wilcox",method.args=list(correct=FALSE),
   size=sigsize,bracket.size=bracketsize)+
  ylab("SHM frequency")+
  xlab("Sorted B cell population")+
  facet_grid(Tissue~Donor)+
  theme_bw()+
  theme(text=element_text(size=fontsize))+
  theme(strip.background =element_rect(fill="white"))+
  ylim(0,max_shm+0.04)+
  labs(color="Isotype")+
  guides(colour = guide_legend(override.aes = list(size=2)))+
  scale_color_manual(values=isotype_palette)

my_comparisons <- list(c("Tonsil", "Adenoid"), c("Adenoid","PBMC"),
  c("Tonsil", "PBMC"))

cloned$Tissue = factor(cloned$Tissue, levels=c("Adenoid","Tonsil","PBMC"))

p3 = ggboxplot(filter(cloned, !is.na(mu_freq)),
    x = "Tissue", y = "mu_freq", outlier.shape=NA,size=boxsize) +
  geom_jitter(width = 0.1, size= pointsize, height=0, aes(color=major_c))+
  stat_compare_means(comparisons = my_comparisons,
   method="wilcox",method.args=list(correct=FALSE),
   size=sigsize,bracket.size=bracketsize)+
  ylab("SHM frequency")+
  xlab("Sorted B cell population")+
  facet_grid(SType~Donor)+
  theme_bw()+
  theme(text=element_text(size=fontsize))+
  theme(strip.background =element_rect(fill="white"))+
  ylim(0,max_shm+0.1)+
  labs(color="Isotype")+
  guides(colour = guide_legend(override.aes = list(size=2)))+
  scale_color_manual(values=isotype_palette)+
  theme(axis.text.x = element_text(angle = 90,vjust = 1, hjust=1))

p4 = ggboxplot(filter(cloned, major_c=="IGHG"),
    x = "Tissue", y = "mu_freq", outlier.shape=NA,size=boxsize) +
  geom_jitter(width = 0.1, size= pointsize, height=0, aes(color=major_c))+
  stat_compare_means(comparisons = my_comparisons,
   method="wilcox",method.args=list(correct=FALSE),
   size=sigsize,bracket.size=bracketsize)+
  ylab("SHM frequency")+
  xlab("Sorted B cell population")+
  facet_grid(SType~Donor)+
  theme_bw()+
  theme(text=element_text(size=fontsize))+
  theme(strip.background =element_rect(fill="white"))+
  ylim(0,max_shm+0.1)+
  labs(color="Isotype")+
  guides(colour = guide_legend(override.aes = list(size=2)))+
  scale_color_manual(values=isotype_palette)+
  theme(axis.text.x = element_text(angle = 90,vjust = 1, hjust=1))

pdf("results/shm.pdf",width=4,height=3,useDingbats=FALSE)
print(p1)
print(p2)
print(p3)
print(p4)
dev.off()

# plot distribution of isotypes
isotypes = cloned %>%
  ungroup() %>%
  group_by(Donor,Tissue,SType,major_c,.drop=FALSE) %>%
  summarize(n=n()) %>%
  mutate(proportion = n/sum(n))

subisotypes = cloned %>%
  ungroup() %>%
  group_by(Donor,Tissue,SType,c_call,.drop=FALSE) %>%
  summarize(n=n()) %>%
  mutate(proportion = n/sum(n))

labsize=3
subisotypes$plot_n = subisotypes$n
subisotypes$plot_n[subisotypes$proportion < 0.15] = ""
isotypes$plot_n = isotypes$n
isotypes$plot_n[isotypes$proportion < 0.15] = ""
p1 = ggplot(subisotypes, aes(x=SType,y=proportion,fill=c_call,label=plot_n))+
  geom_bar(stat='identity',color="black")+
  facet_grid(Tissue ~Donor)+
  scale_fill_manual(values=subisotype_palette)+
  geom_text(position = position_stack(vjust = 0.5),size=labsize)+
  ylab("Proportion of cells")+
  theme_bw()+
  theme(text=element_text(size=fontsize))+
  theme(strip.background =element_rect(fill="white"))+
  labs(fill="Subisotype")

p2 = ggplot(isotypes, aes(x=SType,y=proportion,fill=major_c,label=plot_n))+
  geom_bar(stat='identity',color="black")+
  facet_grid(Tissue ~Donor)+
  scale_fill_manual(values=isotype_palette)+
  geom_text(position = position_stack(vjust = 0.5),size=labsize)+
  ylab("Proportion of cells")+
  theme_bw()+
  theme(text=element_text(size=fontsize))+
  theme(strip.background =element_rect(fill="white"))+
  labs(fill="Isotype")

pdf("results/isotypes.pdf",width=5,height=4,useDingbats=FALSE)
print(p1)
print(p2)
dev.off()

# Clonal diversity analysis of B cell clones
cloned$Donor_Tissue_SType = paste0(cloned$Donor,"-",cloned$Tissue,"-",cloned$SType)
sample_curve_adult <- alphaDiversity(cloned,
        group="Donor_Tissue_SType", clone="clone_id",
    min_q=1.9, max_q=2.1, step_q=0.1,
    ci=0.95, nboot=1000, uniform=TRUE,min_n=100)

saveRDS(sample_curve_adult,file="results/sample_curve.rds")

aq2 = filter(sample_curve_adult@diversity,q==2)
aq2$Donor = unlist(lapply(strsplit(aq2$Donor_Tissue_SType,split="-"),function(x)x[1]))
aq2$Tissue= unlist(lapply(strsplit(aq2$Donor_Tissue_SType,split="-"),function(x)x[2]))
aq2$SType  = unlist(lapply(strsplit(aq2$Donor_Tissue_SType,split="-"),function(x)x[3]))

aq2$SType[aq2$SType == "S1"] = "S1-"
aq2$COVID = "COVID+"
aq2$COVID[aq2$Donor == "CNMC99"] = "COVID-"
aq2$Population = paste0(aq2$COVID, "\n", aq2$SType)
aq2$Population = factor(aq2$Population, levels=c("COVID-\nS1-", "COVID+\nS1-","COVID+\nS1+"))

my_comparisons <- list( c("COVID+\nS1+", "COVID+\nS1-"), 
  c("COVID+\nS1+", "COVID-\nS1-"),
  c("COVID+\nS1-", "COVID-\nS1-"))

shapes = c("CNMC99"=21, "CNMC71"=24, "CNMC89"=25)
boxsize=0.15
bracketsize=0.15
axissize=0.15
panelsize=0.25
sigsize=2
pointsize=1
fontsize=9

max_d = max(aq2$d)
min_d = min(aq2$d)
p1 = ggboxplot(filter(aq2),
    x = "Population", y = "d", outlier.shape=NA,size=boxsize) +
 geom_jitter(size= pointsize, 
    aes(fill=Tissue),shape=21,
    position=position_jitter(seed=4,width = 0.2, height = 0))+
  stat_compare_means(comparisons = my_comparisons,
   method="wilcox",method.args=list(correct=FALSE),
   size=sigsize,bracket.size=bracketsize,step.increase=0.11)+
  ylab("Simpson's diversity")+
  xlab("Sorted B cell population")+
  theme_bw()+
  theme(text=element_text(size=fontsize))+
  theme(strip.background =element_rect(fill="white"))+
  ylim(min_d-0.5,max_d+5)+
  scale_fill_manual(values=tissue_palette)+
  guides(fill = guide_legend(override.aes = list(size=1)))

pdf("results/clonal_diversity.pdf", width=3,height=2, useDingbats=FALSE)
print(p1)
dev.off()

# Calculate clonal overlap among donors
counts = cloned %>%
  group_by(Donor, clone_id) %>%
  mutate(positive = sum(Type=="SpositiveB") > 0, seqs=n()) %>%
  group_by(Donor,positive,clone_id,Tissue,.drop=FALSE) %>%
  summarize(n=n())

results= tibble()
for(d in unique(counts$Donor)){
  for(p in c(TRUE,FALSE)){
    for(a in unique(counts$Tissue)){
        for(b in unique(counts$Tissue)){
          clone_counts = counts %>%
            filter(Donor==d & positive==p & Tissue %in% c(a,b)) %>%
            group_by(positive,clone_id) %>%
            summarize(distinct = n_distinct(Tissue)) 
    
          intersect_clones = clone_counts %>%
            filter(distinct == 2) %>%
            pull(clone_id)
    
          union_clones = clone_counts %>%
            pull(clone_id)
    
          results = bind_rows(results,
            tibble(a=a,b=b,Donor=d,Spositive=p,
              intersect=n_distinct(intersect_clones),
              union=n_distinct(union_clones),
              jaccard=intersect/union))
       }
     }
  }
}

results[results$a==results$b,]$jaccard = NA
results[results$a==results$b,]$intersect = results[results$a==results$b,]$union

results$a = factor(results$a, levels=c("Tonsil","Adenoid","PBMC"))
results$b = factor(results$b, levels=c("Tonsil","Adenoid","PBMC"))

results$SType = ""
results$SType[results$Spositive == TRUE] = "S1+"
results$SType[results$Spositive == FALSE] = "S1-"

boxsize=0.15
bracketsize=0.15
axissize=0.15
panelsize=0.25
sigsize=2
pointsize=1
fontsize=9

pdf("results/clonal_overlap.pdf",width=4,height=4,useDingbats=FALSE)
ggplot(results,aes(x=a,y=b,fill=jaccard))+geom_tile()+
  facet_grid(Donor ~ SType)+
  scale_fill_distiller(palette="RdYlBu")+
  geom_text(aes(label=signif(intersect)),size=3)+
  labs(fill="Jaccard\nindex") + #ggtitle("Normalized clonal overlap")+
  xlab("")+ylab("")+
  ylab("")+
  xlab("")+
  theme_bw()+
  theme(text=element_text(size=fontsize))+
  theme(strip.background =element_rect(fill="white"))+
  theme(axis.text.x = element_text(angle = 90,vjust = 1, hjust=1))
dev.off()

# Read in shared cell cluster information
ts = read.csv("shared_s1_clone_cells_withCITEclust.csv")
m = match(cloned$cell_id, ts$cell_id)
cloned$shared_cluster = ts$Cluster[m]
cloned$shared_cluster[is.na(cloned$shared_cluster)] = ""

# build phylogenetic trees
clones = formatClones(cloned, 
  trait=c("Type","Tissue","c_call","shared_cluster"),
  columns=c("Donor","Tissue","Type"),
  minseq=2)

treesd1 = getTrees(filter(clones, Donor=="CNMC71"),nproc=6,build="igphyml",
  exec="/home/kenneth/Dropbox/Projects/IgPhyML_development/igphyml/src/igphyml")

treesd2 = getTrees(filter(clones, Donor=="CNMC89"),nproc=6,build="igphyml",
  exec="/home/kenneth/Dropbox/Projects/IgPhyML_development/igphyml/src/igphyml")

treesd3 = getTrees(filter(clones, Donor=="CNMC99"),nproc=6,build="igphyml",
  exec="/home/kenneth/Dropbox/Projects/IgPhyML_development/igphyml/src/igphyml")

trees = bind_rows(treesd1, treesd2, treesd3)
trees = trees[order(trees$seqs,decreasing=TRUE),]

saveRDS(trees, file="results/trees.rds")
trees = readRDS("results/trees.rds")

# Figure out which trees have sequences that match with the COV-AbDab database
trees$matches = ""
trees$epitopes = ""
for(i in 1:nrow(trees)){
  matches = which(matches_seq$sequence_id %in% 
    trees$data[[i]]@data$sequence_id)
  if(length(matches) > 0){
    print(i)
    abmatches = unique(unlist(strsplit(matches_seq[matches,]$matches,split=",")))
    epitopes = unique(unlist(strsplit(matches_seq[matches,]$epitopes,split=",")))
    trees$matches[i] = paste0(abmatches,collapse=",")
    trees$epitopes[i] = paste0(epitopes,collapse=",")
  }
}

# AA palette for plotting
aapalette = c(
"H" = "#a6cee3",
"K" = "#a6cee3",
"R" = "#a6cee3",
"D" = "#fb9a99",
"E" = "#fb9a99",
"S" = "#1f78b4",
"T" = "#1f78b4",
"N" = "#1f78b4",
"Q" = "#1f78b4",
"A" = "white",
"V" = "white",
"L" = "white",
"I" = "white",
"M" = "white",
"F" = "#cab2d6",
"Y" = "#cab2d6",
"W" = "#cab2d6",
"P" = "#fb9a99",
"G" = "#fb9a99",
"C" = "#fdbf6f",
"B" = "grey",
"Z" = "grey",
"X" = "grey")

# Make multiple sequence alignment plots 
filter(trees, matches != "")$epitopes

f = filter(cloned, clone_id == "CNMC71-2781" & locus=="IGH")
fasta = paste0(">",f$sequence_id,"\n",
  translateDNA(gsub("\\.","",f$sequence_alignment)))

AA = strsplit(translateDNA(gsub("\\.","",f$sequence_alignment)),split="")
heavy_aa = bind_rows(lapply(1:length(AA),function(x){
  tibble(sequence_id=f$cell_id[x],position=1:length(AA[[x]]),AA=AA[[x]])
}))

ids = strsplit(filter(trees, clone_id == "CNMC71-2781")$matches,split=",")[[1]]

f = filter(covid, Name %in% ids & VH.or.VHH != "ND")
fasta = c(fasta,
  paste0(">",f$Name,"\n",
  f$VH.or.VHH))
writeLines(fasta, con="results/convergent_CNMC71-2781.fasta")

AA = strsplit(f$VH.or.VHH,split="")
heavy_aa = bind_rows(heavy_aa,lapply(1:length(AA),function(x){
  tibble(sequence_id=f$Name[x],position=1:length(AA[[x]]),AA=AA[[x]])
}))

cells = filter(cloned, clone_id == "CNMC71-2781" & locus=="IGH")$cell_id
alldata = readChangeoDb("processed/all_cloned_data.tsv")
f = filter(alldata, cell_id %in% cells & locus != "IGH")
f$sequence_alignment[f$sequence_id == "GTGTTAGGTAGCGTAG-1_contig_1-2712"] = 
  paste0("NN",
    substr(f$sequence_alignment[f$sequence_id == "GTGTTAGGTAGCGTAG-1_contig_1-2712"], 3, 1000))

AA = strsplit(translateDNA(gsub("\\.","",f$sequence_alignment)),split="")
light_aa = bind_rows(lapply(1:length(AA),function(x){
  tibble(sequence_id=f$cell_id[x],position=1:length(AA[[x]]),AA=AA[[x]])
}))

fasta = paste0(">",f$sequence_id,"\n",
  translateDNA(gsub("\\.","",f$sequence_alignment)))
f = filter(covid, Name %in% ids & VL != "ND")
fasta = c(fasta,
  paste0(">",f$Name,"\n",
  f$VL))
writeLines(fasta, con="results/convergent_CNMC71-2781_light.fasta")

AA = strsplit(f$VL,split="")
light_aa = bind_rows(light_aa,lapply(1:length(AA),function(x){
  tibble(sequence_id=f$Name[x],position=1:length(AA[[x]]),AA=AA[[x]])
}))

heavy_aa$sequence_id = factor(heavy_aa$sequence_id,
  levels=unique(heavy_aa$sequence_id))
light_aa$sequence_id = factor(light_aa$sequence_id,
  levels=unique(light_aa$sequence_id))

pdf("results/convergent_clone_aa.pdf",width=9.5,height=5,useDingbats=FALSE)
print(
  ggplot(heavy_aa, aes(x=position,y=sequence_id, fill=AA, label=AA)) + geom_tile() + 
  geom_text(size=2) + 
  theme_bw() + theme(legend.position = "none") + xlab("Position") + ylab("") +
  scale_x_continuous(expand=c(0,0.1)) +
  theme(text = element_text(size = 9), axis.text =element_text(size = 9)) +
  scale_fill_manual(values=aapalette) +
  scale_y_discrete(limits=rev)+
  ggtitle("Clone CNMC71-2781 heavy chain + matches")
  )
print(
  ggplot(light_aa, aes(x=position,y=sequence_id, fill=AA, label=AA)) + geom_tile() + 
  geom_text(size=2) + 
  theme_bw() + theme(legend.position = "none") + xlab("Position") + ylab("") +
  scale_x_continuous(expand=c(0,0.1)) +
  theme(text = element_text(size = 9), axis.text =element_text(size = 9)) +
  scale_fill_manual(values=aapalette) +
  scale_y_discrete(limits=rev)+
  ggtitle("Clone CNMC71-2781 light chain + matches")
  )
dev.off()

# Make tree plots
pal = c(tissue_palette, "Germline"="black")
p = plotTrees(trees)
for(i in 1:length(p)){
  print(i)
  p[[i]] = p[[i]] + geom_tippoint(aes(fill=Tissue),pch=21,size=2)+
    geom_tiplab(aes(subset=(Type=="SpositiveB"),label="+"))+
    scale_fill_manual(values=pal)
}
treesToPDF(p,file="results/alltrees.pdf")

trees$data = lapply(trees$data, function(x){
  x@data$iso_cluster = ""
  x@data$iso_cluster[x@data$shared_cluster != ""] = paste0(x@data$c_call,",",x@data$shared_cluster)
  x@data$iso_cluster[x@data$shared_cluster == ""] = paste0(x@data$c_call)
  x
})

# make S1+ and S1- tree plots
pos = filter(trees, grepl("positive",Type))
neg = filter(trees, !grepl("positive",Type))
maxdivp = max(unlist(lapply(pos[1:3,]$trees,function(x)max(ape::cophenetic.phylo(x)["Germline",]))))
maxdivn = max(unlist(lapply(neg[1:3,]$trees,function(x)max(ape::cophenetic.phylo(x)["Germline",]))))
maxdiv = max(maxdivp, maxdivn)
shapes = c("SpositiveB"=24, "SnegativeB"=21)
ppos = plotTrees(pos[1:3,],scale=0.1)
pneg = plotTrees(neg[1:3,],scale=0.1)
xlim = maxdiv + 0.07
size = 2
offset = 0.015
for(i in 1:3){
  ppos[[i]] = ppos[[i]] + 
  geom_tippoint(aes(fill=Tissue,shape=Type),size=2)+
    scale_shape_manual(values=shapes)+
    scale_fill_manual(values=pal)+
    geom_tiplab(aes(label=iso_cluster),offset=offset,size=size) +
    xlim(0,xlim)

  pneg[[i]] = pneg[[i]] + 
  geom_tippoint(aes(fill=Tissue,shape=Type),size=2)+
    scale_shape_manual(values=shapes)+
    scale_fill_manual(values=pal)+
    geom_tiplab(aes(label=iso_cluster),offset=offset,size=size) +
    xlim(0,xlim)
}

# Combine to make cool figure
treesToPDF(ppos,ncol=1,nrow=3,
  width=3,height=6,useDingbats=FALSE,file="results/positivetrees.pdf")
treesToPDF(pneg,ncol=1,nrow=3,
  width=3,height=6,useDingbats=FALSE,file="results/negativetrees.pdf")
treesToPDF(list(pneg[[1]] + geom_tippoint(aes(fill=Tissue),pch=21,size=2)),ncol=1,nrow=3,
  width=3,height=6,useDingbats=FALSE,file="results/guidetrees.pdf")

maxdiv = max(unlist(lapply(trees$trees,function(x)max(ape::cophenetic.phylo(x)["Germline",]))))
shapes = c("SpositiveB"=24, "SnegativeB"=21)
ppos = plotTrees(trees,scale=0.1)
xlim = maxdiv + 0.07
size = 2
offset = 0.015
for(i in 1:nrow(trees)){
  ppos[[i]] = ppos[[i]] + 
  geom_tippoint(aes(fill=Tissue,shape=Type),size=2)+
    scale_shape_manual(values=shapes)+
    scale_fill_manual(values=pal)+
    geom_tiplab(aes(label=iso_cluster),offset=offset,size=size) +
    xlim(0,xlim)
}
treesToPDF(ppos, file="results/alltrees.pdf")

# plot trees with matches to COV-AbDab database
mtrees = filter(trees, matches != "")
p = plotTrees(mtrees,scale=0.1)
p = lapply(1:length(p), function(x){
 p[[x]] + geom_tippoint(aes(fill=Tissue,shape=Type),size=2)+
    scale_shape_manual(values=shapes)+
    scale_fill_manual(values=pal)+
    geom_tiplab(aes(label=iso_cluster),offset=offset,size=size) +
    xlim(0,xlim)+
    ggtitle(paste0(
      mtrees$clone_id[x],"\nMatches:",
      n_distinct(strsplit(mtrees$matches[x],split=",")[[1]]),"\nEpitopes:",
      paste0(unique(strsplit(mtrees$epitopes[x],split=",")[[1]]),collapse=","))) +
    theme(plot.title = element_text(size=9))+
    theme(text = element_text(size=9))+
    guides(fill = guide_legend(override.aes = list(pch = 24) ))
})
treesToPDF(p,ncol=1,nrow=1,
  width=4,height=3,useDingbats=FALSE,file="results/matchtrees.pdf")

# Analyze the properties of shared clones
shared_clones = counts %>%
  group_by(clone_id) %>%
  filter(sum(Tissue == "Tonsil") > 0 & sum(Tissue == "Adenoid") > 0) %>%
  filter(positive) %>%
  pull(clone_id) %>%
  unique(.)

shared_cells = filter(cloned, clone_id %in% shared_clones)

shared_clone_cells = shared_cells %>%
  select(Donor, clone_id, Tissue, SType, cell_id) %>%
  arrange(Donor, clone_id, Tissue, cell_id, SType) %>% data.frame()

write.csv(shared_clone_cells, file="results/shared_s1_clone_cells.csv", row.names=FALSE, quote=FALSE)

# Compare mutation frequency across donors, tissues, and cell types
cloned$shared = "Not shared"
cloned$shared[cloned$cell_id %in% shared_clone_cells$cell_id] = "Shared"
my_comparisons <- list( unique(cloned$shared))

max_shm = max(filter(cloned, SType == "S1+")$mu_freq)

p1 = ggboxplot(filter(cloned, !is.na(mu_freq) & SType == "S1+"),
    x = "shared", y = "mu_freq", outlier.shape=NA,size=boxsize) +
  geom_jitter(width = 0.1, size= pointsize, height=0, aes(color=major_c))+
  stat_compare_means(comparisons = my_comparisons,
   method="wilcox",method.args=list(correct=FALSE),
   size=sigsize,bracket.size=bracketsize)+
  ylab("SHM frequency, S1+ cells")+
  xlab("Tonsil/adenoid shared clone?")+
  facet_grid(.~Donor)+
  theme_bw()+
  theme(text=element_text(size=fontsize))+
  theme(strip.background =element_rect(fill="white"))+
  ylim(0,max_shm+0.04)+
  labs(color="Isotype")+
  guides(colour = guide_legend(override.aes = list(size=2)))

p2 = ggboxplot(filter(cloned, !is.na(mu_freq) & SType == "S1+" & major_c == "IGHG"),
    x = "shared", y = "mu_freq", outlier.shape=NA,size=boxsize) +
  geom_jitter(width = 0.1, size= pointsize, height=0, aes(color=major_c))+
  stat_compare_means(comparisons = my_comparisons,
   method="wilcox",method.args=list(correct=FALSE),
   size=sigsize,bracket.size=bracketsize)+
  ylab("SHM frequency, S1+ cells")+
  xlab("Tonsil/adenoid shared clone?")+
  facet_grid(.~Donor)+
  theme_bw()+
  theme(text=element_text(size=fontsize))+
  theme(strip.background =element_rect(fill="white"))+
  ylim(0,max_shm+0.04)+
  labs(color="Isotype")+
  guides(colour = guide_legend(override.aes = list(size=2)))

pdf("results/shm_shared.pdf",width=4.5,height=2,useDingbats=FALSE)
print(p1)
print(p2)
dev.off()

# isotype distribution of shared clones across multiple cell types
isotypes = cloned %>%
  ungroup() %>%
  filter(SType == "S1+") %>%
  group_by(Donor,shared,major_c,.drop=FALSE) %>%
  summarize(n=n()) %>%
  mutate(proportion = n/sum(n))

subisotypes = cloned %>%
  ungroup() %>%
  filter(SType == "S1+") %>%
  group_by(Donor,shared,c_call,.drop=FALSE) %>%
  summarize(n=n()) %>%
  mutate(proportion = n/sum(n))

labsize=3
subisotypes$plot_n = subisotypes$n
subisotypes$plot_n[subisotypes$proportion < 0.15] = ""
isotypes$plot_n = isotypes$n
isotypes$plot_n[isotypes$proportion < 0.15] = ""
p1 = ggplot(subisotypes, aes(x=shared,y=proportion,fill=c_call,label=plot_n))+
  geom_bar(stat='identity',color="black")+
  facet_grid(. ~Donor)+
  scale_fill_brewer(palette="Paired")+
  geom_text(position = position_stack(vjust = 0.5),size=labsize)+
  ylab("Proportion of S1+ cells")+
  xlab("")+
  theme_bw()+
  theme(text=element_text(size=fontsize))+
  theme(strip.background =element_rect(fill="white"))+
  labs(fill="Subisotype")

p2 = ggplot(isotypes, aes(x=shared,y=proportion,fill=major_c,label=plot_n))+
  geom_bar(stat='identity',color="black")+
  facet_grid(. ~Donor)+
  scale_fill_brewer(palette="Dark2")+
  geom_text(position = position_stack(vjust = 0.5),size=labsize)+
  ylab("Proportion of S1+ cells")+
  xlab("")+
  theme_bw()+
  theme(text=element_text(size=fontsize))+
  theme(strip.background =element_rect(fill="white"))+
  labs(fill="Isotype")

pdf("results/isotypes_shared.pdf",width=5,height=2.5,useDingbats=FALSE)
print(p1)
print(p2)
dev.off()
