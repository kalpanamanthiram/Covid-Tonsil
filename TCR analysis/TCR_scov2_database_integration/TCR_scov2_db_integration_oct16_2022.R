#####CIHAN OGUZ (cihan.oguz@nih.gov). October 2022.
#####NIAID Collaborative Bioinformatics Resource (NCBR), NIAID/NIH

setwd("/Users/oguzc/Downloads/CITEseq_TCR_final")

library(scRepertoire)
library(Seurat)
library(readxl)
library(stringr)
library(writexl)


Human_CD8_CD4_integrated<-readRDS("Human_CD8_CD4_integrated_cDNA_VDJ.rds")


df_metadata<-Human_CD8_CD4_integrated@meta.data


##########
##########
df_metadata$TRA<-stringr::str_split(df_metadata$CTaa, "_", simplify =T)[,1]

df_metadata$TRB<-stringr::str_split(df_metadata$CTaa, "_", simplify =T)[,2]

df_metadata$TRB_1<-gsub('(.*);\\w+', '\\1', df_metadata$TRB)

df_metadata$TRB_2<-stringr::str_split(df_metadata$TRB, ";", simplify =T)[,2]
##########
##########
epitope_antigen_CDR3_trios_class1<-readRDS("ImmuneCODE_epitope_antigen_CDR3_class1.rds")

epitope_antigen_CDR3_trios_class2<-readRDS("ImmuneCODE_epitope_antigen_CDR3_class2.rds")

epitope_antigen_CDR3_trios_class1and2_combined<-rbind(epitope_antigen_CDR3_trios_class1,epitope_antigen_CDR3_trios_class2)

epitope_antigen_CDR3_trios_class1and2_combined<-unique(epitope_antigen_CDR3_trios_class1and2_combined)

########
########
immunecode_uniq_cdr3<-readRDS("ImmuneCODE_CDR3_all.rds")
######
######
vdjdb_all_scov2<-readRDS("VDJdb_epitope_antigen_CDR3.rds")
colnames(vdjdb_all_scov2)<-make.names(colnames(vdjdb_all_scov2))
vdjdb_TRB_all_scov2<-vdjdb_all_scov2[which(vdjdb_all_scov2$Gene=="TRB"),]
vdjdb_TRA_all_scov2<-vdjdb_all_scov2[which(vdjdb_all_scov2$Gene=="TRA"),]
######
######

###########
###########
ind_ic_TRB_1<-which(df_metadata$TRB_1 %in% immunecode_uniq_cdr3)

df_metadata$IC_hit_TRB_1<-"-"
df_metadata$IC_hit_TRB_1[ind_ic_TRB_1]<-"Yes"

ind_vdjdb_TRB_1<-which(df_metadata$TRB_1 %in% vdjdb_TRB_all_scov2$CDR3)

df_metadata$VDJDB_hit_TRB_1<-"-"
df_metadata$VDJDB_hit_TRB_1[ind_vdjdb_TRB_1]<-"Yes"


ind_vdjdb_TRA<-which(df_metadata$TRA %in% vdjdb_TRA_all_scov2$CDR3)

df_metadata$VDJDB_hit_TRA<-"-"
df_metadata$VDJDB_hit_TRA[ind_vdjdb_TRA]<-"Yes"

###########
###########

ind_ic_TRB_2<-which(df_metadata$TRB_2 %in% immunecode_uniq_cdr3)

df_metadata$IC_hit_TRB_2<-"-"
df_metadata$IC_hit_TRB_2[ind_ic_TRB_2]<-"Yes"

ind_vdjdb_TRB_2<-which(df_metadata$TRB_2 %in% vdjdb_TRB_all_scov2$CDR3)

df_metadata$VDJDB_hit_TRB_2<-"-"
df_metadata$VDJDB_hit_TRB_2[ind_vdjdb_TRB_2]<-"Yes"

###########
###########
############
############
###############

df_metadata_TRA<-df_metadata[-union(which(is.na(df_metadata$TRA)==TRUE),which((df_metadata$TRA=="NA"))),]
abundance_TRA<-df_metadata_TRA %>% dplyr::group_by(TRA,Type_Sample) %>% dplyr::summarise(TRA_Abundance_Type_Sample= dplyr::n())
abundance_TRA<-abundance_TRA %>% dplyr::group_by(Type_Sample) %>%
dplyr::mutate(TRA_Frequency_Type_Sample = TRA_Abundance_Type_Sample/sum(TRA_Abundance_Type_Sample))
abundance_TRA<-unique(abundance_TRA)
abundance_TRA<-abundance_TRA[order(-abundance_TRA$TRA_Frequency_Type_Sample),]


df_metadata_TRB<-df_metadata[-union(which(df_metadata$TRB==""),which((df_metadata$TRB=="NA"))),]
abundance_TRB<-df_metadata_TRB %>% dplyr::group_by(TRB,Type_Sample) %>% dplyr::summarise(TRB_Abundance_Type_Sample= dplyr::n())
abundance_TRB<-abundance_TRB %>% dplyr::group_by(Type_Sample) %>%
dplyr::mutate(TRB_Frequency_Type_Sample = TRB_Abundance_Type_Sample/sum(TRB_Abundance_Type_Sample))
abundance_TRB<-unique(abundance_TRB)
abundance_TRB<-abundance_TRB[order(-abundance_TRB$TRB_Frequency_Type_Sample),]
############
############
IC_TRB_hits<-unique(c(df_metadata$TRB_1[which(df_metadata$IC_hit_TRB_1=="Yes")],df_metadata$TRB_2[which(df_metadata$IC_hit_TRB_2=="Yes")]))
#########
#########
VDJDB_TRB_hits<-unique(c(df_metadata$TRB_1[which(df_metadata$VDJDB_hit_TRB_1=="Yes")],df_metadata$TRB_2[which(df_metadata$VDJDB_hit_TRB_2=="Yes")]))
#########
#########
VDJDB_TRA_hits<-unique(c(df_metadata$TRA[which(df_metadata$VDJDB_hit_TRA=="Yes")]))
#########
#########

#########
#########
df_metadata$IC_TRB_hit<-"No"
df_metadata$IC_TRB_epitope<-"None"
df_metadata$IC_TRB_antigen<-"None"
#########
#########


for (IC_TRB_hit_no in 1:length(IC_TRB_hits)){

print(IC_TRB_hit_no)

dum1<-NULL
dum1<-df_metadata$CTaa[which(df_metadata$TRB_1==IC_TRB_hits[IC_TRB_hit_no])]

dum2<-NULL
dum2<-df_metadata$CTaa[which(df_metadata$TRB_2==IC_TRB_hits[IC_TRB_hit_no])]

ind_relev<-NULL
ind_relev<-which(df_metadata$CTaa %in% c(dum1,dum2))

if (length(ind_relev)>0) {
  df_metadata$IC_TRB_hit[ind_relev]<-"Yes"

  uniq_epitope_hits<-NULL

  uniq_epitope_hits<-unique(epitope_antigen_CDR3_trios_class1and2_combined$Epitope[which(epitope_antigen_CDR3_trios_class1and2_combined$CDR3 %in% IC_TRB_hits[IC_TRB_hit_no])])

  df_metadata$IC_TRB_epitope[ind_relev]<-paste(uniq_epitope_hits, collapse = ',')

  uniq_antigen_hits<-NULL

  uniq_antigen_hits<-unique(epitope_antigen_CDR3_trios_class1and2_combined$Antigen[which(epitope_antigen_CDR3_trios_class1and2_combined$CDR3 %in% IC_TRB_hits[IC_TRB_hit_no])])

  df_metadata$IC_TRB_antigen[ind_relev]<-paste(uniq_antigen_hits, collapse = ',')
}


}

#########
#########

#########
#########
df_metadata$VDJDB_TRB_hit<-"No"
df_metadata$VDJDB_TRB_epitope<-"None"
df_metadata$VDJDB_TRB_antigen<-"None"

df_metadata$VDJDB_TRB_ref<-"None"
df_metadata$VDJDB_TRB_MHC.class<-"None"
df_metadata$VDJDB_TRB_MHC.A<-"None"
df_metadata$VDJDB_TRB_MHC.B<-"None"

#########
#########


for (VDJDB_TRB_hit_no in 1:length(VDJDB_TRB_hits)){

print(VDJDB_TRB_hit_no)

dum1<-NULL
dum1<-df_metadata$CTaa[which(df_metadata$TRB_1==VDJDB_TRB_hits[VDJDB_TRB_hit_no])]

dum2<-NULL
dum2<-df_metadata$CTaa[which(df_metadata$TRB_2==VDJDB_TRB_hits[VDJDB_TRB_hit_no])]

ind_relev<-NULL
ind_relev<-which(df_metadata$CTaa %in% c(dum1,dum2))

ind_vdjdb<-union(which(vdjdb_TRB_all_scov2$CDR3 %in% df_metadata$TRB_1[ind_relev]),which(vdjdb_TRB_all_scov2$CDR3 %in% df_metadata$TRB_2[ind_relev]))

if (length(ind_relev)>0) {
  df_metadata$VDJDB_TRB_hit[ind_relev]<-"Yes"

  uniq_epitope_hits<-NULL

  uniq_epitope_hits<-unique(vdjdb_TRB_all_scov2$Epitope[which(vdjdb_TRB_all_scov2$CDR3 %in% VDJDB_TRB_hits[VDJDB_TRB_hit_no])])

  df_metadata$VDJDB_TRB_epitope[ind_relev]<-paste(uniq_epitope_hits, collapse = ',')

  uniq_antigen_hits<-NULL

  uniq_antigen_hits<-unique(vdjdb_TRB_all_scov2$Epitope.gene[which(vdjdb_TRB_all_scov2$CDR3 %in% VDJDB_TRB_hits[VDJDB_TRB_hit_no])])

  df_metadata$VDJDB_TRB_antigen[ind_relev]<-paste(uniq_antigen_hits, collapse = ',')

###############
###############
  df_metadata$VDJDB_TRB_MHC.class[ind_relev]<-paste(unique(vdjdb_TRB_all_scov2$MHC.class[ind_vdjdb]), collapse = ',')

  df_metadata$VDJDB_TRB_MHC.A[ind_relev]<-paste(unique(vdjdb_TRB_all_scov2$MHC.A[ind_vdjdb]), collapse = ',')

  df_metadata$VDJDB_TRB_MHC.B[ind_relev]<-paste(unique(vdjdb_TRB_all_scov2$MHC.B[ind_vdjdb]), collapse = ',')

  df_metadata$VDJDB_TRB_ref[ind_relev]<-paste(unique(vdjdb_TRB_all_scov2$Reference[ind_vdjdb]), collapse = ',')
###############
###############


}


}

#########
#########

#########
#########
df_metadata$VDJDB_TRA_hit<-"No"
df_metadata$VDJDB_TRA_epitope<-"None"
df_metadata$VDJDB_TRA_antigen<-"None"

df_metadata$VDJDB_TRA_ref<-"None"
df_metadata$VDJDB_TRA_MHC.class<-"None"
df_metadata$VDJDB_TRA_MHC.A<-"None"
df_metadata$VDJDB_TRA_MHC.B<-"None"
#########
#########


for (VDJDB_TRA_hit_no in 1:length(VDJDB_TRA_hits)){

print(VDJDB_TRA_hit_no)

dum1<-NULL
dum1<-df_metadata$CTaa[which(df_metadata$TRA==VDJDB_TRA_hits[VDJDB_TRA_hit_no])]


ind_relev<-NULL
ind_relev<-which(df_metadata$CTaa %in% c(dum1))

ind_vdjdb<-which(vdjdb_TRA_all_scov2$CDR3 %in% df_metadata$TRA[ind_relev])

if (length(ind_relev)>0) {
  df_metadata$VDJDB_TRA_hit[ind_relev]<-"Yes"

  uniq_epitope_hits<-NULL

  uniq_epitope_hits<-unique(vdjdb_TRA_all_scov2$Epitope[which(vdjdb_TRA_all_scov2$CDR3 %in% VDJDB_TRA_hits[VDJDB_TRA_hit_no])])

  df_metadata$VDJDB_TRA_epitope[ind_relev]<-paste(uniq_epitope_hits, collapse = ',')

  uniq_antigen_hits<-NULL

  uniq_antigen_hits<-unique(vdjdb_TRA_all_scov2$Epitope.gene[which(vdjdb_TRA_all_scov2$CDR3 %in% VDJDB_TRA_hits[VDJDB_TRA_hit_no])])

  df_metadata$VDJDB_TRA_antigen[ind_relev]<-paste(uniq_antigen_hits, collapse = ',')

  ###############
  ###############
    df_metadata$VDJDB_TRA_MHC.class[ind_relev]<-paste(unique(vdjdb_TRA_all_scov2$MHC.class[ind_vdjdb]), collapse = ',')

    df_metadata$VDJDB_TRA_MHC.A[ind_relev]<-paste(unique(vdjdb_TRA_all_scov2$MHC.A[ind_vdjdb]), collapse = ',')

    df_metadata$VDJDB_TRA_MHC.B[ind_relev]<-paste(unique(vdjdb_TRA_all_scov2$MHC.B[ind_vdjdb]), collapse = ',')

    df_metadata$VDJDB_TRA_ref[ind_relev]<-paste(unique(vdjdb_TRA_all_scov2$Reference[ind_vdjdb]), collapse = ',')
  ###############
  ###############
}


}

#########
#########

df_metadata$cell_barcode<-rownames(df_metadata)

df_metadata$Frequency<-NULL
df_metadata$cloneType<-NULL
df_metadata$Type_Group<-NULL


df_metadata_final<-merge(abundance_TRA,df_metadata,by=c("Type_Sample","TRA"),all = TRUE)
df_metadata_final<-merge(df_metadata_final,abundance_TRB,by=c("Type_Sample","TRB"),all = TRUE)


df_metadata_final<-df_metadata_final[,c(59,1:3,6:28,30:58,4:5,60:61)]

df_metadata_final<-df_metadata_final[,-c(28)]

df_metadata_final<-df_metadata_final[order(-df_metadata_final$TRB_Abundance_Type_Sample),]

write_xlsx(list(metadata_with_TRA_TRB=df_metadata_final,TRA_abundance_vals=abundance_TRA,TRB_abundance_vals=abundance_TRB),"TCR_VDJDB_TRA_TRB_metadata_abundance_pediatric_samples_oct16_2022.xlsx")
