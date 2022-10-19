##read the excel file of assembled data
library(tidyverse)
library(readxl)
library(stringi)
library(stringr)
library(dplyr)

####Ref1: ""PMID:35026743; Ref2: PMID:33296686,Ref3: PMID:33096020; Ref4: PMID:34647971.
###Path of all the reference tables
Ref1 <- "/Users/xuq6/Documents/Tonsil and adenoid Project/Cihan_TCR_analysis/Paired TCR reference papers/4. Cell_mRNA_vaccine/table S2.xlsx"
Ref2 <-"/Users/xuq6/Documents/Tonsil and adenoid Project/Cihan_TCR_analysis/Paired TCR reference papers/5.IM low avidity CD4,2020/Table S3.xlsx"
Ref3 <-"/Users/xuq6/Documents/Tonsil and adenoid Project/Cihan_TCR_analysis/Paired TCR reference papers/6.Cell regulatory and cytotoxic CD4,2020/Table S4A.xlsx"
Ref4a <-"/Users/xuq6/Documents/Tonsil and adenoid Project/Cihan_TCR_analysis/Paired TCR reference papers/7 JEM expand public cTfh/jem_20211327_tables3 cTfh.xlsx"
Ref4b <- "/Users/xuq6/Documents/Tonsil and adenoid Project/Cihan_TCR_analysis/Paired TCR reference papers/7 JEM expand public cTfh/jem_20211327_tables4_public clone type.xlsx"
###Read the tables
dat1 <- read_excel(path = Ref1)
dat2 <- read_excel(path = Ref2,sheet =2)

###read dat3, only keep reactive to SARS-C0V-2 clones
dat3 <- read_excel(path = Ref3, sheet =2,  col_names = FALSE)
dat3 <- dat3[-c(1:3),]%>%select(, c(1:4))
colnames(dat3) <- dat3[1,]
dat3 <- dat3[-1,]
names(dat3)[names(dat3) == 'Virus-reactivity'] <- "Virus_reactivity"
names(dat3)[names(dat3) == 'CDR3 Aminoacid Sequences'] <- "cdr3aa"
dat3 <- dat3 %>% dplyr::filter(grepl("SARS-CoV-2",dat3$Virus_reactivity))

##read dat4, including tfh and public clones
dat4a <- read_excel(path = Ref4a,sheet =2)
dat4b <- read_excel(path = Ref4b,sheet =2)


#######Further clean the tables
##dat1aa
dat1aa <- select(dat1,c("cdr3a","cdr3b"))

##dat2aa
##remove the healthy cases, covid14-17, covid20-21.There is no covid20 in the table.
dat2_filted <- filter(dat2, !sample %in% c("covid_14","covid_15","covid_16","covid_17","covid_21" ))

dat2<- dat2_filted
##https://stringi.gagolewski.com/rapi/stri_split.html
dat2_TRA_TRB_matrix <- stri_list2matrix(stri_split_fixed(dat2$v_gene, '.',n=3, omit_empty=TRUE))%>%t()
dat2_cdr3_matrix <- stri_list2matrix(stri_split_fixed(dat2$cdr3, '.',n=3, omit_empty=TRUE))%>%t()
dat2_matrix <- cbind(dat2_TRA_TRB_matrix,dat2_cdr3_matrix)

TCR1 <- cbind(substr(dat2_matrix[,1],1,4),dat2_matrix[,4])%>%as.data.frame()
TCR2 <-cbind(substr(dat2_matrix[,2],1,4),dat2_matrix[,5])%>%as.data.frame()
TCR3 <-cbind(substr(dat2_matrix[,3],1,4),dat2_matrix[,6])%>%as.data.frame()

TCR1 <- paste (TCR1$V1,TCR1$V2, sep = "_", collapse = NULL)%>%as.data.frame()
TCR2 <- paste (TCR2$V1,TCR2$V2, sep = "_", collapse = NULL)%>%as.data.frame()
TCR3 <- paste (TCR3$V1,TCR3$V2, sep = "_", collapse = NULL)%>%as.data.frame()

dat2_TCR <- as.data.frame(cbind(TCR1, TCR2, TCR3))
colnames(dat2_TCR) <- c("TCR1","TCR2","TCR3")

dat2_TCR_AA_TRA1 <- str_replace_all(dat2_TCR$TCR1, "TRB", NA_character_)
dat2_TCR_AA_TRB1 <- str_replace_all(dat2_TCR$TCR1, "TRA", NA_character_)
dat2_TCR_AA_TRA2 <- str_replace_all(dat2_TCR$TCR2, "TRB", NA_character_)
dat2_TCR_AA_TRB2 <- str_replace_all(dat2_TCR$TCR2, "TRA", NA_character_)
dat2_TCR_AA_TRA3 <- str_replace_all(dat2_TCR$TCR3, "TRB", NA_character_)
dat2_TCR_AA_TRB3 <- str_replace_all(dat2_TCR$TCR3, "TRA", NA_character_)

dat2_TCR_AA_TRA1 <- str_replace(dat2_TCR_AA_TRA1, "TRAV_", "")
dat2_TCR_AA_TRB1 <- str_replace(dat2_TCR_AA_TRB1, "TRBV_", "")
dat2_TCR_AA_TRA2 <- str_replace(dat2_TCR_AA_TRA2, "TRAV_", "")
dat2_TCR_AA_TRB2 <- str_replace(dat2_TCR_AA_TRB2, "TRBV_", "")
dat2_TCR_AA_TRA3 <- str_replace(dat2_TCR_AA_TRA3, "TRAV_", "")
dat2_TCR_AA_TRB3 <- str_replace(dat2_TCR_AA_TRB3, "TRBV_", "")



dat2aa <- as.data.frame(cbind(dat2_TCR_AA_TRA1,dat2_TCR_AA_TRB1,
                                   dat2_TCR_AA_TRA2,dat2_TCR_AA_TRB2,
                                   dat2_TCR_AA_TRA3,dat2_TCR_AA_TRB3))

##Unique TRAaa in Ref2
dat2aa_TRA <-unique(c(na.omit(dat2_TCR_AA_TRA1),
                                  na.omit(dat2_TCR_AA_TRA2),
                                  na.omit(dat2_TCR_AA_TRA3)))%>%as.data.frame()
colnames(dat2aa_TRA) <- "TRA"

##dat3aa
dat3_TCR_AA <- str_split_fixed(dat3$cdr3aa, ';', n=3)%>%as.data.frame()
## V1 and V2 contains both TRA and TRB, be careful of NA
colnames(dat3_TCR_AA) <- c("TRA_TRB_1", "TRA_TRB_2","TRB")
dat3_TCR_AA <- as.data.frame(dat3_TCR_AA)

dat3_TCR_AA_TRA1 <- str_replace_all(dat3_TCR_AA$TRA_TRB_1, "TRB", NA_character_)
dat3_TCR_AA_TRB1 <- str_replace_all(dat3_TCR_AA$TRA_TRB_1, "TRA", NA_character_)
dat3_TCR_AA_TRA2 <- str_replace_all(dat3_TCR_AA$TRA_TRB_2, "TRB", NA_character_)
dat3_TCR_AA_TRB2 <- str_replace_all(dat3_TCR_AA$TRA_TRB_2, "TRA", NA_character_)

dat3_TCR_AA_TRA1 <- str_replace(dat3_TCR_AA_TRA1, "TRA:", "")
dat3_TCR_AA_TRB1 <- str_replace(dat3_TCR_AA_TRB1, "TRB:", "")
dat3_TCR_AA_TRA2 <- str_replace(dat3_TCR_AA_TRA2, "TRA:", "")
dat3_TCR_AA_TRB2 <- str_replace(dat3_TCR_AA_TRB2, "TRB:", "")
dat3_TCR_AA_TRB3 <- str_replace(dat3_TCR_AA$TRB, "TRB:", "")

dat3_TCR_AA_TCR <-as.data.frame(cbind(dat3_TCR_AA_TRA1,dat3_TCR_AA_TRA2,
                        dat3_TCR_AA_TRB1,  dat3_TCR_AA_TRB2,
                        dat3_TCR_AA_TRB3) )
##Unique TRAaa in Ref3
dat3aa_TRA <-unique(c(na.omit(dat3_TCR_AA_TRA1),
                      na.omit(dat3_TCR_AA_TRA2)))%>%as.data.frame()

sum(is.na(dat3aa_TRA))
colnames(dat3aa_TRA) <- "TRA"
##dat4aa
dat4aa_tfh <-select(dat4a, c("CDR3.alpha","CDR3.beta"))
dat4aa_public <-dat4b[,c(3,6)]
dat4aa_public <- dat4aa_public[-2,]
colnames(dat4aa_public) <- dat4aa_public[1,]
dat4aa_public <- dat4aa_public[-1,]

#####Unique TRAaa in Ref4
dat4aa_TRA <- unique(c(dat4aa_tfh$CDR3.alpha%>% na.omit(),dat4aa_public$CDR3α%>% na.omit()))%>%as.data.frame()
sum(is.na(dat4aa_TRA))
colnames(dat4aa_TRA) <- "TRA"

##read CD4_remove_NA meta.data
###path:/Users/xuq6/Documents/Tonsil and adenoid Project/Cihan_TCR_analysis/TCR_06072022/CD4_TRA_nonNA_meta.csv 

###Get CD4_meta_na_TRA from 10032022_CD4_T_Cell_Clustering_Plotting
Clone_abundance <- table(CD4_meta_na_TRA$VDJDB_hit_TRA)%>%as.data.frame()
Clone_abundance_TRA1 <- table(CD4_TRA_nonNA_meta$TRA1)%>%as.data.frame()
##73 of VDJDB_hit_TRA, none of them are expanded CD4

##############################################################################################################
#################Ref1_hit_TRA
Ref1_hit_TRA <- CD4_meta_na_TRA[CD4_meta_na_TRA$TRA1 %in% dat1aa$cdr3a|CD4_meta_na_TRA$TRA2 %in% dat1aa$cdr3a,]

###############Ref2_hit_TRA
###Test if there are matched clones
Ref2_hit_TRA <- CD4_meta_na_TRA[CD4_meta_na_TRA$TRA1 %in% dat2aa_TRA$TRA|CD4_meta_na_TRA$TRA2 %in% dat2aa_TRA$TRA,]

#################Ref3_hit_TRA
Ref3_hit_TRA <- CD4_meta_na_TRA[CD4_meta_na_TRA$TRA1 %in% dat3aa_TRA$TRA|CD4_meta_na_TRA$TRA2 %in% dat3aa_TRA$TRA,]

#################Ref4_hit_TRA
Ref4_hit_TRA <- CD4_meta_na_TRA[CD4_meta_na_TRA$TRA1 %in% dat4aa_TRA$TRA|CD4_meta_na_TRA$TRA2 %in% dat4aa_TRA$TRA,]

####Merge all ref
Ref_all_hit_TRA <- rbind(Ref1_hit_TRA,Ref2_hit_TRA ,Ref3_hit_TRA ,Ref4_hit_TRA )

Ref_all_hit_TRA <- Ref_all_hit_TRA %>% distinct(.keep_all = TRUE)

###Get ID for the matched
ID_Ref1_matched<- Ref1_hit_TRA[,1]
ID_Ref2_matched<- Ref2_hit_TRA[,1]
ID_Ref3_matched<- Ref3_hit_TRA[,1]
ID_Ref4_matched<- Ref4_hit_TRA[,1]
ID_all_Ref_matched <- Ref_all_hit_TRA[,1]

###Add hits into the meta.data
Ref1_hit_TRAaa <- ifelse(CD4_meta_na_TRA$X %in% ID_Ref1_matched,"pos","neg")
CD4_meta_na_TRA$Ref1_hit_TRAaa <- Ref1_hit_TRAaa

Ref2_hit_TRAaa <- ifelse(CD4_meta_na_TRA$X %in% ID_Ref2_matched,"pos","neg")
CD4_meta_na_TRA$Ref2_hit_TRAaa <- Ref2_hit_TRAaa

Ref3_hit_TRAaa <- ifelse(CD4_meta_na_TRA$X %in% ID_Ref3_matched,"pos","neg")
CD4_meta_na_TRA$Ref3_hit_TRAaa <- Ref3_hit_TRAaa

Ref4_hit_TRAaa <- ifelse(CD4_meta_na_TRA$X %in% ID_Ref4_matched,"pos","neg")
CD4_meta_na_TRA$Ref4_hit_TRAaa <- Ref4_hit_TRAaa


ID_all_Ref_matched_TRAaa <- ifelse(CD4_meta_na_TRA$X %in% ID_all_Ref_matched,"pos","neg")
CD4_meta_na_TRA$ID_all_Ref_matched_TRAaa <- ID_all_Ref_matched_TRAaa

write.csv(CD4_meta_na_TRA,"CD4_new_meta_ref.csv")

