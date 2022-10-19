library(readxl)
library(scales)
library(tidyr)
library(dplyr)
file_name="Cihan_TCR_analysis/20220628_Cihan/Clonotype_abundances_table_oropharyngeal_pediatric_samples_jun28_2022.xlsx"
epitope <- read_excel(file.choose())

##filter table, select data from convalescence CD8 T cells and with single antigen recognition
COVID <- c("CD8_CNMC_71_Adenoid","CD8_CNMC_71_PBMC","CD8_CNMC_71_Tonsil", "CD8_CNMC_89_Adenoid","CD8_CNMC_89_PBMC","CD8_CNMC_89_Tonsil")
epitope_COVID_CD8 <- filter(epitope, Type_Sample %in% COVID)
epitope_COVID_CD8_noNA <- filter(epitope_COVID_CD8, IC_TRB_antigen != "None")
epitope_COVID_CD8_single <- filter(epitope_COVID_CD8_noNA, !grepl(",",IC_TRB_antigen)) 


jColors <- c("#00008B", "#46008B", "#8B008B", "#8B0046" ,"#8B0000", "#8B4500" ,"#8B8B00",
             "#468B00" , "#008B00", "#008B45", "#008B8B" ,"#00468B")

##Calculate the number cells of cells recognize the same antigen
###Remove the package of plyr.It causes problems of using dplyr.
epitope_COVID_CD8_single_among_cells <- epitope_COVID_CD8_single %>% dplyr::group_by(IC_TRB_antigen)%>%summarise(
  Sum = sum(Abundance))

##Calculate the frequencies of cells recognize the same antigen
epitope_COVID_CD8_single_among_cells$Frequency <- epitope_COVID_CD8_single_among_cells$Sum/sum(epitope_COVID_CD8_single_among_cells$Sum)

##plot pie char
ggplot(epitope_COVID_CD8_single_among_cells, aes(x="", y=Frequency,fill=IC_TRB_antigen))+
  scale_fill_manual(values=jColors)+
  geom_bar(stat = "identity",width =1,color="white")+
  coord_polar("y",start = 0)+
  geom_text(aes(x = 1.6,label = percent(Frequency,accuracy = .1)),
            position = position_stack(vjust = .5))+
  theme_void()

ggsave("Distribution of antigen among convalescence CD8 recognize a single antigen.pdf",
       width = 5,
       height =4,
       dpi=300)


##select data from convalescence CD8 T cells that recognize single antigen and large expanded
epitope_COVID_CD8_single_expanded <- filter(epitope_COVID_CD8_single, cloneType_Type_Sample=="Large (0.01 < X <= 0.1)") 
epitope_COVID_CD8_single_expanded_merge<- epitope_COVID_CD8_single_expanded %>% group_by(TRB_1)%>% summarise(
  Sum = sum(Abundance))


##Calculate the frequencies of cells recognize the same antigen
epitope_COVID_CD8_single_expanded_merge$Frequency <- epitope_COVID_CD8_single_expanded_merge$Sum/sum(epitope_COVID_CD8_single_expanded_merge$Sum)
epitope_COVID_CD8_single_expanded_merge$Antigen<- c("S", "ORF1ab","ORF1ab", "S")
ggplot(epitope_COVID_CD8_single_expanded_merge, aes(x="", y=Frequency,fill=Antigen))+
  scale_fill_manual(values=c("#8B0000","#008B45"))+
  geom_bar(stat = "identity",width =1,color="white")+
  coord_polar("y",start = 0)+
  geom_text(aes(x = 1.6,label = percent(Frequency,accuracy = .1)),
             position = position_stack(vjust = .5))+
  theme_void()

ggsave("Proportion of cells from four expanded CD8+ T cell clones which recognize each antigen.pdf",
       width = 5,
       height =4,
       dpi=300)


