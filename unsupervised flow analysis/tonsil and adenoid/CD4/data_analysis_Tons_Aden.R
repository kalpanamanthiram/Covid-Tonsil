library(nlme)
library(lme4)
library(FlowSOM)
library(flowCore)
library(CytoML)
library(flowWorkspace)
library(data.table)
library(tidyverse)
library(readxl)
library(uwot)
library(Seurat)
#library(leiden)
library(parallel)
library(gridExtra)
library(vroom)
library(ggcyto)

#####
##### !!!!! CONCATENATING IN FLOWJO LEADS TO HUGE BATCH EFFECT (with or without metacells) !!!!!
#####

#####
##### !!!!! AFTER cytoqc (MAYBE NOT BC OF cytoqc THOUGH) CELL COUNTS ARE VERY DIFFERENT FROM FLOWJO !!!!!
##### !!!!! THE ABOVE IS TRUE EVEN IF NEW WSP FILE IS USED !!!!!
##### !!!!! THE ABOVE HAPPENS BC THE FSC-A vs TIME GATE IS MISINTERPRETED IN A FEW SAMPLES by flowjo_to_gatingset!!!!!
##### !!!!! (READING ONE FCS AT TIME DOES NOT SOLVE THE PROBLEM) !!!!!
#####

##### 
##### !!!!! CREATE PLUGIN THAT EXPORTS EventNumberDP FROM FLOWJO !!!!!
##### !!!!! EventNumberDP WILL BE USED TO SELECT CELLS INSTEAD OF GatingSet GATES  !!!!!
#####

##### 
##### !!!!! SURPRISINGLY, THE ABOVE WAS HELPED TO SOLVE THE BATCH EFFECT PROBLEM !!!!!
##### !!!!! SINCE DATA WERE TRANSFORMED IN R (flowWorkspace/CytoML) INSTEAD OF FLOWJO !!!!!
#####

# select relevant features
# see EMAIL: 210509 and "Markers for 37c colour.xlsx"
cluster.features <- c("CCR6", "CXCR5", "CXCR3", "CCR7", "CD45RA",
                      "CD161", "CD28", "PD-1", "CD57", "CD25",
                      "CD95", "CD27", "CD127", "HLA-DR", "CD38",
                      "ICOS", "CD11c", "CD24", "CD1c", "CD123",
                      "CD21") %>%
  unique()

# open FlowJo workspace
wsFile_Tonsil <- "~/Box/My Box Notes/4.Flowjo manual gating/Tonsil_with_labels_exclusion.wsp"
ws_Tonsil <- CytoML::open_flowjo_xml(file = wsFile_Tonsil,
                                     sample_names_from = "sampleNode")
# find path to FCS files
sampleFCS_path_Tonsil <- XML::xmlParse(wsFile_Tonsil) %>%
  XML::xpathApply(.,
                  file.path("/Workspace/SampleList/Sample",
                            "DataSet"),
                  function(x)
                    XML::xmlGetAttr(x,
                                    "uri") %>%
                    gsub(pattern = "%20", 
                         replacement = " ", 
                         x = .) %>%
                    gsub(pattern = "file:",
                         replacement = "",
                         x = .)) %>%
  unlist

# TRYING cytoqc PACKAGE (https://github.com/RGLab/cytoqc)
# DUE TO:
### Error in fcs_to_cytoset(sapply(files, normalizePath), list(which.lines = which.lines,  : 
### Found channel inconsistency across samples. 'SSC-W' is missing from XXX.fcs
cqc_data_Tonsil <- normalizePath(path = sampleFCS_path_Tonsil) %>%
  cytoqc::cqc_load_fcs(files = .)
cqc_data_Tonsil

check_results_Tonsil <- cytoqc::cqc_check(x = cqc_data_Tonsil, 
                                          type = "channel")
check_results_Tonsil

res_Tonsil <- cytoqc::cqc_match(x = check_results_Tonsil, 
                                ref = 2)
res_Tonsil

cytoqc::cqc_fix(x = res_Tonsil)

check_results_Tonsil <- cytoqc::cqc_check(x = cqc_data_Tonsil, 
                                          type = "channel")
check_results_Tonsil

cqc_data_Tonsil <- cytoqc::cqc_get_data(check_results_Tonsil)
cqc_data_Tonsil

cs_Tonsil <- flowWorkspace::cytoset(x = cqc_data_Tonsil)

(cs_Tonsil[[1]] %>% flowCore::parameters() %>% flowCore::pData())

# make GatingSet
gs_Tonsil <- flowjo_to_gatingset(
  ws = ws_Tonsil,
  name = 1,
  subset = basename(sampleFCS_path_Tonsil),
  execute = TRUE,
  path = "",
  cytoset = cs_Tonsil,
  backend_dir = tempdir(),
  backend = get_default_backend(),
  includeGates = TRUE,
  additional.keys = "$TOT",
  additional.sampleID = TRUE,
  keywords = character(),
  keywords.source = "XML",
  keyword.ignore.case = FALSE,
  extend_val = -Inf,
  extend_to = -4000,
  channel.ignore.case = FALSE,
  leaf.bool = TRUE,
  include_empty_tree = FALSE,
  skip_faulty_gate = FALSE,
  compensation = NULL,
  transform = TRUE,
  fcs_file_extension = ".fcs",
  greedy_match = FALSE,
  mc.cores = 1) 

# get file paths of FJ exported files containing idx of gated cells
# THIS WAS EXPORTED FROM FLOWJO USING NEW PLUGIN
EventNumberDP.path_Tonsil <- sub(pattern = ".wsp$",
                                 replacement = "", 
                                 x = wsFile_Tonsil,
                                 fixed = FALSE) %>%
  paste0(.,
         "/EventNumberDP") %>%
  list.files(full.names = TRUE) %>%
  grep(pattern = "CD4+", 
       x = ., 
       ignore.case = FALSE, 
       perl = FALSE, 
       value = TRUE,
       fixed = TRUE, 
       useBytes = FALSE, 
       invert = FALSE)

# read csv files and store idx in object
popOfInt.idx_Tonsil <- sampleNames(gs_Tonsil) %>%
  gsub(pattern = "\\..*",
       replacement = "",
       x = .,
       fixed = FALSE) %>%
  setNames(nm = sampleNames(gs_Tonsil)) %>%
  lapply(X = .,
         FUN = function(fcs){
           EventNumberDP.path_Tonsil[
             basename(EventNumberDP.path_Tonsil) %>%
               startsWith(fcs)
           ] %>%
             readr::read_csv(file = ., 
                             col_types = cols(),
                             progress = FALSE) %>%
             .$EventNumberDP
         })

# sample level metacells
metacells_Tonsil <- sampleNames(gs_Tonsil) %>%
  gsub(pattern = "\\..*",
       replacement = "",
       x = .,
       fixed = FALSE) %>%
  paste0(.,
         ".fcs") %>%
  setNames(object = sampleNames(gs_Tonsil),
           nm = .) %>%
  parallel::mclapply(FUN = function(ff){
    # get cytoframe for sample ff
    cs_ff <- flowWorkspace::gs_pop_get_data(obj = gs_Tonsil[[ff]],
                                            y = "root",
                                            inverse.transform = FALSE)[[ff]]
    # get flowset for sample ff
    fs_ff <- flowWorkspace::cytoframe_to_flowFrame(cs_ff[popOfInt.idx_Tonsil[[ff]],])
    # get protein expression for sample ff
    exprs_ff <- flowCore::exprs(object = fs_ff)
    # change colnames of protein expression matrix from channel to antigen name
    colnames(exprs_ff) <- is.na(flowCore::parameters(gh_pop_get_data(gs_Tonsil[[1]])) %>%
                                  flowCore::pData() %>%
                                  .$desc) %>%
      ifelse(yes = flowCore::parameters(gh_pop_get_data(gs_Tonsil[[1]])) %>%
               flowCore::pData() %>%
               .$name,
             no = flowCore::parameters(gh_pop_get_data(gs_Tonsil[[1]])) %>%
               flowCore::pData() %>%
               .$desc)
    irrel.par <- grepl(pattern = "^fsc|^ssc|^time|^viability|^live|^dead|^zombie|^comp-af-a",
                       x = colnames(exprs_ff),
                       ignore.case = TRUE,
                       fixed = FALSE)
    
    paste0(ff,
           ": excluding ",
           paste0(colnames(exprs_ff)[irrel.par], collapse = ", "),
           "\n") %>%
      message()    
    exprs_ff <- exprs_ff[,!irrel.par]
    ## keep cells gated in FlowJo and parameters of interest
    exprs_ff <- exprs_ff[,cluster.features]
    
    paste0(ff,
           ": using only ",
           paste0(colnames(exprs_ff), collapse = ", "),
           "\n") %>%
      message()    
    flush.console()
    # trim values of extreme outliers (0.1 percentile)
    exprs_ff <- apply(X = exprs_ff,
                      MARGIN = 2,
                      FUN = function(marker){
                        small <- quantile(x = marker,
                                          probs = 1e-3)
                        large <- quantile(x = marker,
                                          probs = (1 - 1e-3))
                        marker[marker < small] <- small
                        marker[marker > large] <- large
                        return(marker)
                      })
    # get and return metacells
    set.seed(1985)
    kmeans.res <- suppressWarnings(stats::kmeans(x = exprs_ff,
                                                 centers = 500,
                                                 iter.max = 10,
                                                 nstart = 1))
    if(kmeans.res$ifault == 4) { # https://stackoverflow.com/a/30055776
      set.seed(1985)
      kmeans.res <- suppressWarnings(stats::kmeans(x = exprs_ff,
                                                   centers = kmeans.res$centers,
                                                   iter.max = 10,
                                                   algorithm = "MacQueen"))
    }
    return(kmeans.res)
  },
  mc.cores = 8)

saveRDS(object = metacells_Tonsil,
        file = "metacells_Tonsil.rds")

# open FlowJo workspace
wsFile_Adenoid <- "~/Box/My Box Notes/4.Flowjo manual gating/37c_Adenoids_T_Cells_after_shutdown_exclusion-35,39,52.wsp"
ws_Adenoid <- CytoML::open_flowjo_xml(file = wsFile_Adenoid,
                                      sample_names_from = "sampleNode")
# find path to FCS files
sampleFCS_path_Adenoid <- XML::xmlParse(wsFile_Adenoid) %>%
  XML::xpathApply(.,
                  file.path("/Workspace/SampleList/Sample",
                            "DataSet"),
                  function(x)
                    XML::xmlGetAttr(x,
                                    "uri") %>%
                    gsub(pattern = "%20", 
                         replacement = " ", 
                         x = .) %>%
                    gsub(pattern = "file:",
                         replacement = "",
                         x = .)) %>%
  unlist

# TRYING cytoqc PACKAGE (https://github.com/RGLab/cytoqc)
# DUE TO:
### Error in fcs_to_cytoset(sapply(files, normalizePath), list(which.lines = which.lines,  : 
### Found channel inconsistency across samples. 'SSC-W' is missing from XXX.fcs
cqc_data_Adenoid <- normalizePath(path = sampleFCS_path_Adenoid) %>%
  cytoqc::cqc_load_fcs(files = .)
cqc_data_Adenoid

check_results_Adenoid <- cytoqc::cqc_check(x = cqc_data_Adenoid, 
                                           type = "channel")
check_results_Adenoid

res_Adenoid <- cytoqc::cqc_match(x = check_results_Adenoid, 
                                 ref = 2)
res_Adenoid

cytoqc::cqc_fix(x = res_Adenoid)

check_results_Adenoid <- cytoqc::cqc_check(x = cqc_data_Adenoid, 
                                           type = "channel")
check_results_Adenoid

cqc_data_Adenoid <- cytoqc::cqc_get_data(check_results_Adenoid)
cqc_data_Adenoid

cs_Adenoid <- flowWorkspace::cytoset(x = cqc_data_Adenoid)

(cs_Adenoid[[1]] %>% flowCore::parameters() %>% flowCore::pData())

# make GatingSet
gs_Adenoid <- flowjo_to_gatingset(
  ws = ws_Adenoid,
  name = 1,
  subset = basename(sampleFCS_path_Adenoid),
  execute = TRUE,
  path = "",
  cytoset = cs_Adenoid,
  backend_dir = tempdir(),
  backend = get_default_backend(),
  includeGates = TRUE,
  additional.keys = "$TOT",
  additional.sampleID = TRUE,
  keywords = character(),
  keywords.source = "XML",
  keyword.ignore.case = FALSE,
  extend_val = -Inf,
  extend_to = -4000,
  channel.ignore.case = FALSE,
  leaf.bool = TRUE,
  include_empty_tree = FALSE,
  skip_faulty_gate = FALSE,
  compensation = NULL,
  transform = TRUE,
  fcs_file_extension = ".fcs",
  greedy_match = FALSE,
  mc.cores = 1) 

# get file paths of FJ exported files containing idx of gated cells
# THIS WAS EXPORTED FROM FLOWJO USING NEW PLUGIN
EventNumberDP.path_Adenoid <- sub(pattern = ".wsp$",
                                  replacement = "", 
                                  x = wsFile_Adenoid,
                                  fixed = FALSE) %>%
  paste0(.,
         "/EventNumberDP") %>%
  list.files(full.names = TRUE) %>%
  grep(pattern = "CD4+", 
       x = ., 
       ignore.case = FALSE, 
       perl = FALSE, 
       value = TRUE,
       fixed = TRUE, 
       useBytes = FALSE, 
       invert = FALSE)

# read csv files and store idx in object
popOfInt.idx_Adenoid <- sampleNames(gs_Adenoid) %>%
  gsub(pattern = "\\..*",
       replacement = "",
       x = .,
       fixed = FALSE) %>%
  setNames(nm = sampleNames(gs_Adenoid)) %>%
  lapply(X = .,
         FUN = function(fcs){
           EventNumberDP.path_Adenoid[
             basename(EventNumberDP.path_Adenoid) %>%
               startsWith(fcs)
           ] %>%
             readr::read_csv(file = ., 
                             col_types = cols(),
                             progress = FALSE) %>%
             .$EventNumberDP
         })

# sample level metacells
metacells_Adenoid <- sampleNames(gs_Adenoid) %>%
  gsub(pattern = "\\..*",
       replacement = "",
       x = .,
       fixed = FALSE) %>%
  paste0(.,
         ".fcs") %>%
  setNames(object = sampleNames(gs_Adenoid),
           nm = .) %>%
  parallel::mclapply(FUN = function(ff){
    # get cytoframe for sample ff
    cs_ff <- flowWorkspace::gs_pop_get_data(obj = gs_Adenoid[[ff]],
                                            y = "root",
                                            inverse.transform = FALSE)[[ff]]
    # get flowset for sample ff
    fs_ff <- flowWorkspace::cytoframe_to_flowFrame(cs_ff[popOfInt.idx_Adenoid[[ff]],])
    # get protein expression for sample ff
    exprs_ff <- flowCore::exprs(object = fs_ff)
    # change colnames of protein expression matrix from channel to antigen name
    colnames(exprs_ff) <- is.na(flowCore::parameters(gh_pop_get_data(gs_Adenoid[[1]])) %>%
                                  flowCore::pData() %>%
                                  .$desc) %>%
      ifelse(yes = flowCore::parameters(gh_pop_get_data(gs_Adenoid[[1]])) %>%
               flowCore::pData() %>%
               .$name,
             no = flowCore::parameters(gh_pop_get_data(gs_Adenoid[[1]])) %>%
               flowCore::pData() %>%
               .$desc)
    irrel.par <- grepl(pattern = "^fsc|^ssc|^time|^viability|^live|^dead|^zombie|^comp-af-a",
                       x = colnames(exprs_ff),
                       ignore.case = TRUE,
                       fixed = FALSE)
    
    paste0(ff,
           ": excluding ",
           paste0(colnames(exprs_ff)[irrel.par], collapse = ", "),
           "\n") %>%
      message()    
    exprs_ff <- exprs_ff[,!irrel.par]
    ## keep cells gated in FlowJo and parameters of interest
    exprs_ff <- exprs_ff[,cluster.features]
    
    paste0(ff,
           ": using only ",
           paste0(colnames(exprs_ff), collapse = ", "),
           "\n") %>%
      message()    
    # trim values of extreme outliers (0.1 percentile)
    exprs_ff <- apply(X = exprs_ff,
                      MARGIN = 2,
                      FUN = function(marker){
                        small <- quantile(x = marker,
                                          probs = 1e-3)
                        large <- quantile(x = marker,
                                          probs = (1 - 1e-3))
                        marker[marker < small] <- small
                        marker[marker > large] <- large
                        return(marker)
                      })
    # get and return metacells
    set.seed(1985)
    kmeans.res <- suppressWarnings(stats::kmeans(x = exprs_ff,
                                                 centers = 500,
                                                 iter.max = 10,
                                                 nstart = 1))
    if(kmeans.res$ifault == 4) { # https://stackoverflow.com/a/30055776
      set.seed(1985)
      kmeans.res <- suppressWarnings(stats::kmeans(x = exprs_ff,
                                                   centers = kmeans.res$centers,
                                                   iter.max = 10,
                                                   algorithm = "MacQueen"))
    }
    return(kmeans.res)
  },
  mc.cores = 8)

saveRDS(object = metacells_Adenoid,
        file = "metacells_Adenoid.rds")

rm(list = ls())

metacells_Tonsil <- readRDS(file = "metacells_Tonsil.rds")

conc_Tonsil <- lapply(X = metacells_Tonsil,
                      FUN = "[",
                      "centers") %>%
  unlist(recursive = FALSE,
         use.names = FALSE) %>%
  do.call(rbind,
          .)

metaSampleID_Tonsil <- names(x = metacells_Tonsil) %>%
  rep(each = 500)

metaSampleID_Tonsil_clean <- names(x = metacells_Tonsil) %>% 
  gsub(pattern = "^Tonsil ",
       replacement = "",
       x = .,
       ignore.case = FALSE,
       perl = FALSE,
       fixed = FALSE,
       useBytes = FALSE) %>%
  gsub(pattern = ".fcs$",
       replacement = "",
       x = .,
       ignore.case = FALSE,
       perl = FALSE,
       fixed = FALSE,
       useBytes = FALSE) %>%
  gsub(pattern = "\\-.*|\\_.*",
       replacement = "",
       x = .,
       ignore.case = FALSE,
       perl = FALSE,
       fixed = FALSE,
       useBytes = FALSE) %>%
  gsub(pattern = "[^[:alnum:]]",
       replacement = "",
       x = .,
       fixed = FALSE) %>%
  rep(each = 500)

metacells_Adenoid <- readRDS(file = "metacells_Adenoid.rds")

conc_Adenoid <- lapply(X = metacells_Adenoid,
                       FUN = "[",
                       "centers") %>%
  unlist(recursive = FALSE,
         use.names = FALSE) %>%
  do.call(rbind,
          .)

metaSampleID_Adenoid <- names(x = metacells_Adenoid) %>%
  rep(each = 500)

metaSampleID_Adenoid_clean <- names(x = metacells_Adenoid) %>% 
  gsub(pattern = "^Adenoid ",
       replacement = "",
       x = .,
       ignore.case = FALSE,
       perl = FALSE,
       fixed = FALSE,
       useBytes = FALSE) %>%
  gsub(pattern = ".fcs$",
       replacement = "",
       x = .,
       ignore.case = FALSE,
       perl = FALSE,
       fixed = FALSE,
       useBytes = FALSE) %>%
  gsub(pattern = "\\-.*|\\_.*",
       replacement = "",
       x = .,
       ignore.case = FALSE,
       perl = FALSE,
       fixed = FALSE,
       useBytes = FALSE) %>%
  gsub(pattern = "[^[:alnum:]]",
       replacement = "",
       x = .,
       fixed = FALSE) %>%
  rep(each = 500)

# concatenate organs
conc <- rbind(conc_Tonsil,
              conc_Adenoid)

metaSampleID <- c(metaSampleID_Tonsil,
                  metaSampleID_Adenoid)

metaSampleID_clean <- c(metaSampleID_Tonsil_clean,
                        metaSampleID_Adenoid_clean)

metaOrganID <- c(rep(x = "Tonsil",
                     times = length(metaSampleID_Tonsil)),
                 rep(x = "Adenoid",
                     times = length(metaSampleID_Adenoid)))

# load metadata
meta <- readxl::read_excel("~/Box/My Box Notes/2.Matadata/Full metadata_061021.xlsx")

relevant_meta <- meta[,c("Participant_ID",
                         "Sex",
                         "Age",
                         "COVID",
                         "Time since infection")]

colnames(relevant_meta) <- c("patient_id",
                             "sex",
                             "age",
                             "status",
                             "dpi")

relevant_meta$patient_id <- relevant_meta$patient_id %>% 
  gsub(pattern = "^CNMC ",
       replacement = "",
       x = .,
       ignore.case = FALSE,
       perl = FALSE,
       fixed = FALSE,
       useBytes = FALSE) %>%
  gsub(pattern = "\\-.*|\\_.*",
       replacement = "",
       x = .,
       ignore.case = FALSE,
       perl = FALSE,
       fixed = FALSE,
       useBytes = FALSE) %>%
  gsub(pattern = "[^[:alnum:]]",
       replacement = "",
       x = .,
       fixed = FALSE) %>%
  gsub(pattern = "^0",
       replacement = "",
       x = .,
       fixed = FALSE)

relevant_meta$status <- ifelse(test = relevant_meta$status == "No",
                               yes = "CONTROL",
                               no = "COVID")

relevant_meta <- relevant_meta[relevant_meta$patient_id %in%
                                 unique(metaSampleID_clean),]

#kmeans.res <- readRDS(file = "kmeans.res_Tons_Aden.rds")

conc_pca <- conc %>%
  irlba::prcomp_irlba(x = .,
                      center = FALSE,
                      scale. = FALSE,
                      n = 2) %>%
  .$x %>%
  as.data.frame()

((conc_pca %>%
    cbind(metaSampleID = metaSampleID)) %>%
    ggplot(mapping = aes(x = PC1,
                         y = PC2,
                         color = as.factor(metaSampleID))) +
    labs(title = "sampleID") +
    theme_bw() +
    theme(legend.position = "none") +
    geom_point()) %>%
  egg::set_panel_size(width = unit(3, "in"),
                      height = unit(3, "in")) %>%
  gridExtra::grid.arrange()

((conc_pca %>%
    cbind(metaOrganID = metaOrganID)) %>%
    ggplot(mapping = aes(x = PC1,
                         y = PC2,
                         color = as.factor(metaOrganID))) +
    labs(color = "organID") +
    theme_bw() +
    geom_point()) %>%
  egg::set_panel_size(width = unit(3, "in"),
                      height = unit(3, "in")) %>%
  gridExtra::grid.arrange()

# load metadata
meta <- readxl::read_excel("~/Box/My Box Notes/2.Matadata/Full metadata_061021.xlsx")

relevant_meta <- meta[,c("Participant_ID",
                         "Sex",
                         "Age",
                         "COVID",
                         "Time since infection")]

colnames(relevant_meta) <- c("patient_id",
                             "sex",
                             "age",
                             "status",
                             "dpi")

relevant_meta$patient_id <- relevant_meta$patient_id %>% 
  gsub(pattern = "^CNMC ",
       replacement = "",
       x = .,
       ignore.case = FALSE,
       perl = FALSE,
       fixed = FALSE,
       useBytes = FALSE) %>%
  gsub(pattern = "\\-.*|\\_.*",
       replacement = "",
       x = .,
       ignore.case = FALSE,
       perl = FALSE,
       fixed = FALSE,
       useBytes = FALSE) %>%
  gsub(pattern = "[^[:alnum:]]",
       replacement = "",
       x = .,
       fixed = FALSE) %>%
  gsub(pattern = "^0",
       replacement = "",
       x = .,
       fixed = FALSE)

relevant_meta$status <- ifelse(test = relevant_meta$status == "No",
                               yes = "CONTROL",
                               no = "COVID")

relevant_meta <- relevant_meta[relevant_meta$patient_id %in%
                                 unique(metaSampleID_clean),]

#relevant_meta <- relevant_meta[relevant_meta$patient_id %in%
#                                 cells_dat$patient_id,]

# meta-metacells
set.seed(1985)
kmeans.res <- suppressWarnings(stats::kmeans(x = conc,
                                             centers = 500,
                                             iter.max = 10))
if(kmeans.res$ifault == 4) { # https://stackoverflow.com/a/30055776
  kmeans.res <- suppressWarnings(stats::kmeans(x = conc,
                                               centers = kmeans.res$centers,
                                               iter.max = 10,
                                               algorithm = "MacQueen"))
  rownames(kmeans.res$centers) <- dim(kmeans.res$centers)[1] %>% 
    seq()
}

n_neighbors <- 5
set.seed(1985)
kmeans.res$graphs <- Seurat::FindNeighbors(
  object = kmeans.res$centers, #cosNorm(kmeans.res$centers), # cosine normalized
  k.param = n_neighbors,
  compute.SNN = TRUE,
  prune.SNN = 1/15,
  nn.method = "rann",
  nn.eps = 0,
  verbose = FALSE,
  force.recalc = FALSE)

kmeans.res$leiden <- leiden::leiden(object = kmeans.res$graphs$snn,
                                    resolution_parameter = 1,
                                    seed = 1985,
                                    n_iterations = 10L) %>%
  as.factor()

table(kmeans.res$leiden) 

kmeans.res$idents <- kmeans.res$leiden[kmeans.res$cluster]

names(kmeans.res$idents) <- metaSampleID

table(kmeans.res$idents) 

sub.pops <- levels(kmeans.res$leiden)

# cluster vs sample
table(kmeans.res$idents,
      metaSampleID_clean)

table(kmeans.res$idents,
      metaSampleID_clean) %>%
  apply(MARGIN = 2,
        FUN = function(msID){
          msID * 100 / sum(msID)
        }) %>%
  as.data.frame()

# outliers
freq.PCA <- table(kmeans.res$idents,
                  metaSampleID) %>%
  apply(MARGIN = 2,
        FUN = function(msID){
          msID * 100 / sum(msID)
        }) %>%
  as.data.frame() %>%
  t() %>%
  irlba::prcomp_irlba(x = .,
                      center = FALSE,
                      scale. = FALSE,
                      n = 2) %>%
  .$x %>%
  as.data.frame()

freq.PCA %>%
  ggplot(mapping = aes(x = PC1,
                       y = PC2)) +
  theme_bw() +
  geom_point()

(freq.PCA %>%
    cbind(sampleID = (table(kmeans.res$idents,
                            metaSampleID) %>% t() %>% 
                        rownames() %>%
                        gsub(pattern = ".fcs",
                             replacement = "",
                             x = .,
                             fixed = TRUE)))) %>%
  ggplot(mapping = aes(x = PC1,
                       y = PC2,
                       label = sampleID)) +
  theme_bw() +
  geom_point() +
  geom_text()

# not major ouliers

(c(names(metacells_Tonsil),
   names(metacells_Adenoid)) %>%
    sort() ==
    colnames(table(kmeans.res$idents,
                   metaSampleID) %>%
               apply(MARGIN = 2,
                     FUN = function(msID){
                       msID * 100 / sum(msID)
                     }) %>%
               as.data.frame())) %>%
  all

(table(kmeans.res$idents,
       metaSampleID) %>%
    apply(MARGIN = 2,
          FUN = function(msID){
            msID * 100 / sum(msID)
          }) %>%
    as.data.frame() %>%
    t() %>%
    reshape::melt() %>%
    "colnames<-"(c("sample",
                   "cluster",
                   "percentage"))) %>%
  ggplot(mapping = aes(y = percentage)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, 
                                   vjust = 1, 
                                   hjust = 1)) +
  geom_bar(mapping = aes(x = as.factor(sample), 
                         fill = as.factor(cluster)), 
           position = "fill",
           stat = "identity")

(table(kmeans.res$idents,
       metaSampleID) %>%
    apply(MARGIN = 2,
          FUN = function(msID){
            msID * 100 / sum(msID)
          }) %>%
    as.data.frame() %>%
    t() %>%
    reshape::melt() %>%
    "colnames<-"(c("sample",
                   "cluster",
                   "percentage"))) %>%
  ggplot(mapping = aes(y = percentage)) +
  theme_bw() +
  geom_bar(mapping = aes(x = as.factor(cluster), 
                         fill = as.factor(sample)), 
           position = "fill",
           stat = "identity")

(table(kmeans.res$idents,
       metaOrganID) %>%
    apply(MARGIN = 2,
          FUN = function(moID){
            moID * 100 / sum(moID)
          }) %>%
    as.data.frame() %>%
    t() %>%
    reshape::melt() %>%
    "colnames<-"(c("organ",
                   "cluster",
                   "percentage"))) %>%
  ggplot(mapping = aes(y = percentage)) +
  theme_bw() +
  geom_bar(mapping = aes(x = as.factor(cluster), 
                         fill = as.factor(organ)), 
           position = "fill",
           stat = "identity")

# umap
n_trees <- 100
set.seed(1985)
kmeans.res$umap.res <- uwot::tumap(
  X = kmeans.res$centers,
  n_neighbors = n_neighbors,
  n_components = 2,
  metric = "euclidean",
  n_epochs = NULL,
  learning_rate = 1,
  scale = FALSE,
  init = "spectral",
  init_sdev = 1e-4,
  set_op_mix_ratio = 1,
  local_connectivity = 1,
  bandwidth = 1,
  repulsion_strength = 1,
  negative_sample_rate = 5,
  nn_method = "annoy",
  n_trees = n_trees,
  search_k = 2 * n_neighbors * n_trees,
  n_threads = NULL,
  n_sgd_threads = 0,
  y = NULL,
  pca = NULL,
  pcg_rand = TRUE,
  fast_sgd = FALSE,
  ret_model = TRUE,
  ret_nn = FALSE,
  ret_extra = c(),
  tmpdir = tempdir(),
  verbose = FALSE)

kmeans.res$umap.res <- uwot::umap_transform(
  X = conc,
  model = kmeans.res$umap.res,
  init_weighted = TRUE,
  search_k = NULL,
  tmpdir = tempdir(),
  n_epochs = NULL,
  n_threads = NULL,
  n_sgd_threads = 0,
  grain_size = 1,
  verbose = FALSE
) %>%
  "colnames<-"(c("UMAP_1",
                 "UMAP_2"))

saveRDS(object = kmeans.res,
        file = "kmeans.res_Tons_Aden.rds")

# umap colored by SampleID
(((kmeans.res$umap.res %>%
     as.data.frame() %>%
     "colnames<-"(c("UMAP_1",
                    "UMAP_2"))) %>%
    cbind(metaSampleID = metaSampleID) %>%
    ggplot(data = .,
           mapping = aes(x = UMAP_1,
                         y = UMAP_2,
                         color = as.factor(metaSampleID))) +
    theme_bw() +
    theme(legend.position = "none") +
    geom_point(size = I(0.01)) +
    ggtitle("color: patient ID") #+
  #guides(color = guide_legend(override.aes = list(size = I(5),
  #                                               alpha = 1)))
) %>%
    egg::set_panel_size(width = unit(3, "in"),
                        height = unit(3, "in")) %>%
    gridExtra::grid.arrange()) %>%
  ggsave(filename = "210614_UMAP_metaSampleID_Tons_Aden.png",
         width = 4,
         height = 4,
         dpi = 300)

(((kmeans.res$umap.res %>%
     as.data.frame() %>%
     "colnames<-"(c("UMAP_1",
                    "UMAP_2"))) %>%
    cbind(metaSampleID = metaOrganID) %>%
    ggplot(data = .,
           mapping = aes(x = UMAP_1,
                         y = UMAP_2,
                         color = as.factor(metaSampleID))) +
    theme_bw() +
    theme(legend.position = "right") +
    labs(color = "organ") +
    geom_point(size = I(0.01)) +
    ggtitle("") +
    guides(color = guide_legend(override.aes = list(size = I(5),
                                                    alpha = 1)))
) %>%
    egg::set_panel_size(width = unit(3, "in"),
                        height = unit(3, "in")) %>%
    gridExtra::grid.arrange()) %>%
  ggsave(filename = "210614_UMAP_metaOrganID_Tons_Aden.png",
         width = 4,
         height = 4,
         dpi = 300)

(((kmeans.res$umap.res %>%
     as.data.frame() %>%
     "colnames<-"(c("UMAP_1",
                    "UMAP_2"))) %>%
    cbind(metaClusters = kmeans.res$idents) %>%
    ggplot(data = .,
           mapping = aes(x = UMAP_1,
                         y = UMAP_2,
                         color = as.factor(metaClusters))) +
    theme_bw() +
    theme(legend.position = "right") +
    geom_point(size = I(0.01)) +
    labs(color = "clusters") +
    guides(color = guide_legend(override.aes = list(size = I(5),
                                                    alpha = 1)))) %>%
    egg::set_panel_size(width = unit(3, "in"),
                        height = unit(3, "in")) %>%
    gridExtra::grid.arrange()) %>%
  ggsave(filename = "210614_UMAP_metaCluster_Tons_Aden.png",
         width = 5,
         height = 4,
         dpi = 300)

metadata.umap <- (kmeans.res$umap.res %>%
                    as.data.frame() %>%
                    "colnames<-"(c("UMAP_1",
                                   "UMAP_2"))) %>%
  cbind(status = setNames(object = relevant_meta$status,
                          nm = relevant_meta$patient_id)[metaSampleID_clean],
        cluster = kmeans.res$idents,
        organ = metaOrganID)

set.seed(1985)
((which(metadata.umap$status == "CONTROL") %>%
    sample(x = .,
           size = sum(metadata.umap$status != "CONTROL"),
           replace = FALSE) %>%
    c(which(metadata.umap$status != "CONTROL"),
      .) %>%
    metadata.umap[.,] %>%
    ggplot(data = .,
           mapping = aes(x = UMAP_1,
                         y = UMAP_2)) +
    theme_bw() +
    theme(legend.position = "right") +
    geom_hex(bins = 128,
             mapping = aes(fill = stat(log(count)))) +
    scale_fill_viridis_c(option = "magma") +
    guides(color = guide_legend(override.aes = list(size = I(5),
                                                    alpha = 1))) +
    facet_wrap(facets = ~ status + organ)) %>%
    egg::set_panel_size(width = unit(3, "in"),
                        height = unit(3, "in")) %>%
    gridExtra::grid.arrange()) %>%
  ggsave(filename = "210614_UMAP_COVID_Tons_Aden_hex.png",
         width = 8,
         height = 8,
         dpi = 300)

set.seed(1985)
((which(metadata.umap$status == "CONTROL") %>%
    sample(x = .,
           size = sum(metadata.umap$status != "CONTROL"),
           replace = FALSE) %>%
    c(which(metadata.umap$status != "CONTROL"),
      .) %>%
    metadata.umap[.,] %>%
    ggplot(data = .,
           mapping = aes(x = UMAP_1,
                         y = UMAP_2,
                         color = cluster)) +
    theme_bw() +
    theme(legend.position = "right") +
    geom_point(size = I(0.01)) +
    #scale_fill_viridis_c(option = "magma") +
    guides(color = guide_legend(override.aes = list(size = I(5),
                                                    alpha = 1))) +
    facet_wrap(facets = ~ status + organ)) %>%
    egg::set_panel_size(width = unit(3, "in"),
                        height = unit(3, "in")) %>%
    gridExtra::grid.arrange()) %>%
  ggsave(filename = "210614_UMAP_COVID_Tons_Aden_cluster.png",
         width = 8,
         height = 8,
         dpi = 300)

(((kmeans.res$umap.res %>%
     as.data.frame() %>%
     "colnames<-"(c("UMAP_1",
                    "UMAP_2"))) %>%
    cbind(status = setNames(object = relevant_meta$sex,
                            nm = relevant_meta$patient_id)[metaSampleID_clean]) %>%
    ggplot(data = .,
           mapping = aes(x = UMAP_1,
                         y = UMAP_2)) +
    theme_bw() +
    #theme(legend.position = "top") +
    geom_point(size = I(0.01)) +
    #labs(color = "status") +
    guides(color = guide_legend(override.aes = list(size = I(5),
                                                    alpha = 1))) +
    facet_wrap(facets = ~ status)) %>%
    egg::set_panel_size(width = unit(3, "in"),
                        height = unit(3, "in")) %>%
    gridExtra::grid.arrange()) %>%
  ggsave(filename = "210614_UMAP_gender_Tons_Aden.pdf",
         width = 8,
         height = 4,
         useDingbats = FALSE)

(((kmeans.res$umap.res %>%
     as.data.frame() %>%
     "colnames<-"(c("UMAP_1",
                    "UMAP_2"))) %>%
    cbind(status = setNames(object = relevant_meta$age,
                            nm = relevant_meta$patient_id)[metaSampleID_clean]) %>%
    ggplot(data = .,
           mapping = aes(x = UMAP_1,
                         y = UMAP_2,
                         color = status)) +
    theme_bw() +
    theme(legend.position = "top") +
    geom_point(size = I(0.01)) +
    labs(color = "Age") +
    #guides(color = guide_legend(override.aes = list(size = I(5),
    #                                               alpha = 1))) +
    scale_color_gradientn(colours = c("#000000",
                                      "#ffff00",
                                      "#ff0000"))
) %>%
    egg::set_panel_size(width = unit(3, "in"),
                        height = unit(3, "in")) %>%
    gridExtra::grid.arrange()) %>%
  ggsave(filename = "210614_UMAP_age_Tons_Aden.pdf",
         width = 4,
         height = 4.5,
         useDingbats = FALSE)

top.prot.exp <- conc %>%
  dplyr::as_tibble() %>%
  cbind(cell.pop = as.character(kmeans.res$idents)) %>%
  dplyr::group_by(cell.pop) %>%
  dplyr::summarize_all(.funs = mean)

top.prot.exp <- data.frame(top.prot.exp[,-1],
                           row.names = top.prot.exp$cell.pop) %>%
  as.matrix() #%>%
#t() %>% #SEE BELOW
#scale(center = TRUE,
#       scale = TRUE) #%>%
#t() #SEE BELOW

top.prot.HM <- pheatmap::pheatmap(mat = top.prot.exp,
                                  color = grDevices::colorRampPalette(c("#000000",
                                                                        "#ffff00",
                                                                        "#ff0000"))(100),
                                  kmeans_k = NA,
                                  breaks = NA,
                                  border_color = NA,
                                  scale = "none", ### CHANGE t() ABOVE TO CONTROL BETTER ###
                                  angle_col = 90,
                                  cluster_rows = TRUE,
                                  cellwidth = 10,
                                  cellheight = 10,
                                  treeheight_row = 10,
                                  treeheight_col = 10,
                                  silent = FALSE,
                                  main = "rows: clusters")

top.prot.HM[[4]] %>%
  gridExtra::grid.arrange() %>%
  ggsave(filename = "210614_HM_markers_Tons_Aden.png",
         width = 4,
         height = 3,
         dpi = 300)

# plot cluster marker
markers <- colnames(conc)

markers %>%
  setNames(nm = .) %>%
  lapply(FUN = function(marker) {
    ((kmeans.res$umap.res %>%
        as.data.frame() %>%
        "colnames<-"(c("UMAP_1",
                       "UMAP_2")) %>%
        cbind(marker = conc[,marker])) %>%
       ggplot(data = .,
              mapping = aes(x = UMAP_1,
                            y = UMAP_2,
                            color = marker)) +
       theme_bw() +
       theme(legend.position = "none") +
       geom_point(size = I(0.01)) +
       labs(title = marker) +
       scale_color_gradientn(colours = c("#000000",
                                         "#ffff00",
                                         "#ff0000"))
    ) %>%
      egg::set_panel_size(width = unit(2, "in"),
                          height = unit(2, "in"))
  }) %>%
  c(ncol = ceiling(sqrt(length(markers)))) %>%
  do.call(gridExtra::grid.arrange,
          .) %>%
  ggsave(filename = "210614_UMAP_markers_Tons_Aden.png",
         width = 15,
         height = 15,
         dpi = 300)

# map metacell clusters to cells
cells_dat <- c(names(metacells_Tonsil),
               names(metacells_Adenoid)) %>%
  lapply(FUN = function(fcs){
    if(grepl(pattern = "^Tonsil ",
             x = fcs,
             ignore.case = FALSE,
             perl = FALSE, 
             fixed = FALSE,
             useBytes = FALSE)) {
      res <- kmeans.res$idents[names(kmeans.res$idents) == fcs][metacells_Tonsil[[fcs]][["cluster"]]]  
    } else {
      res <- kmeans.res$idents[names(kmeans.res$idents) == fcs][metacells_Adenoid[[fcs]][["cluster"]]]  
    }
    
    total.n <- length(res)
    
    res <- res %>%
      paste0("cluster_",
             .) %>%
      table() %>%
      as.matrix() %>%
      t() %>%
      "rownames<-"(fcs)
    
    organ <- fcs %>%
      strsplit(split = " ",
               fixed = TRUE,
               perl = FALSE,
               useBytes = FALSE) %>%
      lapply(FUN = "[[",
             1) %>%
      unlist(use.names = FALSE)
    
    patient_id <- fcs %>%
      gsub(pattern = "^Adenoid |^Tonsil ",
           replacement = "",
           x = .,
           ignore.case = FALSE,
           perl = FALSE,
           fixed = FALSE,
           useBytes = FALSE) %>%
      gsub(pattern = ".fcs$",
           replacement = "",
           x = .,
           ignore.case = FALSE,
           perl = FALSE,
           fixed = FALSE,
           useBytes = FALSE) %>%
      gsub(pattern = "\\-.*|\\_.*",
           replacement = "",
           x = .,
           ignore.case = FALSE,
           perl = FALSE,
           fixed = FALSE,
           useBytes = FALSE) %>%
      gsub(pattern = "[^[:alnum:]]",
           replacement = "",
           x = .,
           fixed = FALSE)
    
    status <- relevant_meta$status[relevant_meta$patient_id == unique(patient_id)]
    
    age <- relevant_meta$age[relevant_meta$patient_id == unique(patient_id)]
    
    sex <- relevant_meta$sex[relevant_meta$patient_id == unique(patient_id)]
    
    dpi <- relevant_meta$dpi[relevant_meta$patient_id == unique(patient_id)]
    
    res <- data.table::data.table(patient_id = patient_id,
                                  organ = organ,
                                  status = status,
                                  age = age,
                                  sex = sex,
                                  dpi = dpi,
                                  total.n = total.n,
                                  res)
    
    return(res)
  }) %>%
  data.table::rbindlist(use.names = TRUE,
                        fill = TRUE) %>%
  dplyr::as_tibble()

cells_dat <- cells_dat %>%
  dplyr::mutate(across(.cols = starts_with(match = "cluster_"),
                       .fns = ~ tidyr::replace_na(data = .x,
                                                  replace = 0)),
                organ = factor(x = organ,
                               levels = c("Adenoid",
                                          "Tonsil")))

#saveRDS(object = cells_dat,
#      file = "cells_dat.rds")

cells_dat <- readRDS(file = "cells_dat.rds")

# statistical models and plots of results
mod_Adenoid <- levels(kmeans.res$idents) %>%
  paste0("cluster_",
         .) %>%
  setNames(nm = .) %>%
  lapply(FUN = function(cluster){
    paste0("rank(",
           cluster,
           "/total.n) ~ 
                   sex + 
                   age + 
                   status") %>%
      as.formula() %>%
      lm(formula = .,
         data = cells_dat[cells_dat$organ == "Adenoid",])
  })

source("~/Box/TsangLab/users/PedroMilanezAlmeida/code/custom_functions/plotModelRes.R")

plotModelRes(mod.list = mod_Adenoid) %>%
  ggsave(filename = "stats.Adenoid.png",
         width = 6,
         height = 3,
         dpi = 300)

mod_Tonsil <- levels(kmeans.res$idents) %>%
  paste0("cluster_",
         .) %>%
  setNames(nm = .) %>%
  lapply(FUN = function(cluster){
    paste0("rank(",
           cluster,
           "/total.n) ~ 
                   sex + 
                   age + 
                   status") %>%
      as.formula() %>%
      lm(formula = .,
         data = cells_dat[cells_dat$organ != "Adenoid",])
  })

plotModelRes(mod.list = mod_Tonsil) %>%
  ggsave(filename = "stats.Tonsil.png",
         width = 6,
         height = 3,
         dpi = 300)

mod_Interact <- levels(kmeans.res$idents) %>%
  paste0("cluster_",
         .) %>%
  setNames(nm = .) %>%
  lapply(FUN = function(cluster){
    paste0("rank(",
           cluster,
           "/total.n) ~ 
                   sex + 
                   age + 
                   status * 
                   organ") %>%
      as.formula() %>%
      nlme::lme(fixed = .,
                data = cells_dat,
                random = ~ 1 | patient_id)
  })

plotModelRes(mod.list = mod_Interact) %>%
  ggsave(filename = "stats.interaction.png",
         width = 9,
         height = 3,
         dpi = 300)

# neg binomial
negbin_Adenoid <- levels(kmeans.res$idents) %>%
  paste0("cluster_",
         .) %>%
  setNames(nm = .) %>%
  lapply(FUN = function(cluster){
    paste0(cluster,
           " ~
                   offset(log(total.n)) + 
                   sex + 
                   age + 
                   status") %>%
      as.formula() %>%
      MASS::glm.nb(formula = .,
                   data = cells_dat[cells_dat$organ == "Adenoid",],
                   control = glm.control(maxit = 500),
                   init.theta = 1
      )
  })

plotModelRes(mod.list = negbin_Adenoid)

negbin_Tonsil <- levels(kmeans.res$idents) %>%
  paste0("cluster_",
         .) %>%
  setNames(nm = .) %>%
  lapply(FUN = function(cluster){
    paste0(cluster,
           " ~
                   offset(log(total.n)) + 
                   sex + 
                   age + 
                   status") %>%
      as.formula() %>%
      MASS::glm.nb(formula = .,
                   data = cells_dat[cells_dat$organ != "Adenoid",],
                   control = glm.control(maxit = 500),
                   init.theta = 1
      )
  })

plotModelRes(mod.list = negbin_Tonsil)

negbin_Interact <- levels(kmeans.res$idents) %>%
  paste0("cluster_",
         .) %>%
  setNames(nm = .) %>%
  lapply(FUN = function(cluster){
    paste0(cluster,
           " ~
                   offset(log(total.n)) + 
                   sex + 
                   age + 
                   status * organ +
                   (1 | patient_id)") %>%
      as.formula() %>%
      lme4::glmer.nb(formula = .,
                     data = cells_dat,
                     nAGQ = 0,
                     initCtrl = list(theta = 1)
      )
  })

plotModelRes(mod.list = negbin_Interact)

# model time since PCR
mod_Adenoid_dpi <- levels(kmeans.res$idents) %>%
  paste0("cluster_",
         .) %>%
  setNames(nm = .) %>%
  lapply(FUN = function(cluster){
    paste0("rank(",
           cluster,
           "/total.n) ~ 
                   dpi") %>%
      as.formula() %>%
      lm(formula = .,
         data = cells_dat[cells_dat$organ == "Adenoid",])
  })

plotModelRes(mod.list = mod_Adenoid_dpi)  %>%
  ggsave(filename = "stats.dpi.Adenoid.png",
         width = 2.2,
         height = 3,
         dpi = 300)

mod_Tonsil_dpi <- levels(kmeans.res$idents) %>%
  paste0("cluster_",
         .) %>%
  setNames(nm = .) %>%
  lapply(FUN = function(cluster){
    paste0("rank(",
           cluster,
           "/total.n) ~ 
                   dpi") %>%
      as.formula() %>%
      lm(formula = .,
         data = cells_dat[cells_dat$organ != "Adenoid",])
  })

plotModelRes(mod.list = mod_Tonsil_dpi) %>%
  ggsave(filename = "stats.dpi.Tonsil.png",
         width = 2.2,
         height = 3,
         dpi = 300)

# update on 06/23/21: save frequencies as excel file
library(tidyverse)

cells_dat <- readRDS(file = "cells_dat.rds")

cells_dat_freq <- cells_dat %>%
  dplyr::mutate(dplyr::across(.cols = tidyselect::starts_with(match = "cluster_"),
                              .fns = ~ .x/total.n)) %>%
  dplyr::rename_with(~ paste0(.x,
                              "_freq"), 
                     starts_with("cluster_")) %>%
  ungroup() %>%
  dplyr::select(starts_with("cluster_"))

cbind(cells_dat,
      cells_dat_freq) %>%
  write_csv(file = "cells_dat_freq.csv")

# update of "model time since PCR": 09/20/21
# when using rank(dpi), lm ranks all rows prior to omitting NAs,
# resulting in a ranking that includes the non-COVID patients!
# re-running below but now with data that already excludes NA in dpi
# prior to ranking
library(tidyverse)

cells_dat <- readRDS(file = "cells_dat.rds")

mod_Adenoid_dpi <- colnames(x = cells_dat) %>%
  grep(pattern = "cluster",
       x = .,
       value = TRUE,
       fixed = TRUE) %>%
  setNames(nm = .) %>%
  lapply(FUN = function(cluster){
    paste0("rank(",
           cluster,
           "/total.n) ~ 
                   dpi") %>%
      as.formula() %>%
      lm(formula = .,
         data = cells_dat[cells_dat$organ == "Adenoid" &
                            !is.na(cells_dat$dpi),])
  })

plotModelRes(mod.list = mod_Adenoid_dpi)  %>%
  ggsave(filename = "stats.dpi.Adenoid.UPDATE.png",
         width = 2.2,
         height = 3,
         dpi = 300)

mod_Tonsil_dpi <- colnames(x = cells_dat) %>%
  grep(pattern = "cluster",
       x = .,
       value = TRUE,
       fixed = TRUE) %>%
  setNames(nm = .) %>%
  lapply(FUN = function(cluster){
    paste0("rank(",
           cluster,
           "/total.n) ~ 
                   dpi") %>%
      as.formula() %>%
      lm(formula = .,
         data = cells_dat[cells_dat$organ != "Adenoid" &
                            !is.na(cells_dat$dpi),])
  })

plotModelRes(mod.list = mod_Tonsil_dpi) %>%
  ggsave(filename = "stats.dpi.Tonsil.UPDATE.png",
         width = 2.2,
         height = 3,
         dpi = 300)

# updated on 121220: save plots as pdf for editting in illustrator and change pop names
library(ggrepel)
#library(scattermore)
library(tidyverse)

kmeans.res <- readRDS(file = "kmeans.res_Tons_Aden.rds")

dat <- kmeans.res$umap.res %>%
  "colnames<-"(c("UMAP_1",
                 "UMAP_2")) %>%
  dplyr::as_tibble() %>%
  dplyr::bind_cols(idents = kmeans.res$idents)
#names(pop.names[kmeans.res$idents]))

(ggplot(data = dat,
        mapping = aes(x = !!sym(colnames(dat)[1]),
                      y = !!sym(colnames(dat)[2]),
                      color = idents)) +
    theme_void(base_size = 22) +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5)) +
    geom_point() +
    geom_text_repel(data = dat %>%
                      dplyr::group_by(idents) %>%
                      dplyr::summarize(dplyr::across(dplyr::everything(),
                                                     median)) %>%
                      dplyr::ungroup() %>%
                      as.data.frame(),
                    mapping = aes(x = !!sym(colnames(dat)[1]),
                                  y = !!sym(colnames(dat)[2]),
                                  label = idents),
                    color = "black", 
                    min.segment.length = unit(0, 
                                              'lines'),
                    box.padding = 1,
                    size = I(6)) +
    guides(color = guide_legend(override.aes = list(size = I(5),
                                                    alpha = 1))) +
    labs(color = "cluster",
         title = "CD4 T cells",
         subtitle = "adenoids/tonsils")
) %>%
  egg::set_panel_size(width = unit(4, "in"),
                      height = unit(4, "in")) %>%
  gridExtra::grid.arrange()  %>%
  ggsave(filename = "211220_UMAP_metaCluster_Tons_Aden.pdf",
         width = 5.5,
         height = 5.5,
         useDingbats = FALSE)

#names(pop.names[kmeans.res$idents]))

(ggplot(data = dat,
        mapping = aes(x = !!sym(colnames(dat)[1]),
                      y = !!sym(colnames(dat)[2]),
                      color = idents == "12")) +
    theme_void(base_size = 22) +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5)) +
    geom_point() +
    guides(color = guide_legend(override.aes = list(size = I(5),
                                                    alpha = 1))) +
    labs(color = "cluster 12",
         title = "CD4 T cells",
         subtitle = "adenoids/tonsils")
) %>%
  egg::set_panel_size(width = unit(4, "in"),
                      height = unit(4, "in")) %>%
  gridExtra::grid.arrange()  %>%
  ggsave(filename = "211220_UMAP_metaCluster_Tons_Aden_CD4_12.pdf",
         width = 5.5,
         height = 5.5,
         useDingbats = FALSE)

(ggplot(data = dat,
        mapping = aes(x = !!sym(colnames(dat)[1]),
                      y = !!sym(colnames(dat)[2]),
                      color = idents == "11")) +
    theme_void(base_size = 22) +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5)) +
    geom_point() +
    guides(color = guide_legend(override.aes = list(size = I(5),
                                                    alpha = 1))) +
    labs(color = "cluster 11",
         title = "CD4 T cells",
         subtitle = "adenoids/tonsils")
) %>%
  egg::set_panel_size(width = unit(4, "in"),
                      height = unit(4, "in")) %>%
  gridExtra::grid.arrange()  %>%
  ggsave(filename = "211220_UMAP_metaCluster_Tons_Aden_CD4_11.pdf",
         width = 5.5,
         height = 5.5,
         useDingbats = FALSE)

# updated on 121223: save plots as pdf for editing in illustrator and change bar names
library(ggrepel)
#library(scattermore)
library(tidyverse)

kmeans.res <- readRDS(file = "kmeans.res_Tons_Aden.rds")

cells_dat <- readRDS(file = "cells_dat.rds")

source("~/Box/TsangLab/users/PedroMilanezAlmeida/code/custom_functions/plotModelRes.R")

# new facet label names for variables in models
var.labs <- c(sexMale = paste("(female)","    ","sex","    ","(male)"),
              age = paste("(younger)","    ","age","    ","(older)"),
              statusCOVID = paste("(neg)","    ","COVID status","    ","(pos)"))

# statistical models and plots of results
mod_Adenoid <- levels(kmeans.res$idents) %>%
  paste0("cluster_",
         .) %>%
  setNames(nm = .) %>%
  lapply(FUN = function(cluster){
    paste0("rank(",
           cluster,
           "/total.n) ~ 
                   sex + 
                   age + 
                   status") %>%
      as.formula() %>%
      lm(formula = .,
         data = cells_dat[cells_dat$organ == "Adenoid",])
  })

plotModelRes(mod.list = mod_Adenoid,
             labeller = ggplot2::labeller(variable = var.labs),
             pvalue.text.size = 5,
             base_size = 22) %>%
  egg::set_panel_size(width = unit(4, "in"),
                      height = unit(4, "in")) %>%
  gridExtra::grid.arrange()  %>%
  ggsave(filename = "stats.Adenoid_CD4.pdf",
         width = 14,
         height = 6,
         useDingbats = FALSE)

mod_Tonsil <- levels(kmeans.res$idents) %>%
  paste0("cluster_",
         .) %>%
  setNames(nm = .) %>%
  lapply(FUN = function(cluster){
    paste0("rank(",
           cluster,
           "/total.n) ~ 
                   sex + 
                   age + 
                   status") %>%
      as.formula() %>%
      lm(formula = .,
         data = cells_dat[cells_dat$organ != "Adenoid",])
  })

plotModelRes(mod.list = mod_Tonsil,
             labeller = ggplot2::labeller(variable = var.labs),
             pvalue.text.size = 5,
             base_size = 22) %>%
  egg::set_panel_size(width = unit(4, "in"),
                      height = unit(4, "in")) %>%
  gridExtra::grid.arrange()  %>%
  ggsave(filename = "stats.Tonsil_CD4.pdf",
         width = 14,
         height = 6,
         useDingbats = FALSE)

# updated on 121223: save plots as pdf for editing in illustrator required updating code to use meta-centroids rather than centroids to avoid recalculation
top.prot.exp <- kmeans.res$centers %>%
  dplyr::as_tibble() %>%
  cbind(cell.pop = as.character(kmeans.res$leiden)) %>%
  dplyr::group_by(cell.pop) %>%
  dplyr::summarize_all(.funs = mean)

top.prot.exp <- data.frame(top.prot.exp[,-1],
                           row.names = top.prot.exp$cell.pop) %>%
  as.matrix() #%>%
#t() %>% #SEE BELOW
#scale(center = TRUE,
#       scale = TRUE) #%>%
#t() #SEE BELOW

top.prot.HM <- pheatmap::pheatmap(mat = top.prot.exp,
                                  color = grDevices::colorRampPalette(c("#000000",
                                                                        "#ffff00",
                                                                        "#ff0000"))(100),
                                  kmeans_k = NA,
                                  breaks = NA,
                                  border_color = NA,
                                  scale = "none", ### CHANGE t() ABOVE TO CONTROL BETTER ###
                                  angle_col = 90,
                                  cluster_rows = TRUE,
                                  cellwidth = 10,
                                  cellheight = 10,
                                  treeheight_row = 10,
                                  treeheight_col = 10,
                                  silent = FALSE,
                                  main = "rows: clusters")

top.prot.HM[[4]] %>%
  gridExtra::grid.arrange() %>%
  ggsave(filename = "211223_HM_markers_Tons_Aden_CD4.pdf",
         width = 5,
         height = 3,
         useDingbats = FALSE)

# update on 211226: more pdfs
library(tidyverse)

kmeans.res <- readRDS(file = "kmeans.res_Tons_Aden.rds")

# plot cluster marker
markers <- colnames(kmeans.res$centers)

markers %>%
  setNames(nm = .) %>%
  lapply(FUN = function(marker) {
    ((kmeans.res$umap.res %>%
        as.data.frame() %>%
        "colnames<-"(c("UMAP_1",
                       "UMAP_2")) %>%
        cbind(marker = kmeans.res$centers[,marker][kmeans.res$cluster])) %>%
       ggplot(data = .,
              mapping = aes(x = UMAP_1,
                            y = UMAP_2,
                            color = marker)) +
       theme_bw() +
       theme(legend.position = "none") +
       geom_point(size = I(0.01)) +
       labs(title = marker) +
       scale_color_gradientn(colours = c("#000000",
                                         "#ffff00",
                                         "#ff0000"))
    ) %>%
      egg::set_panel_size(width = unit(2, "in"),
                          height = unit(2, "in"))
  }) %>%
  c(ncol = ceiling(sqrt(length(markers)))) %>%
  do.call(gridExtra::grid.arrange,
          .) %>%
  ggsave(filename = "211226_UMAP_markers_Tons_Aden_CD4.pdf",
         width = 15,
         height = 15,
         useDingbats = FALSE)

