rm(list = ls())

install_packages <- TRUE

if(install_packages){
  BiocManager::install("DOSE")
  BiocManager::install("clusterProfiler")
  BiocManager::install("enrichplot")
  BiocManager::install("ggupset")
  BiocManager::install("gmt")
  BiocManager::install("hypeR")
  BiocManager::install("gprofiler2")
  BiocManager::install("msigdbr")
  
}

rm(list = ls())

library(DOSE)
library(clusterProfiler)
library(tidyverse)
library(readr)
library(DESeq2)
# library(vsn)
library(pheatmap)
library(RColorBrewer)
library(fgsea)
library(hypeR)
library(gprofiler2)
library(magrittr)
library(data.table)
# library(biomaRt)
library(readxl)
library(enrichplot)
# library(org.Hs.eg.db)
library(ggupset)
library(gmt)
library(msigdbr)
library(ggplot2)

# Set PrimaryDirectory where this script is located
dirname(rstudioapi::getActiveDocumentContext()$path)  
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
PrimaryDirectory <- getwd()
PrimaryDirectory

# set working directory
workingDir <- "221010_WorkingDirectory"
dirPath <- file.path(PrimaryDirectory, workingDir)
setwd(dirPath)

# create other suite of directories

resPath <- file.path(dirPath, "results")
resDir <- file.path(dirPath, "results")
dir.create(resPath)
filePath <- file.path(resPath, "csv_files")
dir.create(filePath)

setwd(filePath)

gseaPath <- file.path(filePath, "enrichment_results")
dir.create(gseaPath)

gene_setsPath <- file.path(filePath, "gene_sets")
dir.create(gene_setsPath)

# grab DEG files
files <- list.files(filePath, pattern = ".tsv$", full = FALSE)

#### NEW STUFF ####

# import file with combined counts, contrasts
my_tibble_GSEA <- read_tsv("DA_results_BM_TUM_3contrasts_DESEQ2.tsv")

# list of contrasts
list_of_contrasts <- c(
  TPEX_DIFF = "tissue_cellBM_TPEX - tissue_cellTUM_TPEX",
  TEFF_DIFF = "tissue_cellBM_TEFF - tissue_cellTUM_TEFF",
  TISSUE_EFFECT = "(tissue_cellBM_TPEX + tissue_cellBM_TEFF) / 2 - (tissue_cellTUM_TPEX + tissue_cellTUM_TEFF + tissue_cellTUM_TTERM) / 3"
)

# create a loop to iterate over each contrast and spit out ranked list

## USEFUL PIECES
# counts_DA %>% select(contains("adj.P.Val"))

# get ranked gene list for TPEX differences (BM-TUM)
DA_TPEX_DIFF <- my_tibble_GSEA %>% select(symbol, timepoint, .abundant, contains("TPEX_DIFF"))

rnk_TPEX <- DA_TPEX_DIFF %>%
  dplyr::filter(adj.P.Val___TPEX_DIFF < 0.05) %>%
  dplyr::select(symbol, t___TPEX_DIFF) %>%
  na.omit() %>%
  distinct() %>%
  group_by(symbol) %>%
  summarise(stat = mean(t___TPEX_DIFF))
ranks <- rnk_TPEX %>% arrange(desc(stat))
# # ranks$stat <- rnk$stat * -1
ranks_TPEX <- deframe(ranks)

# get ranked gene list for TEFF differences (BM-TUM)
DA_TEFF_DIFF <- my_tibble_GSEA %>% select(symbol, timepoint, .abundant, contains("TEFF_DIFF"))

rnk_TEFF <- DA_TEFF_DIFF %>%
  dplyr::filter(adj.P.Val___TEFF_DIFF < 0.05) %>%
  dplyr::select(symbol, t___TEFF_DIFF) %>%
  na.omit() %>%
  distinct() %>%
  group_by(symbol) %>%
  summarise(stat = mean(t___TEFF_DIFF))
ranks <- rnk_TEFF %>% arrange(desc(stat))
# # ranks$stat <- rnk$stat * -1
ranks_TEFF <- deframe(ranks)

# get ranked gene list for tissue effect (all BM - all TUM)
DA_TISSUE_EFFECT <- my_tibble_GSEA %>% select(symbol, timepoint, .abundant, contains("TISSUE_EFFECT"))

rnk_TISSUE <- DA_TISSUE_EFFECT %>%
  dplyr::filter(adj.P.Val___TISSUE_EFFECT < 0.05) %>%
  dplyr::select(symbol, t___TISSUE_EFFECT) %>%
  na.omit() %>%
  distinct() %>%
  group_by(symbol) %>%
  summarise(stat = mean(t___TISSUE_EFFECT))
ranks <- rnk_TISSUE %>% arrange(desc(stat))
# # ranks$stat <- rnk$stat * -1
ranks_TISSUE_EFFECT <- deframe(ranks)

# at this point, have 3 objects that contain separate ranked gene lists per condition:
# ranks_TPEX
# ranks_TEFF
# ranks_TISSUE_EFFECT

# de_rank_all is a list of the ranks for each contrast
de_rank_all <- list(ranks_TEFF, ranks_TISSUE_EFFECT, ranks_TPEX)
names(de_rank_all) <- c("ranks_TEFF", "ranks_TISSUE_EFFECT", "ranks_TPEX")

#################################################################
# implement loop to make a list of collections


# TPEX BM v TUM D28

# BMvTUM_TPEX <- fread(files[2])
# BMvTUM_TEFF <- fread(files[1])
# rnk <- BMvTUM_TEFF %>%
#   dplyr::filter(padj < 0.05) %>%
#   dplyr::select(symbol, stat) %>%
#   na.omit() %>%
#   distinct() %>%
#   group_by(symbol) %>%
#   summarise(stat = mean(stat))
# ranks <- rnk %>% arrange(desc(stat))
# # # ranks$stat <- rnk$stat * -1
# ranks <- deframe(ranks)

## push to section with individual plots

# # see all gene sets for mouse
# all_gene_sets = msigdbr(species = "Mus musculus")
# head(all_gene_sets)
# 
# # view available collections
# msigdbr_collections() %>% print(n=23)
# 
# # pull specific collection of gene sets
# h_gene_sets = msigdbr(species = "mouse", category = "H")
# head(h_gene_sets)
# 
# c7_gene_sets = msigdbr(species = "mouse", category = "C7")
# head(c7_gene_sets)
# 
# immune_gene_sets <- c7_gene_sets %>%
#   dplyr::filter(gs_subcat == "IMMUNESIGDB")
# immune_gene_sets


###################################################################################################

### pull in collections from msigdbr
### Claudio looping function to create a collection

create_collection <- function(categories, ...) {
  collection <- lapply(categories, msigdbr::msigdbr, ...)
  names(collection) <- categories
  
  if ("C5" %in% names(collection)) {
    collection$C5 <- collection$C5 %>%
      dplyr::filter(gs_subcat == "GO:BP")
  }
  collection
}
# define categories of interest
categories <- c("C1", "C2", "C3", "C5", "C7", "C8", "H")

collections <- create_collection(categories, species = "Mus musculus")



### Claudio function to correctly format collections 
#' make_fgsea_collections_format
#'
#' fgsea wants pathways organized in list of lists. This is a format different than the one of msigdb, so we do some splitting.
#' We also want to do the enrichment for each subcategory.
#'
#' @param collections List of collections of pathways
#'
#' @return nested list of gs category and subcategories
make_fgsea_collections_format <- function(collections) {
  collections <- lapply(collections, function(x) {
    x[, "gs_cat_subcat"] <- paste0(x$gs_cat, "_", x$gs_subcat)
    split.data.frame(x, x$gs_cat_subcat)
  })
  collections <- flatten(collections)
  names(collections) <- sub("_$", "", names(collections))
  
  
  collections_split <- lapply(collections, function(db) split(x = db$gene_symbol, f = db$gs_name))
  names(collections_split) <- names(collections)
  collections_split
}


# Collection of gene lists in format liked by fgsea
collections_fgsea <- make_fgsea_collections_format(collections)


# de-ranked list of gene lists
de_rank_all


## generate enrichment results and bind into dataframe with collection ID

# generate enrichment results for all gene sets within collection
enr_results_ranks_TEFF <- lapply(collections_fgsea, FUN = fgseaMultilevel, stats = de_rank_all$ranks_TEFF, eps = 0)

#' Bind the enrichment results of a given contrast in a dataframe with .id collection
df_enr_TEFF <- bind_rows(enr_results_ranks_TEFF, .id = "collection")

enr_results_ranks_TPEX <- lapply(collections_fgsea, FUN = fgseaMultilevel, stats = de_rank_all$ranks_TPEX, eps = 0)
df_enr_TPEX <- bind_rows(enr_results_ranks_TPEX, .id = "collection")

enr_results_ranks_TISSUE <- lapply(collections_fgsea, FUN = fgseaMultilevel, stats = de_rank_all$ranks_TISSUE_EFFECT, eps = 0)
df_enr_TISSUE <- bind_rows(enr_results_ranks_TISSUE, .id = "collection")

# from here can access the significantly enriched pathways within the df objects

# take enrichment results and put into a dataframe, identifiable by contrast
enr_results_fl <- bind_rows(list("TPEX" = df_enr_TPEX, 
                                 "TEFF" = df_enr_TEFF, 
                                 "TISSUE" = df_enr_TISSUE), 
                            .id = "contrast")

#' Bind  results of collections
#' #'
#' #' Bind the enrichment results of a contrast in a dataframe with .id collection
#' #' @param enr_result Enriched result
#' bind_collections <- function(enr_result) {
#'   bind_rows(enr_result, .id = "collection")
#' }
#' 



### Claudio's code to make df for heatmap
#'
#' dataframe for heatmatp
#'
#' Create a dataframe with collection, pathway, contrast and l_idx (Luigi's index) (`-log10(pval) * sign(NES)`) to be used for the heatmap.
#' If `ntop` is set it will keep top negative and positive patways with min `padj`. It also calculate the l_idx (Luigi's index) (`-log10(pval) * sign(NES)`) .`C1` is filtered out.
#'
#' @param enr_results_fl A dataframe of flatten enriched results.
#' @param contrasts_keep List of contrast names to keep.
#' @param ntop Number of top pathway for each NES sign to keep. If NULL no filtering.
#' @param p_adj Padj to filter for. Default is 0.05
make_df_heatmap <- function(enr_results_fl,
                            contrasts_keep = c("contrast1", "contrast2", "contrast3"),
                            ntop = NULL,
                            p_adj = 0.05) {
  
  
  # var_plot <- sym(var_plot)
  
  dat <- enr_results_fl %>%
    filter(contrast %in% !!contrasts_keep) %>%
    mutate(contrast = factor(contrast, levels = contrasts_keep)) %>%
    mutate(sign_nes = sign(NES))
  
  # We want to keep the pathways if at least one among the contrast is less then
  # p_adj. The following code is used to  create a filtered dataframe and to pull the name of the pathwyas (pat)
  # that will  be used to filter the dataframe.
  
  # Find name of pathways (pat) with at least one significant and slice for ntop if needed
  dat_filt <- dat %>%
    filter(padj <= !!p_adj)
  
  if (!is.null(ntop)) {
    dat_filt <- dat_filt %>%
      group_by(collection, contrast, sign_nes) %>%
      slice_min(.data$padj, n = ntop)
  }
  
  pat <- dat_filt %>%
    pull(pathway)
  
  dat %>%
    filter(pathway %in% pat) %>%
    # mutate(pathway = droplevels(pathway)) %>%
    mutate(l_idx = -log10(pval) * sign(NES)) %>%
    mutate(sign_nes = if_else(sign_nes == 1, "Up", "Down")) %>%
    mutate(sign_nes = factor(sign_nes, levels = c("Up", "Down"))) %>%
    ungroup()
}


df_l_idx <- make_df_heatmap(
  enr_results_fl,
  contrasts_keep = c(
    "TPEX",
    "TEFF",
    "TISSUE"
  ),
  ntop = NULL,
  p_adj = 0.05
)

# write to file
# df_l_idx %>% write_csv("df_l_idx.csv")

# read from file


# create heatmap of all results

heatmap_allenriched <- df_l_idx %>%
 # filter(padj < 0.05) %>%
 filter(abs(NES) >= 2) %>%
 filter(collection == "C7_IMMUNESIGDB") %>%
  na.omit() %>%
  as_tibble() %>%
#  keep_variable(.abundance = NES, top = 500) %>%
  #as_tibble() %>%
  heatmap(
    .column = contrast,
    .row = pathway,
    .value = NES,
    #transform = log1p,
   scale = "none",
    # name = "TPEX BM v TUM",
    palette_value = c("red", "white", "blue"),
    #column_title = "BM TPEX vs TUM TPEX, D28 TOP 500 Genes",
    #  column_title_gp = gpar(fill = c("red", "blue", "green"), font = 1:2),
  #  row_km = 2,
    #   row_title_gp = gpar(col = c("blue", "red", "green", "orange"), font = 1:4),
   # column_km = 2, 
   # border = TRUE, 
    #width = unit(8, "cm"),
    # height = unit(8, "cm")
  )

heatmap_allenriched


df <- df[is.finite(rowSums(df)),]

#####################################################################################################

# break up given collection by gene sets
msigdbr_list = split(x = immune_gene_sets$gene_symbol, f = immune_gene_sets$gs_name)



# GSEA
fgseaRes <- fgsea(pathways = msigdbr_list, ranks, maxSize = 500, minSize = 20, eps = 0.0)

# write results to file
fgseaRes %>% write_tsv("fgseaResults_immunesets_TEFF.tsv")

# rank results by NES (high to low)
fgseaRes_ranked <- fgseaRes %>% arrange(desc(NES))

# write results to file
fgseaRes_ranked %>% write_tsv("fgseaResults_immunesets_NESranked_TEFF.tsv")


# plot enrichment of given pathway

sets_of_interest <- c("GOLDRATH_EFF_VS_MEMORY_CD8_TCELL_DN",
                      "GOLDRATH_EFF_VS_MEMORY_CD8_TCELL_UP",
                      
                      "GOLDRATH_NAIVE_VS_MEMORY_CD8_TCELL_UP",
                      "GOLDRATH_NAIVE_VS_MEMORY_CD8_TCELL_DN",
                      
                      "GOLDRATH_NAIVE_VS_EFF_CD8_TCELL_UP",
                      "GOLDRATH_NAIVE_VS_EFF_CD8_TCELL_DN",
                      
                      "GSE41867_DAY6_EFFECTOR_VS_DAY30_MEMORY_CD8_TCELL_LCMV_ARMSTRONG_UP",
                      "GSE41867_DAY6_EFFECTOR_VS_DAY30_MEMORY_CD8_TCELL_LCMV_ARMSTRONG_DN",
                      
                      "GSE9650_EXHAUSTED_VS_MEMORY_CD8_TCELL_DN",
                      "GSE9650_EXHAUSTED_VS_MEMORY_CD8_TCELL_UP",
                      
                      "GSE9650_EFFECTOR_VS_MEMORY_CD8_TCELL_DN",
                      "GSE9650_EFFECTOR_VS_MEMORY_CD8_TCELL_UP",
                      
                      "GSE9650_EFFECTOR_VS_EXHAUSTED_CD8_TCELL_UP",
                      "GSE9650_EFFECTOR_VS_EXHAUSTED_CD8_TCELL_DN",
                      
                      "GSE10239_MEMORY_VS_DAY4.5_EFF_CD8_TCELL_UP",
                      "GSE10239_MEMORY_VS_DAY4.5_EFF_CD8_TCELL_DN",
                      
                      "KAECH_DAY8_EFF_VS_MEMORY_CD8_TCELL_DN",
                      "KAECH_DAY8_EFF_VS_MEMORY_CD8_TCELL_UP",
                      
                      "KAECH_DAY15_EFF_VS_MEMORY_CD8_TCELL_DN",
                      "KAECH_DAY15_EFF_VS_MEMORY_CD8_TCELL_UP",
                      
                      "KAECH_NAIVE_VS_MEMORY_CD8_TCELL_UP", 
                      "KAECH_NAIVE_VS_MEMORY_CD8_TCELL_DN", 
            
                      "GSE41867_DAY8_EFFECTOR_VS_DAY30_EXHAUSTED_CD8_TCELL_LCMV_CLONE13_UP",
                      "GSE41867_DAY8_EFFECTOR_VS_DAY30_EXHAUSTED_CD8_TCELL_LCMV_CLONE13_DN"
                      )

for (i in 1:length(sets_of_interest)){
# multi_plot <-
  pdf(paste0("enrichplot_", sets_of_interest[i], ".pdf"))
   plotEnrichment(collections_fgsea$C7_IMMUNESIGDB[[sets_of_interest[i]]], 
                              de_rank_all$ranks_TISSUE_EFFECT) + 
                              labs(title = sets_of_interest[i])
   dev.off()
 #return(multi_plot)
# print(sets_of_interest[i])
# ggsave(paste0("enrichplot_", sets_of_interest[i], ".pdf"), plot = last_plot(), device = "pdf")
# return(multi_plot)
 print(i)

}


for(i in 1:length(sets_of_interest)){
  tiff(paste0("enrichplot_", sets_of_interest[i], ".tiff"), width = 14, height = 14, units = 'cm', res = 400)
  print(plotEnrichment(collections_fgsea$C7_IMMUNESIGDB[[sets_of_interest[i]]], 
                 de_rank_all$ranks_TISSUE_EFFECT) + 
    labs(title = sets_of_interest[i]) +
    theme(axis.title=element_text(size=12,face="bold"), plot.title=element_text(size = 12, face = "bold")))
  dev.off()

}


# for individual plots
# plotEnrichment(msigdbr_list[["GOLDRATH_NAIVE_VS_EFF_CD8_TCELL_UP"]], ranks) + 
#   labs(title = "Naive v EFFCD8, Goldrath")
# 
# plotEnrichment(msigdbr_list[["GOLDRATH_EFF_VS_MEMORY_CD8_TCELL_UP"]], ranks) + 
#   labs(title = "eff v mem CD8_UP, Goldrath")
# 
# plotEnrichment(msigdbr_list[["GSE9650_EXHAUSTED_VS_MEMORY_CD8_TCELL_UP"]], ranks) + 
#   labs(title = "exh vs mem CD8_UP")
# 
# plotEnrichment(msigdbr_list[["GSE9650_EFFECTOR_VS_EXHAUSTED_CD8_TCELL_UP"]], ranks) + 
#   labs(title = "eff vs exh CD8_UP")


############# Make the DF_L_IDX object more digestible
# select immune sets
tibble_immune_sets <- df_l_idx %>%
              as_tibble() %>%
              filter(padj < 0.05) %>%
              filter(abs(NES) >= 2) %>%
              filter(collection == "C7_IMMUNESIGDB")

# break up by contrast
immune_fgseaRes_ranked_TISSUE <- tibble_immune_sets %>%
                                  filter(contrast == "TISSUE") %>%
                                  arrange(desc(NES))

immune_fgseaRes_ranked_TPEX <- tibble_immune_sets %>%
                                filter(contrast == "TPEX") %>%
                                arrange(desc(NES))

immune_fgseaRes_ranked_TEFF <- tibble_immune_sets %>%
                                filter(contrast == "TEFF") %>%
                                arrange(desc(NES))

# pathways containing "CD8", should be filtered by significance 
CD8_pathways <- immune_fgseaRes_ranked_TISSUE %>% 
  filter(grepl("CD8", pathway)) %>%
  pull(pathway)

CD8_pathways_ofinterest <- c(
  "GSE30962_ACUTE_VS_CHRONIC_LCMV_PRIMARY_INF_CD8_TCELL_UP",
  "GSE30962_ACUTE_VS_CHRONIC_LCMV_SECONDARY_INF_CD8_TCELL_UP",
  "GSE39110_DAY3_VS_DAY6_POST_IMMUNIZATION_CD8_TCELL_UP",
  "GSE9650_EXHAUSTED_VS_MEMORY_CD8_TCELL_DN",
  "GOLDRATH_EFF_VS_MEMORY_CD8_TCELL_DN",
  "GSE41867_NAIVE_VS_DAY8_LCMV_EFFECTOR_CD8_TCELL_DN",
  "GSE15324_NAIVE_VS_ACTIVATED_CD8_TCELL_UP",
  "GSE10239_NAIVE_VS_KLRG1INT_EFF_CD8_TCELL_UP",
  "GSE9650_NAIVE_VS_MEMORY_CD8_TCELL_DN",
  "GSE10239_MEMORY_VS_DAY4.5_EFF_CD8_TCELL_UP",
  "GSE41867_NAIVE_VS_DAY15_LCMV_EFFECTOR_CD8_TCELL_DN",
  "GSE9650_EFFECTOR_VS_EXHAUSTED_CD8_TCELL_UP",
  "KAECH_NAIVE_VS_MEMORY_CD8_TCELL_DN", 
  "GSE10239_NAIVE_VS_KLRG1HIGH_EFF_CD8_TCELL_UP",
  "KAECH_DAY8_EFF_VS_MEMORY_CD8_TCELL_DN",
  "GSE23568_ID3_TRANSDUCED_VS_ID3_KO_CD8_TCELL_UP",
  "KAECH_NAIVE_VS_DAY15_EFF_CD8_TCELL_DN",
  "GSE41867_DAY6_EFFECTOR_VS_DAY30_MEMORY_CD8_TCELL_LCMV_ARMSTRONG_UP",
  "GSE14699_NAIVE_VS_ACT_CD8_TCELL_DN",
  "GOLDRATH_NAIVE_VS_EFF_CD8_TCELL_UP",
  "KAECH_DAY8_EFF_VS_DAY15_EFF_CD8_TCELL_DN",
  "GSE9650_EFFECTOR_VS_MEMORY_CD8_TCELL_DN",
  "GSE9650_NAIVE_VS_EXHAUSTED_CD8_TCELL_UP",
  "GSE41867_NAIVE_VS_DAY6_LCMV_EFFECTOR_CD8_TCELL_UP",
  "GSE9650_NAIVE_VS_MEMORY_CD8_TCELL_UP", 
  "GSE41867_DAY8_EFFECTOR_VS_DAY30_EXHAUSTED_CD8_TCELL_LCMV_CLONE13_UP",
  "GSE9650_NAIVE_VS_EXHAUSTED_CD8_TCELL_DN",
  "KAECH_DAY15_EFF_VS_MEMORY_CD8_TCELL_UP",
  "GSE19825_NAIVE_VS_DAY3_EFF_CD8_TCELL_DN",
  "GSE23321_CD8_STEM_CELL_MEMORY_VS_NAIVE_CD8_TCELL_UP", 
  "GSE14699_NAIVE_VS_ACT_CD8_TCELL_UP",
  "GSE10239_MEMORY_VS_KLRG1INT_EFF_CD8_TCELL_DN",
  "GSE41867_DAY6_VS_DAY8_LCMV_ARMSTRONG_EFFECTOR_CD8_TCELL_DN", 
  "GSE41867_DAY8_VS_DAY15_LCMV_CLONE13_EFFECTOR_CD8_TCELL_DN", 
  "GSE23321_CENTRAL_VS_EFFECTOR_MEMORY_CD8_TCELL_UP",
  "GSE41867_LCMV_ARMSTRONG_VS_CLONE13_DAY6_EFFECTOR_CD8_TCELL_UP",
  "GSE10239_KLRG1INT_VS_KLRG1HIGH_EFF_CD8_TCELL_UP",
  "GSE10239_MEMORY_VS_KLRG1HIGH_EFF_CD8_TCELL_DN",
  "GSE41978_KLRG1_HIGH_VS_LOW_EFFECTOR_CD8_TCELL_DN",
  "GSE15324_NAIVE_VS_ACTIVATED_CD8_TCELL_DN",
  "KAECH_NAIVE_VS_DAY8_EFF_CD8_TCELL_DN",
  "GSE41867_NAIVE_VS_EFFECTOR_CD8_TCELL_UP", 
  "GSE30962_PRIMARY_VS_SECONDARY_CHRONIC_LCMV_INF_CD8_TCELL_UP",
  "GSE9650_EXHAUSTED_VS_MEMORY_CD8_TCELL_UP", 
  "GSE10239_MEMORY_VS_DAY4.5_EFF_CD8_TCELL_DN",
  "GSE10239_NAIVE_VS_KLRG1INT_EFF_CD8_TCELL_DN", 
  "GSE9650_EFFECTOR_VS_MEMORY_CD8_TCELL_UP",
  "KAECH_DAY8_EFF_VS_MEMORY_CD8_TCELL_UP",
  "GSE23568_ID3_KO_VS_WT_CD8_TCELL_UP",
  "KAECH_DAY8_EFF_VS_DAY15_EFF_CD8_TCELL_UP",
  "GSE10239_NAIVE_VS_KLRG1HIGH_EFF_CD8_TCELL_DN",
  "GOLDRATH_NAIVE_VS_EFF_CD8_TCELL_DN", 
  "GSE23568_CTRL_VS_ID3_TRANSDUCED_CD8_TCELL_DN",
  "GSE39110_DAY3_VS_DAY6_POST_IMMUNIZATION_CD8_TCELL_DN", 
  "GSE30962_ACUTE_VS_CHRONIC_LCMV_PRIMARY_INF_CD8_TCELL_DN", 
  "GOLDRATH_EFF_VS_MEMORY_CD8_TCELL_UP", 
  "GSE30962_ACUTE_VS_CHRONIC_LCMV_SECONDARY_INF_CD8_TCELL_DN", 
  "GSE10239_NAIVE_VS_DAY4.5_EFF_CD8_TCELL_DN", 
  "GSE30962_PRIMARY_VS_SECONDARY_ACUTE_LCMV_INF_CD8_TCELL_UP"
)


# write to csv to look at 
immune_fgseaRes_ranked_TISSUE %>% write_tsv("immune_fgseaRes_ranked_TISSUE.tsv")

# look at top and bottom pathways
# topPathwaysUp <- immune_fgseaRes_ranked_TISSUE[NES > 0][head(order(pval), n = 50), pathway]

topPathwaysUp <- df_enr_TISSUE[NES > 0][head(order(padj), n=20), pathway]

topCD8pathways <- df_enr_TISSUE[pathway = CD8_pathways]

# topPathwaysUp <-
  #immune_fgseaRes_ranked_TISSUE %>% filter(NES > 0) %>% arrange(padj) %>% pull(pathway) %>% head(20)

# filtered_results <- immune_fgseaRes_ranked_TISSUE %>% select(pathway, pval, padj, NES)


## THIS ONE WORKS
plotGseaTable(collections_fgsea$C7_IMMUNESIGDB[sets_of_interest], 
              stats = de_rank_all$ranks_TISSUE_EFFECT, 
              df_enr_TISSUE, 
              gseaParam = 0.5,
              colwidths = c(5, 2, 1, 1, 1))

#CD8_pathways_ofinterest

# export
pdf("GSEA_CD8pathways_handfulof interest_TISSUEDIFF.pdf", height = 8, width = 16)
plotGseaTable(collections_fgsea$C7_IMMUNESIGDB[sets_of_interest], 
              stats = de_rank_all$ranks_TISSUE_EFFECT, 
              df_enr_TISSUE, 
              gseaParam = 0.5,
              colwidths = c(5, 2, 1, 1, 1))
dev.off()




topPathwaysDown <- immune_fgseaRes_ranked_TISSUE[NES < 0][head(order(pval), n = 50), pathway]


# examine top and bottom n pathways

topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=20), pathway]
topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=20), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(msigdbr_list[topPathways], ranks, fgseaRes, 
              gseaParam=0.5)


# export plot
pdf("GSEA_immune_top20_bottom20_TPEX_allgenes.pdf", height = 10, width = 16)
plotGseaTable(msigdbr_list[topPathways], ranks, fgseaRes, 
              gseaParam=0.5)
dev.off()



# make it tidy
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES)) %>%
  filter(padj < 0.01) %>%
  filter(abs(NES) >= 2) 
fgseaResTidy











# plot NES spectrum - need to pare down selected pathways to use this plot
ggplot(immune_fgseaRes_ranked_TISSUE, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = padj < 0.01)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()

# get database
# getMsigdb(
#   org = c("mm"),
#   id = c("SYM", "EZID"),
#   version = getMsigdbversion()
# )

# sample pathways

# wherryexh.path <- paste(filePath, "gene_sets", 
#                         "GSE41867_MEMORY_VS_EXHAUSTED_CD8_TCELL_DAY30_LCMV_UP.v2022.1.Hs.gmt", sep = "/")
# 
# geneset.wherryexh.DN <- gmtPathways(wherryexh.path)

# hallmark.path <- paste(filePath, "gene_sets", "h.all.v7.1.symbols.gmt", sep = "/")
# pathway.HALLMARK <- gmtPathways(hallmark.path)



# code from Claudio to loop over collections
# looping over human collections bc they have mouse orthologs annotated? or mouse ones?



# GSEA



# KEGG pathway enrichment





#################################### Comparing to an External Dataset ############################################################

dirExternalCounts <- "Philip_norm_counts"
dirExternalCountsPath <- paste(PrimaryDirectory, dirExternalCounts, sep = "/")
dir.create(dirExternalCountsPath)

CountFiles <- list.files(dirExternalCountsPath, pattern = "*.txt", full = FALSE)
CountFiles

# list of files with counts
data_list <- CountFiles %>% lapply(read_tsv)

map_dfr(left_join, by = "entrez_id")
data_list

names(data_list) <- c("N1", "N2", "N3", "E7_1", "E7_2", "E7_3", "M1", "M2", "M3", "L35_1", "L35_2", "L35_3")

names(merged_tib) <- c("entrez_id", "N1", "N2", "N3", "E7_1", "E7_2", "E7_3", "M1", "M2", "M3", "L35_1", "L35_2", "L35_3")

library(annotables)

grcm38

geneIDs <- grcm38 %>% dplyr::select(entrez, symbol)

merged_tib_with_gene_names <- left_join(merged_tib, geneIDs, by = c("entrez_id" = "entrez"))
merged_tib_with_gene_names

unique_T <- merged_tib_with_gene_names %>% distinct(symbol, entrez_id, .keep_all = TRUE)
unique_T
