rm(list = ls())
# MacOS
# install xquartz using homebrew (brew cask install xquartz)
# install cairo (brew install cairo)

install_packages <- TRUE

install.packages("BiocManager")
install.packages("devtools")

if(!require(devtools)){
  BiocManager::install("devtools")
}

if(install_packages){
  Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = TRUE)
  BiocManager::install(c( 
    "tidyverse",
    "tibble",
    "dplyr",
    "tidyr",
    "readr",
    "stringr",
    "ggplot2",
    "purrr",
    "tidyHeatmap",
    "tidybulk",
    "ggrepel",
    "plotly",
    "GGally",
    "DESeq2",
    "limma",
    "edgeR",
    "data.table",
    "cluster",
    "factoextra",
    "tidymodels",
    "remotes",
    "pheatmap",
    "ComplexHeatmap"
  ))
}

devtools::install_github("jbengler/tidyheatmaps")
remotes::install_github("stemangiola/nanny")
devtools::install_github("stephenturner/annotables")


rm(list = ls())

# package install

# load packages
library(tidyverse)
library(tibble)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(ggplot2)
library(purrr)
library(tidyHeatmap)
library(tidybulk)
library(ggrepel)
library(plotly)
library(GGally)
library(DESeq2)
library(limma)
library(edgeR)
library(data.table)
library(cluster)
library(factoextra)
library(tidymodels)
library(nanny)
library(pheatmap)
library(tidyheatmap)
library(ComplexHeatmap)
library(annotables)


# what others?

# set Primary Directory
dirname(rstudioapi::getActiveDocumentContext()$path)  
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
PrimaryDirectory <- getwd()
PrimaryDirectory

# create folders//directories
dirRawCounts <- "raw_counts"
dirRawCountsPath <- paste(PrimaryDirectory, dirRawCounts, sep = "/")
dir.create(dirRawCountsPath)

workingDir <- "221010_WorkingDirectory"
workingDirPath <- paste(PrimaryDirectory, workingDir, sep = "/")
dir.create(workingDirPath)

# list count files
CountFiles <- list.files(dirRawCountsPath, pattern = "*.csv", full = FALSE)
CountFiles

setwd(dirRawCounts)

# read in dataset files
# list of files with counts
data_list <- CountFiles %>% lapply(read_csv)

# %>% map_dfr(left_join, by = "gene_ID")
data_list



 
# create tibble with merged count files
 merged_tibble <- data_list %>% purrr::reduce(., inner_join)

 # above function broken down
 # path.to.files <- file.path(PrimaryDirectory, "RNASeq_rawcounts", files)
 # list.of.tibbles <- map(path.to.files, read_csv)
 # myTibble <- list.of.tibbles %>% purrr::reduce(., inner_join)
 
 # map(file.path(PrimaryDirectory, "RNASeq_rawcounts", files), read_csv) %>% 
 #   purrr::reduce(., inner_join)
 
 
# transform the dimensions of the table, keeping gene_ID column as is
merged_tibble_longer <- pivot_longer(merged_tibble, !gene_ID, names_to = "condition", values_to = "count")

# save long form of table with gene IDs and counts for later access
# merged_tibble_longer %>% write_csv("merged_counts.csv")


merged_tibble_longer <- read_csv("merged_counts.csv")
################# Import More Data Into Our Tibble #########################

# try with annotables package 
# install.packages("devtools")

# mouse table with all info
grcm38

# mouse table with only enstxp and ensgene
# grcm38_tx2gene

geneIDs <- grcm38 %>% dplyr::select(ensgene, symbol)

# append gene names to tibble
with_gene_names <- left_join(merged_tibble_longer, geneIDs, by = c("gene_ID" = "ensgene"))
with_gene_names

# check for NAs in symbols
nanas <- with_gene_names %>% filter(is.na(symbol))
nanas

# create metadata file
exp_md <- fread(file = paste(PrimaryDirectory, "Exp32_metadata.csv", sep = "/"), header = TRUE)
exp_md

# add metadata to the tibble
big_tibble <- left_join(with_gene_names, exp_md, by = c("condition" = "sample_id"))
big_tibble

# remove duplicates and filter unique things by condition and by symbol? gene ID?
unique_T <-
  big_tibble %>% distinct(symbol, condition, .keep_all = TRUE)
unique_T

# check for empty values
nanas2 <- unique_T %>% filter(is.na(symbol))
nanas2

# MAKE THIS MORE SPECIFIC
unique_T_emptycheck <- unique_T %>%
  dplyr::mutate(symbol = replace_na(symbol, "96126"))

unique_T <- unique_T_emptycheck
unique_T
# can also create a unique identifier to add to the tibble, then filter by the combination of symbol and condition


unique_T <-
  unique_T %>% rename("condition" = "sample_ID")

# check uniqueness?

# tidy the data 
# create tidy_counts which renames sample ID, takes the parameters of interest out of the tibble
# feeds parameters into tidybulk

tidy_counts <-
  unique_T %>% mutate(sampleName = sample_ID) %>%
  select(symbol, sampleName, count, tissue, timepoint, mouse, cell_type) %>%
  tidybulk(.sample = sampleName, .transcript = symbol, .abundance = count)
tidy_counts

# change working directory
# setwd(workingDirPath)
# getwd()

# # write my_tibble to external file
# tidy_counts %>% write_csv("tidy_counts.csv")

# my_tibble has abundance identified based on factor of interest and scaled
# TMM is the default scaling method can use EdgeR or Limma as well

my_tibble <- tidy_counts %>%
  identify_abundant(factor_of_interest = tissue) %>%
  scale_abundance(method = "TMM")
my_tibble

# __________________________________________________

# read in scaled values

tidy_counts2 <- fread(file = paste(workingDirPath, "tidy_counts.csv", sep = "/"), header = TRUE)
  
tidy_counts2 <- tidy_counts2 %>%
  as_tibble() %>%
  tidybulk(.sample = sampleName, 
           .transcript = symbol, 
           .abundance = count)
tidy_counts2

my_tibble2 <- tidy_counts2 %>%
  identify_abundant(factor_of_interest = tissue) %>%
  scale_abundance(method = "TMM")
my_tibble2


# make it a tibble 
my_tibble <- my_tibble2


## FILTERING 
# filter selected samples only
my_tibble <- my_tibble %>% filter(grepl("late", sampleName))
my_tibble <- my_tibble %>% filter(grepl("TEFF", sampleName))
my_tibble <- my_tibble %>% filter(!grepl("TTERM", sampleName))

my_tibble

# create a combination term to add to the table which accounts for conditions included
# helps with building contrasts later, e.g. what is the combined effect of tissue and cell type
my_tibble <- my_tibble %>% mutate(tissue_cell = paste0(tissue, "_", cell_type))


# visualize scaling
my_tibble %>%
  filter(.abundant) %>%
  pivot_longer(cols = c("count", "count_scaled"), names_to = "source", values_to = "abundance") %>%
  ggplot(aes(x = abundance + 1, color = sampleName)) +
  geom_density() +
  facet_wrap(~source) +
  scale_x_log10() +
  theme_bw()

# validation # FIX ME!!
# check_if_counts_is_na(.data = my_tibble, .abundance = count_scaled)

# dimensionality reduction
# PCA or MDS plots reduce dimensions of the data to identify greatest sources of variation
# PCA: unsupervised, don't need to specify groups
# if your exp worked well, greatest source of variation should be the treatments / groups you're interested in
# good for QC and checking outliers
# note that this joins the result to the counts object

counts_scal_PCA <-
  my_tibble %>%
  filter(.abundant) %>%
  tidybulk::reduce_dimensions(method = "PCA",
                              top = 500,
                              log_transform = TRUE)
counts_scal_PCA

# pivot_sample takes a tbl as input and returns one with only sample-related columns

counts_scal_PCA %>% pivot_sample()

# plot reduced dimensions
# PCA plot

my_pca_plot <- counts_scal_PCA %>%
  pivot_sample() %>%
  ggplot(aes(x = PC1, y = PC2, colour = cell_type, shape = timepoint)) + 
  scale_color_manual(values = c("darkturquoise", "red")) +
  geom_point() +
  geom_text_repel(aes(label = sampleName), show.legend = FALSE) + 
  theme_bw() +
  labs(title = "PCA: Transcriptomic Profiling of BM TPEX \n D7 v D28, Top 500 Genes") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

my_pca_plot


# heirarchical clustering with heatmaps
# keep_variable extracts the most variable genes which can then be plotted
# check clustering to make sure samples cluster by intended group (e.g. tissue or condition)
heatmap_clustercheck <- my_tibble %>%
  filter(.abundant) %>%
  tidybulk::keep_variable(.abundance = count_scaled, top = 500) %>%
  as_tibble() %>%
  heatmap(
    .column = sampleName,
    .row = symbol,
    .value = count_scaled,
    transform = log1p,
    scale = "row",
   # name = "TPEX BM v TUM",
    palette_value = c("red", "white", "blue"),
    #column_title = "BM TPEX vs TUM TPEX, D28 TOP 500 Genes",
  #  column_title_gp = gpar(fill = c("red", "blue", "green"), font = 1:2),
    row_km = 2,
 #   row_title_gp = gpar(col = c("blue", "red", "green", "orange"), font = 1:4),
    column_km = 2, 
    border = TRUE, 
    #width = unit(8, "cm"),
   # height = unit(8, "cm")
 ) %>%
  # ) %>%
  # split_columns(2) %>%
  # split_rows(4) %>%
   add_tile(timepoint, palette = c("darkturquoise", "red"))# %>%
 # add_tile(timepoint) %>%
  # add_tile(cell_type, palette = c("green", "blue")) 
#add_tile(cell_type, palette = c("#F8766D", "#00BA38", "#619CFF")) %>%
 # packLegend(direction = "vertical")

heatmap_clustercheck



############## DE TESTING ###################

# test_differential_abundance currently uses edgeR to analyze for differential expression
# can specify different method e.g. DESeq2 or limma voom
# give our tidybulk counts object and a formula to specify the column that contains our groups for comparison
# 0 + dex if from same cell line
# 0 + dex + cell since each sample within the treated and untreated groups is from a different cell line
# add additional factor to accommodate this
# also provide the contrasts we want
  # tidybulk says: All testing methods use raw counts, irrespective of if scale_abundance 
  # or adjust_abundance have been calculated. Therefore, it is essential to add covariates 
  # such as batch effects (if applicable) in the formula.

## edgeR and DeSeq2 implementation for normalization and quantification
# (should get 91 to 99% identity between the two)
# how to check identity btwn edgeR and DeSeq2?

## READ MORE ABOUT DESIGN FORMULA AND CONTRASTS

# formula can include cell type, tissue or timepoint depending on analysis
# eventually need to figure out how to implement an interaction term

# DE counts with default (EdgeR)
counts_de <- counts_scal_PCA %>%
  test_differential_abundance(
    .formula = ~ 0 + tissue, #+ cell_type,
    .contrasts = c("tissueBM-tissueTUM"),
    omit_contrast_in_colnames = TRUE
    )

# DE counts with DESeq2
counts_de <- counts_scal_PCA %>%
  test_differential_abundance(
    .formula = ~ 0 + tissue, #+ cell_type,
    .contrasts = list(c("tissue", "BM", "TUM")),
    omit_contrast_in_colnames = FALSE,
    method = "DESeq2"
  )


# look into more detailed contrasts
# .contrasts = c(("timepointD28-timepointD7")-("cell_typeTPEX - cell_typeTEFF"))#, 

# see DE results joined to our tibble
counts_de

# select for DEGs only
counts_de %>% pivot_transcript()

# # write our DEGs into a file that can be loaded into excel
# counts_DE %>%
#   pivot_transcript() %>%
#   write_tsv("DA_results_BM_TUM_3contrasts_DESEQ2.tsv")

############## Differential Abundance Testing, No Subsetting, For GSEA ########################

# DA testing, not scaled nor subsetted, DESeq2
# counts_DA <- my_tibble %>%
#   test_differential_abundance(.formula = ~0 + tissue + cell_type + tissue:cell_type,
#                               contrasts = list(c("tissue", "BM", "TUM"), 
#                                                c("cell_type", "TPEX", "TEFF"), 
#                                                # c("tissue:cell_type", "BMTPEX", "TUMTPEX"),
#                                                # c("tissue:cell_type", "BMTEFF", "TUMTEFF")
#                                                ),
#                               omit_contrast_in_colnames = FALSE,
#                               method = "DESeq2"
#   )
# 
# 
# counts_DA_1 <- my_tibble %>%
#   test_differential_abundance(.formula = ~0 + tissue,
#                               contrasts = list(c("tissue", "BM", "TUM")),
#                               omit_contrast_in_colnames = TRUE,
#                               method = "DESeq2"
#   )

list_of_contrasts <- c(
  TPEX_DIFF = "tissue_cellBM_TPEX - tissue_cellTUM_TPEX",
  TEFF_DIFF = "tissue_cellBM_TEFF - tissue_cellTUM_TEFF",
  TISSUE_EFFECT = "(tissue_cellBM_TPEX + tissue_cellBM_TEFF) / 2 - (tissue_cellTUM_TPEX + tissue_cellTUM_TEFF + tissue_cellTUM_TTERM) / 3"
)

counts_DA <- my_tibble %>%
  test_differential_abundance(.formula = ~0 + tissue_cell,
                              # contrasts = list(c((tissue_cellBM_TPEX + tissue_cellBM_TEFF) / 2 - 
                              #                    (tissue_cellTUM_TPEX + tissue_cellTUM_TEFF + tissue_cellTUM_TTERM) / 3),
                              #                 c(tissue_cellBM_TPEX - tissue_cellTUM_TPEX),
                              #                 c(tissue_cellBM_TEFF - tissue_cellTUM_TEFF)
                              #                 ),
                              # contrasts = list(c("tissue_cellBM_TPEX - tissue_cellTUM_TPEX"),
                              #                  c("tissue_cellBM_TEFF - tissue_cellTUM_TEFF"),
                              #                  c("(tissue_cellBM_TPEX + tissue_cellBM_TEFF) / 2 - (tissue_cellTUM_TPEX + tissue_cellTUM_TEFF) / 2")
                              #                  ),
                              contrasts = list_of_contrasts,
                              omit_contrast_in_colnames = FALSE,
                              method = "limma_voom"
  )


counts_DA_names_adjusted <- counts_DA %>% 
                            rename_with(~ gsub("tissue_cellBM_TPEX - tissue_cellTUM_TPEX", 
                                               "TPEX_DIFF", .x, fixed = TRUE)) %>%
                            rename_with(~ gsub("tissue_cellBM_TEFF - tissue_cellTUM_TEFF", 
                                               "TEFF_DIFF", .x, fixed = TRUE)) %>%
                            rename_with(~ gsub("(tissue_cellBM_TPEX + tissue_cellBM_TEFF) / 2 - (tissue_cellTUM_TPEX + tissue_cellTUM_TEFF + tissue_cellTUM_TTERM) / 3", 
                                               "TISSUE_EFFECT", .x, fixed = TRUE))

# counts_DA_TPEXdifference <- counts_DA_names_adjusted %>% select(contains(""))
# write new counts into a file with compiled results
# do we filter for abundant before GSEA???
counts_DA_names_adjusted %>%
  pivot_transcript() %>%
  write_tsv("DA_results_BM_TUM_3contrasts_DESEQ2.tsv")

# rename again
counts_DA <- counts_DA_names_adjusted

# DA_ALL_GENE_RANK <- counts_DA %>%
#                     filter(padj %>% is.na %>% `!`) %>%
#                     tidybulk::test_gene_rank(
#                       .entrez = symbol,
#                       .arrange_desc = log2FoldChange,
#                       species = "Mus musculus",
#                       gene_sets = c("IMMUNESIGDB")
#                     )

# 
# # write DA tested counts into separate file
# counts_DA_meep %>%
#   pivot_transcript() %>%
#   write_tsv("DA_results_BM_TUM_lateALL_TWOCONTRASTS_DESEQ2.tsv")

##### LOOK AT DEGS, PLOTS ######

# count DEGs, EDGR method
DEGcount <- counts_de %>%
  filter(FDR < 0.05) %>%
  summarise(num_de = n_distinct(symbol))
DEGcount

# count DEGs, DESeq2
DEGcount <- counts_de %>%
  filter(padj < 0.05) %>%
  summarise(num_de = n_distinct(symbol))
DEGcount

# grab a couple of top DEGs
topgenes <- 
  counts_de %>%
  pivot_transcript() %>%
  arrange(PValue) %>%
  head(200)
topgenes

topgenes_symbols <- topgenes %>% pull(symbol)
topgenes_symbols

# to plot specific gene names below
topgenes_of_interest <- c("Lef1", "Rgs16", "Sell", "Mt1", "Lag3", "Pdcd1", "Irf8", "Gzmm", "Ldha", "Xcl1", "CD200", "Bcl2")
# # topgenes_of_interest <- c("Ly6a2", "Gzma", "Rasgrp2", "Klrg1", "Klrb1c",
#                           "S1pr5", "Xcl1","CD200", "Cd7", "Tox", "Zeb2",
#                           "S1pr1", "Myc", "Rgs16", "Pdcd1", "Lag3",
#                           "Gzmm", "Ccr7", "Tigit", "Lef1", "S1pr4")

# Myb paper related gene list
# topgenes_of_interest = c("Lef1", "Sell", "Gzma","Gzmm", "Myb", "Ccr7", "Satb1", "Bach2", "Slamf6", 
#                          "Id3","Tcf7", "Il7r", "Xcl1", "Mki67", "Zeb2", "Cx3cr1", "CD160", "Cd7", 
#                          "Lag3", "Pdcd1", "Bcl2", "S1pr5", "Klrg1", "Klrb1c", "Rgs16")

topgenes_of_interest = c("Rgs16", "Lag3", "Pdcd1", "Cd200", "Lef1", "Gzma", "Mki67", 
                         "Bcl2", "Gzmm", "Lef1", "Xcl1", "Sell", "Nod1", "Tox", "Cxcr4","Cd70", "Il2ra", "Il17ra")

topgenes_of_interest = c("Gzma", "Klrg1", "Klrb1c", "S1pr5", "Xcl1", "Cd200", "Cd7", "Tox", 
                         "Zeb2", "Rasa3", "S1pr1", "Pdcd1", "Rgs16", "Gzmm", "Lef1", "Ctla4", "Havcr2", "Fasl")


# list from BM v TUM late timepoint, all subsets (presumably tissue driven)
topgenes_of_interest = c("Lag3", "Cd7", "Gzmm", "Sell", "Cd200", "Xcl1", "Tox", "S1pr5", "Il2ra", "Klrb1c")

# list from BM D28 v D7
topgenes_of_interest = c("Cd69", "Fasl", "Fosb", "Jun", "Hdac1", "H2ax", "Ifng")

# list from BM D28 v D7 TPEX or TEFF
topgenes_of_interest = c("Gzma", "Gzmb", "CD69", "Bcl2", "Mki67", "Pdcd1", "Mybl2", "Il7r", "Klrg1", "Nfkbia", "Cxcr3")

## PLOTS AFTER DEG TESTING

# # volcano plots, minimal
# counts_de %>%
#   filter(.abundant) %>%
#   ggplot(aes(x = logFC, y = PValue, colour = FDR < 0.05)) + 
#   geom_point() + 
#   scale_y_continuous(trans = "log10_reverse") + 
#   theme_bw()

# # simple heatmap with significant genes of interest
# 
# minimap <- counts_de %>%
# 
#   # subset data
#   filter(.abundant) %>%
#   mutate(significant = FDR < 0.05 & abs(logFC) >= 2) %>%
#   mutate(symbol = ifelse(symbol %in% topgenes_of_interest, as.character(symbol), "")) %>%
#   # arrange(PValue) %>%
#   # head(20) %>%
#   
#   
#     tidy_heatmap(., rows = symbol, columns = sampleName, values = count_scaled
#   
#   
#     )

# more informative volcano plot
my_volcano_plot <- counts_de %>%
  pivot_transcript() %>%
  
  # subset data
  filter(.abundant) %>%
  mutate(significant = FDR < 0.05 & abs(logFC) >= 2) %>%
  mutate(expression = case_when(logFC >= log(2) & FDR <= 0.05 ~ "Up-regulated",
                       logFC <= -log(2) & FDR <= 0.05 ~ "Down-regulated",
                       TRUE ~ "Unchanged")
  ) %>%
 # mutate(Expression = upregulated if (FDR < 0.05 & logFC >= 2)) %>%
  #mutate(symbol = ifelse(symbol %in% topgenes_symbols, as.character(symbol), "")) %>%
   mutate(symbol = ifelse(symbol %in% topgenes_of_interest, as.character(symbol), "")) %>%
  
  
  # plot
  ggplot(aes(x = logFC, y = PValue, label = symbol)) +
  geom_point(aes(color = expression, size = significant, alpha=significant)) +
  geom_text_repel(max.overlaps = Inf) +
  labs(title = "BM D28 v D7 TPEX Subsets, Select Genes of Interest") +
  
  # custom scales
  scale_y_continuous(trans = "log10_reverse") +
  scale_color_manual(values=c("red", "black", "darkturquoise")) +
  scale_size_discrete(range = c(0, 2)) +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
my_volcano_plot

## STRIPCHARTS

# to look at exp level of genes of interest for individual samples
# help show if exp is consistent amongst replicates in groups

strip_chart <-
  my_tibble %>%
  
  # extract counts for top DEGs
 # filter(symbol %in% topgenes_symbols) %>%
   filter(symbol %in% topgenes_of_interest) %>%
  
  # make stripchart
  ggplot(aes(x = tissue, y = count_scaled + 1, fill = tissue, label = "")) +
  geom_boxplot() +
  geom_jitter() +
  facet_wrap(~symbol) +
  scale_y_log10() +
  theme_bw()

strip_chart


# CREATE PLOTS AND SAVE TO WD
setwd(workingDirPath)
getwd()

ggsave("stripchart_top200select_BM_TUM_lateTEFF.svg", plot = strip_chart, width = 6, height = 6)
ggsave("volcano_top200_ofinterest_BM_earlylate_TPEX.pdf", plot = my_volcano_plot, width = 8, height = 6, dpi = 400)
ggsave("PCA_BM_earlylate_top500.svg", plot = my_pca_plot, width = 6, height = 4)
save_pdf(.heatmap = heatmap_clustercheck, "clusteringcheck_BM_TUM_late_allsubsetstop500.svg", width = 6, height = 7)




# nrows <- 200; ncols <- 6
# counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
# rowRanges <- GRanges(rep(c("chr1", "chr2"), c(50, 150)),
#                      IRanges(floor(runif(200, 1e5, 1e6)), width=100),
#                      strand=sample(c("+", "-"), 200, TRUE),
#                      feature_id=sprintf("ID%03d", 1:200))
# colData <- DataFrame(Treatment=rep(c("ChIP", "Input"), 3),
#                      row.names=LETTERS[1:6])
# rse <- SummarizedExperiment(assays=SimpleList(counts=counts),
#                             rowRanges=rowRanges, colData=colData)


################### K MEANS CLUSTERING ########################

## k means clustering of top 500 DEGs
seed <- 123456
set.seed(seed)


# different means of determining cluster amount

# SSE
# wss <- (nrow(scaledata)-1)*sum(apply(scaledata,2,var))
# for (i in 2:20) wss[i] <- sum(kmeans(scaledata,
#                                      centers=i)$withinss)
# plot(1:20, wss, type="b", xlab="Number of Clusters",
#      ylab="Within groups sum of squares")

# # Avg Silhouette Width
# library(cluster)
# sil <- rep(0, 20)
# #repeat k-means for 1:20 and extract silhouette:
# for(i in 2:20){
#   k1to20 <- kmeans(scaledata, centers = i, nstart = 25, iter.max = 20)
#   ss <- silhouette(k1to20$cluster, dist(scaledata))
#   sil[i] <- mean(ss[, 3])
# }
# 
# # Plot the  average silhouette width
# plot(1:20, sil, type = "b", pch = 19, xlab = "Number of clusters k", ylab="Average silhouette width")
# abline(v = which.max(sil), lty = 2)



# pull cluster number from the cluster vector generated by kmeans

# x <- rbind(matrix(rnorm(100, sd = 0.3), ncol = 2),
#            matrix(rnorm(100, mean = 1, sd = 0.3), ncol = 2))
# colnames(x) <- c("x", "y")
# (cl <- kmeans(x, 2))
# 
# #Access the cluster vector
# cl$cluster
# # Map cluster vector to original data with cbind or similar
# out <- cbind(x, clusterNum = cl$cluster)


# clustering, add to tibble
tibble_k5 <- my_tibble %>% 
  
  # subset data - want top 500 abundant genes by scaled count
  filter(.abundant) %>% 
  tidybulk::keep_variable(.abundance = count_scaled, top = 500) %>%
  as_tibble() %>%
  
  # input into clustering
  tidybulk::cluster_elements(
    .element = sampleName,
    .feature = symbol,
    .abundance = count_scaled,
    method = "kmeans",
    log_transform = TRUE,
    centers = 5
  ) 

k5 <- tibble_k5 %>%
  
  # build heatmap
    heatmap(
      .column = sampleName,
      .row = symbol,
      .value = count_scaled,
      transform = log1p,
      scale = "row",
      palette_value = c("red", "white", "darkblue"),
      row_km = 5,# splits rows
      column_km = 3 # splits columns
    ) %>%

    add_tile(tissue) %>%
    add_tile(timepoint) %>%
    add_tile(cell_type) %>%
  add_tile('cluster kmeans') %>%

  as_ComplexHeatmap() %>%
  ComplexHeatmap::rowAnnotation(foo = anno_block(gp = gpar(fill = 1:5),
                                                                                labels = c("group1", "group2", "group3", "4", "5"), 
                                                                                labels_gp = gpar(col = "white", fontsize = 10)))





# add tile cell type (TPEX TEFF ETC)

# try with pheatmap / tidyheatmap
my_tibble %>% 
  
  # subset data - want top 500 abundant genes by scaled count
  filter(.abundant) %>% 
  tidybulk::keep_variable(.abundance = count_scaled, top = 100) %>%
  as_tibble() %>%
  
  tidybulk::cluster_elements(
    .element = sampleName,
    .feature = symbol,
    .abundance = count_scaled,
    method = "kmeans",
    log_transform = TRUE,
    centers = 5
  ) %>%
  tidy_heatmap(
    columns = sampleName,
    rows = symbol,
    values = count_scaled,
    scale = "row",
    colors = c("red", "white", "blue"),
    annotation_col = c(tissue, timepoint, cell_type, 'cluster kmeans'),
  #  annotation_row = c('cluster kmeans'),
    clustering_method = "complete",
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    clustering_distance_cols = "euclidean",
    clustering_distance_rows = "euclidean",
    gaps_row = NULL,
    gaps_col = NULL,
    #  cellwidth = 7,
    # cellheight = 7,
    #   kmeans_k = 5
  )
  
  tidy_heatmap(., rows = symbol, columns = sampleName, values = count_scaled, scale = "row",
               annotation_col = c(tissue, sampleName, cell_type), annotation_row = c(symbol))

tidyheatmap <- my_tibble %>%
  filter(.abundant) %>%
  tidybulk::keep_variable(.abundance = count_scaled, top = 100) %>%
  as_tibble() %>%
  tidy_heatmap(
    columns = sampleName,
    rows = symbol,
    values = count_scaled,
    scale = "row",
    colors = c("red", "white", "blue"),
    annotation_col = c(tissue, timepoint, cell_type),
    clustering_method = "complete",
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    clustering_distance_cols = "euclidean",
    clustering_distance_rows = "euclidean",
    gaps_row = NULL,
    gaps_col = NULL,
    #  cellwidth = 7,
    # cellheight = 7,
    #   kmeans_k = 5
  )
