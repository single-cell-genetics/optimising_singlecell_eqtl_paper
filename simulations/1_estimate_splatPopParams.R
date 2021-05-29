#' Estimate splatPopParams for iPSC Smart-Seq2 and FPP 10X data for the best
#' practices manuscript

library(splatter)
library(data.table)
library(tidyverse)
library(fitdistrplus)
library(SingleCellExperiment)
library(scater)

source("code/splatPop_helper_functions.R")

save_out_data <- "/mnt/mcscratch/cazodi/Projects/EUUI_2019_sceQTL-Workflow/data/sims_best-practices/"
data_path <- "/mnt/mcscratch/cazodi/Projects/EUUI_2019_sceQTL-Workflow/data/"
nCells <- 1255


############################
###### iPSC SmartSeq2 ###### 
############################
# Cuomo et al 2020
# Download from Zenodo: https://zenodo.org/record/3625024#.Xil-0y2cZ0s 

meta <- read.csv("/mnt/mcfiles/Datasets/single-cell-data/cuomo_iPSC_endoderm/cell_metadata_cols.tsv", sep="\t", header=TRUE)
meta$X <- NULL
counts <- as.data.frame(fread("/mnt/mcfiles/Datasets/single-cell-data/cuomo_iPSC_endoderm/counts.tsv", sep="\t", header=TRUE))

row.names(counts) <- counts$GeneID
counts$GeneID <- NULL
sce_ss <- SingleCellExperiment(assays = list(counts = counts),
                               colData = meta)
sce_ss <- subset(sce_ss, , (day == "day0"))

# Get the batch structure 
meta_summary <- meta %>% group_by(donor, experiment) %>% dplyr::count() 
length(unique(meta_summary$experiment)) # 28 batches
length(unique(meta_summary$donor))
mean(table(meta_summary$experiment)) # 5.7 donors per batch
sum(table(meta_summary$donor)-1) / length(unique(meta_summary$donor)) # 28% replicated

# Get the distribution of nCells per donor-experiment
sce_ss$donor_run <- paste(sce_ss$donor, sce_ss$experiment, sep="_")
n <- unlist(table(sce_ss$donor_run), use.names = FALSE)
outliers <- boxplot(n, plot=FALSE)$out
n <- n[-which(n %in% outliers)]
fit.gamma <- fitdist(as.numeric(n), distr = "gamma", method = "mle")
summary(fit.gamma) # shape: 1.90748601, rate: 0.02684891

# Generate plots
sce_ss <- logNormCounts(sce_ss)
sce_ss <- runPCA(sce_ss)
varExpMat <- getVarianceExplained(sce_ss, variables=c("donor", "experiment"))
plot_varExpMat(varExpMat, c("donor", "experiment"))
ggsave(filename = paste0(save_out_data, "iPSC-ss2_varExp.pdf"),
       width = 6, height = 4)

sce_sub <- subset(sce_ss, , experiment %in% c("expt_23", "expt_37"))
plotPCA(sce_sub, shape_by = "experiment", colour_by = "donor",
                theme_size=14, point_size=3) + 
  scale_color_jcolors(palette = "pal8") 
ggsave(filename = paste0(save_out_data, "iPSC-ss2_pca.pdf"), 
       width = 5, height = 4)


# Get mean aggregated matrix for estimating population-scale parameters 
colData(sce_ss)$donor_exp <- paste(sce_ss$donor, sce_ss$experiment, sep="_")
agg <- counts(aggregateAcrossCells(sce_ss, ids= sce_ss$donor_exp, statistics="mean"))
                                   
# Get SCE of individual with most cells (joxm=383) for estimating SC parameters.
# Sampling nCells I want to simulate to get accurate library size estimates
sce.joxm <- subset(sce_ss, , (donor == "joxm"))
sce.joxm <- sce.joxm[sample(1:nrow(sce.joxm), nCells), ]
sce.joxm.counts <- as.matrix(counts(sce.joxm))

# Estimate splatPopParams
params <- splatPopEstimate(counts=sce.joxm.counts, means = agg)
params <- setParams(params, pop.quant.norm=FALSE)
saveRDS(params, file = paste0(save_out_data, "params_smartseq_use.rds"))







#############
## FPP 10x ##
#############

# https://zenodo.org/record/4072909#.X7PUJy-l3dc

# Get subset of data to use for estimating params
# sce <- readRDS("data/sims_best-practices/2020-12-07_d11-sce.Rds")
# gene_list <- read.table("/mnt/mcfiles/Datasets/single-cell-data/cuomo_iPSC_dopaminergic-nerons/eqtl_summary_stats_renamed/D11_unique_genes.txt")
# 
# library("org.Hs.eg.db") # remember to install it if you don't have it already
# gene.symbols <- mapIds(org.Hs.eg.db, keys = gene_list$V1, keytype = "ENSEMBL", column="SYMBOL")
# gene.symbols <- gene.symbols[!is.na(gene.symbols)]
# genes.keep <- intersect(gene.symbols, rownames(sce))
# sce <- sce[genes.keep, ]
# sce <- subset(sce, , celltype == "FPP")
# 
# samples_keep <- table(sce$donor_pool)[table(sce$donor_pool) > 5]
# samples_keep <- names(samples_keep)
# sce <- subset(sce, , donor_pool %in% samples_keep)
# sce$donor_pool <- paste(sce$donor_id, sce$pool_id, sep=":")
# 
# saveRDS(sce, "data/sims_best-practices/2020-12-14_d11-FPP-sce-gene-subset.Rds")
sce <- readRDS("data/sims_best-practices/2020-12-14_d11-FPP-sce-gene-subset.Rds")


# Generate plots
sce <- logNormCounts(sce)
sce <- runPCA(sce)
varExpMat <- getVarianceExplained(sce, variables=c("donor_id", "pool_id"))
plot_varExpMat(varExpMat, c("donor_id", "pool_id"))
ggsave(filename = paste0(save_out_data, "FPP-10x_varExp.pdf"),
       width = 6, height = 4)

sce_sub <- subset(sce_ss, , experiment %in% c("expt_23", "expt_37"))
plotPCA(sce_sub, shape_by = "experiment", colour_by = "donor",
        theme_size=14, point_size=3) + 
  scale_color_jcolors(palette = "pal8") 
ggsave(filename = paste0(save_out_data, "FPP-10x_pca.pdf"), 
       width = 5, height = 4)

# Get the batch structure 
meta_summary <- data.frame(colData(sce)) %>% 
  group_by(donor_id, pool_id) %>% dplyr::count() 
length(unique(meta_summary$pool_id)) # 12 batches
mean(table(meta_summary$pool_id)) # 16.5 donors per batch
sum(table(meta_summary$donor_id)-1) / length(unique(meta_summary$donor_id)) # 13.8% replicated

# Get the distribution of nCells per donor
n <- meta_summary$n
outliers <- boxplot(n, plot=FALSE)$out
n <- n[-which(n %in% outliers)]
hist(n, breaks = 20)
fit.gamma <- fitdist(n, distr = "gamma", method = "mle", lower = c(0, 0))
summary(fit.gamma) # shape: 1.232613761, rate: 0.002342392


# Get mean aggregated matrix for estimating population-scale parameters 
colData(sce)$donor_exp <- paste(sce$donor_id, sce$pool_id, sep="_")
agg <- counts(aggregateAcrossCells(sce, ids= sce$donor_exp, statistics="mean"))

# Get SCE of individual with most cells (HPSI0414i-mita_1 = 14640) for estimating SC parameters.
# Sampling nCells I want to simulate to get accurate library size estimates
sce.mita1 <- subset(sce, , (donor_id == "HPSI0414i-mita_1"))
sce.mita1 <- sce.mita1[sample(1:nrow(sce.mita1), nCells), ]
sce.mita1.counts <- as.matrix(counts(sce.mita1))

# Estimate splatPopParams
params <- splatPopEstimate(counts=sce.mita1.counts, means = agg)
params <- setParams(params, pop.quant.norm=FALSE)
saveRDS(params, file = paste0(save_out_data, "params_10x_use.rds"))


