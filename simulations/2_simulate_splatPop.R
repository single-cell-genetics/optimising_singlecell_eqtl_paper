#' Use splatPop to simulated scRNA-seq data using iPSC Smart-Seq2 and FPP 10X 
#' data as references. 
#' splatPopParams estimated in 1_estimate_splatPopParams.R

#install.packages("/mnt/mcscratch/cazodi/Software/splatter_1.15.1.tar.gz",
#                 repos = NULL, type="source")

library(splatter)
library(VariantAnnotation)
library(SingleCellExperiment)
library(scater)

save_out_data <- "/mnt/mcscratch/cazodi/Projects/EUUI_2019_sceQTL-Workflow/data/sims_best-practices/"
gff <- read.table("/mnt/mcscratch/cazodi/Datasets/references/human/Homo_sapiens.GRCh38.98.chromosome.2.genes.gff3", sep="\t", header=TRUE)


############################
###### iPSC SmartSeq2 ###### 
############################
vcf <- readVcf("/mnt/mcscratch/cazodi/Projects/EUUI_2019_sceQTL-Workflow/output/best-prac_2021-02-15_ss2/01_genotypes/geno_87_01.vcf", "hg38")
params <- readRDS(save_out_data, "params_smartseq_use.rds")

params <- setParams(params, 
                    eqtl.n = 0.35, 
                    similarity.scale = 1.5,
                    batchCells = rep(1, 24),
                    batch.size = 5,
                    batch.facLoc = args$batch_facLoc, 
                    batch.facScale = args$batch_facScale,
                    nCells.sample = TRUE, 
                    nCells.shape = 1.387517805, 
                    nCells.rate = 0.006467842,
                    eqtl.dist = 1000000 ,
                    eqtl.maf.min = 0.1,
                    eqtl.maf.max = 0.5,
                    eqtl.ES.shape = 2.053687,
                    eqtl.ES.rate = 9.358627) 

# Simulate dataset with complex batch structure
sim <- splatPopSimulate(vcf = vcf, params = params, gff=gff)
sim <- logNormCounts(sim)
sim <- runPCA(sim)
# Generate plots
varExpMat <- getVarianceExplained(sim, variables=c("Batch", "Sample"))
plot_varExpMat(varExpMat, c("Batch", "Sample"))

sim_sub <- subset(sim, , Batch %in% c("Batch1", "Batch2"))
plotPCA(sim_sub, shape_by = "Batch", colour_by = "Sample", theme_size=14, point_size=3) + 
  scale_color_jcolors(palette = "pal8") 


#############
## FPP 10x ##
#############

vcf <- readVcf("/mnt/mcscratch/cazodi/Projects/EUUI_2019_sceQTL-Workflow/output/best-prac_2021-02-15_ss2/01_genotypes/geno_87_01.vcf", "hg38")
params <- readRDS(save_out_data, "params_smartseq_use.rds")

params <- setParams(params, 
                    eqtl.n = 0.35, 
                    similarity.scale = 1.5,
                    batchCells = rep(1, 24),
                    batch.size = 5,
                    batch.facLoc = args$batch_facLoc, 
                    batch.facScale = args$batch_facScale,
                    nCells.sample = TRUE, 
                    nCells.shape = 1.387517805, 
                    nCells.rate = 0.006467842,
                    eqtl.dist = 1000000 ,
                    eqtl.maf.min = 0.1,
                    eqtl.maf.max = 0.5,
                    eqtl.ES.shape = 2.053687,
                    eqtl.ES.rate = 9.358627) 

# Simulate dataset with complex batch structure
sim <- splatPopSimulate(vcf = vcf, params = params, gff=gff)
sim <- logNormCounts(sim)
sim <- runPCA(sim)
# Generate plots
varExpMat <- getVarianceExplained(sim, variables=c("Batch", "Sample"))
plot_varExpMat(varExpMat, c("Batch", "Sample"))

sim_sub <- subset(sim, , Batch %in% c("Batch1", "Batch2"))
plotPCA(sim_sub, shape_by = "Batch", colour_by = "Sample", theme_size=14, point_size=3) + 
  scale_color_jcolors(palette = "pal8") 
