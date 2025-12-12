####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Pediatric HLHS Compendium
####  2022-09-01 by Xiao LI (Texas Heart Institute, US)
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Initiate  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
Ver <- '2'
Step <- 'PART21_Milo_Enrichment'
Project <- '2022_hlhs_dturaga'
Code_dir <- paste0('~/Documents/Bioinformatics/project/', Project, '/code/human_v', Ver, '/')

source(Sys.readlink(paste0(Code_dir, 'src/bioinformatics.R')))
source(Sys.readlink(paste0(Code_dir, 'src/scRNAseq.R')))
source(Sys.readlink(paste0(Code_dir, 'src/geneset_human.R')))
source(paste0(Code_dir, 'human_v', Ver, '.helper_functions.R'))

library('miloR')
library('igraph')
InitiateProject('Rivendell', Ver, Step, 'human', Project, 'bree')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Load data  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
full.srt <- readRDS('integrated/PART19.annotated.srt.rds')
scVI.dr <- readRDS('integrated/PART19.all_reductions.srt_dimreducs.rds')
scVI.dr <- scVI.dr$scVI
full.srt@reductions$scVI <- full.srt@reductions$sub_umap
identical(as.vector(full.srt$orig.name), str_split(rownames(scVI.dr@cell.embeddings), pattern = ':', simplify = T)[,3])
full.srt@reductions$scVI@cell.embeddings <- scVI.dr@cell.embeddings
full.srt@reductions$scVI@key <- 'SCVI_'
colnames(full.srt@reductions$scVI@cell.embeddings) <- paste0('SCVI_', 1:50)
rownames(full.srt@reductions$scVI@cell.embeddings) <- Cells(full.srt)
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  MiloR Process  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
RunMilo <- function(seurat, k = 20, d = 50, alpha = 0.05, prop = 0.2){
        sce <- as.SingleCellExperiment(seurat)
        milo.meta <- seurat@meta.data
        milo.obj <- Milo(sce)
        reducedDim(milo.obj, "UMAP") <- reducedDim(sce, "SUB_UMAP")
        milo.obj <- buildGraph(milo.obj, k = k, d = d, reduced.dim = 'SCVI')
        milo.obj <- makeNhoods(milo.obj, k = k, d = d, refined = TRUE, prop = prop, reduced_dims = 'SCVI')
        milo.obj <- calcNhoodDistance(milo.obj, d = d, reduced.dim = 'SCVI')
        milo.obj <- countCells(milo.obj, samples = "group1", meta.data = milo.meta)

        milo.design <- as.data.frame(xtabs(~ group2 + group1, data = milo.meta))
        milo.design <- milo.design[milo.design$Freq > 0, ]
        milo.design <- distinct(milo.design)
        rownames(milo.design) <- milo.design$group1

        milo.res <- testNhoods(milo.obj, design = ~ group2, design.df = milo.design, reduced.dim = 'SCVI')
        milo.obj <- buildNhoodGraph(milo.obj)
        return(list(milo.obj, milo.res))
}
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Milo per cell type HLHS vs Controls  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
tmp.srt <- DropMetaLevels(full.srt[, full.srt$Non_ambiguous & full.srt$group2 != 'VAD'])
Table(tmp.srt$group1, tmp.srt$group2)

srt.list <- list(
        tmp.srt[, tmp.srt$Cell_type == 'Lym'],
        tmp.srt[, tmp.srt$Cell_type == 'Mye'],
        tmp.srt[, tmp.srt$Cell_type == 'EC'],
        tmp.srt[, tmp.srt$Cell_type == 'CM'],
        tmp.srt[, tmp.srt$Cell_type == 'FB'],
        tmp.srt[, tmp.srt$Cell_type %in% c('PC', 'SMC')]
)
name.list <- c('Lym', 'Mye', 'EC', 'CM', 'FB', 'Mural')
k.list <- c(50, 20, 10, 5, 10, 10)
prop.list <- c(1, 0.2, 0.1, 0.2, 0.2, 0.2)
milo_out <- list()
p.list <- list()
for(i in 1:6){
        message(i)
        milo_out[[i]] <- RunMilo(srt.list[[i]], k = k.list[i], prop = prop.list[i])
        p.list[[i]] <- plotNhoodGraphDA(milo_out[[i]][[1]], milo_res = milo_out[[i]][[2]]) +
                scale_fill_distiller(palette = 'RdBu')
        gc()
        PlotPDF(paste0('1.', i, '.umap.', name.list[i], '_milo_result'), 8, 6)
        print(p.list[[i]])
        dev.off()
        ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        saveRDS(milo_out[[i]], paste0('analysis/PART21.', name.list[[i]], '_milo_obj_result_list.rds'))
        ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
}
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Milo per cell type Post-VAD vs Pre-VAD  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
tmp.srt <- DropMetaLevels(full.srt[, full.srt$Non_ambiguous & full.srt$group2 != 'Control'])
Table(tmp.srt$group1, tmp.srt$group2)

srt.list <- list(
        tmp.srt[, tmp.srt$Cell_type == 'Lym'],
        tmp.srt[, tmp.srt$Cell_type == 'Mye'],
        tmp.srt[, tmp.srt$Cell_type == 'EC'],
        tmp.srt[, tmp.srt$Cell_type == 'CM'],
        tmp.srt[, tmp.srt$Cell_type == 'FB'],
        tmp.srt[, tmp.srt$Cell_type %in% c('PC', 'SMC')]
)
name.list <- c('Lym', 'Mye', 'EC', 'CM', 'FB', 'Mural')
k.list <- c(50, 20, 10, 5, 10, 10)
prop.list <- c(1, 0.2, 0.1, 0.2, 0.2, 0.2)
milo_out <- list()
p.list <- list()
for(i in 1:6){
        message(i)
        milo_out[[i]] <- RunMilo(srt.list[[i]], k = k.list[i], prop = prop.list[i])
        p.list[[i]] <- plotNhoodGraphDA(milo_out[[i]][[1]], milo_res = milo_out[[i]][[2]]) +
                scale_fill_distiller(palette = 'RdBu')
        gc()
        PlotPDF(paste0('2.', i, '.umap.', name.list[i], '_milo_result'), 8, 6)
        print(p.list[[i]])
        dev.off()
        ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        saveRDS(milo_out[[i]], paste0('analysis/PART21.pre_vs_postvad_', name.list[[i]], '_milo_obj_result_list.rds'))
        ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
}
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Milo global HLHS vs Control  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
RunMilo <- function(seurat, k = 20, d = 50, alpha = 0.05, prop = 0.2){
        sce <- as.SingleCellExperiment(seurat)
        milo.meta <- seurat@meta.data
        milo.obj <- Milo(sce)
        reducedDim(milo.obj, "UMAP") <- reducedDim(sce, "CLEAN_UMAP")
        milo.obj <- buildGraph(milo.obj, k = k, d = d, reduced.dim = 'SCVI')
        milo.obj <- makeNhoods(milo.obj, k = k, d = d, refined = TRUE, prop = prop, reduced_dims = 'SCVI')
        milo.obj <- calcNhoodDistance(milo.obj, d = d, reduced.dim = 'SCVI')
        milo.obj <- countCells(milo.obj, samples = "group1", meta.data = milo.meta)

        milo.design <- as.data.frame(xtabs(~ group2 + group1, data = milo.meta))
        milo.design <- milo.design[milo.design$Freq > 0, ]
        milo.design <- distinct(milo.design)
        rownames(milo.design) <- milo.design$group1

        milo.res <- testNhoods(milo.obj, design = ~ group2, design.df = milo.design, reduced.dim = 'SCVI')
        milo.obj <- buildNhoodGraph(milo.obj)
        return(list(milo.obj, milo.res))
}

tmp.srt <- full.srt[, full.srt$Non_ambiguous & full.srt$group2 != 'VAD']

tmp.srt <- tmp.srt[, seq(1, ncol(tmp.srt), 5)]

tmp.srt$group2 <- droplevels(tmp.srt$group2)
tmp.srt$group1 <- droplevels(tmp.srt$group1)
Table(tmp.srt$group1, tmp.srt$group2)

milo_out <- RunMilo(tmp.srt, k = 20, d = 50, prop = 0.5)
p <- plotNhoodGraphDA.my(milo_out[[1]], milo_res = milo_out[[2]]) +
        scale_fill_distiller(palette = 'RdBu') +
        theme(aspect.ratio = 1, axis.line = element_line())
p
PlotPDF('3.1.umap.global_hlhs_vs_control', 8, 8)
p
dev.off()

saveRDS(milo_out, paste0('analysis/PART21.global_hlhs_vs_control_milo_obj_result_list.rds'))
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Milo global Post-VAD vs Pre-VAD  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
tmp.srt <- full.srt[, full.srt$Non_ambiguous & full.srt$group2 != 'Control']

tmp.srt <- tmp.srt[, seq(1, ncol(tmp.srt), 5)]

tmp.srt$group2 <- droplevels(tmp.srt$group2)
tmp.srt$group1 <- droplevels(tmp.srt$group1)
Table(tmp.srt$group1, tmp.srt$group2)

tmp.srt$group2 <- factor(tmp.srt$group2, levels = c('VAD', 'HLHS'))

milo_out <- RunMilo(tmp.srt, k = 20, d = 50, prop = 0.5)
p <- plotNhoodGraphDA.my(milo_out[[1]], milo_res = milo_out[[2]])  +
        scale_fill_distiller(palette = 'RdBu') +
        theme(aspect.ratio = 1, axis.line = element_line())
p
PlotPDF('3.2.umap.global_postvad_vs_prevad', 8, 8)
p
dev.off()

saveRDS(milo_out, paste0('analysis/PART21.global_prevad_vs_postvad_milo_obj_result_list.rds'))
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
