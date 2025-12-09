####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Pediatric HLHS Compendium
####  2022-09-01 by Xiao LI (Texas Heart Institute, US)
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Initiate  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
Ver <- '2'
Step <- 'PART02_Merge'
Project <- '2022_hlhs_dturaga'
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Global Functions  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
Preprocess <- function(srt_obj, assay = 'RNA', ...) {
        srt.out <- srt_obj |>
                NormalizeData(verbose = F) |>
                CellCycleScoring(s.features = s.genes, g2m.features = g2m.genes, set.ident = F) |>
                PercentageFeatureSet(pattern = '^MT-', col.name = paste0('pct_mito_', assay), assay = assay)
        srt.out@meta.data[, paste0('pct_mito_', assay)][is.nan(srt.out@meta.data[, paste0('pct_mito_', assay)])] <- 0
        return(srt.out)
}
RunUMAP.my <- function(srt_obj, var.toal = 0.75, haromize.by, assay, ...){
        srt.out <- RunPCA(srt_obj, seed.use = 505, assay = assay)
        dimN <- FindDimNumber(srt_obj = srt.out, var.toal = var.toal, reduction = 'pca')
        srt.out <- RunUMAP(srt.out, dims = 1:dimN, seed.use = 505) |>
                RunHarmony(group.by.vars = haromize.by, assay.use = assay)
        dimN <- FindDimNumber(srt_obj = srt.out, var.toal = var.toal, reduction = 'harmony')
        srt.out <- RunUMAP(srt.out, dims = 1:dimN, seed.use = 505, reduction = 'harmony',
                           reduction.name = 'hmn_umap', reduction.key = 'hmnumap_', ...)
        return(srt.out)
}
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Load sample metadata  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
sample_meta.df <- read.csv(paste0(Docu_dir, 'pediatric_sample_meta.csv'))
studies <- U(sample_meta.df$Study)
names(studies) <- U(sample_meta.df$Study_id)
studies_cellbender <- studies[studies %in% sample_meta.df$Study[sample_meta.df$platform %in% c('10X', 'Drop-seq')]]
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Merge CellRanger   ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
rds_name <- paste0(names(studies), '.', studies, '.raw.srt.rds')
data.list <- list()
ncells <- 0
for(i in 1:L(rds_name)) {
        data.list[[i]] <- readRDS(paste0('individual/', rds_name[i]))
        message(data.list[[i]]$study[1], ' ', ncol(data.list[[i]]))
        ncells <- ncells + ncol(data.list[[i]])
}
merged.raw.srt <- merge(data.list[[1]], data.list[2:L(data.list)])
c(ncells, ncol(merged.raw.srt)) ## 489948 489948
rm(data.list)
gc()
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Merge CellBender    ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
rds_name <- paste0(names(studies), '.', studies, '.cbn.srt.rds')
data.list <- list()
ncells <- 0
for(i in 1:L(rds_name)) {
        data.list[[i]] <- readRDS(paste0('individual/', rds_name[i]))
        message(data.list[[i]]$study[1], ' ', ncol(data.list[[i]]))
        ncells <- ncells + ncol(data.list[[i]])
}
merged.cbn.srt <- merge(data.list[[1]], data.list[2:L(data.list)])
c(ncells, ncol(merged.cbn.srt)) ## 522465 522465
rm(data.list)
gc()

## Intersect data sets
O(Cells(merged.raw.srt), Cells(merged.cbn.srt)) ## common 461146
consensus_cell <- intersect(Cells(merged.raw.srt), Cells(merged.cbn.srt))
merged.cbn.srt <- merged.cbn.srt[, consensus_cell]
consensus_gene <- intersect(rownames(merged.raw.srt), rownames(merged.cbn.srt))
merged.cbn.srt <- DietSeurat(merged.cbn.srt, features = consensus_gene)

rm(merged.raw.srt)
gc()

## Preprocess ## Slow for large data set
merged.cbn.srt$pct_mito <- NULL
merged.cbn.srt <- RenameAssays(merged.cbn.srt, 'RNA' = 'CBN')
assay <- 'CBN'
merged.cbn.srt <- Preprocess(merged.cbn.srt,  assay)
merged.cbn.srt <- FindVariableFeatures(merged.cbn.srt) |>
        ScaleData(vars.to.regress = c('pct_mito_CBN', 'nCount_CBN', 'nFeature_CBN'))
merged.cbn.srt <- RunUMAP.my(merged.cbn.srt, assay = 'CBN', var.toal = 0.75, haromize.by = 'sample')
gc()


merged.cbn.srt <- DietSeurat(merged.cbn.srt, assays = 'CBN', dimreducs = names(merged.cbn.srt@reductions))
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
saveRDS(merged.cbn.srt, 'integrated/PART02.merged.cbn.srt.rds')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
