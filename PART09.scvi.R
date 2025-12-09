####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Pediatric HLHS Compendium
####  2022-09-01 by Xiao LI (Texas Heart Institute, US)
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Initiate  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
Ver <- '2'
Step <- 'PART09_scVI'
Project <- '2022_hlhs_dturaga'
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Load data  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
srt <- readRDS('integrated/PART05.merged.clean_cbn.srt.rds')
DefaultAssay(srt) <- 'CBN'
srt <- DietSeurat(srt, dimreducs = 'hmn_umap')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Follow jupyter notebook PART07 and PART08  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Re-embed based on scVI latent space  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
obsm <- read.csv(file = 'integrated/PART08.scanvi_integration.csvs/obsm.csv', sep = ',', header = T)
obs <- read.csv(file = 'integrated/PART08.scanvi_integration.csvs/obs.csv', sep = ',', header = T)
rownames(obsm) <- str_remove(obs$X, pattern = '-0$')
O(Cells(srt), rownames(obsm))
identical(Cells(srt), rownames(obsm))

obsm2 <- read.csv(file = 'integrated/PART07.scvi_integration.csvs/obsm.csv', sep = ',', header = T)
obs2 <- read.csv(file = 'integrated/PART07.scvi_integration.csvs/obs.csv', sep = ',', header = T)
rownames(obsm2) <- str_remove(obs2$X, pattern = '-0$')
O(Cells(srt), rownames(obsm2))
identical(Cells(srt), rownames(obsm2))

####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
DefaultAssay(srt) <- 'CBN'
## Save reduced dimensions from scVI (no ref data)
srt@reductions$scVI_umap <- srt@reductions$umap
srt@reductions$scVI_umap@cell.embeddings <- as.matrix(obsm2[, paste0('X_scVIumap', 1:2)])
colnames(srt@reductions$scVI_umap@cell.embeddings) <- c('SCVIUMAP_1', 'SCVIUMAP_2')
srt@reductions$scVI_umap@assay.used <- 'CBN'
srt@reductions$scVI_umap@key <- 'SCVIUMAP_'
DimPlot2(srt, reduction = 'scVI_umap', group.by = 'group2')

srt@reductions$scVI <- srt@reductions$umap
srt@reductions$scVI@cell.embeddings <- as.matrix(obsm2[, paste0('X_scVI', 1:50)])
rownames(srt@reductions$scVI@cell.embeddings) <- Cells(srt)
colnames(srt@reductions$scVI@cell.embeddings) <- paste0('SCVI_', 1:50)
srt@reductions$scVI@assay.used <- 'CBN'
srt@reductions$scVI@key <- 'SCVI_'

## Save reduced dimensions from scANVI
srt@reductions$scANVI_umap <- srt@reductions$umap
srt@reductions$scANVI_umap@cell.embeddings <- as.matrix(obsm[, paste0('X_scANVIumap', 1:2)])
colnames(srt@reductions$scANVI_umap@cell.embeddings) <- c('SCANVIUMAP_1', 'SCANVIUMAP_2')
srt@reductions$scANVI_umap@assay.used <- 'CBN'
srt@reductions$scANVI_umap@key <- 'SCANVIUMAP_'
DimPlot2(srt, reduction = 'scANVI_umap', group.by = 'group2')

srt@reductions$scANVI <- srt@reductions$umap
srt@reductions$scANVI@cell.embeddings <- as.matrix(obsm[, paste0('X_scANVI', 1:50)])
rownames(srt@reductions$scANVI@cell.embeddings) <- Cells(srt)
colnames(srt@reductions$scANVI@cell.embeddings) <- paste0('SCANVI_', 1:50)
srt@reductions$scANVI@assay.used <- 'CBN'
srt@reductions$scANVI@key <- 'SCANVI_'

## Save metadata
srt$Cell_type_scANVI <- obs$scANVI_predict_celltype
srt$Cell_state_scANVI <- obs$scANVI_predict_cellstate

####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
saveRDS(srt, 'integrated/PART09.scvi_integration.srt.rds')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
