####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Pediatric HLHS Compendium
####  2022-09-01 by Xiao LI (Texas Heart Institute, US)
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Initiate  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Ver <- '2'
Step <- 'PART50_Visium_ST_Data'
Project <- '2022_hlhs_dturaga'
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Load snRNA data  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
full.srt <- readRDS('integrated/PART19.annotated.srt.rds')
sn.meta <- U(full.srt@meta.data[, c('tissue', 'sex', 'age', 'donor', 'vad_pair', 'vad_duration', 'palliation',
                                    'rv_depression', 'avv_regurg', 'BNP', 'log10BNP', 'group1', 'group2')])
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Collecting HLHS data  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Sample 1 ####
LoadST <- function(name, sample){
        srt <- Load10X_Spatial(
                data.dir = paste0('/Volumes/shire/data/visium/2022_Pediatric_frozen_DTuraga/matrix/', name, '/outs/'),
                filename = "filtered_feature_bc_matrix.h5",
                assay = "ST",
                slice = sample,
                filter.matrix = T)
}
srt <- LoadST('2022_Pediatric_frozen_DTuraga_FC3CB', 'Ctrl_02')
srt@images$Ctrl_02@coordinates$tissue <- as.integer(srt@images$Ctrl_02@coordinates$tissue)
srt@images$Ctrl_02@coordinates$row <- as.integer(srt@images$Ctrl_02@coordinates$row)
srt@images$Ctrl_02@coordinates$col <- as.integer(srt@images$Ctrl_02@coordinates$col)
srt@images$Ctrl_02@coordinates$imagerow <- as.integer(srt@images$Ctrl_02@coordinates$imagerow)
srt@images$Ctrl_02@coordinates$imagecol <- as.integer(srt@images$Ctrl_02@coordinates$imagecol)
srt$Orig_name <- Cells(srt)
srt$Sample <- 'Ctrl_02'
srt$Library <- '2022_Pediatric_frozen_DTuraga_FC3CB'
srt$orig.ident <- NULL
srt$Coord_x_slide <- srt@images$Ctrl_02@coordinates$imagerow
srt$Coord_y_slide <- srt@images$Ctrl_02@coordinates$imagecol
srt <- RenameCells(srt, new.names = paste0(srt$Library, ':', as.vector(srt$Orig_name)))
Ctrl_02 <- srt

#### Sample 2 ####
LoadST <- function(name, sample){
        srt <- Load10X_Spatial(
                data.dir = paste0('/Volumes/shire/data/visium/2023_HLHS_FFPE_DTuraga/matrix/', name, '/outs/'),
                filename = "filtered_feature_bc_matrix.h5",
                assay = "ST",
                slice = sample,
                filter.matrix = T)
}
srt <- LoadST('2023_HLHS_FFPE_DTuraga_P184', 'HLHS_03')
srt@images$HLHS_03@coordinates$tissue <- as.integer(srt@images$HLHS_03@coordinates$tissue)
srt@images$HLHS_03@coordinates$row <- as.integer(srt@images$HLHS_03@coordinates$row)
srt@images$HLHS_03@coordinates$col <- as.integer(srt@images$HLHS_03@coordinates$col)
srt@images$HLHS_03@coordinates$imagerow <- as.integer(srt@images$HLHS_03@coordinates$imagerow)
srt@images$HLHS_03@coordinates$imagecol <- as.integer(srt@images$HLHS_03@coordinates$imagecol)
srt$Orig_name <- Cells(srt)
srt$Sample <- 'HLHS_03'
srt$Library <- '2023_HLHS_FFPE_DTuraga_P184'
srt$orig.ident <- NULL
srt$Coord_x_slide <- srt@images$HLHS_03@coordinates$imagerow
srt$Coord_y_slide <- srt@images$HLHS_03@coordinates$imagecol
srt <- RenameCells(srt, new.names = paste0(srt$Library, ':', as.vector(srt$Orig_name)))
HLHS_03 <- srt

#### Sample 3 ####
srt <- LoadST('2023_HLHS_FFPE_DTuraga_P181', 'HLHS_13')
srt@images$HLHS_13@coordinates$tissue <- as.integer(srt@images$HLHS_13@coordinates$tissue)
srt@images$HLHS_13@coordinates$row <- as.integer(srt@images$HLHS_13@coordinates$row)
srt@images$HLHS_13@coordinates$col <- as.integer(srt@images$HLHS_13@coordinates$col)
srt@images$HLHS_13@coordinates$imagerow <- as.integer(srt@images$HLHS_13@coordinates$imagerow)
srt@images$HLHS_13@coordinates$imagecol <- as.integer(srt@images$HLHS_13@coordinates$imagecol)
srt$Orig_name <- Cells(srt)
srt$Sample <- 'HLHS_13'
srt$Library <- '2023_HLHS_FFPE_DTuraga_P181'
srt$orig.ident <- NULL
srt$Coord_x_slide <- srt@images$HLHS_13@coordinates$imagerow
srt$Coord_y_slide <- srt@images$HLHS_13@coordinates$imagecol
srt <- RenameCells(srt, new.names = paste0(srt$Library, ':', as.vector(srt$Orig_name)))
HLHS_13 <- srt

#### Sample 5  ####
srt <- LoadST('2023_HLHS_FFPE_DTuraga_P125_P151', 'HLHS_14')
srt@images$HLHS_14@coordinates$tissue <- as.integer(srt@images$HLHS_14@coordinates$tissue)
srt@images$HLHS_14@coordinates$row <- as.integer(srt@images$HLHS_14@coordinates$row)
srt@images$HLHS_14@coordinates$col <- as.integer(srt@images$HLHS_14@coordinates$col)
srt@images$HLHS_14@coordinates$imagerow <- as.integer(srt@images$HLHS_14@coordinates$imagerow)
srt@images$HLHS_14@coordinates$imagecol <- as.integer(srt@images$HLHS_14@coordinates$imagecol)
srt$Orig_name <- Cells(srt)
srt$Coord_x_slide <- srt@images$HLHS_14@coordinates$imagerow
srt$Coord_y_slide <- srt@images$HLHS_14@coordinates$imagecol
srt$orig.ident <- NULL
srt$Library <- '2023_HLHS_FFPE_DTuraga_P125_P151'
srt <- RenameCells(srt, new.names = paste0(srt$Library, ':', as.vector(srt$Orig_name)))
SpatialDimPlot(srt, cells.highlight = Cells(srt)[srt$Coord_x_slide < 950])

HLHS_14 <- srt[, srt$Coord_x_slide < 950]
HLHS_14$Sample <- 'HLHS_14'

#### Sample 7  ####
srt <- LoadST('2023_HLHS_FFPE_DTuraga_P125_P151', 'VAD_05')
srt@images$VAD_05@coordinates$tissue <- as.integer(srt@images$VAD_05@coordinates$tissue)
srt@images$VAD_05@coordinates$row <- as.integer(srt@images$VAD_05@coordinates$row)
srt@images$VAD_05@coordinates$col <- as.integer(srt@images$VAD_05@coordinates$col)
srt@images$VAD_05@coordinates$imagerow <- as.integer(srt@images$VAD_05@coordinates$imagerow)
srt@images$VAD_05@coordinates$imagecol <- as.integer(srt@images$VAD_05@coordinates$imagecol)
srt$Orig_name <- Cells(srt)
srt$Coord_x_slide <- srt@images$VAD_05@coordinates$imagerow
srt$Coord_y_slide <- srt@images$VAD_05@coordinates$imagecol
srt$orig.ident <- NULL
srt$Library <- '2023_HLHS_FFPE_DTuraga_P125_P151'
srt <- RenameCells(srt, new.names = paste0(srt$Library, ':', as.vector(srt$Orig_name)))
SpatialDimPlot(srt, cells.highlight = Cells(srt)[srt$Coord_x_slide > 950])

VAD_05 <- srt[, srt$Coord_x_slide > 950]
VAD_05$Sample <- 'VAD_05'

#### Sample 6 ####
LoadST <- function(name, sample){
        srt <- Load10X_Spatial(
                data.dir = paste0('/Volumes/shire/data/visium/2023_HLHS_FFPE2_DTuraga/matrix/', name, '/outs/'),
                filename = "filtered_feature_bc_matrix.h5",
                assay = "ST",
                slice = sample,
                filter.matrix = T)
}
srt <- LoadST('2023_HLHS_FFPE2_DTuraga_P197', 'VAD_04')
srt@images$VAD_04@coordinates$tissue <- as.integer(srt@images$VAD_04@coordinates$tissue)
srt@images$VAD_04@coordinates$row <- as.integer(srt@images$VAD_04@coordinates$row)
srt@images$VAD_04@coordinates$col <- as.integer(srt@images$VAD_04@coordinates$col)
srt@images$VAD_04@coordinates$imagerow <- as.integer(srt@images$VAD_04@coordinates$imagerow)
srt@images$VAD_04@coordinates$imagecol <- as.integer(srt@images$VAD_04@coordinates$imagecol)
srt$Orig_name <- Cells(srt)
srt$Sample <- 'VAD_04'
srt$Library <- '2023_HLHS_FFPE2_DTuraga_P197'
srt$orig.ident <- NULL
srt$Coord_x_slide <- srt@images$VAD_04@coordinates$imagerow
srt$Coord_y_slide <- srt@images$VAD_04@coordinates$imagecol
srt <- RenameCells(srt, new.names = paste0(srt$Library, ':', as.vector(srt$Orig_name)))
VAD_04 <- srt

#### SCT normalize each slide ####
srt.list <- list(Ctrl_02,
                 HLHS_03,
                 HLHS_13,
                 HLHS_14,
                 VAD_04,
                 VAD_05)
for(i in 1:L(srt.list)){
        srt.list[[i]] <- SCTransform(srt.list[[i]],
                                     assay = 'ST',
                                     method = "glmGamPoi",
                                     vars.to.regress = c('nCount_ST', 'nFeature_ST'),
                                     do.correct.umi = T,
                                     seed.use = 505,
                                     return.only.var.genes = F)
}
#### Merge dataset ####
merged <- merge(srt.list[[1]], srt.list[2:L(srt.list)])
VariableFeatures(merged) <- U(c(VariableFeatures(srt.list[[1]], assay = 'SCT'),
                                VariableFeatures(srt.list[[2]], assay = 'SCT'),
                                VariableFeatures(srt.list[[3]], assay = 'SCT'),
                                VariableFeatures(srt.list[[4]], assay = 'SCT'),
                                VariableFeatures(srt.list[[5]], assay = 'SCT'),
                                VariableFeatures(srt.list[[6]], assay = 'SCT')))


#### Store spatial coordinates as a reduced dimension ####
merged@reductions$spatial <- CreateDimReducObject(
        embeddings = as.matrix(data.frame(Spatial_1 = merged$Coord_x_slide,
                                          Spatial_2 = merged$Coord_y_slide,
                                          row.names = Cells(merged))),
        assay = 'ST',
        key = 'Spatial_')
merged@misc$spot_scale <- rep(1, 10)
merged$Sample <- factor(merged$Sample)
Idents(merged) <- 'Sample'

DefaultAssay(merged) <- 'ST'
merged <- NormalizeData(merged)
merged <- JoinLayers(merged)
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  QC spots  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
VlnPlot2(merged, features = c('nFeature_ST', 'nCount_ST'))
## Set filter for max nCount
cutoff_upper <- c()
spot_toss <- c()
for(i in 1:L(levels(merged$Sample))){
        print(i)
        cells_in_group <- Cells(merged)[merged$Sample == levels(merged$Sample)[i]]
        cutoff_upper[i] <- GetOutlier(merged$nCount_ST[cells_in_group], iqr_multiplier = 1.5)[2]
        spot_toss <- c(spot_toss,
                       cells_in_group[merged$nCount_ST[cells_in_group] > cutoff_upper[i]])
}
L(spot_toss) ## 237

spot_toss2 <- Cells(merged)[merged$nCount_ST < 200]

merged <- merged[, ! Cells(merged) %in% c(spot_toss2, spot_toss)]
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
saveRDS(merged, 'integrated/PART50.visium_st_v2.srt.rds')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
