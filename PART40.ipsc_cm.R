####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Pediatric HLHS Compendium
####  2022-09-01 by Xiao LI (Texas Heart Institute, US)
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Initiate  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Ver <- '2'
Step <- 'PART40_iPSC_CM'
Project <- '2022_hlhs_dturaga'
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  iPSC-CM  GSE146341  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
data_path <- '/Volumes/shire/data/scrnaseq/2022_CellStemCell_CLo/matrix/'
ips.srt.list1 <- list(
    CreateSeuratObject(
        Read10X_h5(paste0(data_path, '2022_CellStemCell_CLo_Healthy/outs/filtered_feature_bc_matrix.h5')),
        project = 'iPSCCM_Healthy'),
    CreateSeuratObject(
        Read10X_h5(paste0(data_path, '2022_CellStemCell_CLo_HLHS-Good-A/outs/filtered_feature_bc_matrix.h5')),
        project = 'iPSCCM_HLHS_Good_1'),
    CreateSeuratObject(
        Read10X_h5(paste0(data_path, '2022_CellStemCell_CLo_HLHS-Good-B//outs/filtered_feature_bc_matrix.h5')),
        project = 'iPSCCM_HLHS_Good_2'),
    CreateSeuratObject(
        Read10X_h5(paste0(data_path, '2022_CellStemCell_CLo_HLHS-Poor/outs/filtered_feature_bc_matrix.h5')),
        project = 'iPSCCM_HLHS_Poor'))
for(i in 1:4) {
    ips.srt.list1[[i]] <- RenameCells(ips.srt.list1[[i]],
                                      new.names = paste0('2022_CellStemCell_CLo:',
                                                         ips.srt.list1[[i]]$orig.ident, ':',
                                                         Cells(ips.srt.list1[[i]]))
    )
}
ips.srt1 <- merge(ips.srt.list1[[1]], ips.srt.list1[2:4])

data_path <- '/Volumes/shire/data/scrnaseq/2022_CellStemCell_CLo/cellbender_v1/'
ips.srt.list2 <- list(
    CreateSeuratObject(
        Read10X_h5(paste0(data_path, '2022_CellStemCell_CLo_Healthy/cellbender_filtered.h5')),
        project = 'iPSCCM_Healthy'),
    CreateSeuratObject(
        Read10X_h5(paste0(data_path, '2022_CellStemCell_CLo_HLHS-Good-A/cellbender_filtered.h5')),
        project = 'iPSCCM_HLHS_Good_1'),
    CreateSeuratObject(
        Read10X_h5(paste0(data_path, '2022_CellStemCell_CLo_HLHS-Good-B//cellbender_filtered.h5')),
        project = 'iPSCCM_HLHS_Good_2'),
    CreateSeuratObject(
        Read10X_h5(paste0(data_path, '2022_CellStemCell_CLo_HLHS-Poor/cellbender_filtered.h5')),
        project = 'iPSCCM_HLHS_Poor'))
for(i in 1:4) {
    ips.srt.list2[[i]] <- RenameCells(ips.srt.list2[[i]],
                                      new.names = paste0('2022_CellStemCell_CLo:',
                                                         ips.srt.list2[[i]]$orig.ident, ':',
                                                         Cells(ips.srt.list2[[i]]))
    )
}
ips.srt2 <- merge(ips.srt.list2[[1]], ips.srt.list2[2:4])

O(Cells(ips.srt1), Cells(ips.srt2))
ips.srt <- ips.srt2[, intersect(Cells(ips.srt2), Cells(ips.srt1))]
ips.srt <- RenameAssays(ips.srt, 'RNA' = 'CBN')

head(ips.srt)
ips.srt$Group1 <- ips.srt$orig.ident
ips.srt$Group2 <- revalue(ips.srt$orig.ident, replace = c('iPSCCM_Healthy' = 'Control',
                                                          'iPSCCM_HLHS_Good_1' = 'HLHS_Good',
                                                          'iPSCCM_HLHS_Good_2' = 'HLHS_Good',
                                                          'iPSCCM_HLHS_Poor' = 'HLHS_Poor'))
Table(ips.srt$Group1, ips.srt$Group2)

####  Quality controls  ####
ips.srt <- PercentageFeatureSet(ips.srt, pattern = '^MT-', col.name = 'pct_mito_CBN')
ips.srt$pct_mito_CBN[is.na(ips.srt$pct_mito_CBN)] <- 0
VlnPlot2(ips.srt, group.by = 'Group1', features = c('nCount_CBN', 'nFeature_CBN', 'pct_mito_CBN'))

####  DR  ####
ips.srt <- NormalizeData(ips.srt) |>
    FindVariableFeatures() |>
    ScaleData(vars.to.regress = c('nCount_CBN', 'nFeature_CBN', 'pct_mito_CBN')) |>
    RunPCA() |>
    RunHarmony(group.by.vars = 'Group1')
ips.srt <- RunUMAP(ips.srt, reduction = 'harmony', dims = 1:30)

####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  iPSC-CM  compare with HLHS  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sen_ageatlas <- list(as.vector(read.table('external/senescence_list_from_AgingAtlas.txt', header = F)[, 1]))
sasp_react <- list(as.vector(read.table('external/SASP_reactome_R-HSA-2559582.tsv')[,1]))
sasp_senmayo <- list(as.vector(read.table('external/sasp_senmayo.csv')[,1]))

## Compute Senescence Scores
ips.srt <- AddModuleScore2(ips.srt, features = c(sen_ageatlas, sasp_react, sasp_senmayo), return_z = T,
                           names = c('sen_ageatlas', 'sasp_react', 'sasp_senmayo'))

## Compute HLHS CM State Signature Scores
cm_mk <- readRDS('analysis/PART28.cell_state_marker.srt_mk.rds')
cm_mk <- cm_mk$CM
cm_mk <- split(cm_mk$gene, cm_mk$cluster)
ips.srt <- AddModuleScore2(ips.srt, features = cm_mk, return_z = T, names = names(cm_mk))
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ips.srt <- DietSeurat(ips.srt, dimreducs = names(ips.srt@reductions))
saveRDS(ips.srt, 'individual/PART40.2022_CellStemCell_CLo.srt.rds')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----



####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  iPSC-CM  GSE146763  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
data_path <- 'external/GSE146763/'
ips.srt.list1 <- list(
    CreateSeuratObject(
        Read10X(paste0(data_path, 'GSM4405513_Control/')),
        project = 'iPSCCM_Healthy'),
    CreateSeuratObject(
        Read10X(paste0(data_path, 'GSM4405514_HLHS_146')),
        project = 'iPSCCM_HLHS'))
for(i in 1:2) {
    ips.srt.list1[[i]] <- RenameCells(ips.srt.list1[[i]],
                                      new.names = paste0('2020_Circulation_SWu:',
                                                         ips.srt.list1[[i]]$orig.ident, ':',
                                                         Cells(ips.srt.list1[[i]]))
    )
}
ips.srt <- merge(ips.srt.list1[[1]], ips.srt.list1[[2]])
ips.srt <- JoinLayers(ips.srt)
ips.srt <- RenameAssays(ips.srt, 'RNA' = 'CBN')

ips.srt$Group1 <- ips.srt$orig.ident
ips.srt$Group2 <- revalue(ips.srt$orig.ident, replace = c('iPSCCM_Healthy' = 'Control',
                                                          'iPSCCM_HLHS' = 'HLHS'))
Table(ips.srt$Group1, ips.srt$Group2)

####  DR  ####
ips.srt <- NormalizeData(ips.srt) |>
    FindVariableFeatures() |>
    ScaleData(vars.to.regress = c('nCount_CBN', 'nFeature_CBN', 'pct_mito_CBN')) |>
    RunPCA() |>
    RunHarmony(group.by.vars = 'Group1')
ips.srt <- RunUMAP(ips.srt, reduction = 'harmony', dims = 1:30)


ips.srt <- ips.srt[, ips.srt$CM_identity == 'CM']

####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  iPSC-CM  compare with HLHS  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Compute Senescence Scores
ips.srt <- AddModuleScore2(ips.srt, features = c(sen_ageatlas, sasp_react, sasp_senmayo), return_z = T,
                           names = c('sen_ageatlas', 'sasp_react', 'sasp_senmayo'))

## Compute HLHS CM State Signature Scores
cm_mk <- readRDS('analysis/PART28.cell_state_marker.srt_mk.rds')
cm_mk <- cm_mk$CM
cm_mk <- split(cm_mk$gene, cm_mk$cluster)
ips.srt <- AddModuleScore2(ips.srt, features = cm_mk, return_z = T, names = names(cm_mk))

####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ips.srt <- DietSeurat(ips.srt, dimreducs = names(ips.srt@reductions))
saveRDS(ips.srt, 'individual/PART40.2020_Circulation_SWu.srt.rds')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----



####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  iPSC-CM  GSE135411  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
data_path <- 'external/GSE135411/'
ips.srt.list3 <- list(
    CreateSeuratObject(counts =
                           read.table(gzfile(paste0(data_path, 'GSM4008017_CTRL_1_only_analysis_TPM_matrix.txt.gz'))),
                       project = 'iPSCCM_Ctrl1'),
    CreateSeuratObject(counts =
                           read.table(gzfile(paste0(data_path, 'GSM4008018_CTRL_2_only_analysis_TPM_matrix.txt.gz'))),
                       project = 'iPSCCM_Ctrl2'),
    CreateSeuratObject(counts =
                           read.table(gzfile(paste0(data_path, 'GSM4008019_HLHS_1_only_analysis_TPM_matrix.txt.gz'))),
                       project = 'iPSCCM_HLHS1'),
    CreateSeuratObject(counts =
                           read.table(gzfile(paste0(data_path, 'GSM4008020_HLHS_2_only_analysis_TPM_matrix.txt.gz'))),
                       project = 'iPSCCM_HLHS2')
    )
for(i in 1:4) {
    ips.srt.list3[[i]] <- RenameCells(ips.srt.list3[[i]],
                                      new.names = paste0('2021_Circulation_AMoretti:',
                                                         ips.srt.list3[[i]]$orig.ident, ':',
                                                         Cells(ips.srt.list3[[i]]))
    )
}
ips.srt <- merge(ips.srt.list3[[1]], ips.srt.list3[2:4])
ips.srt <- JoinLayers(ips.srt)

head(ips.srt)
tail(ips.srt)
ips.srt$Group1 <- paste0(str_split(str_split(Cells(ips.srt), ':', simplify = T)[, 3], '_', simplify = T)[, 1], '_',
                         str_split(str_split(Cells(ips.srt), ':', simplify = T)[, 3], '_', simplify = T)[, 2])
ips.srt$Group2 <- str_split(str_split(Cells(ips.srt), ':', simplify = T)[, 3], '_', simplify = T)[, 1]
Table(ips.srt$Group1, ips.srt$Group2)


####  DR  ####
ips.srt <- NormalizeData(ips.srt) |>
    FindVariableFeatures() |>
    ScaleData(vars.to.regress = c('nCount_RNA', 'nFeature_RNA', 'pct_mito_RNA')) |>
    RunPCA() |>
    RunHarmony(group.by.vars = 'Group1')
ips.srt <- RunUMAP(ips.srt, reduction = 'harmony', dims = 1:30)

ips.srt <- ips.srt[, ips.srt$CM_identity == 'CM']

####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  iPSC-CM  compare with HLHS  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Compute Senescence Scores
ips.srt <- AddModuleScore2(ips.srt, features = c(sen_ageatlas, sasp_react, sasp_senmayo), return_z = T,
                           names = c('sen_ageatlas', 'sasp_react', 'sasp_senmayo'))

## Compute HLHS CM State Signature Scores
cm_mk <- readRDS('analysis/PART28.cell_state_marker.srt_mk.rds')
cm_mk <- cm_mk$CM
cm_mk <- split(cm_mk$gene, cm_mk$cluster)
ips.srt <- AddModuleScore2(ips.srt, features = cm_mk, return_z = T, names = names(cm_mk))

####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ips.srt <- DietSeurat(ips.srt, dimreducs = names(ips.srt@reductions))
saveRDS(ips.srt, 'individual/PART40.2021_Circulation_AMoretti.srt.rds')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
