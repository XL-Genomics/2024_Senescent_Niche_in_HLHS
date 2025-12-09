####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Pediatric HLHS Compendium
####  2022-09-01 by Xiao LI (Texas Heart Institute, US)
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Initiate  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
Ver <- '2'
Step <- 'PART32_Pan_Disease_Mapping'
Project <- '2022_hlhs_dturaga'
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Senescence in Adult MI + Lavine Adult DCM  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Adult MI Dataset
x <- c('Cardiomyocyte_snRNA_snATAC.Rds',
       'Endothelial_snRNA_snATAC.Rds',
       'Fibroblast_snRNA_snATAC.Rds',
       'Lymphoid_snRNA_snATAC.Rds',
       'Myeloid_snRNA_snATAC.Rds',
       'Neuronal_snRNA_snATAC.Rds',
       'Pericyte_snRNA_snATAC.Rds',
       'vSMCs_snRNA_snATAC.Rds')
y <- list()
for(i in 1:L(x)){
        message(x[i])
        y[[i]] <- readRDS(paste0('/Volumes/shire/data/scrnaseq/2022_Nature_RKramann/matrix_public/', x[i]))
        y[[i]] <- y[[i]][, y[[i]]$tech == 'RNA']
        y[[i]] <- DietSeurat(y[[i]], assays = 'RNA',
                             dimreducs = names(y[[i]]@reductions), layers = c('counts', 'data'))
}
mi.srt <- merge(y[[1]], y[2:L(x)], merge.dr = T)
mi.srt$group1 <- as.character(mi.srt$patient_region_id)
mi.srt$group2 <- paste0('Kramann_', mi.srt$region)
mi.srt$donor <- mi.srt$patient
mi.srt$Cell_type <- revalue(mi.srt$cell_type, replace = c(
        Endo = 'EC',
        Fib = 'FB',
        Lymphoid = 'Lym',
        Myeloid = 'Mye',
        Neuronal = 'NC',
        vSMCs = 'SMC'
))
mi.srt$Cell_state <- mi.srt$annotation
mi.srt$orig.name <- str_split(Cells(mi.srt), '-', simplify = T)[, 1]
mi.srt$library <- mi.srt$orig.ident
mi.srt$study <- 'RKramann_MI'
mi.srt$patient_group <- revalue(mi.srt$patient_group, replace = c(
        'group_1' = 'Myogenic',
        'group_2' = 'Ischemic',
        'group_3' = 'Fibrotic'))
mi.srt@meta.data <- mi.srt@meta.data[, c(
        'library', 'orig.name', 'study', 'Cell_state', 'Cell_type', 'group2', 'group1', 'donor', 'patient_group'
)]
mi.srt@reductions$sub_umap <- mi.srt@reductions$umap_harmony_v2
mi.srt@reductions$sub_umap@key <- 'cSUBUMAP_'
colnames(mi.srt@reductions$sub_umap@cell.embeddings) <- paste0('cSUBUMAP_', 1:2)

## Adult DCM Dataset
tadros.dcm <- readRDS('../../../2022_pediatric_general/rdata/human_v0/external/HJTadros_Peds_DCM//Final.srt.rds')
tadros.dcm <- DropMetaLevels(tadros.dcm[, tadros.dcm$group2 %in% c('Adult_Ctrl', 'pre-LVAD') &
                                                tadros.dcm$Cell_type %in% c('Cardiomyocyte',
                                                                            'Endothelial',
                                                                            'Endocardium',
                                                                            'Fibroblast',
                                                                            'Myeloid')])
tadros.dcm$Cell_type <- revalue(tadros.dcm$Cell_type, replace = c(
        'Cardiomyocyte' = 'CM',
        'Endothelial' = 'EC',
        'Endocardium' = 'EC',
        'Fibroblast' = 'FB',
        'Myeloid' = 'Mye'
))
tadros.dcm <- RenameAssays(tadros.dcm, 'CBN' = 'RNA')
tadros.dcm <- DietSeurat(tadros.dcm, layers = c('counts', 'data'))
tadros.dcm <- UpdateSeuratObject(tadros.dcm)
tadros.dcm <- NormalizeData(tadros.dcm)

## Current Pediatric HLHS Dataset
hlhs.srt <- RenameAssays(full.srt, 'CBN', 'RNA')
hlhs.srt$study <- 'Current'
pan_hf_v2.srt <- merge(hlhs.srt, list(mi.srt, tadros.dcm), merge.dr = T)
pan_hf_v2.srt <- DietSeurat(pan_hf_v2.srt, dimreducs = 'sub_umap')

pan_hf_v2.srt <- JoinLayers(pan_hf_v2.srt)
pan_hf_v2.srt <- merge(hlhs.srt, list(tadros.dcm, mi.srt))
pan_hf_v2.srt$group3 <- factor(pan_hf_v2.srt$group3, levels = c(
        'Ped_Control',
        'Ped_HLHS',
        'Kramann_Myogenic',
        'Kramann_Ischemic',
        'Kramann_Fibrotic',
        'Adult_Ctrl',
        'Adult_DCM'
))

####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
saveRDS(pan_hf_v2.srt, 'integrated/PART32.ped_hlhs_with_adult_mi_adult_dcm.srt.rds')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
