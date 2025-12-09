####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Pediatric HLHS Compendium
####  2022-09-01 by Xiao LI (Texas Heart Institute, US)
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Initiate  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
Ver <- '2'
Step <- 'PART19_Consolidation'
Project <- '2022_hlhs_dturaga'

Color_cell_type <- c(mycol_40[c(3, 6, 11, 14, 35, 17, 22, 26, 30, 33)], 'grey85', 'grey')
Color_cell_state <- c(mycol_40[1:36], 'grey40', 'grey85')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Load data  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
srt <- readRDS('integrated/PART08.scvi_integration.srt.rds')

srt$Cell_state <- NA
srt$Cell_type <- as.vector(srt$Cell_type)

Table(srt$Cell_state, srt$Cell_type)
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Organize annotation order  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
srt.bkp <- srt
srt <- srt.bkp
Table(srt$Cell_state, srt$Cell_type)

srt$Cell_type <- factor(srt$Cell_type,
                        levels = c('CM', 'FB', 'EC', 'PC', 'SMC',
                                   'Lym', 'Mye', 'EpiC', 'NC', 'AC',
                                   'MitoHi', 'Doublet'))
srt$Cell_state <- factor(srt$Cell_state, levels = c(
        'CM1', 'CM2', 'CM3', 'CM4', 'CM5',
        'FB1', 'FB2', 'FB3', 'FB4', 'FB5',
        'EC1', 'EC2', 'EC3', 'EC4', 'LEC', 'EndoC',
        'PC1', 'PC2', 'SMC1', 'SMC2',
        'TC', 'Prol TC', 'NK',
        'Mono', 'MP1', 'MP2', 'MP3', 'MP4', 'Prol MP', 'DC', 'Mast',
        'EpiC', 'NC', 'AC',
        'MitoHi', 'Doublet'
))
Table(srt$Cell_state, srt$Cell_type)

srt$Cell_type[srt$Cell_state == 'Doublet'] <- 'Doublet'
srt$Cell_type[srt$Cell_state %in% c('PC1', 'PC2')] <- 'PC'
srt$Cell_type[srt$Cell_state %in% c('SMC1', 'SMC2')] <- 'SMC'
srt$Cell_state[srt$Cell_type %in% c('EpiC', 'NC', 'AC')] <- srt$Cell_type[srt$Cell_type %in% c('EpiC', 'NC', 'AC')]
srt$Cell_state[srt$Cell_type == 'Doublet'] <- 'Doublet'
srt$Cell_type[srt$Cell_state %in% c('CM5')] <- 'MitoHi'
srt$Cell_state[srt$Cell_state %in% c('CM5')] <- 'MitoHi'
srt$Cell_state[is.na(srt$Cell_state)] <- 'MitoHi'
srt$Cell_type[srt$Cell_state == 'MitoHi'] <- 'MitoHi'

Table(srt$Cell_state, srt$Cell_type)


## EpiC and MitoHi are ambiguous cells that likely originated from single donor or are technical artifacts
srt$Non_ambiguous <- T
srt$Non_ambiguous[srt$Cell_type %in% c('Doublet', 'MitoHi')] <- F
srt$Non_ambiguous[srt$Cell_state %in% c('Doublet', 'MitoHi', 'EpiC')] <- F
Table(srt$Non_ambiguous, srt$Cell_type)
Table(srt$Non_ambiguous, srt$Cell_state)
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Re-embed global UMAP without ambiguous cells  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
clean.srt <- srt[, srt$Non_ambiguous]
clean.srt <- RunUMAP(clean.srt, reduction = 'scVI', dims = 1:50, n.neighbors = 50, min.dist = 0.5,
                     reduction.name = 'clean_umap', reduction.key = 'cUMAP_')

## Clean umap without ambiguouse cells
srt@reductions$clean_umap <- srt@reductions$scVI2_umap
srt@reductions$clean_umap@cell.embeddings[,1] <- NA
srt@reductions$clean_umap@cell.embeddings[,2] <- NA
srt@reductions$clean_umap@cell.embeddings[Cells(clean.srt), 1:2] <-
        clean.srt@reductions$clean_umap@cell.embeddings[,1:2]
srt@reductions$clean_umap@key <- 'cUMAP_'
colnames(srt@reductions$clean_umap@cell.embeddings) <- c('cUMAP_1', 'cUMAP_2')

## Full umap with doublets and mito-high cells
srt@reductions$full_umap <- srt@reductions$scVI2_umap
srt@reductions$full_umap@key <- 'fUMAP_'
colnames(srt@reductions$full_umap@cell.embeddings) <- c('fUMAP_1', 'fUMAP_2')

## Sub-umap for all cell types
srt@reductions$sub_full_umap <- srt@reductions$sub_scVI_umap
srt@reductions$sub_full_umap@key <- 'fSUBUMAP_'
colnames(srt@reductions$sub_full_umap@cell.embeddings) <- c('fSUBUMAP_1', 'fSUBUMAP_2')

## Sub-umap for all non-ambiguous cell types
srt@reductions$sub_umap <- srt@reductions$sub_clean_scVI_umap
srt@reductions$sub_umap@key <- 'cSUBUMAP_'
colnames(srt@reductions$sub_umap@cell.embeddings) <- c('cSUBUMAP_1', 'cSUBUMAP_2')

## Sub-umap for Main Mono/MP and VEC non-ambiguous cells
srt@reductions$sub_main_umap <- srt@reductions$sub_main_scVI_umap
srt@reductions$sub_main_umap@key <- 'mSUBUMAP_'
colnames(srt@reductions$sub_main_umap@cell.embeddings) <- c('mSUBUMAP_1', 'mSUBUMAP_2')


all_reductions <- srt@reductions

srt <- DietSeurat(srt, counts = T, data = T, scale.data = F,
                  assays = 'CBN',
                  dimreducs = c('sub_umap', 'sub_full_umap', 'sub_main_umap', 'full_umap', 'clean_umap'))
srt@reductions <- srt@reductions[c('clean_umap', 'full_umap', 'sub_full_umap', 'sub_main_umap', 'sub_umap')]
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
saveRDS(all_reductions, 'integrated/PART19.all_reductions.srt_dimreducs.rds')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Adding Clinical Data into Meta  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
srt$vad_pair <- sample_meta[srt$library, 'VAD_pair']
srt$vad_pair[srt$vad_pair == 'Pair_1' & srt$group2 == 'HLHS'] <- 'Pair 1 Pre-VAD'
srt$vad_pair[srt$vad_pair == 'Pair_2' & srt$group2 == 'HLHS'] <- 'Pair 2 Pre-VAD'
srt$vad_pair[srt$vad_pair == 'Pair_3' & srt$group2 == 'HLHS'] <- 'Pair 3 Pre-VAD'
srt$vad_pair[srt$vad_pair == 'Pair_4' & srt$group2 == 'HLHS'] <- 'Pair 4 Pre-VAD'
srt$vad_pair[srt$vad_pair == 'Pair_1' & srt$group2 == 'VAD'] <- 'Pair 1 Post-VAD'
srt$vad_pair[srt$vad_pair == 'Pair_2' & srt$group2 == 'VAD'] <- 'Pair 2 Post-VAD'
srt$vad_pair[srt$vad_pair == 'Pair_3' & srt$group2 == 'VAD'] <- 'Pair 3 Post-VAD'
srt$vad_pair[srt$vad_pair == 'Pair_4' & srt$group2 == 'VAD'] <- 'Pair 4 Post-VAD'
srt$vad_pair <- factor(srt$vad_pair, levels = c(
        'Pair 1 Pre-VAD', 'Pair 1 Post-VAD',
        'Pair 2 Pre-VAD', 'Pair 2 Post-VAD',
        'Pair 3 Pre-VAD', 'Pair 3 Post-VAD',
        'Pair 4 Pre-VAD', 'Pair 4 Post-VAD'
))
srt$vad_duration <- as.numeric(sample_meta[srt$library, 'VAD_duration'])
srt$palliation <- factor(sample_meta[srt$library, 'Palliation'], levels = c('None', 'I', 'II', 'III', 'VAD'))
srt$palliation[is.na(srt$palliation)] <- 'None'
srt$rv_depression <- as.numeric(sample_meta[srt$library, 'RV_function'])
srt$avv_regurg <- as.numeric(sample_meta[srt$library, 'AVV_regurg'])
srt$log10BNP <- log10(sample_meta[srt$library, 'BNP'])

meta <- srt@meta.data
meta <- meta[, c(
        'sample',
        'library',
        'orig.name',
        'study',
        'method',
        'platform',
        'protocol',
        'data_process',
        'tissue',
        'enrichment',
        'preparation',
        'diagnosis',
        'sex',
        'age',
        'donor',
        'replicate',
        'nCount_CBN',
        'nFeature_CBN',
        'pct_mito_CBN',
        'Doublet_SC',
        'Doublet_SC_score',
        'Cell_type',
        'Cell_state',
        'Cell_type_prev',
        'Cell_state_prev',
        'Non_ambiguous',
        'vad_pair',
        'vad_duration',
        'palliation',
        'rv_depression',
        'avv_regurg',
        'BNP',
        'log10BNP',
        'group1',
        'group2'
)]
srt@meta.data <- meta
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Markers  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
clean.srt <- srt[, srt$Non_ambiguous]
clean.srt$Cell_type <- droplevels(clean.srt$Cell_type)
Table(clean.srt$Cell_type)

## Compute Cell_type markers
Idents(clean.srt) <- 'Cell_type'
Cell_type_marker <- FindAllMarkers(clean.srt,
                                   assay = 'CBN',
                                   test.use = 'wilcox',
                                   only.pos = T,
                                   max.cells.per.ident = 1000,
                                   logfc.threshold = 1,
                                   random.seed = 505)

Cell_type_marker.dedup <- Cell_type_marker |> group_by(cluster) |> top_n(25, wt = avg_log2FC)
Cell_type_marker.dedup <- Cell_type_marker.dedup[!duplicated(Cell_type_marker.dedup$gene),]
Table(Cell_type_marker.dedup$cluster)


## Compute Cell_state markers
clean.srt$Cell_state <- droplevels(clean.srt$Cell_state)
Table(clean.srt$Cell_state)
Idents(clean.srt) <- 'Cell_state'
Cell_state_marker <- FindAllMarkers(clean.srt,
                                    assay = 'CBN',
                                    test.use = 'wilcox',
                                    only.pos = T,
                                    max.cells.per.ident = 1000,
                                    logfc.threshold = 0.25,
                                    random.seed = 505)

Cell_state_marker.dedup <- Cell_state_marker |> group_by(cluster) |> top_n(20,  wt = avg_log2FC)
Cell_state_marker.dedup <- Cell_state_marker.dedup[!duplicated(Cell_state_marker.dedup$gene),]
sort(Table(Cell_state_marker.dedup$cluster))

srt@misc$marker <- list()
srt@misc$marker$Cell_type_marker <- Cell_type_marker
# srt@misc$marker$Cell_state_marker <- Cell_state_marker
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
saveRDS(srt, 'integrated/PART19.annotated.srt.rds')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

