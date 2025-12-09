####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Pediatric HLHS Compendium
####  2022-09-01 by Xiao LI (Texas Heart Institute, US)
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Initiate  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
Ver <- '2'
Step <- 'PART26_CellChat'
Project <- '2022_hlhs_dturaga'
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Load Seurat  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
full.srt <- readRDS('integrated/PART19.annotated.srt.rds')
full.srt <- DropMetaLevels(full.srt[, full.srt$Non_ambiguous & full.srt$group2 != 'VAD'])
Table(full.srt$Cell_state, full.srt$Cell_type)
Idents(full.srt) <- 'Cell_state'
ctrl.srt <- DropMetaLevels(full.srt[, full.srt$group2 == 'Control'])
hlhs.srt <- DropMetaLevels(full.srt[, full.srt$group2 == 'HLHS'])

tmp.srt <- RenameAssays(full.srt, 'CBN', 'RNA')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Prepare CellChat database  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~---
DoCellChat <- function(srt, group.by, assay, secreted_only = F, trim = 0.1){
        CellChatDB <- CellChatDB.my
        if(secreted_only){CellChatDB <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation")}
        ####  Build CellChat Object  ####
        cch <- createCellChat(object = srt, group.by = group.by, assay = assay)
        cch <- addMeta(cch, meta = srt@meta.data)
        groupSize <- as.numeric(table(cch@idents))
        cch@DB <- CellChatDB
        ####  Preprocessing the expression data for cell-cell communication analysis  ####
        cch <- subsetData(cch) # subset the expression data of signaling genes for saving computation cost
        cch <- identifyOverExpressedGenes(cch)
        cch <- identifyOverExpressedInteractions(cch)
        cch <- projectData(cch, PPI.human)
        ####  Inference of cell-cell communication network  ####
        cch <- computeCommunProb(cch, type = "truncatedMean", trim = trim)
        # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
        cch <- filterCommunication(cch)
        # CellChat computes the communication probability on signaling pathway level
        cch <- computeCommunProbPathway(cch)
        cch <- aggregateNet(cch)
        return(cch)
}
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  CellChat workflow for CBN assay  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
ctrl_cell_type.cch <- DoCellChat(ctrl.srt, group.by = 'Cell_type', assay = 'CBN', secreted_only = T, trim = 0.01)
hlhs_cell_type.cch <- DoCellChat(hlhs.srt, group.by = 'Cell_type', assay = 'CBN', secreted_only = T, trim = 0.01)

full_cell_type.cch <- mergeCellChat(list(Control = ctrl_cell_type.cch,
                                         HLHS = hlhs_cell_type.cch),
                                    add.names = c('Control', 'HLHS'))

ctrl_cell_state.cch <- DoCellChat(ctrl.srt, group.by = 'Cell_state', assay = 'CBN', secreted_only = T, trim = 0.01)
hlhs_cell_state.cch <- DoCellChat(hlhs.srt, group.by = 'Cell_state', assay = 'CBN', secreted_only = T, trim = 0.01)

full_cell_state.cch <- mergeCellChat(list(Control = ctrl_cell_state.cch,
                                          HLHS = hlhs_cell_state.cch),
                                     add.names = c('Control', 'HLHS'))
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
saveRDS(ctrl_cell_type.cch, 'analysis/PART26.ctrl_cell_type.cellchat.rds')
saveRDS(hlhs_cell_type.cch, 'analysis/PART26.hlhs_cell_type.cellchat.rds')
saveRDS(ctrl_cell_state.cch, 'analysis/PART26.ctrl_cell_state.cellchat.rds')
saveRDS(hlhs_cell_state.cch, 'analysis/PART26.hlhs_cell_state.cellchat.rds')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Re-calculate CellChat focusing on senescent niche cells as senders   ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
full.srt <- readRDS('integrated/PART19.annotated_v2.srt.rds')
tmp.srt <- DropMetaLevels(full.srt[, full.srt$Non_ambiguous & full.srt$group2 == 'HLHS'])
tmp.srt$tmp <- factor(tmp.srt$Cell_type, levels = c(levels(tmp.srt$Cell_type), 'FB5', 'EC4', 'PC2', 'MP3'))
tmp.srt$tmp[tmp.srt$Cell_state %in% c('FB5', 'EC4', 'PC2', 'MP3')] <-
        tmp.srt$Cell_state[tmp.srt$Cell_state %in% c('FB5', 'EC4', 'PC2', 'MP3')]
Table(tmp.srt$tmp, tmp.srt$Cell_type)

tmp.srt <- DropMetaLevels(tmp.srt[, ! tmp.srt$tmp %in% c('NC', 'AC')])

hlhs_cell_type_sen.cch <- DoCellChat(tmp.srt, group.by = 'tmp', assay = 'CBN', secreted_only = T, trim = 0.01)
saveRDS(hlhs_cell_type_sen.cch, 'analysis/PART26.hlhs_cell_type_sen.cellchat.rds')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Filter HLHS-up ligand based on senescence niche cells   ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
full.srt <- readRDS('integrated/PART19.annotated.srt.rds')
tmp.srt <- DropMetaLevels(full.srt[, full.srt$Non_ambiguous & full.srt$group2 == 'HLHS'])
all_lg <- U(hlhs_cell_state.cch@LR$LRsig$ligand)
psb <- AverageExpression(tmp.srt[, tmp.srt$Cell_state %in% c('FB5', "EC4", 'MP3', 'PC2') ], 
                         features = all_lg, group.by = 'Cell_state')$CBN
sen_lig <- data.frame(FB5 = rep(NA, L(all_lg)), EC4 = NA, MP3 = NA, PC2 = NA)
rownames(sen_lig) <- all_lg
for(i in 1:L(all_lg)) {
        sen_lig[i, ] <- c(all_lg[i] %in% cs_deg$FB$gene[cs_deg$FB$cluster == 'FB5'],
                          all_lg[i] %in% cs_deg$EC$gene[cs_deg$EC$cluster == 'EC4'],
                          all_lg[i] %in% cs_deg$Mye$gene[cs_deg$Mye$cluster == 'MP3'],
                          all_lg[i] %in% cs_deg$Mural$gene[cs_deg$Mural$cluster == 'PC2'])
}
sen_lig <- sen_lig[rowSums(sen_lig) > 0, ]

cs_deg <- readRDS('analysis/PART28.cell_state_marker.srt_mk.rds')
fb5_lig <- intersect(cs_deg$FB$gene[cs_deg$FB$cluster == 'FB5'],  all_lg)
ec4_lig <- intersect(cs_deg$EC$gene[cs_deg$EC$cluster == 'EC4'],  all_lg)
mp3_lig <- intersect(cs_deg$Mye$gene[cs_deg$Mye$cluster == 'MP3'],  all_lg)
pc2_lig <- intersect(cs_deg$Mural$gene[cs_deg$Mural$cluster == 'PC2'],  all_lg)

psb <- psb[rownames(sen_lig), ]

ph_result <- pheatmap(as.matrix(psb), scale = 'row', 
                      clustering_method = 'ward.D2', clustering_distance_rows = 'euclidean',
                      cellwidth = 10, cellheight = 10)
psb <- psb[ph_result$tree_row$order,
           ph_result$tree_col$order]
sen_lig <- sen_lig[rownames(psb), colnames(psb)]
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
saveRDS(psb, 'analysis/PART26.sen_niche_ligands.dataframe.rds')
saveRDS(sen_lig, 'analysis/PART26.sen_niche_ligands_state_marker_bin.dataframe.rds')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
