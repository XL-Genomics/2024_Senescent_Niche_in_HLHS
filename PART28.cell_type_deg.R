####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Pediatric HLHS Compendium
####  2022-09-01 by Xiao LI (Texas Heart Institute, US)
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Initiate  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
Ver <- '2'
Step <- 'PART28_Cell_Type_DEG'
Project <- '2022_hlhs_dturaga'
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Load data  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
full.srt <- readRDS('integrated/PART19.annotated.srt.rds')
full.srt <- full.srt[, full.srt$Non_ambiguous & full.srt$group2 != 'VAD']
full.srt$group2 <- droplevels(full.srt$group2)
srt_list <- SplitObject(full.srt, 'Cell_type')
ct <- names(srt_list)
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Compute sc HLHS vs Control DEG  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
sc_deg_list <- list()
for(i in 1:L(ct)){
        message(ct[i])
        Idents(srt_list[[i]]) <- 'group2'
        deg <- FindMarkers(srt_list[[i]], ident.1 = 'HLHS', ident.2 = 'Control')
        deg <- deg[deg$p_val_adj < 0.05, ]
        deg$gene <- rownames(deg)
        sc_deg_list[[i]] <- deg
        names(sc_deg_list)[i] <- ct[i]
        ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        saveRDS(sc_deg_list, 'analysis/PART28.sc_deg_hlhs_vs_control.srt_mk.rds')
        ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
}
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Compute pseudobulk DEG  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
psb_deg_list <- list()
for(i in 1:L(ct)){
        message(ct[i])
        srt_list[[i]] <- RenameAssays(srt_list[[i]], 'CBN' = 'RNA')
        psb <- AggregateExpression(srt_list[[i]], assays = 'RNA', slot = 'count', group.by = 'group1')$RNA
        meta <- srt_list[[i]]@meta.data[!duplicated(srt_list[[i]]$group1), ]
        rownames(meta) <- meta$group1
        meta <- meta[colnames(psb), ]
        dds <- DESeqDataSetFromMatrix(psb, colData = meta, design = ~ group2)
        dds <- estimateSizeFactors(dds)
        dds <- DESeq(dds)
        res <- results(dds, alpha = 0.05, independentFiltering = T, contrast = c('group2', 'HLHS', 'Control'))
        res <- lfcShrink(dds, contrast = c('group2', 'HLHS', 'Control'), res = res, type = "ashr")
        deg_up <- rownames(res)[res$padj < 0.05 & res$log2FoldChange > 0]
        deg_dn <- rownames(res)[res$padj < 0.05 & res$log2FoldChange < 0]
        deg_up <- deg_up[!is.na(deg_up)]
        deg_dn <- deg_dn[!is.na(deg_dn)]
        psb_deg_list[[i]] <- list('UP' = deg_up, 'DN' = deg_dn)
        names(psb_deg_list)[i] <- ct[i]
        ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        saveRDS(psb_deg_list, 'analysis/PART28.psb_deg_hlhs_vs_control.list.rds')
        ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
}
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Intersect sc and psb DEGs  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
final_deg_list <- list()
for(i in 1:L(ct)){
        message(ct[i])
        sc_deg <- sc_deg_list[[i]]
        final_deg_list[[i]] <- sc_deg[(sc_deg$gene %in% psb_deg_list[[i]][['UP']] & sc_deg$avg_log2FC > 0) |
                                              (sc_deg$gene %in% psb_deg_list[[i]][['DN']] & sc_deg$avg_log2FC < 0), ]
        gl <- split(final_deg_list[[i]]$gene, final_deg_list[[i]]$avg_log2FC > 0)
        print(str(gl))
        names(final_deg_list)[i] <- ct[i]
        ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        saveRDS(final_deg_list, 'analysis/PART28.final_deg_hlhs_vs_control.srt_mk.rds')
        ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
}
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Calculate cell state markers for each cell type  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
plan("multisession", workers = 8)
cell_state_mk <- list()
ct <- ct[! ct %in% c('NC', 'AC', 'PC', 'SMC')]
ct[6] <- 'Mural'
srt_list2 <- srt_list[! names(srt_list) %in% c('NC', 'AC', 'PC', 'SMC')]
srt_list2[[6]] <- full.srt[, full.srt$Cell_type %in% c('PC', 'SMC')]
names(srt_list2)[6] <- 'Mural'

PlotPDF('1.heat.cell_state_marker_v2', 10, 10)
for(i in 1:L(ct)){
        message(ct[i])
        Idents(srt_list2[[i]]) <- 'Cell_state'
        mk <- FindAllMarkers(srt_list2[[i]], logfc.threshold = 0.25, only.pos = T, return.thresh = 0.05)
        mk <- mk[mk$p_val_adj < 0.05, ]
        print(MarkerHeatmap(srt_list2[[i]], mk, top = 10, disp.min = 0, raster = T))
        cell_state_mk[[i]] <- mk
        names(cell_state_mk)[i] <- ct[i]
        ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        saveRDS(cell_state_mk, 'analysis/PART28.cell_state_marker_v2.srt_mk.rds')
        ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
}
dev.off()
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

