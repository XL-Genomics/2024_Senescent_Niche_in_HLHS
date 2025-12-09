####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Pediatric HLHS Compendium
####  2022-09-01 by Xiao LI (Texas Heart Institute, US)
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Initiate  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
Ver <- '2'
Step <- 'PART93_Plot_R_V3'
Project <- '2022_hlhs_dturaga'
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Load and edit previously saved objects  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
#### Library meta sheet, all libraries -- lib_meta.df
lib_meta.df <- read.csv(paste0(Docu_dir, 'pediatric_sample_meta_v2.csv'))
sample_meta.df <- lib_meta.df[!duplicated(lib_meta.df$Group1), ]
sample_meta.df <- sample_meta.df[order(sample_meta.df$Group1), ]
rownames(sample_meta.df) <- sample_meta.df$label

#### All cell types, all samples -- full.srt  ####
full.srt <- readRDS('integrated/PART19.annotated.srt.rds')
Idents(full.srt) <- 'Cell_type'
DefaultAssay(full.srt) <- 'CBN'

#### All cell types, all samples, no ambiguous -- full_clean.srt  ####
full_clean.srt <- DropMetaLevels(full.srt[, full.srt$Non_ambiguous])

#### All cell types, paired only, no ambiguous -- paired.srt  ####
paired.srt <- DropMetaLevels(full.srt[, full.srt$Non_ambiguous & !is.na(full.srt$vad_pair)])

#### Clean CM cells -- cm.srt  ####
cm.srt <- DropMetaLevels(full_clean.srt[, full_clean.srt$Cell_type == 'CM'])
Idents(cm.srt) <- 'Cell_state'

#### Clean CM cells with Bulk Controls -- cm_bulk.srt  ####
cm_bulk.srt <- readRDS('integrated/PART30.cm_and_bulk_integrated.srt.rds')

#### Clean Fibro cells -- fb.srt  ####
fb.srt <- DropMetaLevels(full_clean.srt[, full_clean.srt$Cell_type == 'FB'])
Idents(fb.srt) <- 'Cell_state'

#### Clean Myeloid cells -- mye.srt  ####
mye.srt <- full_clean.srt[, full_clean.srt$Cell_type == 'Mye']
Idents(mye.srt) <- 'Cell_state'

#### Clean EC cells -- ec.srt  ####
ec.srt <- full_clean.srt[, full_clean.srt$Cell_type == 'EC']
Idents(ec.srt) <- 'Cell_state'

#### Clean Mural cells -- mu.srt  ####
mu.srt <- full_clean.srt[, full_clean.srt$Cell_type %in% c('SMC', 'PC')]
Idents(mu.srt) <- 'Cell_state'

#### Clean Cell State Markers -- state_mk.list  ####
state_mk.list <- readRDS('analysis/PART28.cell_state_marker_v2.srt_mk.rds')

#### HLHS vs Control DEG -- hlhs_deg.list  ####
hlhs_deg.list <- readRDS('analysis/PART28.final_deg_hlhs_vs_control.srt_mk.rds')

#### CellChat  ####
ctrl.cch <- readRDS('analysis/PART26.ctrl_cell_type.cellchat.rds')
hlhs.cch <- readRDS('analysis/PART26.hlhs_cell_type.cellchat.rds')
ctrl.cs.cch <- readRDS('analysis/PART26.ctrl_cell_state.cellchat.rds')
hlhs.cs.cch <- readRDS('analysis/PART26.hlhs_cell_state.cellchat.rds')
hlhs.sen.cch <- readRDS('analysis/PART26.hlhs_cell_type_sen.cellchat.rds')

#### Integrated Visium ST -- st.srt  ####
st.srt <- readRDS('integrated/PART51.visium_st_niche.srt.rds')
st.srt$Sample <- factor(st.srt$Sample,
                        levels = c('Ctrl_02', 'HLHS_03', 'HLHS_14', 'HLHS_13', 'VAD_04', 'VAD_05'))
st.srt@images <- st.srt@images[levels(st.srt$Sample)]
spot_sizes <- c(0.1, 0.5, 0.6, 0.2, 0.2, 0.8)
asp_ratio <- c(1.123, 1.394, 1.029, 1.024, 0.836, 1.666)
DefaultAssay(st.srt) <- 'ST'


#### Integrated published iPSC-CM scRNA-seq -- ips.srt  ####
ips1.srt <- readRDS('individual/PART40.2022_CellStemCell_CLo.srt.rds')
ips1.srt <- RenameAssays(ips1.srt, 'CBN', 'RNA')
ips1.srt$Study <- 'CLo'
ips1.srt$Group3 <- factor(str_split(ips1.srt$Group2, pattern = '_', simplify = T)[, 1], levels = c('Control', 'HLHS'))
ips2.srt <- readRDS('individual/PART40.2020_Circulation_SWu.srt.rds')
ips2.srt <- RenameAssays(ips2.srt, 'CBN', 'RNA')
ips2.srt[["RNA"]] <- as(object = ips2.srt[["RNA"]], Class = "Assay")
ips2.srt <- NormalizeData(ips2.srt)
ips2.srt$Study <- 'SWu'
ips2.srt$Group3 <- factor(revalue(ips2.srt$Group1, replace = c('iPSCCM_Healthy' = 'Control',  'iPSCCM_HLHS' = 'HLHS')),
                          levels = c('Control', 'HLHS'))
ips3.srt <- readRDS('individual/PART40.2021_Circulation_AMoretti.srt.rds')
ips3.srt$Study <- 'AMoretti'
ips3.srt$Group3 <- factor(revalue(ips3.srt$Group2, replace = c('CTRL' = 'Control')), levels = c('Control', 'HLHS'))
ips.srt <- merge(ips1.srt, list(ips2.srt, ips3.srt))
ips.srt$Group1 <- paste(ips.srt$Study, ips.srt$Group1)
table(ips.srt$Group1, ips.srt$Group3)
ips.srt$Group4 <- paste(ips.srt$Study, ips.srt$Group3)
ips.srt <- JoinLayers(ips.srt)

#### Xenuim -- xen.srt  ####
xen.srt <- readRDS('integrated/PART65.xenium_annotated_no_image.srt.rds')
xen.srt <- DropMetaLevels(xen.srt[, xen.srt$Cell_state != 'Doublet'])
xen.srt$Sen_cell <- ifelse(xen.srt$Cell_state %in% c('CM3', 'FB5', 'EC4', 'MP3', 'PC2'), 'Sen', 'NonSen')
xen.srt@reductions$spatial <- xen.srt@reductions$clean_umap
xen.srt@reductions$spatial@cell.embeddings[,1] <- xen.srt$Coord_x
xen.srt@reductions$spatial@cell.embeddings[,2] <- xen.srt$Coord_y
colnames(xen.srt@reductions$spatial@cell.embeddings) <- paste0('spatial_', 1:2)
xen.srt@reductions$spatial@key <- 'spatial_'
xen_img.srt <- readRDS('integrated/PART65.xenium_annotated_with_image.srt.rds')



#### Gene Sets ####
sen_ageatlas <- list(as.vector(read.table('external/senescence_list_from_AgingAtlas.txt', header = F)[, 1]))
sasp_react <- list(as.vector(read.table('external/SASP_reactome_R-HSA-2559582.tsv')[,1]))
sasp_senmayo <- list(as.vector(read.table('external/sasp_senmayo.csv')[,1]))
sasp <- list(U(c(unlist(sasp_react), unlist(sasp_senmayo))))
sen_combine <- list(U(c(unlist(sen_ageatlas), unlist(sasp))))
hif1a_target <- list(GetKeggGene('hsa04066'))
vegfa_target <- read.table('external/ABE_VEGFA_TARGETS.v2022.1.Hs.grp', comment.char = '#', header = T)[,1]
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Figure 1   ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
##~~  Plot 1  Bubble - sample distribution  ####
df <- U(full_clean.srt@meta.data[, c('group1', 'group2', 'sex', 'age', 'palliation', 'diagnosis', 'vad_pair', 'donor')])
df$type <- factor('snRNA-seq', levels = c('snRNA-seq', 'CM Bulk RNA-seq', 'snRNA-seq + ST'))
df <- rbind(df, data.frame(group1 = paste('Bctrl', 1:4),
                           group2 = 'Control',
                           sex = c('F', 'M')[c(2,1,2,2)],
                           age = c(0.05, 0.2, 2, 4),
                           palliation = 'None',
                           diagnosis = 'Normal',
                           vad_pair = NA,
                           donor = paste('Bctrl', 1:4),
                           type = 'CM Bulk RNA-seq'))
rownames(df) <- df$group1
df$age <- as.numeric(as.vector(df$age))
df$type[df$donor %in% st.srt$donor] <- 'snRNA-seq + ST'
df$type[df$donor %in% xen.srt$Donor] <- 'snRNA-seq + ST'
p <- ggplot(df) +
        geom_point(aes(x = age, y = palliation, shape = sex, color = type),
                   position = position_jitter(width = 0, height = 0.2), size = 3) +
        scale_color_manual(values = c('#D93A3D', '#B2E069', '#6ABCCE')) +
        scale_shape_manual(values = c(16, 15)) +
        # scale_size(breaks = c(1e4, 4e4, 7e4), limits = c(1e2, 7e4)) +
        scale_y_discrete(limits = rev) +
        scale_x_continuous(breaks = seq(0, 17, 2), limits = c(0, 17)) +
        labs(x = ' ', y = 'Palliation Stage', shape = 'Data Type', color = 'Sex') +
        theme_classic() +
        theme(aspect.ratio = 1,
              axis.ticks = element_blank(),
              axis.line = element_line(size = 0.5, color = 'black'),
              legend.background = element_blank(),
              panel.background = element_blank(),
              panel.grid.major = element_line(size = 0.5, color = 'grey75'))
p
PlotPDF('1.01.bubble.stage_sex_age', 10, 10)
p
dev.off()



##~~  Plot 2  UMAP - global HLHS status  ####
tmp.srt <- full_clean.srt[, unlist(DownsampleByMeta(full_clean.srt, meta_var = 'palliation', down_to_min_group = T))]
p <- DimPlot2(tmp.srt, group.by = 'palliation', cols = Color_palliation,
              reduction = 'clean_umap', raster = F, pt.size = 0.1) +
        labs(x = 'UMAP 1', y = 'UMAP 2', title = 'Integration of Control and HLHS')
p
PlotPDF('1.02.umap.global_hlhs_status', 15, 15)
p
dev.off()


##~~  Plot 3  UMAP - global cell type  ####
p <- DimPlot2(full_clean.srt, group.by = 'Cell_type', cols = Color_cell_type, reduction = 'clean_umap',
              label = F, raster = F, pt.size = 0.1) +
        labs(x = 'UMAP 1', y = 'UMAP 2', title = 'Cell Types')
p
PlotPDF('1.03.umap.global_cell_type', 15, 15)
p
dev.off()


##~~  Plot 4  Scatter - Global PCA by gender and age  ####
df <- readRDS('analysis/PART27.ctrl_hlhs_all_cell_pseudobulk_pca.df.rds')
df <- df[df$name %in% full_clean.srt$group1, ]
df2 <- U(full_clean.srt@meta.data[, c('age', 'sex', 'group1')])
df2 <- df2[order(df2$group1), ]
df$Age <- as.numeric(as.vector(df2$age))
df$Sex <- df2$sex
df$Count <- as.vector(Table(full_clean.srt$group1))
p1 <- ggplot(df) +
        geom_point(aes(x = PC1, y = PC2, color = group)) +
        scale_color_manual(values = Color_palliation) +
        labs(color = '', x = '', y = '') +
        labs(x = 'PC1 16% variance', y = 'PC2 12% variance')
p2 <- ggplot(df) +
        geom_point(aes(x = PC1, y = PC2, color = Age, shape = Sex)) +
        scale_shape_manual(values = c(16, 15)) +
        scale_color_distiller(palette = 'Spectral', breaks = seq(0, 17, 3)) +
        labs(x = '', y = '') +
        labs(x = 'PC1 16% variance', y = 'PC2 12% variance')
p <- wrap_plots(p1, p2, ncol = 2) &
        theme(aspect.ratio = 1,
              axis.line = element_line(size = 0.5, color = 'black'),
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              legend.background = element_blank(),
              panel.background = element_blank())
p
PlotPDF('1.04.pca.sample_metadata', 6, 3)
p
dev.off()



##~~  Plot 5  Bar - DEG counts per cell type  ####
x <- readRDS('analysis/PART28.psb_deg_hlhs_vs_control.list.rds')
data <- data.frame(Cell_type = factor(rep(names(x), 2),
                                      levels = c('CM','FB','EC','PC','SMC','Lym','Mye','NC','AC')),
                   N = NA, Direction = rep(c('HLHS Up', 'HLHS Down'), each = L(names(x))))
for(i in 1:L(x)){
        data[i, 2] <- L(x[[i]]$UP)
        data[i+L(x), 2] <- L(x[[i]]$DN)
}
data <- data |> group_by(Cell_type) |> mutate(Total = sum(N))
data$Cell_type <- factor(data$Cell_type, levels = as.vector(U(data$Cell_type[order(data$Total, decreasing = T)])))
p <- ggplot(data) +
        geom_bar(aes(y = N, x = Cell_type, fill = Direction),
                 position = "stack", stat = "identity") +
        scale_fill_manual(values = mycol_10) +
        scale_y_continuous(breaks = seq(0, 2000, 500), limits = c(0, 2000)) +
        scale_x_discrete(limits = rev) +
        labs(x = '', y = 'Number of DEGs') +
        coord_flip() +
        theme_classic() +
        theme(aspect.ratio = 1)
p
PlotPDF('1.05.bar.num_psb_deg_per_cell_type', 4, 4)
p
dev.off()



##~~  Plot 6  Visium Feature - HLHS DEG Score    ####
gl <- readRDS('analysis/PART28.psb_deg_hlhs_vs_control.list.rds')
up <- c()
for(i in 1:L(gl)){up <- c(up, gl[[i]]$UP)}
tmp.srt <- AddModuleScore2(st.srt, features = list(U(up)), names = 'HLHSup')
tmp.srt$HLHSup <- SmoothFeature(tmp.srt, feature = 'HLHSup')
p <- FeaturePlotST_Dark(tmp.srt, features = 'HLHSup', minvals = 0, maxvals = 2.5, pt.sizes = spot_sizes, ncol = 6) &
        scale_colour_viridis_c(option = 'magma', begin = 0.1, limits = c(0, 2.5), values = c(0, 0.25, 1))
p <- wrap_plots(p[[1]], p[[2]], p[[3]], p[[4]], p[[5]], p[[6]], ncol = 6)
p <- p & theme(plot.background = element_rect(fill = 'black'),
               panel.background = element_rect(fill = 'black'))
for(i in 1){for(j in 1:6){p[[6*(i-1)+j]] <- p[[6*(i-1)+j]] + theme(aspect.ratio = asp_ratio[j]) + NoLegend()}}
p
PlotPDF('1.06.st_feat.hlhs_up_deg', 15, 2.2)
p
dev.off()



##~~  Plot 3  Box - HLHS-up gene score in ST - NCVR Revision V1    ####
gl <- readRDS('analysis/PART28.psb_deg_hlhs_vs_control.list.rds')
up <- c()
for(i in 1:L(gl)){up <- c(up, gl[[i]]$UP)}
tmp.srt <- AddModuleScore2(st.srt, features = list(U(up)), names = 'HLHSup')
p1 <- BoxPlot2(tmp.srt, feature = 'HLHSup', group.by = 'group1', cols = Color_palliation[c(1:5, 5)],
               ref.group = 'Ctrl_02') +
        labs(x = '', y = 'Expression Z-score')
p1
data <- tmp.srt@meta.data[, c('HLHSup', 'group1')]
data <- split(data$HLHSup, data$group1)
p <- data.frame(pval = c(0.00000000000000022204460492503131,
                         0.00000000000000022204460492503131,
                         0.00000000000000022204460492503131,
                         0.00000000000000022204460492503131,
                         0.00000000000000022204460492503131),
                pair = 1:5)
p$adj_p <- p.adjust(p$pval, method = "BH")

tmp.srt <- AddModuleScore2(xen.srt, features = list(U(up)), names = 'HLHSup', assay = 'Xenium', nbin = 20, ctrl = 20)
p2 <- BoxPlot2(tmp.srt, feature = 'HLHSup', group.by = 'Group1', cols = Color_palliation[c(1:5, 5)],
               ref.group = 'Ctrl_01') +
        labs(x = '', y = 'Expression Z-score')
p2
PlotPDF('1.03.box.hlhs_up_deg_xenium', 4, 3)
p1+p2
dev.off()


##~~  Plot 8  UMAP - Xenium Cell Type  - NCVR Revision V1  ####
p <- DimPlot2(xen.srt, group.by = 'Cell_type', pt.size = 0.1, cols = Color_cell_type,
              raster = F, reduction = 'clean_umap') +
        theme(panel.border = element_rect(fill = NA))
p
PlotPDF('2.08.umap.global_xenium_cell_type', 12, 10)
p
dev.off()


##~~  Plot 4  ST Xenium - Cell type across sections   ####
# tmp.srt <- xen_img.srt
# tmp.srt$Cell_type[tmp.srt$Cell_type %in% c('Doublet', 'LowQual')] <- NA
plist <- list()
sizes <- c(0.5, 0.7, 0.8, 1, 0.6, 1.0)
for(i in 1:6) {
        plist[[i]] <- ImageDimPlot(xen_img.srt, fov = names(xen_img.srt@images)[i], axes = T, size = sizes[i],
                                   cols = Color_Xen,  border.size = NA,
                                   group.by = 'Cell_type', na.value = NA,
                                   cells = Cells(xen_img.srt)[!xen_img.srt$Cell_type %in% c('Doublet', 'LowQual')],
                                   coord.fixed = T)
}
p <- wrap_plots(plist, nrow = 1) &
        theme_Publication(grid = F) &
        theme(panel.background = element_rect(fill = 'black')) & NoLegend()
p[[6]] <- p[[6]] + RestoreLegend()
p
PlotPDF('3.04.st_xenium.spatial_centroid_color_by_cell_type', 40, 10)
p
dev.off()




####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
####  Extended Data Figure 1   ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
##~~  Plot 8  Baton -  BNP level in paired pre/post VAD  ####
data <- sample_meta.df[! is.na(sample_meta.df$VAD_pair),
                       c('BNP', 'RV_function', 'Palliation', "VAD_pair", 'Diagnosis', "VAD_duration")]
data$Delta <- log10(data$BNP[1:4])-log10(data$BNP[5:8]) > 0
data <- data[data$VAD_pair != 'Pair_0', ]
p <- ggplot(data, aes(x = Diagnosis, y = log10(BNP), color = VAD_pair, group = VAD_pair)) +
        geom_point() +
        geom_line() +
        theme_classic() +
        theme(aspect.ratio = 2) +
        scale_color_manual(values = Color_pair)
p
PlotPDF('1.08.baton.paired_bnp_level_in_postvad', 4, 6)
p
dev.off()



##~~  Plot 9  Heatmap Mapping with HCA annotation  ####
p <- MappingHeatmap(full_clean.srt, que_var = 'Cell_type', ref_var = 'Cell_type_scANVI',
                      percentage = T, ref.disp.min = 0.1, center = T,
                      que_order = levels(full_clean.srt$Cell_type)[1:10],
                      ref_order = c(
                              'Ventricular_Cardiomyocyte',
                              'Atrial_Cardiomyocyte',
                              'Fibroblast',
                              'Endothelial',
                              'Pericytes',
                              'Smooth_muscle_cells',
                              'Lymphoid',
                              'Myeloid',
                              'Mesothelial',
                              'Neuronal',
                              'Adipocytes'
                      )) +
        labs(x = 'Adult HCA Annotation', y = 'Integrated Pediatric Data', title = 'Mapping with HCA') +
        theme(aspect.ratio = 10/11)
p
PlotPDF('1.09.heat.cell_type_mapping_with_HCA', 5, 5)
p
dev.off()


##~~  Plot 10  Heatmap - Cell type markers  ####
mk <- full_clean.srt@misc$marker$Cell_type_marker
p <- MarkerHeatmap(full_clean.srt, marker.df = mk, n_cells = 500, top = 30,
                    disp.min = 0.5, disp.max = 2.5, raster = T)+
        theme(aspect.ratio = 1) +
        theme(axis.text.y = element_blank()) +
        scale_fill_viridis() +
        NoLegend()
p
PlotPDF('1.10.heat.cell_type_marker', 10, 10)
p
dev.off()
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
####  Figure 2   ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
##~~  Plot 1  Box - CM young and adult signature scores   ####
mean_exp <- readRDS('analysis/PART30.cm_and_bulk_young_adult_score.srt_meta.rds')
mean_exp <- mean_exp[mean_exp$variable == 'Mean_Score', ]
mean_exp$age_group <- factor(NA, levels = c('0-1', '2-9', '>10'))
mean_exp$age_group[mean_exp$age <= 2] <- '0-1'
mean_exp$age_group[mean_exp$age > 2 & mean_exp$age < 10] <- '2-9'
mean_exp$age_group[mean_exp$age >= 10] <- '>10'
p <- ggplot(mean_exp) +
        geom_boxplot(aes(x = age_group, y = value, fill = group2), outliers = F) +
        geom_quasirandom(aes(x = age_group, y = value)) +
        stat_compare_means(aes(x = age_group, y = value, fill = group2),
                           method = 'wilcox.test', label = 'p.signif') +
        labs(x = 'Chronological age group', y = 'Normalized aging score', color = '', fill = '') +
        scale_fill_manual(values = Color_palliation[1:2]) +
        scale_color_manual(values = Color_palliation[1:2]) +
        theme_classic() +
        theme(aspect.ratio = 1)
p
PlotPDF('2.01.box.young_adult_signature_expression_in_all_cm_with_bulk', 5, 4)
p
dev.off()



##~~  Plot 3  Box - CM EC FB Mye Senescence score HLHS MI DCM   - NCVR Revision V1  ####
mix.srt <- readRDS('integrated/PART32.ped_hlhs_with_adult_mi_adult_dcm.srt.rds')
Table(mix.srt$group4)
mix.srt <- DropMetaLevels(mix.srt[, mix.srt$group4 != 'Kramann_Fibrotic'])
Table(mix.srt$Cell_type, mix.srt$group4)

srt.list <- SplitObject(mix.srt, split.by = 'Cell_type')
p.list <- list()
for(i in 1:4){srt.list[[i]] <- AddModuleScore2(srt.list[[i]], features = sen_combine, names = 'Sen_sig', return_z = T)}
for(i in 1:4){
        p.list[[i]] <- BoxPlot2(srt.list[[i]], features = 'Sen_sig', group.by = 'group4', min = 0.1, max = 0.9,
                                cols = c('grey80', mycol_10[2], 'grey80', mycol_10[5], 'grey80', mycol_10[1])) +
                labs(title = names(srt.list)[i]) +
                scale_y_continuous(limits = c(-1.5, 1.5)) +
                theme(aspect.ratio = 1)
}
p <- wrap_plots(p.list, ncol = 2)
p
PlotPDF('2.03.box.sen_sig_across_4_cell_types_adult_mi_dcm', 6, 6)
p
dev.off()



##~~  Plot 3  UMAP - CM Cell State    ####
p <- DimPlot2(cm.srt, group.by = 'Cell_state', cols = Color_cell_type, reduction = 'sub_umap', label = F, raster = F) +
        labs(x = 'UMAP 1', y = 'UMAP 2', title = 'CM States')
p
PlotPDF('2.03.umap.cm_states', 10, 10)
p
dev.off()



##~~  Plot 4  Feature UMAP - CM Kramman MI Signatures    ####
p <- FeaturePlot2(cm.srt, features = c('Healthy_CM', 'Pre-stressed CM', 'Stressed CM'),
                  min.cutoff = 0, max.cutoff = 3, raster = F, pt.size = 0.1, ncol = 1)
p
PlotPDF('2.04.feature.mi_cm_signatures_score', 4, 12)
p
dev.off()



##~~  Plot 5  Density UMAP - CM enrichemnet over palliation  ####
tmp.srt <- paired.srt
tmp.srt <- cm.srt[, Cells(cm.srt) %in% Cells(tmp.srt) | cm.srt$group2 == 'Control']
tmp.srt <- tmp.srt[, unlist(DownsampleByMeta(tmp.srt, meta_var = 'group2',
                                             down_to_min_group = T, random = T))]
Idents(tmp.srt) <- 'group2'
df2 <- data.frame(tmp.srt@reductions$sub_umap@cell.embeddings)
df2$group2 <- tmp.srt$group2
df3 <- data.frame(tmp.srt@reductions$sub_umap@cell.embeddings)
df3$group2 <- tmp.srt$group2
alpha <- 0.3
smooth <- 3
bins <- 15
p_list <- list()
for(i in 1:3){
        p_list[[i]] <- ggplot(df2, aes(x = cSUBUMAP_1, y = cSUBUMAP_2)) +
                geom_point(data = df3, color = 'black') +
                theme_classic() +
                stat_density_2d(data = df2[df2$group2 == levels(df2$group2)[i],],
                                aes(fill = after_stat(level)),
                                geom = "polygon", contour = T, bins = bins, h = c(smooth, smooth),
                                alpha = alpha, size = 0) +
                scale_fill_viridis() +
                theme(aspect.ratio = 1, axis.text = element_blank(), axis.ticks = element_blank(),
                      axis.title = element_blank()) +
                labs(title = levels(tmp.srt$group2)[i]) +
                NoLegend()
}
p <- wrap_plots(p_list, nrow = 1)
p[[3]] <- p[[3]]+ RestoreLegend()
p
PlotPDF('2.05.density_umap.cm_transition_across_palliation', 13, 3)
p
dev.off()



##~~  Plot 6  Bar - CM2/CM3 state composition  ####
df <- melt(as.matrix(Table(cm.srt$Cell_state, cm.srt$group2)))
df$pct <- df$value/(Table(cm.srt$group2)[df$Var2])
p1 <- ggplot(df[df$Var1 %in% c('CM2', 'CM3'), ], aes(x = Var1, y = pct, fill = Var2)) +
        geom_bar(stat = "identity", position = "dodge") +
        scale_fill_manual(values = Color_palliation) +
        labs(x = '', y = 'Percentage of Total CMs (%)', title = 'Aggregated Abundance') +
        theme_Publication(aspect.ratio = 1) +
        NoLegend()
tmp.srt <- DropMetaLevels(xen.srt[, xen.srt$Cell_type == 'CM'])
df <- reshape2::melt(as.matrix(Table(tmp.srt$Cell_state, tmp.srt$Group2)))
df$pct <- df$value/(Table(tmp.srt$Group2)[df$Var2])
p2 <- ggplot(df[df$Var1 %in% c('CM2', 'CM3'), ], aes(x = Var1, y = pct, fill = Var2)) +
        geom_bar(stat = "identity", position = "dodge") +
        scale_fill_manual(values = Color_palliation) +
        labs(x = '', y = 'Percentage of Total CMs (%)', title = 'Aggregated Abundance') +
        theme_Publication(aspect.ratio = 1) +
        NoLegend()
p <- p1 + p2
PlotPDF('2.06.bar.cm2_cm3_abundance_across_disease', 3, 3)
p
dev.off()



##~~  Plot 09  Box - CM senescence in conditions  - NCVR Revision V1  ####
tmp.srt <- AddModuleScore2(cm.srt[, cm.srt$Cell_state != 'CM4'],
                           features = c(sen_ageatlas, sasp_react, sasp_senmayo),
                           names = c('sen_ageatlas', 'sasp_react', 'sasp_senmayo'))
p <- BoxPlot2(tmp.srt, feature = 'sen_ageatlas', group.by = 'Cell_state') +
        BoxPlot2(tmp.srt, feature = 'sasp_react', group.by = 'Cell_state') +
        BoxPlot2(tmp.srt, feature = 'sasp_senmayo', group.by = 'Cell_state')
p
PlotPDF('1.09.box.cm_state_3_senescence_geneset_score', 6, 3)
p
dev.off()




##~~  Plot 8  Heatmap - CM SASP example    ####
p <- MarkerHeatmapMean(cm.srt,
                       features= c('CDKN2A', 'CDKN1A',  'JUN', 'MMP2', 'PDK4', 'TGFB2',
                                   'GDF15', 'CCN1', 'CCL2', 'CSF1'),
                       group.by = 'Cell_state') +
        scale_fill_distiller(palette = 'RdBu')
p
PlotPDF('2.08.heatmap.cm_state_sasp_example', 6, 6)
p
dev.off()



##~~  Plot 9  Box - Published iPSC-CM Expression of CM markers and Sen Genes   ####
ips.srt <- NormalizeData(ips.srt) |> FindVariableFeatures(nfeatures = 1e4)
ips.srt <- AddModuleScore2(ips.srt, features = c(sen_combine, CM1, CM2, CM3),
                           names = c('Sen_Signature', 'CM1', 'CM2', 'CM3'), return_z = T)
cmp <- list(c('AMoretti Control', 'AMoretti HLHS'),
            c('CLo Control', 'CLo HLHS'),
            c('SWu Control', 'SWu HLHS'))
p1 <- BoxPlot2(ips.srt, feature = 'CM1', group.by = 'Group4', cols = Color_ipsc,
               compare_list = cmp, test = 't.test') +
        NoLegend()
p2 <- BoxPlot2(ips.srt, feature = 'CM2', group.by = 'Group4', cols = Color_ipsc,
               compare_list = cmp, test = 't.test') +
        NoLegend()
p3 <- BoxPlot2(ips.srt, feature = 'CM3', group.by = 'Group4', cols = Color_ipsc,
               compare_list = cmp, test = 't.test')+
        NoLegend()
p4 <- BoxPlot2(ips.srt, feature = 'Sen_Signature', group.by = 'Group4', cols = Color_ipsc,
               compare_list = cmp, test = 't.test') +
        NoLegend()
p <- wrap_plots(p1, p2, p3, p4, ncol = 4) &
        labs(x = '', color = '') &
        scale_y_continuous(limits = c(-2, 3)) &
        theme_Publication(aspect.ratio = 2) &
        RotatedAxis() &
        NoLegend()
p
PlotPDF('2.09.box.ipsc_cm_state_signature_and_senescence', 8, 4)
p
dev.off()


##~~  Plot 4  Spatial scatter - Xenium CM1,2,3 Single Marker   - NCVR Revision V1 ####
x <- LoadXenium('~/Temporary/xenium/output-XETG00106__0004175__sample-2__20240926__181354') ## HLHS16
options(future.globals.maxSize = 10 * 1024^3)
cropped.coords <- Crop(x[["fov"]], x = c(1560, 1800), y = c(1860, 2240), coords = "plot")
x[["zoom"]] <- cropped.coords
DefaultBoundary(x[["zoom"]]) <- "centroids"
p <- ImageDimPlot(
        object = x,
        border.color = NA,
        fov = "zoom",
        molecules = c('NPPB', 'FHL2', 'TCAP'), mols.cols = c('#ffde17', '#00aeef', '#ec008c'),
        mols.size = 0.001, mols.alpha = 1, nmols = 1e5
)
p
PlotPDF('2.04.scatter.xenium_hlhs16_cm_state_single_marker', 8, 10)
p
dev.off()



####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
####  Extended Data Figure 2   ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
##~~  Plot 16  Heatmap - CM pseudobulk HLHS vs Control DEGs   ####
deg <- hlhs_deg.list$CM
up <- deg$gene[deg$avg_log2FC > 0]
dn <- deg$gene[deg$avg_log2FC < 0]
psb <- AverageExpression(cm.srt, group.by = 'group1', features = c(up, dn), assays = 'CBN', return.seurat = T)
psb <- ScaleData(psb[, !grepl('^VAD', Cells(psb))], features = rownames(psb))
df <- melt(GetAssayData(psb, layer = 'scale.data')[c(deg$gene[order(deg$avg_log2FC, decreasing = T)][1:50],
                                                     deg$gene[order(deg$avg_log2FC, decreasing = F)][1:50]), ])
p <- ggplot(df) +
        geom_tile(aes(x = Var2, y = Var1, fill = value)) +
        scale_fill_distiller(palette = 'RdBu') +
        labs(x = '', y = '', fill = 'Scaled expr.') +
        theme_classic() +
        theme(aspect.ratio = 5/2,
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              axis.line = element_blank())
p
PlotPDF('2.16.heat.cm_global_psudobulk_deg', 5, 10)
p
dev.off()


##~~  Plot 17  Heatmap Sim et al Adult Young Signautres  ####
df <- reshape2::melt(GetAssayData(cm_bulk.srt, layer = 'scale.data')[unlist(cm_bulk_young_adult_deg), ])
cutoff <- 2
df$value[df$value > cutoff] <- cutoff
df$value[df$value < -cutoff] <- -cutoff
p <- ggplot(df) +
        geom_tile(aes(x = Var2, y = Var1, fill = value)) +
        scale_fill_distiller(palette = 'RdBu') +
        labs(x = '', y = '', fill = 'Scaled expr.') +
        scale_y_discrete(limit = rev) +
        theme_classic() +
        theme(aspect.ratio = 2,
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              axis.line = element_blank())
p
PlotPDF('2.17.heat.sim_bulk_cm_young_adult_genes', 5, 10)
p
dev.off()


##~~  Plot 11   - CM1/CM2/CM3 state composition in snRNA  ####
mtx <- Table(full_clean.srt$group1, full_clean.srt$Cell_state)
total <- Table(full_clean.srt$group1)
df <- as.data.frame(mtx)
df$Norm <- df$Freq/total[df$Var1]
df$group <- str_split(df$Var1, '_', simplify = T)[, 1]
df <- df[df$Var2 %in% c('CM1', 'CM2', 'CM3'), ]
df$pairs <- NA
df$pairs[df$Var1 %in% c('HLHS_07', 'VAD_02')] <- 'Pair1'
df$pairs[df$Var1 %in% c('HLHS_12', 'VAD_04')] <- 'Pair2'
df$pairs[df$Var1 %in% c('HLHS_16', 'VAD_05')] <- 'Pair3'
p <- ggplot(df, aes(x = group, y = Norm, group = group)) +
        geom_boxplot(outliers = F, aes(fill = NULL)) +
        geom_beeswarm(aes(color = pairs), cex = 6) +
        stat_summary(
                fun = mean,
                geom = "crossbar",
                width = 0.5,
                color = "red",
                aes(group = group)
        ) +
        facet_wrap(~Var2, scale = 'free_y') +
        theme_Publication(aspect.ratio = 1)
p
PlotPDF('1.11.beeswarm.cm1_cm2_cm3_abundance_per_individual', 6, 2)
p
dev.off()
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
####  Figure 3   ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
##~~  Plot 6  Baton - SASP in paired data   - NCVR Revision V1   ####
tmp.srt <- AddModuleScore2(paired.srt, features = sen_combine, names = 'Sen')
tmp.srt$pair <- str_split(tmp.srt$vad_pair, pattern = ' ', simplify = T)[, 2]
data <- tmp.srt@meta.data[, c('pair', 'vad_pair', 'Sen', 'group2')] |>
        group_by(vad_pair) |>
        mutate(Sen_mean = mean(Sen))
data <- U(data[, c('pair', 'vad_pair', 'Sen_mean', 'group2')])
p <- ggplot(data, aes(x = group2, y = Sen_mean, color = pair, group = pair)) +
        geom_point() +
        geom_line() +
        # scale_y_continuous(limits = c(-0.3, 0.3), n.breaks = 6) +
        scale_color_manual(values = Color_pair) +
        theme_Publication(aspect.ratio = 2)
p
PlotPDF('1.06.baton.sen_combine_in_paired_vad', 3, 3)
p
dev.off()



##~~  Plot 5 ST Feature - SASP in visium ST  - NCVR Revision V1  ####
tmp.srt <- AddModuleScore2(st.srt, features = sen_combine, names = 'Sen')
p <- FeaturePlotST_Dark(tmp.srt,
                        features = 'Sen',
                        minvals = 0,
                        maxvals = 1.5,
                        pt.sizes = spot_sizes,
                        ncol = 6) &
        scale_colour_viridis_c(option = 'magma', begin = 0.1, limits = c(0, 1.5), values = c(0, 0.5, 1))
p
PlotPDF('1.05.st_feature.sen_combine_in_st', 15, 2.3)
p
dev.off()


##~~  Plot 4 Dot - Sen in Visium & Xenium ST - NCVR Revision V1   ####
vi <- AddModuleScore2(st.srt, features = sen_combine, names = 'Sen')
vi$Sen_smooth <- SmoothFeature(vi, feature = 'Sen')
vi$Group2 <- vi$group2
vi$Group1 <- vi$group1
vi$Type <- 'Visium'
vi_mean <- vi@meta.data |> group_by(Group1) |> mutate(Sen = mean(Sen_smooth))
vi_mean <- U(vi_mean[, c('Sen', 'Group1', 'Group2', 'Type')])

xe <- AddModuleScore2(xen.srt, features = sen_combine, names = 'Sen', assay = 'Xenium', nbin = 20, ctrl = 20)
xe$Group2
xe$Type <- 'Xenium'
xe_mean <- xe@meta.data |> group_by(Group1) |> mutate(Sen = mean(Sen))
xe_mean <- U(xe_mean[, c('Sen', 'Group1', 'Group2', 'Type')])

data <- rbind(vi_mean, xe_mean)
p1 <- ggplot(data, aes(x = Group2, y = Sen, color = Type)) +
        geom_beeswarm() +
        ggrepel::geom_label_repel(aes(label = Group1)) +
        scale_y_continuous(limits = c(-1,2)) +
        theme_Publication(aspect.ratio = 1)
p1
PlotPDF('1.04.dot.sen_combined_visium_xenium', 3, 3)
p1
dev.off()



##~~  Plot 1  UMAP - Cell states  ####
p <- DimPlot2(full_clean.srt, group.by = 'Cell_state', cols = Color_cell_state, reduction = 'clean_umap',
               label = F, raster = F, pt.size = 0.2, alpha = 1) +
        labs(x = 'UMAP 1', y = 'UMAP 2', title = 'Cell States')
p
PlotPDF('3.01.umap.global_cell_states', 20, 20)
p
dev.off()


##~~  Plot 2  Heatmap - Cell states enrichment in HLHS  ####
include <- c('CM1', 'CM2', 'CM3', 'FB1', 'FB2', 'FB3', 'FB4', 'FB5', 'EC1', 'EC2', 'EC3', 'EC4', 'PC1', 'PC2',
             'SMC1', 'SMC2', 'TC', 'Mono', 'MP1', 'MP2', 'MP3')
tmp.srt <- DropMetaLevels(full_clean.srt[, full_clean.srt$Cell_state %in% include])
mtx <- as.data.frame(Table(tmp.srt$group2, tmp.srt$Cell_state))
total <- as.vector(Table(tmp.srt$Cell_state))
names(total) <- names(Table(tmp.srt$Cell_state))
total2 <- as.vector(Table(tmp.srt$group2))
names(total2) <- names(Table(tmp.srt$group2))
mtx$Freq <- mtx$Freq/total[mtx$Var2]
mtx$Freq <- mtx$Freq*1e6/total2[mtx$Var1]
data <- reshape2::acast(mtx, formula = Var1~Var2)
data <- melt(data)
p <- ggplot(data) +
        geom_tile(aes(y = Var1, x = Var2, fill = value)) +
        # geom_bar(aes(x = Var2, y = value)) +
        scale_fill_viridis_c(values = c(0, 0.8, 1)) +
        scale_y_discrete(limits = rev) +
        theme_classic() +
        labs(x = '', y = '', fill = 'Relative Abundance') +
        theme(aspect.ratio = 3/L(include)) +
        RotatedAxis()
p
PlotPDF('3.02.heat.cell_type_vs_disease_distribution', 7, 3)
p
dev.off()



##~~  Plot 3  Milo UMAP - Global HLHS vs Control, PreVAD vs PostVAD   ####
milo <- readRDS('analysis/PART21.global_hlhs_vs_control_milo_obj_result_list.rds')
p1 <- plotNhoodGraphDA.my(milo[[1]], milo_res = milo[[2]], alpha = 0.5, edge_alpha = 0.1, edge_color = 'grey90') +
        scale_fill_distiller(palette = 'RdBu') +
        theme(aspect.ratio = 1, axis.line = element_line())
milo <- readRDS('analysis/PART21.global_prevad_vs_postvad_milo_obj_result_list.rds')
p2 <- plotNhoodGraphDA.my(milo[[1]], milo_res = milo[[2]], alpha = 0.5, edge_alpha = 0.1, edge_color = 'grey90') +
        scale_fill_distiller(palette = 'RdBu') +
        theme(aspect.ratio = 1, axis.line = element_line())
PlotPDF('3.03.milo.global_hlhs_vs_control_and_prevad_vs_postvad', 15, 7)
p1 + p2
dev.off()



##~~  Plot 6  Bar - Xenium Sen Cell States across groups   ####
mtx <- as.data.frame(Table(xen.srt$Cell_state, xen.srt$Group2))
total <- table(xen.srt$Group2)
mtx$Freq <- mtx$Freq*100/total[mtx$Var2]
p <- ggplot(mtx[mtx$Var1 %in% c('CM3', 'FB5', 'EC4', 'MP3', 'PC2'), ]) +
        geom_bar(aes(fill = Var1, y = Freq, x = Var2), position = "stack", stat = "identity") +
        scale_fill_manual(values = Color_cell_state_xen) +
        facet_wrap(~Var1, nrow = 1) +
        labs(x = '', y = 'Percentage of Total Cells (%)') +
        theme_Publication(aspect.ratio = 2)
p
PlotPDF('3.06.bar.xenium_sen_cell_state_composition', 6, 2)
p
dev.off()



##~~  Plot 5  UMAP - Global Xenium Cell States   ####
names(Color_cell_state) <- levels(full_clean.srt$Cell_state)
names(Color_cell_state) <- revalue(names(Color_cell_state), replace = c('TC' = 'Lym'))
Color_cell_state_xen <- Color_cell_state[levels(xen.srt$Cell_state)]
p <- DimPlot2(xen.srt, group.by = 'Cell_state', pt.size = 0.1, cols = Color_cell_state_xen,
              raster = F, reduction = 'clean_umap') +
        theme(panel.border = element_rect(fill = NA))
p
PlotPDF('3.05.umap.global_xenium_cell_states', 12, 10)
p
dev.off()



####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
####  Extended Data Figure 3   ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
##~~  Plot 7  Heatmap - Global cell state markers  ####
ngene <- 20
pct <- 0.2
lfc <- 1
p1 <- MarkerHeatmap(cm.srt,
                      marker.df = state_mk.list$CM[state_mk.list$CM$pct.1 > pct &
                                                           state_mk.list$CM$avg_log2FC > lfc, ],
                      disp.min =  0.5, disp.max = 2.5, top = ngene, raster = T, group.cols = mycol_10)
p2 <- MarkerHeatmap(fb.srt,
                      marker.df = state_mk.list$FB[state_mk.list$FB$pct.1 > pct &
                                                           state_mk.list$FB$avg_log2FC > lfc, ],
                      disp.min = 0.5, disp.max = 2.5, top = ngene, raster = T, group.cols = mycol_10)
p3 <- MarkerHeatmap(ec.srt,
                      marker.df = state_mk.list$EC[state_mk.list$EC$pct.1 > pct &
                                                           state_mk.list$EC$avg_log2FC > lfc, ],
                      disp.min = 0.5, disp.max = 2.5, top = ngene, raster = T, group.cols = mycol_10)
p4 <- MarkerHeatmap(mu.srt,
                      marker.df = state_mk.list$Mural[state_mk.list$Mural$pct.1 > pct &
                                                              state_mk.list$Mural$avg_log2FC > lfc, ],
                      disp.min = 0.5, disp.max = 2.5, top = ngene, raster = T, group.cols = mycol_10)
p5 <- MarkerHeatmap(mye.srt,
                      marker.df = state_mk.list$Mye[state_mk.list$Mye$pct.1 > pct &
                                                            state_mk.list$Mye$avg_log2FC > lfc, ],
                      disp.min = 0.5, disp.max = 2.5, top = ngene, raster = T, group.cols = mycol_10)
p6 <- MarkerHeatmap(lym.srt,
                      marker.df = state_mk.list$Lym[state_mk.list$Lym$pct.1 > pct &
                                                            state_mk.list$Lym$avg_log2FC > lfc, ],
                      disp.min = 0.5, disp.max = 2.5, top = ngene, raster = T, group.cols = mycol_10)

p <- wrap_plots(p1, p2, p3, p4, p5, p6, ncol = 2) &
        theme(aspect.ratio = 1) &
        theme(axis.text.y = element_blank()) &
        scale_fill_viridis() &
        NoLegend()
p
PlotPDF('3.07.heat.global_cell_state_marker', 20, 30)
p
dev.off()



##~~  Plot 8  Milo beeswarm - Global HLHS vs Control, PreVAD vs PostVAD   ####
milo <- readRDS('analysis/PART21.global_hlhs_vs_control_milo_obj_result_list.rds')
da_results <- annotateNhoods(milo[[1]], milo[[2]], coldata_col = "Cell_state")
da_results$logFC2 <- da_results$logFC
da_results$logFC2[da_results$SpatialFDR > 0.1] <- NA
da_results$Cell_state <- factor(da_results$Cell_state, levels = levels(full_clean.srt$Cell_state))
p1 <- ggplot(da_results) +
        geom_quasirandom(aes(x = Cell_state, y = logFC, color = logFC2), size = 0.1)
milo <- readRDS('analysis/PART21.global_prevad_vs_postvad_milo_obj_result_list.rds')
da_results <- annotateNhoods(milo[[1]], milo[[2]], coldata_col = "Cell_state")
da_results$logFC2 <- da_results$logFC
da_results$logFC2[da_results$SpatialFDR > 0.1] <- NA
da_results$Cell_state <- factor(da_results$Cell_state, levels = levels(full_clean.srt$Cell_state))
p2 <- ggplot(da_results) +
        geom_quasirandom(aes(x = Cell_state, y = logFC, color = logFC2), size = 0.1)
p <- p1/p2 &
        coord_flip() &
        scale_color_distiller(palette = 'RdBu') &
        scale_x_discrete(limit = rev) &
        theme_classic() &
        theme(aspect.ratio = 3, axis.line = element_line())
p
PlotPDF('3.08.milo_beeswarm.global_hlhs_vs_control', 7, 10)
p
dev.off()


##~~  Plot 1  Bar - Sen vs Non-Sen cell composition in sn and xenium  - NCVR Revision V1   ####
tmp.srt <- full_clean.srt
tmp.srt$Cell_state_sen <- 'Non_Senescent'
tmp.srt$Cell_state_sen[tmp.srt$Cell_state %in% c('CM3', 'MP3', 'PC2', 'EC4', 'FB5')] <- 'Senescent'
p1 <- CompositionBarPlot(tmp.srt, group.by = 'group2', stack.by = 'Cell_state_sen') +
        scale_fill_manual(values = c('grey90', mycol_10[2])) +
        labs(title = 'snRNA-seq')

tmp.srt <- xen.srt
tmp.srt$Cell_state_sen <- 'Non_Senescent'
tmp.srt$Cell_state_sen[tmp.srt$Cell_state %in% c('CM3', 'MP3', 'PC2', 'EC4', 'FB5')] <- 'Senescent'
p2 <- CompositionBarPlot(tmp.srt, group.by = 'Group2', stack.by = 'Cell_state_sen') +
        scale_fill_manual(values = c('grey90', mycol_10[2])) +
        labs(title = 'Xenium ST')
p1 + p2
PlotPDF('1.01.bar.composition_of_senescent_cells_sn_xenium', 10, 10)
p1 + p2
dev.off()



##~~  Plot 1  Radar - snRNA per patient senescence score  - NCVR Revision V1 ####
srt_list <- list(cm.srt, fb.srt, ec.srt, mu.srt, mye.srt)
mtx_list <- list()
for(x in 1:5){
        data <- Table(srt_list[[x]]$group1, srt_list[[x]]$Cell_state)
        total <- Table(srt_list[[x]]$group1)
        norm <- as.matrix(data)
        for(i in 1:ncol(data)){norm[, i] <- norm[, i]/total[rownames(data)]}
        mtx_list[[x]] <- norm
}
final <- cbind(mtx_list[[1]],
               mtx_list[[2]],
               mtx_list[[3]],
               mtx_list[[4]],
               mtx_list[[5]])
sen_final <- final[, c('CM3', 'FB5', 'EC4', 'PC2', 'MP3')]
df_long <- melt(sen_final)
p <- ggplot(df_long, aes(x = Var2, y = value*100)) +
        geom_bar(stat = "identity", fill = "steelblue") +
        coord_polar(start = 0) +
        theme_classic() +
        theme_Publication(aspect.ratio = 1, grid = T) +
        facet_wrap(~Var1)
p
PlotPDF('2.01.radar.sen_niche_cells', 10, 10)
p
dev.off()



##~~  Plot 2  Radar - Xenium per patient senescence score  - NCVR Revision V1  ####
srt_list <- SplitObject(xen.srt, split.by = 'Cell_type')
srt_list <- srt_list[c('CM', 'FB', 'EC', 'PC', 'Mye')]
mtx_list <- list()
for(x in 1:5){
        data <- Table(srt_list[[x]]$Group1, srt_list[[x]]$Cell_state)
        total <- Table(srt_list[[x]]$Group1)
        norm <- as.matrix(data)
        for(i in 1:ncol(data)){norm[, i] <- norm[, i]/total[rownames(data)]}
        mtx_list[[x]] <- norm
}
final <- cbind(mtx_list[[1]],
               mtx_list[[2]],
               mtx_list[[3]],
               mtx_list[[4]],
               mtx_list[[5]])
sen_final <- final[, c('CM3', 'FB5', 'EC4', 'PC2', 'MP3')]
df_long <- melt(sen_final)
p <- ggplot(df_long, aes(x = Var2, y = value*100)) +
        geom_bar(stat = "identity", fill = "steelblue") +
        coord_polar(start = 0) +
        theme_classic() +
        theme_Publication(aspect.ratio = 1, grid = T) +
        facet_wrap(~Var1, nrow = 3)
p
PlotPDF('2.02.radar.xenium_sen_niche_cells', 6, 6)
p
dev.off()



####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
####  Figure 4    ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
##~~  Plot 1  Box - Cell states SASP expression  ####
include <- c('CM1', 'CM2', 'CM3',
             'FB1', 'FB2', 'FB3', 'FB4', 'FB5',
             'EC1', 'EC2', 'EC3', 'EC4',
             'Mono', 'MP1', 'MP2', 'MP3',
             'PC1', 'PC2')
tmp.srt <- DropMetaLevels(full_clean.srt[, full_clean.srt$Cell_state %in% include])
tmp.srt <- AddModuleScore2(tmp.srt, features = list(sen_combine), names = 'Sen', return_z = T)
x <- c('grey80', 'pink')[c(1, 2, 2,
                           1, 1, 1, 1, 2,
                           1, 1, 1, 2,
                           1, 2,
                           1, 1, 1, 2)]
p <- BoxPlot2(tmp.srt, feature = 'Sen', group.by = 'Cell_state', cols = x) +
        labs(x = '') +
        RotatedAxis() +
        theme(aspect.ratio = 1)
p
PlotPDF('4.01.box.fb_ec_mye_pc_sasp', 9, 9)
p
dev.off()



##~~  Plot 8  Bar plot - Sen Niche shared modules  ####
x <- c(state_mk.list$FB$gene[state_mk.list$FB$cluster == 'FB5'],
       state_mk.list$EC$gene[state_mk.list$EC$cluster == 'EC4'],
       state_mk.list$Mye$gene[state_mk.list$Mye$cluster == 'MP3'],
       state_mk.list$Mural$gene[state_mk.list$Mural$cluster == 'PC2'],
       state_mk.list$CM$gene[state_mk.list$CM$cluster == 'CM3'])
sort(Table(x), decreasing = T)
shared_mod <- names(Table(x))[as.vector(Table(x)) == 5]
shared_mod_enr <- ModuleEnrichment(list(shared_mod), human_or_mouse = 'human')
data1 <- shared_mod_enr$Reactome$Reactome_
data1 <- data1[c(data1$Description %in% c(
        'Interleukin-6 signaling',
        'Transcriptional Regulation by TP53'
)), c('Description', 'p.adjust')]
data2 <- shared_mod_enr$KEGG$KEGG_
data2 <- data2[c(data2$Description %in% c(
        'JAK-STAT signaling pathway',
        'MAPK signaling pathway',
        'Cellular senescence'
)), c('Description', 'p.adjust') ]
data <- rbind(data1, data2)
data$p.adjust <- -log10(data$p.adjust)
data <- data[order(data$p.adjust, decreasing = T), ]
data$Description <- factor(data$Description, levels = data$Description)
p <- ggplot(data, aes(x = p.adjust, y = Description))+
        geom_bar(stat = 'identity', position = "dodge", fill = mycol_10[1], alpha = 0.4) +
        scale_y_discrete(limit = rev) +
        labs(x = '-Log10 p value', y = 'Sen Niche Signatures') +
        theme_Publication(aspect.ratio = 0.5)
p
PlotPDF('5.08.bar.sen_niche_shared_sig_enrichment', 4, 2)
p
dev.off()



##~~  Plot 3  Heatmap - HLHS cell state pairwise colocalization  ####
mk1 <- readRDS('analysis/PART28.cell_state_marker.srt_mk.rds')
mk1 <- rbind(mk1$CM, mk1$FB, mk1$EC, mk1$ECm, mk1$Mural, mk1$Mye, mk1$Lym)
unique_gene <- names(Table(mk1$gene))[Table(mk1$gene) == 1]
mk_cs <- mk1[mk1$avg_log2FC > 0.25 & mk1$p_val_adj < 1e-3 & mk1$gene %in% unique_gene,]
Table(mk_cs$cluster)
st.srt <- AddModuleScore2(st.srt, features = split(mk_cs$gene, mk_cs$cluster),
                          names = paste0('CS_Score_', levels(mk_cs$cluster)), return_z = T)
tmp.srt <- st.srt[, st.srt$group2 == 'HLHS']
p_cutoff <- 0.001
ct <- levels(mk_cs$cluster)
mtx <- matrix(NA, L(ct), L(ct))
rownames(mtx) <- ct
colnames(mtx) <- ct
for(i in 1:L(ct)){
        message(i)
        for(j in 1:L(ct)){
                if(i < j){
                        x <- GetColocalProb(tmp.srt, meta_features = paste0('CS_Score_', c(ct[i], ct[j])))
                        y <- x > -log10(p_cutoff)
                        mtx[i, j] <- sum(y)/L(y)
                }
        }
}
data <- reshape2::melt(mtx)
p <- ggplot(data) +
        geom_tile(aes(x = Var1, y = Var2, fill = value)) +
        scale_fill_viridis_c(na.value = 0, values = c(0, 0.8, 1)) +
        theme_classic() +
        theme(aspect.ratio = 1) +
        RotatedAxis()
p
PlotPDF('03.1.heat.pairwise_cell_state_colocalization', 10, 10)
p
dev.off()



##~~  Plot 4  ST Feature - Colocalization p value   ####
p <- FeaturePlotST_Dark(st.srt,
                        features = 'Colocal_pval',
                        minvals = 1,
                        maxvals = 3,
                        pt.sizes = spot_sizes*1.8,
                        ncol = 6) &
        scale_colour_viridis_c(option = 'magma', begin = 0.1, limits = c(1, 3), values = c(0, 0.75, 1))
p
PlotPDF('4.04.st_feat.colocalization_pval', 22, 3)
p
dev.off()



##~~  Plot 5  ST - Niche and neighbor   ####
p <- DimPlotST(st.srt, group.by = 'dNiche', pt.sizes = spot_sizes*1.8, ncol = 6,
               cols = c(Color_dniche, 'grey80'), legend = 6)
p
PlotPDF('4.05.st_dim.niche_neighbor_binary', 18, 3)
p
dev.off()



##~~  Plot 6  Violin Colocalization p value across disease group  ####
p <- VlnPlot2(st.srt, feature = 'Colocal_pval', group.by = 'group1', pt.size = 0.1,
               cols = c('grey', mycol_10), adjust = 5) +
        labs(caption = paste0('Wilcox P < ', pval)) +
        geom_hline(yintercept = -log10(0.001), color = 'red', linetype = 'dashed')
pval <- wilcox.test(split(st.srt$Colocal_pval, st.srt$group1)[[2]],
                    split(st.srt$Colocal_pval, st.srt$group1)[[3]])$p.value
PlotPDF('4.06.vln.colocal_pval_across_groups', 4, 4)
p
dev.off()



##~~  Plot 7  ST Neighbor Trend - Niche Neighbor Senescence Score and SASP Score   - NCVR Revision V1  ####
tmp.srt <- AddModuleScore2(st.srt, features = sen_combine, names = 'Sen', return_z = T, assay = 'SCT')
data <- tmp.srt@meta.data |> group_by(dNiche, group2) |> mutate(Mean_SASP = mean(Sen))
data <- U(data[, c('group1', 'group2', 'dNiche', 'Mean_SASP')])
p <- ggplot(data[data$group2 != 'Control', ]) +
        geom_point(aes(x = dNiche, y = Mean_SASP, color = dNiche)) +
        geom_line(aes(x = dNiche, y = Mean_SASP, group = group2)) +
        # scale_y_continuous(limits = c(0, 1.5)) +
        scale_color_manual(values = Color_dniche) +
        labs(x='', y = 'Mean Expression Z-score', title = 'SASP') +
        theme_classic() +
        theme(aspect.ratio = 2) +
        RotatedAxis() +
        facet_wrap(~group2)
p
PlotPDF('2.07.visium_trend.senenscens_in_dniche_neighbors', 5, 5)
p
dev.off()




##~~  Plot 6  Feature Xenium - Sen sig expression on Xenium slides   ####
sub.srt <- DropMetaLevels(xen.srt[, xen.srt$Group1 == 'HLHS_15'])
sub.srt <- AddModuleScore2(sub.srt, features = sen_combine,
                           names = 'Sen', return_z = T, assay = 'Xenium',
                           nbin = 20, ctrl = 20)
sub.srt@reductions$spatial <- sub.srt@reductions$clean_umap
sub.srt@reductions$spatial@cell.embeddings[,1] <- sub.srt$Coord_x
sub.srt@reductions$spatial@cell.embeddings[,2] <- sub.srt$Coord_y
colnames(sub.srt@reductions$spatial@cell.embeddings) <- paste0('spatial_', 1:2)
sub.srt@reductions$spatial@key <- 'spatial_'
p <- FeaturePlot2(sub.srt, features = 'Sen', min.cutoff = 'q5', max.cutoff = 'q95', pt.size = 1) +
        scale_color_viridis_c(option = 'B', begin = 0.1) +
        theme(aspect.ratio = max(xen.srt$Coord_y)/ max(xen.srt$Coord_x),
              plot.background = element_rect(fill = "black"),
              panel.background = element_rect(fill = "black"),
              legend.background = element_rect(fill = "black"),
              legend.box.background = element_rect(fill = "black", size = 0),
              legend.key = element_rect(fill = "black", size = 0),
              strip.background = element_rect(fill = "grey50", colour = NA),
              axis.line.x = element_line(colour = "white"),
              axis.line.y = element_line(colour = "white"),
              panel.grid = element_blank(),
              panel.grid.minor = element_blank(),
              axis.title = element_blank(),
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              text = element_text(colour = "white")) +
        NoLegend()
p
PlotPDF('2.06.feature.sen_combine_expression', 6, 10)
p
dev.off()




##~~  Plot 13  Line trend - EC4 surrounding SASP/SenSig level   ####
samples <- c('HLHS_16', 'HLHS_12', 'HLHS_15', 'VAD_04', 'VAD_05')
tmp.srt <- AddModuleScore2(xen.srt, features = sen_combine,
                           names = 'Sen', return_z = T, assay = 'Xenium',
                           nbin = 20, ctrl = 20)
max_dist <- 150
bin <- 5
mean_val.mtx <- matrix(NA, L(samples), max_dist/bin)
for(x in 1:L(samples)){
        message(x)
        sub.srt <- DropMetaLevels(tmp.srt[, tmp.srt$Group1 == samples[x]])
        sub.srt$Cell_state[sub.srt$Cell_state %in% c('FB5', 'MP3', 'PC2', 'EC4')] <- 'SenNiche' ## consider all niche cells
        mtx <- matrix(0, Table(sub.srt$Cell_state)['SenNiche'], max_dist/bin,
                      dimnames = list(Cells(sub.srt)[sub.srt$Cell_state == 'SenNiche'], 1:(max_dist/bin)))
        coords_que <- as.matrix(sub.srt@meta.data[, c("Coord_x", "Coord_y")])
        for(j in 1:Table(sub.srt$Cell_state)['SenNiche']) {
                coords_viewpoint <- as.matrix(sub.srt@meta.data[Cells(sub.srt)[sub.srt$Cell_state == 'SenNiche'][j],
                                                                c("Coord_x", "Coord_y")])
                distance_matrix <- proxy::dist(coords_que, coords_viewpoint, method = "Euclidean")
                nn_list <- list()
                mean_sasp <- c()
                for(i in 1:(max_dist/bin)){
                        nn_list[[i]] <- sub.srt@meta.data[
                                rownames(distance_matrix)[as.vector(distance_matrix) <= bin*i &
                                                                  as.vector(distance_matrix) > bin*(i-1)],
                                'Sen']
                        mean_sasp[i] <- mean(nn_list[[i]])
                }
                mtx[j,] <- mean_sasp
        }
        mean_val.mtx[x, ] <- colMedians(mtx, na.rm = T)
}
rownames(mean_val.mtx) <- samples
mean_val.mtx2 <- rbind(colMeans2(mean_val.mtx[1:3, ]), colMeans2(mean_val.mtx[4:5, ]))
df <- data.frame(t(mean_val.mtx2))
colnames(df) <- c('HLHS', 'VAD')
df$Distance <- paste(1:(max_dist/bin)*bin, 'µm')
df <- reshape2::melt(df)
df$Distance <- as.numeric(str_split(df$Distance, pattern = ' ', simplify = T)[, 1])
p <- ggplot(df, aes(x = Distance, y = value)) +
        geom_point() +
        geom_smooth(method = 'loess', span = 0.8) +
        labs(x = 'Distance from a Sen. niche cell (µm)', y = 'Sen Sig Mean Expr.') +
        theme_Publication(aspect.ratio = 2) +
        facet_wrap(~variable, nrow = 1)
p
PlotPDF('2.09.spatial_trend_line.xenium_sen_combine_decline_from_sen_niche_cell', 8, 4)
p
dev.off()



##~~  Plot 9  Heatmap - Sen Cell State Colocalization   ####
mtx <- read.csv(paste0('~/Documents/Bioinformatics/project/2022_hlhs_dturaga/meta/human_v2/',
                'PART69.xenium_analysis_squidpy/HLHS15_nhood_enrichment_zscore.csv'),
                header = T, row.names = 1)
for(i in 1:ncol(mtx)){for (j in 1:nrow(mtx)) {if(i<=j){mtx[j, i] <- NA}}}
mtx$Cell_state <- rownames(mtx)
df <- melt(mtx)
df$value[df$value < 1] <- 1
df$value[df$value > 30] <- 30
df$Cell_state <- factor(df$Cell_state, levels = rownames(mtx))
df$variable <- factor(df$variable, levels = rownames(mtx))
p <- ggplot(df) +
        geom_tile(aes(x = Cell_state, y = variable, fill = value)) +
        scale_y_discrete(limit = rev) +
        scale_fill_viridis_c(values = c(0, 0.25, 1), na.value = 'white') +
        labs(x = '', y = '', fill = 'Z-score', title = 'Neighborhood enrichment') +
        theme_Publication(aspect.ratio = 1) +
        RotatedAxis()
p
PlotPDF('4.09.heatmap.sen_cell_state_neighborhood_enrichment', 4, 3)
p
dev.off()




##~~  Plot 10  Heatmap + UMAP - H15 spatial cluster cell state enrichment   ####
H15 <- readRDS('analysis/PART68.xenium_niche_clustering_HLHS_15')
obj <- DropMetaLevels(xen.srt[, xen.srt$Group1 == 'HLHS_15'])
obj@reductions$spatial <- obj@reductions$clean_umap
coor <- xen_img.srt@images$HLHS_15$centroids@coords[, 1:2]
rownames(coor) <- Cells(xen_img.srt)[xen_img.srt$Group1 == 'HLHS_15']
obj@reductions$spatial@cell.embeddings <- coor[Cells(obj), 1:2]
colnames(obj@reductions$spatial@cell.embeddings) <- c('spatial_1', 'spatial_2')
obj@reductions$spatial@key <- 'spatial_'
obj@meta.data <- H15
p <- MetaMatchingHeatmap(obj,  ref_var = 'Cell_state', que_var = 'niche_label', percentage = T,
                         center = F, ref.disp.min = 0.02) +
        scale_fill_viridis_c(option = 'D')
data <- p$data
data$Value[data$Value > 0.25] <- 0.25

obj$niche_bin <- 'Niche1'
obj$niche_bin[obj$niche_label%in% c(13,15,12,4)] <- 'Niche2'
obj$niche_bin[obj$niche_label%in% c(9,10,0,11,5)] <- 'Sen_niche'

p1 <- ggplot(data)+
        geom_tile(aes(y = Reference, x = Query, fill = Value)) +
        scale_fill_viridis_c(limits = c(0, 0.25), option = 'D') +
        scale_y_discrete(limits = rev) +
        theme(aspect.ratio = LU(data$Reference)/LU(data$Query))
p2 <- DimPlot2(obj, reduction = 'spatial', group.by = 'niche_bin',
               cols = c(mycol_10[1], mycol_10[1], mycol_10[2]), pt.size = 0.1) +
        theme(aspect.ratio = (max(obj$Coord_x)-min(obj$Coord_x))/(max(obj$Coord_y)-min(obj$Coord_y)))
p1 + p2
PlotPDF('2.10.spatial_niche_xenium_HLHS15', 16, 8)
p1 + p2
dev.off()



##~~  Plot 10  Line trend - EC4 Colocalization with Sen Cell States   ####
mtx <- read.csv('PART68.xenium_analysis_squidpy/HLHS_EC4_co_occurrence.csv', header = T, row.names = 1)
mtx$Cell_state <- rownames(mtx)
df <- melt(mtx)
df$variable <- as.numeric(str_remove(df$variable, pattern = 'X'))
df$Cell_state <- factor(df$Cell_state, levels = levels(xen.srt$Cell_state))
color <- c('grey', 'red')[as.numeric(levels(xen.srt$Cell_state) %in% c('EC4', 'FB5', 'MP3', 'PC2')) + 1]
df <- df[df$variable <= 150, ]
df$thick <- as.numeric(df$Cell_state %in% c('EC4', 'FB5', 'MP3', 'PC2'))+1
p <- ggplot() +
        geom_smooth(data = df, aes(x = variable, y = value, color = Cell_state, linewidth = thick),
                    se = F, method = 'loess', span = 1) +
        ggrepel::geom_label_repel(data = df[df$variable == min(df$variable), ],
                                  aes(label = Cell_state, x = variable, y = value) )+
        scale_color_manual(values = Color_cell_state) +
        theme_Publication(aspect.ratio = 1)
PlotPDF('4.10.line.ec4_co_occur_cell_state', 8, 7)
p
dev.off()



####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
####  Extended Data Figure 4   ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
##~~  Plot 15  Scatter  - Senescence Marker Expression   ####
include <- c('CM1', 'CM2', 'CM3',
             'FB1', 'FB2', 'FB3', 'FB4', 'FB5',
             'EC1', 'EC2', 'EC3', 'EC4',
             'PC1', 'PC2',
             'Mono', 'MP1', 'MP2', 'MP3')
ct <- c('CM', 'EC', 'FB', 'Mye', 'PC')
sen_cs <- c('CM2', 'CM3', 'FB5', 'EC4', 'PC2', 'MP3')
cols <- Color_cell_state[levels(full_clean.srt$Cell_state) %in% include]

meta <- data.frame()
for(i in 1:L(ct)){
        tmp.srt <- DropMetaLevels(full_clean.srt[, full_clean.srt$Cell_state %in% include &
                                                         full_clean.srt$Cell_type %in% ct[i]])
        tmp.srt <- AddModuleScore2(tmp.srt,
                                   features = c(gs5,
                                                'CDKN1A', 'CDKN2A', 'ATF6', 'BHLHE40', 'TNFRSF10D', 'ARF1',
                                                'SERPINE1', 'HIF1A', 'MMP10', 'MAPK14', 'LMNB1', 'GATA4'),
                                   names = c('Sen',
                                             paste0(c('CDKN1A', 'CDKN2A', 'ATF6', 'BHLHE40', 'TNFRSF10D', 'ARF1',
                                                      'SERPINE1', 'HIF1A', 'MMP10', 'MAPK14', 'LMNB1', 'GATA4'), '_Score')),
                                   return_z = T)
        meta <- rbind(meta, tmp.srt@meta.data)
}
data <- meta[, c('Cell_state', 'Cell_type', 'Sen',
                 paste0(c('CDKN1A', 'CDKN2A', 'ATF6', 'BHLHE40', 'TNFRSF10D', 'ARF1',
                          'SERPINE1', 'HIF1A', 'MMP10', 'MAPK14', 'LMNB1', 'GATA4'), '_Score'))] |>
        group_by(Cell_state) |>
        mutate(Mean_Sen = mean(Sen),
               Mean_CDKN1A= mean(CDKN1A_Score),
               Mean_CDKN2A = mean(CDKN2A_Score),
               Mean_BHLHE40 = mean(BHLHE40_Score),
               Mean_TNFRSF10D = mean(TNFRSF10D_Score),
               Mean_ARF1 = mean(ARF1_Score),
               Mean_SERPINE1 = mean(SERPINE1_Score),
               Mean_HIF1A = mean(HIF1A_Score),
               Mean_LMNB1 = mean(LMNB1_Score),
               Mean_MMP10 = mean(MMP10_Score),
        )
data <- U(data[, c(grep('Mean_', colnames(data), value = T), 'Cell_state', 'Cell_type')])
data$Cell_state <- factor(data$Cell_state, levels = include)
data <- reshape2::melt(data)
data$SASP <- data$value[data$variable == 'Mean_Sen']
data <- data[data$variable != 'Mean_Sen', ]
data$highlight <- ifelse(data$Cell_state %in% sen_cs, yes = 'Sen', no = 'Non-Sen')
data <- data |> group_by(variable) |> mutate(value2 =  Range01(value))
p <- ggplot(data) +
        geom_label(aes(x = SASP, y = value2, label = Cell_state, color = highlight)) +
        scale_color_manual(values = c('black', 'red')) +
        theme_classic() +
        theme(aspect.ratio = 1) +
        facet_wrap(~variable, scales = 'free')
p
PlotPDF('4.15.scatter.sen_marker_expression_in_cell_states', 15, 15)
p
dev.off()


##~~  Plot 10  Heatmap + UMAP - V4 spatial cluster cell state enrichment   ####
V4 <- readRDS('analysis/PART68.xenium_niche_clustering_VAD_04')
obj <- DropMetaLevels(xen.srt[, xen.srt$Group1 == 'VAD_04'])
obj@reductions$spatial <- obj@reductions$clean_umap
coor <- xen_img.srt@images$VAD_04$centroids@coords[, 1:2]
rownames(coor) <- Cells(xen_img.srt)[xen_img.srt$Group1 == 'VAD_04']
obj@reductions$spatial@cell.embeddings <- coor[Cells(obj), 1:2]
colnames(obj@reductions$spatial@cell.embeddings) <- c('spatial_1', 'spatial_2')
obj@reductions$spatial@key <- 'spatial_'
obj@meta.data <- V4
obj <- DropMetaLevels(obj[, !obj$Cell_state %in% c('AC', 'NC', 'LEC')])
p <- MetaMatchingHeatmap(obj,  ref_var = 'Cell_state', que_var = 'niche_label', percentage = T) +
        scale_fill_viridis_c( option = 'D')
p
data <- p$data
p1 <- ggplot(data)+
        geom_tile(aes(y = Reference, x = Query, fill = Value)) +
        scale_fill_viridis_c(limits = c(0, .25), option = 'D') +
        scale_y_discrete(limits = rev) +
        theme(aspect.ratio = LU(data$Reference)/LU(data$Query))
p1
PlotPDF('2.11.spatial_niche_xenium_VAD5', 16, 8)
p1
dev.off()
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
####  Figure 5    ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
##~~  Plot 1  Heatmap - Hif1a in all cell states  ####
cs_exclude <- c('Prol MP', 'NK', 'Prol TC', 'AC', 'NC', 'LEC', 'EndoC', 'MP4', 'DC', 'Mast') ## small clusters
tmp.srt <- AddModuleScore2(full_clean.srt[, ! full_clean.srt$Cell_state %in% cs_exclude],
                           features = c(hif1a, 'HIF1A', 'IL6'),
                           names = c('HIF1A_Targets', 'HIF1A_Expr', 'IL6_Expr'))
tmp.srt <- tmp.srt[, tmp.srt$group2 != 'Control']
p <- MarkerHeatmapMean(tmp.srt, features = c('IL6_Expr', 'HIF1A_Targets', 'HIF1A_Expr'), group.by = 'Cell_state',
                         col.max = 1.8, col.min = 0) +
        scale_fill_viridis_c() +
        scale_x_discrete(limits = rev) +
        RotatedAxis() +
        coord_flip() +
        theme(aspect.ratio = 21/3)
p
PlotPDF('5.01.heat.hif1a_il6_cell_states', 3, 10)
p
dev.off()



##~~  Plot 2 UMAP - EC Cell States   ####
p <- DimPlot2(ec.srt, group.by = 'Cell_state', reduction = 'sub_umap', cols = mycol_10, pt.size = 0.1) +
        theme(aspect.ratio = 1)
p
PlotPDF('5.02.umap.ec_ell_state', 8, 8)
p
dev.off()



##~~  Plot 3 Box - EC4 Abundance across Disease Groups   ####
meta <- U(ec.srt@meta.data[, c('group1', 'group2')])
meta <- meta[order(meta$group1),]
mtx <- Table(ec.srt$group1, ec.srt$Cell_state) + 1
mtx <- as.matrix(mtx[, colSums(mtx)>0])
mtx <- mtx/as.vector(table(ec.srt$group1))
df <- data.frame(Cell = mtx[, 'EC4'], group2 = as.factor(meta$group2), value = as.numeric(meta$group2))
wilcox.test(split(df$Cell, df$group2)[[1]], split(df$Cell, df$group2)[[2]])$p.value < 0.05
wilcox.test(split(df$Cell, df$group2)[[3]], split(df$Cell, df$group2)[[2]])$p.value < 0.05
p <- ggplot(df) +
        geom_boxplot(aes(y = Cell, x = group2, fill = group2), outlier.shape = NA) +
        geom_beeswarm(aes(y = Cell, x = group2, color = group2), size = 2, cex = 3) +
        labs(title = 'EC4', y = 'Fraction of Total ECs', x = 'group2') +
        scale_fill_manual(values = Color_palliation[c(1, 2, 5)]) &
        scale_color_manual(values = Color_palliation[c(1, 2, 5)]) &
        theme_classic() &
        theme_Publication(aspect.ratio = 1) &
        NoLegend()
p
PlotPDF('5.03.box.ec4_across_disease_group', 6, 6)
p
dev.off()



##~~  Plot 4 Bar - EC4 Abundance across pairs    ####
tmp.srt <- paired.srt[, paired.srt$Cell_state %in% c(paste0('EC', 1:4))]
meta <- U(tmp.srt@meta.data[, c('palliation', 'group2', "vad_pair")])
meta <- meta[order(meta$vad_pair),]
mtx <- Table(tmp.srt$vad_pair, tmp.srt$Cell_state)
mtx <- as.matrix(mtx[, colMaxs(mtx ) > 0])
total <- table(paired.srt$vad_pair)
mtx <- mtx/as.vector(total[rownames(mtx)])
df <- data.frame(Frac = mtx[, 'EC4'],
                 group2 = as.factor(meta$group2),
                 value = as.numeric(as.factor(meta$group2)),
                 pair = str_split(meta$vad_pair, ' ', simplify = T)[,2])
data2 <- data.frame(Group = c('1', '2', '3'),
                    Fold_reduction = log2(df[df$group2 == 'VAD', ]$Frac/df[df$group2 == 'HLHS', ]$Frac))
data2$Cell_type <- 'EC'
p <- ggplot(data2) +
        geom_bar(aes(x = Group, y = Fold_reduction, fill = Group), stat = 'identity') +
        scale_fill_manual(values = Color_pair) +
        labs(title = 'EC4', y = 'Log2 FC (Pre vs Post-VAD)', x = 'Pre/Post-VAD Pair') +
        theme_classic() +
        theme(aspect.ratio = 2) +
        NoLegend()
p
PlotPDF('5.04.bar.ec4_abundance_in_pre_post_pair', 4, 6)
p
dev.off()



##~~  Plot 5 Heatmap - EC4 hypoxia and senescence    ####
tmp.srt <- DropMetaLevels(ec.srt[, ec.srt$Cell_state %in% paste0('EC', 1:4)])
tmp.srt$tmp <- 'Other ECs'
tmp.srt$tmp[tmp.srt$Cell_state == 'EC4'] <- 'EC4'
tmp.srt$tmp <- paste(tmp.srt$tmp, tmp.srt$group2)
p <- MarkerHeatmapMean(tmp.srt, features = c('HIF1A', 'HIF1_Score', 'VEGFA_Score',
                                             'CDKN1A', 'IL6', 'CSF3'), group.by = 'tmp',
                       col.max = 2, col.min = 0) +
        scale_fill_distiller(palette = 'RdBu', values = c(0, .5, 1)) +
        coord_flip() +
        RotatedAxis() +
        theme(aspect.ratio = 6/6)
p
PlotPDF('5.05.heatmap.ec_latent_gene_post_vad', 6, 6)
p
dev.off()



##~~  Plot 6  Vortex - HLHS EC4 secrets IL6  ####
PlotPDF('5.06.vortex.ec4_il6_cellchat', 5, 5)
netVisual_individual(hlhs.cs.cch, 
                     remove.isolate = T,
                     pairLR.use = 'IL6_IL6R_IL6ST',
                     signaling = 'IL6',
                     color.use = Color_cell_state,
                     layout = "circle")
dev.off()




##~~  Plot 7  ST  - IL6 LR in ST  ####
tmp.srt <- AddModuleScore2(st.srt, features = c(list('IL6'), list(c('IL6R', 'IL6ST'))),
                           names = c('IL6_score', 'IL6R_score'), return_z = T)
tmp.srt$IL6_score_smooth <- SmoothFeature(tmp.srt, feature = 'IL6_score')
tmp.srt$IL6R_score_smooth <- SmoothFeature(tmp.srt, feature = 'IL6R_score')
tmp.srt$IL6_Coloc <- GetColocalProb(tmp.srt, meta_features = c('IL6_score_smooth', 'IL6R_score_smooth'))
minvals <- 0.5
maxvals <- 2.5
feat <- c('IL6_score_smooth', 'IL6R_score_smooth')
p1 <- FeaturePlotST_Dark(tmp.srt, features = feat ,
                            minvals = rep(minvals, 2),
                            maxvals = rep(maxvals, 2), pt.sizes = spot_sizes*1.8, ncol = 6)
for(i in 1:L(feat)){for(j in 1:6){p1[[6*(i-1)+j]] <- p1[[6*(i-1)+j]] + theme(aspect.ratio = asp_ratio[j])}}
minvals <- 1
maxvals <- 4
feat <- c('IL6_Coloc')
p2 <- FeaturePlotST_Dark(tmp.srt, features = feat ,
                            minvals = minvals,
                            maxvals = maxvals, pt.sizes = spot_sizes*1.8, ncol = 6) &
        scale_colour_viridis_c(option = 'magma', begin = 0.1, limits = c(minvals, maxvals), values = c(0, 0.3, 1))
for(i in 1:L(feat)){for(j in 1:6){p2[[6*(i-1)+j]] <- p2[[6*(i-1)+j]] + theme(aspect.ratio = asp_ratio[j])}}
p <- wrap_plots(p1[[1]], p1[[2]], p1[[3]], p1[[4]], p1[[5]], p1[[6]],
                  p1[[7]], p1[[8]], p1[[9]], p1[[10]], p1[[11]], p1[[12]],
                  p2[[1]], p2[[2]], p2[[3]], p2[[4]], p2[[5]], p2[[6]], ncol = 6)
p
PlotPDF('5.07.st_feat.il6_lr_interaction_probability', 22, 8)
p
dev.off()


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
####  Extended Data Figure 5   ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
##~~  Plot 10  Heatmap EC/Mural Identity Prediction  ####
p <- MappingHeatmap(ec.srt, que_var = 'Cell_state', ref_var = 'scANVI_predict_cellstate',
                      percentage = T, ref.disp.min = 0.1,
                      que_order = levels(droplevels(ec.srt$Cell_state)),
                      ref_order = c('EC5_art', 'EC2_cap', 'EC6_ven', 'EC4_immune', 'EC8_ln', 'EC7_atria'),
                      center = F) +
        scale_fill_distiller(palette = 'RdYlBu', limits = c(0, 1)) +
        theme(aspect.ratio = 1)
p
PlotPDF('5.10.heatmap.ec_identity_mapping', 5, 5)
p
dev.off()



##~~  Plot 11 Feature UMAP - EC4 pseudotime   ####
tmp.srt <- DropMetaLevels(ec.srt[, ec.srt$Cell_state %in% paste0('EC', 2:4) & !is.na(ec.srt$Latent_time)])
tmp.srt <- tmp.srt[, tmp.srt@reductions$sub_main_umap@cell.embeddings[, 1] > 0 &
                           tmp.srt@reductions$sub_main_umap@cell.embeddings[, 2] < 0 ]
p <- FeaturePlot2(tmp.srt, features = 'Latent_time', reduction = 'sub_main_umap', order = F) +
        scale_color_viridis_c() +
        theme(aspect.ratio = 1.5)
p
PlotPDF('5.11.feature_umap.ec4_pseudotime', 5, 10)
p
dev.off()



##~~  Plot 12 Gene Trend - EC latent time module   ####
tmp.srt <- DropMetaLevels(ec.srt[, ec.srt$Cell_state %in% paste0('EC', 2:4) & !is.na(ec.srt$Latent_time)])
tmp.srt <- tmp.srt[, tmp.srt@reductions$sub_main_umap@cell.embeddings[, 1] > 0 &
                           tmp.srt@reductions$sub_main_umap@cell.embeddings[, 2] < 0 ]
tmp.srt <- ScaleData(tmp.srt, features = rownames(tmp.srt), assay = 'CBN')
ext.mtx <- tmp.srt@assays$CBN@scale.data
gl <- data$Gene[data$PCC > 0.05 & data$AdjPval < 1e-20]
L(gl) ##233
ext.mtx <- tmp.srt@assays$CBN@scale.data[as.vector(gl), colnames(tmp.srt)[order(tmp.srt$Latent_time_rank)]]
for(i in 1:nrow(ext.mtx)){ext.mtx[i, ] <- smooth.spline(ext.mtx[i, ], spar = 1.1)$y}
sub_ext.mtx <- ext.mtx[, round(ncol(ext.mtx)/10):(ncol(ext.mtx)-round(ncol(ext.mtx)/10)) ]
par(pty = "s")
PlotPDF('5.12.line_trend.ec4_latent_time_correlated_gene_module', 4, 4.5)
plot(c(1, ncol(sub_ext.mtx)), c(-.3, .4), type = 'l', col = alpha('white', 0))
for(i in 1:nrow(sub_ext.mtx)){
        points(1:ncol(sub_ext.mtx), sub_ext.mtx[i, ], type = 'l', col = alpha('black', 0.2))}
dev.off()



##~~  Plot 12 Trend EC and Mural correlation in Hypoxia  ####
gl <- c('ARNT', 'ASPH', 'COPS5', 'CREB1', 'EDN1', 'EP300', 'EPO', 
        'HIF1A', 'HSP90AA1', 'JUN', 'LDHA', 'NOS3', 'P4HB', 'UBE2A', 'VEGFA', 'VHL')
# Hypoxia-Inducible Factor In The Cardivascular System - Biocarta Pathways
pc.srt <- mu.srt[, mu.srt$Cell_type %in% c('PC')]
ec.score <- split(ec.srt$Hypoxia_Score, ec.srt$group1)
mu.score <- split(pc.srt$Hypoxia_Score, pc.srt$group1)
pal <- U(ec.srt@meta.data[, c('group1', 'palliation')])[order(U(ec.srt@meta.data[, c('group1', 'palliation')])[,1]), 2]
pal <- pal[pal != 'VAD']
data <- data.frame(
        Palliation = pal,
        EC_hyoxia_target = rep(NA, 23),
        MU_hyoxia_target = rep(NA, 23)
)
for(i in 1:23){
        data$EC_hyoxia_target[i] = mean(ec.score[[i]])
        data$MU_hyoxia_target[i] = mean(mu.score[[i]])
}
data <- data[7:23, ]
k <- 16
rsq <- summary(mgcv::gam(data = data, formula = EC_hyoxia_target ~ s(MU_hyoxia_target, k = k)))$r.sq
p <- ggplot(data, aes(y = MU_hyoxia_target, x = EC_hyoxia_target, color = Palliation)) +
        geom_point() +
        geom_smooth(method = "gam", formula = y ~ s(x, k = k), fullrange = T, se = T, color = 'grey30') +
        labs(title = 'Hypoxia Response', y = 'Mural Cell Hypoxia Response', x = 'EC Hypoxia Response',
             caption = paste0('R2 = ', round(rsq, 3))) +
        scale_color_manual(values = Color_palliation[2:5]) +
        theme_classic() +
        theme(aspect.ratio = 1)
p
PlotPDF('5.13.trend.ec_mural_hypoxia_response', 3, 3)
p
dev.off()



##~~  Plot 14 UMAP - Mural Cell States   ####
p <- DimPlot2(mu.srt, group.by = 'Cell_state', reduction = 'sub_umap', cols = mycol_10, pt.size = 0.7) +
        theme(aspect.ratio = 1)
p
PlotPDF('5.14.umap.mu_ell_state', 8, 8)
p
dev.off()



##~~  Plot 15 Box - PC2 Abundance across disease groups   ####
meta <- U(mu.srt@meta.data[, c('group1', 'group2')])
meta <- meta[order(meta$group1),]
mtx <- Table(mu.srt$group1, mu.srt$Cell_state)
mtx <- as.matrix(mtx[, colSums(mtx)>0])
mtx <- mtx/as.vector(table(mu.srt$group1))
df <- data.frame(Cell = mtx[, 'PC2'], group2 = as.factor(meta$group2), value = as.numeric(meta$group2))
wilcox.test(split(df$Cell, df$group2)[[1]], split(df$Cell, df$group2)[[2]])$p.value > 0.05
wilcox.test(split(df$Cell, df$group2)[[3]], split(df$Cell, df$group2)[[2]])$p.value < 0.05
p <- ggplot(df) +
        geom_boxplot(aes(y = Cell, x = group2, fill = group2), outlier.shape = NA) +
        geom_beeswarm(aes(y = Cell, x = group2, color = group2), size = 2, cex = 3) +
        labs(title = 'PC2', y = 'Fraction of Total ECs', x = 'group2') +
        scale_fill_manual(values = Color_palliation[c(1, 2, 5)]) &
        scale_color_manual(values = Color_palliation[c(1, 2, 5)]) &
        theme_classic() &
        theme_Publication(aspect.ratio = 1) &
        NoLegend()
p
PlotPDF('5.15.box.pc2_across_disease_group', 6, 6)
p
dev.off()




##~~  Plot 16 Bar - PC2 abundance across pairs    ####
tmp.srt <- paired.srt[, paired.srt$Cell_type %in% c('PC', 'SMC')]
meta <- U(tmp.srt@meta.data[, c('palliation', 'group2', "vad_pair")])
meta <- meta[order(meta$vad_pair),]
mtx <- Table(tmp.srt$vad_pair, tmp.srt$Cell_state)
mtx <- as.matrix(mtx[, colMaxs(mtx ) > 0])
total <- table(tmp.srt$vad_pair)
mtx <- mtx/as.vector(total[rownames(mtx)])
df <- data.frame(Frac = mtx[, 'PC2'],
                 group2 = as.factor(meta$group2),
                 value = as.numeric(as.factor(meta$group2)),
                 pair = str_split(meta$vad_pair, ' ', simplify = T)[,2])
data3 <- data.frame(Group = c('1', '2', '3'),
                    Fold_reduction = log2(df[df$group2 == 'VAD', ]$Frac/df[df$group2 == 'HLHS', ]$Frac))
data3$Cell_type <- 'Mural'
p <- ggplot(data3) +
        geom_bar(aes(x = Group, y = Fold_reduction, fill = Group), stat = 'identity') +
        scale_fill_manual(values = Color_pair) +
        labs(title = 'PC2', y = 'Log2 FC (Pre vs Post-VAD)', x = 'Pre/Post-VAD Pair') +
        theme_classic() +
        theme(aspect.ratio = 2) +
        NoLegend()
p
PlotPDF('5.16.bar.pc2_abundance_in_pre_post_pair', 4, 6)
p
dev.off()




##~~  Plot 17  Bar  - IL6 LR Sig Spot Composition in ST  ####
st.srt$IL6_Coloc <- GetColocalProb(st.srt, meta_features = c('IL6_score_smooth', 'IL6R_score_smooth'))
st.srt$IL6_Coloc_bin <- st.srt$IL6_Coloc > 2 ## p < 0.01
p <- CountCellBarPlot(st.srt, stack.by = 'IL6_Coloc_bin', group.by = 'group2', percentage = T, cols = mycol_10)
p
PlotPDF('5.17.bar.il6_lr_interaction_probability_in_visium', 4, 4)
p
dev.off()
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
####  Figure 6   ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
##~~  Plot 1 UMAP - FB Cell States   ####
p <- DimPlot2(fb.srt, group.by = 'Cell_state', reduction = 'sub_umap', cols = mycol_10, pt.size = 0.1)
p
PlotPDF('6.01.umap.fb_ell_state', 8, 8)
p
dev.off()



##~~  Plot 2  Feature - FB pseudobulk deg  ####
cutoff <- 2
tmp.srt <- AddModuleScore2(fb.srt, features = list(up, dn), return_z = T, names = c('HLHS_Up', 'HLHS_Down'))
p <- FeaturePlot2(tmp.srt, reduction = 'sub_umap', features = c('HLHS_Up', 'HLHS_Down'), ncol = 2, raster = F,
                     min.cutoff = -cutoff, max.cutoff = cutoff, order = F, pt.size = 0.1) &
        scale_color_distiller(palette = 'RdBu', limits = c(-cutoff, cutoff))
p
PlotPDF('6.02.feature.fb_global_deg', 10, 10)
p
dev.off()



##~~  Plot 3  Scatter - FB PCA Disease Mapping  ####
pc <- readRDS('analysis/PART32.fb_compare_mi_cm.pca_df.rds')
pc$Disease <- factor(pc$Disease, levels = c('Ped_Ctrl', 'Ped_HLHS',
                                            'Adult_Ctrl_1', 'MI_Ischemic', 'MI_Myogenic', 'MI_Fibrotic',
                                            'Adult_Ctrl_2', 'NCCM', 'DCM', 'ARVCM'))
pc$group <- c(rep('Current', 2), rep('MI', 4),  rep('CM', 4))
p <- ggplot(pc) +
        geom_point(aes(x = PC1, y = PC2, color = Disease, shape = group), size = 5) +
        scale_color_manual(values = c(mycol_10[1:9], 'grey50')) +
        labs(color = '', x = '', y = '') +
        labs(x = 'PC1 41% variance', y = 'PC2 20% variance') +
        theme(aspect.ratio = 1,
              axis.line = element_line(size = 0.5, color = 'black'),
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              legend.background = element_blank(),
              panel.background = element_blank())
p
PlotPDF('6.03.pca.mapping_heart_disease', 4, 4)
p
dev.off()



##~~  Plot 4  Box -  Cell states Pan Yap Target expression  ####
include <- c('CM1', 'CM2', 'CM3',
             'FB1', 'FB2', 'FB3', 'FB4', 'FB5',
             'EC1', 'EC2', 'EC3', 'EC4',
             'Mono', 'MP1', 'MP2', 'MP3',
             'PC1', 'PC2')
tmp.srt <- DropMetaLevels(full_clean.srt[, full_clean.srt$Cell_state %in% include])
tmp.srt$Cell_state <- factor(tmp.srt$Cell_state, levels = include)
tmp.srt <- AddModuleScore2(tmp.srt, features = c(Yap_target),
                           names = c('Pan_Yap_target'), return_z = T)
cols <- c('grey80', 'pink')[c(1, 1, 2,
                              1, 1, 1, 1, 2,
                              1, 1, 1, 2,
                              1, 1, 1, 2,
                              1, 2)]
p <- BoxPlot2(tmp.srt, feature = 'Pan_Yap_target', group.by = 'Cell_state', cols = cols) &
        theme(aspect.ratio = 0.3)
p
PlotPDF('6.04.box.yap_target', 9, 3)
p
dev.off()



##~~  Plot 5  Scatter - Trend SASP, Yap   ####
include <- c('CM1', 'CM2', 'CM3',
             'FB1', 'FB2', 'FB3', 'FB4', 'FB5',
             'EC1', 'EC2', 'EC3', 'EC4',
             'PC1', 'PC2',
             'Mono', 'MP1', 'MP2', 'MP3')
ct <- c('CM', 'EC', 'FB', 'Mye', 'PC')
cols <- Color_cell_state[levels(full_clean.srt$Cell_state) %in% include]
meta <- data.frame()
for(i in 1:L(ct)){
        tmp.srt <- DropMetaLevels(full_clean.srt[, full_clean.srt$Cell_state %in% include &
                                                         full_clean.srt$Cell_type %in% ct[i]])
        tmp.srt <- AddModuleScore2(tmp.srt, features = c(sasp, yap_target),
                                   names = c('SASP', 'Pan_Yap_target'), return_z = T)
        meta <- rbind(meta, tmp.srt@meta.data)
}
data <- meta[, c('SASP', 'Pan_Yap_target', 'Cell_state', 'Cell_type')] |>
        group_by(Cell_state) |>
        mutate(Mean_SASP = mean(SASP),
               Mean_YAP = mean(Pan_Yap_target),
               Mean_mTORC1 = mean(mTORC1_Acticity))
data <- U(data[, c('Mean_YAP', 'Mean_SASP', 'Cell_state', 'Cell_type')])
data$Cell_state <- factor(data$Cell_state, levels = include)
p <- ggplot(data) +
        geom_point(aes(x = Mean_SASP, y = Mean_YAP, color = Cell_state)) +
        scale_color_manual(values = cols) +
        theme_Publication(aspect.ratio = 1) +
        scale_y_continuous(limits = c(-0.5, 1.5), breaks = seq(-0.5, 1.5, 0.5)) &
        scale_x_continuous(limits = c(-0.2, 0.8), breaks = seq(-0.2, 0.8, 0.2))
p
PlotPDF('6.05.scatter.sasp_yap_correlation', 3, 3)
p
dev.off()



##~~  Plot 6  Dot - FB4/FB5 Yap target gene examples  ####
p <- DotPlot2(fb.srt, features = c('ACTA2', 'ACTN1', 'SHROOM3', 'SERPINB6', 'ANXA1', 'PTX3', 'CCN1', 'MYC'),
              col.min = 0, col.max = 1.5, cols = 'RdYlBu', group.by = 'Cell_state')
p
PlotPDF('6.06.dot.fb4_fb5_yap_targets', 5, 5)
p
dev.off()


##~~  Plot 8 Box - FB5 Abundance across Disease Groups   ####
meta <- U(fb.srt@meta.data[, c('group1', 'group2')])
meta <- meta[order(meta$group1),]
mtx <- Table(fb.srt$group1, fb.srt$Cell_state)
mtx <- as.matrix(mtx[, colSums(mtx) > 0])
mtx <- mtx/as.vector(table(fb.srt$group1))
df <- data.frame(Cell = mtx[, 'FB5'], group2 = as.factor(meta$group2), value = as.numeric(meta$group2))
wilcox.test(split(df$Cell, df$group2)[[1]], split(df$Cell, df$group2)[[2]])$p.value < 0.05
wilcox.test(split(df$Cell, df$group2)[[3]], split(df$Cell, df$group2)[[2]])$p.value < 0.05
p <- ggplot(df) +
        geom_boxplot(aes(y = Cell, x = group2, fill = group2), outlier.shape = NA) +
        geom_beeswarm(aes(y = Cell, x = group2, color = group2), size = 2, cex = 3) +
        labs(title = 'FB5', y = 'Fraction of Total ECs', x = 'group2') +
        scale_fill_manual(values = Color_palliation[c(1, 2, 5)]) &
        scale_color_manual(values = Color_palliation[c(1, 2, 5)]) &
        theme_classic() &
        theme_Publication(aspect.ratio = 1) &
        NoLegend()
p
PlotPDF('6.08.box.fb5_across_disease_group', 6, 6)
p
dev.off()



##~~  Plot 9 Bar - FB5 Abundance across pairs    ####
tmp.srt <- paired.srt[, paired.srt$Cell_state %in% c(paste0('FB', 1:5))]
meta <- U(tmp.srt@meta.data[, c('palliation', 'group2', "vad_pair")])
meta <- meta[order(meta$vad_pair),]
mtx <- Table(tmp.srt$vad_pair, tmp.srt$Cell_state)
mtx <- as.matrix(mtx[, colMaxs(mtx ) > 0])
total <- table(paired.srt$vad_pair)
mtx <- mtx/as.vector(total[rownames(mtx)])
df <- data.frame(Frac = mtx[, 'FB5'],
                 group2 = as.factor(meta$group2),
                 value = as.numeric(as.factor(meta$group2)),
                 pair = str_split(meta$vad_pair, ' ', simplify = T)[,2])
data2 <- data.frame(Group = c('1', '2', '3'),
                    Fold_reduction = log2(df[df$group2 == 'VAD', ]$Frac/df[df$group2 == 'HLHS', ]$Frac))
data2$Cell_type <- 'FB'
p <- ggplot(data2) +
        geom_bar(aes(x = Group, y = Fold_reduction, fill = Group), stat = 'identity') +
        scale_fill_manual(values = Color_pair) +
        labs(title = 'FB5', y = 'Log2 FC (Pre vs Post-VAD)', x = 'Pre/Post-VAD Pair') +
        theme_classic() +
        theme(aspect.ratio = 2) +
        NoLegend()
p
PlotPDF('6.09.bar.fb5_abundance_in_pre_post_pair', 4, 6)
p
dev.off()



##~~  Plot 10  Feature UMAP - ITGAV and SASP across palliation  ####
tmp.srt <- fb.srt[, unlist(DownsampleByMeta(fb.srt, meta_var = 'palliation', down_to_min_group = T, random = T))]
tmp.srt <- RunALRA(tmp.srt, genes.use = rownames(tmp.srt))
tmp.srt <- AddModuleScore2(tmp.srt, features = sasp, names = 'SASP')
p1 <- FeaturePlot2(tmp.srt, features = 'ITGAV', reduction = 'sub_umap',
                     raster = F, min.cutoff = 3, max.cutoff = 5,  pt.size = 1, split.by = 'palliation') &
        scale_color_distiller(palette = 'RdYlBu', limits = c(3, 5), values = c(0, 0.33, 1)) &
        # scale_color_viridis_c() &
        scale_x_continuous(limits = c(min(tmp.srt@reductions$sub_umap@cell.embeddings[,1]),
                                      max(tmp.srt@reductions$sub_umap@cell.embeddings[,1]))) &
        scale_y_continuous(limits = c(min(tmp.srt@reductions$sub_umap@cell.embeddings[,2]),
                                      max(tmp.srt@reductions$sub_umap@cell.embeddings[,2]))) &
        theme(aspect.ratio = 1, axis.line = element_blank())
p1[[5]] <- p1[[5]] + RestoreLegend()
p2 <- FeaturePlot2(tmp.srt, features = 'SASP', reduction = 'sub_umap',
                     raster = F, min.cutoff = 2, max.cutoff = 4,  pt.size = 1, split.by = 'palliation') &
        scale_color_distiller(palette = 'RdYlBu', limits = c(2, 4)) &
        # scale_color_viridis_c() &
        scale_x_continuous(limits = c(min(tmp.srt@reductions$sub_umap@cell.embeddings[,1]),
                                      max(tmp.srt@reductions$sub_umap@cell.embeddings[,1]))) &
        scale_y_continuous(limits = c(min(tmp.srt@reductions$sub_umap@cell.embeddings[,2]),
                                      max(tmp.srt@reductions$sub_umap@cell.embeddings[,2]))) &
        theme(aspect.ratio = 1, axis.line = element_blank())
p2[[5]] <- p2[[5]] + RestoreLegend()
p <- p1/p2
p
PlotPDF('6.10.feat.fb_ITGAV_score', 10, 5)
p
dev.off()



####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
####  Extended Data Figure 6   ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
##~~  Plot 11  Heatmap - FB pseudobulk deg  ####
deg <- hlhs_deg.list$FB
up <- deg$gene[deg$avg_log2FC > 0]
dn <- deg$gene[deg$avg_log2FC < 0]
psb <- AverageExpression(fb.srt, group.by = 'group1', features = c(up, dn), assays = 'CBN', return.seurat = T)
psb <- ScaleData(psb[, !grepl('^VAD', Cells(psb))], features = rownames(psb))
df <- melt(GetAssayData(psb, layer = 'scale.data')[c(deg$gene[order(deg$avg_log2FC, decreasing = T)][1:50],
                                       deg$gene[order(deg$avg_log2FC, decreasing = F)][1:50]), ])
cutoff <- 2
df$value[df$value > cutoff] <- cutoff
df$value[df$value < -cutoff] <- -cutoff
p <- ggplot(df) +
        geom_tile(aes(x = Var1, y = Var2, fill = value)) +
        scale_fill_distiller(palette = 'RdBu') +
        labs(x = '', y = '', fill = 'Scaled expr.') +
        theme_classic() +
        theme(aspect.ratio = 19/50,
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              axis.line = element_blank())
p
PlotPDF('6.11.heatmap.fb_deg', 10, 3)
p
dev.off()


##~~  Plot 12  Feature PC1/2 gene expression  ####
pc_genes <- readRDS('analysis/PART32.fb_compare_mi_cm_pc1_pc2_top_genes.list.rds')
p <- FeaturePlot3(fb.srt, features = 'Pan_Compare_PC2', adjust = 1.5, rescale_neg = T, pt.size = 1.5)
p
PlotPDF('6.12.feature_density.fb_mapping_heart_disease_pc1_pc2_expression', 4, 3)
p
dev.off()
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
####  Figure 7    ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
##~~  Plot 1 UMAP - FB Cell States   ####
p <- DimPlot2(mye.srt, group.by = 'Cell_state', reduction = 'sub_umap', cols = mycol_10, pt.size = 0.2)
p
PlotPDF('7.01.umap.mye_ell_state', 8, 8)
p
dev.off()



##~~  Plot 2 Box - MP3 Abundance across Disease Groups   ####
meta <- U(mye.srt@meta.data[, c('group1', 'group2')])
meta <- meta[order(meta$group1),]
mtx <- Table(mye.srt$group1, mye.srt$Cell_state)
mtx <- as.matrix(mtx[, colSums(mtx)>0])
mtx <- mtx/as.vector(table(mye.srt$group1))
df <- data.frame(Cell = mtx[, 'MP3'], group2 = as.factor(meta$group2), value = as.numeric(meta$group2))
wilcox.test(split(df$Cell, df$group2)[[1]], split(df$Cell, df$group2)[[2]])$p.value < 0.05
wilcox.test(split(df$Cell, df$group2)[[3]], split(df$Cell, df$group2)[[2]])$p.value < 0.05
p <- ggplot(df) +
        geom_boxplot(aes(y = Cell, x = group2, fill = group2), outlier.shape = NA) +
        geom_beeswarm(aes(y = Cell, x = group2, color = group2), size = 2, cex = 3) +
        labs(title = 'MP3', y = 'Fraction of Total Myeloid', x = 'group2') +
        scale_fill_manual(values = Color_palliation[c(1, 2, 5)]) &
        scale_color_manual(values = Color_palliation[c(1, 2, 5)]) &
        theme_classic() &
        theme_Publication(aspect.ratio = 1) &
        NoLegend()
p
PlotPDF('7.02.box.mp3_across_disease_group', 6, 6)
p
dev.off()



##~~  Plot 3 Bar - MP3 Abundance across pairs    ####
tmp.srt <- paired.srt[, paired.srt$Cell_type == 'Mye']
meta <- U(tmp.srt@meta.data[, c('palliation', 'group2', "vad_pair")])
meta <- meta[order(meta$vad_pair),]
mtx <- Table(tmp.srt$vad_pair, tmp.srt$Cell_state)
mtx <- as.matrix(mtx[, colMaxs(mtx ) > 0])
total <- table(tmp.srt$vad_pair)
mtx <- mtx/as.vector(total[rownames(mtx)])
df <- data.frame(Frac = mtx[, 'MP3'],
                 group2 = as.factor(meta$group2),
                 value = as.numeric(as.factor(meta$group2)),
                 pair = str_split(meta$vad_pair, ' ', simplify = T)[,2])
data4 <- data.frame(Group = c('1', '2', '3'),
                    Fold_reduction = log2(df[df$group2 == 'VAD', ]$Frac/df[df$group2 == 'HLHS', ]$Frac))
data4$Cell_type <- 'Mye'
p <- ggplot(data4) +
        geom_bar(aes(x = Group, y = Fold_reduction, fill = Group), stat = 'identity') +
        scale_fill_manual(values = Color_pair) +
        labs(title = 'MP3', y = 'Log2 FC (Pre vs Post-VAD)', x = 'Pre/Post-VAD Pair') +
        theme_classic() +
        theme(aspect.ratio = 2) +
        NoLegend()
p
PlotPDF('7.03.bar.mp3_abundance_in_pre_post_pair', 4, 6)
p
dev.off()



##~~  Plot 6 Bar - MP2 signature enrichment   ####
deg <- FindMarkers(momp.srt, ident.1 = 'MP2')
deg <- deg[deg$p_val_adj < 0.05, ]
deg$gene <- rownames(deg)
deg$direction <- ifelse(deg$avg_log2FC > 0, 'MP2', 'Other')
gl <- split(deg$gene, deg$direction)
enrich <- ModuleEnrichment(gl, 'human')
p <- EnrichBarPlot(enrich_df = enrich$GO$GO_MP2, terms = c(
        'lamellipodium organization',
        'filopodium assembly',
        'phagocytosis, engulfment',
        'phagocytosis',
        'macroautophagy',
        'endosomal transport',
        'lysosome organization'
))
p
PlotPDF('7.06.bar.mp2_function_enrich', 4, 4)
p
dev.off()



##~~  Plot 7  Density UMAP - MP enrichemnet over palliation  ####
p_list <- list()
tmp.srt <- momp.srt[, unlist(DownsampleByMeta(momp.srt, meta_var = 'palliation', down_to_min_group = T, random = T))]
Idents(tmp.srt) <- 'palliation'
df2 <- data.frame(tmp.srt@reductions$sub_main_umap@cell.embeddings)
df2$palliation <- tmp.srt$palliation
df3 <- data.frame(momp.srt@reductions$sub_main_umap@cell.embeddings)
df3$palliation <- momp.srt$palliation
alpha <- 0.3
smooth <- 3
bins <- 15
for(i in 1:5){
        p_list[[i]] <- ggplot(df2, aes(x = mSUBUMAP_1, y = mSUBUMAP_2)) +
                geom_point(data = df3,
                           color = 'black', alpha = 0.05, size = 0.1) +
                theme_classic() +
                scale_x_continuous(limits = c(-9, 8)) +
                scale_y_continuous(limits = c(-8, 6)) +
                stat_density_2d(data = df2[df2$palliation == levels(df2$palliation)[i],],
                                aes(fill = after_stat(level)),
                                geom = "polygon", contour = T, bins = bins, h = c(smooth, smooth),
                                alpha = alpha, size = 0) +
                scale_fill_viridis() +
                theme(aspect.ratio = 1, axis.text = element_blank(), axis.ticks = element_blank(),
                      axis.title = element_blank()) +
                labs(title = levels(tmp.srt$palliation)[i]) +
                NoLegend()
}
p <- wrap_plots(p_list, nrow = 1)
p
PlotPDF('7.07.density_umap.mp2_mp3_transition_across_palliation', 15, 3)
p
dev.off()



####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
####  Extended Data Figure 7   ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
##~~  Plot 8  Feature - Known MP markers  ####
p <- FeaturePlot3(mye.srt, features = c('Impute_CCR2', 'TFL', 'DC2', 'Mo_CD14',
                                        'LAM', 'MP_LYVE1_hi', 'MP_HLA_hi'),
                  rescale_neg = T, ncol = 7, adjust = 2, pt.size = 0.75)
p
PlotPDF('7.08.feat_density.known_markers', 21, 2.5)
p
dev.off()
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
####  Figure 8  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
##~~  Plot 1  Heatmap - subset of ligands from senescent niche cells  ####
data <- readRDS('analysis/PART26.sen_niche_ligands.dataframe.rds')
data <- data[rownames(data) %in% hlhs.sen.cch@LR$LRsig$ligand, ]
data <- t(log1p(as.matrix(data)))
viridis_colors <- viridis(100, option = 'D')
data2 <- data[, colMaxs(data) > 2 ]
data2 <- cbind(data[, order(data[1,], decreasing = T)][, 1:10],
               data[, order(data[2,], decreasing = T)][, 1:10],
               data[, order(data[3,], decreasing = T)][, 1:10],
               data[, order(data[4,], decreasing = T)][, 1:10])
data2 <- data2[, !duplicated(colnames(data2))]
rst <- pheatmap(t(data2), scale = 'none', cluster_rows = T, cluster_cols = T,
                clustering_method = 'average', clustering_distance_cols = 'manhattan',
                cellwidth = 10, cellheight = 10,
                color = viridis_colors)
intersect(rst$tree_row$labels[rst$tree_row$order], sasp[[1]])
PlotPDF('8.01.heatmap.selective_sen_niche_ligands', 5, 10)
rst
dev.off()


##~~  Plot 2  Chord - Flow of HLHS-enriched signals  ####
data <- readRDS('analysis/PART26.sen_niche_ligands.dataframe.rds')
sig <- readRDS('analysis/PART26.sen_niche_ligands_state_marker_bin.dataframe.rds')
PlotPDF('8.02.chord.sen_niche_signaling', 20, 20)
lig <- U(as.vector(intersect(rownames(data), rownames(sig)[sig$FB5])))
pw.df <- hlhs.sen.cch@LR$LRsig[hlhs.sen.cch@LR$LRsig$ligand %in% rownames(data), ]
pw.df <- pw.df[!pw.df$interaction_name %in% c('NAMPT_CCR5','NAMPT_ITGA5_ITGB1','NAMPT_INSR'), ]
top <- sort(table(U(pw.df[,c('pathway_name', 'ligand')])$pathway_name), decreasing = T)
netVisual_chord_gene(hlhs.sen.cch,
                     targets.use = c('CM', 'FB', 'EC', 'PC', 'SMC', 'Lym', 'Mye'),
                     sources.use = 'FB5',
                     slot.name = "netP",
                     pairLR.use = pw.df,
                     #signaling = names(top),
                     link.target.prop = F,
                     scale = F,
                     reduce = 0.01,
                     title.name = 'FB5 Signaling',
                     color.use = c(Color_cell_type, 'grey', 'black'),
                     thresh = 0.001,
                     annotationTrackHeight = 0.05,
                     link.border= T, lab.cex = 0.5, legend.pos.x = 10, legend.pos.y = 30)

lig <- U(as.vector(intersect(rownames(data), rownames(sig)[sig$EC4])))
pw.df <- hlhs.sen.cch@LR$LRsig[hlhs.sen.cch@LR$LRsig$ligand %in% lig, ]
pw.df <- pw.df[!pw.df$interaction_name %in% c('NAMPT_CCR5','NAMPT_ITGA5_ITGB1','NAMPT_INSR'), ]
top <- sort(table(U(pw.df[,c('pathway_name', 'ligand')])$pathway_name), decreasing = T)
netVisual_chord_gene(hlhs.sen.cch,
                     targets.use = c('CM', 'FB', 'EC', 'PC', 'SMC', 'Lym', 'Mye'),
                     sources.use = 'EC4',
                     slot.name = "netP",
                     pairLR.use = pw.df,
                     #signaling = names(top),
                     link.target.prop = F,
                     scale = F,
                     reduce = 0.01,
                     title.name = 'EC4 Signaling',
                     color.use = c(Color_cell_type, 'grey', 'black'),
                     thresh = 0.001,
                     annotationTrackHeight = 0.05,
                     link.border= T, lab.cex = 0.5, legend.pos.x = 10, legend.pos.y = 30)

lig <- U(as.vector(intersect(rownames(data), rownames(sig)[sig$MP3])))
pw.df <- hlhs.sen.cch@LR$LRsig[hlhs.sen.cch@LR$LRsig$ligand %in% lig, ]
pw.df <- pw.df[!pw.df$interaction_name %in% c('NAMPT_CCR5','NAMPT_ITGA5_ITGB1','NAMPT_INSR'), ]
top <- sort(table(U(pw.df[,c('pathway_name', 'ligand')])$pathway_name), decreasing = T)
netVisual_chord_gene(hlhs.sen.cch,
                     targets.use = c('CM', 'FB', 'EC', 'PC', 'SMC', 'Lym', 'Mye'),
                     sources.use = 'MP3',
                     slot.name = "netP",
                     pairLR.use = pw.df,
                     #signaling = names(top),
                     link.target.prop = F,
                     scale = F,
                     reduce = 0.01,
                     title.name = 'MP3 Signaling',
                     color.use = c(Color_cell_type, 'grey', 'black'),
                     thresh = 0.001,
                     annotationTrackHeight = 0.05,
                     link.border= T, lab.cex = 0.5, legend.pos.x = 10, legend.pos.y = 30)

lig <- U(as.vector(intersect(rownames(data), rownames(sig)[sig$PC2])))
pw.df <- hlhs.sen.cch@LR$LRsig[hlhs.sen.cch@LR$LRsig$ligand %in% lig, ]
pw.df <- pw.df[!pw.df$interaction_name %in% c('NAMPT_CCR5','NAMPT_ITGA5_ITGB1','NAMPT_INSR'), ]
top <- sort(table(U(pw.df[,c('pathway_name', 'ligand')])$pathway_name), decreasing = T)
netVisual_chord_gene(hlhs.sen.cch,
                     targets.use = c('CM', 'FB', 'EC', 'PC', 'SMC', 'Lym', 'Mye'),
                     sources.use = 'PC2',
                     slot.name = "netP",
                     pairLR.use = pw.df,
                     #signaling = names(top),
                     link.target.prop = F,
                     scale = F,
                     reduce = 0.01,
                     title.name = 'PC2 Signaling',
                     color.use = c(Color_cell_type, 'grey', 'black'),
                     thresh = 0.001,
                     annotationTrackHeight = 0.05,
                     link.border= T, lab.cex = 0.5, legend.pos.x = 10, legend.pos.y = 30)
dev.off()



##~~  Plot 3  Vortex - Postn Wnt9 Nampt signaling  ####
PlotPDF('8.03.vortex.individual_signaling', 5, 5)
netVisual_individual(hlhs.sen.cch,
                     targets.use = c('CM', 'FB', 'EC', 'PC', 'SMC', 'Lym', 'Mye'),
                     sources.use = c('FB5', 'MP3', 'EC4', 'PC2'),
                     pairLR.use = 'POSTN_ITGAV_ITGB5',
                     signaling = 'PERIOSTIN',
                     color.use = c(Color_cell_type, 'grey', 'black'),
                     layout = "circle")
netVisual_individual(hlhs.sen.cch,
                     targets.use = c('CM', 'FB', 'EC', 'PC', 'SMC', 'Lym', 'Mye'),
                     sources.use = c('FB5', 'MP3', 'EC4', 'PC2'),
                     pairLR.use = 'NAMPT_TLR4',
                     signaling = 'VISFATIN',
                     color.use = c(Color_cell_type, 'grey', 'black'),
                     layout = "circle")
dev.off()



##~~  Plot 2  Dot - Per patient NAMPT/POSTN expression in ST   - NCVR Revision V1   ####
vi_mean <- AverageExpression(st.srt, features = c('NAMPT', 'POSTN'), group.by = 'group1')$ST
vi_mean[1,] <- scale(vi_mean[1,])
vi_mean[2,] <- scale(vi_mean[2,])
xe_mean <- AverageExpression(xen.srt, features = c('NAMPT', 'POSTN'), group.by = 'Group1')$Xenium
xe_mean[1,] <- scale(xe_mean[1,])
xe_mean[2,] <- scale(xe_mean[2,])

data1 <- reshape2::melt(t(vi_mean), value.name = c('NAMPT', 'POSTN'))
data1$Group2 <- str_split(rownames(data1), pattern = '-', simplify = T)[, 1]
data1$Group1 <- rownames(data1)
data1$Type <- 'Visium'
data2 <- reshape2::melt(t(xe_mean), value.name = c('NAMPT', 'POSTN'))
data2$Group2 <- str_split(rownames(data2), pattern = '-', simplify = T)[, 1]
data2$Group1 <- rownames(data2)
data2$Type <- 'Xenium'
data <- rbind(data1, data2)
p1 <- ggplot(data, aes(x = Group2, y = POSTN, color = Type)) +
        geom_beeswarm() +
        ggrepel::geom_label_repel(aes(label = Group1)) +
        scale_y_continuous(limits = c(-1,2)) +
        theme_Publication(aspect.ratio = 1)
p2 <- ggplot(data, aes(x = Group2, y = NAMPT, color = Type)) +
        geom_beeswarm() +
        ggrepel::geom_label_repel(aes(label = Group1)) +
        scale_y_continuous(limits = c(-1,2)) +
        theme_Publication(aspect.ratio = 1)
p1 + p2
PlotPDF('1.02.dot.postn_nampt_expression_in_visium_and_xenium', 8, 4)
p1 + p2
dev.off()



##~~  Plot 6  ST Trend - Sen Niche Neighbor NAMPT, POSTN  ####
data <- st.srt@meta.data
data$NAMPT_score <- st.srt@assays$ST@data['NAMPT',]
data$POSTN_score <- st.srt@assays$ST@data['POSTN',]
data1 <- data |> group_by(dNiche, group2) |> mutate(MEAN = mean(NAMPT_score), Gene = 'NAMPT')
data2 <- data |> group_by(dNiche, group2) |> mutate(MEAN = mean(POSTN_score), Gene = 'POSTN')
data <- rbind(data1,  data2)
data <- U(data[, c('group2', 'dNiche', 'MEAN', 'Gene')])
p <- ggplot(data) +
        geom_point(aes(x = dNiche, y = MEAN, color = dNiche)) +
        geom_line(aes(x = dNiche, y = MEAN, group = group2)) +
        scale_color_manual(values = Color_dniche) +
        scale_y_continuous(limits = c(0, 0.6)) +
        facet_wrap(Gene~group2) +
        theme_classic() +
        theme(aspect.ratio = 2)
p
PlotPDF('8.06.trend.ligand_in_sen_niche_neighbors', 6, 12)
p
dev.off()


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
####  Extended Data Figure 8   ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
##~~  Plot 7  Heatmap - full ligands from senescent niche cells  ####
data <- readRDS('analysis/PART26.sen_niche_ligands.dataframe.rds')
data <- data[rownames(data) %in% hlhs.sen.cch@LR$LRsig$ligand, ]
data <- t(log1p(as.matrix(data)))
viridis_colors <- viridis(100, option = 'D')
rst <- pheatmap(data, scale = 'none', cluster_rows = T, cluster_cols = T,
                clustering_method = 'average', clustering_distance_cols = 'manhattan',
                cellwidth = 10, cellheight = 10,
                color = viridis_colors)
intersect(rst$tree_col$labels[rst$tree_col$order], sasp[[1]])
PlotPDF('8.07.heatmap.full_sen_niche_ligands', 15, 5)
rst
dev.off()


##~~  Plot 8  Heatmap - Ligand receptor expression of enriched signals  ####
tmp.srt <- DropMetaLevels(full_clean.srt[, ! full_clean.srt$Cell_type %in% c('AC', 'NC')])
full_clean.srt$tmp <- paste(full_clean.srt$group2, full_clean.srt$Cell_type)
full_clean.srt$tmp <- factor(full_clean.srt$tmp,
                      levels = c(paste('Control', levels(full_clean.srt$Cell_type)),
                                 paste('HLHS', levels(full_clean.srt$Cell_type)),
                                 paste('VAD', levels(full_clean.srt$Cell_type))))
p1 <- MarkerHeatmapMean(full_clean.srt, features = 'POSTN', group.by = 'tmp', col.max = 2.4, col.min = 1.2) +
        scale_fill_viridis(limits = c(1.2, 2.4)) +
        theme(aspect.ratio = 1/21)
p2 <- MarkerHeatmapMean(full_clean.srt, features = c('ITGAV', 'ITGB5'), group.by = 'tmp', col.max = 1.5, col.min = -1.5) +
        scale_fill_viridis(limits = c(-1.5, 1.5)) +
        theme(aspect.ratio = 2/21)
p3 <- MarkerHeatmapMean(full_clean.srt, features = 'NAMPT', group.by = 'tmp', col.max = 2.8, col.min = -0.8) +
        scale_fill_viridis(limits = c(-0.8, 2.8)) +
        theme(aspect.ratio = 1/21)
p4 <- MarkerHeatmapMean(full_clean.srt, features = c('TLR4'), group.by = 'tmp', col.max = 2.5, col.min = -1) +
        scale_fill_viridis(limits = c(-1, 2.5)) +
        theme(aspect.ratio = 1/21)
p <- p1/p2/p3/p4 &
        RotatedAxis()
p
PlotPDF('8.08.heat.example_enriched_ligand_receptor_expression', 12, 10)
p
dev.off()
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
####  Method Figure | Data QC   ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
##~~  Plot 1  UMAP - Ambiguous cells before filter cell type  ####
pct <- (table(full.srt$Cell_type)[['MitoHi']]+table(full.srt$Cell_type)[['Doublet']])/ncol(full.srt)
p <- DimPlot2(full.srt, group.by = 'Cell_type', cols = c(Color_cell_type, 'blue2', 'black'),
               reduction = 'full_umap', raster = F) +
        labs(x = 'UMAP 1', y = 'UMAP 2', title = paste0('Ambiguous Cells (', round(pct*100, 1),'%)'))
p
PlotPDF('9.01.umap.global_cell_type_with_ambiguous', 15, 15)
p
dev.off()



##~~  Plot 2  UMAP - Scrublet doublets  ####
pct <- table(full.srt$Doublet_SC)[[2]]/ncol(full.srt)
p <- DimPlot2(full.srt, group.by = 'Doublet_SC', cols = c('grey80', 'black'),
               reduction = 'full_umap', raster = F) +
        labs(x = 'UMAP 1', y = 'UMAP 2', title = paste0('Scrublet Doublets (', round(pct*100, 1),'%)'))
p
PlotPDF('9.02.umap.global_srublet_doublets', 15, 15)
p
dev.off()



###~~  Plot 3  Violin - Per cell UMI nGene and mito% QC  ####
tmp.df <- readRDS('integrated/PART02.merged.cbn.srt_meta.rds')
cutoffs <- readRDS('analysis/PART03.cell_filtered.low_quality_cutoffs.df.rds')
tmp.df <- inner_join(x = tmp.df, y = lib_meta.df, by = c('sample' = 'Sample_id'))
tmp.df <- inner_join(x = tmp.df, y = cutoffs, by = c('sample' = 'name'))
tmp.df$Palliation[is.na(tmp.df$Palliation)] <- 'Control'
tmp.df$Palliation <- factor(tmp.df$Palliation, levels = c('Control', 'I', 'II', 'III', 'VAD'))
tmp.df$group3 <- paste(tmp.df$group1, tmp.df$Replicate)
tmp.df$group3 <- str_replace(tmp.df$group3, pattern = '_0', replacement = ' ')
tmp.df$group3 <- str_replace(tmp.df$group3, pattern = '_', replacement = ' ')
tmp.df$group3 <- factor(tmp.df$group3, labels = str_sort(U(tmp.df$group3), numeric = T))
p1 <- ggplot(tmp.df) +
        geom_violin(aes(x = group3, y = log10(nCount_CBN), fill = group2)) +
        labs(x = '', y = 'Log10 UMI count per nucleus') +
        geom_point(data = tmp.df[!duplicated(tmp.df$sample), ],
                   aes(x = group3, y = log10(nCount_max)), shape = 3, size = 2) +
        geom_point(data = tmp.df[!duplicated(tmp.df$sample), ],
                   aes(x = group3, y = log10(nCount_min)), shape = 3, size = 2)
p2 <- ggplot(tmp.df) +
        geom_violin(aes(x = group3, y = nFeature_CBN/1e3, fill = group2)) +
        scale_y_continuous(limits = c(0, 15)) +
        labs(x = '', y = 'Gene count per nucleus') +
        geom_point(data = tmp.df[!duplicated(tmp.df$sample), ],
                   aes(x = group3, y = nFeature_max/1e3), shape = 3, size = 2) +
        geom_point(data = tmp.df[!duplicated(tmp.df$sample), ],
                   aes(x = group3, y = nFeature_min/1e3), shape = 3, size = 2)
p3 <- ggplot(tmp.df) +
        geom_violin(aes(x = group3, y = pct_mito_CBN, fill = group2), scale = 'width') +
        scale_y_continuous(limits = c(0, 10)) +
        labs(x = '', y = 'Mitochondrial UMI % per nucleus') +
        geom_point(data = tmp.df[!duplicated(tmp.df$sample), ],
                   aes(x = group3, y = 5), shape = 3, size = 2)

p <- p1/p2/p3 &
        scale_fill_manual(values = Color_disease[c(1, 2, 2)]) &
        theme_classic() +
        theme(aspect.ratio = 0.33) &
        RotatedAxis()
p
PlotPDF('9.03.vln.ncount_nfeature', 10, 10)
p
dev.off()



###~~  Plot 4  Bar - Cell cluster distribution across samples  ####
p4.1 <- MappingHeatmap(full.srt, que_var = 'Cell_type', ref_var = 'group1',
                       percentage = T, log10_scale = T, center = F,
                       que_order = levels(full.srt$Cell_type),
                       ref_order = levels(full.srt$group1), ref.disp.min = 0)
p4.2 <- MappingHeatmap(full.srt, que_var = 'Cell_state', ref_var = 'group1',
                       percentage = T, log10_scale = T, center = F,
                       que_order = levels(full.srt$Cell_state),
                       ref_order = levels(full.srt$group1), ref.disp.min = 0)
data1 <- p4.1$data |> group_by(Query) |> transmute(max(Value)) |> distinct()
data2 <- p4.2$data |> group_by(Query) |> transmute(max(Value)) |> distinct()
p1 <- ggplot(data1) +
        geom_bar(stat = 'identity', aes(x = Query, y = `max(Value)`, fill = `max(Value)` < 0.5)) +
        geom_hline(yintercept = 0.5) +
        theme_classic() +
        theme(aspect.ratio = 0.3) +
        NoLegend()
p2 <- ggplot(data2) +
        geom_bar(stat = 'identity', aes(x = Query, y = `max(Value)`, fill = `max(Value)` < 0.5)) +
        geom_hline(yintercept = 0.5) +
        theme_classic() +
        theme(aspect.ratio = 0.3/(35/12))
p <- p1 + p2 &
        RotatedAxis()
p
PlotPDF('9.04.bar.cell_type_vs_donor_distribution', 10, 3)
p
dev.off()



##~~  Plot 8  Bar - Silhouette score for clusters   ####
umap_df <- as.data.frame(full_clean.srt@reductions$clean_umap@cell.embeddings)
umap_df$Cell_type <- full_clean.srt$Cell_type
umap_df$Cell_state <- full_clean.srt$Cell_state
umap_df$Group1 <- full_clean.srt$group1
umap_df$subUMAP_1 <- full_clean.srt@reductions$sub_umap@cell.embeddings[,1]
umap_df$subUMAP_2 <- full_clean.srt@reductions$sub_umap@cell.embeddings[,2]

## Global downsampled:
sub <- umap_df[sample(1:nrow(umap_df), size = 30e3), ]
dist_mat <- Rfast::Dist(as.matrix(sub[, c("cUMAP_1", "cUMAP_2")]))
sil <- silhouette(as.numeric(sub$Cell_type), dist_mat)
global_sil_score <- mean(sil[, "sil_width"])  # Higher value (~1) indicates better separation
## 0.5893597

## Cell State Level:
ct <- list('CM', 'FB', 'EC', 'Mye', c('PC', 'SMC'))
Sil_score <- c()
for(i in 1:L(ct)) {
        df <- sub[sub$Cell_type %in% ct[[i]], ]
        sil <- silhouette(as.numeric(df$Cell_state),
                          Rfast::Dist(as.matrix(df[, c("subUMAP_1", "subUMAP_2")])))
        Sil_score[i] <- mean(sil[, "sil_width"])  # Higher value (~1) indicates better separation
}
names(Sil_score) <- c('CM', 'FB', 'EC', 'Mye', 'Mural')
Sil_score > 0.25 ## TRUE
# Sil_score
# CM        FB        EC       Mye     Mural
# 0.2502551 0.2627522 0.3410668 0.3415686 0.3251319

data <- data.frame('Cluster Level' = c('Cell type', 'CM state', 'FB state', 'EC state', 'Mye state', 'Mural state'),
                   'Silhouette Score' = c(global_sil_score, Sil_score))
p <- ggplot(data) +
        geom_bar(aes(x = Cluster.Level, y = Silhouette.Score), stat = 'identity') +
        geom_hline(yintercept = 0.25) +
        theme_Publication(aspect.ratio = 1)
p
PlotPDF('1.08.bar.cluster_silhouette_score', 4, 4)
p
dev.off()
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
