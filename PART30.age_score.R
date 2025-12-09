####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Pediatric HLHS Compendium
####  2022-09-01 by Xiao LI (Texas Heart Institute, US)
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Initiate  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
Ver <- '2'
Step <- 'PART30_Age_Score'
Project <- '2022_hlhs_dturaga'
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Load Data  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
full.srt <- readRDS('integrated/PART19.annotated.srt.rds')
full.srt <- full.srt[, full.srt$Non_ambiguous]
cm.srt <- full.srt[, full.srt$Cell_type == 'CM' & full.srt$group2 != 'VAD']
fb.srt <- full.srt[, full.srt$Cell_type == 'FB' & full.srt$group2 != 'VAD']
mye.srt <- full.srt[, full.srt$Cell_type == 'Mye' & full.srt$group2 != 'VAD']

bulk.mtx <- read.table(paste0('/Volumes/shire/data/rnaseq/2021_Circulation_EPorrello/matrix_public/',
                              'GSE156702_hrna_dev_mf_fulllen_se_strrev_q30.mx.all.fix_filt.csv'), header = T, sep = ',')
bulk.mtx <- bulk.mtx[!duplicated(str_split(bulk.mtx$Geneid, pattern = '_', simplify = T)[, 2]), ]
rownames(bulk.mtx) <- str_split(bulk.mtx$Geneid, pattern = '_', simplify = T)[, 2]
bulk.mtx$Geneid <- NULL
gc()
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Find DEG between Young vs Adult  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
## Check data variance
meta <- data.frame(Name = colnames(bulk.mtx),
                   Group = c(rep('Fetal', 6), rep('Young', 4), rep('Adult', 11)),
                   Age = c(rep(0, 6), 0.05, 0.2, 2, 4, 45, 40, 41, 39, 48, 42, 54, 60, 62, 62, 65),
                   Sex = c('F', 'M')[c(1,1,1,2,2,1,  2,1,2,2,  2,2,1,2,1,1,1,2,2,1,1)]
)
rownames(meta) <- meta$Name
all_dds <- DESeqDataSetFromMatrix(bulk.mtx, colData = meta, design = ~ Group)
all_dds <- estimateSizeFactors(all_dds)
all_dds <- DESeq(all_dds)
all_rld <- rlog(all_dds, blind = T)


## Calculate DE
sub.mtx <- bulk.mtx[, 7:ncol(bulk.mtx)]
sub.meta <- meta[7:nrow(meta), ]
dds <- DESeqDataSetFromMatrix(sub.mtx, colData = sub.meta, design = ~ Group)
dds <- estimateSizeFactors(dds)
dds <- DESeq(dds)
res <- results(dds, alpha = 0.05, contrast = c('Group', 'Young', 'Adult'))
res <- lfcShrink(dds, contrast = c('Group', 'Young', 'Adult'), res = res, type = 'normal')
deg_up <- rownames(res)[res$padj < 0.05 & res$log2FoldChange > 0.5]
deg_dn <- rownames(res)[res$padj < 0.05 & res$log2FoldChange < -0.5]
deg_up <- intersect(deg_up, rownames(cm.srt))
deg_dn <- intersect(deg_dn, rownames(cm.srt))
deg_list <- list('Young' = deg_up, 'Adult' = deg_dn)
str(deg_list)

norm_count <- counts(all_dds, normalized=TRUE)
deg.mtx <- as.data.frame(log(norm_count[c(deg_up, deg_dn), ]))
deg.mtx$Gene <- rownames(deg.mtx)
deg.mtx$DEG <- c(rep('Young', L(deg_up)), rep('Adult', L(deg_dn)))
deg.df <- melt(deg.mtx)
deg.df$Age <- as.numeric(as.vector(mapvalues(deg.df$variable, from = sub.meta$Name, to = sub.meta$Age)))

deg_list$Adult <- deg_list$Adult[! grepl('\\.', deg_list$Adult)]
deg_list$Young <- deg_list$Young[! grepl('\\.', deg_list$Young)]
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
saveRDS(deg_list, 'analysis/PART30.deseq2_sim_et_al_young_adult_cm_signautures.list.rds')
WriteCSV(data.frame('Young' = paste0("'", deg_list$Young)), 'PART30.cm_bulk_young_genes')
WriteCSV(data.frame('Adult' = paste0("'", deg_list$Adult)), 'PART30.cm_bulk_adult_genes')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Score CMs with Adult and Young CM signatures in snRNA-seq Data + Bulk Data  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
## Include Bulk RNA data
bulk.srt <- CreateSeuratObject(counts = bulk.mtx[, 7:10],
                               meta.data = data.frame(row.names = c('Young1', 'Young2', 'Young3', 'Young4'),
                                                      group1 = c('BCtrl1', 'BCtrl2', 'BCtrl3', 'BCtrl4'),
                                                      group2 = 'Control',
                                                      age = c(0.05, 0.2, 2, 4)),
                               min.cells = 0, min.features = 0, assay = 'CBN')
bulk.srt <- NormalizeData(bulk.srt)
cm.srt2 <- merge(cm.srt, bulk.srt)
cm.srt2 <- JoinLayers(cm.srt2)
cm.srt2 <- AddModuleScore2(cm.srt2, features = deg_list, names = c('Young_score', 'Adult_score'), return_z = T)

mean_exp <- cm.srt2@meta.data[, c('Young_score', 'Adult_score', 'group1', 'group2', 'age')] |>
        dplyr::group_by(group1, age) |>
        dplyr::mutate(Mean_Young = mean(Young_score), Mean_Adult = mean(Adult_score),
                      Mean_Score = mean(Adult_score)-mean(Young_score))
mean_exp$age <- as.numeric(as.vector(mean_exp$age))
mean_exp <- mean_exp[!duplicated(mean_exp[, 3:ncol(mean_exp)]), 3:ncol(mean_exp)]
mean_exp <- mean_exp[order(mean_exp$age), ]
mean_exp <- reshape2::melt(mean_exp, id = 1:3)

saveRDS(mean_exp, 'analysis/PART30.cm_and_bulk_young_adult_score.srt_meta.rds')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
