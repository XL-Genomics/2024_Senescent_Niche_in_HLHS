####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Pediatric HLHS Compendium
####  2022-09-01 by Xiao LI (Texas Heart Institute, US)
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Initiate  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
Ver <- '2'
Step <- 'PART27_Pseudobulk_PCA'
Project <- '2022_hlhs_dturaga'
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Load HLHS Data  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
full.srt <- readRDS('integrated/PART19.annotated.srt.rds')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Global PCA  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
DefaultAssay(full.srt) <- 'CBN'
tmp.srt <- RenameAssays(full.srt, 'CBN' = 'RNA')
tmp.srt <- tmp.srt[, tmp.srt$Non_ambiguous]
# tmp.srt <- tmp.srt[, tmp.srt$Cell_type != 'EpiC']
tmp.srt <- NormalizeData(tmp.srt)
tmp.srt@meta.data <- tmp.srt@meta.data[, c('group1', 'group2', 'age', 'sex', 'palliation', 'study',
                                           'rv_depression', 'log10BNP', 'avv_regurg')]
tmp.srt <- FindVariableFeatures(tmp.srt, nfeatures = 1000)

psb <- AggregateExpression(tmp.srt, assays = 'RNA', slot = 'count', group.by = 'group1',
                           features = VariableFeatures(tmp.srt))$RNA
meta <- tmp.srt@meta.data[!duplicated(tmp.srt$group1), ]
rownames(meta) <- meta$group1
meta <- meta[colnames(psb), ]
meta$age <- as.numeric(as.vector(meta$age))
dds <- DESeqDataSetFromMatrix(psb, colData = meta, design = ~ group1)
rld <- rlog(dds, blind = T)
plotPCA(rld, intgroup = "palliation", ntop = 1000) +
        scale_color_manual(values = c(mycol_10)) +
        theme_classic() +
        theme(aspect.ratio = 1)
plotPCA(rld, intgroup = "age", ntop = 1000) +
        scale_color_distiller(palette = 'Spectral') +
        theme_classic() +
        theme(aspect.ratio = 1)
plotPCA(rld, intgroup = "sex", ntop = 1000) +
        scale_color_manual(values = c(mycol_10)) +
        theme_classic() +
        theme(aspect.ratio = 1)
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
saveRDS(p$data, 'analysis/PART27.ctrl_hlhs_all_cell_pseudobulk_pca.df.rds')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  CM main source of variance  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
tmp.srt <- RenameAssays(cm.srt, 'CBN' = 'RNA')
tmp.srt <- NormalizeData(tmp.srt)
tmp.srt@meta.data <- tmp.srt@meta.data[, c('group1', 'group2', 'age', 'sex', 'palliation', 'study')]
tmp.srt <- FindVariableFeatures(tmp.srt)

psb <- AggregateExpression(tmp.srt, assays = 'RNA', slot = 'count', group.by = 'group1',
                           features = VariableFeatures(tmp.srt))$RNA
meta <- tmp.srt@meta.data[!duplicated(tmp.srt$group1), ]
rownames(meta) <- meta$group1
meta <- meta[colnames(psb), ]
meta$age <- as.numeric(meta$age)
dds <- DESeqDataSetFromMatrix(psb, colData = meta, design = ~ group2)
dds <- estimateSizeFactors(dds)
cds <- estimateDispersions(dds)
vsd <- varianceStabilizingTransformation(cds)
plotPCA(vsd, intgroup = c("palliation"))

## Get PCA scree plot
## calculate the variance for each gene
rv <- rowVars(assay(vsd))
## select the ntop genes by variance
select <- order(rv, decreasing=TRUE)[seq_len(min(500, length(rv)))]
## perform a PCA on the data in assay(x) for the selected genes
pca <- prcomp(t(assay(vsd)[select, ]))
## the contribution to the total variance for each component
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
## plot the "percentVar"
scree_plot <- data.frame(percentVar[1:10])
scree_plot[, 2]<- c(1:10)

colnames(scree_plot) <- c("variance", 'component_number')

## Get covariates correlation with top PCs
res <- DEGreport::degCovariates(log2(counts(cds)+0.5), metadata = colData(cds), minPC = 2)
df <- res$corMatrix
df <- df[df$covar %in% c('age', 'sex', 'palliation', 'group2'), ]
df$covar[df$covar == 'group2'] <- 'disease'
df$logP <- -log10(df$pvalue)
df$covar <- factor(df$covar, levels = c('disease', 'palliation', 'age', 'sex'))
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
saveRDS(cds, 'analysis/PART30.cm_pseudobulk.deseq2.rds')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
