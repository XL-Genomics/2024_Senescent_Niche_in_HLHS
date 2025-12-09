####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Pediatric HLHS Compendium
####  2022-09-01 by Xiao LI (Texas Heart Institute, US)
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Initiate  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Ver <- '2'
Step <- 'PART65_Xenium_Integration'
Project <- '2022_hlhs_dturaga'
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Load snRNA Data  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### All cell types, all samples -- full.srt  ####
full.srt <- readRDS('integrated/PART19.annotated.srt.rds')
#### All cell types, all samples, no ambiguous -- full_clean.srt  ####
full_clean.srt <- DropMetaLevels(full.srt[, full.srt$Non_ambiguous])

cs_mk <- readRDS('analysis/PART28.cell_state_marker.srt_mk.rds')
ct_mk <- full.srt@misc$marker$Cell_type_marker
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Load Xenium Data  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ReadXenium2 <- function (data.dir, outs = c("matrix", "microns"), type = "centroids",
                        mols.qv.threshold = 20) {
    type <- match.arg(arg = type, choices = c("centroids", "segmentations"),
                      several.ok = TRUE)
    outs <- match.arg(arg = outs, choices = c("matrix", "microns"),
                      several.ok = TRUE)
    outs <- c(outs, type)
    has_dt <- requireNamespace("data.table", quietly = TRUE) &&
        requireNamespace("R.utils", quietly = TRUE)
    data <- sapply(outs, function(otype) {
        switch(EXPR = otype, matrix = {
            matrix <- suppressWarnings(Read10X(data.dir = file.path(data.dir,
                                                                    "cell_feature_matrix/")))
            matrix
        }, centroids = {
            if (has_dt) {
                cell_info <- as.data.frame(data.table::fread(file.path(data.dir,
                                                                       "cells.csv.gz")))
            } else {
                cell_info <- read.csv(file.path(data.dir, "cells.csv.gz"))
            }
            cell_centroid_df <- data.frame(x = cell_info$x_centroid,
                                           y = cell_info$y_centroid, cell = cell_info$cell_id,
                                           stringsAsFactors = FALSE)
            cell_centroid_df
        }, segmentations = {
            if (has_dt) {
                cell_boundaries_df <- as.data.frame(data.table::fread(file.path(data.dir,
                                                                                "cell_boundaries.csv.gz")))
            } else {
                cell_boundaries_df <- read.csv(file.path(data.dir,
                                                         "cell_boundaries.csv.gz"), stringsAsFactors = FALSE)
            }
            names(cell_boundaries_df) <- c("cell", "x", "y")
            cell_boundaries_df
        }, microns = {

            transcripts <- arrow::read_parquet(file.path(data.dir, "transcripts.parquet"))
            transcripts <- subset(transcripts, qv >= mols.qv.threshold)

            df <- data.frame(x = transcripts$x_location, y = transcripts$y_location,
                             gene = transcripts$feature_name, stringsAsFactors = FALSE)
            df
        }, stop("Unknown Xenium input type: ", otype))
    }, USE.NAMES = TRUE)
    return(data)
}

path <- '/Volumes/dale/data/xenium/2024_HLHS_EMB/K2Bio_10_2_24/'

## UK2_3B62D ####
data <- ReadXenium2(paste0(path, 'output-XETG00106__0004672__3B-62D__20240926__181354/'))
segmentations.data <- list(centroids = CreateCentroids(data$centroids),
                           segmentation = CreateSegmentation(data$microns))
coords <- CreateFOV(coords = segmentations.data, type = c("segmentation", "centroids"),
                    molecules = data$microns,
                    assay = "Xenium")
Ctrl_01 <- CreateSeuratObject(counts = data$matrix[["Gene Expression"]], assay = "Xenium")
Ctrl_01[["BlankCodeword"]] <- CreateAssayObject(counts = data$matrix[["Unassigned Codeword"]])
Ctrl_01[["ControlCodeword"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Codeword"]])
Ctrl_01[["ControlProbe"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Probe"]])
Ctrl_01[["Ctrl_01"]] <- coords
Ctrl_01$donor <- 'UK2_3B62D'
Ctrl_01$group1 <- 'Ctrl_01'
Ctrl_01$group2 <- 'Control'
saveRDS(Ctrl_01, 'individual/PART65.xenium_Ctrl_01_3B62D.raw.srt.rds')


## P175 ####
data <- ReadXenium2(paste0(path, 'output-XETG00106__0004672__P175-RV__20240926__181354/'))
segmentations.data <- list(centroids = CreateCentroids(data$centroids),
                           segmentation = CreateSegmentation(data$microns))
coords <- CreateFOV(coords = segmentations.data, type = c("segmentation", "centroids"),
                    molecules = data$microns,
                    assay = "Xenium")
HLHS_12 <- CreateSeuratObject(counts = data$matrix[["Gene Expression"]], assay = "Xenium")
HLHS_12[["BlankCodeword"]] <- CreateAssayObject(counts = data$matrix[["Unassigned Codeword"]])
HLHS_12[["ControlCodeword"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Codeword"]])
HLHS_12[["ControlProbe"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Probe"]])
HLHS_12[["HLHS_12"]] <- coords
HLHS_12$donor <- 'P175'
HLHS_12$group1 <- 'HLHS_12'
HLHS_12$group2 <- 'HLHS'
saveRDS(HLHS_12, 'individual/PART65.xenium_HLHS_12_P175.raw.srt.rds')

## P157 ####
data <- ReadXenium2(paste0(path, 'output-XETG00106__0004175__sample-4__20240926__181354/'))
segmentations.data <- list(centroids = CreateCentroids(data$centroids),
                           segmentation = CreateSegmentation(data$microns))
coords <- CreateFOV(coords = segmentations.data, type = c("segmentation", "centroids"),
                    molecules = data$microns,
                    assay = "Xenium")
HLHS_15 <- CreateSeuratObject(counts = data$matrix[["Gene Expression"]], assay = "Xenium")
HLHS_15[["BlankCodeword"]] <- CreateAssayObject(counts = data$matrix[["Unassigned Codeword"]])
HLHS_15[["ControlCodeword"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Codeword"]])
HLHS_15[["ControlProbe"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Probe"]])
HLHS_15[["HLHS_15"]] <- coords
HLHS_15$donor <- 'P157'
HLHS_15$group1 <- 'HLHS_15'
HLHS_15$group2 <- 'HLHS'
saveRDS(HLHS_15, 'individual/PART65.xenium_HLHS_15_P157.raw.srt.rds')



## P93 ####
data <- ReadXenium2(paste0(path, 'output-XETG00106__0004175__sample-2__20240926__181354/'))
segmentations.data <- list(centroids = CreateCentroids(data$centroids),
                           segmentation = CreateSegmentation(data$microns))
coords <- CreateFOV(coords = segmentations.data, type = c("segmentation", "centroids"),
                    molecules = data$microns,
                    assay = "Xenium")
HLHS_16 <- CreateSeuratObject(counts = data$matrix[["Gene Expression"]], assay = "Xenium")
HLHS_16[["BlankCodeword"]] <- CreateAssayObject(counts = data$matrix[["Unassigned Codeword"]])
HLHS_16[["ControlCodeword"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Codeword"]])
HLHS_16[["ControlProbe"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Probe"]])
HLHS_16[["HLHS_16"]] <- coords
HLHS_16$donor <- 'P93'
HLHS_16$group1 <- 'HLHS_16'
HLHS_16$group2 <- 'HLHS'
saveRDS(HLHS_16, 'individual/PART65.xenium_HLHS_16_P93.raw.srt.rds')

## P197 ####
data <- ReadXenium2(paste0(path, 'output-XETG00106__0004672__P197-new__20240926__181354/'))
segmentations.data <- list(centroids = CreateCentroids(data$centroids),
                           segmentation = CreateSegmentation(data$microns))
coords <- CreateFOV(coords = segmentations.data, type = c("segmentation", "centroids"),
                    molecules = data$microns,
                    assay = "Xenium")
VAD_04 <- CreateSeuratObject(counts = data$matrix[["Gene Expression"]], assay = "Xenium")
VAD_04[["BlankCodeword"]] <- CreateAssayObject(counts = data$matrix[["Unassigned Codeword"]])
VAD_04[["ControlCodeword"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Codeword"]])
VAD_04[["ControlProbe"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Probe"]])
VAD_04[["VAD_04"]] <- coords
VAD_04$donor <- 'P197'
VAD_04$group1 <- 'VAD_04'
VAD_04$group2 <- 'VAD'
saveRDS(VAD_04, 'individual/PART65.xenium_VAD_04_P197.raw.srt.rds')

## P125 ####
data <- ReadXenium2(paste0(path, 'output-XETG00106__0004175__sample-3__20240926__181354/'))
segmentations.data <- list(centroids = CreateCentroids(data$centroids),
                           segmentation = CreateSegmentation(data$microns))
coords <- CreateFOV(coords = segmentations.data, type = c("segmentation", "centroids"),
                    molecules = data$microns,
                    assay = "Xenium")
VAD_05 <- CreateSeuratObject(counts = data$matrix[["Gene Expression"]], assay = "Xenium")
VAD_05[["BlankCodeword"]] <- CreateAssayObject(counts = data$matrix[["Unassigned Codeword"]])
VAD_05[["ControlCodeword"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Codeword"]])
VAD_05[["ControlProbe"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Probe"]])
VAD_05[["VAD_05"]] <- coords
VAD_05$donor <- 'P125'
VAD_05$group1 <- 'VAD_05'
VAD_05$group2 <- 'VAD'
saveRDS(VAD_05, 'individual/PART65.xenium_VAD_05_P125.raw.srt.rds')



## Merge:
Ctrl_01 <- readRDS('individual/PART65.xenium_Ctrl_01_3B62D.raw.srt.rds')
HLHS_12 <- readRDS('individual/PART65.xenium_HLHS_12_P175.raw.srt.rds')
HLHS_15 <- readRDS('individual/PART65.xenium_HLHS_15_P157.raw.srt.rds')
HLHS_16 <- readRDS('individual/PART65.xenium_HLHS_16_P93.raw.srt.rds')
VAD_04 <- readRDS('individual/PART65.xenium_VAD_04_P197.raw.srt.rds')
VAD_05 <- readRDS('individual/PART65.xenium_VAD_05_P125.raw.srt.rds')

merged.srt <- merge(Ctrl_01, list(HLHS_12, HLHS_15, HLHS_16, VAD_04, VAD_05))
merged.srt <- NormalizeData(merged.srt)
merged.srt <- JoinLayers(merged.srt)

merged.srt$Sample <- revalue(merged.srt$group1, replace = c('Ctrl_01' = 'Xenium_01',
                                                            'HLHS_12' = 'Xenium_02',
                                                            'HLHS_15' = 'Xenium_03',
                                                            'HLHS_16' = 'Xenium_04',
                                                            'VAD_04' = 'Xenium_05',
                                                            'VAD_05' = 'Xenium_06'))
merged.srt$Donor <- merged.srt$donor
merged.srt$Group1 <- merged.srt$group1
merged.srt$Group2 <- merged.srt$group2
merged.srt$orig.ident <- NULL
merged.srt$group1 <- NULL
merged.srt$group2 <- NULL
merged.srt$donor <- NULL

merged.srt$Sex <- revalue(merged.srt$Donor, replace = c(
    'UK2_3B62D' = 'F',
    'P175' = 'M',
    'P157' = 'M',
    'P93' = 'M',
    'P197' = 'M',
    'P125' = 'M'
))

merged.srt$Age <- revalue(merged.srt$Donor, replace = c(
    'UK2_3B62D' = 3.78,
    'P175' = 7,
    'P157' = 11,
    'P93' = 15,
    'P197' = 8,
    'P125' = 17
))
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
saveRDS(merged.srt, 'integrated/PART65.xenium_merged.raw.srt.rds')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Xenium QC and Preprocessing  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Set hard quality filters  ####
merged.srt$LowQual <- ifelse(merged.srt$nCount_Xenium < 10 |
                                 merged.srt$nCount_Xenium > 2e3 |
                                 merged.srt$nFeature_Xenium > 200 |
                                 merged.srt$nCount_Xenium < 10, yes = T, no = F)
merged.srt$LowQual <- factor(merged.srt$LowQual, levels = c(F, T))

## To avoid bug when subsetting, we create a Seurat obj without the images.
srt <- CreateSeuratObject(counts = GetAssayData(merged.srt, assay = 'Xenium', layer = 'count'),
                          assay = 'Xenium',
                          meta.data = merged.srt@meta.data)
srt <- srt[, srt$LowQual == F]
srt <- NormalizeData(srt)


## Identify doublets
library('reticulate')
scr <- import('scrublet')
plt <- import("matplotlib.pyplot")

GetDoublet <- function(srt_obj, doublet_rate, dimN.var.toal){
    ## Scrublet (run via reticulate)
    mtx <- GetAssayData(srt_obj, layer = "counts", assay = 'Xenium')
    mtx <- t(mtx)
    scrub_model <- scr$Scrublet(mtx, expected_doublet_rate = doublet_rate)
    rst <- scrub_model$scrub_doublets(min_gene_variability_pctl = dimN.var.toal*100,
                                      n_prin_comps = 30L,
                                      min_counts = 2, min_cells = 3)
    rst[[2]] <- scrub_model$call_doublets(threshold = 0.3) ## adjusted based on histogram
    sc_doublets <- Cells(srt_obj)[rst[[2]]]
    sc_singlets <- Cells(srt_obj)[!rst[[2]]]
    srt_obj$Scrublet_doublet <- 'Singlet'
    srt_obj$Scrublet_doublet[rst[[2]]] <- 'Doublet'
    Scrublet <- rst[[1]]
    names(Scrublet) <- Cells(srt_obj)

    # p2 <- DimPlotSplit(srt_obj, split_by = 'Scrublet_doublet', split_order = c('Singlet', 'Doublet'),
    #                    cols.highlight = mycol_14[c(2, 1)], ncol = 2)
    # p2[[1]] <- p2[[1]] + labs(title = paste0('Srub Singlet: ', L(sc_singlets), ' Cells'))
    # p2[[2]] <- p2[[2]] + labs(title = paste0('Srub Doublet: ', L(sc_doublets), ' Cells'))
    # p <- wrap_plots(
    #     p2[[1]],
    #     p2[[2]],
    #     ncol = 2)
    return(list(
        sc_doublets,
        # p,
        Scrublet,
        scrub_model
    ))
}

doublet_rate <- 0.15 ## Assuming 15% doublet formation rate same as snRNA-seq
all_Doublet_SC <- c()
all_Scrublet <- c()
all_samples <- U(srt$Group1)
for(i in 1:L(all_samples)){
    gc()
    message(paste0('Processing ', all_samples[i], ' ...'))
    tmp.srt <- srt[, srt$Group1 == all_samples[i]]
    results <- GetDoublet(srt_obj = tmp.srt, doublet_rate = doublet_rate, dimN.var.toal = 0.8)
    all_Doublet_SC <- c(all_Doublet_SC, results[[1]])
    all_Scrublet <- c(all_Scrublet, results[[2]])
    ## plot scrublet histogram
    PlotPDF(paste0('2.', str_pad(i, pad = 0, width = 2), '.', all_samples[i], '.scrublet_hist'), 8, 4)
    print(plt$show(results[[3]]$plot_histogram()[[1]]))
    dev.off()
}

####  Evaluate doublets  ####
## save doublets to the main Seurat
srt$Doublet_SC <- F
srt$Doublet_SC[all_Doublet_SC] <- T
srt$Doublet_SC_score <- NA
srt$Doublet_SC_score[names(all_Scrublet)] <- all_Scrublet
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Xenium Multi-sample Integration  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Harmony-batch correction RNA integration  ####
srt$orig.ident <- NULL
genes_include <- rownames(srt@assays$Xenium)

####  FastMNN-batch correction RNA integration  ####
tmp.srt <- DietSeurat(srt, assays = 'Xenium', features = genes_include)
tmp.srt <- NormalizeData(tmp.srt)
tmp.srt[['Xenium']] <- as(tmp.srt[['Xenium']], Class = 'Assay') ## Switch to v3 assay

srt.list <- SplitObject(tmp.srt, split.by = 'Group1')
for(i in 1:L(srt.list)){
    srt.list[[i]] <- SCTransform(object = srt.list[[i]], method ="glmGamPoi", assay = 'Xenium',
                                 variable.features.n = 400)
}

tmp.srt <- RunFastMNN(object.list = srt.list, assay = 'SCT', features = 400)
tmp.srt <- RunUMAP(tmp.srt, reduction = 'mnn', dims = 1:50,
                   reduction.name = 'RNA_FastMNN_umap', reduction.key = 'mmnrnaumap_',
                   min.dist = 0.5, n.neighbors = 50)
tmp.srt$Group1 <- factor(tmp.srt$Group1, levels = levels(srt$Group1))

srt@assays$SCT <- tmp.srt@assays$SCT
srt@reductions$mnn <- tmp.srt@reductions$mnn
srt@reductions$RNA_FastMNN_umap <- tmp.srt@reductions$RNA_FastMNN_umap
srt@reductions$mnn@assay.used <- 'Xenium'
srt@reductions$RNA_FastMNN_umap@assay.used <- 'Xenium'


####  Clustering ####
srt <- FindNeighbors(srt, dims = 1:50, reduction = 'mnn')
srt <- FindClusters(srt, resolution = seq(0.1, 0.9, 0.1))
srt <- FindClusters(srt, resolution = c(1, 2))
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
saveRDS(srt, 'integrated/PART65.xenium_integrated_no_images.srt.rds')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----



####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Convert H5AD for scANVI and CellTypist  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
x <- RenameAssays(srt, 'Xenium', 'RNA')
x <- DietSeurat(x, dimreducs = c('RNA_FastMNN_umap', 'mnn'), misc = F, layers = c('counts', 'data'), assays = 'RNA')
x[['RNA']] <- as(object = x[['RNA']], Class = "Assay")

SaveH5ad(x, path = 'integrated/', name = 'PART65.xenium_integrated_no_images.ann',
         assay = 'RNA', raw_count_only = T, verbose = T)
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----



####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Cell type annotation  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Try scANVI  ####
## Load scANVI cell state label ####
pred.df <- read.csv('integrated/PART66.xenium_scanvi_integration.csvs/obs.csv')
pred.df <- pred.df[pred.df$batch == 0, ]
identical(str_split(pred.df$X, '-0', simplify = T)[,1], Cells(srt))

srt$Cell_type_scANVI <- pred.df$scANVI_predict_celltype
srt$Cell_state_scANVI <- pred.df$scANVI_predict_cellstate


#### Try CellTypist  ####
## Load CellTypist cell state label ####
pred.df <- read.csv('analysis/part53.celltypist_adata_obs.csv') ## csv from Yi
identical(pred.df$X, Cells(srt))

srt$Cell_state_celltypist <- pred.df$majority_voting


#### Try SingleR  ####
## Prepare reference
library('SingleR')
library('BiocParallel')
ref.srt <- DropMetaLevels(full.srt[, ! full.srt$Cell_type  %in% c('MitoHi', 'Doublet')])
ref.sce <- as.SingleCellExperiment(DietSeurat(ref.srt))

que.sce <- as.SingleCellExperiment(DietSeurat(srt))

system.time({
    pred_celltype <- SingleR(test = que.sce,
                         ref = ref.sce,
                         labels = factor(ref.sce$Cell_type),
                         aggr.ref = T,
                         BPPARAM = SnowParam(16))
}) ## 31.853
srt$Cell_type_SingleR <- pred_celltype$pruned.labels
srt$Cell_type_SingleR[is.na(srt$Cell_type_SingleR)] <- 'Unmapped'

system.time({
    pred_cellstate <- SingleR(test = que.sce,
                              ref = ref.sce,
                              labels = factor(ref.sce$Cell_state),
                              aggr.ref = T,
                              BPPARAM = SnowParam(16))
}) ## 42.205
srt$Cell_state_SingleR <- pred_cellstate$pruned.labels
srt$Cell_state_SingleR[is.na(srt$Cell_state_SingleR)] <- 'Unmapped'
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
saveRDS(srt, 'integrated/PART65.xenium_integrated_no_images.srt.rds') ## Update
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


#### Consolidate scANVI, CellTypist and SingleR annotations  ####
srt$Cell_state[is.na(srt$Cell_state)] <- as.vector(srt$Cell_type[is.na(srt$Cell_state)])
MappingHeatmap(srt, 'Cell_state', 'Cell_type')
srt$Cell_state <- factor(srt$Cell_state, levels = c('CM1', 'CM2', 'CM3',
                                                    'FB1', 'FB2', 'FB3', 'FB4', 'FB5',
                                                    'EC1', 'EC2', 'EC3', 'EC4', 'LEC', 'EndoC',
                                                    'PC1', 'PC2', 'SMC1', 'SMC2',
                                                    'Lym', 'Mono', 'MP1', 'MP2', 'MP3',
                                                    'NC', 'AC', 'Doublet'))
srt$Cell_state[is.na(srt$Cell_state)] <- 'Doublet'

all_meta <- srt@meta.data
all_reduction <- srt@reductions

srt@meta.data <- srt@meta.data[, c(
    'nCount_Xenium',
    'nFeature_Xenium',
    'nCount_BlankCodeword',
    'nFeature_BlankCodeword',
    'nCount_ControlCodeword',
    'nFeature_ControlCodeword',
    'nCount_ControlProbe',
    'nFeature_ControlProbe',
    'Sample',
    'Donor',
    'Group1',
    'Group2',
    'Sex',
    'Age',
    'Doublet_SC',
    'Doublet_SC_score',
    'HLHS_Sen_Niche',
    'Cell_type',
    'Cell_state'
)]
srt@reductions$full_umap <- srt@reductions$RNA_FastMNN_umap
colnames(srt@reductions$full_umap@cell.embeddings) <- paste0('fullumap_', 1:2)
srt@reductions$full_umap@key <- 'fullumap_'
srt <- DietSeurat(srt, dimreducs = c('mnn', 'full_umap', 'sub_umap'), assays = 'Xenium')

srt_clean <- DropMetaLevels(srt[, srt$Cell_state != 'Doublet'])
srt_clean <- RunUMAP(srt_clean, dims = 1:50, reduction = 'mnn',  min.dist = 0.8,
                     reduction.name = 'clean_umap', reduction.key = 'cleanumap_')

srt@reductions$clean_umap <- srt@reductions$full_umap
srt@reductions$clean_umap@cell.embeddings[, 1:2] <- NA
srt@reductions$clean_umap@cell.embeddings[Cells(srt_clean), 1:2] <-
    srt_clean@reductions$clean_umap@cell.embeddings[, 1:2]
colnames(srt@reductions$clean_umap@cell.embeddings) <- paste0('cleanumap_', 1:2)
srt@reductions$clean_umap@key <- 'cleanumap_'

####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
saveRDS(all_meta, 'integrated/PART65.xenium_annotated_no_image.srt_meta.rds')
saveRDS(all_reduction, 'integrated/PART65.xenium_annotated_no_image.srt_dr.rds')
saveRDS(srt, 'integrated/PART65.xenium_annotated_no_image.srt.rds')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Export Cell type label for Xenium explorer:
for(i in 1:6){
    df <- srt@meta.data[srt$Donor == U(srt$Donor)[i], ]
    df$cell_id <- str_split(rownames(df), '_', simplify = T)[, 1]
    df$group <- df$Cell_type
    WriteCSV(df[, c('cell_id', 'group')], title = paste0(U(srt$Donor)[i], '_Cell_type'))
    df$group <- df$Cell_state
    WriteCSV(df[, c('cell_id', 'group')], title = paste0(U(srt$Donor)[i], '_Cell_state'))
}
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
