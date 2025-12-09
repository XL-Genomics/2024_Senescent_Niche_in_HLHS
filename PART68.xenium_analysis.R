####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Pediatric HLHS Compendium
####  2022-09-01 by Xiao LI (Texas Heart Institute, US)
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Initiate  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Ver <- '2'
Step <- 'PART68_Xenium_Analysis'
Project <- '2022_hlhs_dturaga'
Color_palliation <- c('grey85', '#D02E30FF', '#FBB157', '#6FB9A2FF', '#65A8DC')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Load snRNA Data  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
xen.srt <- readRDS('integrated/PART65.xenium_annotated_no_image.srt.rds')
xen_clean.srt <- DropMetaLevels(xen.srt[, xen.srt$Cell_state != 'Doublet'])
xen_clean.srt$Sen_cell <- ifelse(xen_clean.srt$Cell_state %in% c('CM3', 'FB5', 'EC4', 'MP3', 'PC2'), 'Sen', 'NonSen')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Spatial NN  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
xen_image.srt <- readRDS('integrated/PART65.xenium_annotated_with_image.srt.rds')
max_dist <- 100

## Use the Squidpy algo for spatial NN clustering:
# sq.gr.spatial_neighbors(
#         adata,
#         coord_type="generic",     # (use "visium" for Visium)
#         n_neighs=50,                # or set 'radius=...' (in the same units as coordinates)
#         delaunay=False
# )
# from scipy import sparse
# import numpy as np
#
# # Get adjacency
# A = adata.obsp["spatial_connectivities"].tocsr()
#
# # One-hot encoding of cell types
# ct = adata.obs["Cell_state"].astype("category")
# M = sparse.csr_matrix(
#         (np.ones(ct.size), (np.arange(ct.size), ct.cat.codes.values)),
#         shape=(ct.size, len(ct.cat.categories))
# )
#
# # Neighborhood composition per cell = A * M  (counts of each type in its neighborhood)
# N = A @ M                                # shape: n_cells x n_cell_types
# N = N.toarray()
# # Row-normalize to proportions (optional)
# N = N / (N.sum(axis=1, keepdims=True) + 1e-9)
#
# # Store as obsm for downstream clustering/UMAP
# adata.obsm["niche_comp"] = N

## R version:
DoSpatialClustering <- function(obj, res = 0.8){
    # Compute pairwise Euclidean distances
    dist_mat <- as.matrix(dist(obj@reductions$spatial@cell.embeddings))
    # Build adjacency (1 if within radius)
    adj <- (dist_mat <= max_dist) * 1
    diag(adj) <- 0
    # Convert to sparse matrix
    adj <- Matrix(adj, sparse = TRUE)
    # Create a Seurat graph
    obj@graphs$spat <- adj
    # Assume we already built the spatial graph in (A) with graph.name = "spat"
    A <- obj@graphs$spat              # sparse adjacency (spots x spots)
    A <- (A > 0) + 0                      # binarize
    # Provide/compute discrete labels (e.g., cell_type or deconvolved cell class)
    lab <- factor(obj$Cell_state)          # or any discrete annotation in obj@meta.data
    # One-hot of labels (spots x levels)
    M <- sparse.model.matrix(~ lab - 1)   # columns = levels(lab)
    # Neighborhood composition = A %*% M  (counts of each label in each spot’s neighborhood)
    N <- A %*% M
    # Convert to proportions per spot
    rs <- Matrix::rowSums(N); rs[rs == 0] <- 1
    N <- N / rs
    # Embed & cluster in “niche space”
    npc <- 20
    pcs <- prcomp(as.matrix(N), center = TRUE, scale. = TRUE)$x[, 1:min(npc, ncol(N)), drop = FALSE]
    obj[["niche"]] <- CreateDimReducObject(embeddings = as.matrix(pcs), key = "niche_", assay = DefaultAssay(obj))

    obj <- FindNeighbors(obj, reduction = "niche", dims = 1:ncol(pcs), k.param = 50, graph.name = "niche_nn")
    obj <- FindClusters(obj, graph.name = "niche_nn", resolution = res)
    obj$niche_label <- Idents(obj)
    return(obj)
}

## HLHS_15 ####
obj <- DropMetaLevels(xen_clean.srt[, xen_clean.srt$Group1 == 'HLHS_15'])
obj@reductions$spatial <- obj@reductions$clean_umap
coor <- xen_image.srt@images$HLHS_15$centroids@coords[, 1:2]
rownames(coor) <- Cells(xen_image.srt)[xen_image.srt$Group1 == 'HLHS_15']
obj@reductions$spatial@cell.embeddings <- coor[Cells(obj), 1:2]
colnames(obj@reductions$spatial@cell.embeddings) <- c('spatial_1', 'spatial_2')
obj@reductions$spatial@key <- 'spatial_'
obj <- UpdateSeuratObject(obj)
obj <- DoSpatialClustering(obj)
# Visualize
p1 <- DimPlot2(obj, group.by = "niche_label",  reduction = 'spatial')
p2 <- MetaMatchingHeatmap(obj,  que_var = 'Cell_state', ref_var = 'niche_label', percentage = T)
p3 <- MetaMatchingHeatmap(obj,  ref_var = 'Cell_state', que_var = 'niche_label', percentage = T)
PlotPDF('3.1.xenium_HLHS_15_spatial_niche_umap_cell_state_enrichment', 24, 8)
p1 + p2 + p3
dev.off()
data <- obj@meta.data[, c('Group1', 'Cell_type', 'Cell_state',  'Coord_x',  'Coord_y', 'Sen_cell', 'niche_label')]
saveRDS(data, 'analysis/PART68.xenium_niche_clustering_HLHS_15')

## HLHS_16 ####
obj <- DropMetaLevels(xen_clean.srt[, xen_clean.srt$Group1 == 'HLHS_16'])
obj@reductions$spatial <- obj@reductions$clean_umap
coor <- xen_image.srt@images$HLHS_16$centroids@coords[, 1:2]
rownames(coor) <- Cells(xen_image.srt)[xen_image.srt$Group1 == 'HLHS_16']
obj@reductions$spatial@cell.embeddings <- coor[Cells(obj), 1:2]
colnames(obj@reductions$spatial@cell.embeddings) <- c('spatial_1', 'spatial_2')
obj@reductions$spatial@key <- 'spatial_'
obj <- UpdateSeuratObject(obj)
obj <- DoSpatialClustering(obj)
# Visualize
p1 <- DimPlot2(obj, group.by = "niche_label",  reduction = 'spatial')
p2 <- MetaMatchingHeatmap(obj,  que_var = 'Cell_state', ref_var = 'niche_label', percentage = T)
p3 <- MetaMatchingHeatmap(obj,  ref_var = 'Cell_state', que_var = 'niche_label', percentage = T)
PlotPDF('3.2.xenium_HLHS_16_spatial_niche_umap_cell_state_enrichment', 24, 8)
p1 + p2 + p3
dev.off()
data <- obj@meta.data[, c('Group1', 'Cell_type', 'Cell_state',  'Coord_x',  'Coord_y', 'Sen_cell', 'niche_label')]
saveRDS(data, 'analysis/PART68.xenium_niche_clustering_HLHS_16')


## HLHS_12 ####
obj <- DropMetaLevels(xen_clean.srt[, xen_clean.srt$Group1 == 'HLHS_12'])
obj@reductions$spatial <- obj@reductions$clean_umap
coor <- xen_image.srt@images$HLHS_12$centroids@coords[, 1:2]
rownames(coor) <- Cells(xen_image.srt)[xen_image.srt$Group1 == 'HLHS_12']
obj@reductions$spatial@cell.embeddings <- coor[Cells(obj), 1:2]
colnames(obj@reductions$spatial@cell.embeddings) <- c('spatial_1', 'spatial_2')
obj@reductions$spatial@key <- 'spatial_'
obj <- UpdateSeuratObject(obj)
obj <- DoSpatialClustering(obj, res = 0.6)
# Visualize
p1 <- DimPlot2(obj, group.by = "niche_label",  reduction = 'spatial')
p2 <- MetaMatchingHeatmap(obj,  que_var = 'Cell_state', ref_var = 'niche_label', percentage = T)
p3 <- MetaMatchingHeatmap(obj,  ref_var = 'Cell_state', que_var = 'niche_label', percentage = T)
PlotPDF('3.3.xenium_HLHS_12_spatial_niche_umap_cell_state_enrichment', 24, 8)
p1 + p2 + p3
dev.off()
data <- obj@meta.data[, c('Group1', 'Cell_type', 'Cell_state',  'Coord_x',  'Coord_y', 'Sen_cell', 'niche_label')]
saveRDS(data, 'analysis/PART68.xenium_niche_clustering_HLHS_12')


## Donor Ctrl_01 ####
obj <- DropMetaLevels(xen_clean.srt[, xen_clean.srt$Group1 == 'Ctrl_01'])
obj@reductions$spatial <- obj@reductions$clean_umap
coor <- xen_image.srt@images$Ctrl_01$centroids@coords[, 1:2]
rownames(coor) <- Cells(xen_image.srt)[xen_image.srt$Group1 == 'Ctrl_01']
obj@reductions$spatial@cell.embeddings <- coor[Cells(obj), 1:2]
colnames(obj@reductions$spatial@cell.embeddings) <- c('spatial_1', 'spatial_2')
obj@reductions$spatial@key <- 'spatial_'
obj <- UpdateSeuratObject(obj)
obj <- DoSpatialClustering(obj, res = 0.5)
# Visualize
p1 <- DimPlot2(obj, group.by = "niche_label",  reduction = 'spatial')
p2 <- MetaMatchingHeatmap(obj,  que_var = 'Cell_state', ref_var = 'niche_label', percentage = T)
p3 <- MetaMatchingHeatmap(obj,  ref_var = 'Cell_state', que_var = 'niche_label', percentage = T)
PlotPDF('3.4.xenium_Ctrl_01_spatial_niche_umap_cell_state_enrichment', 24, 8)
p1 + p2 + p3
dev.off()
data <- obj@meta.data[, c('Group1', 'Cell_type', 'Cell_state',  'Coord_x',  'Coord_y', 'Sen_cell', 'niche_label')]
saveRDS(data, 'analysis/PART68.xenium_niche_clustering_Ctrl_01')


## VAD_04 ####
obj <- DropMetaLevels(xen_clean.srt[, xen_clean.srt$Group1 == 'VAD_04'])
obj@reductions$spatial <- obj@reductions$clean_umap
coor <- xen_image.srt@images$VAD_04$centroids@coords[, 1:2]
rownames(coor) <- Cells(xen_image.srt)[xen_image.srt$Group1 == 'VAD_04']
obj@reductions$spatial@cell.embeddings <- coor[Cells(obj), 1:2]
colnames(obj@reductions$spatial@cell.embeddings) <- c('spatial_1', 'spatial_2')
obj@reductions$spatial@key <- 'spatial_'
obj <- UpdateSeuratObject(obj)
obj <- DoSpatialClustering(obj, res = 0.5)
# Visualize
p1 <- DimPlot2(obj, group.by = "niche_label",  reduction = 'spatial')
p2 <- MetaMatchingHeatmap(obj,  que_var = 'Cell_state', ref_var = 'niche_label', percentage = T)
p3 <- MetaMatchingHeatmap(obj,  ref_var = 'Cell_state', que_var = 'niche_label', percentage = T)
PlotPDF('3.5.xenium_VAD_04_spatial_niche_umap_cell_state_enrichment', 24, 8)
p1 + p2 + p3
dev.off()
data <- obj@meta.data[, c('Group1', 'Cell_type', 'Cell_state',  'Coord_x',  'Coord_y', 'Sen_cell', 'niche_label')]
saveRDS(data, 'analysis/PART68.xenium_niche_clustering_VAD_04')


## VAD_05 ####
obj <- DropMetaLevels(xen_clean.srt[, xen_clean.srt$Group1 == 'VAD_05'])
obj@reductions$spatial <- obj@reductions$clean_umap
coor <- xen_image.srt@images$VAD_05$centroids@coords[, 1:2]
rownames(coor) <- Cells(xen_image.srt)[xen_image.srt$Group1 == 'VAD_05']
obj@reductions$spatial@cell.embeddings <- coor[Cells(obj), 1:2]
colnames(obj@reductions$spatial@cell.embeddings) <- c('spatial_1', 'spatial_2')
obj@reductions$spatial@key <- 'spatial_'
obj <- UpdateSeuratObject(obj)
obj <- DoSpatialClustering(obj, res = 0.5)
# Visualize
p1 <- DimPlot2(obj, group.by = "niche_label",  reduction = 'spatial')
p2 <- MetaMatchingHeatmap(obj,  que_var = 'Cell_state', ref_var = 'niche_label', percentage = T)
p3 <- MetaMatchingHeatmap(obj,  ref_var = 'Cell_state', que_var = 'niche_label', percentage = T)
PlotPDF('3.6.xenium_VAD_05_spatial_niche_umap_cell_state_enrichment', 24, 8)
p1 + p2 + p3
dev.off()
data <- obj@meta.data[, c('Group1', 'Cell_type', 'Cell_state',  'Coord_x',  'Coord_y', 'Sen_cell', 'niche_label')]
saveRDS(data, 'analysis/PART68.xenium_niche_clustering_VAD_05')
