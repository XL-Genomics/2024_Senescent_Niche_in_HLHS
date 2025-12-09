####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Pediatric HLHS Compendium
####  2022-09-01 by Xiao LI (Texas Heart Institute, US)
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Initiate  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
Ver <- '2'
Step <- 'PART04_Doublet'
Project <- '2022_hlhs_dturaga'
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Load sample metadata  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
merged.dlt.srt <- readRDS('integrated/PART03.merged.flt.srt.rds')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Global Functions  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
GetDoublet <- function(srt_obj, doublet_rate, dimN.var.toal){
        ## Scrublet (run via reticulate)
        mtx <- srt_obj@assays$CBN@counts
        mtx <- t(mtx)
        scrub_model <- scr$Scrublet(mtx, expected_doublet_rate = doublet_rate)
        rst <- scrub_model$scrub_doublets(min_gene_variability_pctl = dimN.var.toal*100,
                                          n_prin_comps = 30L,
                                          min_counts = 2, min_cells = 3)
        rst[[2]] <- scrub_model$call_doublets(threshold = 0.25) ## adjusted based on histogram
        sc_doublets <- Cells(srt_obj)[rst[[2]]]
        sc_singlets <- Cells(srt_obj)[!rst[[2]]]
        srt_obj$Scrublet_doublet <- 'Singlet'
        srt_obj$Scrublet_doublet[rst[[2]]] <- 'Doublet'
        Scrublet <- rst[[1]]
        names(Scrublet) <- Cells(srt_obj)

        p2 <- DimPlotSplit(srt_obj, split_by = 'Scrublet_doublet', split_order = c('Singlet', 'Doublet'),
                           cols.highlight = mycol_14[c(2, 1)], ncol = 2)
        p2[[1]] <- p2[[1]] + labs(title = paste0('Srub Singlet: ', L(sc_singlets), ' Cells'))
        p2[[2]] <- p2[[2]] + labs(title = paste0('Srub Doublet: ', L(sc_doublets), ' Cells'))
        p <- wrap_plots(
                p2[[1]],
                p2[[2]],
                ncol = 2)
        return(list(
                sc_doublets,
                p,
                Scrublet,
                scrub_model
        ))
}
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
###  Identify doublets for each dataset  (linear processing) ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
doublet_rate <- 0.15 ## Assuming 15% doublet formation rate based on 10X estimation
all_Doublet_SC <- c()
all_Scrublet <- c()

all_samples <- levels(merged.dlt.srt$sample)
for(i in 1:L(all_samples)){
        gc()
        message(paste0('Processing ', all_samples[i], ' ...'))
        tmp.srt <- merged.dlt.srt[, merged.dlt.srt$sample == all_samples[i]]
        results <- GetDoublet(srt_obj = tmp.srt, doublet_rate = doublet_rate, dimN.var.toal = 0.85)
        all_Doublet_SC <- c(all_Doublet_SC, results[[1]])
        ## plot umap
        PlotPDF(paste0('1.', str_pad(i, pad = 0, width = 2), '.', all_samples[i], '.doublets_found'), 10, 5)
        print(results[[2]])
        dev.off()
        all_Scrublet <- c(all_Scrublet, results[[3]])
        ## plot scrublet histogram
        PlotPDF(paste0('1.', str_pad(i, pad = 0, width = 2), '.', all_samples[i], '.scrublet_hist'), 8, 4)
        print(plt$show(results[[4]]$plot_histogram()[[1]]))
        dev.off()
}

####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Evaluate doublets  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
## save doublets to the main Seurat
merged.dlt.srt$Doublet_SC <- F
merged.dlt.srt$Doublet_SC[all_Doublet_SC] <- T
merged.dlt.srt$Doublet_SC_score <- NA
merged.dlt.srt$Doublet_SC_score[names(all_Scrublet)] <- all_Scrublet

####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
saveRDS(merged.dlt.srt@meta.data, 'integrated/PART04.merged.dlt.srt_meta.rds')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
