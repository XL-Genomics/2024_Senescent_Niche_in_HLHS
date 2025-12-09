####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Pediatric HLHS Compendium
####  2022-09-01 by Xiao LI (Texas Heart Institute, US)
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Initiate  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
Ver <- '2'
Step <- 'PART03_Filter_Cells'
Project <- '2022_hlhs_dturaga'
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Load data  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
merged.flt.srt <- readRDS('integrated/PART02.merged.cbn.srt.rds')
studies <- levels(merged.flt.srt$study)
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Set hard quality filters  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
## Filter cells with less than 200 umi(counts)/cell
cutoff1 <- 200
pct_filtered <- rep(NA, L(studies))
names(pct_filtered) <- studies

for(i in 1:L(studies)){
        n_filtered <- sum(merged.flt.srt$nCount_CBN < cutoff1 & merged.flt.srt$study == studies[i])
        pct_filtered[i] <- n_filtered*100/sum(merged.flt.srt$study == studies[i])
}
## Filter cells with less than 150 genes/cell
cutoff2 <- 150
pct_filtered <- rep(NA, L(studies))
names(pct_filtered) <- studies

for(i in 1:L(studies)){
        n_filtered <- sum(merged.flt.srt$nFeature_CBN < cutoff2 & merged.flt.srt$study == studies[i])
        pct_filtered[i] <- n_filtered*100/sum(merged.flt.srt$study == studies[i])
}
## Filter cells with higher than 5% mitochondrial counts/cell
cutoff3 <- 5
pct_filtered <- rep(NA, L(studies))
names(pct_filtered) <- studies

for(i in 1:L(studies)){
        n_filtered <- sum(merged.flt.srt$pct_mito_CBN > cutoff3 & merged.flt.srt$study == studies[i])
        pct_filtered[i] <- n_filtered*100/sum(merged.flt.srt$study == studies[i])
}

cells_toss_1 <- Cells(merged.flt.srt)[
        merged.flt.srt$nCount_CBN < cutoff1 |
                merged.flt.srt$nFeature_CBN < cutoff2 |
                merged.flt.srt$pct_mito_CBN > cutoff3
]
L(cells_toss_1) ## 34556
L(cells_toss_1)*100/ncol(merged.flt.srt) ## 7.404265%
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Set dynamic quality filters  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
L(cells_toss_1) ## 34556

## Remove low quality cells that did not pass the hard filters before setting dynamic filters
merged.flt.srt2 <- merged.flt.srt[, !Cells(merged.flt.srt) %in% cells_toss_1]
all_samples <- U(merged.flt.srt2$sample)
L(all_samples) ## 39

cutoffs <- data.frame(name = all_samples)

## Set filter for max nCount
cutoff_upper <- c()
cutoff_lower <- c()
cells_toss_2 <- c()
for(i in 1:L(all_samples)){
        print(i)
        cells_in_group <- Cells(merged.flt.srt2)[merged.flt.srt2$sample == all_samples[i]]
        GetOutlier(merged.flt.srt2$nCount_CBN[cells_in_group])
        cutoff_upper[i] <- GetOutlier(merged.flt.srt2$nCount_CBN[cells_in_group])[2]
        cutoff_lower[i] <- GetOutlier(merged.flt.srt2$nCount_CBN[cells_in_group])[1]
        if(cutoff_lower[i] < cutoff1){cutoff_lower[i] <- cutoff1}
        cells_toss_2 <- c(cells_toss_2,
                          cells_in_group[merged.flt.srt2$nCount_CBN[cells_in_group] > cutoff_upper[i] |
                                                 merged.flt.srt2$nCount_CBN[cells_in_group] < cutoff_lower[i]])
}
L(cells_toss_2) ## 24012

cutoffs$'nCount_max' <- cutoff_upper
cutoffs$'nCount_min' <- cutoff_lower

## Set filter for max nFeature
cutoff_upper <- c()
cutoff_lower <- c()
cells_toss_3 <- c()
for(i in 1:L(all_samples)){
        print(i)
        cells_in_group <- Cells(merged.flt.srt2)[merged.flt.srt2$sample == all_samples[i]]
        GetOutlier(merged.flt.srt2$nCount_CBN[cells_in_group])
        cutoff_upper[i] <- GetOutlier(merged.flt.srt2$nFeature_CBN[cells_in_group])[2]
        cutoff_lower[i] <- GetOutlier(merged.flt.srt2$nFeature_CBN[cells_in_group])[1]
        if(cutoff_lower[i] < cutoff2){cutoff_lower[i] <- cutoff2}
        cells_toss_3 <- c(cells_toss_3,
                          cells_in_group[merged.flt.srt2$nFeature_CBN[cells_in_group] > cutoff_upper[i] |
                                                 merged.flt.srt2$nFeature_CBN[cells_in_group] < cutoff_lower[i]])
}
L(cells_toss_3) ## 11578


cutoffs$'nFeature_max' <- cutoff_upper
cutoffs$'nFeature_min' <- cutoff_lower
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Filter cells  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
cells_toss <- U(c(cells_toss_1, cells_toss_2, cells_toss_3))
L(cells_toss) ## 59059

merged.flt.srt$LowQual <- F
merged.flt.srt$LowQual[cells_toss_3] <- 'cells_toss_3'
merged.flt.srt$LowQual[cells_toss_2] <- 'cells_toss_2'
merged.flt.srt$LowQual[cells_toss_1] <- 'cells_toss_1'

cell_toss.list <- list(cells_toss_1, cells_toss_2, cells_toss_3)

rm(merged.flt.srt2, p, p1, p2, p3, p4)
gc()
merged.flt.srt <- merged.flt.srt[, !Cells(merged.flt.srt) %in% cells_toss]
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
saveRDS(merged.flt.srt, 'integrated/PART03.merged.flt.srt.rds')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
