####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Pediatric HLHS Compendium
####  2022-09-01 by Xiao LI (Texas Heart Institute, US)
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Initiate  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
Ver <- '2'
Step <- 'PART01_Data_Collection'
Project <- '2022_hlhs_dturaga'
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Global Functions  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
MakeSrt <- function(mode, matrix_dir, sample, study, method, platform, protocol, data_process, tissue, enrichment,
                    preparation, diagnosis, sex, age, donor, replicate) {
        if (mode == '10x') {
                matrix <- paste0(matrix_dir, '/outs/filtered_feature_bc_matrix/')
                srt <- CreateSeuratObject(counts = Read10X(data.dir = matrix),
                                          min.cells = 1, min.features = 1, project = study)
        }
        else if (mode == 'cellbender') {
                matrix <- paste0(matrix_dir, '/cellbender_filtered.h5')
                srt <- CreateSeuratObject(counts = ReadCB_h5(matrix), min.cells = 1, min.features = 1, project = study)
        }
        else if (mode == 'matrix') {
                matrix <- read.table(gzfile(paste0(matrix_dir[1], '.matrix.csv.gz')), header = T, sep = ',')
                srt <- CreateSeuratObject(counts = matrix, min.cells = 1, min.features = 1, project = study)
        }
        srt$sample <- sample
        srt$orig.name <- Cells(srt)
        srt$study <- study
        srt$method <- method
        srt$platform <- platform
        srt$protocol <- protocol
        srt$data_process <- data_process
        srt$tissue <- tissue
        srt$enrichment <- enrichment
        srt$preparation <- preparation
        srt$diagnosis <- diagnosis
        srt$sex <- sex
        srt$age <- age
        srt$donor <- donor
        srt$replicate <- replicate
        srt <- RenameCells(srt, new.names = paste(srt$study,
                                                  srt$sample,
                                                  srt$orig.name,
                                                  sep = ':'), for.merge = F)
        srt <- PercentageFeatureSet(srt, pattern = '^MT-', col.name = 'pct_mito', assay = 'RNA')
        srt$pct_mito[is.nan(srt$pct_mito)] <- 0
        return(srt)
}
MakeDataset <- function(study, study_id, sample_name, mode, matrix_dir, starting_sample = 1){
        srt.list <- list()
        for(i in 1:L(sample_name)) {
                sample_id = paste0(study_id, '_S', str_pad(starting_sample - 1 + i, 3, pad = '0'))
                message('Processing sample:', sample_id)
                sample_meta_sub.df <- sample_meta.df[sample_meta.df$Study == study &
                                                             sample_meta.df$Sample_id == sample_id &
                                                             sample_meta.df$Name_on_disk == sample_name[i], ]
                message('Sample metadata found')
                srt.list[[i]] <- MakeSrt(mode = mode,
                                         matrix_dir = matrix_dir[i],
                                         study = study,
                                         sample = sample_id,
                                         method = sample_meta_sub.df$Method,
                                         platform = sample_meta_sub.df$Platform,
                                         protocol = sample_meta_sub.df$Protocol,
                                         data_process = sample_meta_sub.df$Data_process,
                                         tissue = sample_meta_sub.df$Tissue,
                                         enrichment = sample_meta_sub.df$Enrichment,
                                         preparation = sample_meta_sub.df$Preparation,
                                         diagnosis = sample_meta_sub.df$Diagnosis,
                                         sex = sample_meta_sub.df$Sex,
                                         age = sample_meta_sub.df$Age,
                                         donor = sample_meta_sub.df$Donor,
                                         replicate = sample_meta_sub.df$Replicate
                )
                print('Seurat generated...')
                # print(srt.list[[i]])
                # cat('\n_____________________________________________________\n')
        }
        if(L(srt.list) > 1) {
                merge.srt <- merge(srt.list[[1]], srt.list[2:L(srt.list)])
        } else {
                merge.srt <- srt.list[[1]]
        }
        return(merge.srt)
}
MakeRawDataset <- function(study, sample_name, raw_matrix_type, starting_sample = 1){
        no. <- names(studies[studies==study])
        message('Collecting Raw Data...')
        merge.srt <- MakeDataset(study = study,
                                 study_id = no.,
                                 sample_name = sample_name,
                                 mode = raw_matrix_type,
                                 matrix_dir = paste0('/Volumes/shire/data/scrnaseq/',
                                                     study, '/matrix/', sample_name),
                                 starting_sample = starting_sample)
        message('Processing Raw Seurat Object...')
        # merge.srt <- Process(merge.srt, assay = 'RNA')
        message('Saving Raw Seurat Object...')
        saveRDS(merge.srt, paste0('individual/', no., '.', study, '.raw.srt.rds'))
        # SaveH5ad(merge.srt, path = 'individual/', name = paste0(no., '.', study, '.raw'),
        #          assay = 'RNA', raw_count_only = F, verbose = F)
        rm(merge.srt)
        gc()
}
MakeCbnDataset <- function(study, sample_name, raw_matrix_type, starting_sample = 1, cb_folder = 'cellbender_v1'){
        no. <- names(studies[studies==study])
        message('Collecting CellBender Data...')
        merge.srt <- MakeDataset(study = study,
                                 study_id = no.,
                                 sample_name = sample_name,
                                 mode = 'cellbender',
                                 matrix_dir = paste0('/Volumes/shire/data/scrnaseq/',
                                                     study, '/matrix/', cb_folder, '/', sample_name),
                                 starting_sample = starting_sample)
        message('Processing CellBender Seurat Object...')
        # merge.srt <- Process(merge.srt, assay = 'RNA')
        message('Saving CellBender Seurat Object...')
        saveRDS(merge.srt, paste0('individual/', no., '.', study, '.cbn.srt.rds'))
        # SaveH5ad(merge.srt, path = 'individual/', name = paste0(no., '.', study, '.cbn'),
        #          assay = 'RNA', raw_count_only = F, verbose = F)
        rm(merge.srt)
        gc()
}
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Load sample metadata  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
sample_meta.df <- read.csv(paste0(Docu_dir, 'pediatric_sample_meta.csv'))
studies <- U(sample_meta.df$Study)
names(studies) <- U(sample_meta.df$Study_id)
studies_cellbender <- studies[studies %in% sample_meta.df$Study[sample_meta.df$Platform %in% c('10X', 'Drop-seq')]]
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Dataset #1: 2022_Nature_JMartin (previously 2019_Unpub7_JMartin)  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
study <- '2022_Nature_JMartin'
sample_name <- c(
        '2019_Unpub7_JMartin_ctrl_s1',
        '2019_Unpub7_JMartin_ctrl_s2',
        '2019_Unpub7_JMartin_ctrl_s3',
        '2019_Unpub7_JMartin_ctrl_s8',
        '2019_Unpub7_JMartin_ctrl_s9',
        '2019_Unpub7_JMartin_ctrl_s10',
        '2019_Unpub7_JMartin_ctrl_s4',
        '2019_Unpub7_JMartin_ctrl_s5',
        '2019_Unpub7_JMartin_ctrl_s11',
        '2019_Unpub7_JMartin_ctrl_s12',
        '2019_Unpub7_JMartin_ctrl_s6',
        '2019_Unpub7_JMartin_ctrl_s7',
        '2019_Unpub7_JMartin_hlhs_s2',
        '2019_Unpub7_JMartin_hlhs_s3',
        '2019_Unpub7_JMartin_hlhs_s4',
        '2019_Unpub7_JMartin_htx_s1',
        '2019_Unpub7_JMartin_htx_s2'
        )
matrix_type <- '10x'
MakeRawDataset(study = study, sample_name = sample_name, raw_matrix_type = matrix_type)
MakeCbnDataset(study = study, sample_name = sample_name, raw_matrix_type = matrix_type, cb_folder = 'cellbender')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Dataset #2: 2022_Unpub3_DTuraga  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
study <- '2022_Unpub3_DTuraga'
sample_name <- c(
        '2022_Unpub3_DTuraga_P184',
        '2022_Unpub3_DTuraga_P179',
        '2022_Unpub3_DTuraga_P181'
)
matrix_type <- '10x'
MakeRawDataset(study = study, sample_name = sample_name, raw_matrix_type = matrix_type)
MakeCbnDataset(study = study, sample_name = sample_name, raw_matrix_type = matrix_type)
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Dataset #3: 2022_Unpub4_DTuraga  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
study <- '2022_Unpub4_DTuraga'
sample_name <- c(
        '2022_Unpub4_DTuraga_P81',
        '2022_Unpub4_DTuraga_P84'
)
matrix_type <- '10x'
MakeRawDataset(study = study, sample_name = sample_name, raw_matrix_type = matrix_type)
MakeCbnDataset(study = study, sample_name = sample_name, raw_matrix_type = matrix_type)
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Dataset #4: 2022_Unpub5_DTuraga  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
study <- '2022_Unpub5_DTuraga'
sample_name <- c(
       '2022_Unpub5_DTuraga_P70',
       '2022_Unpub5_DTuraga_P175',
       '2022_Unpub5_DTuraga_P151',
       '2022_Unpub5_DTuraga_P157',
       '2022_Unpub5_DTuraga_P163'
)
matrix_type <- '10x'
MakeRawDataset(study = study, sample_name = sample_name, raw_matrix_type = matrix_type)
MakeCbnDataset(study = study, sample_name = sample_name, raw_matrix_type = matrix_type)
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Dataset #5: 2021_Circulation_EPorrello  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
study <- '2021_Circulation_EPorrello'
sample_name <- c(
        '2021_Circulation_EPorrello_Y2',
        '2021_Circulation_EPorrello_Y3'
)
matrix_type <- '10x'
MakeRawDataset(study = study, sample_name = sample_name, raw_matrix_type = matrix_type)
MakeCbnDataset(study = study, sample_name = sample_name, raw_matrix_type = matrix_type, cb_folder = 'cellbender')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Dataset #6: 2021_Unpub1_DTuraga  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
study <- '2021_Unpub1_DTuraga'
sample_name <- c(
        '2021_Unpub1_DTuraga_P93_s1',
        '2021_Unpub1_DTuraga_P125_s1'
        )
matrix_type <- '10x'
MakeRawDataset(study = study, sample_name = sample_name, raw_matrix_type = matrix_type)
MakeCbnDataset(study = study, sample_name = sample_name, raw_matrix_type = matrix_type)
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Dataset #7: 2022_Unpub6_DTuraga - Pre/Post VAD Samples  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
study <- '2022_Unpub6_DTuraga'
sample_name <- c(
        '2022_Unpub6_DTuraga_P91'
)
matrix_type <- '10x'
MakeRawDataset(study = study, sample_name = sample_name, raw_matrix_type = matrix_type)
MakeCbnDataset(study = study, sample_name = sample_name, raw_matrix_type = matrix_type)
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Dataset #8: 2023_Pediatric1_DTuraga - Pre/Post VAD Samples  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
study <- '2023_Pediatric1_DTuraga'
sample_name <- c(
        '2023_Pediatric1_DTuraga_P195',
        '2023_Pediatric1_DTuraga_P197'
)
matrix_type <- '10x'
MakeRawDataset(study = study, sample_name = sample_name, raw_matrix_type = matrix_type)
MakeCbnDataset(study = study, sample_name = sample_name, raw_matrix_type = matrix_type)
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Dataset #9: 2023_Pediatric2_DTuraga - Pre/Post VAD Samples  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
study <- '2023_Pediatric2_DTuraga'
sample_name <- c(
        '2023_Pediatric2_DTuraga_P162',
        '2023_Pediatric2_DTuraga_P138',
        '2023_Pediatric2_DTuraga_P89',
        '2023_Pediatric2_DTuraga_P203'
)
matrix_type <- '10x'
MakeRawDataset(study = study, sample_name = sample_name, raw_matrix_type = matrix_type)
MakeCbnDataset(study = study, sample_name = sample_name, raw_matrix_type = matrix_type)
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
