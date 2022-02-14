###########    ---------   Get data to start    ----------    ###########
###########    -- for the COVID-19 CD4 follow-up story  --    ###########

# Version.
# Full version: 0.1


### -------------------------- Description -------------------------- ###
# Script to get the basic data to get the final figures for the COVID-19 CD4 follow-up story.

# NOTE: If any of the objects previously processed need any further processing, you may as well load the already processed object and do so for this session.


### --------------------------- Libraries --------------------------- ###
library(Seurat)
library(stringr)
library(data.table)
library(tidyr)
library(gtools)


### --------------------------- Functions --------------------------- ###
coll.version <- 0.8
coll.file <- paste0('/mnt/BioAdHoc/Groups/vd-vijay/vfajardo/COVID-19/paper_developments/COVID-19_CD4-FollowUp/final_figures/jobs_scripts/functions_collection.', coll.version, '.R')
file.exists(coll.file)
source(coll.file)

### ----------------------- General Arguments ----------------------- ###

# ---> General definitions.
# @ Seed.
set.seed(seed=1)
# @ Cluster labels for each dataset.
this.figure <- 'figure_stack'
pt.six.all.sf.clusts.lab <- 'RNA_snn_res.0.3'
pt.tf.all.sf.clusts.lab <- 'RNA_snn_res.0.1'
pt.six.conv.sfce.clusts.lab <- 'RNA_snn_res.0.1'
hlty.six.expa.fce.clusts.lab <- 'RNA_snn_res.0.1'
hlty.six.expb.sfc.clusts.lab <- 'RNA_snn_res.0.3'
hlty.six.expc.all.clusts.lab <- 'RNA_snn_res.0.3'
hlty.six.expd.s3.clusts.lab <- 'RNA_snn_res.0.1'
# ---> Path definitions.
final.figs.path <- '/mnt/beegfs/vfajardo/COVID-19/paper_developments/COVID-19_CD4-FollowUp/final_figures'
gen.data.path <- paste0(final.figs.path, '/data_to_start')
module.sigs.path <- paste0(gen.data.path, '/module_signatures')
# ---> File definitions.
# Seurat objects.
pt.six.all.sf.seurat.obj.file <- '/mnt/beegfs/vfajardo/COVID-19/seurat_analysis/COVID_GenEx_CD4/COVID_GenEx_Batches-1-to-7_CD4_STIM_6_PT/COVID_GenEx_Batches-1-to-7_CD4_STIM_6_PT_07-26-2021_qc-std_var-20_pc-30_hto-class_stim-stim_virus-all_stimtime-six/seurat_objects/SeuratObjectForPrjCOVID_GenEx_Batches-1-to-7_CD4_STIM_6_PT_WithArgs_NoPCs_30_WithTCRTags_ChangedOn_2021-09-30.RDS'
pt.tf.all.sf.seurat.obj.file <- '/mnt/beegfs/vfajardo/COVID-19/seurat_analysis/COVID_GenEx_CD4/COVID_GenEx_Batches-1-to-7_CD4_STIM_24_PT/COVID_GenEx_Batches-1-to-7_CD4_STIM_24_PT_07-26-2021_qc-std_var-30_pc-30_hto-class_stim-stim_virus-all_stimtime-twentyfour/seurat_objects/SeuratObjectForPrjCOVID_GenEx_Batches-1-to-7_CD4_STIM_24_PT_WithArgs_NoPCs_30_WithTCRTags_ChangedOn_2021-11-16.RDS' # NOTE
pt.six.conv.sfce.seurat.obj.file <- '/mnt/beegfs/vfajardo/COVID-19/seurat_analysis/COVID_GenEx_CD4/COVID_GenEx_Batches-4-5-6-7-8_CD4_STIM_6_PT/COVID_GenEx_Batches-4-5-6-7-8_CD4_STIM_6_PT_04-13-2021_qc-std_var-18_pc-16_hto-all_stim-stim_virus-all_stimtime-six/seurat_objects/SeuratObjectForPrjCOVID_GenEx_Batches-4-5-6-7-8_CD4_STIM_6_PT_WithArgs_NoPCs_16_WithTCRTags_ChangedOn_2021-04-20.RDS'
hlty.six.expa.fce.seurat.obj.file <- '/mnt/beegfs/vfajardo/COVID-19/seurat_analysis/COVID_GenEx_CD4/COVID_GenEx_Batch0_CD4_STIM_6_HLTY/COVID_GenEx_Batch0_CD4_STIM_6_HLTY_05-13-2021_qc-std_var-30_pc-30_hto-all_stim-stim_virus-all_stimtime-six/seurat_objects/SeuratObjectForPrjCOVID_GenEx_Batch0_CD4_STIM_6_HLTY_WithArgs_NoPCs_30_WithTCRTags_ChangedOn_2021-05-26.RDS'
hlty.six.expb.sfc.seurat.obj.file <- '/mnt/beegfs/vfajardo/COVID-19/seurat_analysis/COVID_GenEx_CD4/COVID-CD4_ExpVal_Batches-1_Exp-2/COVID-CD4_ExpVal_Batches-1_Exp-2_09-26-2021_qc-std_var-30_pc-30_hto-class/seurat_objects/SeuratObjectForPrjCOVID-CD4_ExpVal_Batches-1_Exp-2_WithArgs_NoPCs_30_WithTCRTags_ChangedOn_2021-09-30.RDS'
hlty.six.expc.all.seurat.obj.file <- '/mnt/beegfs/vfajardo/COVID-19/seurat_analysis/COVID_GenEx_CD4/COVID-CD4_ExpVal_Batches-1_Exp-1/COVID-CD4_ExpVal_Batches-1_Exp-1_09-26-2021_qc-std_var-20_pc-20_hto-class/seurat_objects/SeuratObjectForPrjCOVID-CD4_ExpVal_Batches-1_Exp-1_WithArgs_NoPCs_20_WithTCRTags_ChangedOn_2021-09-30.RDS'
hlty.six.expd.s3.seurat.obj.file <- '/mnt/beegfs/vfajardo/COVID-19/seurat_analysis/COVID_GenEx_CD4/COVID_GenEx_Batches-4-5-6-7_CD4_STIM_6_HLTY/COVID_GenEx_Batches-4-5-6-7_CD4_STIM_6_HLTY_05-26-2021_qc-std_var-30_pc-30_hto-all_stim-stim_virus-all_stimtime-six/seurat_objects/SeuratObjectForPrjCOVID_GenEx_Batches-4-5-6-7_CD4_STIM_6_HLTY_WithArgs_NoPCs_30_WithTCRTags_ChangedOn_2021-05-26.RDS'
# ---> Check directories and files.
if(!all(dir.exists(gen.data.path), dir.exists(module.sigs.path))) stop(paste0('Following paths must be already defined:\n', gen.data.path, '\n', module.sigs.path))
reports.path <- gen.data.path
essential.files <- c(
  pt.six.all.sf.seurat.obj.file,
  pt.tf.all.sf.seurat.obj.file,
  pt.six.conv.sfce.seurat.obj.file,
  hlty.six.expa.fce.seurat.obj.file,
  hlty.six.expb.sfc.seurat.obj.file,
  hlty.six.expc.all.seurat.obj.file,
  hlty.six.expd.s3.seurat.obj.file
)
essential.files <- essential.files[!unlist(lapply(X=essential.files, FUN=file.exists))]
if(length(essential.files) > 0) stop(paste0('Next essential files are not appropriately defined -make sure you\'ve got their names right-:\n', paste0(essential.files, collapse='\n'), '\n'))


### ------------------------- Data Loading -------------------------- ###

# ---> Seurat obejcts.
pt.six.all.sf.obj <- readRDS(file=pt.six.all.sf.seurat.obj.file)
pt.tf.all.sf.obj <- readRDS(file=pt.tf.all.sf.seurat.obj.file)
pt.six.conv.sfce.obj <- readRDS(file=pt.six.conv.sfce.seurat.obj.file)
hlty.six.expa.fce.obj <- readRDS(file=hlty.six.expa.fce.seurat.obj.file)
hlty.six.expb.sfc.obj <- readRDS(file=hlty.six.expb.sfc.seurat.obj.file)
hlty.six.expc.all.obj <- readRDS(file=hlty.six.expc.all.seurat.obj.file)
hlty.six.expd.s3.obj <- readRDS(file=hlty.six.expd.s3.seurat.obj.file)


### ---------------------- Data preprocessing ----------------------- ###

# ---> PT-6-ALL-SF
pt.six.all.sf.obj <- preprocess.obj(seurat.obj=pt.six.all.sf.obj, set.lab='PT-6-ALL-SF', clusts.lab=pt.six.all.sf.clusts.lab, pr.infer=TRUE, get.metadata=FALSE)

# ---> PT-24-ALL-SF
pt.tf.all.sf.obj <- preprocess.obj(seurat.obj=pt.tf.all.sf.obj, set.lab='PT-24-ALL-SF', clusts.lab=pt.tf.all.sf.clusts.lab, pr.infer=TRUE, get.metadata=FALSE)

# ---> PT-6-Conv-SFCE
pt.six.conv.sfce.obj <- preprocess.obj(seurat.obj=pt.six.conv.sfce.obj, set.lab='PT-6-Conv-SFCE', clusts.lab=pt.six.conv.sfce.clusts.lab, pr.infer=TRUE, get.metadata=FALSE)

# ---> HLTY-6-EXPA-FCE
hlty.six.expa.fce.obj <- preprocess.obj(seurat.obj=hlty.six.expa.fce.obj, set.lab='HLTY-6-ExpA-FCE', clusts.lab=hlty.six.expa.fce.clusts.lab, pr.infer=TRUE, get.metadata=FALSE)

# ---> HLTY-6-EXPB-SFC
hlty.six.expb.sfc.obj <- preprocess.obj(seurat.obj=hlty.six.expb.sfc.obj, set.lab='HLTY-6-ExpB-SFC', clusts.lab=hlty.six.expb.sfc.clusts.lab, pr.infer=TRUE, get.metadata=FALSE)

# ---> HLTY-6-EXPC-ALL
hlty.six.expc.all.obj <- preprocess.obj(seurat.obj=hlty.six.expc.all.obj, set.lab='HLTY-6-ExpC-ALL', clusts.lab=hlty.six.expc.all.clusts.lab, pr.infer=TRUE, get.metadata=FALSE)

# ---> HLTY-6-EXPD-S3
hlty.six.expd.s3.obj <- preprocess.obj(seurat.obj=hlty.six.expd.s3.obj, set.lab='HLTY-6-ExpD-S3', clusts.lab=hlty.six.expd.s3.clusts.lab, pr.infer=TRUE, get.metadata=FALSE)


### ------------------------- Main program -------------------------- ###

# ---> PT-6-ALL-SF
tmp.file.name <- paste0(reports.path, '/SeuratObj_PT-6-ALL-SF_Processed.RDS')
saveRDS(file=tmp.file.name, object=pt.six.all.sf.obj)

# ---> PT-24-ALL-SF
tmp.file.name <- paste0(reports.path, '/SeuratObj_PT-24-ALL-SF_Processed.RDS')
saveRDS(file=tmp.file.name, object=pt.tf.all.sf.obj)

# ---> PT-6-Conv-SFCE
tmp.file.name <- paste0(reports.path, '/SeuratObj_PT-6-Conv-SFCE_Processed.RDS')
saveRDS(file=tmp.file.name, object=pt.six.conv.sfce.obj)

# ---> HLTY-6-EXPA-FCE
tmp.file.name <- paste0(reports.path, '/SeuratObj_HLTY-6-ExpA-FCE_Processed.RDS')
saveRDS(file=tmp.file.name, object=hlty.six.expa.fce.obj)

# ---> HLTY-6-EXPB-SFC
tmp.file.name <- paste0(reports.path, '/SeuratObj_HLTY-6-ExpB-SFC_Processed.RDS')
saveRDS(file=tmp.file.name, object=hlty.six.expb.sfc.obj)

# ---> HLTY-6-EXPC-ALL
tmp.file.name <- paste0(reports.path, '/SeuratObj_HLTY-6-ExpC-ALL_Processed.RDS')
saveRDS(file=tmp.file.name, object=hlty.six.expc.all.obj)

# ---> HLTY-6-EXPD-S3
tmp.file.name <- paste0(reports.path, '/SeuratObj_HLTY-6-ExpD-S3_Processed.RDS')
saveRDS(file=tmp.file.name, object=hlty.six.expd.s3.obj)
