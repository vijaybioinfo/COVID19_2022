############    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^    ############
############    -----   Preliminary tables stack    -----    ############
############    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^    ############

# ---> About the script.
# Version: 0
# Subversion: 1
# Script to get all bioinformatics-based panel figures for any of the supplementary tables for the paper known as COVID CD4 II ('follow-up') from the general COVID project of Vijay's lab.


############    -----------------------------------------    ############
### -------------------------- Description -------------------------- ###
############    -----------------------------------------    ############

# Copy description from DICE's CD4STIM project's script (figure 1), but if possible, include in a different file.


############    -----------------------------------------    ############
### --------------------------- Libraries --------------------------- ###
############    -----------------------------------------    ############
library(Seurat)
library(stringr)
library(data.table)
library(tidyr)
library(gtools)
library(rjson)


############    -----------------------------------------    ############
### --------------------------- Functions --------------------------- ###
############    -----------------------------------------    ############

coll.version <- 0.8
coll.file <- paste0('/mnt/BioAdHoc/Groups/vd-vijay/vfajardo/COVID-19/paper_developments/COVID-19_CD4-FollowUp/final_figures/jobs_scripts/functions_collection.', coll.version, '.R')
file.exists(coll.file)
source(coll.file)
total.workers <- try(expr=as.integer(system(command='echo $PBS_NUM_PPN', intern=TRUE)))


############    -----------------------------------------    ############
### ----------------------- General Arguments ----------------------- ###
############    -----------------------------------------    ############

# ---> General definitions.
# @ Paper.
this.prj <- 'COVID-19_CD4-FollowUp'
this.table <- 'table_stack_v1'
# @ Seed.
set.seed(seed=1)
# @ Dataset labels.
pt.six.all.sf <- 'PT-6-ALL-SF'
pt.tf.all.sf <- 'PT-24-ALL-SF'
pt.six.conv.sfce <- 'PT-6-Conv-SFCE'
hlty.six.expa.fce <- 'HLTY-6-ExpA-FCE'
hlty.six.expb.sfc <- 'HLTY-6-ExpB-SFC'
hlty.six.expc.all <- 'HLTY-6-ExpC-ALL'
hlty.six.expd.s3 <- 'HLTY-6-ExpD-S3'
dataset.labs <- c(pt.six.all.sf, pt.tf.all.sf, pt.six.conv.sfce, hlty.six.expa.fce, hlty.six.expb.sfc, hlty.six.expc.all, hlty.six.expd.s3)
# @ Cluster labels for each dataset.
pt.six.all.sf.clusts.lab <- 'RNA_snn_res.0.3'
pt.tf.all.sf.clusts.lab <- 'RNA_snn_res.0.1'
pt.six.conv.sfce.clusts.lab <- 'RNA_snn_res.0.1'
hlty.six.expa.fce.clusts.lab <- 'RNA_snn_res.0.1'
hlty.six.expb.sfc.clusts.lab <- 'RNA_snn_res.0.3'
hlty.six.expc.all.clusts.lab <- 'RNA_snn_res.0.3'
hlty.six.expd.s3.clusts.lab <- 'RNA_snn_res.0.1'
clusts.lab <- c(pt.six.all.sf.clusts.lab, pt.tf.all.sf.clusts.lab, pt.six.conv.sfce.clusts.lab, hlty.six.expa.fce.clusts.lab, hlty.six.expb.sfc.clusts.lab, hlty.six.expc.all.clusts.lab, hlty.six.expd.s3.clusts.lab)
names(clusts.lab) <- dataset.labs
# @ Gene signature definitions.
cd4.ctl.sign <- 'cell.cytotoxicity.patil.score'
# @ For subsets.
lane.tags <- c(
  `PT-6-ALL-SF`='cv.batch.tag;species.tag;arte.pep.pool.tag;stim.time.tag;cell.type.tag;tenx.chemistry.tag;donors.no.tag;virus.tag;stim.tag;library.tag;cv.arte.pep.pool.tag;donors.kind.tag;blood.draw.phase.tag',
  `PT-24-ALL-SF`='cv.batch.tag;species.tag;arte.pep.pool.tag;stim.time.tag;cell.type.tag;tenx.chemistry.tag;donors.no.tag;virus.tag;stim.tag;library.tag;cv.arte.pep.pool.tag;donors.kind.tag;blood.draw.phase.tag',
  `PT-6-Conv-SFCE`='cv.batch.tag;species.tag;arte.pep.pool.tag;stim.time.tag;cell.type.tag;tenx.chemistry.tag;donors.no.tag;virus.tag;stim.tag;library.tag;cv.arte.pep.pool.tag;donors.kind.tag',
  `HLTY-6-ExpA-FCE`='cv.batch.tag;virus.tag;stim.tag;stim.time.tag;species.tag;cell.type.tag;donors.no.tag;arte.pep.pool.tag;donors.kind.tag',
  `HLTY-6-ExpB-SFC`='cell.type.tag;exp.val.tag;orig.chrom.batch.tag',
  `HLTY-6-ExpC-ALL`='cell.type.tag;exp.val.tag;orig.chrom.batch.tag',
  `HLTY-6-ExpD-S3`='cv.batch.tag;species.tag;arte.pep.pool.tag;stim.time.tag;cell.type.tag;tenx.chemistry.tag;donors.no.tag;virus.tag;stim.tag;library.tag;cv.arte.pep.pool.tag;donors.kind.tag'
)
# ---> Path definitions.
final.tables.path <- paste0('/mnt/BioAdHoc/Groups/vd-vijay/vfajardo/COVID-19/paper_developments/', this.prj, '/final_tables')
gen.data.path <- paste0(final.tables.path, '/data_to_start')
this.table.path <- paste0(final.tables.path, '/', this.table)
reports.path <- this.table.path
# ---> File definitions.
# Seurat objects.
objs.files.list <- paste0(gen.data.path, '/', 'SeuratObj_', dataset.labs, '_Processed.RDS')
# Cell counts previous to filtering.
cell.count.files <- paste0(gen.data.path, '/CellCountSummary_', dataset.labs, '.csv')
# Gene sets' files (for FGSE analyses).
module.feats.dir <- paste0(gen.data.path, '/module_signatures_for_fgsea')
# Cluster markers.
marker.files.list <- paste0(gen.data.path, '/ClusterMarkers_', dataset.labs, '.csv')
# Specific files.
#     SARS-CoV-2-reactive cell count per donor.
sr.cell.counts.file <- paste0(gen.data.path, '/SortedCellsPerMillionPBMCs.csv')
#     DEA, HLTY-6-ExpD-S3, HCoV vs SARS-CoV-2 cells.
# hlty.six.expd.s3.dea.file <- paste0(gen.data.path, '/HLTY-6-ExpD-S3_DEA_HCoV-vs-SARS.csv')
#     CMV seropositivity test
cmv.igg.res.file <- paste0(gen.data.path, '/CMVSeropositivityTestResults.csv')
#     HLA typing resuls.
hla.data.file <- paste0(gen.data.path, '/DataFromHLATypingServiceProvider.csv')
#     COVID-19 donor order.
donor.order.file <- paste0(gen.data.path, '/COVID-19DonorOrder.csv')
# ---> Check directories and files.
if(!all(dir.exists(gen.data.path), dir.exists(this.table.path))) stop(paste0('Following paths must be already defined:\n', gen.data.path, '\n', this.table.path))
if(!dir.exists(reports.path)) dir.create(reports.path)
essential.files <- c(
  lapply(X=objs.files.list, FUN=function(x) return(x))
)
essential.files <- essential.files[!unlist(lapply(X=essential.files, FUN=file.exists))]
if(length(essential.files) > 0) stop(paste0('Next essential files are not appropriately defined -make sure you\'ve got their names right-:\n', paste0(essential.files, collapse='\n'), '\n'))
# ---> Present arguments (to do)


############    -----------------------------------------    ############
### ------------------------- Data Loading -------------------------- ###
############    -----------------------------------------    ############

# ---> Seurat obejcts.
srt.objs.list <- lapply(X=objs.files.list, FUN=readRDS)
names(srt.objs.list) <- dataset.labs

# ---> CLuster markers.
markers.list <- lapply(X=marker.files.list, FUN=fread)
names(markers.list) <- dataset.labs

# ---> Specific files.
#     Cell counts summary.
cell.count.summs <- lapply(X=cell.count.files, FUN=fread)
names(cell.count.summs) <- dataset.labs
#     SARS-CoV-2-reactive cell count per donor.
sr.cell.counts <- fread(file=sr.cell.counts.file)
sr.cell.counts <- sr.cell.counts[time.point=='6', .(patient.id, cells.acute, cells.memory, interval.acute, interval.memory)]
#     DEA, HLTY-6-ExpD-S3, HCoV vs SARS-CoV-2 cells.
# hlty.six.expd.s3.dea <- fread(file=hlty.six.expd.s3.dea.file)
#     CMV seropositivity test
cmv.igg.res <- fread(file=cmv.igg.res.file)
#     HLA typing resuls.
hla.data <- fread(file=hla.data.file)
#     COVID-19 donor order.
donor.order <- fread(file=donor.order.file)

# ---> Gene sets' files (for FGSE analyses).
# module.feats.files <- list.files(path=module.feats.dir, recursive=FALSE, include.dirs=FALSE, pattern='signature.csv$', full.names=TRUE)
# modules.feats <- lapply(X=module.feats.files, FUN=function(tmp.file){
#   tmp.feats <- read.csv(file=tmp.file, stringsAsFactors=FALSE)
#   if(!'feature' %in% colnames(tmp.feats)) stop(paste0('File ', tmp.file, ' not appropriately defined. Column listing the gene features should be named \'feature\'.\n\n'))
#   tmp.feats <- tmp.feats$feature
#   # Replace underscores with points in the gene names just as Seurat does when one creates a seurat object.
#   tmp.check <- any(str_detect(string=tmp.feats, pattern='_'))
#   if(tmp.check){
#     tmp.warning <- paste0('Feature names cannot have underscores (\'_\'), replacing with dashes (\'-\'). This for file: ', basename(tmp.file))
#     warning(tmp.warning)
#     tmp.feats <- str_replace_all(string=tmp.feats, pattern='_', replacement='-')
#   }
#   # Return.
#   return(tmp.feats)
# })
# names(modules.feats) <- list.files(path=module.feats.dir, recursive=FALSE, include.dirs=FALSE, pattern='signature.csv$', full.names=FALSE)
# names(modules.feats) <- str_replace(string=names(modules.feats), pattern='_signature.csv$', replacement='')
# names(modules.feats) <- str_replace_all(string=names(modules.feats), pattern='_', replacement='.')
# modules.names <- names(modules.feats)
# names(modules.names) <- str_replace_all(string=modules.names, pattern='\\.', replacement='_')


############    -----------------------------------------    ############
### ---------------------- Data preprocessing ----------------------- ###
############    -----------------------------------------    ############

# ---> PT-6-ALL-SF
# @ Clusters definition.
# Merge clusters 1 and 7, the CD4-CTL clusters, into a single cluster (cluster 1).
new.vals <- as.character(srt.objs.list[[pt.six.all.sf]]@meta.data[, clusts.lab[pt.six.all.sf]])
new.vals[new.vals=='7'] <- '1'
new.lvls <- mixedsort(unique(new.vals))
srt.objs.list[[pt.six.all.sf]]@meta.data[, clusts.lab[pt.six.all.sf]] <- factor(x=new.vals, levels=new.lvls)

# ---> PT-24-ALL-SF
# @ Clusters definition.
# Remove cluster 5 (around ~150 cells that seem not to be T cells).
cells.to.keep <- as.character(srt.objs.list[[pt.tf.all.sf]]@meta.data[, clusts.lab[pt.tf.all.sf]])!='5'
cells.to.keep <- Cells(srt.objs.list[[pt.tf.all.sf]])[cells.to.keep]
srt.objs.list[[pt.tf.all.sf]] <- subset(x=srt.objs.list[[pt.tf.all.sf]], cells=cells.to.keep)

# ---> PT-6-Conv-SFCE


# ---> HLTY-6-EXPA-FCE


# ---> HLTY-6-EXPB-SFC
# @ Clusters definition.
# Remove clusters 7 and 8 (around ~150 cells in each cluster and together they account for less than 1% of the whole dataset).
cells.to.keep <- !as.character(srt.objs.list[[hlty.six.expb.sfc]]@meta.data[, clusts.lab[hlty.six.expb.sfc]])%in%c('7', '8')
cells.to.keep <- Cells(srt.objs.list[[hlty.six.expb.sfc]])[cells.to.keep]
srt.objs.list[[hlty.six.expb.sfc]] <- subset(x=srt.objs.list[[hlty.six.expb.sfc]], cells=cells.to.keep)
# Tag defined by virus reactivity and ARTE peptide pool concentration.
srt.objs.list[[hlty.six.expb.sfc]]@meta.data[, 'virus.arte.con.tag'] <- paste(srt.objs.list[[hlty.six.expb.sfc]]@meta.data[, 'virus.tag'], srt.objs.list[[hlty.six.expb.sfc]]@meta.data[, 'arte.concentration.tag'], sep=';')

# ---> HLTY-6-EXPC-ALL
# @ Clusters definition.
# Remove clusters 6, 7 and 8 (less than 1% of the whole dataset).
cells.to.keep <- !as.character(srt.objs.list[[hlty.six.expc.all]]@meta.data[, clusts.lab[hlty.six.expb.sfc]])%in%c('6', '7', '8')
cells.to.keep <- Cells(srt.objs.list[[hlty.six.expc.all]])[cells.to.keep]
srt.objs.list[[hlty.six.expc.all]] <- subset(x=srt.objs.list[[hlty.six.expc.all]], cells=cells.to.keep)
# @ SLAMF7 population tag.
srt.objs.list[[hlty.six.expc.all]]@meta.data[, 'slamf7.pop.tag'] <- ifelse(test=srt.objs.list[[hlty.six.expc.all]]@meta.data[, 'cell.pop.tag']=='SLAMF7+', yes='Positive', no='Negative')

# ---> Specific files.
#     DEA, HLTY-6-ExpD-S3, HCoV vs SARS-CoV-2 cells.
# hlty.six.expd.s3.dea <- fread(file=hlty.six.expd.s3.dea.file)

#     CMV seropositivity test
# Harmonize donor ID.
colnames(cmv.igg.res) <- c('sample.id', 'assay', 'results')
cmv.igg.res[, donor.id:=str_extract(string=sample.id, pattern='^[0-9]+')]
cmv.igg.res[, donor.id:=str_replace(string=donor.id, pattern='000$', replacement='')]
# Test results.
cmv.igg.res[, donor.id:=paste0('P', ifelse(test=as.numeric(donor.id)<10, yes='0', no=''), donor.id)]
cmv.igg.res[, cmv.detected:=str_detect(string=results, pattern='\\(')]
cmv.igg.res[, test.result:=ifelse(test=cmv.detected==TRUE, yes='Detected', no='Not detected')]
# Ag amount.
cmv.igg.res[, ig.g:=str_extract(string=results, pattern='\\(.+\\)$')]
cmv.igg.res[, ig.g:=str_replace_all(string=ig.g, pattern='\\(|\\)', replacement='')]
cmv.igg.res[, ig.g:=str_replace(string=ig.g, pattern='\\s+U/ml|\\s+U/ML', replacement='')]
cmv.igg.res[, ig.g:=as.numeric(ig.g)]
cmv.igg.res[, results:=NULL]
# Take consensus.
cmv.igg.res <- cmv.igg.res[, .(donor.id, cmv.detected, test.result, ig.g)]
cmv.igg.res <- cmv.igg.res[,
  .(
    cmv.detected=sum(cmv.detected)>0,
    ig.g=min(ig.g, na.rm=TRUE)
  ),
  by=donor.id
]
cmv.igg.res[ig.g==Inf, ig.g:=NA]
cmv.igg.res[, test.result:=ifelse(cmv.detected==TRUE, yes='Detected', no='Not detected')]
# cmv.igg.res[cmv.detected==FALSE, unique(ig.g)]


############    -----------------------------------------    ############
### ------------------------ On introduction ------------------------ ###
############    -----------------------------------------    ############

# ----> Output directory.
table.1.path <- paste0(reports.path, '/table_on_introduction')
if(!dir.exists(table.1.path)) dir.create(table.1.path)

# ----> Demographic and clinical characteristics of donors recovered from COVID-19 illness.
# ---> Merge individual pieces of data.
tmp.data <- merge(x=cmv.igg.res[, .(donor.id.tag=donor.id, test.result, ig.g)], y=hla.data, by='donor.id.tag', all=TRUE)
tmp.data <- merge(x=sr.cell.counts[, .(donor.id.tag=patient.id, interval.acute, interval.memory)], y=tmp.data, by='donor.id.tag', all=TRUE)
tmp.data <- merge(x=donor.order, y=tmp.data, by='donor.id.tag', sort=FALSE)
tmp.file.name <- paste0(table.1.path, '/CovidPatientsInformation.csv')
fwrite(file=tmp.file.name, x=tmp.data, quote=FALSE, na=NA)


############    -----------------------------------------    ############
### ---------------- On convalescent Covid patients ----------------- ###
### -------------- Preprocessing and GEx data details --------------- ###
############    -----------------------------------------    ############

# ----> Output directory.
table.2.path <- paste0(reports.path, '/table_on_conv_covid')
if(!dir.exists(table.2.path)) dir.create(table.2.path)

# ---> General definitions.
# Seurat objects of interest.
table.2.objs <- c('PT-6-ALL-SF', 'PT-24-ALL-SF', 'PT-6-Conv-SFCE')

# ---> Aggregation QC summaries.
# @ Load summaries.
aggr.qc.files <- paste0(gen.data.path, '/AggrSummary_', table.2.objs, '.json')
aggr.qcs <- lapply(X=aggr.qc.files, FUN=function(tmp.file) fromJSON(file=tmp.file))
names(aggr.qcs) <- table.2.objs
# @ Get relevant data in a tidy format per aggregation of interest.
aggr.qcs <- lapply(X=names(aggr.qcs), FUN=process.aggr.json)
names(aggr.qcs) <- table.2.objs
# @ Get general aggregation metrics only.
tmp.data <- lapply(X=aggr.qcs, FUN=function(x) return(x[['gen.qc']]))
tmp.data <- rbindlist(l=tmp.data, use.names=TRUE, idcol='aggr')
tmp.data <- spread(data=tmp.data, key=group, value=value)
cols.order <- c('aggr', 'pre.total', 'post.total', 'pre.mean', 'post.mean')
tmp.data <- tmp.data[, ..cols.order]
colnames(tmp.data) <- c(
  'Aggregation ID',
  'Pre-Normalization Total Number of Reads in unstimulated data',
  'Post-Normalization Total Number of Reads in unstimulated data',
  'Pre-Normalization Mean Reads per Cell in unstimulated data',
  'Post-Normalization Mean Reads per Cell in unstimulated data'
)
# @ Output general aggregation metrics only.
tmp.file.name <- paste0(table.2.path, '/QCSumm_Aggr.csv')
fwrite(file=tmp.file.name, x=tmp.data, na=NA, quote=TRUE)
# @ Merge per-library metrics.
ind.aggr.qcs <- lapply(X=aggr.qcs, FUN=function(x) return(x[['libs.qc']]))
ind.aggr.qcs <- Reduce(x=ind.aggr.qcs, f=function(x, y){
  to.output <- merge(x=x, y=y, by='sample.id', all=TRUE)
})
colnames(ind.aggr.qcs) <- c('sample.id', paste0(c('cmb.reads.per.bc.', 'frac.reads.kept.', 'raw.reads.per.bc.'), rep(x=names(aggr.qcs), each=3)))

# ---> Per-library QC summaries.
# @ General QC information for GEx and HTO data.
sample.ids <- lapply(X=table.2.objs, FUN=list.samples.info)
sample.ids <- unique(rbindlist(l=sample.ids, use.names=TRUE, fill=TRUE))
sample.ids[is.na(blood.draw.phase.tag), blood.draw.phase.tag:='convalescent']
lane.info <- unique(sample.ids[, .(sample.id, facs.batch=cv.batch.tag, seq.batch=seq.batch.tag, blood.draw.phase=blood.draw.phase.tag, virus=virus.tag, stim.time=stim.time.tag, donors.no=donors.no.tag)])
lane.info[, donors.no:=str_replace(string=donors.no, pattern='D$', replacement='')]
lane.info[, blood.draw.phase:=ifelse(test=blood.draw.phase=='convalescent', yes='Convalescent', no='Acute')]

# @ GEx summries.
# Define metric summary files.
uniq.sample.ids <- unique(sample.ids[, sample.id])
sample.qc.files <- unique(paste0(uniq.sample.ids, '/metrics_summary.csv'))
to.check <- all(file.exists(sample.qc.files))
# Load data as a single list keepin track of ID.
tmp.data <- lapply(X=sample.qc.files, FUN=fread)
names(tmp.data) <- uniq.sample.ids
tmp.data <- rbindlist(l=tmp.data, use.names=TRUE, idcol='sample.id')
# Merge metadata with metrics.
tmp.data <- merge(x=lane.info, y=tmp.data, by='sample.id')
tmp.data <- merge(x=tmp.data, y=unique(sample.ids[, .(sample.id, library_id)]), by='sample.id')
# Include per-library aggr information.
tmp.data <- merge(x=tmp.data, y=ind.aggr.qcs, by.x='library_id', by.y='sample.id', all=TRUE)
# Redefine sample ID.
tmp.data[, sample.id:=paste0('Covid-CD4_', facs.batch, '_', ifelse(test=as.integer(stim.time)==6, yes='06', no='24'), '_', virus, '_', donors.no, '_GEx')]
tmp.data[, sample.id:=str_replace(string=sample.id, pattern='SARS-CoV-2', replacement='SARS')]
# Final formatting.
cols.order <- c('Estimated Number of Cells',	'Mean Reads per Cell',	'Median Genes per Cell',	'Number of Reads',	'Valid Barcodes',	'Fraction Reads in Cells',	'Total Genes Detected',	'Median UMI Counts per Cell', 'Reads Mapped to Genome',	paste0('Reads Mapped Confidently to ', c('Genome',	'Intergenic Regions', 'Intronic Regions',	'Exonic Regions',	'Transcriptome')), 'Reads Mapped Antisense to Gene', 'Sequencing Saturation',	'Q30 Bases in Barcode',	'Q30 Bases in RNA Read',	'Q30 Bases in RNA Read 2',	'Q30 Bases in Sample Index',	'Q30 Bases in UMI', colnames(ind.aggr.qcs)[2:length(ind.aggr.qcs)])
# all(cols.order %in% colnames(tmp.data))
cols.order <- c(setdiff(x=colnames(tmp.data), y=cols.order), cols.order)
tmp.data <- tmp.data[, ..cols.order]
setorderv(x=tmp.data, cols=c('stim.time', 'blood.draw.phase', 'virus', 'seq.batch', 'sample.id'), order=c(-1, 1, -1, 1, 1))
# Output table.
tmp.file.name <- paste0(table.2.path, '/QCSumm_GEx.csv')
fwrite(file=tmp.file.name, x=tmp.data, na=NA, quote=TRUE)

# @ HTO summary
# Define metric summary files.
uniq.sample.ids <- unique(sample.ids[, .(sample.id, hto.tag)])
for.files <- uniq.sample.ids[, hto.tag]
for.files <- str_replace(string=for.files, pattern='/mnt/BioAdHoc/Groups/vd-vijay/vfajardo/COVID-19/deconvolution/HTO_based/', replacement='')
tmp.dates <- str_extract(string=for.files, pattern='^\\d{2}-\\d{2}-\\d{4}')
for.files <- str_replace(string=for.files, pattern='^\\d{2}-\\d{2}-\\d{4}', replacement=paste0(tmp.dates, '/count'))
for.files <- str_extract(string=for.files, pattern='^([^/]+/){4}')
sample.qc.files <- paste0('/mnt/BioAdHoc/Groups/vd-vijay/vfajardo/sequencing_data/', for.files, '/outs/metrics_summary.csv')
to.check <- all(file.exists(sample.qc.files))
# Load data as a single list keepin track of ID.
tmp.data <- lapply(X=sample.qc.files, FUN=fread)
names(tmp.data) <- uniq.sample.ids[, sample.id]
tmp.data <- rbindlist(l=tmp.data, use.names=TRUE, idcol='sample.id')
# Merge metadata with metrics.
tmp.data <- merge(x=lane.info, y=tmp.data, by='sample.id')
# Redefine sample ID.
tmp.data[, sample.id:=paste0('Covid-CD4_', facs.batch, '_', ifelse(test=as.integer(stim.time)==6, yes='06', no='24'), '_', virus, '_', donors.no, '_CITE')]
tmp.data[, sample.id:=str_replace(string=sample.id, pattern='SARS-CoV-2', replacement='SARS')]
# Final formatting.
colnames(tmp.data) <- str_replace(string=colnames(tmp.data), pattern='^Antibody:\\s', replacement='')
cols.order <- c('Mean Reads per Cell',	'Valid Barcodes', 'Fraction Antibody Reads',	'Fraction Antibody Reads Usable',	'Antibody Reads Usable per Cell',	'Fraction Reads in Barcodes with High UMI Counts',	'Fraction Unrecognized Antibody',	'Antibody Reads in Cells',	'Median UMIs per Cell (summed over all recognized antibody barcodes)', 'Sequencing Saturation',	'Q30 Bases in Barcode',	'Q30 Bases in Antibody Read',	'Q30 Bases in Sample Index',	'Q30 Bases in UMI')
# all(cols.order %in% colnames(tmp.data))
cols.order <- c(setdiff(x=colnames(tmp.data), y=cols.order), cols.order)
tmp.data <- tmp.data[, ..cols.order]
setorderv(x=tmp.data, cols=c('blood.draw.phase', 'stim.time', 'virus', 'seq.batch', 'sample.id'), order=c(1, -1, -1, 1, 1))
# Output table.
tmp.file.name <- paste0(table.2.path, '/QCSumm_HTO.csv')
fwrite(file=tmp.file.name, x=tmp.data, na=NA, quote=TRUE)

# @ TCR summary
# Define metric summary files.
uniq.sample.ids <- unique(sample.ids[, .(sample.id, vdj.path)])
sample.qc.files <- paste0(uniq.sample.ids[, vdj.path], '/metrics_summary.csv')
to.check <- all(file.exists(sample.qc.files))
# Load data as a single list keepin track of ID.
tmp.data <- lapply(X=sample.qc.files, FUN=fread)
names(tmp.data) <- uniq.sample.ids[, sample.id]
tmp.data <- rbindlist(l=tmp.data, use.names=TRUE, idcol='sample.id')
# Merge metadata with metrics.
tmp.data <- merge(x=lane.info, y=tmp.data, by='sample.id')
# Redefine sample ID.
tmp.data[, sample.id:=paste0('Covid-CD4_', facs.batch, '_', ifelse(test=as.integer(stim.time)==6, yes='06', no='24'), '_', virus, '_', donors.no, '_TCR')]
tmp.data[, sample.id:=str_replace(string=sample.id, pattern='SARS-CoV-2', replacement='SARS')]
# Final formatting.
colnames(tmp.data) <- str_replace(string=colnames(tmp.data), pattern='^Antibody:\\s', replacement='')
cols.order <- c('Estimated Number of Cells',	'Mean Read Pairs per Cell',	'Number of Cells With Productive V-J Spanning Pair',	'Number of Read Pairs',	'Valid Barcodes',	'Fraction Reads in Cells',	'Median TRA UMIs per Cell',	'Median TRB UMIs per Cell', 'Reads Mapped to Any V(D)J Gene',	'Reads Mapped to TRA',	'Reads Mapped to TRB',	'Mean Used Read Pairs per Cell', 'Q30 Bases in Barcode',	'Q30 Bases in RNA Read 1',	'Q30 Bases in RNA Read 2',	'Q30 Bases in Sample Index',	'Q30 Bases in UMI', 'Cells With Productive V-J Spanning Pair',	'Cells With Productive V-J Spanning (TRA, TRB) Pair',	'Paired Clonotype Diversity',	'Cells With TRA Contig',	'Cells With TRB Contig',	'Cells With CDR3-annotated TRA Contig',	'Cells With CDR3-annotated TRB Contig',	'Cells With V-J Spanning TRA Contig',	'Cells With V-J Spanning TRB Contig',	'Cells With Productive TRA Contig',	'Cells With Productive TRB Contig')
# all(cols.order %in% colnames(tmp.data))
cols.order <- c(setdiff(x=colnames(tmp.data), y=cols.order), cols.order)
tmp.data <- tmp.data[, ..cols.order]
setorderv(x=tmp.data, cols=c('blood.draw.phase', 'stim.time', 'virus', 'seq.batch', 'sample.id'), order=c(1, -1, -1, 1, 1))
# Output table.
tmp.file.name <- paste0(table.2.path, '/QCSumm_TCR.csv')
fwrite(file=tmp.file.name, x=tmp.data, na=NA, quote=TRUE)

# ---> Cell counts per condition per donor.
# @ Obtain each table.
sep.conv.tables <- lapply(X=table.2.objs, FUN=get.sep.conv.table)
names(sep.conv.tables) <- table.2.objs
# @ Output each table.
lapply(X=names(sep.conv.tables), FUN=function(table.obj){
  # Retrieve data.
  tmp.data <- sep.conv.tables[[table.obj]]
  # Clusters order.
  cluster.cols <- grep(x=colnames(tmp.data), pattern='^\\d+$', value=TRUE)
  other.cols <- setdiff(x=colnames(tmp.data), y=cluster.cols)
  cols.order <- c(other.cols, mixedsort(cluster.cols))
  tmp.data <- tmp.data[, ..cols.order]
  # Sort
  setorderv(x=tmp.data, cols=c('Virus peptide pool', 'Hospitalization status', 'Sample collention timepoint', 'Donor'), order=c(-1, 1, 1, 1))
  # Output
  tmp.file.name <- paste0(table.2.path, '/CellCountsPerCond_Set-', table.obj, '.csv')
  fwrite(file=tmp.file.name, x=tmp.data, na=NA, quote=FALSE)
})

# ---> Clusters' markers.
# Output definition as a function.
output.markers <- function(obj.of.int, output.path){
  # Retrieve and pre-process table.
  tmp.data <- markers.list[[obj.of.int]]
  tmp.data <- tmp.data[p.adj<=0.05 & avg.lfc>=0.25]
  # Formatting.
  cols.order <- c('gene.id', 'element', 'pct.1', 'pct.2', 'avg.lfc', 'p.val', 'p.adj', grep(x=colnames(tmp.data), pattern='mean.resolution.|mean.cell.pop', value=TRUE), grep(x=colnames(tmp.data), pattern='prop.resolution.|prop.cell.pop', value=TRUE))
  tmp.data <- tmp.data[, ..cols.order]
  setorderv(x=tmp.data, cols=c('element', 'p.adj', 'avg.lfc'), order=c(1, 1, -1))
  # Output.
  tmp.file.name <- paste0(output.path, '/ClusterMarkers_Set-', obj.of.int, '.csv')
  fwrite(file=tmp.file.name, x=tmp.data, na=NA, quote=FALSE)
  return(NA)
}
# Output each table.
for(table.obj in table.2.objs){
  output.markers(obj.of.int=table.obj, output.path=table.2.path)
}


############    -----------------------------------------    ############
### ---------------- On convalescent Covid patients ----------------- ###
### ---------------------- scTCR data details ----------------------- ###
############    -----------------------------------------    ############

# ----> Output directory.
table.3.path <- paste0(reports.path, '/table_on_conv_covid_tcr')
if(!dir.exists(table.3.path)) dir.create(table.3.path)

# ---> General definitions.
# Seurat objects of interest.
table.3.objs <- c('PT-6-ALL-SF', 'PT-24-ALL-SF', 'PT-6-Conv-SFCE')

# ---> TCR extended information.

# Custom cross-reactivity info and phase tag info availability per dataset.
sets.cr.info <- rep(x=TRUE, times=length(table.3.objs)); names(sets.cr.info) <- table.3.objs
sets.phase.info <- c(
  `PT-6-ALL-SF`=TRUE,
  `PT-24-ALL-SF`=TRUE,
  `PT-6-Conv-SFCE`=FALSE
)
# Do process per dataset.
for(table.obj in table.3.objs){
  cat('Process for object:', table.obj, '\n')
  # Obtain vdj gene usage information.
  vdj.gene.info <- get.vdj.info(obj.of.int=table.obj)
  # Obtain data separately for clusters and virus reactivities and then merge.
  cr.info <- sets.cr.info[table.obj]; phase.info <- sets.phase.info[table.obj]
  tcr.by.cluster <- get.tcr.by.tag.table(obj.of.int=table.obj, vdj.gene.info=vdj.gene.info, tag.of.int=clusts.lab[table.obj], phase.info=phase.info, cr.info=cr.info)
  tcr.by.virus <- get.tcr.by.tag.table(obj.of.int=table.obj, vdj.gene.info=vdj.gene.info, tag.of.int='virus.tag', phase.info=phase.info, cr.info=cr.info)
  to.check <- tcr.by.cluster[, .N]==tcr.by.virus[, .N] & tcr.by.cluster[, uniqueN(clonotype)]==tcr.by.virus[, uniqueN(clonotype)]
  if(!to.check) stop('Merging of virus and cluster data unsuccessful, error type 1.\n')
  tmp.data <- merge(x=tcr.by.virus, y=tcr.by.cluster, by=c('clonotype', 'trb.aa', 'tra.aa', 'trb.nt', 'tra.nt', 'trb.v', 'tra.v', 'trb.j', 'tra.j', 'pr.tag', 'donor', 'phase', 'global.size', 'local.size'), suffixes=c('.virus', '.cluster'))
  to.check <- tmp.data[, .N]==tcr.by.cluster[, .N] & tmp.data[, .N]==tcr.by.virus[, .N]
  if(!to.check) stop('Merging of virus and cluster data unsuccessful, error type 2.\n')
  # Final formatting.
  tmp.data[, to.sort:=as.integer(str_replace(string=clonotype, pattern='clonotype', replacement=''))]
  setorderv(x=tmp.data, cols=c('phase', 'to.sort', 'donor'), order=c(-1, 1, 1)); tmp.data[, to.sort:=NULL]
  if(!cr.info){
    tmp.data[, pr.tag:=NULL]
  }else{
    tmp.data[, pr.tag:=ifelse(test=pr.tag=='pR', yes='Cross-reactive', no='Other')]
  }
  if(!phase.info){
    tmp.data[, phase:=NULL]
  }else{
    tmp.data[, phase:=ifelse(test=phase=='acute', yes='Acute', no='Convalescent')]
  }
  # Output.
  tmp.file.name <- paste0(table.3.path, '/TCRExtendedData_Set-', table.obj, '.csv')
  fwrite(file=tmp.file.name, x=tmp.data, na=NA, quote=FALSE)
}


############    -----------------------------------------    ############
### ----------------- On cross-reactive CD4 T cells ----------------- ###
### ----------------------- in healthy donors ----------------------- ###
############    -----------------------------------------    ############

# ----> Output directory.
table.4.path <- paste0(reports.path, '/table_on_cr_in_hlty')
if(!dir.exists(table.4.path)) dir.create(table.4.path)

# ---> General definitions.
# Seurat objects of interest.
table.4.objs <- c('HLTY-6-ExpA-FCE', 'HLTY-6-ExpD-S3')

# ---> Aggregation QC summaries.
# @ Load summaries.
aggr.qc.files <- paste0(gen.data.path, '/AggrSummary_', table.4.objs, '.json')
aggr.qcs <- lapply(X=aggr.qc.files, FUN=function(tmp.file) fromJSON(file=tmp.file))
names(aggr.qcs) <- table.4.objs
# @ Get relevant data in a tidy format per aggregation of interest.
aggr.qcs <- lapply(X=names(aggr.qcs), FUN=process.aggr.json)
names(aggr.qcs) <- table.4.objs
# @ Get general aggregation metrics only.
tmp.data <- lapply(X=aggr.qcs, FUN=function(x) return(x[['gen.qc']]))
tmp.data <- rbindlist(l=tmp.data, use.names=TRUE, idcol='aggr')
tmp.data <- spread(data=tmp.data, key=group, value=value)
cols.order <- c('aggr', 'pre.total', 'post.total', 'pre.mean', 'post.mean')
tmp.data <- tmp.data[, ..cols.order]
colnames(tmp.data) <- c(
  'Aggregation ID',
  'Pre-Normalization Total Number of Reads in unstimulated data',
  'Post-Normalization Total Number of Reads in unstimulated data',
  'Pre-Normalization Mean Reads per Cell in unstimulated data',
  'Post-Normalization Mean Reads per Cell in unstimulated data'
)
# @ Output general aggregation metrics only.
tmp.file.name <- paste0(table.4.path, '/QCSumm_Aggr.csv')
fwrite(file=tmp.file.name, x=tmp.data, na=NA, quote=TRUE)
# @ Merge per-library metrics.
ind.aggr.qcs <- lapply(X=aggr.qcs, FUN=function(x) return(x[['libs.qc']]))
ind.aggr.qcs <- Reduce(x=ind.aggr.qcs, f=function(x, y){
  to.output <- merge(x=x, y=y, by='sample.id', all=TRUE)
})
colnames(ind.aggr.qcs) <- c('sample.id', paste0(c('cmb.reads.per.bc.', 'frac.reads.kept.', 'raw.reads.per.bc.'), rep(x=names(aggr.qcs), each=3)))

# ---> Per-library QC summaries.
# @ General QC information for GEx and HTO data.
sample.ids <- lapply(X=table.4.objs, FUN=list.samples.info)
sample.ids <- unique(rbindlist(l=sample.ids, use.names=TRUE, fill=TRUE))
lane.info <- unique(sample.ids[, .(sample.id, facs.batch=cv.batch.tag, seq.batch=seq.batch.tag, virus=virus.tag, stim.time=stim.time.tag, donors.no=donors.no.tag)])
lane.info[, donors.no:=str_replace(string=donors.no, pattern='D$', replacement='')]

# @ GEx summries.
# Define metric summary files.
uniq.sample.ids <- unique(sample.ids[, sample.id])
sample.qc.files <- unique(paste0(uniq.sample.ids, '/metrics_summary.csv'))
to.check <- all(file.exists(sample.qc.files))
# Load data as a single list keepin track of ID.
tmp.data <- lapply(X=sample.qc.files, FUN=fread)
names(tmp.data) <- uniq.sample.ids
tmp.data <- rbindlist(l=tmp.data, use.names=TRUE, idcol='sample.id')
# Merge metadata with metrics.
tmp.data <- merge(x=lane.info, y=tmp.data, by='sample.id')
tmp.data <- merge(x=tmp.data, y=unique(sample.ids[, .(sample.id, library_id)]), by='sample.id')
# Include per-library aggr information.
tmp.data <- merge(x=tmp.data, y=ind.aggr.qcs, by.x='library_id', by.y='sample.id', all=TRUE)
# Redefine sample ID.
tmp.data[, sample.id:=paste0('Covid-CD4_', facs.batch, '_', ifelse(test=as.integer(stim.time)==6, yes='06', no='24'), '_', virus, '_', donors.no, '_GEx')]
tmp.data[, sample.id:=str_replace(string=sample.id, pattern='SARS-CoV-2', replacement='SARS')]
# Final formatting.
cols.order <- c('Estimated Number of Cells',	'Mean Reads per Cell',	'Median Genes per Cell',	'Number of Reads',	'Valid Barcodes',	'Fraction Reads in Cells',	'Total Genes Detected',	'Median UMI Counts per Cell', 'Reads Mapped to Genome',	paste0('Reads Mapped Confidently to ', c('Genome',	'Intergenic Regions', 'Intronic Regions',	'Exonic Regions',	'Transcriptome')), 'Reads Mapped Antisense to Gene', 'Sequencing Saturation',	'Q30 Bases in Barcode',	'Q30 Bases in RNA Read',	'Q30 Bases in RNA Read 2',	'Q30 Bases in Sample Index',	'Q30 Bases in UMI', colnames(ind.aggr.qcs)[2:length(ind.aggr.qcs)])
# all(cols.order %in% colnames(tmp.data))
cols.order <- c(setdiff(x=colnames(tmp.data), y=cols.order), cols.order)
tmp.data <- tmp.data[, ..cols.order]
setorderv(x=tmp.data, cols=c('stim.time', 'virus', 'seq.batch', 'sample.id'), order=c(-1, -1, 1, 1))
# Output table.
tmp.file.name <- paste0(table.4.path, '/QCSumm_GEx.csv')
fwrite(file=tmp.file.name, x=tmp.data, na=NA, quote=TRUE)

# @ HTO summary
# Define metric summary files.
uniq.sample.ids <- unique(sample.ids[, .(sample.id, hto.tag)])
for.files <- uniq.sample.ids[, hto.tag]
for.files <- str_replace(string=for.files, pattern='/mnt/BioAdHoc/Groups/vd-vijay/vfajardo/COVID-19/deconvolution/HTO_based/', replacement='')
tmp.dates <- str_extract(string=for.files, pattern='^\\d{2}-\\d{2}-\\d{4}')
for.files <- str_replace(string=for.files, pattern='^\\d{2}-\\d{2}-\\d{4}', replacement=paste0(tmp.dates, '/count'))
for.files <- str_extract(string=for.files, pattern='^([^/]+/){4}')
to.track <- paste0('/mnt/BioAdHoc/Groups/vd-vijay/vfajardo/sequencing_data/', for.files, '/outs/metrics_summary.csv')
sample.qc.files <- to.track[file.exists(to.track)]
# Load data as a single list keepin track of ID.
tmp.data <- lapply(X=sample.qc.files, FUN=fread)
names(tmp.data) <- uniq.sample.ids[file.exists(to.track), sample.id]
tmp.data <- rbindlist(l=tmp.data, use.names=TRUE, idcol='sample.id')
# Merge metadata with metrics.
tmp.data <- merge(x=lane.info, y=tmp.data, by='sample.id')
# Redefine sample ID.
tmp.data[, sample.id:=paste0('Covid-CD4_', facs.batch, '_', ifelse(test=as.integer(stim.time)==6, yes='06', no='24'), '_', virus, '_', donors.no, '_CITE')]
tmp.data[, sample.id:=str_replace(string=sample.id, pattern='SARS-CoV-2', replacement='SARS')]
# Final formatting.
colnames(tmp.data) <- str_replace(string=colnames(tmp.data), pattern='^Antibody:\\s', replacement='')
cols.order <- c('Mean Reads per Cell',	'Valid Barcodes', 'Fraction Antibody Reads',	'Fraction Antibody Reads Usable',	'Antibody Reads Usable per Cell',	'Fraction Reads in Barcodes with High UMI Counts',	'Fraction Unrecognized Antibody',	'Antibody Reads in Cells',	'Median UMIs per Cell (summed over all recognized antibody barcodes)', 'Sequencing Saturation',	'Q30 Bases in Barcode',	'Q30 Bases in Antibody Read',	'Q30 Bases in Sample Index',	'Q30 Bases in UMI')
# all(cols.order %in% colnames(tmp.data))
cols.order <- c(setdiff(x=colnames(tmp.data), y=cols.order), cols.order)
tmp.data <- tmp.data[, ..cols.order]
setorderv(x=tmp.data, cols=c('stim.time', 'virus', 'seq.batch', 'sample.id'), order=c(-1, -1, 1, 1))
# Output table.
tmp.file.name <- paste0(table.4.path, '/QCSumm_HTO.csv')
fwrite(file=tmp.file.name, x=tmp.data, na=NA, quote=TRUE)

# @ TCR summary
# Define metric summary files.
uniq.sample.ids <- unique(sample.ids[, .(sample.id, vdj.path)])
sample.qc.files <- paste0(uniq.sample.ids[, vdj.path], '/metrics_summary.csv')
to.check <- all(file.exists(sample.qc.files))
# Load data as a single list keepin track of ID.
tmp.data <- lapply(X=sample.qc.files, FUN=fread)
names(tmp.data) <- uniq.sample.ids[, sample.id]
tmp.data <- rbindlist(l=tmp.data, use.names=TRUE, idcol='sample.id', fill=TRUE)
# Merge metadata with metrics.
tmp.data <- merge(x=lane.info, y=tmp.data, by='sample.id')
# Redefine sample ID.
tmp.data[, sample.id:=paste0('Covid-CD4_', facs.batch, '_', ifelse(test=as.integer(stim.time)==6, yes='06', no='24'), '_', virus, '_', donors.no, '_TCR')]
tmp.data[, sample.id:=str_replace(string=sample.id, pattern='SARS-CoV-2', replacement='SARS')]
# Final formatting.
cols.order <- c('Estimated Number of Cells',	'Mean Read Pairs per Cell',	'Number of Cells With Productive V-J Spanning Pair',	'Number of Read Pairs',	'Valid Barcodes',	'Fraction Reads in Cells',	'Median TRA UMIs per Cell',	'Median TRB UMIs per Cell', 'Reads Mapped to Any V(D)J Gene',	'Reads Mapped to TRA',	'Reads Mapped to TRB',	'Mean Used Read Pairs per Cell', 'Q30 Bases in Barcode',	'Q30 Bases in RNA Read 1',	'Q30 Bases in RNA Read 2',	'Q30 Bases in Sample Index',	'Q30 Bases in UMI', 'Cells With Productive V-J Spanning Pair',	'Cells With Productive V-J Spanning (TRA, TRB) Pair',	'Paired Clonotype Diversity',	'Cells With TRA Contig',	'Cells With TRB Contig',	'Cells With CDR3-annotated TRA Contig',	'Cells With CDR3-annotated TRB Contig',	'Cells With V-J Spanning TRA Contig',	'Cells With V-J Spanning TRB Contig',	'Cells With Productive TRA Contig',	'Cells With Productive TRB Contig')
# all(cols.order %in% colnames(tmp.data))
cols.order <- c(setdiff(x=colnames(tmp.data), y=cols.order), cols.order)
tmp.data <- tmp.data[, ..cols.order]
setorderv(x=tmp.data, cols=c('stim.time', 'virus', 'seq.batch', 'sample.id'), order=c(-1, -1, 1, 1))
# Output table.
tmp.file.name <- paste0(table.4.path, '/QCSumm_TCR.csv')
fwrite(file=tmp.file.name, x=tmp.data, na=NA, quote=TRUE)

# ---> Cell counts per condition per donor.
# @ Obtain each table.
sep.conv.tables <- lapply(X=table.4.objs, FUN=get.sep.conv.table)
names(sep.conv.tables) <- table.4.objs
# @ Output each table.
lapply(X=names(sep.conv.tables), FUN=function(table.obj){
  # Retrieve data.
  tmp.data <- sep.conv.tables[[table.obj]]
  # Clusters order.
  cluster.cols <- grep(x=colnames(tmp.data), pattern='^\\d+$', value=TRUE)
  other.cols <- setdiff(x=colnames(tmp.data), y=cluster.cols)
  cols.order <- c(other.cols, mixedsort(cluster.cols))
  tmp.data <- tmp.data[, ..cols.order]
  # Sort
  setorderv(x=tmp.data, cols=c('Virus peptide pool', 'Hospitalization status', 'Sample collention timepoint', 'Donor'), order=c(-1, 1, 1, 1))
  # Output
  tmp.file.name <- paste0(table.4.path, '/CellCountsPerCond_Set-', table.obj, '.csv')
  fwrite(file=tmp.file.name, x=tmp.data, na=NA, quote=FALSE)
})

# ---> Clusters' markers.
# Output each table.
for(table.obj in table.4.objs){
  output.markers(obj.of.int=table.obj, output.path=table.4.path)
}

# ---> TCR extended information.
# Custom cross-reactivity info and phase tag info availability per dataset.
sets.cr.info <- rep(x=TRUE, times=length(table.4.objs)); names(sets.cr.info) <- table.4.objs
sets.phase.info <- rep(x=FALSE, times=length(table.4.objs)); names(sets.phase.info) <- table.4.objs
# Do process per dataset.
for(table.obj in table.4.objs){
  cat('Process for object:', table.obj, '\n')
  # Obtain vdj gene usage information.
  vdj.gene.info <- get.vdj.info(obj.of.int=table.obj)
  # Obtain data separately for clusters and virus reactivities and then merge.
  cr.info <- sets.cr.info[table.obj]; phase.info <- sets.phase.info[table.obj]
  tcr.by.cluster <- get.tcr.by.tag.table(obj.of.int=table.obj, vdj.gene.info=vdj.gene.info, tag.of.int=clusts.lab[table.obj], phase.info=phase.info, cr.info=cr.info)
  tcr.by.virus <- get.tcr.by.tag.table(obj.of.int=table.obj, vdj.gene.info=vdj.gene.info, tag.of.int='virus.tag', phase.info=phase.info, cr.info=cr.info)
  to.check <- tcr.by.cluster[, .N]==tcr.by.virus[, .N] & tcr.by.cluster[, uniqueN(clonotype)]==tcr.by.virus[, uniqueN(clonotype)]
  if(!to.check) stop('Merging of virus and cluster data unsuccessful, error type 1.\n')
  tmp.data <- merge(x=tcr.by.virus, y=tcr.by.cluster, by=c('clonotype', 'trb.aa', 'tra.aa', 'trb.nt', 'tra.nt', 'trb.v', 'tra.v', 'trb.j', 'tra.j', 'pr.tag', 'donor', 'phase', 'global.size', 'local.size'), suffixes=c('.virus', '.cluster'))
  to.check <- tmp.data[, .N]==tcr.by.cluster[, .N] & tmp.data[, .N]==tcr.by.virus[, .N]
  if(!to.check) stop('Merging of virus and cluster data unsuccessful, error type 2.\n')
  # Final formatting.
  tmp.data[, to.sort:=as.integer(str_replace(string=clonotype, pattern='clonotype', replacement=''))]
  setorderv(x=tmp.data, cols=c('phase', 'to.sort', 'donor'), order=c(-1, 1, 1)); tmp.data[, to.sort:=NULL]
  if(!cr.info){
    tmp.data[, pr.tag:=NULL]
  }else{
    tmp.data[, pr.tag:=ifelse(test=pr.tag=='pR', yes='Cross-reactive', no='Other')]
  }
  if(!phase.info){
    tmp.data[, phase:=NULL]
  }else{
    tmp.data[, phase:=ifelse(test=phase=='acute', yes='Acute', no='Convalescent')]
  }
  # Output.
  tmp.file.name <- paste0(table.4.path, '/TCRExtendedData_Set-', table.obj, '.csv')
  fwrite(file=tmp.file.name, x=tmp.data, na=NA, quote=FALSE)
}


############    -----------------------------------------    ############
### --------------------- On further validation --------------------- ###
### ----------------------- in healthy donors ----------------------- ###
############    -----------------------------------------    ############

# ----> Output directory.
table.5.path <- paste0(reports.path, '/table_on_further_val')
if(!dir.exists(table.5.path)) dir.create(table.5.path)

# ---> General definitions.
# Seurat objects of interest.
table.5.objs <- c('HLTY-6-ExpB-SFC', 'HLTY-6-ExpC-ALL')

# ---> Aggregation QC summaries.
# @ Load summaries.
aggr.qc.files <- paste0(gen.data.path, '/AggrSummary_', table.5.objs, '.json')
aggr.qcs <- lapply(X=aggr.qc.files, FUN=function(tmp.file) fromJSON(file=tmp.file))
names(aggr.qcs) <- table.5.objs
# @ Get relevant data in a tidy format per aggregation of interest.
aggr.qcs <- lapply(X=names(aggr.qcs), FUN=process.aggr.json)
names(aggr.qcs) <- table.5.objs
# @ Get general aggregation metrics only.
tmp.data <- lapply(X=aggr.qcs, FUN=function(x) return(x[['gen.qc']]))
tmp.data <- rbindlist(l=tmp.data, use.names=TRUE, idcol='aggr')
tmp.data <- spread(data=tmp.data, key=group, value=value)
cols.order <- c('aggr', 'pre.total', 'post.total', 'pre.mean', 'post.mean')
tmp.data <- tmp.data[, ..cols.order]
colnames(tmp.data) <- c(
  'Aggregation ID',
  'Pre-Normalization Total Number of Reads in unstimulated data',
  'Post-Normalization Total Number of Reads in unstimulated data',
  'Pre-Normalization Mean Reads per Cell in unstimulated data',
  'Post-Normalization Mean Reads per Cell in unstimulated data'
)
# @ Output general aggregation metrics only.
tmp.file.name <- paste0(table.5.path, '/QCSumm_Aggr.csv')
fwrite(file=tmp.file.name, x=tmp.data, na=NA, quote=TRUE)
# @ Merge per-library metrics.
ind.aggr.qcs <- lapply(X=aggr.qcs, FUN=function(x) return(x[['libs.qc']]))
ind.aggr.qcs <- Reduce(x=ind.aggr.qcs, f=function(x, y){
  to.output <- merge(x=x, y=y, by='sample.id', all=TRUE)
})
colnames(ind.aggr.qcs) <- c('sample.id', paste0(c('cmb.reads.per.bc.', 'frac.reads.kept.', 'raw.reads.per.bc.'), rep(x=names(aggr.qcs), each=3)))

# ---> Per-library QC summaries.
# @ General QC information for GEx and HTO data.
sample.ids <- lapply(X=table.5.objs, FUN=list.samples.info)
sample.ids <- unique(rbindlist(l=sample.ids, use.names=TRUE, fill=TRUE))
lane.info <- unique(sample.ids[, .(sample.id, facs.batch=orig.chrom.batch.tag, seq.batch=seq.batch.tag)])

# @ GEx summries.
# Define metric summary files.
uniq.sample.ids <- unique(sample.ids[, sample.id])
sample.qc.files <- unique(paste0(uniq.sample.ids, '/metrics_summary.csv'))
to.check <- all(file.exists(sample.qc.files))
# Load data as a single list keepin track of ID.
tmp.data <- lapply(X=sample.qc.files, FUN=fread)
names(tmp.data) <- uniq.sample.ids
tmp.data <- rbindlist(l=tmp.data, use.names=TRUE, idcol='sample.id')
# Merge metadata with metrics.
tmp.data <- merge(x=lane.info, y=tmp.data, by='sample.id')
tmp.data <- merge(x=tmp.data, y=unique(sample.ids[, .(sample.id, library_id)]), by='sample.id')
# Include per-library aggr information.
tmp.data <- merge(x=tmp.data, y=ind.aggr.qcs, by.x='library_id', by.y='sample.id', all=TRUE)
# Redefine sample ID.
tmp.data[, sample.id:=paste0('Covid-CD4_', facs.batch, '_', ifelse(test=facs.batch!='BeMe09', yes='06', no='00'), '_', ifelse(test=facs.batch!='BeMe09', yes='SFC', no='None'), '_', '2', '_GEx')]
# Final formatting.
cols.order <- c('Estimated Number of Cells',	'Mean Reads per Cell',	'Median Genes per Cell',	'Number of Reads',	'Valid Barcodes',	'Fraction Reads in Cells',	'Total Genes Detected',	'Median UMI Counts per Cell', 'Reads Mapped to Genome',	paste0('Reads Mapped Confidently to ', c('Genome',	'Intergenic Regions', 'Intronic Regions',	'Exonic Regions',	'Transcriptome')), 'Reads Mapped Antisense to Gene', 'Sequencing Saturation',	'Q30 Bases in Barcode',	'Q30 Bases in RNA Read',	'Q30 Bases in RNA Read 2',	'Q30 Bases in Sample Index',	'Q30 Bases in UMI', colnames(ind.aggr.qcs)[2:length(ind.aggr.qcs)])
# all(cols.order %in% colnames(tmp.data))
cols.order <- c(setdiff(x=colnames(tmp.data), y=cols.order), cols.order)
tmp.data <- tmp.data[, ..cols.order]
setorderv(x=tmp.data, cols=c('seq.batch', 'sample.id'), order=c(1, 1))
# Output table.
tmp.file.name <- paste0(table.5.path, '/QCSumm_GEx.csv')
fwrite(file=tmp.file.name, x=tmp.data, na=NA, quote=TRUE)

# @ HTO summary
# Define metric summary files.
uniq.sample.ids <- unique(sample.ids[, .(sample.id, hto.tag)])
for.files <- uniq.sample.ids[, hto.tag]
for.files <- str_replace(string=for.files, pattern='/mnt/BioAdHoc/Groups/vd-vijay/vfajardo/COVID-19/deconvolution/HTO_based/', replacement='')
tmp.dates <- str_extract(string=for.files, pattern='^\\d{2}-\\d{2}-\\d{4}')
for.files <- str_replace(string=for.files, pattern='^\\d{2}-\\d{2}-\\d{4}', replacement=paste0(tmp.dates, '/count'))
for.files <- str_extract(string=for.files, pattern='^([^/]+/){4}')
to.track <- paste0('/mnt/BioAdHoc/Groups/vd-vijay/vfajardo/sequencing_data/', for.files, '/outs/metrics_summary.csv')
sample.qc.files <- to.track[file.exists(to.track)]
# Load data as a single list keepin track of ID.
tmp.data <- lapply(X=sample.qc.files, FUN=fread)
names(tmp.data) <- uniq.sample.ids[file.exists(to.track), sample.id]
tmp.data <- rbindlist(l=tmp.data, use.names=TRUE, idcol='sample.id')
# Merge metadata with metrics.
tmp.data <- merge(x=lane.info, y=tmp.data, by='sample.id')
# Redefine sample ID.
tmp.data[, sample.id:=paste0('Covid-CD4_', facs.batch, '_', ifelse(test=facs.batch!='BeMe09', yes='06', no='00'), '_', ifelse(test=facs.batch!='BeMe09', yes='SFC', no='None'), '_', '2', '_CITE')]
# Final formatting.
colnames(tmp.data) <- str_replace(string=colnames(tmp.data), pattern='^Antibody:\\s', replacement='')
cols.order <- c('Mean Reads per Cell',	'Valid Barcodes', 'Fraction Antibody Reads',	'Fraction Antibody Reads Usable',	'Antibody Reads Usable per Cell',	'Fraction Reads in Barcodes with High UMI Counts',	'Fraction Unrecognized Antibody',	'Antibody Reads in Cells',	'Median UMIs per Cell (summed over all recognized antibody barcodes)', 'Sequencing Saturation',	'Q30 Bases in Barcode',	'Q30 Bases in Antibody Read',	'Q30 Bases in Sample Index',	'Q30 Bases in UMI')
# all(cols.order %in% colnames(tmp.data))
cols.order <- c(setdiff(x=colnames(tmp.data), y=cols.order), cols.order)
tmp.data <- tmp.data[, ..cols.order]
setorderv(x=tmp.data, cols=c('seq.batch', 'sample.id'), order=c(1, 1))
# Output table.
tmp.file.name <- paste0(table.5.path, '/QCSumm_HTO.csv')
fwrite(file=tmp.file.name, x=tmp.data, na=NA, quote=TRUE)

# @ TCR summary
# Define metric summary files.
uniq.sample.ids <- unique(sample.ids[, .(sample.id, vdj.path)])
sample.qc.files <- paste0(uniq.sample.ids[, vdj.path], '/metrics_summary.csv')
to.check <- all(file.exists(sample.qc.files))
# Load data as a single list keepin track of ID.
tmp.data <- lapply(X=sample.qc.files, FUN=fread)
names(tmp.data) <- uniq.sample.ids[, sample.id]
tmp.data <- rbindlist(l=tmp.data, use.names=TRUE, idcol='sample.id', fill=TRUE)
# Merge metadata with metrics.
tmp.data <- merge(x=lane.info, y=tmp.data, by='sample.id')
# Redefine sample ID.
tmp.data[, sample.id:=paste0('Covid-CD4_', facs.batch, '_', ifelse(test=facs.batch!='BeMe09', yes='06', no='00'), '_', ifelse(test=facs.batch!='BeMe09', yes='SFC', no='None'), '_', '2', '_TCR')]
# Final formatting.
cols.order <- c('Estimated Number of Cells',	'Mean Read Pairs per Cell',	'Number of Cells With Productive V-J Spanning Pair',	'Number of Read Pairs',	'Valid Barcodes',	'Fraction Reads in Cells',	'Median TRA UMIs per Cell',	'Median TRB UMIs per Cell', 'Reads Mapped to Any V(D)J Gene',	'Reads Mapped to TRA',	'Reads Mapped to TRB',	'Mean Used Read Pairs per Cell', 'Q30 Bases in Barcode',	'Q30 Bases in RNA Read 1',	'Q30 Bases in RNA Read 2',	'Q30 Bases in Sample Index',	'Q30 Bases in UMI', 'Cells With Productive V-J Spanning Pair',	'Cells With Productive V-J Spanning (TRA, TRB) Pair',	'Paired Clonotype Diversity',	'Cells With TRA Contig',	'Cells With TRB Contig',	'Cells With CDR3-annotated TRA Contig',	'Cells With CDR3-annotated TRB Contig',	'Cells With V-J Spanning TRA Contig',	'Cells With V-J Spanning TRB Contig',	'Cells With Productive TRA Contig',	'Cells With Productive TRB Contig')
# all(cols.order %in% colnames(tmp.data))
cols.order <- c(setdiff(x=colnames(tmp.data), y=cols.order), cols.order)
tmp.data <- tmp.data[, ..cols.order]
setorderv(x=tmp.data, cols=c('seq.batch', 'sample.id'), order=c(1, 1))
# Output table.
tmp.file.name <- paste0(table.5.path, '/QCSumm_TCR.csv')
fwrite(file=tmp.file.name, x=tmp.data, na=NA, quote=TRUE)

# ---> Cell counts per condition per donor.
# @ Obtain each table.
extra.tags <- c(`HLTY-6-ExpB-SFC`='arte.concentration.tag', `HLTY-6-ExpC-ALL`='cell.pop.tag')
sep.conv.tables <- lapply(X=table.5.objs, FUN=function(table.obj){
  get.sep.conv.table(table.obj, extra.tag=extra.tags[table.obj])
})
names(sep.conv.tables) <- table.5.objs
# @ Output each table.
lapply(X=names(sep.conv.tables), FUN=function(table.obj){
  # Retrieve data.
  tmp.data <- sep.conv.tables[[table.obj]]
  # Clusters order.
  cluster.cols <- grep(x=colnames(tmp.data), pattern='^\\d+$', value=TRUE)
  other.cols <- setdiff(x=colnames(tmp.data), y=cluster.cols)
  cols.order <- c(other.cols, mixedsort(cluster.cols))
  tmp.data <- tmp.data[, ..cols.order]
  # Sort
  setorderv(x=tmp.data, cols=c('Virus peptide pool', 'Extra tag', 'Hospitalization status', 'Sample collention timepoint', 'Donor'), order=c(-1, 1, 1, 1, 1))
  # Output
  tmp.file.name <- paste0(table.5.path, '/CellCountsPerCond_Set-', table.obj, '.csv')
  fwrite(file=tmp.file.name, x=tmp.data, na=NA, quote=FALSE)
})

# ---> Clusters' markers.
# Output each table.
for(table.obj in table.5.objs){
  output.markers(obj.of.int=table.obj, output.path=table.5.path)
}

# ---> TCR extended information.
# Custom cross-reactivity info and phase tag info availability per dataset.
sets.cr.info <- c(TRUE, FALSE); names(sets.cr.info) <- table.5.objs
sets.phase.info <- rep(x=FALSE, times=length(table.5.objs)); names(sets.phase.info) <- table.5.objs
# Do process per dataset.
for(table.obj in table.5.objs){
  cat('Process for object:', table.obj, '\n')
  # Obtain vdj gene usage information.
  vdj.gene.info <- get.vdj.info(obj.of.int=table.obj)
  # Obtain data separately for clusters and virus reactivities and then merge.
  cr.info <- sets.cr.info[table.obj]; phase.info <- sets.phase.info[table.obj]
  tcr.by.cluster <- get.tcr.by.tag.table(obj.of.int=table.obj, vdj.gene.info=vdj.gene.info, tag.of.int=clusts.lab[table.obj], phase.info=phase.info, cr.info=cr.info)
  second.tag <- if(table.obj=='HLTY-6-ExpC-ALL') 'cell.pop.tag' else 'virus.arte.con.tag'
  tcr.by.second <- get.tcr.by.tag.table(obj.of.int=table.obj, vdj.gene.info=vdj.gene.info, tag.of.int=second.tag, phase.info=phase.info, cr.info=cr.info)
  to.check <- tcr.by.cluster[, .N]==tcr.by.second[, .N] & tcr.by.cluster[, uniqueN(clonotype)]==tcr.by.second[, uniqueN(clonotype)]
  if(!to.check) stop('Merging of virus and cluster data unsuccessful, error type 1.\n')
  tmp.data <- merge(x=tcr.by.second, y=tcr.by.cluster, by=c('clonotype', 'trb.aa', 'tra.aa', 'trb.nt', 'tra.nt', 'trb.v', 'tra.v', 'trb.j', 'tra.j', 'pr.tag', 'donor', 'phase', 'global.size', 'local.size'), suffixes=c('.virus', '.cluster'))
  to.check <- tmp.data[, .N]==tcr.by.cluster[, .N] & tmp.data[, .N]==tcr.by.second[, .N]
  if(!to.check) stop('Merging of virus and cluster data unsuccessful, error type 2.\n')
  # Final formatting.
  tmp.data[, to.sort:=as.integer(str_replace(string=clonotype, pattern='clonotype', replacement=''))]
  setorderv(x=tmp.data, cols=c('phase', 'to.sort', 'donor'), order=c(-1, 1, 1)); tmp.data[, to.sort:=NULL]
  if(!cr.info){
    tmp.data[, pr.tag:=NULL]
  }else{
    tmp.data[, pr.tag:=ifelse(test=pr.tag=='pR', yes='Cross-reactive', no='Other')]
  }
  if(!phase.info){
    tmp.data[, phase:=NULL]
  }else{
    tmp.data[, phase:=ifelse(test=phase=='acute', yes='Acute', no='Convalescent')]
  }
  # Output.
  tmp.file.name <- paste0(table.5.path, '/TCRExtendedData_Set-', table.obj, '.csv')
  fwrite(file=tmp.file.name, x=tmp.data, na=NA, quote=FALSE)
}
