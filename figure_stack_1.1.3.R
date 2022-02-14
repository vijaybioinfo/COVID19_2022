############    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^    ############
############    -----   Preliminary figure stack    -----    ############
############    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^    ############

# ---> About the script.
# Version: 1
# Subversion: 0
# Script to get all bioinformatics-based panel figures for any of the final main and supplementary figures for the paper known as COVID CD4 II ('follow-up') from the general COVID project of Vijay's lab.


############    -----------------------------------------    ############
### -------------------------- Description -------------------------- ###
############    -----------------------------------------    ############

# Copy description from DICE's CD4STIM project's script (figure 1), but if possible, include in a different file.


############    -----------------------------------------    ############
### --------------------------- Libraries --------------------------- ###
############    -----------------------------------------    ############
library(Seurat)
library(ggplot2)
library(ggpubr)
library(VennDiagram)
library(ggnewscale)
library(pheatmap)
library(stringr)
library(data.table)
library(tidyr)
library(gtools)
library(mefa)
library(UpSetR)
library(fgsea)
# library(NMF)


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
this.figure <- 'figure_stack_v1'
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
# @ For plotting purposes.
blank.complement.1 <- theme(line=element_blank(), text=element_blank(), panel.background=element_blank(), panel.border=element_blank(), legend.position='none', axis.ticks=element_line(size=1.2), axis.line=element_line(size=1.2), axis.ticks.length=unit(0.4, "cm")) # Default blank.
blank.complement.1.1 <- theme(line=element_blank(), text=element_blank(), panel.background=element_blank(), panel.border=element_blank(), legend.position='none', axis.ticks=element_line(size=1.2), axis.line=element_line(size=1.2), axis.ticks.length=unit(0.4, "cm"), axis.ticks.x=element_blank()) # Alternative to 1 with no ticks in the x axis.
blank.complement.2 <- theme(line=element_blank(), text=element_blank(), panel.background=element_blank(), panel.border=element_blank(), legend.position='none', axis.ticks=element_line(size=0.4), axis.line=element_line(size=0.4), axis.ticks.length=unit(0.15, "cm")) # For the dot plot.
blank.complement.3 <- theme(line=element_blank(), text=element_blank(), panel.background=element_blank(), panel.border=element_blank(), legend.position='none', axis.ticks=element_blank(), axis.line=element_blank()) # For UMAP plots (axes are not drawn).
blank.complement.4 <- theme(line=element_blank(), text=element_blank(), panel.background=element_blank(), panel.border=element_blank(), legend.position='none') # For donut plots.
# @ Gene signature definitions.
cd4.ctl.sign <- 'cell.cytotoxicity.patil.score'
# @ For tag analysis.
min.no.cells.per.group <- 50 # @ Internal parameter for tag specific analysis (part of figure ?).
colors.table <- NULL
# ---> Color defintions.
# @ Color Qper cluster per dataset.
clusts.cols <- list(
  `PT-6-ALL-SF`=c(
    `0`='#a6b0e3', # TFH
    `1`='#cc0000', # CD4-CTL
    `2`='#e1ba7f', # TH17/polyfunctional
    `3`='#ffd700', # THIFNR
    `4`='#754600', # TH1/polyfunctional
    `5`='#b95872', # Cytotoxic TFH
    `6`='#cdc495' # Proliferative
  ),
  `PT-24-ALL-SF`=c(
    `0`='#319272', # TREG
    `1`='#cc0000', # CD4-CTL
    `2`='#a6b0e3', # TFH
    `3`='#83bdaa', # TREG
    `4`='#ff5656' # CD4-CTL
  ),
  `PT-6-Conv-SFCE`=c(
    `0`='#cc0000',
    `1`='#a6b0e3',
    `2`='#ffd700'
  ),
  `HLTY-6-ExpA-FCE`=c(
    `0`='#808080', # Unimportant
    `1`='#cc0000', # CD4-CTL
    `2`='#b3b3b3', # Unimportant
    `3`='#5d5d5d' # Unimportant
  ),
  `HLTY-6-ExpB-SFC`=c(
    `0`='#a6b0e3', # TFH
    `1`='#cc0000', # NK-like CD4-CTL.
    `2`='#800000', # CD4-CTL
    `3`='#ff5656', # CD4-CTL
    `4`='#e1ba7f', # TH17/polyfunctional # NOTE: This may be a TFH population because it's got CXCL10 and IL32 as markers.
    `5`='#ffd700', # THIFNR
    `6`='#ba0000' # CD4-CTL
  ),
  `HLTY-6-ExpC-ALL`=c(
    `0`='#00bfff', # Naive CD4 T cells (markers: CD27, IL7R, and TCF7)
    `1`='#cc0000', # CD4-CTL.
    `2`='#800000', # CD4-CTL
    `3`='#ff5656', # CD4-CTL
    `4`='#1a9fff', # Naive CD4 T cells (markers: CCR7, SELL, and TCF7) or more like TCM.
    `5`='#319272' # Memory TREGs
  ),
  `HLTY-6-ExpD-S3`=c(
    `0`='#808080', # Unimportant 1
    `1`='#b3b3b3', # Unimportant 2
    `2`='#cc0000' # CD4-CTL
  )
)
# @ Other terms' colors.
virus.cols <- c('SARS-CoV-2'='#FE0002', 'HCoV'='#a76400', 'FLU'='#00947F', 'HCoV'='#EEA101', 'EBV'='#005746', 'CMV'='#42a2cc', 'HIV'='#720000', 'FLU|CMV'='#219ba6', 'FLU|CMV|EBV'='#b39700')
pr.groups.cols <- c('pR'='#cccc00', 'non-pR'='#0074e6')
# pr.cluster.class.cols <- c(`Inter-cluster`='#cc9900', `Intra-cluster`='#cccc00')
blood.phases.cols <- c('Acute'='#cb6400', 'Convalescent'='#0067cb')
sharing.cols <- c('Shared'='#cccc00', 'Unique'='#0074e6')
slamf7.cols <- c(Positive='#cc0000', Negative='#00bfff')
# pr.degree.cols <- c('1'='#00008b', '2'='#556b2f', '3'='#ffff00', '4'='#8b0000')
concentration.cols <- c(`0.001`='#d3d3d3', `0.01`='#a9a9a9', `0.1`='#808080', `1`='#000000')
signatures.col.scale <- c('#ffffff', '#ffffe0', '#ffffad', '#ffdf32', '#ff9a00', '#ff5a00', '#ff5719','#EE0000','#b30000', '#670000')
hosp.cols <- c(No='#9AD34C', Yes='#F45C54')
# signatures.alt.col.scale <- c('#021893', '#3537ae', '#740699', '#b70b0b', '#990623')
# featplot.col.scale <- c('#ffdf32', '#ff9a00', '#ff5a00', '#ff5719','#EE0000','#b30000', '#670000')
# ---> Path definitions.
final.figs.path <- paste0('/mnt/BioAdHoc/Groups/vd-vijay/vfajardo/COVID-19/paper_developments/', this.prj, '/final_figures')
gen.data.path <- paste0(final.figs.path, '/data_to_start')
this.fig.path <- paste0(final.figs.path, '/', this.figure)
reports.path <- paste0(this.fig.path, '/', this.figure, '_panels')
# ---> File definitions.
# Seurat objects.
objs.files.list <- paste0(gen.data.path, '/', 'SeuratObj_', dataset.labs, '_Processed.RDS')
# Gene sets' files (for FGSE analyses).
module.feats.dir <- paste0(gen.data.path, '/module_signatures_for_fgsea')
# Specific files.
#     SARS-CoV-2-reactive cell count per donor.
sr.cell.counts.file <- paste0(gen.data.path, '/SortedCellsPerMillionPBMCs.csv')
#     HLA typing resuls.
hla.data.file <- paste0(gen.data.path, '/DataFromHLATypingServiceProvider.csv')
#     SARS-CoV-2- and HCoV-reactive cells per 10e6 PBMCs.
hlty.hcov.prof.file <- paste0(gen.data.path, '/HltyDonor_HCoV_and_SARS_Cells.csv')
#     CD4-CTL profile of healthy donors.
hlty.cyt.prof.file <- paste0(gen.data.path, '/HltyDonor_CD4-CTL-and-Virus-ReactiveTCellProfile.csv')
#     Titration experiments data.
titration.exp.file <- paste0(gen.data.path, '/HltyDonor_TitrationExpData.csv')
#     DEA, HLTY-6-ExpD-S3, HCoV vs SARS-CoV-2 cells.
hlty.six.expd.s3.dea.file <- paste0(gen.data.path, '/HLTY-6-ExpD-S3_DEA_HCoV-vs-SARS.csv')
# ---> Check directories and files.
if(!all(dir.exists(gen.data.path), dir.exists(this.fig.path))) stop(paste0('Following paths must be already defined:\n', gen.data.path, '\n', this.fig.path))
if(!dir.exists(reports.path)) dir.create(reports.path)
essential.files <- c(
  lapply(X=objs.files.list, FUN=function(x) return(x)),
  hlty.six.expd.s3.dea.file
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

# ---> Specific files.
#     SARS-CoV-2-reactive cell count per donor.
sr.cell.counts <- fread(file=sr.cell.counts.file)
#     HLA typing resuls.
hla.data <- fread(file=hla.data.file)
#     SARS-CoV-2- and HCoV-reactive cells per 10e6 PBMCs.
hlty.hcov.prof <- fread(file=hlty.hcov.prof.file)
#     CD4-CTL profile of healthy donors.
hlty.cyt.prof <- fread(file=hlty.cyt.prof.file)
#     Titration experiments data.
titration.exp <- fread(file=titration.exp.file)
#     DEA, HLTY-6-ExpD-S3, HCoV vs SARS-CoV-2 cells.
hlty.six.expd.s3.dea <- fread(file=hlty.six.expd.s3.dea.file)

# ---> Gene sets' files (for FGSE analyses).
module.feats.files <- list.files(path=module.feats.dir, recursive=FALSE, include.dirs=FALSE, pattern='signature.csv$', full.names=TRUE)
modules.feats <- lapply(X=module.feats.files, FUN=function(tmp.file){
  tmp.feats <- read.csv(file=tmp.file, stringsAsFactors=FALSE)
  if(!'feature' %in% colnames(tmp.feats)) stop(paste0('File ', tmp.file, ' not appropriately defined. Column listing the gene features should be named \'feature\'.\n\n'))
  tmp.feats <- tmp.feats$feature
  # Replace underscores with points in the gene names just as Seurat does when one creates a seurat object.
  tmp.check <- any(str_detect(string=tmp.feats, pattern='_'))
  if(tmp.check){
    tmp.warning <- paste0('Feature names cannot have underscores (\'_\'), replacing with dashes (\'-\'). This for file: ', basename(tmp.file))
    warning(tmp.warning)
    tmp.feats <- str_replace_all(string=tmp.feats, pattern='_', replacement='-')
  }
  # Return.
  return(tmp.feats)
})
names(modules.feats) <- list.files(path=module.feats.dir, recursive=FALSE, include.dirs=FALSE, pattern='signature.csv$', full.names=FALSE)
names(modules.feats) <- str_replace(string=names(modules.feats), pattern='_signature.csv$', replacement='')
names(modules.feats) <- str_replace_all(string=names(modules.feats), pattern='_', replacement='.')
modules.names <- names(modules.feats)
names(modules.names) <- str_replace_all(string=modules.names, pattern='\\.', replacement='_')


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

# ---> HLTY-6-EXPC-ALL
# @ Clusters definition.
# Remove clusters 6, 7 and 8 (less than 1% of the whole dataset).
cells.to.keep <- !as.character(srt.objs.list[[hlty.six.expc.all]]@meta.data[, clusts.lab[hlty.six.expb.sfc]])%in%c('6', '7', '8')
cells.to.keep <- Cells(srt.objs.list[[hlty.six.expc.all]])[cells.to.keep]
srt.objs.list[[hlty.six.expc.all]] <- subset(x=srt.objs.list[[hlty.six.expc.all]], cells=cells.to.keep)
# @ SLAMF7 population tag.
srt.objs.list[[hlty.six.expc.all]]@meta.data[, 'slamf7.pop.tag'] <- ifelse(test=srt.objs.list[[hlty.six.expc.all]]@meta.data[, 'cell.pop.tag']=='SLAMF7+', yes='Positive', no='Negative')

# ---> SARS-CoV-2- and HCoV-reactive cells per 10e6 PBMCs.
hlty.hcov.prof <- as.data.table(gather(data=hlty.hcov.prof, key='virus', value='cells', -`donor.id.tag`))
hlty.hcov.prof[, virus:=str_replace(string=virus, pattern='\\.cells$', replacement='')]
hlty.hcov.prof[, virus:=ifelse(test=virus=='hcov', yes='HCoV', no='SARS-CoV-2')]

# ---> Titration experiments data.
# Tidy data up.
titration.exp <- as.data.table(gather(data=titration.exp, `020`, `064`, `067`, `085`, key='donor.id', value='count'))
# Set concentration as a factor.
new.lvls <- as.character(mixedsort(titration.exp[, unique(concentration)], decreasing=FALSE))
titration.exp$concentration <- factor(x=titration.exp[, as.character(concentration)], levels=new.lvls)


############    -----------------------------------------    ############
### --------------------------- Figure 1 ---------------------------- ###
############    -----------------------------------------    ############

# ---> General definitions.
# Seurat object of interest.
main.obj.1 <- 'PT-6-ALL-SF'
sec.obj.1 <- 'PT-24-ALL-SF'
# Seurat object's meta data
meta.data <- cbind(
  data.table(cell=Cells(srt.objs.list[[main.obj.1]])),
  srt.objs.list[[main.obj.1]]@meta.data,
  srt.objs.list[[main.obj.1]]@reductions$umap@cell.embeddings
)
sec.meta.data <- cbind(
  data.table(cell=Cells(srt.objs.list[[sec.obj.1]])),
  srt.objs.list[[sec.obj.1]]@meta.data,
  srt.objs.list[[sec.obj.1]]@reductions$umap@cell.embeddings
)

# ---> Preliminary: Provide list of cohort for either dataset.
# Because of the issue: Age does not match with the original age.
donor.info.file <- '/mnt/BioAdHoc/Groups/vd-vijay/vfajardo/COVID-19/paper_developments/Batches-4-5-6-7/donors_metadata/GeneralMetaDataTable2ForCOVIDGenExBatches-4-5-6-7.1.1.csv'
donor.info <- fread(file=donor.info.file)
donor.info <- donor.info[gender.tag %in% c('Male', 'Female'),
  .(
    donor.id=paste0(
      'P', ifelse(test=as.integer(donor.alt.id.tag)>9, yes='', no='0'), donor.alt.id.tag
    ),
    age=age.tag,
    sex=gender.tag,
    ethnicity=population.tag,
    hosp.status=ifelse(test=hospitalization.tag==TRUE, yes='Yes', no='No')
  )
]
#
tmp.data.1 <- meta.data[, .(age=unique(age.tag), sex=unique(gender.tag), ethnicity=unique(population.tag), hosp.status=unique(hospitalization.tag), `6.hr`=TRUE), by=.(donor.id=as.character(donor.id.tag))]
tmp.data.2 <- sec.meta.data[, .(age=unique(age.tag), sex=unique(gender.tag), ethnicity=unique(population.tag), hosp.status=unique(hospitalization.tag), `24.hr`=TRUE), by=.(donor.id=as.character(donor.id.tag))]
tmp.data <- merge(x=tmp.data.1, y=tmp.data.2, by=c('donor.id', 'age', 'sex', 'ethnicity', 'hosp.status'), fill=FALSE, sort=TRUE, all=TRUE)
weird.donors <- tmp.data[, .N, by=donor.id][N>1, donor.id]# Sanity check
to.check <- tmp.data[, uniqueN(donor.id)]
tmp.data <- merge(x=tmp.data, y=donor.info, by=c('donor.id', 'age', 'sex', 'ethnicity', 'hosp.status'), fill=FALSE)
tmp.data[, uniqueN(donor.id)]==to.check
tmp.file.name <- paste0(fig.1.path, '/CohortDescription_CD4-FollowUp.csv')
fwrite(file=tmp.file.name, x=tmp.data, na=NA, quote=FALSE)
# TO DO: EVRYTHING IS FINE EXCEPT FOR THE AGE TAG. SOMETHING MAY GET SCREWED UP WHILE ANNOTATING WITH INTEGERS. THIS IS WORRYING BECAUSE MAY POTENTIALLY AFFECT THE RESULTS OF OTHER PROJECTS. STRANGELY, I DON'T REMEMBER SEEING THIS PROBLEM FOR ANY OTHER PROJECT: E.G., THE R24 PROJECT.

# ---> Preliminary. Stats for COVID-19 patients enrolled in our study.
tmp.data <- unique(rbind(unique(meta.data[!is.na(donor.id.tag), .(donor=donor.id.tag, hosp.tag=hospitalization.tag)]), unique(sec.meta.data[!is.na(donor.id.tag), .(donor=donor.id.tag, hosp.tag=hospitalization.tag)])))
unique.donors <- tmp.data[, sort(unique(as.character(donor)))]
# c('P04', 'P37') %in% meta.data[, as.character(donor.id.tag)]
sr.cell.counts[patient.id %in% unique.donors, range(interval.memory)]
sr.cell.counts[patient.id %in% unique.donors, median(interval.memory)]


### ------------------------- Main Figure 1 ------------------------- ###
# ----> Output directory.
fig.1.path <- paste0(reports.path, '/figure_1')
if(!dir.exists(fig.1.path)) dir.create(fig.1.path)
# ---> Comments: Main dataset for this section, PT-6-ALL-SF

# ---> Numbers needed for the study overview.
meta.data[virus.tag=='SARS-CoV-2', .N]
meta.data[virus.tag=='SARS-CoV-2' & !is.na(clonotype.tag), .N]
meta.data[virus.tag=='FLU', .N]
meta.data[virus.tag=='FLU' & !is.na(clonotype.tag), .N]
sec.meta.data[virus.tag=='SARS-CoV-2', .N]
sec.meta.data[virus.tag=='SARS-CoV-2' & !is.na(clonotype.tag), .N]
sec.meta.data[virus.tag=='FLU', .N]
sec.meta.data[virus.tag=='FLU' & !is.na(clonotype.tag), .N]

# ---> Clonotype sharing between virus reactivities.
tmp.cols <- virus.cols[c('SARS-CoV-2', 'FLU')]
depict.clone.sharing(meta.data=meta.data, tag.of.int='virus.tag', filter.tag='blood.draw.phase.tag', groups.to.filter='convalescent', reports.path=fig.1.path, file.preffix='CovidPt_TCR-Repertoire_ClonotypeSharingBetweenReactivities', fill.cols=tmp.cols)

# ---> Contrasting clonal expansion between crossreactive and non-crossreactive TCRs seen in the convalescent phase.
# Plot and output.
tmp.cols <- sharing.cols; names(tmp.cols) <- c('pR', 'non-pR')
srt.objs.list[[main.obj.1]]@meta.data[, 'log.clon.size.tag'] <- log2(srt.objs.list[[main.obj.1]]@meta.data[, 'clon.size.tag'])
tmp.ggplot <-
vln.plot(seurat.obj=srt.objs.list[[main.obj.1]], feature='log.clon.size.tag', groups.tag='pr.tag', na.rm=TRUE, groups.of.int=NULL, groups.label=NULL, filter.tags=c('blood.draw.phase.tag', 'virus.tag'), groups.to.filter=list(c('convalescent'), c('SARS-CoV-2')), keep=TRUE, feature.thold=1, color='median', vln.type='density', size.thold=0, file.name=NULL, adjust.val=2, trim.val=TRUE, plot.limits=NULL, add.points=FALSE, this.color.scale=tmp.cols)
tmp.ggplot <- tmp.ggplot + scale_x_continuous(expand=expansion(add=0), breaks=scales::pretty_breaks(n=3)) + scale_y_continuous(expand=expansion(add=0), breaks=scales::pretty_breaks(n=2), limits=c(0, 0.5))
publish.plot(tmp.ggplot=tmp.ggplot, output.path=fig.1.path, file.name='CovidPt_TCR-Repertoire_SizeBetween-pR-and-nonpR_ExpA', type='pdf', blank.comp=blank.complement.1, do.legend=TRUE, height=4, width=4)

# ---> Cross-reactive TCR representation per donor.
# Get data.
tmp.data <- meta.data[
  !is.na(clonotype.tag) & blood.draw.phase.tag=='convalescent' & virus.tag=='SARS-CoV-2',
  .(
    pr.prop=.SD[pr.tag=='pR', .N]/.N
  ),
  by=.(donor=donor.id.tag, hosp.status=hospitalization.tag)
]
setorderv(x=tmp.data, cols='pr.prop', order=-1)
tmp.data$donor <- factor(x=as.character(tmp.data$donor), levels=as.character(tmp.data$donor))
# Plot and output.
tmp.ggplot <- ggplot(data=tmp.data, aes(x=donor, y=pr.prop)) +
  geom_bar(stat='identity', fill='lightgray', width=0.7) +
  scale_y_continuous(limits=c(0, 1), expand=expansion(add=0), breaks=scales::pretty_breaks(n=3)) + scale_x_discrete(expand=expansion(add=0.5)) +
  labs(x='Donors ranked according to crossreactivity representation', y='Percent of total cells') +
  theme(axis.ticks.x=element_blank())
publish.plot(tmp.ggplot=tmp.ggplot, output.path=fig.1.path, file.name='CovidPt_TCR-Repertoire_pR-RepresentationPerDonor_ExpA', type='pdf', blank.comp=blank.complement.1, do.legend=TRUE, width=10)

# ---> pR profile for shared TCRs between phases.
# Obtain plots.
tmp.cols <- sharing.cols; names(tmp.cols) <- c('pR', 'non-pR')
tmp.ggplot.1 <- plot.props(seurat.obj=srt.objs.list[[main.obj.1]], plot.type='pie', groups.tag='pr.tag', na.rm=TRUE, groups.of.int=NULL, groups.label=NULL, filter.tags=c('blood.draw.phase.tag', 'virus.tag'), groups.to.filter=list('convalescent', 'SARS-CoV-2'), keep=TRUE, color.vals=tmp.cols)
# Output.
publish.plot(tmp.ggplot=tmp.ggplot.1, output.path=fig.1.path, file.name='CovidPt-TCR-Repertoire_pR-RepresentationOverall', blank.comp=blank.complement.4, do.legend=TRUE, height=5, width=5)
# Report proportions per pR status:
meta.data[blood.draw.phase.tag=='convalescent'&virus.tag=='SARS-CoV-2'&!is.na(clonotype.tag), .(cell.frac=.N/meta.data[blood.draw.phase.tag=='convalescent'&virus.tag=='SARS-CoV-2'&!is.na(clonotype.tag), .N]), by=pr.tag]
# pR 0.3891782
# non-pR 0.6108218
meta.data[blood.draw.phase.tag=='convalescent'&virus.tag=='SARS-CoV-2'&!is.na(clonotype.tag), .N]

# ---> Contrasting clonal expansion between crossreactive and non-crossreactive TCRs seen in the convalescent phase for selected donors.
# Plot and output.
donors.of.int <- c('P07', 'P43', 'P03')
tmp.cols <- sharing.cols; names(tmp.cols) <- c('pR', 'non-pR')
srt.objs.list[[main.obj.1]]@meta.data[, 'log.clon.size.tag'] <- log2(srt.objs.list[[main.obj.1]]@meta.data[, 'clon.size.tag'])
for(tmp.donor in donors.of.int){
  tmp.ggplot <- vln.plot(seurat.obj=srt.objs.list[[main.obj.1]], feature='log.clon.size.tag', groups.tag='pr.tag', na.rm=TRUE, groups.of.int=NULL, groups.label=NULL, filter.tags=c('blood.draw.phase.tag', 'donor.id.tag'), groups.to.filter=list(c('convalescent'), c(tmp.donor)), keep=TRUE, feature.thold=1, color='median', vln.type='density', size.thold=0, file.name=NULL, adjust.val=2, trim.val=TRUE, plot.limits=NULL, add.points=FALSE, this.color.scale=tmp.cols)
  tmp.ggplot <- if(tmp.donor=='P03') tmp.ggplot + scale_x_continuous(expand=expansion(add=0), breaks=c(5, 8)) + scale_y_continuous(expand=expansion(add=0), breaks=scales::pretty_breaks(n=2)) else tmp.ggplot + scale_x_continuous(expand=expansion(add=0), breaks=c(5, 10)) + scale_y_continuous(expand=expansion(add=0), breaks=scales::pretty_breaks(n=2))
  publish.plot(tmp.ggplot=tmp.ggplot, output.path=fig.1.path, file.name=paste0('CovidPt_TCR-Repertoire_SizeBetween-pR-and-nonpR_ExpA_ForDonor-', tmp.donor), type='pdf', blank.comp=blank.complement.1, do.legend=TRUE, height=4, width=4)
}

# ---> Comparison between SARS-reactive cell count and cross-reactive potential -crossreactive cell count- (at the donor level).
# Get data.
tmp.data.1 <- meta.data[
  !is.na(clonotype.tag) & blood.draw.phase.tag=='convalescent' & virus.tag=='SARS-CoV-2',
  .(
    pr.freq=.SD[pr.tag=='pR', .N],
    pr.prop=.SD[pr.tag=='pR', .N]/.N
  ),
  by=.(donor=donor.id.tag, hosp.status=hospitalization.tag)
]
tmp.data.2 <- sr.cell.counts[time.point=='6']
tmp.data.2[, hosp.status:=ifelse(test=hosp.status=='Hosp', yes='Yes', no='No')]
tmp.data <- merge(x=tmp.data.1, y=tmp.data.2, by.x=c('donor', 'hosp.status'), by.y=c('patient.id', 'hosp.status'))
# tmp.data[, all(hosp.status.x==hosp.status.y)] # Sanity check. Annotations may maintain stable across files.
# Plot and output.
tmp.ggplot <- ggplot(data=tmp.data, aes(x=cells.memory, y=pr.prop)) +
  geom_point(size=7, shape=1) + geom_smooth(formula=y~x, method='lm', se=FALSE, color='black') +
  scale_y_continuous(expand=expansion(add=0.015), breaks=scales::pretty_breaks(n=3)) + scale_x_continuous(expand=expansion(add=38), breaks=scales::pretty_breaks(n=3)) +
  labs(x='SARS-CoV-2-reactive cell count per 10e6 PBMC', y='Crossreactive cell fraction')
publish.plot(tmp.ggplot=tmp.ggplot, output.path=fig.1.path, file.name='CovidPt_TCR-Repertoire_SARSCellCount-vs-pR-RepresentationPerDonor_ExpA', type='pdf', blank.comp=blank.complement.1, do.legend=FALSE, stat.cor=TRUE)

# ---> pR profile on UMAP along clone size depiction
tmp.cells <- meta.data[blood.draw.phase.tag=='convalescent' & !is.na(clonotype.tag), cell]
meta.data[, tmp.clon.size.tag:=clon.size.tag]; meta.data[clon.size.tag>=1500, tmp.clon.size.tag:=1500]
tmp.ggplot <- ggplot(data=meta.data[cell %in% tmp.cells], aes(x=UMAP_1, y=UMAP_2, fill=pr.tag, size=tmp.clon.size.tag)) +
  geom_point(stroke=0.4, shape=21, col='black', alpha=0.8) +
  scale_fill_manual(values=pr.groups.cols) +
  scale_radius(breaks=c(10, 100, 1000), limits=c(1, meta.data[, max(tmp.clon.size.tag, na.rm=TRUE)]), range=c(3, 13)) +
  scale_y_continuous(expand=expansion(add=0)) + scale_x_continuous(expand=expansion(add=0)) +
  labs(x='UMAP 1', y='UMAP 2', col='')
# Output.
publish.plot(tmp.ggplot=tmp.ggplot, output.path=fig.1.path, file.name='CovidPt_SharingStatusOnUMAP', type='tiff', blank.comp=blank.complement.3, do.legend=TRUE)

# ---> UMAP depicting clusters.
# Plot.
tmp.cells <- meta.data[blood.draw.phase.tag=='convalescent', cell]
tmp.ggplot <- ggplot(data=meta.data[cell %in% tmp.cells], aes_string(x='UMAP_1', y='UMAP_2', col=clusts.lab[main.obj.1])) +
  geom_point(size=4, alpha=0.4, shape=19) +
  scale_color_manual(values=clusts.cols[[main.obj.1]]) +
  scale_y_continuous(expand=expansion(add=0)) + scale_x_continuous(expand=expansion(add=0)) +
  labs(x='UMAP 1', y='UMAP 2', col='Cluster')
# Output.
publish.plot(tmp.ggplot=tmp.ggplot, output.path=fig.1.path, file.name='CovidPt_ClustersOnUMAP', type='tiff', blank.comp=blank.complement.3, do.legend=TRUE)

# ---> Gene signatures.
# Plot.
tmp.cells <- meta.data[blood.draw.phase.tag=='convalescent', cell]
tmp.ggplot <- ggplot(data=meta.data[cell %in% tmp.cells], aes_string(x='UMAP_1', y='UMAP_2', col='cell.cytotoxicity.patil.score', fill='cell.cytotoxicity.patil.score')) +
  geom_point(size=4, alpha=0.8, shape=19) +
  scale_fill_gradientn(colors=signatures.col.scale) +
  scale_color_gradientn(colors=signatures.col.scale) +
  scale_y_continuous(expand=expansion(add=0)) + scale_x_continuous(expand=expansion(add=0)) +
  labs(x='UMAP 1', y='UMAP 2', col='Score')
# Output.
publish.plot(tmp.ggplot=tmp.ggplot, output.path=fig.1.path, file.name='CovidPt_Signature-Cytotoxicity', type='tiff', blank.comp=blank.complement.3, do.legend=TRUE)

# ---> CTL-related phenotype (fractions) per cross-reactive group (cross-reactive cells and others).
# Define CD4-CTL tag in seurat object.
tmp.corrs <- c('CD4-CTL', rep(x='Rest', times=6)); names(tmp.corrs) <- c('1', '0', 2:6)
tmp.data <- as.character(srt.objs.list[[main.obj.1]]@meta.data[, clusts.lab[main.obj.1]])
tmp.data <- tmp.corrs[tmp.data]
srt.objs.list[[main.obj.1]]@meta.data[, 'ctl.tag'] <- tmp.data
# seurat.obj@meta.data[, 'ctl.tag'] <- tmp.data
tmp.cols <- c('Rest'='#808080', 'CD4-CTL'='#cc0000')
# Obtain plots.
tmp.ggplot.1 <- plot.props(seurat.obj=srt.objs.list[[main.obj.1]], plot.type='donut', groups.tag='ctl.tag', na.rm=TRUE, groups.of.int=NULL, groups.label=NULL, filter.tags=c('pr.tag', 'blood.draw.phase.tag'), groups.to.filter=list('pR', 'convalescent'), keep=TRUE, color.vals=tmp.cols)
tmp.ggplot.2 <- plot.props(seurat.obj=srt.objs.list[[main.obj.1]], plot.type='donut', groups.tag='ctl.tag', na.rm=TRUE, groups.of.int=NULL, groups.label=NULL, filter.tags=c('pr.tag', 'blood.draw.phase.tag'), groups.to.filter=list('non-pR', 'convalescent'), keep=TRUE, color.vals=tmp.cols)
# Output.
publish.plot(tmp.ggplot=tmp.ggplot.1, output.path=fig.1.path, file.name='CovidPt_CD4-CTLs-and-pRAssessment-pR', blank.comp=blank.complement.4, do.legend=TRUE, height=5, width=5)
publish.plot(tmp.ggplot=tmp.ggplot.2, output.path=fig.1.path, file.name='CovidPt_CD4-CTLs-and-pRAssessment-nonpR', blank.comp=blank.complement.4, do.legend=TRUE, height=5, width=5)

# ---> pR profile (fraction) for CD4-CTLs and TFHs
# Define CD4-CTL tag in seurat object.
tmp.corrs <- c('CTL', rep(x='Rest', times=6)); names(tmp.corrs) <- c('1', '0', 2:6)
tmp.data <- as.character(srt.objs.list[[main.obj.1]]@meta.data[, clusts.lab[main.obj.1]])
tmp.data <- tmp.corrs[tmp.data]
srt.objs.list[[main.obj.1]]@meta.data[, 'ctl.tag'] <- tmp.data
# Define TFH tag in seurat object.
tmp.corrs <- c(rep(x='TFH', times=2), rep(x='Rest', times=5)); names(tmp.corrs) <- c('0', '5', 1:4, '6')
tmp.data <- as.character(srt.objs.list[[main.obj.1]]@meta.data[, clusts.lab[main.obj.1]])
tmp.data <- tmp.corrs[tmp.data]
srt.objs.list[[main.obj.1]]@meta.data[, 'tfh.tag'] <- tmp.data
# seurat.obj@meta.data[, 'ctl.tag'] <- tmp.data
tmp.cols <- sharing.cols; names(tmp.cols) <- c('pR', 'non-pR')
# Obtain plots.
tmp.ggplot.1 <- plot.props(seurat.obj=srt.objs.list[[main.obj.1]], plot.type='donut', groups.tag='pr.tag', na.rm=TRUE, groups.of.int=NULL, groups.label=NULL, filter.tags=c('ctl.tag', 'blood.draw.phase.tag'), groups.to.filter=list('CTL', 'convalescent'), keep=TRUE, color.vals=tmp.cols)
tmp.ggplot.2 <- plot.props(seurat.obj=srt.objs.list[[main.obj.1]], plot.type='donut', groups.tag='pr.tag', na.rm=TRUE, groups.of.int=NULL, groups.label=NULL, filter.tags=c('tfh.tag', 'blood.draw.phase.tag'), groups.to.filter=list('TFH', 'convalescent'), keep=TRUE, color.vals=tmp.cols)
# Output.
publish.plot(tmp.ggplot=tmp.ggplot.1, output.path=fig.1.path, file.name='CovidPt_pRFractions-CD4-CTL', blank.comp=blank.complement.4, do.legend=TRUE, height=5, width=5)
publish.plot(tmp.ggplot=tmp.ggplot.2, output.path=fig.1.path, file.name='CovidPt_pRFractions-TFH', blank.comp=blank.complement.4, do.legend=TRUE, height=5, width=5)


### -------------------- Supplementary Figure 1 --------------------- ###
# ----> Output directory.
sup.1.path <- paste0(reports.path, '/sup_figure_1')
if(!dir.exists(sup.1.path)) dir.create(sup.1.path)
# ---> Comments: Main dataset for this section, PT-24-ALL-SF

# ---> SARS-CoV-2-specific T-cell response stability
# Get data.
tmp.data <- gather(data=sr.cell.counts, key='point', value='cells', `cells.acute`, `cells.memory`)
tmp.data$hosp.status <- ifelse(test=tmp.data$hosp.status=='Hosp', yes='Yes', no='No')
tmp.data$point <- str_replace(string=tmp.data$point, pattern='cells\\.', replacement='')
tmp.data$interval <- unlist(lapply(X=1:nrow(tmp.data), FUN=function(idx) return(tmp.data[idx, paste0('interval.', tmp.data[idx, 'point'])])))
tmp.data[, 'interval.memory'] <- NULL; tmp.data[, 'interval.acute'] <- NULL
tmp.data <- tmp.data[tmp.data$time.point==6, ] # Keep data for 6h stimulation only.
tmp.data <- tmp.data[tmp.data$patient.id!='P04', ] # Remove data for patient with NA values from the acute data.
# Plot
for(tmp.status in c('Yes', 'No')){
  tmp.ggplot <- ggplot(data=tmp.data[tmp.data$hosp.status==tmp.status, ], aes(x=interval, y=cells, col=hosp.status)) +
    geom_point(size=5) + geom_line(aes(group=patient.id)) +
    scale_fill_manual(values=hosp.cols[names(hosp.cols)==tmp.status]) + scale_color_manual(values=hosp.cols) +
    scale_y_continuous(limits=c(NA, max(tmp.data$cells)), breaks=c(0, 1000, 2000)) +
    labs(x='Time interval', y='Cell count per 10e6 PBMCs (log10)', col='Status')
  publish.plot(tmp.ggplot=tmp.ggplot, output.path=sup.1.path, file.name=paste0('CovidPt_CellStability_Hosp-', tmp.status), type='pdf', blank.comp=blank.complement.1, do.legend=TRUE)
}

# ---> QC-related figures.
# @ Distribution of the number of genes per cell across sequencing batches.
new.lvls <- as.character(mixedsort(unique(srt.objs.list[[main.obj.1]]@meta.data[, 'seq.batch.tag'])))
srt.objs.list[[main.obj.1]]@meta.data[, 'seq.batch.tag'] <- factor(x=as.character(srt.objs.list[[main.obj.1]]@meta.data[, 'seq.batch.tag']), levels=new.lvls)
tmp.ggplot <- vln.plot(seurat.obj=srt.objs.list[[main.obj.1]], feature='nFeature_RNA', groups.tag='seq.batch.tag', na.rm=TRUE, groups.of.int=NULL, groups.label=NULL, filter.tags=NULL, groups.to.filter=NULL, keep=TRUE, feature.thold=NULL, color='color', size.thold=0, file.name=NULL, adjust.val=1, trim.val=TRUE, plot.limits=NULL, add.points=FALSE, this.color.scale=NULL)
tmp.cols <- rep(x='lightblue', times=meta.data[, uniqueN(seq.batch.tag)]); names(tmp.cols) <- meta.data[, unique(seq.batch.tag)]
tmp.ggplot <- tmp.ggplot + scale_fill_manual(values=tmp.cols) + scale_y_continuous(limits=c(0, NA), breaks=c(1000, 2500, 4000))
publish.plot(tmp.ggplot=tmp.ggplot, output.path=sup.1.path, file.name='CovidPt_QCs_GenesPerCellDist', type='pdf', blank.comp=blank.complement.1, do.legend=TRUE, height=4, width=5)
# @ Fraction of cell cluster per sequencing batch.
# Get data
tmp.data <- get.tag.analysis(seurat.obj=srt.objs.list[[main.obj.1]], tmp.meta.data=as.data.frame(meta.data), tmp.tag='seq.batch.tag', clusters.tag=clusts.lab[main.obj.1], tag.reports.path=NA, reports.pref=NA, vals.cols=NULL, output.umaps=FALSE)
tmp.data <- as.data.frame(t(tmp.data), stringsAsFactors=TRUE); tmp.data$cluster <- row.names(tmp.data)
tmp.data <- tmp.data[tmp.data$cluster %in% names(clusts.cols[[main.obj.1]]), ]
tmp.data <- gather(data=tmp.data, key='seq.batch', value='fraction', -`cluster`)
# Plot and output.
tmp.ggplot <- ggplot(data=tmp.data, aes(x=seq.batch, y=fraction, fill=cluster)) +
  geom_bar(stat='identity', position='fill', alpha=0.7) +
  scale_fill_manual(values=clusts.cols[[main.obj.1]]) +
  scale_y_continuous(breaks=c(0, 0.5, 1)) + coord_flip() +
  labs(x='Sequencing batch', y='Fraction', fill='Cluster')
publish.plot(tmp.ggplot=tmp.ggplot, output.path=sup.1.path, file.name='CovidPt_QCs_ClusterFractsPerSeqBatch', type='pdf', blank.comp=blank.complement.1, do.legend=TRUE)

# ---> Clonotype sharing between virus reactivities.
tmp.cols <- virus.cols[c('SARS-CoV-2', 'FLU')]
depict.clone.sharing(meta.data=sec.meta.data, tag.of.int='virus.tag', filter.tag='blood.draw.phase.tag', groups.to.filter='convalescent', reports.path=sup.1.path, file.preffix='CovidPt-24_TCR-Repertoire_ClonotypeSharingBetweenReactivities', fill.cols=tmp.cols)

# ---> Contrasting clonal expansion between crossreactive and non-crossreactive TCRs seen in the convalescent phase.
# Plot and output.
tmp.cols <- sharing.cols; names(tmp.cols) <- c('pR', 'non-pR')
srt.objs.list[[sec.obj.1]]@meta.data[, 'log.clon.size.tag'] <- log2(srt.objs.list[[sec.obj.1]]@meta.data[, 'clon.size.tag'])
tmp.ggplot <-
vln.plot(seurat.obj=srt.objs.list[[sec.obj.1]], feature='log.clon.size.tag', groups.tag='pr.tag', na.rm=TRUE, groups.of.int=NULL, groups.label=NULL, filter.tags='blood.draw.phase.tag', groups.to.filter=list(c('convalescent')), keep=TRUE, feature.thold=1, color='median', vln.type='density', size.thold=0, file.name=NULL, adjust.val=2, trim.val=TRUE, plot.limits=NULL, add.points=FALSE, this.color.scale=tmp.cols)
tmp.ggplot <- tmp.ggplot + scale_x_continuous(expand=expansion(add=0), breaks=scales::pretty_breaks(n=3)) + scale_y_continuous(expand=expansion(add=0), breaks=scales::pretty_breaks(n=2))#, limits=c(0, 0.2))
publish.plot(tmp.ggplot=tmp.ggplot, output.path=sup.1.path, file.name='CovidPt-24_TCR-Repertoire_SizeBetween-pR-and-nonpR_ExpA', type='pdf', blank.comp=blank.complement.1, do.legend=TRUE, height=4, width=4)

# ---> Cross-reactive TCR representation per donor.
# Get data.
tmp.data <- sec.meta.data[
  !is.na(clonotype.tag) & blood.draw.phase.tag=='convalescent' & virus.tag=='SARS-CoV-2',
  .(
    pr.prop=.SD[pr.tag=='pR', .N]/.N
  ),
  by=.(donor=donor.id.tag, hosp.status=hospitalization.tag)
]
setorderv(x=tmp.data, cols='pr.prop', order=-1)
tmp.data$donor <- factor(x=as.character(tmp.data$donor), levels=as.character(tmp.data$donor))
# Plot and output.
tmp.ggplot <- ggplot(data=tmp.data, aes(x=donor, y=pr.prop)) +
  geom_bar(stat='identity', fill='lightgray', width=0.7) +
  scale_y_continuous(limits=c(0, 1), expand=expansion(add=0), breaks=scales::pretty_breaks(n=3)) + scale_x_discrete(expand=expansion(add=0.5)) +
  labs(x='Donors ranked according to crossreactivity representation', y='Percent of total cells') +
  theme(axis.ticks.x=element_blank())
publish.plot(tmp.ggplot=tmp.ggplot, output.path=sup.1.path, file.name='CovidPt-24_TCR-Repertoire_pR-RepresentationPerDonor_ExpA', type='pdf', blank.comp=blank.complement.1, do.legend=TRUE, width=10)

# ---> FGSEA, cytotoxicity signature.
do.fgsea(metrics.src=srt.objs.list[[main.obj.1]], tag.of.int='pr.tag', metric='signal.to.noise', output.path=sup.1.path, vals.to.depict='pR', files.preffix='CovidPt')

# ---> pR profile for shared TCRs between phases.
# Obtain plots.
tmp.cols <- sharing.cols; names(tmp.cols) <- c('pR', 'non-pR')
tmp.ggplot.1 <- plot.props(seurat.obj=srt.objs.list[[sec.obj.1]], plot.type='pie', groups.tag='pr.tag', na.rm=TRUE, groups.of.int=NULL, groups.label=NULL, filter.tags=c('blood.draw.phase.tag', 'virus.tag'), groups.to.filter=list('convalescent', 'SARS-CoV-2'), keep=TRUE, color.vals=tmp.cols)
# Output.
publish.plot(tmp.ggplot=tmp.ggplot.1, output.path=sup.1.path, file.name='CovidPt-24-TCR-Repertoire_pR-RepresentationOverall', blank.comp=blank.complement.4, do.legend=TRUE, height=5, width=5)
# Report proportions per pR status:
sec.meta.data[blood.draw.phase.tag=='convalescent'&virus.tag=='SARS-CoV-2'&!is.na(clonotype.tag), .(cell.frac=.N/sec.meta.data[blood.draw.phase.tag=='convalescent'&virus.tag=='SARS-CoV-2'&!is.na(clonotype.tag), .N]), by=pr.tag]
# pR 0.2135541
# non-pR 0.7864459
sec.meta.data[blood.draw.phase.tag=='convalescent'&virus.tag=='SARS-CoV-2'&!is.na(clonotype.tag), .N]

# ---> Contrasting clonal expansion between crossreactive and non-crossreactive TCRs seen in the convalescent phase for selected donors.
# Plot and output.
donors.of.int <- c('P07', 'P12')
tmp.cols <- sharing.cols; names(tmp.cols) <- c('pR', 'non-pR')
srt.objs.list[[sec.obj.1]]@meta.data[, 'log.clon.size.tag'] <- log2(srt.objs.list[[sec.obj.1]]@meta.data[, 'clon.size.tag'])
for(tmp.donor in donors.of.int){
  tmp.ggplot <- vln.plot(seurat.obj=srt.objs.list[[sec.obj.1]], feature='log.clon.size.tag', groups.tag='pr.tag', na.rm=TRUE, groups.of.int=NULL, groups.label=NULL, filter.tags=c('blood.draw.phase.tag', 'donor.id.tag'), groups.to.filter=list(c('convalescent'), c(tmp.donor)), keep=TRUE, feature.thold=1, color='median', vln.type='density', size.thold=0, file.name=NULL, adjust.val=2, trim.val=TRUE, plot.limits=NULL, add.points=FALSE, this.color.scale=tmp.cols)
  tmp.ggplot <- if(tmp.donor=='P03') tmp.ggplot + scale_x_continuous(expand=expansion(add=0), breaks=c(5, 8)) + scale_y_continuous(expand=expansion(add=0), breaks=scales::pretty_breaks(n=2)) else tmp.ggplot + scale_x_continuous(expand=expansion(add=0), breaks=c(5, 10)) + scale_y_continuous(expand=expansion(add=0), breaks=scales::pretty_breaks(n=2))
  publish.plot(tmp.ggplot=tmp.ggplot, output.path=sup.1.path, file.name=paste0('CovidPt-24_TCR-Repertoire_SizeBetween-pR-and-nonpR_ExpA_ForDonor-', tmp.donor), type='pdf', blank.comp=blank.complement.1, do.legend=TRUE, height=4, width=4)
}

# ---> Comparison between SARS-reactive cell count and cross-reactive potential -crossreactive cell count- (at the donor level).
# Get data.
tmp.data.1 <- meta.data[
  !is.na(clonotype.tag) & blood.draw.phase.tag=='convalescent' & virus.tag=='SARS-CoV-2',
  .(
    pr.freq=.SD[pr.tag=='pR', .N],
    pr.prop=.SD[pr.tag=='pR', .N]/.N
  ),
  by=.(donor=donor.id.tag, hosp.status=hospitalization.tag)
]
tmp.data.2 <- sr.cell.counts[time.point=='24']
tmp.data.2[, hosp.status:=ifelse(test=hosp.status=='Hosp', yes='Yes', no='No')]
tmp.data <- merge(x=tmp.data.1, y=tmp.data.2, by.x='donor', by.y='patient.id')
# tmp.data[, all(hosp.status.x==hosp.status.y)] # Sanity check. Annotations may maintain stable across files.
# Plot and output.
tmp.ggplot <- ggplot(data=tmp.data, aes(x=cells.memory, y=pr.prop)) +
  geom_point(size=7, shape=1) + geom_smooth(formula=y~x, method='lm', se=FALSE, color='black') +
  scale_y_continuous(expand=expansion(add=0.015), breaks=scales::pretty_breaks(n=3)) + scale_x_continuous(expand=expansion(add=38), breaks=scales::pretty_breaks(n=3)) +
  labs(x='SARS-CoV-2-reactive cell count per 10e6 PBMC', y='Crossreactive cell fraction')
publish.plot(tmp.ggplot=tmp.ggplot, output.path=sup.1.path, file.name='CovidPt-24_TCR-Repertoire_SARSCellCount-vs-pR-RepresentationPerDonor_ExpA', type='pdf', blank.comp=blank.complement.1, do.legend=FALSE, stat.cor=TRUE)

# ---> pR profile on UMAP along clone size depiction
tmp.cells <- sec.meta.data[blood.draw.phase.tag=='convalescent' & !is.na(clonotype.tag), cell]
sec.meta.data[, tmp.clon.size.tag:=clon.size.tag]; sec.meta.data[clon.size.tag>=1500, tmp.clon.size.tag:=1500]
tmp.ggplot <- ggplot(data=sec.meta.data[cell %in% tmp.cells], aes(x=UMAP_1, y=UMAP_2, fill=pr.tag, size=tmp.clon.size.tag)) +
  geom_point(stroke=0.4, shape=21, col='black', alpha=0.8) +
  scale_fill_manual(values=pr.groups.cols) +
  scale_radius(breaks=c(10, 100, 1000), limits=c(1, sec.meta.data[, max(tmp.clon.size.tag, na.rm=TRUE)]), range=c(3, 13)) +
  scale_y_continuous(expand=expansion(add=0)) + scale_x_continuous(expand=expansion(add=0), limits=c(NA, 5)) +
  labs(x='UMAP 1', y='UMAP 2', col='')
# Output.
publish.plot(tmp.ggplot=tmp.ggplot, output.path=sup.1.path, file.name='CovidPt-24_SharingStatusOnUMAP', type='tiff', blank.comp=blank.complement.3, do.legend=TRUE)

# ---> UMAP depicting clusters.
# Plot.
tmp.cells <- sec.meta.data[blood.draw.phase.tag=='convalescent', cell]
tmp.ggplot <- ggplot(data=sec.meta.data[cell %in% tmp.cells], aes_string(x='UMAP_1', y='UMAP_2', col=clusts.lab[sec.obj.1])) +
  geom_point(size=4, alpha=0.4, shape=19) +
  scale_color_manual(values=clusts.cols[[sec.obj.1]]) +
  scale_y_continuous(expand=expansion(add=0)) + scale_x_continuous(expand=expansion(add=0), limits=c(NA, 5)) +
  labs(x='UMAP 1', y='UMAP 2', col='Cluster')
# Output.
publish.plot(tmp.ggplot=tmp.ggplot, output.path=sup.1.path, file.name='CovidPt-24_ClustersOnUMAP', type='tiff', blank.comp=blank.complement.3, do.legend=TRUE)

# ---> Gene signatures.
# Plot.
tmp.cells <- sec.meta.data[blood.draw.phase.tag=='convalescent', cell]
tmp.ggplot <- ggplot(data=sec.meta.data[cell %in% tmp.cells], aes_string(x='UMAP_1', y='UMAP_2', col='cell.cytotoxicity.patil.score', fill='cell.cytotoxicity.patil.score')) +
  geom_point(size=4, alpha=0.8, shape=19) +
  scale_fill_gradientn(colors=signatures.col.scale) +
  scale_color_gradientn(colors=signatures.col.scale) +
  scale_y_continuous(expand=expansion(add=0)) + scale_x_continuous(expand=expansion(add=0), limits=c(NA, 5)) +
  labs(x='UMAP 1', y='UMAP 2', col='Score')
# Output.
publish.plot(tmp.ggplot=tmp.ggplot, output.path=sup.1.path, file.name='CovidPt-24_Signature-Cytotoxicity', type='tiff', blank.comp=blank.complement.3, do.legend=TRUE)

# ---> CTL-related phenotype (fractions) per cross-reactive group (cross-reactive cells and others).
# Define CD4-CTL tag in seurat object.
tmp.corrs <- c(rep(x='CD4-CTL', times=2), rep(x='Rest', times=3)); names(tmp.corrs) <- c('1', '4', '0', 2:3)
tmp.data <- as.character(srt.objs.list[[sec.obj.1]]@meta.data[, clusts.lab[sec.obj.1]])
tmp.data <- tmp.corrs[tmp.data]
srt.objs.list[[sec.obj.1]]@meta.data[, 'ctl.tag'] <- tmp.data
tmp.cols <- c('Rest'='#808080', 'CD4-CTL'='#cc0000')
# Obtain plots.
tmp.ggplot.1 <- plot.props(seurat.obj=srt.objs.list[[sec.obj.1]], plot.type='donut', groups.tag='ctl.tag', na.rm=TRUE, groups.of.int=NULL, groups.label=NULL, filter.tags=c('pr.tag', 'blood.draw.phase.tag'), groups.to.filter=list('pR', 'convalescent'), keep=TRUE, color.vals=tmp.cols)
tmp.ggplot.2 <- plot.props(seurat.obj=srt.objs.list[[sec.obj.1]], plot.type='donut', groups.tag='ctl.tag', na.rm=TRUE, groups.of.int=NULL, groups.label=NULL, filter.tags=c('pr.tag', 'blood.draw.phase.tag'), groups.to.filter=list('non-pR', 'convalescent'), keep=TRUE, color.vals=tmp.cols)
# Output.
publish.plot(tmp.ggplot=tmp.ggplot.1, output.path=sup.1.path, file.name='CovidPt-24_CD4-CTLs-and-pRAssessment-pR', blank.comp=blank.complement.4, do.legend=TRUE, height=5, width=5)
publish.plot(tmp.ggplot=tmp.ggplot.2, output.path=sup.1.path, file.name='CovidPt-24_CD4-CTLs-and-pRAssessment-nonpR', blank.comp=blank.complement.4, do.legend=TRUE, height=5, width=5)

# ---> pR profile (fraction) for CD4-CTLs and TREGs
# Define CD4-CTL tag in seurat object.
tmp.corrs <- c(rep(x='CD4-CTL', times=2), rep(x='Rest', times=3)); names(tmp.corrs) <- c('1', '4', '0', 2:3)
tmp.data <- as.character(srt.objs.list[[sec.obj.1]]@meta.data[, clusts.lab[sec.obj.1]])
tmp.data <- tmp.corrs[tmp.data]
srt.objs.list[[sec.obj.1]]@meta.data[, 'ctl.tag'] <- tmp.data
# Define TREG tag in seurat object.
tmp.corrs <- c(rep(x='TREG', times=2), rep(x='Rest', times=3)); names(tmp.corrs) <- c('0', '3', 1:2, '4')
tmp.data <- as.character(srt.objs.list[[sec.obj.1]]@meta.data[, clusts.lab[sec.obj.1]])
tmp.data <- tmp.corrs[tmp.data]
srt.objs.list[[sec.obj.1]]@meta.data[, 'treg.tag'] <- tmp.data
# seurat.obj@meta.data[, 'ctl.tag'] <- tmp.data
tmp.cols <- sharing.cols; names(tmp.cols) <- c('pR', 'non-pR')
# Obtain plots.
tmp.ggplot.1 <- plot.props(seurat.obj=srt.objs.list[[sec.obj.1]], plot.type='donut', groups.tag='pr.tag', na.rm=TRUE, groups.of.int=NULL, groups.label=NULL, filter.tags=c('ctl.tag', 'blood.draw.phase.tag'), groups.to.filter=list('CD4-CTL', 'convalescent'), keep=TRUE, color.vals=tmp.cols)
tmp.ggplot.2 <- plot.props(seurat.obj=srt.objs.list[[sec.obj.1]], plot.type='donut', groups.tag='pr.tag', na.rm=TRUE, groups.of.int=NULL, groups.label=NULL, filter.tags=c('treg.tag', 'blood.draw.phase.tag'), groups.to.filter=list('TREG', 'convalescent'), keep=TRUE, color.vals=tmp.cols)
# Output.
publish.plot(tmp.ggplot=tmp.ggplot.1, output.path=sup.1.path, file.name='CovidPt-24_pRFractions-CD4-CTL', blank.comp=blank.complement.4, do.legend=TRUE, height=5, width=5)
publish.plot(tmp.ggplot=tmp.ggplot.2, output.path=sup.1.path, file.name='CovidPt-24_pRFractions-TREG', blank.comp=blank.complement.4, do.legend=TRUE, height=5, width=5)

# ---> Dot plot, 6-h data.
# @ Version A.
# Define genes of interest.
these.markers <- c(
  'CD200', 'BTLA', 'POU2AF1', 'IL21', # TFH
  'GZMB', 'PRF1', 'GNLY', 'NKG7', 'XCL1', 'XCL2', # CD4-CTL
  'IFNG', 'IL2', 'TNF', 'CSF2', # TH1/polyfunctional
  'LTB', 'FLT3LG', 'CCR6', 'CTSH', 'IL4I1', # TH17/polyfunctional
  'ISG15', 'MX1', 'OAS1', 'IFIT3', # THIFNR
  'MKI67', 'CDK1' # Proliferative
)
# Define clusters' order.
these.pops <- c('0', '5', '1', '4', '2', '3', '6')
# Plot.
tmp.ggplot <- dot.plot(
  seurat.obj=srt.objs.list[[main.obj.1]], features=these.markers, slot='data', do.norm=FALSE, ensembl=FALSE,
  groups.tag=clusts.lab[main.obj.1], groups.order=these.pops, groups.of.int=NULL, filter.tag=NULL, groups.to.filter=NULL, keep=TRUE, na.rm=TRUE, feature.thold=NULL,
  this.color.scale=signatures.col.scale, col.min=NULL, col.max=NULL,
  scale.by='radius', dot.scale=6, size.min=NA, size.max=NA,
  file.name=NULL
)
# Output.
publish.plot(tmp.ggplot=tmp.ggplot, output.path=sup.1.path, file.name='CovidPt_Markers_Opt-A', type='pdf', blank.comp=blank.complement.2, do.legend=TRUE, height=7, width=5)

# ---> Dot plot, 24-h data.
# @ Version A.
# Define genes of interest.
these.markers <- c(
  'FOXP3', 'TIGIT', 'CCR8', 'IKZF2',
  'GZMB', 'PRF1', 'GNLY', # CD4-CTL
  'CD200', 'BTLA' # TFH
)
# Define clusters' order.
these.pops <- c('0', '3', '1', '4', '2')
# Plot.
tmp.ggplot <- dot.plot(
  seurat.obj=srt.objs.list[[sec.obj.1]], features=these.markers, slot='data', do.norm=FALSE, ensembl=FALSE,
  groups.tag=clusts.lab[sec.obj.1], groups.order=these.pops, groups.of.int=NULL, filter.tag=NULL, groups.to.filter=NULL, keep=TRUE, na.rm=TRUE, feature.thold=NULL,
  this.color.scale=signatures.col.scale, col.min=NULL, col.max=NULL,
  scale.by='radius', dot.scale=6, size.min=NA, size.max=NA,
  file.name=NULL
)
# Output.
publish.plot(tmp.ggplot=tmp.ggplot, output.path=sup.1.path, file.name='CovidPt-24_Markers_Opt-A', type='pdf', blank.comp=blank.complement.2, do.legend=TRUE, height=7, width=5)


### ------------- Other things that may fit in Figure 1 ------------- ###

# ---> Contrasting clonal expansion between crossreactive and non-crossreactive TCRs seen in the convalescent phase per donor.
# Plot and output.
tmp.cols <- sharing.cols; names(tmp.cols) <- c('pR', 'non-pR')
srt.objs.list[[main.obj.1]]@meta.data[, 'log.clon.size.tag'] <- log2(srt.objs.list[[main.obj.1]]@meta.data[, 'clon.size.tag'])
tmp.ggplot <-
vln.plot(seurat.obj=srt.objs.list[[main.obj.1]], feature='log.clon.size.tag', groups.tag='donor.id.tag', na.rm=TRUE, groups.of.int=NULL, groups.label=NULL, filter.tags='blood.draw.phase.tag', groups.to.filter=list(c('convalescent')), keep=TRUE, feature.thold=NULL, color='pr.tag', vln.type='split', size.thold=0, file.name=NULL, adjust.val=2, trim.val=TRUE, plot.limits=NULL, add.points=FALSE, this.color.scale=tmp.cols)
publish.plot(tmp.ggplot=tmp.ggplot, output.path=fig.1.path, file.name='CovidPt_TCR-Repertoire_SizeBetween-pR-and-nonpR_ExpA_PerDonor', type='pdf', blank.comp=blank.complement.1, do.legend=TRUE, height=5, width=10)


############    -----------------------------------------    ############
### --------------------------- Figure 2 ---------------------------- ###
############    -----------------------------------------    ############

# ---> General definitions.
# Seurat object of interest.
main.obj.2 <- 'PT-6-Conv-SFCE'
# Seurat object's meta data
meta.data <- cbind(
  data.table(cell=Cells(srt.objs.list[[main.obj.2]])),
  srt.objs.list[[main.obj.2]]@meta.data,
  srt.objs.list[[main.obj.2]]@reductions$umap@cell.embeddings
)

### ------------------------- Main Figure 2 ------------------------- ###
# ----> Output directory.
fig.2.path <- paste0(reports.path, '/figure_2')
if(!dir.exists(fig.2.path)) dir.create(fig.2.path)
# ---> Comments: Main dataset for this section, PT-6-Conv-SFCE

# ---> Numbers needed for overview figure (Fig. 1A).
meta.data[virus.tag=='CMV', .N]
meta.data[virus.tag=='CMV' & !is.na(clonotype.tag), .N]
meta.data[virus.tag=='EBV', .N]
meta.data[virus.tag=='EBV' & !is.na(clonotype.tag), .N]

# ---> Clonotype sharing between viruses.
# @ Dataset as a whole.
tmp.cols <- virus.cols[c('SARS-CoV-2', 'FLU', 'CMV', 'EBV')]
depict.clone.sharing(meta.data=meta.data, tag.of.int='virus.tag', reports.path=fig.2.path, file.preffix='CovidPt_TCR-Repertoire_ClonotypeSharingBetweenReactivities', fill.cols=tmp.cols)

# ---> Contribution to the SARS cross-reactive T cell pool per virus (other than SARS).
# @ Define clonotypes that are reactive to each virus.
tmp.meta <- meta.data[!is.na(clonotype.tag)]
tmp.meta[, tag.of.int:=virus.tag]
tmp.vals <- tmp.meta[, as.character(unique(tag.of.int))]
tmp.data <- lapply(X=tmp.vals, FUN=function(tmp.val){
  to.return <- tmp.meta[tag.of.int==tmp.val, .(count=.N), by=.(clonotype=clonotype.tag)][count>1, clonotype] # This may be changed to adjust for the minimum size for a clonotype to be taken into account.
  return(to.return)
})
names(tmp.data) <- tmp.vals
# @ Define clonotypes that are cross-reactive with SARS per virus reactivity.
sars.cr.clones <- lapply(X=tmp.vals[tmp.vals!='SARS-CoV-2'], FUN=function(tmp.virus){
  to.return <- intersect(x=tmp.data[[tmp.virus]], y=tmp.data[['SARS-CoV-2']])
  return(to.return)
})
names(sars.cr.clones) <- tmp.vals[tmp.vals!='SARS-CoV-2']
# @ Define groups of TCRs cross-reactive to SARS (i.e., combinations of different viruses).
cr.groups <- combn(x=tmp.vals[tmp.vals!='SARS-CoV-2'], m=2, simplify=FALSE)
names(cr.groups) <- unlist(lapply(X=cr.groups, FUN=paste0, collapse='|'))
cr.for.all <- Reduce(x=sars.cr.clones, f=intersect, accumulate=FALSE)
cr.groups <- lapply(X=cr.groups, FUN=function(cr.group){
  to.assess <- sars.cr.clones[cr.group]
  to.return <- Reduce(x=to.assess, f=intersect, accumulate=FALSE)
  to.return <- setdiff(x=to.return, y=cr.for.all)
})
cr.groups[['FLU|CMV|EBV']] <- cr.for.all
cr.for.many <- unique(unlist(cr.groups))
sars.cr.clones <- lapply(X=sars.cr.clones, FUN=setdiff, y=cr.for.many)
sars.cr.clones <- c(sars.cr.clones, cr.groups)
# @ Restructure data as a data table.
sars.cr.clones <- lapply(X=names(sars.cr.clones), FUN=function(cr.group){
  if(length(sars.cr.clones[[cr.group]])==0) return(NA)
  to.return <- data.table(cr.group=cr.group, clonotype=sars.cr.clones[[cr.group]])
  return(to.return)
}); sars.cr.clones <- sars.cr.clones[!is.na(sars.cr.clones)]
sars.cr.clones <- rbindlist(l=sars.cr.clones, use.names=TRUE)
# @ Annotate clonotypes with single cross-reactivity with SARS-CoV-2 in the seurat obj.
tmp.groups <- c('FLU', 'CMV', 'FLU|CMV')
tmp.data.1 <- sars.cr.clones#[cr.group %in% tmp.groups]
tmp.data.2 <- data.table(cell=Cells(srt.objs.list[[main.obj.2]]), clonotype=srt.objs.list[[main.obj.2]]@meta.data[, 'clonotype.tag'])
tmp.data <- merge(x=tmp.data.2, y=tmp.data.1, by='clonotype', sort=FALSE, all=TRUE)
# all(tmp.data[, cell]==Cells(srt.objs.list[[main.obj.2]]))
srt.objs.list[[main.obj.2]]@meta.data[, 'sars.cr.group'] <- tmp.data[, cr.group]
tmp.data[!is.na(cr.group), .N, by=cr.group]
# @ Tag-specific analysis.
# Get data
tmp.data <- get.tag.analysis(seurat.obj=srt.objs.list[[main.obj.2]], tmp.meta.data=as.data.frame(meta.data), tmp.tag='sars.cr.group', clusters.tag=clusts.lab[main.obj.2], tag.reports.path=fig.2.path, reports.pref='CovidPt_TagAnalysis_CRGroup', vals.cols=virus.cols[c('FLU', 'CMV', 'EBV', 'FLU|CMV', 'FLU|CMV|EBV' )], output.umaps=FALSE, output.single=TRUE, blank.complement=blank.complement.3, alpha.val=1)
# @ Summary.
tmp.data <- as.data.table(srt.objs.list[[main.obj.2]]@meta.data)
tmp.data[
  !is.na(sars.cr.group),
  .(
    total.cells=.N,
    total.clones=uniqueN(clonotype.tag),
    size.median=median(clon.size.tag)
  ),
  by=sars.cr.group
]
tmp.data.2 <- tmp.data[
  !is.na(clonotype.tag) & virus.tag=='CMV',
  .N,
  by=sars.cr.group
]
tmp.data.2[, sars.reactive:=ifelse(test=!is.na(sars.cr.group), yes='Yes', no='No')]
tmp.data.2[sars.reactive=='Yes', sum(N)]/tmp.data.2[, sum(N)]

# ---> Sharing and expansion relation.
# Get data.
tmp.data <- meta.data[!is.na(clonotype.tag), .(size=.N), by=.(clone.id=clonotype.tag, virus=virus.tag)]
tmp.data <- spread(data=tmp.data, key=virus, value=size, fill=0)
# Save raw data.
tmp.file.name <- paste0(fig.2.path, '/CovidPt_TCR-RepertoireAcrossReactivities.R.pdf')
fwrite(file=tmp.file.name, x=tmp.data, na=NA, quote=FALSE)
# Further process data to calculate z-scores.
tmp.data[, clone.id:=NULL]
tmp.data <- as.matrix(tmp.data)
rownames(tmp.data) <- paste0('row.', 1:nrow(tmp.data))

# @ ---> Option 1.
# For order.
to.order <- tmp.data
to.order[to.order>=1] <- 1
to.order <- to.order[rowSums(to.order)>1, ]
cols.order <- c('SARS-CoV-2', 'FLU', 'CMV', 'EBV')
to.order <- to.order[, cols.order]
for(tmp.col in rev(cols.order)){
  to.order <- to.order[order(to.order[, tmp.col], decreasing=TRUE), ]
}
# Calculate z-scores.
to.display <- tmp.data[row.names(to.order), ]
to.display[to.display==0] <- NA
to.display <- scale(x=to.display, center=TRUE, scale=TRUE)
# Set a maximum z-score value.
to.display[to.display>1.5] <- 1.5
# Set metadata.
col.metadata <- data.frame(
  row.names=colnames(tmp.data),
  virus.reactivity=colnames(tmp.data)
)
# Define colors for metadata tracks.
ann.colors <- list(
  virus.reactivity=virus.cols[colnames(tmp.data)]
)
# Set color scale and breaks for heatmap.
col.breaks <- seq(from=min(range(to.display, na.rm=TRUE)), to=max(range(to.display, na.rm=TRUE)), length.out=100)
mid.point <- which.min(abs(col.breaks - 0))
hmap.col.scale.1 <- colorRampPalette(c('blue', 'mediumblue', 'black'))(mid.point)
hmap.col.scale.2 <- colorRampPalette(c('black', 'gold', 'yellow'))(100-(mid.point+1))
hmap.col.scale <- c(hmap.col.scale.1, hmap.col.scale.2)
# hmap.col.scale <- colorRampPalette(c('blue', 'yellow'))(101)
# Plot 'complete' version
tmp.file.name <- paste0(fig.2.path, '/CovidPt_TCR-RepertoireAcrossReactivities.C.pdf')
pheatmap(mat=to.display, color=hmap.col.scale, breaks=col.breaks, scale='none', cluster_rows=FALSE, cluster_cols=FALSE, annotation_col=col.metadata, annotation_colors=ann.colors, show_rownames=FALSE, show_colnames=FALSE, filename=tmp.file.name, na_col='gray')
# Plot 'blank' version
tmp.file.name <- paste0(fig.2.path, '/CovidPt_TCR-RepertoireAcrossReactivities.B.pdf')
pheatmap(mat=to.display, color=hmap.col.scale, breaks=col.breaks, scale='none', cluster_rows=FALSE, cluster_cols=FALSE, annotation_col=col.metadata, annotation_colors=ann.colors, legend=FALSE, annotation_legend=FALSE, annotation_names_row=FALSE, annotation_names_col=FALSE, show_rownames=FALSE, show_colnames=FALSE, filename=tmp.file.name)

# @ ---> Option 2.
# For order.
to.order <- tmp.data
to.order[to.order>=1] <- 1
to.order <- to.order[rowSums(to.order)>1, ]
cols.order <- c('SARS-CoV-2', 'FLU', 'CMV', 'EBV')
to.order <- to.order[, cols.order]
for(tmp.col in rev(cols.order)){
  to.order <- to.order[order(to.order[, tmp.col], decreasing=TRUE), ]
}
# Set categorical values.
to.display <- tmp.data[row.names(to.order), ]
to.display[to.display>1] <- 2
to.display[to.display==0] <- NA
# Set metadata.
col.metadata <- data.frame(
  row.names=colnames(tmp.data),
  virus.reactivity=colnames(tmp.data)
)
# Define colors for metadata tracks.
ann.colors <- list(
  virus.reactivity=virus.cols[colnames(tmp.data)]
)
# Set color scale and breaks for heatmap.
col.breaks <- seq(from=min(range(to.display, na.rm=T)), to=max(range(to.display, na.rm=T)), length.out=3)
hmap.col.scale <- c('blue', 'yellow')
# Plot 'complete' version
tmp.file.name <- paste0(fig.2.path, '/CovidPt_TCR-RepertoireAcrossReactivities_B.C.pdf')
pheatmap(mat=to.display, color=hmap.col.scale, breaks=col.breaks, scale='none', cluster_rows=FALSE, cluster_cols=FALSE, annotation_col=col.metadata, annotation_colors=ann.colors, show_rownames=FALSE, show_colnames=FALSE, filename=tmp.file.name, na_col='gray')
# Plot 'blank' version
tmp.file.name <- paste0(fig.2.path, '/CovidPt_TCR-RepertoireAcrossReactivities_B.B.pdf')
pheatmap(mat=to.display, color=hmap.col.scale, breaks=col.breaks, scale='none', cluster_rows=FALSE, cluster_cols=FALSE, annotation_col=col.metadata, annotation_colors=ann.colors, legend=FALSE, annotation_legend=FALSE, annotation_names_row=FALSE, annotation_names_col=FALSE, show_rownames=FALSE, show_colnames=FALSE, filename=tmp.file.name)

# ---> pR profile on UMAP along clone size depiction
tmp.cells <- meta.data[!is.na(clonotype.tag), cell]
meta.data[, tmp.clon.size.tag:=clon.size.tag]; meta.data[clon.size.tag>=1500, tmp.clon.size.tag:=1500]
tmp.ggplot <- ggplot(data=meta.data[cell %in% tmp.cells], aes(x=UMAP_1, y=UMAP_2, fill=pr.tag, size=tmp.clon.size.tag)) +
  geom_point(stroke=0.4, shape=21, col='black') +
  scale_fill_manual(values=pr.groups.cols) +
  scale_radius(breaks=c(3, 100, 500, 1000), limits=c(1, meta.data[, max(tmp.clon.size.tag, na.rm=TRUE)]), range=c(3, 13)) +
  scale_y_continuous(expand=expansion(add=0)) + scale_x_continuous(expand=expansion(add=0)) +
  labs(x='UMAP 1', y='UMAP 2', col='')
# Output.
publish.plot(tmp.ggplot=tmp.ggplot, output.path=fig.2.path, file.name='CovidPt_SharingStatusOnUMAP', type='tiff', blank.comp=blank.complement.3, do.legend=TRUE)

# ---> Contrasting clonal expansion between cross-reactive and others.
# Plot and output.
tmp.cols <- sharing.cols; names(tmp.cols) <- c('pR', 'non-pR')
srt.objs.list[[main.obj.2]]@meta.data[, 'log.clon.size.tag'] <- log2(srt.objs.list[[main.obj.2]]@meta.data[, 'clon.size.tag'])
tmp.ggplot <-
vln.plot(seurat.obj=srt.objs.list[[main.obj.2]], feature='log.clon.size.tag', groups.tag='pr.tag', na.rm=TRUE, groups.of.int=NULL, groups.label=NULL, filter.tags=NULL, groups.to.filter=NULL, keep=TRUE, feature.thold=1, color='median', vln.type='density', size.thold=0, file.name=NULL, adjust.val=2, trim.val=TRUE, plot.limits=NULL, add.points=FALSE, this.color.scale=tmp.cols)
tmp.ggplot <- tmp.ggplot + scale_x_continuous(expand=expansion(add=0), breaks=c(5, 10)) + scale_y_continuous(expand=expansion(add=0), breaks=scales::pretty_breaks(n=2))
publish.plot(tmp.ggplot=tmp.ggplot, output.path=fig.2.path, file.name='CovidPt_TCR-Repertoire_SizeBetween-pR-and-nonpR', type='pdf', blank.comp=blank.complement.1, do.legend=TRUE)

# ---> Gene signatures.
# Plot.
tmp.cells <- meta.data[, cell]
tmp.ggplot <- ggplot(data=meta.data[cell %in% tmp.cells], aes_string(x='UMAP_1', y='UMAP_2', col='cell.cytotoxicity.patil.score', fill='cell.cytotoxicity.patil.score')) +
  geom_point(size=4, alpha=0.8, shape=19) +
  scale_fill_gradientn(colors=signatures.col.scale) +
  scale_color_gradientn(colors=signatures.col.scale) +
  scale_y_continuous(expand=expansion(add=0)) + scale_x_continuous(expand=expansion(add=0)) +
  labs(x='UMAP 1', y='UMAP 2', col='Score')
# Output.
publish.plot(tmp.ggplot=tmp.ggplot, output.path=fig.2.path, file.name='CovidPt_Signature-Cytotoxicity', type='tiff', blank.comp=blank.complement.3, do.legend=TRUE)

# ---> FGSEA, cytotoxicity signature.
do.fgsea(metrics.src=srt.objs.list[[main.obj.2]], tag.of.int='pr.tag', metric='signal.to.noise', output.path=fig.2.path, vals.to.depict='pR', files.preffix='CovidPt')

# ---> Cell fractions per cross-reactive status for CT4-CTLs (cluster 0) and non-CD4-CTLs (rest of the clusters).
# Cell fractions per cross-reactive status for CD4-CTLs.
tmp.ggplot <- plot.props(seurat.obj=srt.objs.list[[main.obj.2]], plot.type='donut', groups.tag='pr.tag', na.rm=TRUE, groups.of.int=NULL, groups.label=NULL, filter.tags=clusts.lab[main.obj.2], groups.to.filter='0', keep=TRUE, color.vals=tmp.cols, file.name=NULL)
publish.plot(tmp.ggplot=tmp.ggplot, output.path=fig.2.path, file.name='CovidPt_CellFracPerCrossReactiveStatus_CTLs', type='pdf', blank.comp=blank.complement.4, do.legend=TRUE)
# Cell fractions per cross-reactive status for non-CD4-CTLs.
tmp.ggplot <- plot.props(seurat.obj=srt.objs.list[[main.obj.2]], plot.type='donut', groups.tag='pr.tag', na.rm=TRUE, groups.of.int=NULL, groups.label=NULL, filter.tags=clusts.lab[main.obj.2], groups.to.filter='0', keep=FALSE, color.vals=tmp.cols, file.name=NULL)
publish.plot(tmp.ggplot=tmp.ggplot, output.path=fig.2.path, file.name='CovidPt_CellFracPerCrossReactiveStatus_Others', type='pdf', blank.comp=blank.complement.4, do.legend=TRUE)

# ---> pR profile across CD4-CTLs and rest of the cells.
# Define CD4-CTL tag in seurat object.
tmp.corrs <- c('CD4-CTL', rep(x='Rest', times=2)); names(tmp.corrs) <- as.character(0:2)
tmp.data <- as.character(srt.objs.list[[main.obj.2]]@meta.data[, clusts.lab[main.obj.2]])
tmp.data <- tmp.corrs[tmp.data]
srt.objs.list[[main.obj.2]]@meta.data[, 'ctl.tag'] <- tmp.data
# seurat.obj@meta.data[, 'ctl.tag'] <- tmp.data
tmp.cols <- c('Rest'='#808080', 'CD4-CTL'='#cc0000')
# Obtain plots.
tmp.ggplot.1 <- plot.props(seurat.obj=srt.objs.list[[main.obj.2]], plot.type='donut', groups.tag='ctl.tag', na.rm=TRUE, groups.of.int=NULL, groups.label=NULL, filter.tags='pr.tag', groups.to.filter='pR', keep=TRUE, color.vals=tmp.cols)
tmp.ggplot.2 <- plot.props(seurat.obj=srt.objs.list[[main.obj.2]], plot.type='donut', groups.tag='ctl.tag', na.rm=TRUE, groups.of.int=NULL, groups.label=NULL, filter.tags='pr.tag', groups.to.filter='non-pR', keep=TRUE, color.vals=tmp.cols)
# Output.
publish.plot(tmp.ggplot=tmp.ggplot.1, output.path=fig.2.path, file.name='CovidPt_CD4-CTLs-and-pRAssessment-pR', blank.comp=blank.complement.4, do.legend=TRUE, height=5, width=5)
publish.plot(tmp.ggplot=tmp.ggplot.2, output.path=fig.2.path, file.name='CovidPt_CD4-CTLs-and-pRAssessment-nonpR', blank.comp=blank.complement.4, do.legend=TRUE, height=5, width=5)


### -------------------- Supplementary Figure 2 --------------------- ###
# ----> Output directory.
sup.2.path <- paste0(reports.path, '/sup_figure_2')
if(!dir.exists(sup.2.path)) dir.create(sup.2.path)
# ---> Comments: Main dataset for this section, PT-6-Conv-SFCE

# ---> QC-related figures.
# @ Distribution of the number of genes per cell across sequencing batches.
new.lvls <- as.character(mixedsort(unique(srt.objs.list[[main.obj.2]]@meta.data[, 'seq.batch.tag'])))
srt.objs.list[[main.obj.2]]@meta.data[, 'seq.batch.tag'] <- factor(x=as.character(srt.objs.list[[main.obj.2]]@meta.data[, 'seq.batch.tag']), levels=new.lvls)
tmp.ggplot <- vln.plot(seurat.obj=srt.objs.list[[main.obj.2]], feature='nFeature_RNA', groups.tag='seq.batch.tag', na.rm=TRUE, groups.of.int=NULL, groups.label=NULL, filter.tags=NULL, groups.to.filter=NULL, keep=TRUE, feature.thold=NULL, color='color', size.thold=0, file.name=NULL, adjust.val=1, trim.val=TRUE, plot.limits=NULL, add.points=FALSE, this.color.scale=NULL)
tmp.cols <- rep(x='lightblue', times=meta.data[, uniqueN(seq.batch.tag)]); names(tmp.cols) <- meta.data[, unique(seq.batch.tag)]
tmp.ggplot <- tmp.ggplot + scale_fill_manual(values=tmp.cols) + scale_y_continuous(breaks=c(1000, 2500, 4000))
publish.plot(tmp.ggplot=tmp.ggplot, output.path=sup.2.path, file.name='CovidPt_QCs_GenesPerCellDist', type='pdf', blank.comp=blank.complement.1, do.legend=TRUE, height=4, width=5)
# @ Fraction of cell cluster per sequencing batch.
# Get data
tmp.data <- get.tag.analysis(seurat.obj=srt.objs.list[[main.obj.2]], tmp.meta.data=as.data.frame(meta.data), tmp.tag='seq.batch.tag', clusters.tag=clusts.lab[main.obj.2], tag.reports.path=NA, reports.pref=NA, vals.cols=NULL, output.umaps=FALSE)
tmp.data <- as.data.frame(t(tmp.data), stringsAsFactors=TRUE); tmp.data$cluster <- row.names(tmp.data)
tmp.data <- tmp.data[tmp.data$cluster %in% names(clusts.cols[[main.obj.1]]), ]
tmp.data <- gather(data=tmp.data, key='seq.batch', value='fraction', -`cluster`)
# Plot and output.
tmp.ggplot <- ggplot(data=tmp.data, aes(x=seq.batch, y=fraction, fill=cluster)) +
  geom_bar(stat='identity', position='fill', alpha=0.7) +
  scale_fill_manual(values=clusts.cols[[main.obj.2]]) +
  scale_y_continuous(breaks=c(0, 0.5, 1)) + coord_flip() +
  labs(x='Sequencing batch', y='Fraction', fill='Cluster')
publish.plot(tmp.ggplot=tmp.ggplot, output.path=sup.2.path, file.name='CovidPt_QCs_ClusterFractsPerSeqBatch', type='pdf', blank.comp=blank.complement.1, do.legend=TRUE)

# ---> UMAP depicting clusters.
# Plot.
tmp.cells <- meta.data[, cell]
tmp.ggplot <- ggplot(data=meta.data[cell %in% tmp.cells], aes_string(x='UMAP_1', y='UMAP_2', col=clusts.lab[main.obj.2])) +
  geom_point(size=4, alpha=0.4, shape=19) +
  scale_color_manual(values=clusts.cols[[main.obj.2]]) +
  scale_y_continuous(expand=expansion(add=0)) + scale_x_continuous(expand=expansion(add=0)) +
  labs(x='UMAP 1', y='UMAP 2', col='Cluster')
# Output.
publish.plot(tmp.ggplot=tmp.ggplot, output.path=sup.2.path, file.name='CovidPt_ClustersOnUMAP', type='tiff', blank.comp=blank.complement.3, do.legend=TRUE)

# ---> Dot plot.
# @ Version A.
# Define genes of interest.
these.markers <- c(
  'GZMB', 'PRF1', 'GNLY', 'NKG7', 'XCL1', 'XCL2' # CD4-CTL
)
# Define clusters' order.
these.pops <- c('0', '1', '2')
# Plot.
tmp.ggplot <- dot.plot(
  seurat.obj=srt.objs.list[[main.obj.2]], features=these.markers, slot='data', do.norm=FALSE, ensembl=FALSE,
  groups.tag=clusts.lab[main.obj.2], groups.order=NULL, groups.of.int='0', filter.tag=NULL, groups.to.filter=NULL, keep=TRUE, na.rm=TRUE, feature.thold=NULL,
  this.color.scale=signatures.col.scale, col.min=NULL, col.max=NULL,
  scale.by='radius', dot.scale=6, size.min=NA, size.max=NA,
  file.name=NULL
)
# Output.
publish.plot(tmp.ggplot=tmp.ggplot, output.path=sup.2.path, file.name='CovidPt_Markers_Opt-A', type='pdf', blank.comp=blank.complement.2, do.legend=TRUE, height=4, width=3)

# ---> Tag-specific analysis for virus reactivity.
# Get data
tmp.data <- get.tag.analysis(seurat.obj=srt.objs.list[[main.obj.2]], tmp.meta.data=as.data.frame(meta.data), tmp.tag='virus.tag', clusters.tag=clusts.lab[main.obj.2], tag.reports.path=sup.2.path, reports.pref='CovidPt_TagAnalysis_Virus', vals.cols=virus.cols[c('SARS-CoV-2', 'FLU', 'CMV', 'EBV')], output.umaps=TRUE, downsample=FALSE, blank.complement=blank.complement.3)
tmp.data <- as.data.frame(t(tmp.data), stringsAsFactors=TRUE); tmp.data$cluster <- row.names(tmp.data)
tmp.data <- tmp.data[tmp.data$cluster %in% names(clusts.cols[[main.obj.2]]), ]
tmp.data <- gather(data=tmp.data, key='virus', value='fraction', -`cluster`)
# Plot and output.
tmp.ggplot <- ggplot(data=tmp.data, aes(x=virus, y=fraction, fill=cluster)) +
  geom_bar(stat='identity', position='fill', alpha=0.7) +
  scale_fill_manual(values=clusts.cols[[main.obj.2]]) +
  scale_y_continuous(breaks=c(0, 0.5, 1)) + coord_flip() +
  labs(x='Sequencing batch', y='Fraction', fill='Cluster')
publish.plot(tmp.ggplot=tmp.ggplot, output.path=sup.2.path, file.name='CovidPt_TagAnalysis_GenProps', type='pdf', blank.comp=blank.complement.1, do.legend=TRUE)


### ------------- Other things that may fit in Figure 2 ------------- ###

# ---> Cytotoxicity signature comparison for pR status
# Plot and output.
# srt.objs.list[[main.obj.1]]@meta.data[, 'log.clon.size.tag'] <- log2(srt.objs.list[[main.obj.1]]@meta.data[, 'clon.size.tag'])
tmp.ggplot <- vln.plot(seurat.obj=srt.objs.list[[main.obj.2]], feature='cell.cytotoxicity.patil.score', groups.tag='pr.tag', na.rm=TRUE, groups.of.int=NULL, groups.label=NULL, filter.tags=NULL, groups.to.filter=NULL, keep=TRUE, feature.thold=NULL, color='median', size.thold=0, file.name=NULL, adjust.val=1, trim.val=TRUE, plot.limits=NULL, add.points=FALSE, this.color.scale=NULL)
tmp.ggplot <- tmp.ggplot + scale_y_continuous(expand=expansion(add=0))
publish.plot(tmp.ggplot=tmp.ggplot, output.path=fig.2.path, file.name='CovidPt_TCR-Repertoire_CytotoxicityScoreBetween-pR-and-nonpR', type='pdf', blank.comp=blank.complement.1, do.legend=TRUE, width=5)

# ---> Clonal expansion per sharing status (according to SARS) per antigen.
# Add tag of clonotype sharing according to SARS.
sars.clones <- meta.data[
  !is.na(clonotype.tag) & virus.tag=='SARS-CoV-2',
  .(count=.N),
  by=clonotype.tag
][count>1, clonotype.tag]
viruses.of.int <- c('FLU', 'CMV', 'EBV')
cr.status.info <- lapply(X=viruses.of.int, FUN=function(tmp.virus){
  to.return <- meta.data[
    !is.na(clonotype.tag) & virus.tag==tmp.virus,
    .(
      cr.status.tag=clonotype.tag%chin%sars.clones,
      ag.clon.size.tag=.N
    ),
    by=clonotype.tag
  ]
  return(to.return)
})
names(cr.status.info) <- viruses.of.int
cr.status.info <- rbindlist(l=cr.status.info, use.names=TRUE, idcol='virus.tag')
meta.data <- merge(x=meta.data, y=cr.status.info, by=c('virus.tag', 'clonotype.tag'), sort=FALSE, all.x=TRUE, all.y=FALSE)
# Merge cross-reactivity status info with the ag info and further include this info in the seurat object as well.
meta.data[
  !is.na(cr.status.tag),
  new.cr.status.tag:=paste0(
      virus.tag,
      ifelse(test=cr.status.tag==TRUE, yes='+', no='-')
  )
]
meta.data[, cr.status.tag:=new.cr.status.tag]; meta.data[, new.cr.status.tag:=NULL]
# all(tmp[, cell]==Cells(srt.objs.list[[main.obj.2]]))
new.lvls <- mixedsort(unique(meta.data[!is.na(cr.status.tag), cr.status.tag]), decreasing=TRUE)
srt.objs.list[[main.obj.2]]@meta.data[, 'cr.status.tag'] <- factor(x=meta.data[, cr.status.tag], levels=new.lvls)
srt.objs.list[[main.obj.2]]@meta.data[, 'ag.clon.size.tag'] <- meta.data[, ag.clon.size.tag]
srt.objs.list[[main.obj.2]]@meta.data[, 'log.ag.clon.size.tag'] <- log2(srt.objs.list[[main.obj.2]]@meta.data[, 'ag.clon.size.tag'])
# Plot and output.
tmp.ggplot <- vln.plot(seurat.obj=srt.objs.list[[main.obj.2]], feature='log.ag.clon.size.tag', groups.tag='cr.status.tag', na.rm=TRUE, groups.of.int=NULL, groups.label=NULL, filter.tags=NULL, groups.to.filter=NULL, keep=TRUE, feature.thold=NULL, color='median', size.thold=0, file.name=NULL, adjust.val=4, trim.val=FALSE, plot.limits=NULL, add.points=FALSE, this.color.scale=NULL)
tmp.breaks <- c(1, 2, 50, 100, 500, 1000); tmp.breaks <- log2(tmp.breaks)
tmp.ggplot <- tmp.ggplot + scale_y_continuous(expand=expansion(add=0), limits=c(0, 11), breaks=tmp.breaks)
publish.plot(tmp.ggplot=tmp.ggplot, output.path=fig.2.path, file.name='CovidPt_TCR-Repertoire_CytotoxicityScoreBetween-pR-and-nonpR_AcrossReactivities', type='pdf', blank.comp=blank.complement.1, do.legend=TRUE)


############    -----------------------------------------    ############
### --------------------------- Figure 3 ---------------------------- ###
############    -----------------------------------------    ############

# ---> General definitions.
# Seurat object of interest.
main.obj.3 <- 'PT-6-ALL-SF'
# Seurat object's meta data
meta.data <- cbind(
  data.table(cell=Cells(srt.objs.list[[main.obj.3]])),
  srt.objs.list[[main.obj.3]]@meta.data,
  srt.objs.list[[main.obj.3]]@reductions$umap@cell.embeddings
)

### ------------------------- Main Figure 3 ------------------------- ###
# ----> Output directory.
fig.3.path <- paste0(reports.path, '/figure_3')
if(!dir.exists(fig.3.path)) dir.create(fig.3.path)
# ---> Comments: Main dataset for this section, PT-6-Conv-SFCE

# ---> HLA-DRB3*02:02 profile according to cross-reactivity status.
# Split donors into two groups according to their proportions of cross-reactive CD4 T cells (more or less than 10 percent).
tmp.data.1 <- meta.data[
  !is.na(clonotype.tag) & blood.draw.phase.tag=='convalescent' & virus.tag=='SARS-CoV-2',
  .(
    pr.prop=.SD[pr.tag=='pR', .N]/.N
  ),
  by=.(donor=donor.id.tag, hosp.status=hospitalization.tag)
]
tmp.data.1$donor <- tmp.data.1$donor
tmp.data.1[, pr.group:=ifelse(test=pr.prop<0.05, yes='Poor', no='Rich')]
# Get relevant HLA info.
tmp.data.2 <- hla.data[, .(donor=donor.id.tag, hla.group=`HLA-DRB3`)]
tmp.data <- tmp.data.2[, str_split(string=hla.group, pattern=' \\+ ', simplify=TRUE)]
tmp.data.2$drb3.1 <- tmp.data[, 1]; tmp.data.2$drb3.2 <- tmp.data[, 2]
tmp.data.2[, drb3.status:=ifelse(
  test=str_detect(string=drb3.1, pattern='^02:02') | str_detect(string=drb3.2, pattern='^02:02'),
  yes='DRB3*02:02', no='Other'
)]
# Merge all data of interest and get plot.
tmp.data <- merge(x=tmp.data.1, y=tmp.data.2, by='donor', all=FALSE)
# tmp.data[, .N, by=.(pr.group, drb3.status)]
tmp.cols <- c(`DRB3*02:02`='#cc0000', `Other`='#9a9a9a')
tmp.ggplot <- ggplot(data=tmp.data, aes(x=pr.group, fill=drb3.status)) +
  geom_bar(position='fill', width=0.4) +
  scale_fill_manual(values=tmp.cols) +
  scale_y_continuous(limits=c(0, 1), expand=expansion(add=0), breaks=scales::pretty_breaks(n=2)) + scale_x_discrete(expand=expansion(add=0.5)) +
  labs(x='Group per cross-reactive CD4 T cell percent', y='Donors with DRB3*02:02')
publish.plot(tmp.ggplot=tmp.ggplot, output.path=fig.3.path, file.name='CovidPt-Convalescent_DRB3-02-02-StatusByCrossReactiveGroup.pdf', type='pdf', blank.comp=blank.complement.1, do.legend=TRUE, width=5, height=7)


############    -----------------------------------------    ############
### --------------------------- Figure 4 ---------------------------- ###
############    -----------------------------------------    ############

# ---> General definitions.
# Seurat object of interest.
main.obj.4 <- 'PT-6-ALL-SF'
sec.obj.4 <- 'PT-24-ALL-SF'
# Seurat object's meta data
meta.data <- cbind(
  data.table(cell=Cells(srt.objs.list[[main.obj.4]])),
  srt.objs.list[[main.obj.4]]@meta.data,
  srt.objs.list[[main.obj.4]]@reductions$umap@cell.embeddings
)
sec.meta.data <- cbind(
  data.table(cell=Cells(srt.objs.list[[sec.obj.4]])),
  srt.objs.list[[sec.obj.4]]@meta.data,
  srt.objs.list[[sec.obj.4]]@reductions$umap@cell.embeddings
)

### ------------------------- Main Figure 4 ------------------------- ###
# ----> Output directory.
fig.4.path <- paste0(reports.path, '/figure_4')
if(!dir.exists(fig.4.path)) dir.create(fig.4.path)
# ---> Comments: Main dataset for this section, PT-6-ALL-SF

# ---> Cross-reactive TCR representation per donor.
# Get data.
tmp.data <- meta.data[
  !is.na(clonotype.tag) & blood.draw.phase.tag=='acute' & virus.tag=='SARS-CoV-2',
  .(
    pr.prop=.SD[pr.tag=='pR', .N]/.N
  ),
  by=.(donor=donor.id.tag, hosp.status=hospitalization.tag)
]
setorderv(x=tmp.data, cols='pr.prop', order=-1)
tmp.data$donor <- factor(x=as.character(tmp.data$donor), levels=as.character(tmp.data$donor))
# Plot and output.
tmp.ggplot <- ggplot(data=tmp.data, aes(x=donor, y=pr.prop)) +
  geom_bar(stat='identity', fill='lightgray', width=0.7) +
  scale_y_continuous(limits=c(0, 1), expand=expansion(add=0), breaks=scales::pretty_breaks(n=2)) + scale_x_discrete(expand=expansion(add=0.5)) +
  labs(x='Donors ranked according to crossreactivity representation', y='Percent of total cells') +
  theme(axis.ticks.x=element_blank())
publish.plot(tmp.ggplot=tmp.ggplot, output.path=fig.4.path, file.name='CovidPt-Acute_TCR-Repertoire_pR-RepresentationPerDonor', type='pdf', blank.comp=blank.complement.1, do.legend=TRUE, width=10)

# ---> pR profile for shared TCRs between phases.
# Obtain plots.
tmp.cols <- sharing.cols; names(tmp.cols) <- c('pR', 'non-pR')
tmp.ggplot.1 <- plot.props(seurat.obj=srt.objs.list[[main.obj.4]], plot.type='pie', groups.tag='pr.tag', na.rm=TRUE, groups.of.int=NULL, groups.label=NULL, filter.tags='blood.draw.phase.tag', groups.to.filter='acute', keep=TRUE, color.vals=tmp.cols)
# Output.
publish.plot(tmp.ggplot=tmp.ggplot.1, output.path=fig.4.path, file.name='CovidPt-Acute-pRAssessment', blank.comp=blank.complement.4, do.legend=TRUE, height=5, width=5)
# Report proportions per pR status:
meta.data[blood.draw.phase.tag=='acute'&!is.na(clonotype.tag), .(cell.frac=.N/meta.data[blood.draw.phase.tag=='acute'&!is.na(clonotype.tag), .N]), by=pr.tag]
# pR 0.1963221
# non-pR 0.8036779
meta.data[blood.draw.phase.tag=='acute'&!is.na(clonotype.tag), .N]

# ---> Correlation of pR cell fraction per donor between blood draw phases.
# Get data by keeping fractions only for the donors that were present in both phases of the study.
tmp.data <- meta.data[, .(cell.count=.N), by=.(donor=donor.id.tag, phase=blood.draw.phase.tag)]
tmp.data <- spread(data=tmp.data, key=phase, value=cell.count) # Apparently, we've included only the donors for which we have data frm both phases. We're so smart.
donor.to.keep <- tmp.data[acute>0 & convalescent>0, as.character(donor)]
tmp.data <- meta.data[
  !is.na(clonotype.tag) & donor.id.tag%in%donor.to.keep,
  .(
    acute.freq=.SD[blood.draw.phase.tag=='acute' & pr.tag=='pR', .N],
    acute.prop=.SD[blood.draw.phase.tag=='acute' & pr.tag=='pR', .N]/.SD[blood.draw.phase.tag=='acute', .N],
    conv.freq=.SD[blood.draw.phase.tag=='convalescent' & pr.tag=='pR', .N],
    conv.prop=.SD[blood.draw.phase.tag=='convalescent' & pr.tag=='pR', .N]/.SD[blood.draw.phase.tag=='convalescent', .N]
  ),
  by=.(donor=donor.id.tag, hosp.status=hospitalization.tag)
]
# Proportions.
tmp.ggplot <- ggplot(data=tmp.data, aes(x=conv.prop, y=acute.prop)) +
  geom_jitter(size=7, shape=1, width=0.02, height=0.02) + geom_smooth(formula=y~x, method='lm', se=FALSE, color='black') +
  scale_y_continuous(expand=expansion(add=0.02), breaks=c(0, 0.3, 0.6)) + scale_x_continuous(expand=expansion(add=0.02), breaks=c(0, 0.3, 0.6)) +
  labs(x='Convalescent phase pR cell fraction', y='Acute phase pR cell fraction')
publish.plot(tmp.ggplot=tmp.ggplot, output.path=fig.4.path, file.name='CovidPt-Acute_TCR-Repertoire_pR-RepresentationPerDonor_Acute-vs-Memory_Prop', type='pdf', blank.comp=blank.complement.1, do.legend=FALSE, stat.cor=TRUE)
# Raw numbers.
tmp.ggplot <- ggplot(data=tmp.data, aes(x=conv.freq, y=acute.freq)) +
  geom_point(size=7, shape=1) + geom_smooth(formula=y~x, method='lm', se=FALSE, color='black') +
  scale_y_continuous(expand=expansion(add=0.2), breaks=scales::pretty_breaks(n=3)) + scale_x_continuous(expand=expansion(add=0.2), breaks=scales::pretty_breaks(n=3)) +
  labs(x='Convalescent phase pR cell fraction', y='Acute phase pR cell fraction')
publish.plot(tmp.ggplot=tmp.ggplot, output.path=fig.4.path, file.name='CovidPt-Acute_TCR-Repertoire_pR-RepresentationPerDonor_Acute-vs-Memory_Count', type='pdf', blank.comp=blank.complement.1, do.legend=FALSE, stat.cor=TRUE)
# Logged raw numbers.
tmp.ggplot <- ggplot(data=tmp.data, aes(x=log2(conv.freq+1), y=log2(acute.freq+1))) +
  geom_jitter(size=7, shape=1, width=0.2, height=0.2) + geom_smooth(formula=y~x, method='lm', se=FALSE, color='black') +
  scale_y_continuous(expand=expansion(add=0.2), breaks=c(0, 5, 10)) + scale_x_continuous(expand=expansion(add=0.2), breaks=c(0, 5, 10)) +
  labs(x='Convalescent phase pR cell fraction', y='Acute phase pR cell fraction')
publish.plot(tmp.ggplot=tmp.ggplot, output.path=fig.4.path, file.name='CovidPt-Acute_TCR-Repertoire_pR-RepresentationPerDonor_Acute-vs-Memory_LogCount', type='pdf', blank.comp=blank.complement.1, do.legend=FALSE, stat.cor=TRUE)

# ---> pR profile on UMAP along clone size depiction
tmp.cells <- meta.data[!is.na(clonotype.tag) & blood.draw.phase.tag=='acute', cell]
meta.data[, tmp.clon.size.tag:=clon.size.tag]; meta.data[clon.size.tag>=1500, tmp.clon.size.tag:=1500]
tmp.ggplot <- ggplot(data=meta.data[cell %in% tmp.cells], aes(x=UMAP_1, y=UMAP_2, fill=pr.tag, size=tmp.clon.size.tag)) +
  geom_point(stroke=0.4, shape=21, col='black', alpha=0.8) +
  scale_fill_manual(values=pr.groups.cols) +
  scale_radius(breaks=c(10, 100, 1000), limits=c(1, meta.data[, max(tmp.clon.size.tag, na.rm=TRUE)]), range=c(3, 13)) +
  scale_y_continuous(expand=expansion(add=0)) + scale_x_continuous(expand=expansion(add=0)) +
  labs(x='UMAP 1', y='UMAP 2', col='')
# Output.
publish.plot(tmp.ggplot=tmp.ggplot, output.path=fig.4.path, file.name='CovidPt-Acute_TCR-Repertoire_SharingStatusOnUMAP', type='tiff', blank.comp=blank.complement.3, do.legend=TRUE)

# ---> Contrasting clonal expansion between crossreactive and non-crossreactive TCRs seen in the acute phase.
# Plot and output.
tmp.cols <- sharing.cols; names(tmp.cols) <- c('pR', 'non-pR')
srt.objs.list[[main.obj.4]]@meta.data[, 'log.clon.size.tag'] <- log2(srt.objs.list[[main.obj.4]]@meta.data[, 'clon.size.tag'])
tmp.ggplot <-
vln.plot(seurat.obj=srt.objs.list[[main.obj.4]], feature='log.clon.size.tag', groups.tag='pr.tag', na.rm=TRUE, groups.of.int=NULL, groups.label=NULL, filter.tags='blood.draw.phase.tag', groups.to.filter='acute', keep=TRUE, feature.thold=1, color='median', vln.type='density', size.thold=0, file.name=NULL, adjust.val=2, trim.val=TRUE, plot.limits=NULL, add.points=FALSE, this.color.scale=tmp.cols)
tmp.ggplot <- tmp.ggplot + scale_x_continuous(expand=expansion(add=0), breaks=scales::pretty_breaks(n=3)) + scale_y_continuous(expand=expansion(add=0), breaks=scales::pretty_breaks(n=2))
publish.plot(tmp.ggplot=tmp.ggplot, output.path=fig.4.path, file.name='CovidPt-Acute_TCR-Repertoire_SizeBetween-CrossReactiveStatus', type='pdf', blank.comp=blank.complement.1, do.legend=TRUE)

# ---> Clone size comparison between blood draw phases only for clonotypes intersected between them.
# Keep track of clonotypes that are intersected between blood draw phases (with a given clone size).
clones.of.int <- meta.data[!is.na(donor.id.tag) & !is.na(clonotype.tag),
  .(ovlpd=(.SD[blood.draw.phase.tag=='acute', .N]>=2 & .SD[blood.draw.phase.tag=='convalescent', .N]>=2)),
  by=.(clonotype=clonotype.tag, donor.id.tag)
][ovlpd==TRUE, clonotype]
# Harverst data.
tmp.data <- meta.data[clonotype.tag%chin%clones.of.int & !is.na(clonotype.tag), .(clone.size=.N, pr.tag=unique(pr.tag)), by=.(clonotype=clonotype.tag, draw.phase=blood.draw.phase.tag, donor.id.tag)]
tmp.data <- spread(data=tmp.data, key=draw.phase, value=clone.size, fill=0)
# Remove paris of clonotype-donor when clone size is below the threshold for the corresponding donor.
tmp.data <- tmp.data[acute>=2 & convalescent>=2]
# tmp.clones <- tmp.data[, .N, by=clonotype][N>1, clonotype]; tmp.data[clonotype %in% tmp.clones] # Sanity check. We corroborate there's a single pR class assigned to every clonotype of interest.
# Set log scale for clone size of either phase.
tmp.data[, log.acute:=log2(acute)]; tmp.data[, log.convalescent:=log2(convalescent)]
# Get plot.
tmp.cols <- sharing.cols; names(tmp.cols) <- c('pR', 'non-pR')
tmp.ggplot <- ggplot(data=tmp.data, aes(x=log.convalescent, y=log.acute)) +
  # geom_point(size=3) + #geom_smooth(formula=y~x, method='lm', se=FALSE) +
  geom_jitter(aes(col=pr.tag), size=7, alpha=0.4, width=0.1, height=0.1) +
  geom_smooth(data=tmp.data[pr.tag=='pR'], formula=y~x, method='lm', se=FALSE, color='black') +
  scale_color_manual(values=tmp.cols) + scale_x_continuous(expand=expansion(add=0.2), limits=c(1, tmp.data[, max(log.convalescent)]), breaks=scales::pretty_breaks(n=3)) + scale_y_continuous(expand=expansion(add=0.2), limits=c(1, tmp.data[, max(log.acute)]), breaks=scales::pretty_breaks(n=3)) +
  labs(x='Acute size', y='Convalescent size') + theme_bw()
# Output.
publish.plot(tmp.ggplot=tmp.ggplot, output.path=fig.4.path, file.name='CovidPt-Acute_TCR-Repertoire_SizeStabilityBetweenPhases', blank.comp=blank.complement.1, do.legend=TRUE, stat.cor=TRUE, cor.group='pr.tag')

# ---> CTL-related phenotype (fractions) per cross-reactive group (cross-reactive cells and others).
# Define CD4-CTL tag in seurat object.
tmp.corrs <- c('CD4-CTL', rep(x='Rest', times=6)); names(tmp.corrs) <- c('1', '0', 2:6)
tmp.data <- as.character(srt.objs.list[[main.obj.4]]@meta.data[, clusts.lab[main.obj.4]])
tmp.data <- tmp.corrs[tmp.data]
srt.objs.list[[main.obj.4]]@meta.data[, 'ctl.tag'] <- tmp.data
# seurat.obj@meta.data[, 'ctl.tag'] <- tmp.data
tmp.cols <- c('Rest'='#808080', 'CD4-CTL'='#cc0000')
# Obtain plots.
tmp.ggplot.1 <- plot.props(seurat.obj=srt.objs.list[[main.obj.4]], plot.type='donut', groups.tag='ctl.tag', na.rm=TRUE, groups.of.int=NULL, groups.label=NULL, filter.tags=c('pr.tag', 'blood.draw.phase.tag'), groups.to.filter=list('pR', 'acute'), keep=TRUE, color.vals=tmp.cols)
tmp.ggplot.2 <- plot.props(seurat.obj=srt.objs.list[[main.obj.4]], plot.type='donut', groups.tag='ctl.tag', na.rm=TRUE, groups.of.int=NULL, groups.label=NULL, filter.tags=c('pr.tag', 'blood.draw.phase.tag'), groups.to.filter=list('non-pR', 'acute'), keep=TRUE, color.vals=tmp.cols)
# Output.
publish.plot(tmp.ggplot=tmp.ggplot.1, output.path=fig.4.path, file.name='CovidPt-Acute_CD4-CTLs-and-pRAssessment-pR', blank.comp=blank.complement.4, do.legend=TRUE, height=5, width=5)
publish.plot(tmp.ggplot=tmp.ggplot.2, output.path=fig.4.path, file.name='CovidPt-Acute_CD4-CTLs-and-pRAssessment-nonpR', blank.comp=blank.complement.4, do.legend=TRUE, height=5, width=5)

# ---> pR profile (fraction) for CD4-CTLs and TFHs
# Define CD4-CTL tag in seurat object.
tmp.corrs <- c('CTL', rep(x='Rest', times=6)); names(tmp.corrs) <- c('1', '0', 2:6)
tmp.data <- as.character(srt.objs.list[[main.obj.4]]@meta.data[, clusts.lab[main.obj.4]])
tmp.data <- tmp.corrs[tmp.data]
srt.objs.list[[main.obj.4]]@meta.data[, 'ctl.tag'] <- tmp.data
# Define TFH tag in seurat object.
tmp.corrs <- c(rep(x='TFH', times=2), rep(x='Rest', times=5)); names(tmp.corrs) <- c('0', '5', 1:4, '6')
tmp.data <- as.character(srt.objs.list[[main.obj.4]]@meta.data[, clusts.lab[main.obj.4]])
tmp.data <- tmp.corrs[tmp.data]
srt.objs.list[[main.obj.4]]@meta.data[, 'tfh.tag'] <- tmp.data
# seurat.obj@meta.data[, 'ctl.tag'] <- tmp.data
tmp.cols <- sharing.cols; names(tmp.cols) <- c('pR', 'non-pR')
# Obtain plots.
tmp.ggplot.1 <- plot.props(seurat.obj=srt.objs.list[[main.obj.4]], plot.type='donut', groups.tag='pr.tag', na.rm=TRUE, groups.of.int=NULL, groups.label=NULL, filter.tags=c('ctl.tag', 'blood.draw.phase.tag'), groups.to.filter=list('CTL', 'acute'), keep=TRUE, color.vals=tmp.cols)
tmp.ggplot.2 <- plot.props(seurat.obj=srt.objs.list[[main.obj.4]], plot.type='donut', groups.tag='pr.tag', na.rm=TRUE, groups.of.int=NULL, groups.label=NULL, filter.tags=c('tfh.tag', 'blood.draw.phase.tag'), groups.to.filter=list('TFH', 'acute'), keep=TRUE, color.vals=tmp.cols)
# Output.
publish.plot(tmp.ggplot=tmp.ggplot.1, output.path=fig.4.path, file.name='CovidPt-Acute_pRFractions-CD4-CTL', blank.comp=blank.complement.4, do.legend=TRUE, height=5, width=5)
publish.plot(tmp.ggplot=tmp.ggplot.2, output.path=fig.4.path, file.name='CovidPt-Acute_pRFractions-TFH', blank.comp=blank.complement.4, do.legend=TRUE, height=5, width=5)


### -------------------- Supplementary Figure 4 --------------------- ###
# ----> Output directory.
sup.4.path <- paste0(reports.path, '/sup_figure_4')
if(!dir.exists(sup.4.path)) dir.create(sup.4.path)
# ---> Comments: Main dataset for this section, PT-24-ALL-SF

# ---> Cross-reactive TCR representation per donor.
# Get data.
tmp.data <- sec.meta.data[
  !is.na(clonotype.tag) & blood.draw.phase.tag=='acute',
  .(
    pr.prop=.SD[pr.tag=='pR', .N]/.N
  ),
  by=.(donor=donor.id.tag, hosp.status=hospitalization.tag)
]
setorderv(x=tmp.data, cols='pr.prop', order=-1)
tmp.data$donor <- factor(x=as.character(tmp.data$donor), levels=as.character(tmp.data$donor))
# Plot and output.
# Plot and output.
tmp.ggplot <- ggplot(data=tmp.data, aes(x=donor, y=pr.prop)) +
  geom_bar(stat='identity', fill='lightgray', width=0.7) +
  scale_y_continuous(limits=c(0, 1), expand=expansion(add=0), breaks=scales::pretty_breaks(n=2)) + scale_x_discrete(expand=expansion(add=0.5)) +
  labs(x='Donors ranked according to crossreactivity representation', y='Percent of total cells') +
  theme(axis.ticks.x=element_blank())
publish.plot(tmp.ggplot=tmp.ggplot, output.path=sup.4.path, file.name='CovidPt-Acute-24_TCR-Repertoire_pR-RepresentationPerDonor', type='pdf', blank.comp=blank.complement.1, do.legend=TRUE, width=10)

# ---> pR profile for shared TCRs between phases.
# Obtain plots.
tmp.cols <- sharing.cols; names(tmp.cols) <- c('pR', 'non-pR')
tmp.ggplot.1 <- plot.props(seurat.obj=srt.objs.list[[sec.obj.4]], plot.type='pie', groups.tag='pr.tag', na.rm=TRUE, groups.of.int=NULL, groups.label=NULL, filter.tags='blood.draw.phase.tag', groups.to.filter='acute', keep=TRUE, color.vals=tmp.cols)
# Output.
publish.plot(tmp.ggplot=tmp.ggplot.1, output.path=sup.4.path, file.name='CovidPt-Acute-24-pRAssessment', blank.comp=blank.complement.4, do.legend=TRUE, height=5, width=5)
# Report proportions per pR status:
sec.meta.data[blood.draw.phase.tag=='acute'&!is.na(clonotype.tag), .(cell.frac=.N/sec.meta.data[blood.draw.phase.tag=='acute'&!is.na(clonotype.tag), .N]), by=pr.tag]
# pR 0.2147186
# non-pR 0.7852814

# ---> Correlation of pR cell fraction per donor between blood draw phases.
# Get data by keeping fractions only for the donors that were present in both phases of the study.
tmp.data <- sec.meta.data[, .(cell.count=.N), by=.(donor=donor.id.tag, phase=blood.draw.phase.tag)]
tmp.data <- spread(data=tmp.data, key=phase, value=cell.count) # Apparently, we've included only the donors for which we have data frm both phases. We're so smart.
donor.to.keep <- tmp.data[acute>0 & convalescent>0, as.character(donor)]
tmp.data <- sec.meta.data[
  !is.na(clonotype.tag) & donor.id.tag%in%donor.to.keep,
  .(
    acute.freq=.SD[blood.draw.phase.tag=='acute' & pr.tag=='pR', .N],
    acute.prop=.SD[blood.draw.phase.tag=='acute' & pr.tag=='pR', .N]/.SD[blood.draw.phase.tag=='acute', .N],
    conv.freq=.SD[blood.draw.phase.tag=='convalescent' & pr.tag=='pR', .N],
    conv.prop=.SD[blood.draw.phase.tag=='convalescent' & pr.tag=='pR', .N]/.SD[blood.draw.phase.tag=='convalescent', .N]
  ),
  by=.(donor=donor.id.tag, hosp.status=hospitalization.tag)
]
# Proportions.
tmp.ggplot <- ggplot(data=tmp.data, aes(x=conv.prop, y=acute.prop)) +
  geom_jitter(size=7, shape=1, width=0.02, height=0.02) + geom_smooth(formula=y~x, method='lm', se=FALSE, color='black') +
  scale_y_continuous(expand=expansion(add=0.02), breaks=c(0, 0.3, 0.6)) + scale_x_continuous(expand=expansion(add=0.02), breaks=c(0, 0.3, 0.6)) +
  labs(x='Convalescent phase pR cell fraction', y='Acute phase pR cell fraction')
publish.plot(tmp.ggplot=tmp.ggplot, output.path=sup.4.path, file.name='CovidPt-Acute-24_TCR-Repertoire_pR-RepresentationPerDonor_Acute-vs-Memory_Prop', type='pdf', blank.comp=blank.complement.1, do.legend=FALSE, stat.cor=TRUE)
# Raw numbers.
tmp.ggplot <- ggplot(data=tmp.data, aes(x=conv.freq, y=acute.freq)) +
  geom_point(size=7, shape=1) + geom_smooth(formula=y~x, method='lm', se=FALSE, color='black') +
  scale_y_continuous(expand=expansion(add=0.2), breaks=scales::pretty_breaks(n=3)) + scale_x_continuous(expand=expansion(add=0.2), breaks=scales::pretty_breaks(n=3)) +
  labs(x='Convalescent phase pR cell fraction', y='Acute phase pR cell fraction')
publish.plot(tmp.ggplot=tmp.ggplot, output.path=sup.4.path, file.name='CovidPt-Acute-24_TCR-Repertoire_pR-RepresentationPerDonor_Acute-vs-Memory_Count', type='pdf', blank.comp=blank.complement.1, do.legend=FALSE, stat.cor=TRUE)
# Logged raw numbers.
tmp.ggplot <- ggplot(data=tmp.data, aes(x=log2(conv.freq+1), y=log2(acute.freq+1))) +
  geom_jitter(size=7, shape=1, width=0.2, height=0.2) + geom_smooth(formula=y~x, method='lm', se=FALSE, color='black') +
  scale_y_continuous(expand=expansion(add=0.2), breaks=c(0, 5, 10)) + scale_x_continuous(expand=expansion(add=0.2), breaks=c(0, 5, 10)) +
  labs(x='Convalescent phase pR cell fraction', y='Acute phase pR cell fraction')
publish.plot(tmp.ggplot=tmp.ggplot, output.path=sup.4.path, file.name='CovidPt-Acute-24_TCR-Repertoire_pR-RepresentationPerDonor_Acute-vs-Memory_LogCount', type='pdf', blank.comp=blank.complement.1, do.legend=FALSE, stat.cor=TRUE)

# ---> pR profile on UMAP along clone size depiction
tmp.cells <- sec.meta.data[!is.na(clonotype.tag) & blood.draw.phase.tag=='acute', cell]
sec.meta.data[, tmp.clon.size.tag:=clon.size.tag]; sec.meta.data[clon.size.tag>=1500, tmp.clon.size.tag:=1500]
tmp.ggplot <- ggplot(data=sec.meta.data[cell %in% tmp.cells], aes(x=UMAP_1, y=UMAP_2, fill=pr.tag, size=tmp.clon.size.tag)) +
  geom_point(stroke=0.4, shape=21, col='black', alpha=0.8) +
  scale_fill_manual(values=pr.groups.cols) +
  scale_radius(breaks=c(10, 100, 1000), limits=c(1, sec.meta.data[, max(tmp.clon.size.tag, na.rm=TRUE)]), range=c(3, 13)) +
  scale_y_continuous(expand=expansion(add=0)) + scale_x_continuous(expand=expansion(add=0), limits=c(NA, 5)) +
  labs(x='UMAP 1', y='UMAP 2', col='')
# Output.
publish.plot(tmp.ggplot=tmp.ggplot, output.path=sup.4.path, file.name='CovidPt-Acute-24_TCR-Repertoire_SharingStatusOnUMAP', type='tiff', blank.comp=blank.complement.3, do.legend=TRUE)

# ---> Contrasting clonal expansion between crossreactive and non-crossreactive TCRs seen in the acute phase.
# Plot and output.
tmp.cols <- sharing.cols; names(tmp.cols) <- c('pR', 'non-pR')
srt.objs.list[[sec.obj.4]]@meta.data[, 'log.clon.size.tag'] <- log2(srt.objs.list[[sec.obj.4]]@meta.data[, 'clon.size.tag'])
tmp.ggplot <-
vln.plot(seurat.obj=srt.objs.list[[sec.obj.4]], feature='log.clon.size.tag', groups.tag='pr.tag', na.rm=TRUE, groups.of.int=NULL, groups.label=NULL, filter.tags='blood.draw.phase.tag', groups.to.filter='acute', keep=TRUE, feature.thold=1, color='median', vln.type='density', size.thold=0, file.name=NULL, adjust.val=2, trim.val=TRUE, plot.limits=NULL, add.points=FALSE, this.color.scale=tmp.cols)
tmp.ggplot <- tmp.ggplot + scale_x_continuous(expand=expansion(add=0), breaks=scales::pretty_breaks(n=2)) + scale_y_continuous(expand=expansion(add=0), breaks=scales::pretty_breaks(n=2))
publish.plot(tmp.ggplot=tmp.ggplot, output.path=sup.4.path, file.name='CovidPt-Acute-24_TCR-Repertoire_SizeBetween-CrossReactiveStatus', type='pdf', blank.comp=blank.complement.1, do.legend=TRUE)

# ---> Clone size comparison between blood draw phases only for clonotypes intersected between them.
# Keep track of clonotypes that are intersected between blood draw phases (with a given clone size).
clones.of.int <- sec.meta.data[!is.na(donor.id.tag) & !is.na(clonotype.tag),
  .(ovlpd=(.SD[blood.draw.phase.tag=='acute', .N]>=2 & .SD[blood.draw.phase.tag=='convalescent', .N]>=2)),
  by=.(clonotype=clonotype.tag, donor.id.tag)
][ovlpd==TRUE, clonotype]
# Harverst data.
tmp.data <- sec.meta.data[clonotype.tag%chin%clones.of.int & !is.na(clonotype.tag), .(clone.size=.N, pr.tag=unique(pr.tag)), by=.(clonotype=clonotype.tag, draw.phase=blood.draw.phase.tag, donor.id.tag)]
tmp.data <- spread(data=tmp.data, key=draw.phase, value=clone.size, fill=0)
# Remove paris of clonotype-donor when clone size is below the threshold for the corresponding donor.
tmp.data <- tmp.data[acute>=2 & convalescent>=2]
# tmp.clones <- tmp.data[, .N, by=clonotype][N>1, clonotype]; tmp.data[clonotype %in% tmp.clones] # Sanity check. We corroborate there's a single pR class assigned to every clonotype of interest.
# Set log scale for clone size of either phase.
tmp.data[, log.acute:=log2(acute)]; tmp.data[, log.convalescent:=log2(convalescent)]
# Get plot.
tmp.cols <- sharing.cols; names(tmp.cols) <- c('pR', 'non-pR')
tmp.ggplot <- ggplot(data=tmp.data, aes(x=log.convalescent, y=log.acute)) +
  # geom_point(size=3) + #geom_smooth(formula=y~x, method='lm', se=FALSE) +
  geom_jitter(aes(col=pr.tag), size=7, alpha=0.4, width=0.1, height=0.1) +
  geom_smooth(data=tmp.data[pr.tag=='pR'], formula=y~x, method='lm', se=FALSE, color='black') +
  scale_color_manual(values=tmp.cols) + scale_x_continuous(expand=expansion(add=0.2), limits=c(1, tmp.data[, max(log.convalescent)]), breaks=scales::pretty_breaks(n=3)) + scale_y_continuous(expand=expansion(add=0.2), limits=c(1, tmp.data[, max(log.acute)]), breaks=scales::pretty_breaks(n=3)) +
  labs(x='Acute size', y='Convalescent size') + theme_bw()
# Output.
publish.plot(tmp.ggplot=tmp.ggplot, output.path=sup.4.path, file.name='CovidPt-Acute-24_TCR-Repertoire_SizeStabilityBetweenPhases', blank.comp=blank.complement.1, do.legend=TRUE, stat.cor=TRUE, cor.group='pr.tag')

# ---> CTL-related phenotype (fractions) per cross-reactive group (cross-reactive cells and others).
# Define CD4-CTL tag in seurat object.
tmp.corrs <- c(rep(x='CD4-CTL', times=2), rep(x='Rest', times=3)); names(tmp.corrs) <- c('1', '4', '0', 2:3)
tmp.data <- as.character(srt.objs.list[[sec.obj.4]]@meta.data[, clusts.lab[sec.obj.4]])
tmp.data <- tmp.corrs[tmp.data]
srt.objs.list[[sec.obj.4]]@meta.data[, 'ctl.tag'] <- tmp.data
tmp.cols <- c('Rest'='#808080', 'CD4-CTL'='#cc0000')
# Obtain plots.
tmp.ggplot.1 <- plot.props(seurat.obj=srt.objs.list[[sec.obj.4]], plot.type='donut', groups.tag='ctl.tag', na.rm=TRUE, groups.of.int=NULL, groups.label=NULL, filter.tags=c('pr.tag', 'blood.draw.phase.tag'), groups.to.filter=list('pR', 'acute'), keep=TRUE, color.vals=tmp.cols)
tmp.ggplot.2 <- plot.props(seurat.obj=srt.objs.list[[sec.obj.4]], plot.type='donut', groups.tag='ctl.tag', na.rm=TRUE, groups.of.int=NULL, groups.label=NULL, filter.tags=c('pr.tag', 'blood.draw.phase.tag'), groups.to.filter=list('non-pR', 'acute'), keep=TRUE, color.vals=tmp.cols)
# Output.
publish.plot(tmp.ggplot=tmp.ggplot.1, output.path=sup.4.path, file.name='CovidPt-Acute-24_CD4-CTLs-and-pRAssessment-pR', blank.comp=blank.complement.4, do.legend=TRUE, height=5, width=5)
publish.plot(tmp.ggplot=tmp.ggplot.2, output.path=sup.4.path, file.name='CovidPt-Acute-24_CD4-CTLs-and-pRAssessment-nonpR', blank.comp=blank.complement.4, do.legend=TRUE, height=5, width=5)

# ---> pR profile (fraction) for CD4-CTLs and TFHs
# Define CD4-CTL tag in seurat object.
tmp.corrs <- c(rep(x='CD4-CTL', times=2), rep(x='Rest', times=3)); names(tmp.corrs) <- c('1', '4', '0', 2:3)
tmp.data <- as.character(srt.objs.list[[sec.obj.4]]@meta.data[, clusts.lab[sec.obj.4]])
tmp.data <- tmp.corrs[tmp.data]
srt.objs.list[[sec.obj.4]]@meta.data[, 'ctl.tag'] <- tmp.data
# Define TFH tag in seurat object.
tmp.corrs <- c(rep(x='TREG', times=2), rep(x='Rest', times=3)); names(tmp.corrs) <- c('0', '3', 1:2, '4')
tmp.data <- as.character(srt.objs.list[[sec.obj.4]]@meta.data[, clusts.lab[sec.obj.4]])
tmp.data <- tmp.corrs[tmp.data]
srt.objs.list[[sec.obj.4]]@meta.data[, 'treg.tag'] <- tmp.data
# seurat.obj@meta.data[, 'ctl.tag'] <- tmp.data
tmp.cols <- sharing.cols; names(tmp.cols) <- c('pR', 'non-pR')
# Obtain plots.
tmp.ggplot.1 <- plot.props(seurat.obj=srt.objs.list[[sec.obj.4]], plot.type='donut', groups.tag='pr.tag', na.rm=TRUE, groups.of.int=NULL, groups.label=NULL, filter.tags=c('ctl.tag', 'blood.draw.phase.tag'), groups.to.filter=list('CD4-CTL', 'acute'), keep=TRUE, color.vals=tmp.cols)
tmp.ggplot.2 <- plot.props(seurat.obj=srt.objs.list[[sec.obj.4]], plot.type='donut', groups.tag='pr.tag', na.rm=TRUE, groups.of.int=NULL, groups.label=NULL, filter.tags=c('treg.tag', 'blood.draw.phase.tag'), groups.to.filter=list('TREG', 'acute'), keep=TRUE, color.vals=tmp.cols)
# Output.
publish.plot(tmp.ggplot=tmp.ggplot.1, output.path=sup.4.path, file.name='CovidPt-Acute-24_pRFractions-CD4-CTL', blank.comp=blank.complement.4, do.legend=TRUE, height=5, width=5)
publish.plot(tmp.ggplot=tmp.ggplot.2, output.path=sup.4.path, file.name='CovidPt-Acute-24_pRFractions-TREG', blank.comp=blank.complement.4, do.legend=TRUE, height=5, width=5)


### ------------- Other things that may fit in Figure 4 ------------- ###

# ---> Clonotype sharing between blood draw phases.
tmp.cols <- blood.phases.cols; names(tmp.cols) <- tolower(names(tmp.cols))
depict.clone.sharing(meta.data=meta.data, tag.of.int='blood.draw.phase.tag', reports.path=fig.4.path, file.preffix='CovidPt-Acute_TCR-Repertoire_ClonotypeSharingBetweenPhases', fill.cols=tmp.cols)

# ---> Phenotype stability for intersected clonotypes between phases.
# Define phenotypes of interest to depict in the plot.
clusts.of.int <- c('0'='TFH', '1'='CD4-CTL') # cTFH and CD4-CTL
# Keep track of clonotypes that are intersected between blood draw phases (with a given clone size) taking into account only the populations of interest.
clones.of.int <- meta.data[
  as.character(get(clusts.lab[main.obj.4])) %chin% names(clusts.of.int),
  .(ovlpd=(.SD[blood.draw.phase.tag=='acute', .N]>=2 & .SD[blood.draw.phase.tag=='convalescent', .N]>=2)),
  by=.(clonotype=clonotype.tag)
][ovlpd==TRUE, clonotype]
# Harvest subset of data.
tmp.data <- meta.data[clonotype.tag%chin%clones.of.int & !is.na(clonotype.tag) & as.character(get(clusts.lab[main.obj.1])) %chin% names(clusts.of.int)]
tmp.data <- tmp.data[,
  .(clone.size=.N),
  by=.(
    clonotype=clonotype.tag,
    cluster=as.character(get(clusts.lab[main.obj.1])),
    draw.phase=blood.draw.phase.tag
  )]
tmp.data[, cluster:=clusts.of.int[cluster]]
# Tidy data.
tmp.data <- spread(data=tmp.data, key=cluster, value=clone.size, fill=0)
tmp.data[, clusts.sizes:=paste(`CD4-CTL`, `TFH`, sep=';')]; tmp.data[, `:=`(`CD4-CTL`=NULL, `TFH`=NULL)]
tmp.data <- spread(data=tmp.data, key=draw.phase, value=clusts.sizes, fill='0;0')
for(tmp.val in c('acute', 'convalescent')){
  tmp.data <- separate(data=tmp.data, col=tmp.val, into=paste(c('CD4-CTL', 'TFH'), tmp.val, sep='.'), sep=';', convert=TRUE)
}
# Identify zero values previous to calculating z-scores.
tmp.rows <- tmp.data[, clonotype]; tmp.data[, clonotype:=NULL]
tmp.data <- as.matrix(as.matrix(tmp.data)); row.names(tmp.data) <- tmp.rows
tmp.cols <- sort(colnames(tmp.data)); tmp.data <- tmp.data[, tmp.cols]
zero.track <- tmp.data==0
# Calulate z-scores.
tmp.data <- t(scale(x=t(tmp.data)))
# Set column metadata.
col.metadata <- data.frame(
  row.names=tmp.cols,
  Population=ifelse(test=str_detect(string=tmp.cols, pattern='CD4-CTL'), yes='CD4-CTL', no='TFH'),
  `Blood draw phase`=ifelse(test=str_detect(string=tmp.cols, pattern='acute'), yes='Acute', no='Convalescent')
)
# Set row metadata.
row.metadata.1 <- as.data.frame(meta.data[!is.na(clonotype.tag) & clonotype.tag%in%row.names(tmp.data), .(Sharing.Status=.SD[, .N, by=pr.tag][N==max(N), pr.tag]), by=clonotype.tag], stringsAsFactors=FALSE)
row.metadata.2 <- data.frame(
  row.names=row.names(tmp.data),
  clonotype.tag=row.names(tmp.data)
)
row.metadata <- merge(x=row.metadata.2, y=row.metadata.1, by='clonotype.tag', sort=FALSE)
if(all(row.metadata$clonotype.tag==row.names(tmp.data))){ row.names(row.metadata) <- row.metadata$clonotype.tag; row.metadata$clonotype.tag <- NULL }
# Define colors for metadata tracks.
pop.cols <- clusts.cols[[main.obj.1]][names(clusts.of.int)]; names(pop.cols) <- clusts.of.int[names(pop.cols)]
tmp.cols <- sharing.cols; names(tmp.cols) <- c('pR', 'non-pR')
ann.colors <- list(
  Population=pop.cols,
  Blood.draw.phase=blood.phases.cols,
  Sharing.Status=tmp.cols
)
# Set color scale and breaks for heatmap.
col.breaks <- seq(from=min(range(tmp.data)), to=max(range(tmp.data)), length.out=100)
mid.point <- which.min(abs(col.breaks - 0))
hmap.col.scale.1 <- colorRampPalette(c('blue', 'mediumblue', 'black'))(mid.point)
hmap.col.scale.2 <- colorRampPalette(c('black', 'gold', 'yellow'))(100-(mid.point+1))
hmap.col.scale <- c(hmap.col.scale.1, hmap.col.scale.2)
# Set order of rows according to hierarchical clustering.
dist.mat <- dist(x=tmp.data, method='euclidean')
tmp.dgrm <- hclust(d=dist.mat, method='complete')
# Replace scaled values by NA for those entries that originally were 0.
tmp.data[zero.track] <- NA
# Plot 'complete' version
tmp.file.name <- paste0(fig.4.path, '/CovidPt_TCR-Acute-RepertoirePhenotypeStability.C.pdf')
pheatmap(mat=tmp.data, color=hmap.col.scale, breaks=col.breaks, scale='none', cluster_rows=tmp.dgrm, cluster_cols=FALSE, annotation_col=col.metadata, annotation_row=row.metadata, annotation_colors=ann.colors, na_col='#808080', show_colnames=FALSE, show_rownames=FALSE, filename=tmp.file.name)
# Plot 'blank' version
tmp.file.name <- paste0(fig.4.path, '/CovidPt_TCR-Acute-RepertoirePhenotypeStability.B.pdf')
pheatmap(mat=tmp.data, color=hmap.col.scale, breaks=col.breaks, scale='none', cluster_rows=tmp.dgrm, cluster_cols=FALSE, annotation_col=col.metadata, annotation_row=row.metadata, annotation_colors=ann.colors, legend=FALSE, annotation_legend=FALSE, annotation_names_row=FALSE, annotation_names_col=FALSE, show_rownames=FALSE, show_colnames=FALSE, filename=tmp.file.name)
# Plot legend.
tmp.ggplot <- ggplot(data=as.data.frame(tmp.data), aes(x=`CD4-CTL.acute`, y=`CD4-CTL.acute`, col=`CD4-CTL.acute`)) + geom_point() + scale_color_gradientn(colors=hmap.col.scale, breaks=col.breaks, name=NULL, labels=NULL); tmp.ggplot <- get_legend(p=tmp.ggplot)
tmp.file.name <- paste0(fig.4.path, '/CovidPt_TCR-Acute-RepertoirePhenotypeStability.L.pdf')
pdf(file=tmp.file.name)
print(as_ggplot(tmp.ggplot))
dev.off()

# ---> Comparison between SARS-reactive cell count and cross-reactive potential -crossreactive cell count- (at the donor level).
# Get data.
tmp.data.1 <- meta.data[
  !is.na(clonotype.tag) & blood.draw.phase.tag=='acute',
  .(
    pr.freq=.SD[pr.tag=='pR', .N],
    pr.prop=.SD[pr.tag=='pR', .N]/.N
  ),
  by=.(donor=donor.id.tag, hosp.status=hospitalization.tag)
]
tmp.data.2 <- sr.cell.counts[time.point=='6']
tmp.data.2[, hosp.status:=ifelse(test=hosp.status=='Hosp', yes='Yes', no='No')]
tmp.data <- merge(x=tmp.data.1, y=tmp.data.2, by.x='donor', by.y='patient.id')
# tmp.data[, all(hosp.status.x==hosp.status.y)] # Sanity check. Annotations may maintain stable across files.
# Plot and output.
tmp.ggplot <- ggplot(data=tmp.data, aes(x=cells.acute, y=pr.prop)) +
  geom_point(size=3) + geom_smooth(formula=y~x, method='lm', se=FALSE, color='black') +
  scale_y_continuous(expand=expansion(add=0.01)) + scale_x_continuous(expand=expansion(add=0.01)) +
  labs(x='SARS-CoV-2-reactive cell count per 10e6 PBMC', y='Crossreactive cell fraction')
publish.plot(tmp.ggplot=tmp.ggplot, output.path=fig.4.path, file.name='CovidPt_TCR-Acute-Repertoire_SARSCellCount-vs-pR-RepresentationPerDonor', type='pdf', blank.comp=blank.complement.1, do.legend=FALSE, stat.cor=TRUE)


############    -----------------------------------------    ############
### --------------------------- Figure 5 ---------------------------- ###
############    -----------------------------------------    ############

# ---> General definitions.
# Seurat objects of interest.
main.obj.5 <- 'HLTY-6-ExpA-FCE'
sec.obj.5 <- 'HLTY-6-ExpD-S3'
# Seurat objects' meta data
meta.data <- cbind(
  data.table(cell=Cells(srt.objs.list[[main.obj.5]])),
  srt.objs.list[[main.obj.5]]@meta.data,
  srt.objs.list[[main.obj.5]]@reductions$umap@cell.embeddings
)
sec.meta.data <- cbind(
  data.table(cell=Cells(srt.objs.list[[sec.obj.5]])),
  srt.objs.list[[sec.obj.5]]@meta.data,
  srt.objs.list[[sec.obj.5]]@reductions$umap@cell.embeddings
)

# ---> Whole list of healthy donors considered in scRNA-seq
mixedsort(unique(c(meta.data[, unique(as.character(donor.id.tag))], sec.meta.data[, unique(as.character(donor.id.tag))])), decreasing=TRUE)


### ------------------------- Main Figure 5 ------------------------- ###

# ---> Output directory.
fig.5.path <- paste0(reports.path, '/figure_5')
if(!dir.exists(fig.5.path)) dir.create(fig.5.path)
# ---> Comments: Main dataset for this section, HLTY-6-ExpA-FCE & HLTY-6-ExpD-S3

# ---> Number of SARS-CoV-2- and HCoV-reactive cells per 10e6 PBMCs in healthy donors.
tmp.cols <- virus.cols[c('SARS-CoV-2', 'HCoV')]
tmp.ggplot <- ggplot(data=hlty.hcov.prof, aes(x=virus, y=cells)) +
  geom_jitter(aes(col=virus), size=8, shape=21, alpha=1, width=0.1) + scale_color_manual(values=tmp.cols) +
  geom_line(aes(group=donor.id.tag), size=0.8, color='black') +
  labs(x='', y='Cells per 10e6 PBMCs', col='Virus')
publish.plot(tmp.ggplot=tmp.ggplot, output.path=fig.5.path, file.name='HltyDonor-ExpD_VirusReactiveCellProfile', type='pdf', blank.comp=blank.complement.1, do.legend=TRUE, width=4)

# ---> Clonotype sharing between virus reactivities, experiment D.
tmp.cols <- virus.cols[c('HCoV', 'SARS-CoV-2')]
depict.clone.sharing(meta.data=sec.meta.data, tag.of.int='virus.tag', reports.path=fig.5.path, file.preffix='HltyDonor-ExpD_TCR-Repertoire_ClonotypeSharingBetweenReactivities', fill.cols=tmp.cols)
sec.meta.data[!is.na(clonotype.tag), .N, by=.(virus.tag, clonotype.tag)][N>1][, uniqueN(clonotype.tag), by=virus.tag]

# ---> Overall cross-reactive TCR representation, experiment D.
# Obtain plots.
tmp.cols <- sharing.cols; names(tmp.cols) <- c('pR', 'non-pR')
tmp.ggplot.1 <- plot.props(seurat.obj=srt.objs.list[[sec.obj.5]], plot.type='pie', groups.tag='pr.tag', na.rm=TRUE, groups.of.int=NULL, groups.label=NULL, filter.tags=NULL, groups.to.filter=NULL, keep=TRUE, color.vals=tmp.cols)
# Output.
publish.plot(tmp.ggplot=tmp.ggplot.1, output.path=fig.5.path, file.name='HltyDonor-ExpD_TCR-Repertoire_pR-RepresentationOverall', blank.comp=blank.complement.4, do.legend=TRUE, height=5, width=5)
# Report proportions per pR status:
sec.meta.data[!is.na(clonotype.tag), .(cell.frac=.N/sec.meta.data[!is.na(clonotype.tag), .N]), by=pr.tag]
# pR 0.05
# non-pR 0.95

# ---> Cross-reactive TCR representation per donor, experiment D.
# Get data.
tmp.data <- sec.meta.data[
  !is.na(clonotype.tag) & !is.na(donor.id.tag),
  .(
    pr.prop=.SD[pr.tag=='pR', .N]/.N
  ),
  by=.(donor=donor.id.tag)
]
tmp.data[, donor:=paste0('D', str_extract(string=donor, pattern='\\d{2}$'))]
setorderv(x=tmp.data, cols='pr.prop', order=-1)
tmp.data$donor <- factor(x=as.character(tmp.data$donor), levels=as.character(tmp.data$donor))
# Plot and output.
tmp.ggplot <- ggplot(data=tmp.data, aes(x=donor, y=pr.prop)) +
  geom_bar(stat='identity', fill='lightgray', width=0.7) +
  scale_y_continuous(limits=c(0, 1), expand=expansion(add=0), breaks=scales::pretty_breaks(n=3)) + scale_x_discrete(expand=expansion(add=0.5)) +
  labs(x='Donors ranked according to crossreactivity representation', y='Percent of total cells') +
  theme(axis.ticks.x=element_blank())
publish.plot(tmp.ggplot=tmp.ggplot, output.path=fig.5.path, file.name='HltyDonor-ExpD_TCR-Repertoire_pR-RepresentationPerDonor', type='pdf', blank.comp=blank.complement.1, do.legend=TRUE, width=6)

# ---> Contrasting clonal expansion between crossreactive and non-crossreactive TCRs.
# Plot and output.
tmp.cols <- sharing.cols; names(tmp.cols) <- c('pR', 'non-pR')
srt.objs.list[[sec.obj.5]]@meta.data[, 'log.clon.size.tag'] <- log2(srt.objs.list[[sec.obj.5]]@meta.data[, 'clon.size.tag'])
tmp.ggplot <-
vln.plot(seurat.obj=srt.objs.list[[sec.obj.5]], feature='log.clon.size.tag', groups.tag='pr.tag', na.rm=TRUE, groups.of.int=NULL, groups.label=NULL, filter.tags=NULL, groups.to.filter=NULL, keep=TRUE, feature.thold=1, color='median', vln.type='density', size.thold=0, file.name=NULL, adjust.val=2, trim.val=TRUE, plot.limits=NULL, add.points=FALSE, this.color.scale=tmp.cols)
tmp.ggplot <- tmp.ggplot + scale_x_continuous(expand=expansion(add=0), breaks=scales::pretty_breaks(n=2)) + scale_y_continuous(expand=expansion(add=0), breaks=scales::pretty_breaks(n=2))
publish.plot(tmp.ggplot=tmp.ggplot, output.path=fig.5.path, file.name='HltyDonor-ExpD_TCR-Repertoire_SizeBetween-pR-and-nonpR', type='pdf', blank.comp=blank.complement.1, do.legend=TRUE, width=4, height=4)

# ---> CD4-CTL profile in healthy donors.
# Info of CD4-CTL profile per donor assessed.
tmp.data.1 <- hlty.cyt.prof[, .(donor.id, CMV=mem.cd154.pop.cmv.prop*100, `CMV-reactive`=cd4.ctl.cmv.prop*100, `Total`=cd4.ctl.total.prop*100)]
tmp.data.1 <- gather(data=tmp.data.1, -`donor.id`, -`CMV`, key='group', value='cell.prop')
# Info of donors that were picked for further analyses after assessment.
tmp.data.2 <- data.table(donor.id=paste0('SDBB-', c('088', '085', '020', '039', '065', '070')), picked='Yes')
tmp.data.1 <- merge(x=tmp.data.2, y=tmp.data.1, by='donor.id', all=TRUE)
tmp.data.1[is.na(picked), picked:='No']
# Plot.
# tmp.sizes <- c(Yes=8, No=4)
tmp.cols <- c(Yes='#00008b', No=unname(virus.cols[c('CMV')]))
tmp.alphas <- c(Yes=1, No=0.2)
tmp.ggplot <- ggplot(data=tmp.data.1[group=='CMV-reactive'], aes(x=CMV, y=cell.prop)) +
  geom_point(aes(color=picked, fill=picked, alpha=picked), shape=21, size=8) + geom_smooth(method='lm', se=FALSE, size=0.6) +
  scale_color_manual(values=tmp.cols) +
  scale_fill_manual(values=tmp.cols) +
  scale_alpha_manual(values=tmp.alphas) +
  scale_y_continuous(expand=expansion(add=c(1.5, 1.5)), breaks=c(30, 60)) + scale_x_continuous(expand=expansion(add=c(0.07, 0.07)), breaks=c(1, 3)) +
  labs(x='% CMV-reactive cells', y='% CD4-CTLs', col='')
# Output.
publish.plot(tmp.ggplot=tmp.ggplot, output.path=fig.5.path, file.name='HltyDonor-ExpA_CD4-CTLProfile_CMVReactiveCells', blank.comp=blank.complement.1, do.legend=TRUE, stat.cor=TRUE)

# ---> Clonotype sharing between virus reactivities, experiment A.
tmp.cols <- virus.cols[c('FLU', 'CMV', 'EBV')]
depict.clone.sharing(meta.data=meta.data, tag.of.int='virus.tag', reports.path=fig.5.path, file.preffix='HltyDonor-ExpA_TCR-Repertoire_ClonotypeSharingBetweenReactivities', fill.cols=tmp.cols)
meta.data[!is.na(clonotype.tag), .N, by=.(virus.tag, clonotype.tag)][N>1][, uniqueN(clonotype.tag), by=virus.tag]

# ---> Overall cross-reactive TCR representation, experiment A.
# Obtain plots.
tmp.cols <- sharing.cols; names(tmp.cols) <- c('pR', 'non-pR')
tmp.ggplot.1 <- plot.props(seurat.obj=srt.objs.list[[main.obj.5]], plot.type='pie', groups.tag='pr.tag', na.rm=TRUE, groups.of.int=NULL, groups.label=NULL, filter.tags=NULL, groups.to.filter=NULL, keep=TRUE, color.vals=tmp.cols)
# Output.
publish.plot(tmp.ggplot=tmp.ggplot.1, output.path=fig.5.path, file.name='HltyDonor-ExpA_TCR-Repertoire_pR-RepresentationOverall', blank.comp=blank.complement.4, do.legend=TRUE, height=5, width=5)
# Report proportions per pR status:
meta.data[!is.na(clonotype.tag), .(cell.frac=.N/meta.data[!is.na(clonotype.tag), .N]), by=pr.tag]
# pR 0.1015945
# non-pR 0.8984055

# ---> Cross-reactive TCR representation per donor, experiment A.
# Get data.
tmp.data <- meta.data[
  !is.na(clonotype.tag) & !is.na(donor.id.tag),
  .(
    pr.prop=.SD[pr.tag=='pR', .N]/.N
  ),
  by=.(donor=donor.id.tag)
]
tmp.data[, donor:=paste0('D', str_extract(string=donor, pattern='\\d{2}$'))]
setorderv(x=tmp.data, cols='pr.prop', order=-1)
tmp.data$donor <- factor(x=as.character(tmp.data$donor), levels=as.character(tmp.data$donor))
# Plot and output.
tmp.ggplot <- ggplot(data=tmp.data, aes(x=donor, y=pr.prop)) +
  geom_bar(stat='identity', fill='lightgray', width=0.7) +
  scale_y_continuous(limits=c(0, 1), expand=expansion(add=0), breaks=scales::pretty_breaks(n=3)) + scale_x_discrete(expand=expansion(add=0.5)) +
  labs(x='Donors ranked according to crossreactivity representation', y='Percent of total cells') +
  theme(axis.ticks.x=element_blank())
publish.plot(tmp.ggplot=tmp.ggplot, output.path=fig.5.path, file.name='HltyDonor-ExpA_TCR-Repertoire_pR-RepresentationPerDonor', type='pdf', blank.comp=blank.complement.1, do.legend=TRUE, width=2)

# ---> pR profile on UMAP along clone size depiction, experiment A.
tmp.cells <- meta.data[!is.na(clonotype.tag), cell]
meta.data[, tmp.clon.size.tag:=clon.size.tag]; meta.data[clon.size.tag>=100, tmp.clon.size.tag:=100]
tmp.ggplot <- ggplot(data=meta.data[cell %in% tmp.cells], aes(x=UMAP_1, y=UMAP_2, fill=pr.tag, size=tmp.clon.size.tag)) +
  geom_point(stroke=0.4, shape=21, col='black') +
  scale_fill_manual(values=pr.groups.cols) +
  scale_radius(breaks=c(1, 10, 100), limits=c(1, meta.data[, max(tmp.clon.size.tag, na.rm=TRUE)]), range=c(3, 13)) +
  scale_y_continuous(expand=expansion(add=0)) + scale_x_continuous(expand=expansion(add=0)) +
  labs(x='UMAP 1', y='UMAP 2', col='')
# Output.
publish.plot(tmp.ggplot=tmp.ggplot, output.path=fig.5.path, file.name='HltyDonor-ExpA_SharingStatusOnUMAP', type='tiff', blank.comp=blank.complement.3, do.legend=TRUE)

# ---> Contrasting clonal expansion between crossreactive and non-crossreactive TCRs.
# Plot and output.
tmp.cols <- sharing.cols; names(tmp.cols) <- c('pR', 'non-pR')
srt.objs.list[[main.obj.5]]@meta.data[, 'log.clon.size.tag'] <- log2(srt.objs.list[[main.obj.5]]@meta.data[, 'clon.size.tag'])
tmp.ggplot <-
vln.plot(seurat.obj=srt.objs.list[[main.obj.5]], feature='log.clon.size.tag', groups.tag='pr.tag', na.rm=TRUE, groups.of.int=NULL, groups.label=NULL, filter.tags=NULL, groups.to.filter=NULL, keep=TRUE, feature.thold=1, color='median', vln.type='density', size.thold=0, file.name=NULL, adjust.val=2, trim.val=TRUE, plot.limits=NULL, add.points=FALSE, this.color.scale=tmp.cols)
tmp.ggplot <- tmp.ggplot + scale_x_continuous(expand=expansion(add=0), breaks=c(2, 6)) + scale_y_continuous(expand=expansion(add=0), breaks=scales::pretty_breaks(n=2))
publish.plot(tmp.ggplot=tmp.ggplot, output.path=fig.5.path, file.name='HltyDonor-ExpA_TCR-Repertoire_SizeBetween-pR-and-nonpR', type='pdf', blank.comp=blank.complement.1, do.legend=TRUE, width=4, height=4)

# ---> Gene signatures for experiment A
# Plot.
tmp.cells <- meta.data[, cell]
tmp.ggplot <- ggplot(data=meta.data[cell %in% tmp.cells], aes_string(x='UMAP_1', y='UMAP_2', col='cell.cytotoxicity.patil.score', fill='cell.cytotoxicity.patil.score')) +
  geom_point(size=4, alpha=0.8, shape=19) +
  scale_fill_gradientn(colors=signatures.col.scale) +
  scale_color_gradientn(colors=signatures.col.scale) +
  scale_y_continuous(expand=expansion(add=0)) + scale_x_continuous(expand=expansion(add=0)) +
  labs(x='UMAP 1', y='UMAP 2', col='Score')
# Output.
publish.plot(tmp.ggplot=tmp.ggplot, output.path=fig.5.path, file.name='HltyDonor-ExpA_Signature-Cytotoxicity', type='tiff', blank.comp=blank.complement.3, do.legend=TRUE)

# ---> FGSEA, cytotoxicity signature for experiment A.
do.fgsea(metrics.src=srt.objs.list[[main.obj.5]], tag.of.int='pr.tag', metric='signal.to.noise', output.path=fig.5.path, vals.to.depict='pR', files.preffix='HltyDonor-ExpA')

# ---> pR profile across CD4-CTLs and rest of the cells.
# Define CD4-CTL tag in seurat object.
tmp.corrs <- c('0'='Rest', '1'='CD4-CTL', '2'='Rest', '3'='Rest')
tmp.data <- as.character(srt.objs.list[[main.obj.5]]@meta.data[, clusts.lab[main.obj.5]])
tmp.data <- tmp.corrs[tmp.data]
srt.objs.list[[main.obj.5]]@meta.data[, 'ctl.tag'] <- tmp.data
# seurat.obj@meta.data[, 'ctl.tag'] <- tmp.data
tmp.cols <- c('Rest'='#808080', 'CD4-CTL'='#cc0000')
# Obtain plots.
tmp.ggplot.1 <- plot.props(seurat.obj=srt.objs.list[[main.obj.5]], plot.type='donut', groups.tag='ctl.tag', na.rm=TRUE, groups.of.int=NULL, groups.label=NULL, filter.tags='pr.tag', groups.to.filter='pR', keep=TRUE, color.vals=tmp.cols)
tmp.ggplot.2 <- plot.props(seurat.obj=srt.objs.list[[main.obj.5]], plot.type='donut', groups.tag='ctl.tag', na.rm=TRUE, groups.of.int=NULL, groups.label=NULL, filter.tags='pr.tag', groups.to.filter='non-pR', keep=TRUE, color.vals=tmp.cols)
# Output.
publish.plot(tmp.ggplot=tmp.ggplot.1, output.path=fig.5.path, file.name='HltyDonor-ExpA_CD4-CTLs-and-pRAssessment-pR', blank.comp=blank.complement.4, do.legend=TRUE, height=5, width=5)
publish.plot(tmp.ggplot=tmp.ggplot.2, output.path=fig.5.path, file.name='HltyDonor-ExpA_CD4-CTLs-and-pRAssessment-nonpR', blank.comp=blank.complement.4, do.legend=TRUE, height=5, width=5)


### -------------------- Supplementary Figure 5 --------------------- ###

# ----> Output directory.
sup.5.path <- paste0(reports.path, '/sup_figure_5')
if(!dir.exists(sup.5.path)) dir.create(sup.5.path)
# ---> Comments: Main dataset for this section, The same as in the main figure.

# ---> QC-related figures, experiment D.
# @ Distribution of the number of genes per cell across sequencing batches.
tmp.ggplot <- vln.plot(seurat.obj=srt.objs.list[[sec.obj.5]], feature='nFeature_RNA', groups.tag='library.tag', na.rm=TRUE, groups.of.int=NULL, groups.label=NULL, filter.tags=NULL, groups.to.filter=NULL, keep=TRUE, feature.thold=NULL, color='color', size.thold=0, file.name=NULL, adjust.val=1, trim.val=TRUE, plot.limits=NULL, add.points=FALSE, this.color.scale=NULL)
tmp.cols <- rep(x='lightblue', times=sec.meta.data[, uniqueN(library.tag)]); names(tmp.cols) <- sec.meta.data[, unique(library.tag)]
tmp.ggplot <- tmp.ggplot + scale_fill_manual(values=tmp.cols) + scale_y_continuous(limits=c(0, NA), breaks=c(1000, 2500, 4000))
publish.plot(tmp.ggplot=tmp.ggplot, output.path=sup.5.path, file.name='HltyDonor-ExpD_QCs_GenesPerCellDist', type='pdf', blank.comp=blank.complement.1, do.legend=TRUE, height=5, width=5)
# @ Fraction of cell cluster per sequencing batch.
# Get data
tmp.data <- get.tag.analysis(seurat.obj=srt.objs.list[[sec.obj.5]], tmp.meta.data=as.data.frame(meta.data), tmp.tag='library.tag', clusters.tag=clusts.lab[sec.obj.5], tag.reports.path=NA, reports.pref=NA, vals.cols=NULL, output.umaps=FALSE)
tmp.data <- as.data.frame(t(tmp.data), stringsAsFactors=TRUE); tmp.data$cluster <- row.names(tmp.data)
tmp.data <- tmp.data[tmp.data$cluster %in% names(clusts.cols[[sec.obj.5]]), ]
tmp.data <- gather(data=tmp.data, key='seq.batch', value='fraction', -`cluster`)
# Plot and output.
tmp.ggplot <- ggplot(data=tmp.data, aes(x=seq.batch, y=fraction, fill=cluster)) +
  geom_bar(stat='identity', position='fill', alpha=0.7) +
  scale_fill_manual(values=clusts.cols[[sec.obj.5]]) +
  scale_y_continuous(breaks=c(0, 0.5, 1)) + coord_flip() +
  labs(x='Sequencing batch', y='Fraction', fill='Cluster')
publish.plot(tmp.ggplot=tmp.ggplot, output.path=sup.5.path, file.name='HltyDonor-ExpD_QCs_ClusterFractsPerSeqBatch', type='pdf', blank.comp=blank.complement.1, do.legend=TRUE)

# ---> FGSEA, cytotoxicity signature, experiment D.
do.fgsea(metrics.src=srt.objs.list[[sec.obj.5]], tag.of.int='pr.tag', metric='signal.to.noise', output.path=sup.5.path, vals.to.depict='pR', files.preffix='HltyDonor-ExpD')

# ---> DEA, HCoV- vs SARS-CoV-2-reactive CD4+ T cells.
# Get data.
tmp.data <- hlty.six.expd.s3.dea[,
.(
    gene=gene.id, lfc=lfc, p.adj,
    log.p.adj=-log10(p.adj),
    mean.hcov=mean.virus.HCoV, mean.sars=`mean.virus.SARS-CoV-2`,
    prop.hcov=prop.virus.HCoV, prop.sars=`prop.virus.SARS-CoV-2`
  )
]
# Set color scale.
tmp.data[p.adj<=0.05 & lfc<=0.25, col.scale:=log2(mean.sars)]
tmp.data[p.adj<=0.05 & lfc>=0.25, col.scale:=log2(mean.hcov)]
tmp.data[col.scale<0, col.scale:=0]
# Set size scale.
tmp.data[p.adj<=0.05 & lfc<=0.25, size.scale:=prop.sars]
tmp.data[p.adj<=0.05 & lfc>=0.25, size.scale:=prop.hcov]
tmp.data[is.na(size.scale), size.scale:=0.01]
# Set boundaries in scales for better scale definition after outlier 'limitation'
tmp.data[, size.scale:=size.scale*100] # Color scale.
tmp.data[lfc>4, lfc:=4]
tmp.data[log.p.adj>300, log.p.adj:=300]
# Get selected DEGs to display in the paper.
top.genes <- c('IFNG', 'IL2', 'IL3', 'TBX21', 'CSF2', 'IL21', 'IRF1')
top.genes <- tmp.data[gene %in% top.genes]
# Get plot.
tmp.ggplot <- ggplot(data=tmp.data, aes(x=lfc, y=log.p.adj)) +
  geom_jitter(aes(col=col.scale, size=size.scale)) +
  # ggrepel::geom_text_repel(data=top.genes, aes(label=gene), color='black') +
  geom_vline(xintercept=-0.25, linetype='dashed') + geom_vline(xintercept=0.25, linetype='dashed') +
  geom_hline(yintercept=-log10(0.05), linetype='dashed') +
  scale_color_gradientn(colors=signatures.col.scale, na.value='#a6a6a6') +
  scale_radius(breaks=c(1, 25, 50, 75), range=c(1, 3)) +
  scale_x_continuous(limits=c(-2.1, 4.1), expand=expansion(add=0.2), breaks=c(-2, 0, 2, 4)) +
  scale_y_continuous(expand=expansion(add=5), breaks=c(0, 150, 300)) +
  labs(x='Log fold change', y='Log10 FDR', col='Avg. exp.', size='Prop. Exp.')
publish.plot(tmp.ggplot=tmp.ggplot, output.path=sup.5.path, file.name='HltyDonor-ExpD_DEA_HCoV-vs-SARS', type='pdf', blank.comp=blank.complement.1, repel.data=top.genes, repel.label='gene', do.legend=TRUE)
DONAS

# ---> QC-related figures, experiment A.
tmp.ggplot <- vln.plot(seurat.obj=srt.objs.list[[main.obj.5]], feature='nFeature_RNA', groups.tag='chrom.batch.tag', na.rm=TRUE, groups.of.int=NULL, groups.label=NULL, filter.tags=NULL, groups.to.filter=NULL, keep=TRUE, feature.thold=NULL, color='color', size.thold=0, file.name=NULL, adjust.val=1, trim.val=TRUE, plot.limits=NULL, add.points=FALSE, this.color.scale=NULL)
tmp.cols <- rep(x='lightblue', times=meta.data[, uniqueN(chrom.batch.tag)]); names(tmp.cols) <- meta.data[, unique(chrom.batch.tag)]
tmp.ggplot <- tmp.ggplot + scale_fill_manual(values=tmp.cols) + scale_y_continuous(limits=c(0, NA), breaks=c(1000, 2500, 4000))
publish.plot(tmp.ggplot=tmp.ggplot, output.path=sup.5.path, file.name='HltyDonor-ExpA_QCs_GenesPerCellDist', type='pdf', blank.comp=blank.complement.1, do.legend=TRUE, height=5, width=5)
# @ Fraction of cell cluster per sequencing batch.
# Get data
tmp.data <- get.tag.analysis(seurat.obj=srt.objs.list[[main.obj.5]], tmp.meta.data=as.data.frame(meta.data), tmp.tag='chrom.batch.tag', clusters.tag=clusts.lab[main.obj.5], tag.reports.path=NA, reports.pref=NA, vals.cols=NULL, output.umaps=FALSE)
tmp.data <- as.data.frame(t(tmp.data), stringsAsFactors=TRUE); tmp.data$cluster <- row.names(tmp.data)
tmp.data <- tmp.data[tmp.data$cluster %in% names(clusts.cols[[main.obj.5]]), ]
tmp.data <- gather(data=tmp.data, key='chrom.batch', value='fraction', -`cluster`)
# Plot and output.
tmp.ggplot <- ggplot(data=tmp.data, aes(x=chrom.batch, y=fraction, fill=cluster)) +
  geom_bar(stat='identity', position='fill', alpha=0.7) +
  scale_fill_manual(values=clusts.cols[[main.obj.5]]) +
  scale_y_continuous(breaks=c(0, 0.5, 1)) + coord_flip() +
  labs(x='Library', y='Fraction', fill='Cluster')
publish.plot(tmp.ggplot=tmp.ggplot, output.path=sup.5.path, file.name='HltyDonor-ExpA_QCs_ClusterFractsPerSeqBatch', type='pdf', blank.comp=blank.complement.1, do.legend=TRUE)

# ---> UMAP depicting clusters for experiment A.
# Plot.
tmp.ggplot <- ggplot(data=meta.data, aes_string(x='UMAP_1', y='UMAP_2', col=clusts.lab[main.obj.5])) +
  geom_point(size=4, alpha=0.4, shape=19) +
  scale_color_manual(values=clusts.cols[[main.obj.5]]) +
  scale_y_continuous(expand=expansion(add=0)) + scale_x_continuous(expand=expansion(add=0)) +
  labs(x='UMAP 1', y='UMAP 2', col='Cluster')
# Output.
publish.plot(tmp.ggplot=tmp.ggplot, output.path=sup.5.path, file.name='HltyDonor-ExpA_ClustersOnUMAP', type='tiff', blank.comp=blank.complement.3, do.legend=TRUE)

# ---> Dot plot.
# @ Version A.
# Define genes of interest.
these.markers <- c(
  'GZMB', 'PRF1', 'GNLY', 'NKG7', 'XCL1', 'XCL2' # CD4-CTL
)
# Plot.
tmp.ggplot <- dot.plot(
  seurat.obj=srt.objs.list[[main.obj.5]], features=these.markers, slot='data', do.norm=FALSE, ensembl=FALSE,
  groups.tag=clusts.lab[main.obj.5], groups.order=NULL, groups.of.int=c('1'), filter.tag=NULL, groups.to.filter=NULL, keep=TRUE, na.rm=TRUE, feature.thold=NULL,
  this.color.scale=signatures.col.scale, col.min=NULL, col.max=NULL,
  scale.by='radius', dot.scale=6, size.min=NA, size.max=NA,
  file.name=NULL
)
# Output.
publish.plot(tmp.ggplot=tmp.ggplot, output.path=sup.5.path, file.name='HltyDonor-ExpA_Markers_Opt-A', type='pdf', blank.comp=blank.complement.2, do.legend=TRUE, height=4, width=3)


### ------------- Other things that may fit in Figure 5 ------------- ###

# ---> pR profile on UMAP along clone size depiction, experiment D.
tmp.cells <- sec.meta.data[!is.na(clonotype.tag), cell]
sec.meta.data[, tmp.clon.size.tag:=clon.size.tag]; sec.meta.data[clon.size.tag>=100, tmp.clon.size.tag:=100]
tmp.ggplot <- ggplot(data=sec.meta.data[cell %in% tmp.cells], aes(x=UMAP_1, y=UMAP_2, fill=pr.tag, size=tmp.clon.size.tag)) +
  geom_point(stroke=0.4, shape=21, col='black') +
  scale_fill_manual(values=pr.groups.cols) +
  scale_radius(breaks=c(1, 3, 50, 100), limits=c(1, sec.meta.data[, max(tmp.clon.size.tag, na.rm=TRUE)]), range=c(3, 13)) +
  scale_y_continuous(expand=expansion(add=0)) + scale_x_continuous(expand=expansion(add=0)) +
  labs(x='UMAP 1', y='UMAP 2', col='')
# Output.
publish.plot(tmp.ggplot=tmp.ggplot, output.path=fig.5.path, file.name='HltyDonor-ExpD_SharingStatusOnUMAP', type='tiff', blank.comp=blank.complement.3, do.legend=TRUE)

# ---> Gene signatures for experiment D
# Plot.
tmp.cells <- sec.meta.data[, cell]
tmp.ggplot <- ggplot(data=sec.meta.data[cell %in% tmp.cells], aes_string(x='UMAP_1', y='UMAP_2', col='cell.cytotoxicity.patil.score', fill='cell.cytotoxicity.patil.score')) +
  geom_point(size=4, alpha=0.8, shape=19) +
  scale_fill_gradientn(colors=signatures.col.scale) +
  scale_color_gradientn(colors=signatures.col.scale) +
  scale_y_continuous(expand=expansion(add=0)) + scale_x_continuous(expand=expansion(add=0)) +
  labs(x='UMAP 1', y='UMAP 2', col='Score')
# Output.
publish.plot(tmp.ggplot=tmp.ggplot, output.path=fig.5.path, file.name='HltyDonor-ExpD_Signature-Cytotoxicity', type='tiff', blank.comp=blank.complement.3, do.legend=TRUE)

# ---> pR profile (fraction) for CD4-CTLs and TFHs, experiment D.
# Define CD4-CTL tag in seurat object.
tmp.corrs <- c('CTL', rep(x='Rest', times=2)); names(tmp.corrs) <- c('2', '0', '1')
tmp.data <- as.character(srt.objs.list[[sec.obj.5]]@meta.data[, clusts.lab[sec.obj.5]])
tmp.data <- tmp.corrs[tmp.data]
srt.objs.list[[sec.obj.5]]@meta.data[, 'ctl.tag'] <- tmp.data
tmp.cols <- sharing.cols; names(tmp.cols) <- c('pR', 'non-pR')
# Obtain plots.
tmp.ggplot.1 <- plot.props(seurat.obj=srt.objs.list[[sec.obj.5]], plot.type='donut', groups.tag='pr.tag', na.rm=TRUE, groups.of.int=NULL, groups.label=NULL, filter.tags='ctl.tag', groups.to.filter='CTL', keep=TRUE, color.vals=tmp.cols)
tmp.ggplot.2 <- plot.props(seurat.obj=srt.objs.list[[sec.obj.5]], plot.type='donut', groups.tag='pr.tag', na.rm=TRUE, groups.of.int=NULL, groups.label=NULL, filter.tags='ctl.tag', groups.to.filter='Rest', keep=TRUE, color.vals=tmp.cols)
# Output.
publish.plot(tmp.ggplot=tmp.ggplot.1, output.path=fig.5.path, file.name='HltyDonor-ExpD_pRFractions-CD4-CTL', blank.comp=blank.complement.4, do.legend=TRUE, height=5, width=5)
publish.plot(tmp.ggplot=tmp.ggplot.2, output.path=fig.5.path, file.name='HltyDonor-ExpD_pRFractions-Non-CTL', blank.comp=blank.complement.4, do.legend=TRUE, height=5, width=5)

# ---> pR profile across CD4-CTLs and rest of the cells.
# Define CD4-CTL tag in seurat object.
tmp.corrs <- c('0'='Rest', '1'='Rest', '2'='CD4-CTL')
tmp.data <- as.character(srt.objs.list[[sec.obj.5]]@meta.data[, clusts.lab[sec.obj.5]])
tmp.data <- tmp.corrs[tmp.data]
srt.objs.list[[sec.obj.5]]@meta.data[, 'ctl.tag'] <- tmp.data
# seurat.obj@meta.data[, 'ctl.tag'] <- tmp.data
tmp.cols <- c('Rest'='#808080', 'CD4-CTL'='#cc0000')
# Obtain plots.
tmp.ggplot.1 <- plot.props(seurat.obj=srt.objs.list[[sec.obj.5]], plot.type='donut', groups.tag='ctl.tag', na.rm=TRUE, groups.of.int=NULL, groups.label=NULL, filter.tags='pr.tag', groups.to.filter='pR', keep=TRUE, color.vals=tmp.cols)
tmp.ggplot.2 <- plot.props(seurat.obj=srt.objs.list[[sec.obj.5]], plot.type='donut', groups.tag='ctl.tag', na.rm=TRUE, groups.of.int=NULL, groups.label=NULL, filter.tags='pr.tag', groups.to.filter='non-pR', keep=TRUE, color.vals=tmp.cols)
# Output.
publish.plot(tmp.ggplot=tmp.ggplot.1, output.path=fig.5.path, file.name='HltyDonor-ExpD_CD4-CTLs-and-pRAssessment-pR', blank.comp=blank.complement.4, do.legend=TRUE, height=5, width=5)
publish.plot(tmp.ggplot=tmp.ggplot.2, output.path=fig.5.path, file.name='HltyDonor-ExpD_CD4-CTLs-and-pRAssessment-nonpR', blank.comp=blank.complement.4, do.legend=TRUE, height=5, width=5)

# ---> UMAP depicting clusters for experiment D.
# Plot.
tmp.ggplot <- ggplot(data=sec.meta.data, aes_string(x='UMAP_1', y='UMAP_2', col=clusts.lab[sec.obj.5])) +
  geom_point(size=1.2, alpha=0.6) +
  scale_color_manual(values=clusts.cols[[sec.obj.5]]) +
  scale_y_continuous(expand=expansion(add=0.3)) + scale_x_continuous(expand=expansion(add=0.3)) +
  labs(x='UMAP 1', y='UMAP 2', col='Cluster')
# Output.
publish.plot(tmp.ggplot=tmp.ggplot, output.path=fig.5.path, file.name='HltyDonor-ExpD_ClustersOnUMAP', type='tiff', blank.comp=blank.complement.3, do.legend=TRUE)

# ---> Virus-reactive profile of CD4 T cells in healthy donors.
# Info of virus-reactive profile per donor assessed.
tmp.data.1 <- hlty.cyt.prof[, .(donor.id, CMV=mem.cd154.pop.cmv.prop*100, EBV=mem.cd154.pop.ebv.prop*100)]
setorderv(x=tmp.data.1, cols=c('CMV', 'EBV'), order=-1)
tmp.data.1[, rank:=1:.N]
tmp.data.1 <- gather(data=tmp.data.1, -`donor.id`, -`rank`, key='reactivity', value='cell.prop')
# Info of donors that were picked for further analyses after assessment.
tmp.data.2 <- data.table(donor.id=paste0('SDBB-', c('088', '085', '020', '039', '065', '070')))
tmp.data.2 <- unique(merge(x=tmp.data.2, y=tmp.data.1, by='donor.id')[, .(donor.id, rank)])
# Plot.
tmp.cols <- virus.cols[c('CMV', 'EBV')]
tmp.ggplot <- ggplot(data=tmp.data.1, aes(x=rank, y=cell.prop)) +
  geom_point(aes(color=reactivity), size=3) + geom_line(aes(group=reactivity)) +
  scale_color_manual(values=tmp.cols) +
  scale_y_continuous(expand=expansion(add=c(0.01, 0.1))) + scale_x_continuous(expand=expansion(add=c(0.2, 0.2))) +
  labs(x='Donor Rank', y='CD154+/CD69+ CD4 T cells %', col='Peptide\npool')
for(tmp.rank in tmp.data.2[, rank]) tmp.ggplot <- tmp.ggplot + geom_vline(xintercept=tmp.rank, size=3, alpha=0.3, col='#e6e600')
# Output.
publish.plot(tmp.ggplot=tmp.ggplot, output.path=fig.5.path, file.name='HltyDonor-ExpA_VirusReactivityProfile', blank.comp=blank.complement.1, do.legend=TRUE, width=10)

# ---> CD4-CTL profile in healthy donors.
# Info of CD4-CTL profile per donor assessed.
tmp.data.1 <- hlty.cyt.prof[, .(donor.id, CMV=mem.cd154.pop.cmv.prop*100, `CMV-reactive`=cd4.ctl.cmv.prop*100, `Total`=cd4.ctl.total.prop*100)]
tmp.data.1 <- gather(data=tmp.data.1, -`donor.id`, -`CMV`, key='group', value='cell.prop')
# Info of donors that were picked for further analyses after assessment.
tmp.data.2 <- data.table(donor.id=paste0('SDBB-', c('088', '085', '020', '039', '065', '070')), picked='Yes')
tmp.data.1 <- merge(x=tmp.data.2, y=tmp.data.1, by='donor.id', all=TRUE)
tmp.data.1[is.na(picked), picked:='No']
# Plot.
tmp.cols <- c(`CMV-reactive`=unname(virus.cols[c('CMV')]), `Total`='#000000')
tmp.sizes <- c(Yes=2.5, No=0.8)
tmp.ggplot <- ggplot(data=tmp.data.1, aes(x=CMV, y=cell.prop, color=group)) +
  # geom_point(aes(shape=picked), size=1.5, alpha=0.8) + geom_smooth(method='lm', se=FALSE, size=0.6, alpha=0.8) +
  geom_point(aes(size=picked), alpha=0.8) + geom_smooth(method='lm', se=FALSE, size=0.6, alpha=0.8) +
  # stat_cor() +
  scale_color_manual(values=tmp.cols) +
  scale_size_manual(values=tmp.sizes) +
  scale_y_continuous(expand=expansion(add=c(0.01, 1))) + scale_x_continuous(expand=expansion(add=c(0.01, 0.02))) +
  labs(x='CD154+/CD69+ CD4 T cells %', y='CD4-CTLs %', col='', shape='Chosen')
# Output.
publish.plot(tmp.ggplot=tmp.ggplot, output.path=fig.5.path, file.name='HltyDonor-ExpA_CD4-CTLProfile', blank.comp=blank.complement.1, do.legend=TRUE, width=10, stat.cor=TRUE)

# ---> pR profile (fraction) for CD4-CTLs and TFHs, experiment A.
# Define CD4-CTL tag in seurat object.
# tmp.corrs <- c(rep(x='CTL', times=2), rep(x='Rest', times=2)); names(tmp.corrs) <- c('1', '3', '0', '2')
tmp.corrs <- c('CTL', rep(x='Rest', times=3)); names(tmp.corrs) <- c('1', '3', '0', '2')
tmp.data <- as.character(srt.objs.list[[main.obj.5]]@meta.data[, clusts.lab[main.obj.5]])
tmp.data <- tmp.corrs[tmp.data]
srt.objs.list[[main.obj.5]]@meta.data[, 'ctl.tag'] <- tmp.data
tmp.cols <- sharing.cols; names(tmp.cols) <- c('pR', 'non-pR')
# Obtain plots.
tmp.ggplot.1 <- plot.props(seurat.obj=srt.objs.list[[main.obj.5]], plot.type='donut', groups.tag='pr.tag', na.rm=TRUE, groups.of.int=NULL, groups.label=NULL, filter.tags='ctl.tag', groups.to.filter='CTL', keep=TRUE, color.vals=tmp.cols)
tmp.ggplot.2 <- plot.props(seurat.obj=srt.objs.list[[main.obj.5]], plot.type='donut', groups.tag='pr.tag', na.rm=TRUE, groups.of.int=NULL, groups.label=NULL, filter.tags='ctl.tag', groups.to.filter='Rest', keep=TRUE, color.vals=tmp.cols)
# Output.
publish.plot(tmp.ggplot=tmp.ggplot.1, output.path=fig.5.path, file.name='HltyDonor-ExpA_pRFractions-CD4-CTL', blank.comp=blank.complement.4, do.legend=TRUE, height=5, width=5)
publish.plot(tmp.ggplot=tmp.ggplot.2, output.path=fig.5.path, file.name='HltyDonor-ExpA_pRFractions-Non-CTL', blank.comp=blank.complement.4, do.legend=TRUE, height=5, width=5)


############    -----------------------------------------    ############
### --------------------------- Figure 6 ---------------------------- ###
############    -----------------------------------------    ############

# ---> General definitions.
# Seurat objects of interest.
main.obj.6 <- 'HLTY-6-ExpB-SFC'
sec.obj.6 <- 'HLTY-6-ExpC-ALL'
# Seurat objects' meta data
meta.data <- cbind(
  data.table(cell=Cells(srt.objs.list[[main.obj.6]])),
  srt.objs.list[[main.obj.6]]@meta.data,
  srt.objs.list[[main.obj.6]]@reductions$umap@cell.embeddings
)
sec.meta.data <- cbind(
  data.table(cell=Cells(srt.objs.list[[sec.obj.6]])),
  srt.objs.list[[sec.obj.6]]@meta.data,
  srt.objs.list[[sec.obj.6]]@reductions$umap@cell.embeddings
)

### ------------------------- Main Figure 6 ------------------------- ###

# ----> Output directory.
fig.6.path <- paste0(reports.path, '/figure_6')
if(!dir.exists(fig.6.path)) dir.create(fig.6.path)
# ---> Comments: Main dataset for this section, The same as in the main figure.

# ---> Titration experiments.
donors.of.int <- titration.exp[, unique(donor.id)]
pops.of.int <- c('pmpbmc.cd154.cd69.dp', 'cd154.cd69.slamf7.tp')
for(tmp.donor in donors.of.int){
  for(tmp.pop in pops.of.int){
    # Retrieve data for donor and population of interest.
    tmp.data <- titration.exp[donor.id==tmp.donor & cell.population==tmp.pop]
    tmp.cols <- virus.cols[names(virus.cols) %in% tmp.data[, unique(peptide.pool)]]
    # Get plot.
    tmp.ggplot <- ggplot(data=tmp.data, aes(x=concentration, y=count, col=peptide.pool)) +
      geom_point(size=7) + geom_line(aes(group=peptide.pool), size=2.5) +
      scale_color_manual(values=tmp.cols) +
      scale_y_continuous(breaks=scales::pretty_breaks(n=2)) + scale_x_discrete(expand=expansion(add=0.1)) +
      labs(x='Concentration (ug/mL)', y=ifelse(test=tmp.pop=='pmpbmc.cd154.cd69.dp', yes='Antigen-reactive cells per 10e6 PBMC', no='SLAMF7+ antigen-reactive cells'))
    tmp.file.name <- paste0('HltyDonor-ExpB_TitrationResults_Donor-', tmp.donor, '_Pop-', ifelse(test=tmp.pop=='pmpbmc.cd154.cd69.dp', yes='General', no='SLAMF7'))
    publish.plot(tmp.ggplot=tmp.ggplot, output.path=fig.6.path, file.name=tmp.file.name, type='pdf', blank.comp=blank.complement.1, do.legend=FALSE, height=4, width=4)
  }
}

# ---> Clonotype sharing among virus reactivities, experiment B.
tmp.cols <- virus.cols[c('FLU', 'CMV', 'SARS-CoV-2')]
depict.clone.sharing(meta.data=meta.data, tag.of.int='virus.tag', reports.path=fig.6.path, file.preffix='HltyDonor-ExpB_TCR-Repertoire_ClonotypeSharingBetweenReactivities', fill.cols=tmp.cols)
meta.data[!is.na(clonotype.tag), .N, by=.(virus.tag, clonotype.tag)][N>1][, uniqueN(clonotype.tag), by=virus.tag]

# ---> pR profile on UMAP along clone size depiction, experiment B.
tmp.cells <- meta.data[!is.na(clonotype.tag), cell]
meta.data[, tmp.clon.size.tag:=clon.size.tag]; meta.data[clon.size.tag>=1500, tmp.clon.size.tag:=1500]
tmp.ggplot <- ggplot(data=meta.data[cell %in% tmp.cells], aes(x=UMAP_1, y=UMAP_2, fill=pr.tag, size=tmp.clon.size.tag)) +
  geom_point(stroke=0.4, shape=21, col='black') +
  scale_fill_manual(values=pr.groups.cols) +
  scale_radius(breaks=c(1, 7, 200, 1000), limits=c(1, meta.data[, max(tmp.clon.size.tag, na.rm=TRUE)]), range=c(3, 13)) +
  scale_y_continuous(expand=expansion(add=0)) + scale_x_continuous(expand=expansion(add=0)) +
  labs(x='UMAP 1', y='UMAP 2', col='')
# Output.
publish.plot(tmp.ggplot=tmp.ggplot, output.path=fig.6.path, file.name='HltyDonor-ExpB_SharingStatusOnUMAP', type='tiff', blank.comp=blank.complement.3, do.legend=TRUE)

# ---> Contrasting clonal expansion between crossreactive and non-crossreactive TCRs, for experiment B.
# Plot and output.
tmp.cols <- sharing.cols; names(tmp.cols) <- c('pR', 'non-pR')
srt.objs.list[[main.obj.6]]@meta.data[, 'log.clon.size.tag'] <- log2(srt.objs.list[[main.obj.6]]@meta.data[, 'clon.size.tag'])
tmp.ggplot <-
vln.plot(seurat.obj=srt.objs.list[[main.obj.6]], feature='log.clon.size.tag', groups.tag='pr.tag', na.rm=TRUE, groups.of.int=NULL, groups.label=NULL, filter.tags=NULL, groups.to.filter=NULL, keep=TRUE, feature.thold=1, color='median', vln.type='density', size.thold=0, file.name=NULL, adjust.val=2, trim.val=TRUE, plot.limits=NULL, add.points=FALSE, this.color.scale=tmp.cols)
tmp.ggplot <- tmp.ggplot + scale_x_continuous(expand=expansion(add=0), breaks=scales::pretty_breaks(n=2)) + scale_y_continuous(expand=expansion(add=0), breaks=scales::pretty_breaks(n=2))#, limits=c(0, 0.4))
publish.plot(tmp.ggplot=tmp.ggplot, output.path=fig.6.path, file.name='HltyDonor-ExpB_TCR-Repertoire_SizeBetween-pR-and-nonpR', type='pdf', blank.comp=blank.complement.1, do.legend=TRUE, width=4, height=4)

# ---> Gene signatures for experiment B.
# Plot.
tmp.cells <- meta.data[, cell]
tmp.ggplot <- ggplot(data=meta.data[cell %in% tmp.cells], aes_string(x='UMAP_1', y='UMAP_2', col='cell.cytotoxicity.patil.score', fill='cell.cytotoxicity.patil.score')) +
  geom_point(size=4, alpha=0.8, shape=19) +
  scale_fill_gradientn(colors=signatures.col.scale) +
  scale_color_gradientn(colors=signatures.col.scale) +
  scale_y_continuous(expand=expansion(add=0)) + scale_x_continuous(expand=expansion(add=0)) +
  labs(x='UMAP 1', y='UMAP 2', col='Score')
# Output.
publish.plot(tmp.ggplot=tmp.ggplot, output.path=fig.6.path, file.name='HltyDonor-ExpB_Signature-Cytotoxicity', type='tiff', blank.comp=blank.complement.3, do.legend=TRUE)

# ---> FGSEA, cytotoxicity signature, experiment B
do.fgsea(metrics.src=srt.objs.list[[main.obj.6]], tag.of.int='pr.tag', metric='signal.to.noise', output.path=fig.6.path, vals.to.depict='pR', files.preffix='HltyDonor-ExpB')

# ---> pR profile across virus reativities and ARTE peptide pool concentrations.
# Get data.
tmp.data <- meta.data[!is.na(clonotype.tag), .(cell.count=.N), by=.(clonotype=clonotype.tag, reactivity.and.arte=paste(virus.tag, arte.concentration.tag, sep='-'))]
tmp.data <- spread(data=tmp.data, key=reactivity.and.arte, value=cell.count, fill=NA)
# Check 90% percentile per virus reactivity and ARTE peptide pool concentration. These will the thresholds to call a given clonotype whether as expanded or not.
size.tholds <- apply(X=tmp.data[, -'clonotype'], MARGIN=2, FUN=quantile, probs=0.9, na.rm=TRUE)
# Replace sizes with expansion values.
tmp.data <- lapply(X=names(size.tholds), FUN=function(tmp.val){
  # Retrieve threshold.
  tmp.thold <- size.tholds[tmp.val]
  # Set threshold.
  tmp.data <- tmp.data[!is.na(clonotype)]
  to.return <- as.integer(tmp.data[, get(tmp.val)>=tmp.thold])
  to.return <- data.table(clonotype=tmp.data[, clonotype], val=to.return)
  return(to.return)
})
names(tmp.data) <- names(size.tholds); tmp.data <- rbindlist(l=tmp.data, use.names=TRUE, idcol='reactivity.and.arte')
tmp.data <- spread(data=tmp.data, key=reactivity.and.arte, value=val)
# Keep only the clonotypes that are seen in over three combinations.
tmp.data <- tmp.data[rowSums(x=tmp.data[, -'clonotype'], na.rm=TRUE)>3, ]
# Identify bonafie (high avidity) clonotypes (those seen in all four titrations for at least two viruses).
query.viruses <- c('CMV', 'FLU', 'SARS-CoV-2')
bonafide.cr.tcrs <- lapply(X=query.viruses, FUN=function(tmp.virus){
  query.cols <- paste(tmp.virus, c('0.001', '0.01', '0.1', '1'), sep='-')
  to.return <- rowSums(tmp.data[, ..query.cols], na.rm=TRUE)==4
  to.return <- data.table(clonotype=tmp.data[, clonotype], bonafide=to.return)
  return(to.return)
})
names(bonafide.cr.tcrs) <- query.viruses
bonafide.cr.tcrs <- rbindlist(l=bonafide.cr.tcrs, use.names=TRUE, idcol='virus')
bonafide.cr.tcrs <- spread(data=bonafide.cr.tcrs, key=virus, value=bonafide)
bonafide.cr.tcrs$bonafide <- rowSums(bonafide.cr.tcrs[, -'clonotype'])>1
# Further metadata for bonafide clonotypes.
bonafide.cr.tcrs <- merge(
  x=bonafide.cr.tcrs,
  y=meta.data[,
    .(
      size=unique(clon.size.tag),
      trb.aa=unique(TRB.aa.chains.tag),
      tra.aa=unique(TRA.aa.chains.tag)
    ), by=.(clonotype=clonotype.tag)],
  by='clonotype',
  all.x=TRUE, all.y=FALSE
)
# Set rows' order.
tmp.data <- tmp.data[, -'clonotype']
for(tmp.col in rev(colnames(tmp.data))){
  setorderv(x=tmp.data, cols=tmp.col, order=-1, na.last=TRUE)
}
row.order <- rowSums(x=tmp.data[, -'clonotype'], na.rm=TRUE); names(row.order) <- as.character(1:(length(row.order)))
row.order <- as.integer(names(sort(row.order, decreasing=TRUE)))
tmp.data <- tmp.data[row.order, ]
# Set metadata.
col.metadata <- data.frame(
  row.names=colnames(tmp.data),
  virus.reactivity=str_replace(string=colnames(tmp.data), pattern='-[^-]+$', replacement=''),
  arte.concentration=str_replace(string=colnames(tmp.data), pattern='^.+-', replacement='')
)
# Define colors for metadata tracks.
ann.colors <- list(
  virus.reactivity=virus.cols[unique(as.character(col.metadata$virus.reactivity))],
  arte.concentration=concentration.cols
)
# Set color scale and breaks for heatmap.
col.breaks <- seq(from=min(range(tmp.data, na.rm=T)), to=max(range(tmp.data, na.rm=T)), length.out=3)
hmap.col.scale <- colorRampPalette(c('blue', 'yellow'))(2)
# Plot 'complete' version
tmp.file.name <- paste0(fig.6.path, '/HltyDonor_TCR-RepertoireAcrossReactivitiesAndConcentrations.C.pdf')
pheatmap(mat=tmp.data, color=hmap.col.scale, breaks=col.breaks, scale='none', cluster_rows=FALSE, cluster_cols=FALSE, annotation_col=col.metadata, annotation_colors=ann.colors, show_colnames=FALSE, na_col='gray', filename=tmp.file.name)
# Plot 'blank' version
tmp.file.name <- paste0(fig.6.path, '/HltyDonor_TCR-RepertoireAcrossReactivitiesAndConcentrations.B.pdf')
pheatmap(mat=tmp.data, color=hmap.col.scale, breaks=col.breaks, scale='none', cluster_rows=FALSE, cluster_cols=FALSE, annotation_col=col.metadata, annotation_colors=ann.colors,
legend=FALSE, annotation_legend=FALSE, annotation_names_row=FALSE, annotation_names_col=FALSE, show_rownames=FALSE, show_colnames=FALSE, filename=tmp.file.name)

# ---> UMAP depicting SLAMF7-related groups for experiment C.
# Plot.
tmp.cells <- sec.meta.data[!is.na(clonotype.tag), cell]
sec.meta.data[, tmp.clon.size.tag:=clon.size.tag]; sec.meta.data[clon.size.tag>=500, tmp.clon.size.tag:=500]
tmp.ggplot <- ggplot(data=sec.meta.data[cell %in% tmp.cells], aes(x=UMAP_1, y=UMAP_2, fill=slamf7.pop.tag, size=tmp.clon.size.tag)) +
  geom_point(stroke=0.4, shape=21, col='black') +
  scale_fill_manual(values=slamf7.cols) +
  scale_radius(breaks=c(1, 7, 200, 500), limits=c(1, sec.meta.data[, max(tmp.clon.size.tag, na.rm=TRUE)]), range=c(3, 13)) +
  scale_y_continuous(expand=expansion(add=0.3), limits=c(-3, 4.6)) + scale_x_continuous(expand=expansion(add=0.3)) +
  labs(x='UMAP 1', y='UMAP 2', col='Population')
# Output.
publish.plot(tmp.ggplot=tmp.ggplot, output.path=fig.6.path, file.name='HltyDonor-ExpC_CellPopsOnUMAP', type='tiff', blank.comp=blank.complement.3, do.legend=TRUE)

# ---> Contrasting clonal expansion between SLAMF7- and SLAMF7+ T cells.
# Plot and output.
srt.objs.list[[sec.obj.6]]@meta.data[, 'log.clon.size.tag'] <- log2(srt.objs.list[[sec.obj.6]]@meta.data[, 'clon.size.tag'])
tmp.ggplot <-
vln.plot(seurat.obj=srt.objs.list[[sec.obj.6]], feature='log.clon.size.tag', groups.tag='slamf7.pop.tag', na.rm=TRUE, groups.of.int=NULL, groups.label=NULL, filter.tags=NULL, groups.to.filter=NULL, keep=TRUE, feature.thold=1, color='median', vln.type='density', size.thold=0, file.name=NULL, adjust.val=2, trim.val=TRUE, plot.limits=NULL, add.points=FALSE, this.color.scale=slamf7.cols)
tmp.ggplot <- tmp.ggplot + scale_x_continuous(expand=expansion(add=0), breaks=scales::pretty_breaks(n=3)) + scale_y_continuous(expand=expansion(add=0), breaks=scales::pretty_breaks(n=2))
publish.plot(tmp.ggplot=tmp.ggplot, output.path=fig.6.path, file.name='HltyDonor-ExpC_TCR-Repertoire_SizeBetween-SLAMF7Pos-and-Neg', type='pdf', blank.comp=blank.complement.1, do.legend=TRUE)

# ---> Gene signatures for experiment C
# Plot.
tmp.cells <- sec.meta.data[, cell]
tmp.ggplot <- ggplot(data=sec.meta.data[cell %in% tmp.cells], aes_string(x='UMAP_1', y='UMAP_2', col='cell.cytotoxicity.patil.score', fill='cell.cytotoxicity.patil.score')) +
  geom_point(size=4, alpha=0.8, shape=19) +
  scale_fill_gradientn(colors=signatures.col.scale) +
  scale_color_gradientn(colors=signatures.col.scale) +
  scale_y_continuous(expand=expansion(add=0), limits=c(-3, 4.6)) + scale_x_continuous(expand=expansion(add=0)) +
  labs(x='UMAP 1', y='UMAP 2', col='Score')
# Output.
publish.plot(tmp.ggplot=tmp.ggplot, output.path=fig.6.path, file.name='HltyDonor-ExpC_Signature-Cytotoxicity', type='tiff', blank.comp=blank.complement.3, do.legend=TRUE)

# ---> FGSEA, cytotoxicity signature, experiment C
do.fgsea(metrics.src=srt.objs.list[[sec.obj.6]], tag.of.int='slamf7.pop.tag', metric='signal.to.noise', output.path=fig.6.path, vals.to.depict='Positive', files.preffix='HltyDonor-ExpC')

# ---> High avidity cross-reactive clonotypes (beta chain level) are detected in the SLAMF7+ population and not in the SLAMF7- population.
# @ General code as a function to apply for either population.
comp.exps.fig6 <- function(pop.of.int, plot.cols){
  # Define clonotypes at the beta chain level for the population of interest. For clonotypes defined with multiple beta chains, we split their clone size according to the total number of beta chains (we see only a few with only two chains).
  tmp.data.1 <- sec.meta.data[
    !is.na(TRB.aa.chains.tag) & cell.pop.tag==pop.of.int,
    .(size=.N, chain.count=unique(str_count(string=TRB.aa.chains.tag, pattern=';')+1)),
    by=.(trb.aa=TRB.aa.chains.tag)
  ]
  tmp.data.1 <- as.data.table(separate_rows(data=tmp.data.1, trb.aa, sep=';'))
  tmp.data.1[, size:=size/chain.count]; tmp.data.1[, chain.count:=NULL]
  # Apply the clonotype calling in a similar fashion for the bonafide cross-reactive clonotypes.
  tmp.data.2 <- bonafide.cr.tcrs[bonafide==TRUE, .(trb.aa, size, chain.count=unique(str_count(string=trb.aa, pattern=';')+1))]
  tmp.data.2 <- as.data.table(separate_rows(data=tmp.data.2, trb.aa, sep=';'))
  tmp.data.2[, size:=size/chain.count]; tmp.data.2[, chain.count:=NULL]
  # Merge data from both datasets according to the TCR beta chains
  tmp.data <- merge(x=tmp.data.1, y=tmp.data.2, by='trb.aa', suffixes=c('.slamf7', '.bonafide'), all.x=TRUE, all.y=FALSE)
  tmp.data[, log.size.slamf7:=log2(size.slamf7+1)]
  # Rescale bonafide clonotype size (i.e., size in experiment C) to account for the maximum value in the y scale using next formula:
  #   New bonafide scale=\frac{max(count of geom_histogram -i.e., for the size provided through log.size.slamf7 with binwidth 1-) * previous bonafide scale}{max(previous bonafide scale)}
  #     For that, first calculate the maximum count from geom_histogram and the maximum value from log.size.slamf7
  tmp.data[, hist.val:=ceiling(log.size.slamf7)]
  max.hist <- tmp.data[, .N, by=hist.val][, max(N)]
  max.prev <- tmp.data[, max(size.bonafide, na.rm=TRUE)]
  #     Then, apply formula.
  tmp.data[, trans.size.bonafide:=(max.hist*size.bonafide)/max.prev]
  # Plot and output.
  tmp.ggplot <- ggplot(data=tmp.data, aes(x=log.size.slamf7)) +
    geom_histogram(binwidth=1, alpha=0.5, fill=plot.cols[1]) +
    geom_point(aes(y=trans.size.bonafide), size=7, col=plot.cols[2]) +
    scale_y_continuous(
      name='Histogram count', expand=expansion(add=5),
      sec.axis=sec_axis(trans=~.*1, name='Clone size after enrichment', labels=function(new.size) round(x=(max.prev*new.size)/max.hist, digits=2))
    ) +
    scale_x_continuous(breaks=1:11, labels=function(log.size){(2**log.size)-1}, expand=expansion(add=c(0, 0))) +
    labs(x='Clone size previous to enrichment (log2)')
  publish.plot(tmp.ggplot=tmp.ggplot, output.path=fig.6.path, file.name=paste0('HltyDonor_TCR-Repertoire_ComparisonOfExps_Pop-', ifelse(test=pop.of.int=='SLAMF7+', yes='SLAMF7Pos', no='SLAMF7Neg')), type='pdf', blank.comp=blank.complement.1, do.legend=FALSE, width=10)
}
# @ Do comparison taking into account the SLAMF7+ cells only.
comp.exps.fig6(pop.of.int='SLAMF7+', plot.cols=c('indianred', 'darkred'))
# @ Do comparison taking into account the SLAMF7- cells only.
comp.exps.fig6(pop.of.int='SLAMF7-/CCR7+', plot.cols=c('lightblue', 'darblue'))


### -------------------- Supplementary Figure 6 --------------------- ###

# ----> Output directory.
sup.6.path <- paste0(reports.path, '/sup_figure_6')
if(!dir.exists(sup.6.path)) dir.create(sup.6.path)
# ---> Comments: Main dataset for this section, The same as in the main figure.

# ---> UMAP depicting clusters for experiment B.
# Plot.
tmp.ggplot <- ggplot(data=meta.data, aes_string(x='UMAP_1', y='UMAP_2', col=clusts.lab[main.obj.6])) +
  geom_point(size=4, alpha=0.4, shape=19) +
  scale_color_manual(values=clusts.cols[[main.obj.6]]) +
  scale_y_continuous(expand=expansion(add=0.3)) + scale_x_continuous(expand=expansion(add=0.3)) +
  labs(x='UMAP 1', y='UMAP 2', col='Cluster')
# Output.
publish.plot(tmp.ggplot=tmp.ggplot, output.path=sup.6.path, file.name='HltyDonor-ExpB_ClustersOnUMAP', type='tiff', blank.comp=blank.complement.3, do.legend=TRUE)

# ---> QC-related figures, experiment B.
# @ Distribution of the number of genes per cell across sequencing batches.
tmp.ggplot <- vln.plot(seurat.obj=srt.objs.list[[main.obj.6]], feature='nFeature_RNA', groups.tag='chrom.batch.tag', na.rm=TRUE, groups.of.int=NULL, groups.label=NULL, filter.tags=NULL, groups.to.filter=NULL, keep=TRUE, feature.thold=NULL, color='color', size.thold=0, file.name=NULL, adjust.val=1, trim.val=TRUE, plot.limits=NULL, add.points=FALSE, this.color.scale=NULL)
tmp.cols <- rep(x='lightblue', times=meta.data[, uniqueN(chrom.batch.tag)]); names(tmp.cols) <- meta.data[, unique(chrom.batch.tag)]
tmp.ggplot <- tmp.ggplot + scale_fill_manual(values=tmp.cols) + scale_y_continuous(breaks=c(1000, 2500, 4000))
publish.plot(tmp.ggplot=tmp.ggplot, output.path=sup.6.path, file.name='HltyDonor-ExpB_QCs_GenesPerCellDist', type='pdf', blank.comp=blank.complement.1, do.legend=TRUE, height=4, width=5)
# @ Fraction of cell cluster per sequencing batch.
# Get data
tmp.data <- get.tag.analysis(seurat.obj=srt.objs.list[[main.obj.6]], tmp.meta.data=as.data.frame(meta.data), tmp.tag='chrom.batch.tag', clusters.tag=clusts.lab[main.obj.6], tag.reports.path=NA, reports.pref=NA, vals.cols=NULL, output.umaps=FALSE)
tmp.data <- as.data.frame(t(tmp.data), stringsAsFactors=TRUE); tmp.data$cluster <- row.names(tmp.data)
tmp.data <- tmp.data[tmp.data$cluster %in% names(clusts.cols[[main.obj.6]]), ]
tmp.data <- gather(data=tmp.data, key='library.batch', value='fraction', -`cluster`)
# Plot and output.
tmp.ggplot <- ggplot(data=tmp.data, aes(x=library.batch, y=fraction, fill=cluster)) +
  geom_bar(stat='identity', position='fill', alpha=0.7) +
  scale_fill_manual(values=clusts.cols[[main.obj.6]]) +
  scale_y_continuous(breaks=c(0, 0.5, 1)) + coord_flip() +
  labs(x='Sequencing batch', y='Fraction', fill='Cluster')
publish.plot(tmp.ggplot=tmp.ggplot, output.path=sup.6.path, file.name='HltyDonor-ExpB_QCs_ClusterFractsPerSeqBatch', type='pdf', blank.comp=blank.complement.1, do.legend=TRUE)

# ---> Dot plot, experiment B.
# Define CD4-CTL tag for seurat object of interest.
tmp.corrs <- c(rep(x='CD4-CTL', times=4), rep(x='Rest', times=3)); names(tmp.corrs) <- c(1:3, '6', '0', 4:5)
tmp.data <- as.character(srt.objs.list[[main.obj.6]]@meta.data[, clusts.lab[main.obj.6]])
tmp.data <- tmp.corrs[tmp.data]
srt.objs.list[[main.obj.6]]@meta.data[, 'ctl.tag'] <- tmp.data
# Define genes of interest.
these.markers <- c(
  'GZMB', 'PRF1', 'GNLY', 'NKG7', 'XCL1', 'XCL2' # CD4-CTL
)
# Plot.
tmp.ggplot <- dot.plot(
  seurat.obj=srt.objs.list[[main.obj.6]], features=these.markers, slot='data', do.norm=FALSE, ensembl=FALSE,
  groups.tag='ctl.tag', groups.order=c('CD4-CTL', 'Rest'), groups.of.int=NULL, filter.tag=NULL, groups.to.filter=NULL, keep=TRUE, na.rm=TRUE, feature.thold=NULL,
  this.color.scale=signatures.col.scale, col.min=NULL, col.max=NULL,
  scale.by='radius', dot.scale=6, size.min=NA, size.max=NA,
  file.name=NULL
)
# Output.
publish.plot(tmp.ggplot=tmp.ggplot, output.path=sup.6.path, file.name='HltyDonor-ExpB_Markers_Opt-A', type='pdf', blank.comp=blank.complement.2, do.legend=TRUE, height=4, width=3)

# ---> QC-related figures, experiment C.
# @ Distribution of the number of genes per cell across sequencing batches.
sec.meta.data[, donor.pop.tag:=paste(cell.pop.tag, donor.id.tag, sep=';')]
srt.objs.list[[sec.obj.6]]@meta.data[, 'donor.pop.tag'] <- sec.meta.data[, donor.pop.tag]
tmp.ggplot <- vln.plot(seurat.obj=srt.objs.list[[sec.obj.6]], feature='nFeature_RNA', groups.tag='donor.pop.tag', na.rm=TRUE, groups.of.int=NULL, groups.label=NULL, filter.tags=NULL, groups.to.filter=NULL, keep=TRUE, feature.thold=NULL, color='color', size.thold=0, file.name=NULL, adjust.val=1, trim.val=TRUE, plot.limits=NULL, add.points=FALSE, this.color.scale=NULL)
tmp.cols <- c(`SLAMF7+;SDBB-064`='#cc0000', `SLAMF7+;SDBB-085`='#800000', `SLAMF7-/CCR7+;SDBB-064`='#add8e6', `SLAMF7-/CCR7+;SDBB-085`='#00bfff')
tmp.ggplot <- tmp.ggplot + scale_fill_manual(values=tmp.cols) + scale_y_continuous(breaks=c(1000, 2500, 4000))
publish.plot(tmp.ggplot=tmp.ggplot, output.path=sup.6.path, file.name='HltyDonor-ExpC_QCs_GenesPerCellDist', type='pdf', blank.comp=blank.complement.1.1, do.legend=TRUE, height=4, width=5)

# ---> Dot plot, experiment D.
# Define genes of interest.
these.markers <- c(
  'GZMB', 'PRF1', 'GNLY', 'NKG7', 'XCL1', 'XCL2' # CD4-CTL
)
# Plot.
tmp.ggplot <- dot.plot(
  seurat.obj=srt.objs.list[[sec.obj.6]], features=these.markers, slot='data', do.norm=FALSE, ensembl=FALSE,
  groups.tag='cell.pop.tag', groups.order=NULL, groups.of.int=NULL, filter.tag=NULL, groups.to.filter=NULL, keep=TRUE, na.rm=TRUE, feature.thold=NULL,
  this.color.scale=signatures.col.scale, col.min=NULL, col.max=NULL,
  scale.by='radius', dot.scale=6, size.min=NA, size.max=NA,
  file.name=NULL
)
# Output.
publish.plot(tmp.ggplot=tmp.ggplot, output.path=sup.6.path, file.name='HltyDonor-ExpD_Markers_Opt-A', type='pdf', blank.comp=blank.complement.2, do.legend=TRUE, height=4, width=3)


### ------------- Other things that may fit in Figure 5 ------------- ###

# ---> UMAP depicting clusters for experiment C.
# Plot.
tmp.ggplot <- ggplot(data=sec.meta.data, aes_string(x='UMAP_1', y='UMAP_2', col=clusts.lab[sec.obj.6])) +
  geom_point(size=1.2, alpha=0.6) +
  scale_color_manual(values=clusts.cols[[sec.obj.6]]) +
  scale_y_continuous(expand=expansion(add=0.3), limits=c(-3, 4.6)) + scale_x_continuous(expand=expansion(add=0.3)) +
  labs(x='UMAP 1', y='UMAP 2', col='Cluster')
# Output.
publish.plot(tmp.ggplot=tmp.ggplot, output.path=fig.6.path, file.name='HltyDonor-ExpC_ClustersOnUMAP', type='tiff', blank.comp=blank.complement.3, do.legend=TRUE)

# ---> Cell fractions per SLAMF7 group or per donor (SLAMF7 experiments) for each cluster.
# Get data
tmp.data <- get.tag.analysis(seurat.obj=srt.objs.list[[sec.obj.6]], tmp.meta.data=as.data.frame(sec.meta.data), tmp.tag='cell.pop.tag', clusters.tag=clusts.lab[sec.obj.6], tag.reports.path=NA, reports.pref=NA, vals.cols=NULL, output.umaps=FALSE)
tmp.data <- as.data.frame(t(tmp.data), stringsAsFactors=TRUE); tmp.data$cluster <- row.names(tmp.data)
tmp.data <- tmp.data[tmp.data$cluster %in% names(clusts.cols[[sec.obj.6]]), ]
# Cell fractions per cluster for SLAMF7+ cells.
tmp.ggplot <- ggplot(data=tmp.data, aes(x=0, y=`SLAMF7+`, fill=cluster)) +
  geom_bar(stat='identity', position='fill', alpha=0.7) +
  scale_fill_manual(values=clusts.cols[[sec.obj.6]]) +
  coord_polar(theta = "y") + xlim(c(-2, 0.5)) +
  labs(x='', y='', fill='Cluster') +
  theme_void()
publish.plot(tmp.ggplot=tmp.ggplot, output.path=fig.6.path, file.name='HltyDonor_ExpC_SLAMF7Pos_CellFracPerCluster', type='pdf', blank.comp=blank.complement.4, do.legend=TRUE)
# Cell fractions per cluster for SLAMF7- cells.
tmp.ggplot <- ggplot(data=tmp.data, aes(x=0, y=`SLAMF7-/CCR7+`, fill=cluster)) +
  geom_bar(stat='identity', position='fill', alpha=0.7) +
  scale_fill_manual(values=clusts.cols[[sec.obj.6]]) +
  coord_polar(theta = "y") + xlim(c(-2, 0.5)) +
  labs(x='', y='', fill='Cluster') +
  theme_void()
publish.plot(tmp.ggplot=tmp.ggplot, output.path=fig.6.path, file.name='HltyDonor_ExpC_SLAMF7Neg_CellFracPerCluster', type='pdf', blank.comp=blank.complement.4, do.legend=TRUE)

# ---> Contrasting SLAMF7+ and SLAMF7- cell populations through TCR data.
# Get data.
tmp.data <- sec.meta.data[!is.na(clonotype.tag)]
# Plot and output.
tmp.ggplot <- vln.plot(seurat.obj=srt.objs.list[[sec.obj.6]], feature='clon.size.tag', groups.tag='cell.pop.tag', na.rm=TRUE, groups.of.int=NULL, groups.label=NULL, filter.tags=NULL, groups.to.filter=NULL, keep=TRUE, feature.thold=NULL, color='median', size.thold=0, file.name=NULL, adjust.val=20, trim.val=TRUE, plot.limits=NULL, add.points=FALSE, this.color.scale=NULL)
publish.plot(tmp.ggplot=tmp.ggplot, output.path=fig.6.path, file.name='HltyDonor_TCR-Repertoire_SizeBetweenCellPops_ExpC', type='pdf', blank.comp=blank.complement.1, do.legend=TRUE, width=5)
