###############################################################################
# Example data analysis                                                       #
# Literature dataset (SRA)                                                    #
# Macrophage stimulation timecourses                                          #
###############################################################################
folder <- "macrophage_stimulation_timecourse"
sra.id.vector <- "SRP056625"

hpc.mount <- "/path/to/HPC/partition/"

## Output directory for result website ##
webSiteDir <- "C:/xampp/htdocs/"

###############################################################################
## Parameters                                                                ##

species <- "mm"
project.code <- "p126"
project_id <- "sll126"
lab_id = "example_lab"
user_ids = c(
    "userGroup1", 
    "userGroup2", 
    "userGroup3"
)

original.data.folder <- fastqDir

paired.end <- FALSE
stranded <- TRUE
ref.genome <-  "GRCm38-release-86"

## If read length is unknown probe FASTQ files using ##
# zcat $filename | awk 'NR%2==0' | awk '{print length($1)}' | head
read.length <-  "50bp"


labname <- "Example"

count.table.headline = "TPM Values for all Samples"
count.table.sidelabel = "TPM"
heamap.headline.text = "Heatmap: Row-averaged TPM"


## Setting database parameters ##
prim.data.db <- "example_data"
db.user <- "database_user"
db.pwd <- db.pwd
host <- "database_url"

ref.cat.db = "reference_categories_db_new"
ref.cat.db.table <- "reference_categories_db_new"
cat.ref.db.table = paste0(project_id, "_cat_reference_db_table")
cat.ref.db <- prim.data.db

PCAdbTableName <- paste0(project.code, "_PCA")
rnaseqdbTableName <- paste0(project.code, "_rna_seq_table")
enriched.categories.dbTableName = paste0(project.code, "_enriched_categories_table")
lab.categories.table = "sl_lab_categories"


## Others ##
timecourse.units <- "min"

## Create design file name ##
design.file  <- paste0(project.code, ".design.file.txt")

## Create required directories ##
datadir <- paste0(
    hpc.mount,
    "Projects/",
    folder, 
    "/basedata/"
)

localDataDir <- paste0(
    hpc.mount,
    "Projects/",
    folder, 
    "/basedata/"
)

workdir <- paste0(
    hpc.mount,
    "Projects/",
    folder, 
    "/workdir/"
)

localWorkDir <- paste0(
    hpc.mount,
    "Projects/",
    folder, 
    "/workdir/"
)

fastqDir <- paste0(
    hpc.mount,
    "Projects/",
    folder,
    "/FASTQ_files/"
)

## Set indices ##
ENSG.index.mouse <- paste0(
    hpc.mount,
    "/path/to/genomes/mus_musculus/ensembl/GRCm38/release-86/genome_idx/rsem/star/", 
    read.length,
    "/genome"
)

ENSG.index.human = paste0(
    "/path/to/genomes/homo_sapiens/ensembl/GRCh38/release-86/genome_idx/rsem/star/",
    read.length,
    "/genome"
)

if (species == "hs"){
    ref.index <- ENSG.index.human
    primary.alignment.gene.id <- "ENSG"
    gene.id.column = "hgnc_symbol"
} else if (species == "mm"){
    ref.index <- ENSG.index.mouse 
    primary.alignment.gene.id <- "ENSMUSG"
    gene.id.column = "mgi_symbol"
}

if (species == "hs"){
    GTFfile="/path/to/genomes/homo_sapiens/ensembl/GRCh38/release-86/gtf/Homo_sapiens.GRCh38.86.rnaseqc.gtf"
    genome.fa="/path/to/genomes/homo_sapiens/ensembl/GRCh38/release-86/genome/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa"
    rRNAfile <- "/path/to/genomes/homo_sapiens/ensembl/GRCh38/release-86/gtf/Homo_sapiens.GRCh38.86.rRNA.list"
    gene.id.column <- "hgnc_symbol"
} else if (species == "mm"){
    GTFfile="/path/to/genomes/mus_musculus/ensembl/GRCm38/release-86/gtf/Mus_musculus.GRCm38.86.rnaseqc.gtf"
    genome.fa="/path/to/genomes/mus_musculus/ensembl/GRCm38/release-86/genome/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa"
    rRNAfile <- "/path/to/genomes/mus_musculus/ensembl/GRCm38/release-86/gtf/Mus_musculus.GRCm38.86.rRNA.list"
    gene.id.column <- "mgi_symbol"
}



###############################################################################
# Prepare annotation tables                                                   #
###############################################################################
if (species == "hs"){
  gene.id.table <- paste0(
    hpc.mount,
    "Projects/reference_data/gene_id_annotation_files/",
    "20170202.ENSG.mgi.entrez.uniprot.description.mgi.table.txt"
  )
} else if (species == "mm"){
    gene.id.table <- paste0(
        hpc.mount,
        "Projects/reference_data/gene_id_annotation_files/",
        "20170124.ENSMUSG.mgi.entrez.uniprot.description.hgnc.table.txt"
    )
}


###############################################################################
## Download FASTQ files from SRA                                             ##
library(SBwebtools)

## Create shell script to download FASTQfiles ##
# For this, create EMPTY FASTQ_files folder
setwd(localWorkDir)
createSRAdownloadScript(
    sra.id.vector = sra.id.vector,
    module.load.cmd = "module use path/to/module/library;module load sratoolkit/2.8.2-1",
    fastqDir = fastqDir
)

##  Done downloading FASTQ files from SRA                                    ##
###############################################################################


###############################################################################
## Load   design file                                                        ##

setwd(localWorkDir)
df.design <- read.delim(
    design.file, 
    header = TRUE, 
    sep = "\t",
    stringsAsFactors = FALSE
)

## Reorder df design ##
# order by sample group - sample name
df.design <- df.design[order(
    df.design$dataseries, 
    df.design$timepoint, 
    df.design$sample.group, 
    df.design$sample.id
),]

print("Design file ordered correctly? If in doubt, check!")
## Done creating/loading design file                                         ##
###############################################################################

###############################################################################
## Load libraries                                                            ##

library(SBwebtools)

##                                                                           ##
###############################################################################

###############################################################################
## Create shell script to rename files                                       ##
setwd(localWorkDir)
sink("create.fastq.softlinks.sh")
for (i in 1:nrow(df.design)){
    string <- paste0(
        "ln -s ",
        df.design$original.NGS[i],
        " ",
        fastqDir,
        df.design$NGS[i]
)
    cat(string);cat("\n")
}

sink()

## Run ##
print("sh conv.create.fastq.softlinks.sh")

###############################################################################
## Create strings for documentation                                          ##
## fastq folder
print(paste0("FASTQ file location:", original.data.folder))

print(paste0(
    "FASTQ file names:", 
    paste0(
        gsub(original.data.folder, "",df.design$original.NGS), 
        collapse = ", "
    )
))

print(paste0(
    "Project FASTQ file names:", 
    paste0(
        df.design$NGS, 
        collapse = ", "
    )
))

if (paired.end){
    print("Alignment mode: paired-end")
} else {
    print("Alignment mode: single-end")
}

if (stranded){
    print("This is a stranded dataset")
} else {
    print("This is not a stranded dataset")
}

print(
    paste0(
        "Reference genome/transcriptome: ",
        ref.genome
    )
)


## End create strings for documentation                                      ##
###############################################################################

###############################################################################
# Create trim galore shell script                                             #
###############################################################################
setwd(localWorkDir)
df.design <- create.trim.galore.shell.script(
    df.design, #design file with FASTQ column
    project.code = project.code,
    fastqDir  = fastqDir,
    localWorkDir = localWorkDir,
    paired.end = paired.end, 
    module.to.load = "Trim_Galore/0.4.2-foss-2016b",
    min.length = 25,
    quality = 20
)


# Write df.design to file so it can be re-read once the alignment is done
setwd(localWorkDir)
write.table(
    df.design, 
    design.file, 
    row.names= FALSE, 
    sep="\t"
)


###############################################################################
# Align                                                                       #
###############################################################################
setwd(localWorkDir)

create.alignment.shell.script(
    localWorkDir = localWorkDir,
    df.design = df.design,
    df.design.fastq.column = "FASTQ_trimgalore",
    workDir = workdir,
    fastqDir =fastqDir,
    sampleDir = paste0(
        hpc.mount,
        "Projects/",
        folder,
        "/workdir/RSEM/Ensembl"
    ),
    outputEnsDir = paste0(
        workdir, 
        "RSEM/Ensembl/"
    ),
    FASTQC = paste0(workdir, "FASTQC/"),
    logDir = paste0(workdir, "logs/"),
    msubDir = paste0(workdir, "msub/"),
    localMsubDir = paste0(localWorkDir, "msub/"),
    msubFile = paste0(project.code, ".fastqc.and.rsem.all.samples.sh"),
    sample.order = "",
    forward_prob = "--forward-prob 0.0",
    ENSG.index.human = ENSG.index.human,
    ENSG.index.mouse = ENSG.index.mouse,
    species = species,
    paired.end = paired.end
)


# Run
print(paste0("sh ", "conv.fastqc.and.rsem.all.samples.sh"))

###############################################################################
# Prepare RNAseQC script                                                      #
###############################################################################
setwd(localWorkDir)
create.rnaseqc.script(
    df.design = df.design,
    sample.column = "sample.id",
    project.code = project.code,
    project="SB_RNAseqQC",
    basedataDir=paste0(workdir, "RSEM/Ensembl"),
    bam.suffix = "STAR.genome.bam",
    GTFfile=GTFfile,
    rRNAfile=rRNAfile,
    genome.fa=genome.fa,
    paired.end = paired.end
)



###############################################################################
# Do differential gene expression analysis                                    #
###############################################################################

###############################################################################
# Prepare TPM and FPKM tables                                                 #
###############################################################################
## Make sure df.design is ordered properly ##
samples = as.vector(unique(df.design$sample.id))
files = paste(localWorkDir, "RSEM/Ensembl/", samples, ".genes.results", sep="")

list.tpm.fpkm <- create.tpm.and.fpkm.tables(
    workdir = localWorkDir, 
    samples = samples,
    files = files
)

df.tpm  <- list.tpm.fpkm$df.tpm
df.fpkm <- list.tpm.fpkm$df.fpkm

###############################################################################
## Do differential gene expression analysis                                  ##
###############################################################################
## Review raw counts file column names ##
col.names <- names(df.raw <- read.delim(
    paste0(
        localWorkDir, "RSEM/",
        count.data.file
    ), 
    header=TRUE, 
    sep="\t", 
    stringsAsFactors = FALSE
    )
)


setwd(localWorkDir)

## Do differential gene expression analysis using DESEQ2 ##
df.summary <- do.differential.expression.analyis(
    raw.count.dir = paste0(localWorkDir, "RSEM"),
    count.data.file = count.data.file,             #count data filename
    DEseq2Dir = paste0(localWorkDir, "DESeq2"),    #directory for results
    df.design = df.design,                         #df.design
    gene.id = primary.alignment.gene.id,                           #primary gene id after alignment 
    batch.mode = FALSE, #if true, df.design needs to contain a 'replicate' column 
    string.to.be.deleted.in.raw.counts.columns = paste0("X", gsub("/", ".",paste0(workdir, "RSEM/Ensembl/")))
)


###################################
# Upload pca table to database    #
###################################
setwd(localWorkDir)
pca.dir <- paste0(
    localWorkDir, 
    "DESeq2", 
    "/pca.table.txt"
)

df.pca <- read.delim(
    pca.dir, 
    header = TRUE, 
    sep="\t", 
    stringsAsFactors = FALSE
)

names(df.pca) <- gsub(
    "[.]", 
    "_", 
    names(df.pca)
)

## Upload to database ##
 
upload.pca.table.to.db(
    df.pca = df.pca,
    host = "www.biologic-db.org",
    prim.data.db = prim.data.db,
    password = db.pwd,
    db.user = db.user, 
    PCAdbTableName = PCAdbTableName
)

# Finished creating pca database table

###############################################################################
# Add -logFC p-values to the datatable (for use in Volcano plots)
###############################################################################
# Get padj columns
padj  <- names(df.summary)[grep("_padj_", names(df.summary))]
lg10p <- gsub("padj", "lg10p", padj) 

for (i in 1:length(padj)){
  temp <- -1*log10(as.numeric(df.summary[,padj[i]]))
  temp[temp >= 50] = 50
  df.summary[,lg10p[i]] <- temp
}


###############################################################################
# Add TPM values to data table                                               #
###############################################################################
new.order.tpm <- new.order[grep("norm_counts_", new.order)]
new.order.tpm <- gsub("norm_counts_", "", new.order.tpm)
new.order.tpm <- paste0(new.order.tpm, ".TPM")
new.order.tpm <- append("gene_id", new.order.tpm)
df.tpm <- df.tpm[,new.order.tpm]
names(df.tpm) <- gsub(".TPM", "",names(df.tpm))
names(df.tpm)[2:length(df.tpm)] <- paste0("TPM________",names(df.tpm)[2:length(df.tpm)])
names(df.tpm) <- gsub("gene_id", primary.alignment.gene.id,names(df.tpm))
df.summary <- merge(df.summary, df.tpm, by.x = primary.alignment.gene.id, by.y = primary.alignment.gene.id)

###############################################################################
# Exchange TPM and norm_counts in df.summary                                  #
###############################################################################
names(df.summary) <- gsub("norm_counts_", "Xtemp_", names(df.summary))
names(df.summary) <- gsub("TPM________", "norm_counts_", names(df.summary))
names(df.summary) <- gsub( "Xtemp_", "raw_counts__", names(df.summary))


###############################################################################
# Upload to website                                                           #
###############################################################################
#library(SBwebtools)
setwd(localWorkDir)

###############################################################################
# Prepare database table                                                      #
###############################################################################

## Select for heatmap: abs change of at least 0.5 in any contrast ##
database.table <- datatable.to.website.ptm(
    df.data = df.summary, 
    gene.id.column = primary.alignment.gene.id, 
    n.cluster.genes = 2000, 
    count.data = TRUE, 
    gene.id.table = gene.id.table,
    add.uniprot.column = TRUE, 
    selector4heatmap.cols = "norm_counts",
    heatmap.preprocessing = "lg2.row.avg", # possible: "lg2", "lg2.row.avg", "none"
    hm.cut.off = 10,
    n.hm.cluster = 10
)

## Remove all genes with low expression ##
## Require on average 1 TPM per sample
database.table <- database.table[database.table$count_cut_off >= 1,]

###############################################################################
# Upload to database                                                          #
###############################################################################
cmd.vec <- upload.datatable.to.database(
    host = host, 
    user = db.user,
    password = db.pwd,
    prim.data.db = prim.data.db,
    dbTableName = rnaseqdbTableName,
    df.data = database.table, #[database.table$count_cut_off > 5, ],
    db.col.parameter.list = list(
        "VARCHAR(255) CHARACTER SET latin1 COLLATE latin1_swedish_ci" = c("gene_description"),
        "VARCHAR(50) CHARACTER SET latin1 COLLATE latin1_swedish_ci" = c("ENSG", "ENSMUSG", "hgnc_symbol", "mgi_symbol", "uniprot", "entrezgene","display_ptm", "^sequence_window", "p_site_env","for_GSEA_gene_chip","associated_gene_name", "gene_type"),
        "VARCHAR(1) CHARACTER SET latin1 COLLATE latin1_swedish_ci" = c("ppos", "amino_acid", "charge","known_site"),
        "BIGINT(8) NULL DEFAULT NULL" = c("row_names"),
        "INT(8) NULL DEFAULT NULL" = c("CoVarOrder","row_id", "cluster_order","cluster_id", "count_cut_off", "^position$", "raw_counts"),
        "DECIMAL(6,3) NULL DEFAULT NULL" = c("^CoVar$","NES", "logFC", "lg2_avg", "intensity", "^int", "iBAQ","^localization_prob$", "stat", "lg10p"),
        "DECIMAL(6,5) NULL DEFAULT NULL" = c("padj", "pvalue","^pep$"),
        "DECIMAL(6,1) NULL DEFAULT NULL" = c("norm_counts")
    ),
    new.table = TRUE
)

  
###############################################################################
# Do GSEA                                                                     #
###############################################################################
#Create rnk files for all logFC comparisons
#Remove lowly expressed genes from table. Cut off: Less than 1/counts per sample. Should be changed to TPM cut-off
# In this case, an average of 1 TPM is required for a sample to qualify

create.gsea.rnk.files(
    localWorkDir, 
    df.dataTable = database.table,
    GSEA.colum.type = "logFC",
    gene.symbol.column.name = "mgi_symbol"
)



## Run GSEA on the command line or in JAVA interface ##


## Collect and process GSEA output files ##

gsea.cat.lines <- create.GSEA.table(
    GSEADir = paste0(localWorkDir, "GSEA"),
    logFC.column.name = "logFC",
    host = host,
    refdbname= "reference_categories_db_new",
    refDBTableName = enriched.categories.dbTableName,
    db.user = db.user,
    db.password = db.pwd, 
    tables = tables,
    df.dataTable = database.table
)
            

## Create GSEA plot directory by executing farm command ##
## Check [project.code.enrichment.file and copy to outputs]

###############################################################################
## Create reference category selection for this project                      ##

library(SBwebtools)
df.data <- import.db.table.from.db(
    dbname = prim.data.db,
    dbtable = rnaseqdbTableName,
    password = db.pwd
)

## Get relevant genes = genes that exhibit statistically siginficant alterations
# in at least one logFC comparison
relevant.genes <- vector(mode="character", length=0)
df.sel <- df.data
padj.cols <- names(df.sel)[grep("padj", names(df.data))]
for (i in 1:length(padj.cols)){
    temp.vec <- as.vector(df.sel[df.sel[,padj.cols[i]] < 0.05,gene.id.column])
    relevant.genes <- unique(append(relevant.genes, temp.vec))
}

print(paste0(length(relevant.genes), " relevant genes selected."))
## Further reduce comlexity ##
#relevant.genes <- vector(mode="character", length=0)
#logFC.cols <- names(df.sel)[grep("logFC", names(df.data))]
#cut.off <- 1.0
#for (i in 1:length(logFC.cols)){
#  temp.vec <- as.vector(df.sel[((abs(df.sel[,logFC.cols[i]]) > cut.off) & (df.sel[,padj.cols[i]] < 0.05)),gene.id.column])
#  relevant.genes <- unique(append(relevant.genes, temp.vec))
#}
#print(paste0(length(relevant.genes), " relevant genes selected."))

sel.vec = c(
    "cat_name",  
    "cat_id",
    "hgnc_symbol",
    gene.id.column,
    "cat_item_size"
)

###############################################################################
## Create lab specific reference categories                                  ##

## Lab categories ##
df.ref <- import.db.table.from.db(
    dbname = ref.cat.db,
    dbtable = "es_lab_categories",
    password = db.pwd
)

## Remove all project temp categories from df.ref
df.ref <- df.ref[-grep("temp_", df.ref$cat_type),]

df.temp <- add2labCatSelectionDBtable(
    df.ref = df.ref,
    cat_group_name = paste0(labname, " lab"),
    reference.gene.vector = relevant.genes,
    ref.gene.vec.id = gene.id.column,
    cat_views = NA
)

df.db.table <- df.temp

###############################################################################
## ChIP/ChEA                                                                 ##
###############################################################################
df.ref <- import.db.table.from.db(
    dbname = ref.cat.db,
    dbtable = "ChEA_2016",
    password = db.pwd
)

df.temp <- add2labCatSelectionDBtable(
    df.ref = df.ref,
    cat_group_name = "ChIP",
    reference.gene.vector = reference.gene.vector,
    ref.gene.vec.id = gene.id.column,
    cat_views = NA
)

df.db.table <- rbind(df.temp, df.db.table)

###############################################################################
##mysigdb c3                                                                 ##
###############################################################################
df.ref <- import.db.table.from.db(
    dbname = ref.cat.db,
    dbtable = "mysigdb_c3_TF_targets",
    password = db.pwd
)

df.temp <- add2labCatSelectionDBtable(
    df.ref = df.ref,
    cat_group_name = "TF target collection",
    reference.gene.vector = reference.gene.vector,
    ref.gene.vec.id = gene.id.column,
    cat_views = NA
)

df.db.table <- rbind(df.temp, df.db.table)

###############################################################################
##TRANSFAC/JASPAR                                                            ##
###############################################################################
df.ref <- import.db.table.from.db(
    dbname = ref.cat.db,
    dbtable = "TRANSFAC_and_JASPAR_PWMs",
    password = db.pwd
)

df.temp <- add2labCatSelectionDBtable(
    df.ref = df.ref,
    cat_group_name = "JASPAR TF binding motifs",
    reference.gene.vector = reference.gene.vector,
    ref.gene.vec.id = gene.id.column,
    cat_views = NA
)

df.db.table <- rbind(df.temp, df.db.table)
# to be added

###############################################################################
##GEO single gene pertubations                                               ##
###############################################################################
df.ref <- import.db.table.from.db(
    dbname = ref.cat.db,
    dbtable = "SGP_from_GEO_up_down_combined",
    password = db.pwd
)

df.temp <- add2labCatSelectionDBtable(
    df.ref = df.ref,
    cat_group_name = "Single Gene Pertubations",
    reference.gene.vector = reference.gene.vector,
    ref.gene.vec.id = gene.id.column,
    cat_views = NA
)

df.db.table <- rbind(df.temp, df.db.table)

###############################################################################
## Pathways                                                                  ##
###############################################################################
## MysigDB ##
df.ref <- import.db.table.from.db(
    dbname = "reference_categories_db_new",
    dbtable = "mysigdb_c2_1329_canonical_pathways",
    password = db.pwd
)

df.temp <- add2labCatSelectionDBtable(
    df.ref = df.ref,
    cat_group_name = "Pathways",
    reference.gene.vector = reference.gene.vector,
    ref.gene.vec.id = gene.id.column,
    cat_views = NA
)

df.db.table <- rbind(df.temp, df.db.table)


###############################################################################
## Oncogenic Signatures                                                    ##
###############################################################################
## MysigDB ##
df.ref <- import.db.table.from.db(
    dbname = "reference_categories_db_new",
    dbtable = "mysigdb_c6_oncogenic_signatures",
    password = db.pwd
)

df.temp <- add2labCatSelectionDBtable(
    df.ref = df.ref,
    cat_group_name = "Oncogenic Signatures",
    reference.gene.vector = reference.gene.vector,
    ref.gene.vec.id = gene.id.column,
    cat_views = NA
)

df.db.table <- rbind(df.temp, df.db.table)

###############################################################################
## Protein Complexes                                                         ##
###############################################################################
df.ref <- import.db.table.from.db(
    dbname = "reference_categories_db_new",
    dbtable = "networkcategories",
    password = db.pwd
)

df.temp <- add2labCatSelectionDBtable(
    df.ref = df.ref,
    cat_group_name = "Protein-Protein Interactions",
    reference.gene.vector = reference.gene.vector,
    ref.gene.vec.id = gene.id.column,
    cat_views = NA
)

df.db.table <- rbind(df.temp, df.db.table)
###############################################################################
## Filter cat dataset                                                        ##
###############################################################################
## Remove all datasets with less than 5 proteins matched
dim(df.db.table)
df.db.table <- df.db.table[df.db.table$cat_count >= 5,]
dim(df.db.table)
#df.db.table <- df.db.table[df.db.table$cat_weight > 0.3,]
#dim(df.db.table)


###############################################################################
## Add full temp categories associated with this project                     ##
## Lab categories ##
df.ref <- import.db.table.from.db(
    dbname = ref.cat.db,
    dbtable = lab.categories.table,
    password = db.pwd
)

## Remove all project temp categories from df.ref
df.ref <- df.ref[grep(paste0("temp_", project_id), df.ref$cat_type),]

df.temp <- add2labCatSelectionDBtable(
    df.ref = df.ref,
    cat_group_name = "This project",
    reference.gene.vector = reference.gene.vector,
    ref.gene.vec.id = gene.id.column,
    cat_views = NA
)

df.db.table <- rbind(df.temp, df.db.table)
## Done adding temp categories associated with this project
###############################################################################


###############################################################################
## Upload to database                                                        ##
###############################################################################
upload.datatable.to.database(
    host = host, 
    user = db.user,
    password = db.pwd,
    prim.data.db = cat.ref.db,
    dbTableName = cat.ref.db.table,
    df.data = df.db.table,
    db.col.parameter.list = list(
        "VARCHAR(255) CHARACTER SET utf8 COLLATE utf8_general_ci" = c("comments_1","cat_name"),
        "VARCHAR(50) CHARACTER SET utf8 COLLATE utf8_general_ci" = c("cat_id", "cat_group"),
        "VARCHAR(1) CHARACTER SET utf8 COLLATE utf8_general_ci" = c("ppos", "amino_acid", "charge","known_site"),
        "BIGINT(8) NULL DEFAULT NULL" = c("row_names"),
        "INT(8) NULL DEFAULT NULL" = c("cat_item_size", "cat_count", "cat_views"),
        "DECIMAL(6,3) NULL DEFAULT NULL" = c("cat_weight"),
        "DECIMAL(6,5) NULL DEFAULT NULL" = c("padj", "pvalue","^pep$")
    )
)



##  End creating reference category selection for this project               ##
###############################################################################


###############################################################################
# Create timecourse specifications                                            #

timecourse.cat.lines <- create.timecourse.cat.lines(
    df.design
)


###############################################################################
# Create microwebsite                                                         #
###############################################################################
sample.order <- names(database.table)[grepl("norm_counts_", names(database.table))]
            
library(RColorBrewer)
selcol <- colorRampPalette(brewer.pal(12,"Set3"))
            
color.sel = selcol(10)
color.sel = c(rbind(color.sel, color.sel, color.sel))
            
#library(SBwebtools)
setwd(localWorkDir)

create.website.parameters(
    df.data = database.table,
    gene.id.column = gene.id.column,
    ptm.colum = "",
    lab_id = lab_id , 
    user_ids = user_ids, 
    project_id = project_id, 
    download_result_table = paste0(project.code, ".result.table"),
    download_cat_enrichment_table = paste0(project.code, ".enriched.categories"),
    database = prim.data.db,
    reference_categories_db = "reference_categories_db_new",
    labname = labname,
    rnaseqdbTableName,
    lab.categories.table = lab.categories.table,
    sample.order = sample.order, #set to "" to go with default sorting
    count.sample.colors = "",
    count.table.headline = count.table.headline,
    count.table.sidelabel = count.table.sidelabel,
    webSiteDir = webSiteDir,
    heamap.headline.text = heamap.headline.text,
    upper_heatmap_limit = 3, 
    lower_heatmap_limit = -3,
    slider.selection.name = "logFC",
    presentation.file = paste0(project.code, ".project.presentation.pptx"),
    number_of_slides = 18,
    default.sequence = "",
    use.logFC.columns.for.heatmap = FALSE,
    peptide.view.link = "",
    create.2d.scatterplot.button = TRUE,
    low_highlight = -1,
    high_highlight =1,
    display.qc = TRUE,
    display.pca = TRUE,
    pca.table.name = PCAdbTableName, 
    gsea.cat.lines = gsea.cat.lines,
    timecourse.cat.lines = timecourse.cat.lines,
    venn.slider.selector.strings = c("_logFC", "_padj"),
    plot.selection.strings = c("_logFC", "_lg10p"), #NA, #strings to grep from col names for plot display
    plate.view.db.table = NA,
    plate.view.column.vec = NA,
    cat.seletion.table.vec = c(cat.ref.db, cat.ref.db.table)
)
            
            
###########################################
# Do multiQC and add to QC section 
###########################################

## Run multiqc ##
print("module use /path/to/modules/all;module load multiqc/0.9.1-2016b-Python-2.7.12;multiqc.")


###############################################################################
## Module Add metacore results                                               ##

## Retrieve dataset from database for metacore ##
df.metacore <- import.db.table.from.db(
    password = db.pwd,
    dbname = prim.data.db,
    dbtable = rnaseqdbTableName
)

## select padj and logFC columns only ##
sel.vec <- names(df.metacore)[grep("contrast_", names(df.metacore))]
sel.vec <- sel.vec[-grep("stat", sel.vec)]
sel.vec <- sel.vec[-grep("lg10p", sel.vec)]
## Remove contrast_x_ prefix > remove front 11 characters
sel.vec <- append(gene.id.column, sel.vec)
df.metacore <- unique(df.metacore[,sel.vec])

## Remove contrast_x_ tag from column labels ##
names(df.metacore)[grep("contrast", names(df.metacore))] <- substr(
    names(df.metacore)[grep("contrast", names(df.metacore))],
    12,
    100
)

setwd(localWorkDir)

write.table(
    df.metacore, 
    paste0(
        project.code, 
        "metacore.input.file.txt"
    ), 
    row.names = FALSE, 
    sep="\t"
)

print(
    paste0(
        "Upload ",   
        project.code, 
        ".metacore.input.file.txt", 
        " to metacore"
    )
)

## Open file in excel and save as metacore.input.file.xls (2003)
print("Perform enrichment analysis by subsetting data to logFC cut off: 1, padj 0.05")
## Select Workflows & Reports
## Select Enrichment analysis
## Set threshold: Threshold 1 p-value 0.05 direction both
## Run analysis
## If necessary, repeat with lower theshold
## If successful hit Get report button and safe as 
print(paste0(project.code, ".metacore.results.enrichment.analysis.xls"))

## Next do transcription factor target analysis ##
## Select One-click Analysis > Transcription Factors
## Set FDR threshold to 0.05
print(paste0("Export result table as: ", project.code, ".metacore.TF.analysis.", "logFC_nonAligned_vs_aligned"))

## Save results as p111.metacore.result.

## For selected TF targets, export MC results an curate into project category ##
## Select transcription factor of interest (Object name)
## Limit selection: Direction: Outgoing Effect Activation Mechanism influence on expression and transcription
## regulation >> Aplly >> Build Network
## Select additional options
## Pre filters Interaction types transcription regulation
## Additional options Directions: Downstream
## Hit build network
## Select all >> File >> Export >> safe as [TFname.mc.targets.xls]


## End Module add metacore results                                           ##
###############################################################################

