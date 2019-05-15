###############################################################################
# Ley lab                                                                     #
# Experiment:                                                                 #
# Primary data anlalysis                                                      #
# Ensembl release-89                                                          #
###############################################################################

## The sample sample.ids.txt file needs to have the following columns:
# sampleID	sample.name	sample.id


rm(list = ls())

###############################################################################
## Initialize S4 object                                                      ##


if (dir.exists("/Volumes/babs/working/boeings/")){
    hpc.mount <- "/Volumes/babs/working/boeings/"
} else if (dir.exists("Y:/working/boeings/")){
    hpc.mount <- "Y:/working/boeings/"
} else if (dir.exists("/camp/stp/babs/working/boeings/")){
    hpc.mount <- "/camp/stp/babs/working/boeings/"
} else {
    hpc.mount <- ""
}

source(
    paste0(
        hpc.mount,
        "Stefan/protocol_files/github/boeings/packages/cDev.package.SBwebtools.NEW201905.r"
    )
)



###############################################################################
## Create S4 object for this analysis                                        ##

Obio = new(
    "bioLOGIC",
    parameterList = list(
        "user_ids" = c(
            "project", 
            "sl.lab.all",
            "joan.manils",
            "paul.osullivan",
            "louise.webb",
            "sl.lab.all"
        ),
        "folder" = "leys/joan.manils/299_SLL_JM_bulk_RNA_Seq_psoriasis_normal_treated_GSE74697",
        "sra.id.vector" = "",
        "gse.id.vector" = "",
        "lims.id" = "",
        "machine" = "",
        "experiment.type" = "rna_seq",    # any
        "species" = "mus_musculus", # Species has to be "mus_musculus", "homo_sapiens", "danio_rerio" 
        "release" = "release-89", #release-86, release-89
        "project_id" = "sll299",
        "labname" = "Ley",
        "db.user" = "boeings",
        "host" = "www.biologic-db.org",
        "timecourse.units" = "",
        "count.table.headline" = "TPM Values for all Samples",
        "count.table.sidelabel" = "TPM",
        "heamap.headline.text" = "Heatmap: Row-averaged TPM",
        "paired.end" = TRUE,
        "stranded" = TRUE,
        "loadR" = "module purge;module load R/3.3.1-foss-2016b-bioc-3.3-libX11-1.6.3;R;",
        "ref.cat.db" = "reference_categories_db_new",
        "ref.cat.db.table" = "reference_categories_db_new",
        "read.length" =  "100bp",
        "webSiteDir" = "",
        "pathToSeqStorageFolder" = c(
            "/camp/stp/babs/working/boeings/Projects/leys/joan.manils/299_SLL_JM_bulk_RNA_Seq_psoriasis_normal_treated_GSE74697/FASTQ_files/"
        ),
        "addFullTPMtable" = FALSE,
        "hpcMount" = ""
    )
)

## Done creating S4 object for this analysis                                 ##
###############################################################################

###############################################################################
## Set database password                                                     ##
if (exists("db.pwd")){
    print("Database password is set.")
} else {
    ## Set database password ##
    library(rstudioapi)
    db.pwd <- rstudioapi::askForPassword("Please enter database password")
}

##                                                                           ##
###############################################################################

## Check if S4LvariableListing has been saved before ##

## Determine sequencer for report (from Harshil Patel)
## Crick Sequencing
## more /camp/stp/babs/working/data/crick_runs_instruments.csv | grep "190205_K00102_0309_AH37TNBBXY"
## CRUK Sequencing
## more /camp/stp/babs/working/data/lif_runs_instruments.csv | grep "190205_K00102_0309_AH37TNBBXY"
## Gathering flow cell details
## https://github.com/10XGenomics/supernova/blob/master/tenkit/lib/python/tenkit/illumina_instrument.py#L12-L45

## Determine read length of unknown fastq files
## zcat SRR2926047_pass_1.fastq.gz | awk '{if(NR%4==2) print length($1)}' | head


###############################################################################
## Set HPC access point                                                      ##

Obio <- setMountingPoint(Obio)

## Done setting HPC access point                                             ##
###############################################################################

###############################################################################
## Set Genome and reference table                                            ##

Obio <- setCrickGenomeAndGeneNameTable(Obio)

## Create required folder structure
Obio <- createAnalysisFolders(
    Obio,
    baseDir="/camp/stp/babs/working/boeings/Projects/",
    localBaseDir = "Y:/working/boeings/Projects/"
)


## Done setting genome and reference table                                   ##
###############################################################################


###############################################################################
## Optional module:Download FASTQ files from SRA                             ##

## End: Optional module:Download FASTQ files from SRA                        ##
###############################################################################


###############################################################################
## Organizing FASTQ files                                                    ##

concatenationRequired <- organizeFastqFiles(
    baseMount = "Y:/working/",
    pathToSeqStorageFolder =  Obio@parameterList$pathToSeqStorageFolder,
    fastqOutputDir = Obio@parameterList$fastqDir,
    localWorkDir = Obio@parameterList$localWorkDir
)

## Done organizing FASTQ files                                               ##
###############################################################################

###############################################################################
## Assemble design file                                                      ##

if (concatenationRequired){
    print("Concatenation script executed?")
    Obio@parameterList$pathToSeqStorageFolder = Obio@parameterList$fastqDir
} 

###############################################################################
## Create Design File                                                        ##

dfDesign <- createDesignFileCrickASFsamples(
    pathToSeqStorageFolder = Obio@parameterList$pathToSeqStorageFolder,
    FNsampleAnnotation = paste0(
        Obio@parameterList$localWorkDir,
        "sample.ids.txt"
    ),
    paired.end = Obio@parameterList$paired.end,
    fastqDir = Obio@parameterList$fastqDir,
    baseMount = gsub(
        "boeings/",
        "",
        Obio@parameterList$hpcMount
    )
)

## Add DGE comparisons ##
unique(dfDesign$sample.group)

dfDesign <- addDGEcomparisons2DesignFile(
    dfDesign = dfDesign,
    comparisonList = list(
        "S_SPU_vs_SN" = c("SkinPsoriasis_Unt", "SkinNormal_Unt"),
        "S_STR_vs_SPU" = c("SkinPsoriasis_adalimuMab", "SkinPsoriasis_Unt"),
        "S_STR_vs_SUT" = c("SkinPsoriasis_adalimuMab", "SkinNormal_Unt")
    )
)

## Add design file to the Obio object ##
Obio@dfDesign <- dfDesign

setwd(Obio@parameterList$localWorkDir)
write.table(dfDesign, Obio@parameterList$designFN, row.names = FALSE, sep = "\t")
    

## Done creating design file                                                 ##
###############################################################################
## Done data organization                                                    ##
###############################################################################


###############################################################################
## Prepare Alignment and QC scripts                                          ##

## Done Preparing Alignment and QC scripts                                   ##
###############################################################################

###############################################################################
## Load   design file                                                        ##

setwd(localWorkDir)

dfDesign <- read.delim(
    Obio@parameterList$designFN, 
    header = TRUE, 
    sep = "\t",
    stringsAsFactors = FALSE
)


## Done creating/loading design file                                         ##
###############################################################################



###############################################################################
## Create analysis scripts                                                   ##
setwd(localWorkDir)
returnList <- createbulkRNASeqAnalysisBashScripts(
    dfDesign = dfDesign,
    shellScriptVector = as.vector(NULL, mode = "character")
)

## Run analysis BASH script ##

## Done creating analysis bash scripts                                       ##
###############################################################################


## Save workspace ##
workspaceFN <- paste0(
    localWorkDir,
    project_id,
    ".RData"
)

save.image(workspaceFN)

###############################################################################
## Stopping point                                                            ##
###############################################################################

## Load workspace ##
workspaceFN <- paste0(
    localWorkDir,
    project_id,
    ".RData"
)

load(workspaceFN)

## Save workspace, so you can resume once the RSEM alignment and             ##
## RNASeqQC is done                                                          ##
###############################################################################


###############################################################################
# Do differential gene expression analysis                                    #
###############################################################################

###############################################################################
#Generate input matrix from expression results                                #
###############################################################################
#At this stage human samples from previous alignments will be integrated
#Get full sample list
setwd(localWorkDir)

count.data.file <- paste0(
    localWorkDir,
    "RSEM/",
    project.code, 
    ".count.data.txt"
)

###############################################################################
# Prepare TPM and FPKM tables                                                 #
###############################################################################
## Make sure dfDesign is ordered properly ##
samples <- as.vector(
    unique(
        dfDesign$sample.id
        )
)

files <- paste(
    localWorkDir, 
    "RSEM/Ensembl/", 
    samples, 
    ".genes.results", 
    sep=""
)


#library(SBwebtools)
list.tpm.fpkm <- create.tpm.and.fpkm.tables(
    workdir = localWorkDir, 
    samples = samples,
    files = files
)

df.tpm  <- list.tpm.fpkm$df.tpm
#df.fpkm <- list.tpm.fpkm$df.fpkm
rm(list.tpm.fpkm)

## Remove all 0 only rows ##
rSums <- rowSums(df.tpm[,2:ncol(df.tpm)])
df.tpm <- df.tpm[rSums > 0, ]

###############################################################################
## Do differential gene expression analysis                                  ##
###############################################################################
## Review raw counts file column names ##
col.names <- names(
    read.delim(
        count.data.file, 
        header=TRUE, 
        sep="\t", 
        stringsAsFactors = FALSE
    )
)

print(head(col.names))

## Add replicate column to dfDesign ##
dfDesign[["replicate"]] <- sapply(
    dfDesign$sample.id,
    function(x) unlist(strsplit(x, "_"))[3]
)

dfDesign[["FASTQ"]] <- dfDesign$sample.id

setwd(localWorkDir)

###############################################################################
## Retrieve count matrix (RSEM output)                                       ##
raw.counts.filt <- readAndPrepareCountMatrix(
    count.data.fn = paste0(raw.count.dir, "/",count.data.file),
    string.to.be.deleted.in.raw.counts.columns = paste0("X", gsub("/", ".",paste0(workdir, "RSEM/Ensembl/"))),
    df.design = dfDesign
)

## Done retrieving count matrix                                              ##
###############################################################################

###############################################################################
## Determine variation in the data set                                       ##

determineCoefficientOfVariation <- function(
    raw.counts.filt = 'raw.counts.filt'
    df.design = dfDesign
) {
    ###########################################################################
    ## Load required libraries                                               ##
    library(DESeq2)
    
    ## Done Loading libraries                                                ##
    ###########################################################################
    
    ###########################################################################
    ## Create required folders                                               ##
    
    if (length(grep("DESeq2$", list.files())) == 0){
        dir.create(DEseq2Dir)
    }
    
    setwd(DEseq2Dir)
    ## Done creating folders                                                 ##
    ###########################################################################
    
    ###########################################################################
    ## Format design file                                                    ##
    ## Ensure that replicate is set ##
    
    if (!batch.mode) {
        df.design[["replicate"]] = paste0("B", df.design$sample.id)
    }
    
    ## Done formating design file                                           ##
    ###########################################################################
    
    
    ###########################################################################
    ## Getting started with the PCA                                          ##
    ###########################################################################
    ## Create PCA plot for sample.groups                                     ##
    
    if (batch.mode){
        colData = unique(df.design[, c("sample.id", "sample.group","replicate")])
        rownames(colData) = colData[,"sample.id"]
        colData$sample.id <- NULL
        colnames(colData)[1] = "condition"
        colData$condition <- as.factor(colData$condition)
        colData$replicate <- as.factor(colData$replicate)
        
        
        dds <- DESeqDataSetFromMatrix(
            countData = raw.counts.filt,
            colData   = colData,
            design    = ~ replicate
        )
    } else {
        colData = unique(df.design[, c("sample.id", "sample.group")])
        rownames(colData) = colData[,"sample.id"]
        colData$sample.id <- NULL
        colnames(colData)[1] = "condition"
        colData$condition <- as.factor(colData$condition)
        
        dds <- DESeqDataSetFromMatrix(
            countData = raw.counts.filt,
            colData   = colData,
            design    = ~ condition
        )
    }
    
    
    dds <- DESeq(
        dds,
        test = "Wald",
        betaPrior = FALSE,
        parallel = parallel
    )
    
    ## Extract norm counts ##
    ## RSEM-generate-matrix produces a raw-count matrix
    normCounts = round(counts(dds, normalized=toBeNormalized))
    df.normCounts = data.frame(normCounts)
    #Remove all rows 0 counts for all samples from df.normCounts
    df.normCounts = df.normCounts[rowSums(df.normCounts)!=0,]
    
    df.data["CoVar"] <- 0
    
    ## Ignore low-intesity rows ##
    df.data[df.data$count_cut_off > 1,"CoVar"] <- apply(
        df.data[df.data$count_cut_off > 1, grep("^norm_counts_", names(df.data))],
        1,
        function(x) sd(x)/mean(x)
    )
    
    df.data[is.na(df.data)] <- 0
    df.data[df.data$CoVar == Inf, "CoVar"] <- max(df.data[df.data$CoVar < Inf ,"CoVar"])
    
    
    # Order from highest to lowest CoVar
    df.data <- df.data[order(df.data$CoVar, decreasing = TRUE),]
    df.data[["CoVarOrder"]] <- 1:nrow(df.data)
    
    plot(
        df.data$CoVarOrder,
        df.data$CoVar,
        type = "l"
    )
    abline(h=500, col="red")
    
    return(df.data)
    
}

#dfCorVar <- determineCoefficientOfVariation(
#    raw.counts.filt = raw.counts.filt
#)

## Done determine variation in the data set                                  ##
###############################################################################

###############################################################################
## Perform PCA, MV-analysis, and Clusterdendrogram                           ##
pcaDimensionsToInvestigate <- c(1:5)

df.pca <- do.PcaVarianceSampleClustering(
    raw.counts.filt = raw.counts.filt,             #count data filename
    DEseq2Dir = paste0(localWorkDir, "DESeq2"), #directory for results
    df.design = dfDesign,                      #df.design
    LRTcolTag = "LRT_",
    batch.mode = FALSE, 
    parallel = FALSE,
    timeseries = FALSE,
    Ntop4pca = 500,
    doPCA = TRUE,
    lmFitDim = pcaDimensionsToInvestigate ,
    primary.alignment.gene.id = primary.alignment.gene.id
)


dfPcaRes <- df.pca$dfPCAgeneRes
names(dfPcaRes) <- gsub("intercept.PCA", "intercept_PCA", names(dfPcaRes))
names(dfPcaRes) <- gsub("estimate.PCA", "contrast_P_PCA_estimatePCA", names(dfPcaRes))
names(dfPcaRes) <- gsub("padj.PCA", "contrast_P_padj_PCA", names(dfPcaRes))
names(dfPcaRes) <- gsub("rsquared.PCA", "r2.PCA", names(dfPcaRes))

dfPcaSamples <- df.pca$df.pca

## First, run DESeq2 in normal mode ##
dfResList <- do.differential.expression.analyis(
    raw.counts.filt = raw.counts.filt,             #count data filename
    DEseq2Dir = paste0(localWorkDir, "DESeq2"),    #directory for results
    df.design = dfDesign,                         #dfDesign
    gene.id = primary.alignment.gene.id,                           #primary gene id after alignment 
    batch.mode = FALSE, #if true, dfDesign needs to contain a 'replicate' column 
    timeseries = FALSE
)

df.summary <- dfResList$df.summary
names(df.summary) <- gsub("contrast_D_padj_LRTdataseries" , "contrast_D_padj_LRT_Disease_Status", names(df.summary) )
names(df.summary) <- gsub("contrast_G_padj_LRTsampleGroup" , "contrast_G_padj_Sample_Group", names(df.summary) )
df.summary$LRT_logFC <- NULL


###############################################################################
## Attach the estimate                                                       ##

dfPcaRes <- dfPcaRes[dfPcaRes[,primary.alignment.gene.id] %in% df.summary[,primary.alignment.gene.id],]
dfPcaRes <- dfPcaRes[c(primary.alignment.gene.id, names(dfPcaRes)[grep("contrast_P", names(dfPcaRes))])]

df.summary <- merge(
    df.summary, 
    dfPcaRes, 
    by.x = primary.alignment.gene.id,
    by.y = primary.alignment.gene.id,
    all =TRUE
)

df.summary[is.na(df.summary)] <- 0

## Done attaching                                                            ##
###############################################################################


###############################################################################
## Adding human disease datasets                                             ##

dfAdd <- import.db.table.from.db(
    dbname = prim.data.db, 
    dbtable = "p285_rna_seq_table",
    password = db.pwd
)

selVec <- c(names(dfAdd)[grep("hr_", names(dfAdd))])
selVec <- c("ENSG", selVec[grep("contrast_", selVec)])

dfAdd <- unique(
    dfAdd[,selVec]
)

dfAdd <- dfAdd[dfAdd$ENSG %in% df.summary$ENSG, ]

names(dfAdd) <- gsub("contrast_", "contrast_P", names(dfAdd))

df.summary <- merge(
    df.summary, 
    dfAdd, 
    by.x = "ENSG",
    by.y = "ENSG",
    all = TRUE
)
df.summary[is.na(df.summary)] <- 0

## Run DESeq2 in batch mode, if requiered##
#df.summary2 <- do.differential.expression.analyis(
#    raw.count.dir = paste0(localWorkDir, "RSEM"),
#    count.data.file = count.data.file,             #count data filename
#    DEseq2Dir = paste0(localWorkDir, "DESeq2"),    #directory for results
#    df.design = dfDesign,                         #dfDesign
#    gene.id = primary.alignment.gene.id,                           #primary gene id after alignment 
#    batch.mode = TRUE, #if true, dfDesign needs to contain a 'replicate' column 
#    string.to.be.deleted.in.raw.counts.columns = paste0("X", gsub("/", ".",paste0(workdir, "RSEM/Ensembl/")))
#)
#df.summary2 <- unique(df.summary2[,c("ENSG",names(df.summary2)[grep("contrast_2", names(df.summary2))])])
#names(df.summary2) <- gsub("contrast_2", "contrast_5", names(df.summary2))

#df.summary <- unique(merge(df.summary, df.summary2, by.x = "ENSG", by.y="ENSG", all = TRUE))
#df.summary[is.na(df.summary)] <- 0

## Add comparison 1 to full list ##


###################################
# Upload pca table to database    #
###################################
names(dfPcaSamples) <- gsub("[.]", "_", names(dfPcaSamples))

## Documentation ##
# For each sample_group_[sample_group_name] 
# 

upload.pca.table.to.db(
    df.pca = dfPcaSamples,
    host = host,
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
    preprocess <- as.numeric(df.summary[,padj[i]])
    
    if (length(grep("padj_LRT", padj[i])) > 0){
        preprocess <- as.numeric(df.summary[,padj[i]])
        minNum <- min(preprocess[preprocess != 0])
        preprocess[preprocess == 0] <- minNum
    } else {
        preprocess <- as.numeric(df.summary[,padj[i]])
    }
    
    temp <- -1*log10(preprocess)
    #temp[temp >= 50] = 50
    df.summary[,lg10p[i]] <- temp
}


###############################################################################
# Add TPM values to data table                                               #
###############################################################################
new.order <- names(df.summary)
new.order.tpm <- new.order[grep("norm_counts_", new.order)]
new.order.tpm <- gsub("norm_counts_", "", new.order.tpm)
new.order.tpm <- paste0(new.order.tpm, ".TPM")
new.order.tpm <- append("gene_id", new.order.tpm)
df.tpm <- df.tpm[,new.order.tpm]
names(df.tpm) <- gsub(".TPM", "",names(df.tpm))
names(df.tpm)[2:length(df.tpm)] <- paste0("TPM________",names(df.tpm)[2:length(df.tpm)])
names(df.tpm) <- gsub("gene_id", primary.alignment.gene.id,names(df.tpm))

if (!(exists("addFullTPMtable"))){
    addFullTPMtable <- FALSE
}

df.summary <- merge(
    df.summary, 
    df.tpm, 
    by.x = primary.alignment.gene.id, 
    by.y = primary.alignment.gene.id,
    all = addFullTPMtable
)

df.summary[is.na(df.summary)] <- 0

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

## Ensure all contrast columns are numeric ##
df.summary[,grep("contrast_", names(df.summary))] <- apply(
    df.summary[,grep("contrast_", names(df.summary))],
    2,
    as.numeric
)

df.summary[is.na(df.summary)] <- 0

## Select for heatmap ##
df.summary <- selectHeatmapGenes(
    dfData = df.summary,
    cutOff = 6,
    zeroOneCol = "logFC_cut_off",
    selCol = "contrast_G_lg10p_Sample_Group"
)

## Select for heatmap all genes with a TPM row sum of 2 or higher ##
# df.summary[["logFC_cut_off"]] <- 0
# df.summary[,"logFC_cut_off"] <- rowSums(df.summary[,grep("norm_counts_", names(df.summary))])
# nSamples <- length(unique(dfDesign$sample.id))
# df.summary[,"logFC_cut_off"] <- ifelse(df.summary$logFC_cut_off >= 5*nSamples, 1, 0)

## Select for heatmap: abs change of at least 0.5 in any contrast ##
database.table <- datatable.to.website.ptm(
    df.data = df.summary, 
    gene.id.column = primary.alignment.gene.id, 
    heatmap.genes = "",
    n.cluster.genes = 2000, 
    count.data = TRUE, 
    logFC.cut.off = 1, 
    gene.id.table = gene.id.table,
    add.uniprot.column = TRUE, 
    #use.logFC.columns.for.heatmap = FALSE,
    selector4heatmap.cols = "norm_counts",
    heatmap.preprocessing = "lg2.row.avg", # possible: "lg2", "lg2.row.avg", "none"
    hm.cut.off = 4,
    n.hm.cluster = 10,
    count.cut.off.filter = 0
)

###############################################################################
## Create CoVar plot for documentation                                       ##

dfCoVar <- unique(
    database.table[,c(primary.alignment.gene.id, "CoVar","CoVarOrder")]
)

dfCoVar <- dfCoVar[order(dfCoVar$CoVarOrder, decreasing = FALSE),]
dfCoVar <- dfCoVar[dfCoVar$CoVar != 0,]

dfCoVar50 <- dfCoVar[ (0.5*nrow(dfCoVar)):nrow(dfCoVar),]
fit <- lm(CoVar ~ CoVarOrder, data=dfCoVar50)

png("CoVarPlot.png", type="cairo")
    plot(dfCoVar$CoVarOrder, dfCoVar$CoVar, type="l")
    abline(v=500, col="red")
    abline(fit, col="grey")
    points(dfCoVar$CoVarOrder, dfCoVar$CoVar, type="l")
dev.off()

## Done CoVar for documentation                                              ##
###############################################################################

###############################################################################
## Create Excel output files                                                 ##
createAndFormatExcelOutputFiles(
    database.table = database.table,
    outputDir = outputDir,
    metaCoreCountFilter = 1
)

## Done creating Excel output files                                          ##
###############################################################################


###############################################################################
## Add additional plot columns from database                                 ##
dfAdd <- import.db.table.from.db(
    dbname = prim.data.db,
    dbtable = "p276_rna_seq_table",
    password = db.pwd
)

selVec <- c(
    "hgnc_symbol",
    names(dfAdd)[grep("contrast_", names(dfAdd))]
)

dfAdd <- unique(dfAdd[, selVec])

dfAdd <- dfAdd[dfAdd$hgnc_symbol %in% database.table$hgnc_symbol,]
names(dfAdd) <- gsub("contrast_", "contrast_J", names(dfAdd))

database.table <- merge(
    database.table,
    dfAdd, 
    by.x = "hgnc_symbol",
    by.y = "hgnc_symbol",
    all=TRUE
)

## Done adding additional plot columns                                       ##
###############################################################################

###############################################################################
# Upload to database                                                          #
###############################################################################
database.table$Gene_name <- NULL
database.table$Gene_type <- NULL
database.table$Gene_description <- NULL

cmd.vec <- upload.datatable.to.database(
    host = host, 
    user = db.user,
    password = db.pwd,
    prim.data.db = prim.data.db,
    dbTableName = rnaseqdbTableName,
    df.data = database.table[database.table$count_cut_off > 1, ],
    db.col.parameter.list = list(
        "VARCHAR(255) CHARACTER SET latin1 COLLATE latin1_swedish_ci" = c("Gene_name","gene_description"),
        "VARCHAR(50) CHARACTER SET latin1 COLLATE latin1_swedish_ci" = c("Gene_description","Gene_type","ENSG", "ENSMUSG", "hgnc_symbol", "mgi_symbol", "uniprot", "entrezgene","display_ptm", "^sequence_window", "p_site_env","for_GSEA_gene_chip","associated_gene_name", "gene_type"),
        "VARCHAR(1) CHARACTER SET latin1 COLLATE latin1_swedish_ci" = c("ppos", "amino_acid", "charge","known_site"),
        "BIGINT(8) NULL DEFAULT NULL" = c("row_names"),
        "INT(8) NULL DEFAULT NULL" = c("CoVarOrder","row_id", "cluster_order","cluster_id", "count_cut_off", "^position$", "raw_counts"),
        "DECIMAL(6,3) NULL DEFAULT NULL" = c( "contrast_P_PCA_estimatePCA","^CoVar$","NES", "logFC", "lg2_avg", "intensity", "^int", "iBAQ","^localization_prob$", "stat", "lg10p"),
        "DECIMAL(6,5) NULL DEFAULT NULL" = c("padj", "pvalue","^pep$", "p_value_PC1"),
        "DECIMAL(6,1) NULL DEFAULT NULL" = c("norm_counts")
    ),
    new.table = TRUE
)

killDbConnections()

  
###############################################################################
# Do GSEA                                                                     #
###############################################################################

database.table2 <- database.table

## Remove unnecessary columns, if needed ##
names(database.table2) <- gsub("J1_logFC_", "J1_AAAA_", names(database.table2))
names(database.table2) <- gsub("P1_logFC_", "P1_AAAA_", names(database.table2))
names(database.table2) <- gsub("P2_logFC_", "P2_AAAA_", names(database.table2))
names(database.table2) <- gsub("P3_logFC_", "P3_AAAA_", names(database.table2))
names(database.table2) <- gsub("P4_logFC_", "P4_AAAA_", names(database.table2))

## Create GSEA rank files ##
create.gsea.rnk.files(
    localWorkDir, 
    df.dataTable = database.table2,
    GSEA.colum.type = "_logFC_",
    gene.symbol.column.name = "hgnc_symbol"
)

## Remove last character from file ##
#truncate -s -2 file
#sed '$d' file # remove last line

## Remove last character from file ##
#truncate -s -2 file
#sed '$d' file # remove last line

## Function to create gmt file ##
tables <- c(
    "mysigdb_h_hallmarks"
)

dfRefGmt <- create.gmt.file.from.ref.data.table(
    host = 'www.biologic-db.org',
    dbname = "reference_categories_db_new",
    dataTable = tables,
    pwd = db.pwd,
    user="boeings",
    gene.id.column = "hgnc_symbol"
)

###############################################################################
## Save gmt file                                                             ##
#"/camp/stp/babs/working/boeings/Projects/reference_data/GSEA.gmt.files/20160508.rna.seq.txn.analysis.gmt"

localGmtDir <- paste0(
    localWorkDir,
    "GSEA/"
)

gmtDir<- paste0(
    workdir,
    "GSEA/"
)

gmtFileName <- paste0(
    project_id,
    ".",
    "projectGmtFile.gmt"
)

dfRefGmt <- dfRefGmt[!(duplicated(dfRefGmt[,1])),]

write.table(
    dfRefGmt,
    paste0(localGmtDir, gmtFileName),
    col.names = FALSE,
    row.names = FALSE,
    sep="\t"
)

contrasts <- names(database.table2)
contrasts <- contrasts[grep("logFC", contrasts)]
contrasts <- contrasts[contrasts != "logFC_cut_off"]
contrasts <- contrasts[contrasts != "contrast_J1_logFC_Skin_Card14KI_vs_Skin_WT"]
contrasts <- contrasts[contrasts != "contrast_P1_logFC_CARD14E138A_0hr_vs_CARD14WT_0hr"]
contrasts <- contrasts[contrasts != "contrast_P2_logFC_CARD14E138A_3hr_vs_CARD14WT_3hr"]
contrasts <- contrasts[contrasts != "contrast_P3_logFC_CARD14E138A_6hr_vs_CARD14WT_6hr"]
contrasts <- contrasts[contrasts != "contrast_P4_logFC_CARD14E138A_9hr_vs_CARD14WT_9hr"]



GSEAfn <- paste0(
    localWorkDir,
    "/GSEA/GSEAcommands.sh"
)
sink(GSEAfn)

cat("module load Java/1.8.0_131");cat("\n");cat("\n")
for (i in 1:length(contrasts)){
    gmtFile <- paste0(gmtDir, gmtFileName)
    contrastNo <- unlist(strsplit(contrasts[i], "_"))[2]
    nTopPlots <- 50
    GSEAdir <- paste0(workdir, "GSEA")
    rnkFile <- paste0(GSEAdir, "/",contrasts[i],".rnk")
    
    gseaCMD <- paste0(
        "java -Xmx2512m -cp /camp/stp/babs/working/boeings/Projects/software/gsea-3.0.jar xtools.gsea.GseaPreranked -gmx ",
        gmtFile,
        " -rnk ",
        rnkFile,
        " -rpt_label ",
        "contrast_",
        contrastNo,
        "_rnaSeqTxnTest",
        " -out ",
        GSEAdir,
        " -collapse false -mode Max_probe -norm meandiv -nperm 1000 -scoring_scheme classic -include_only_symbols true -make_sets true -plot_top_x ",
        nTopPlots,
        " -rnd_seed timestamp -set_max 2500 -set_min 10 -zip_report false -gui false"
    )
    cat(gseaCMD);cat("\n");cat("\n");
    
    
}
sink()

## Convert ## 
print(
    paste0(
        "Remove end of line characters: tr -d '\r' <GSEAcommands.sh> ",
        "conv.GSEAcommands.sh"
    )
)

#srun --mem=64G --partition=hmem --x11 xterm


# .rnk files need to be opened and saved in Excel once before being forwarded to the GSEA JAVA application.
#GSEA parameters
#Collapse dataset to gene symbols = FALSE
#Number of permutations = 1000
#C2: Min. size 15, max. size = 500

## load reference gmt files 
## 20160508.rna.seq.txn.analysis.gmt
## pathways.and.complexesfor.rna.seq.gsea.gmt

#Got to GSEA to calculate category enrichments
#Result name convention contrast_X_C2
#Write result folders to .../GSEA

#############################
# Run GSEA in JAVA interface
#############################

###############################################################################
## Optional module hypergeometric test                                       ##

## Done: Optional module hypergeometric test                                 ##
###############################################################################


###############################################################################
## Stopping point                                                            ##
## Save workspace ##
workspaceFN <- paste0(
    localWorkDir,
    project_id,
    ".RData"
)

save.image(workspaceFN)


## Load workspace ##
workspaceFN <- paste0(
    localWorkDir,
    project_id,
    ".RData"
)

load(workspaceFN)

##                                                                           ##
###############################################################################

###############################################################################
## Process GSEA outputs                                                      ##

# Subfunction to select relevant tables
library(RMySQL)
dbDB = dbConnect(
    drv = RMySQL::MySQL(), 
    user = db.user, 
    password = db.pwd, 
    dbname = "reference_categories_db_new",
    host = host
) 

#This will move to the header
availableTables = dbGetQuery(dbDB, "SHOW TABLES in reference_categories_db_new")
availableTables <- as.vector(availableTables$Tables_in_reference_categories_db_new)
dbDisconnect(dbDB)

## Check if all reference tables are present in the db ##
print(
    paste0(
        sum(!(tables %in% availableTables)),
        " tables NOT present in the reference database."
    )
)



#temporary change
#tables <- tables[c(1,2,3,5,6)]
#tables <- append(tables, "vt_lab_categories")
# End of parts to be moved to the header

gsea.cat.lines <- create.GSEA.table(
    GSEADir = paste0(localWorkDir, "GSEA"),
    logFC.column.name = "logFC",
    host = host,
    refdbname= "reference_categories_db_new",
    refDBTableName = enriched.categories.dbTableName,
    db.user = db.user,
    db.password = db.pwd, 
    tables = tables,
    df.dataTable = database.table2,
    outputDir = outputDir
)

## Done processing GSEA outputs                                              ##
###############################################################################


###############################################################################
## Optional Module: Correlation analysis                                     ##

## If required, run the module 
# templateRNAseqAnalysis_pluginModule_correlationAnalysis.R
# here.

## DoneOptional Module: Correlation analysis                                 ##
###############################################################################


###############################################################################
## Select tables for reference listing                                       ##


## List all available reference tables ##
list.db.tables.in.db(
    dbname = "reference_categories_db_new",
    password = db.pwd
)

labCatName <-  paste0(labname, " Lab")

referenceTableList <- list(
    labCatName = lab.categories.table,
    "Hallmark Signatures" = "mysigdb_h_hallmarks",
    "Pathways" = "mysigdb_c2_1329_canonical_pathways",
    "GO-BP" = "mysigdb_c5_BP",                               
    "GO-MF" = "mysigdb_c5_MF",
    "TF Motifs" = "TRANSFAC_and_JASPAR_PWMs",
    "Protein Complexes" = "networkcategories"
)

## Define relevant genes for selection ##
relevant.genes <- as.vector(
    unique(
        database.table[database.table$cluster_order, gene.id.column]
    )
)

length(relevant.genes)

for (i in 1:length(referenceTableList)) {
    df.ref <- import.db.table.from.db(
        dbname = ref.cat.db,
        dbtable = referenceTableList[[i]],
        password = db.pwd
    )
    
    ## Remove temp categories ##
    temPos <- grep("temp_", df.ref$cat_type)
    
    if (length(temPos) > 0) {
        df.ref <- df.ref[-temPos,]
    }
    
    df.temp <- add2labCatSelectionDBtable(
        df.ref = df.ref,
        cat_group_name = names(referenceTableList)[i],
        reference.gene.vector = relevant.genes,
        ref.gene.vec.id = gene.id.column,
        cat_views = NA
    )
    
    if (i == 1) {
        df.db.table <- df.temp
    } else {
        df.db.table <- rbind(df.temp, df.db.table)
    }
}

###############################################################################
## Filter cat dataset                                                        ##
###############################################################################
## Remove all datasets with less than 5 proteins matched
catIdGroupsExemptFromFiltering <- c(
    lab.categories.table
)

dfUnfiltered <- df.db.table[as.vector(unlist(sapply(catIdGroupsExemptFromFiltering, function(x) grep(x, df.db.table$cat_id)))),]    

dim(df.db.table)
df.db.table <- df.db.table[df.db.table$cat_count >= 5,]
dim(df.db.table)

## Add unfiltered categories ##
if (exists("dfUnfiltered") & (nrow(dfUnfiltered) > 0)){
    df.db.table <- unique(
        rbind(
            dfUnfiltered,
            df.db.table
        )
    )
}

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

if (exists("df.ref") & nrow(df.ref) > 0){
    df.temp <- add2labCatSelectionDBtable(
        df.ref = df.ref,
        cat_group_name = "This project",
        reference.gene.vector = reference.gene.vector,
        ref.gene.vec.id = gene.id.column,
        cat_views = NA
    )
    
    df.db.table <- rbind(df.temp, df.db.table)
}
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
    ),
    new.table = TRUE
)

##  End creating reference category selection for this project               ##
###############################################################################
    
###############################################################################
# Create timecourse specifications                                            #
###############################################################################
# Add dataseries_color to dfDesign

#timecourse.cat.lines <- create.timecourse.cat.lines(
#     dfDesign
# )


###############################################################################
# Create microwebsite                                                         #
###############################################################################

## Color samples by sample group ##
# Get sample order from database.table
## oder samples according to PC1
# df.pca <- df.pca[order(df.pca$pca1, decreasing = FALSE),]
# sample.order <- paste0("norm_counts_", df.pca$sample_id)

sample.order <- sort(names(database.table)[grep("norm_counts_", names(database.table))])


## Sort database.table ##
orderVec <- names(database.table)
lg2Order <- gsub("norm_counts_", "lg2_avg_", sample.order)

orderVec <- orderVec[!(orderVec %in% sample.order)]
orderVec <- orderVec[!(orderVec %in% lg2Order)]

orderVec <- c(
    orderVec,
    sample.order,
    lg2Order
)

database.table <- database.table[,orderVec]

## Get relevant colors ##
sample.colors <- substr(sample.order, 1, nchar(sample.order)-2)
sample.colors <- gsub("_$", "", sample.colors)

color.groups <- unique(sample.colors)

#dfDesign[["sample_color"]] <- "orange"
nSampleGroups <- length(unique(color.groups))

library(RColorBrewer)
selcol <- colorRampPalette(brewer.pal(12,"Set3"))

color.sel = selcol(nSampleGroups)

for (i in 1:length(color.groups)){
    sample.colors[grep(color.groups[i], sample.colors)] <- color.sel[i]
}


if (!exists("gsea.cat.lines")){
    gsea.cat.lines <- ""
    downloadCatEnrichmentFNxlsx <- ""
} else {
    downloadCatEnrichmentFNxlsx <- paste0("outputs/", project.code, ".enriched.categories.txt")
}

if (!exists("timecourse.cat.lines")){
    timecourse.cat.lines <- NA
}

#library(SBwebtools)
setwd(localWorkDir)

#library(SBwebtools)
webSiteDir <- localWorkDir

create.website.parameters(
    df.data = database.table,
    gene.id.column = gene.id.column,
    ptm.colum = "",
    lab_id = lab_id , 
    user_ids = user_ids, 
    project_id = project_id, 
    download_result_table = paste0("outputs/", project.code, ".result.table.txt"),
    download_cat_enrichment_table = downloadCatEnrichmentFNxlsx,
    database = prim.data.db,
    reference_categories_db = "reference_categories_db_new",
    labname = labname,
    rnaseqdbTableName,
    lab.categories.table = lab.categories.table,
    sample.order = sample.order, #set to "" to go with default sorting
    count.sample.colors = sample.colors,
    count.table.headline = count.table.headline,
    count.table.sidelabel = count.table.sidelabel,
    webSiteDir = webSiteDir,
    heamap.headline.text = heamap.headline.text,
    upper_heatmap_limit = 3, 
    lower_heatmap_limit = -3,
    slider.selection.name = "logFC",
    presentation.file = paste0("outputs/", project.code, ".project.presentation.pptx"),
    number_of_slides = 7,
    default.sequence = "",
    use.logFC.columns.for.heatmap = FALSE,
    peptide.view.link = "",
    create.2d.scatterplot.button = TRUE,
    low_highlight = -1,
    high_highlight =1,
    display.qc = TRUE,
    display.pca = TRUE,
    display.report = TRUE,
    pca.table.name = PCAdbTableName, 
    gsea.cat.lines = gsea.cat.lines,
    timecourse.cat.lines = timecourse.cat.lines,
    venn.slider.selector.strings = c(
        "_logFC", 
        "contrast_1_padj",
        "contrast_2_padj",
        "contrast_3_padj",
        "contrast_J1_padj",
        "contrast_P1_padj",
        "contrast_P2_padj",
        "contrast_P3_padj",
        "contrast_P4_padj",
        "contrast_D_lg10p_LRT",
        "contrast_G_lg10p_LRT"
        
    ),
    plot.selection.strings = c(
        "_logFC", 
        "_PCA_",
        "_lg10p",
        "contrast_P_PCA_estimatePCA",
        "contrast_P_padj_PCA"
    ), #NA, #strings to grep from col names for plot display
    plate.view.db.table = NA,
    plate.view.column.vec = NA,
    cat.seletion.table.vec = c(cat.ref.db, cat.ref.db.table)
)

###############################################################################
## Assemble microwebsite                                                     ##

mWVec <- c(
    "cp -r ../GSEA/enrichment_plots/ .",
    "cp -r ../RSEM/Ensembl/multiqc_report.html .",
    "cp -r /camp/stp/babs/working/boeings/Stefan/protocol_files/github/boeings/biologic/GD_bulkRNA_seq/* .",
    paste0("mv outputs/p999.project.presentation.pptx outputs/", project.code, ".project.presentation.pptx"),
    "cp ../../outputs/*xl* ./outputs/",
    "cp ../DESeq2/MA* ./outputs/",
    "cp ../DESeq2/PCA* ./outputs/",
    "cp ../heatmap.col.* ./outputs/",
    "cp ../DESeq2/covariance.plot.* ./report_figures"
)

setwd(paste0(localWorkDir, "/", project_id))
sink("assembleMicrowebsite.sh")
for (i in 1:length(mWVec)){
    cat(mWVec[i]);cat("\n");
}
sink()

## Done assembling microwebsite                                              ##
###############################################################################

###############################################################################
## Deploy next generation website                                            ##
createSettingsFile(
    df.data = database.table,
    primDataTable = rnaseqdbTableName,
    sample.order = sample.order, #set to "" to go with default sorting
    heatmapSampleOrder = "",
    sample.names = "", # default is sample.order
    count.sample.colors = sample.colors,
    ptm.colum = "",
    count.table.headline = count.table.headline,
    count.table.sidelabel = count.table.sidelabel,
    venn.slider.selector.strings = c(
        "_logFC", 
        "_PCA_",
        "contrast_1_padj",
        "contrast_2_padj",
        "contrast_3_padj",
        "contrast_J1_padj",
        "contrast_P1_padj",
        "contrast_P2_padj",
        "contrast_P3_padj",
        "contrast_P4_padj",
        "contrast_D_lg10p_LRT",
        "contrast_G_lg10p_LRT"
        
    ),
    plot.selection.strings = c(
        "_logFC", 
        "_lg10p"
    ),
    webSiteDir = paste0(hpc.mount, "Stefan/protocol_files/github/biologic/src/experiments/"),
    upper_heatmap_limit = 3, 
    lower_heatmap_limit = -3,
    heamap.headline.text = heamap.headline.text,
    project_id = project_id
)

## Done deploying next generation website                                    ##
###############################################################################


###########################################
# Do multiQC and add to QC section 
###########################################
# run 
# >> module use /camp/stp/babs/working/software/modules/all
# >> module load multiqc/0.9-2016b-Python-2.7.12

## Copy mga.tsv file into the working directory ##
# print(
#     paste0(
#         "cp ",
#         dataFolder,
#         "/mga.tsv ",
#         workdir
#     )
# )
# 
# # run in top project directory:
# # >> multiqc .
# # multiqc /camp/stp/babs/working/boeings/Projects/103_VTL_ES_RNA_seq_BAFF_timecourse_hs/workdir/RSEM/Ensembl /camp/stp/babs/working/boeings/Projects/103_VTL_ES_RNA_seq_BAFF_timecourse_hs/workdir/RSEM/Ensembl/RNAseQC /camp/stp/babs/working/boeings/Projects/103_VTL_ES_RNA_seq_BAFF_timecourse_hs/workdir/logs /camp/stp/babs/working/boeings/Projects/103_VTL_ES_RNA_seq_BAFF_timecourse_hs/workdir/FASTQC_stranded
# print("module purge; module use /camp/stp/babs/working/software/modules/all; module load multiqc/1.3-2016b-Python-2.7.12")
# print("multiqc --config /camp/stp/babs/working/escudem/bin/multiqc_config.yaml .")
# 
# print(paste0("multiqc . -o ", webSiteDir))
# print(paste0("multiqc . -o ", workdir, "outputs"))
# 
# ## Transfer multiqc_report.html to outputs ##
# 
# ## list directories to include in multiqc search
# #paste0(workdir, "logs")
# paste0(workdir, "RSEM/Ensembl")
# paste0(workdir, "RSEM/Ensembl/RNAseQC")
# paste0(workdir, "FASTQC")


#[61] knitr_1.15.1     


###############################################################################
## Create report for this analysis                                           ##

## Use templateReportRNASeq.html and customize to this analysi               ##
###############################################################################

listExistingProjects(
    user     = "boeings",
    password = db.pwd,
    host     = "www.biologic-db.org",
    dbname = "reference_categories_db_new",
    dbtable = "project_description_table"
)

###############################################################################
## Create new project if necessary                                           ##

createNewProject(
    dbname = "reference_categories_db_new",
    dbtable = "project_description_table",
    password = db.pwd,
    project_name = "Mrtf_Srf_regulated_cytoskeletal_dynamics",
    project_lab = ";Treisman;",
    project_description = "Characterization of Mrtf and Srf transcription programs in the regulation of cytoskeletal dynamics."
)

## Done creating new project                                                 ##
###############################################################################

###############################################################################
## List project in projects table                                            ##

addProject2ProjectTable(
   dbname = "reference_categories_db_new",
   dbtable = "project_db_table",
   password = db.pwd,
   experiment_id= project_id,
   experiment_question = "What are the potential target genes of YAP and TAZ in palate development?",
   experiment_description = "<b>RNA Seq analysis</b>. Examination of Yap Taz cartilage knockout animals at one day prior to birth has demonstrated the presence of a number of cleft palate phenotypes, including a failure of palatal shelf elevation. This process occurs at approximately",
   experiment_owner = ";Hannah Vanyai;",
   experiment_lab = paste0(";",labname,";"),
   experiment_viewers = paste0(";", paste(user_ids, collapse = ";"),";"),
   experiment_project = ";YAP_TAZ_in_Palate_Development;",
   experiment_type = experiment.type,
   experiment_details = "Experimental details will follow here",
   experiment_x_coordinate = ";siCtrl;siTarget;",
   experiment_x_coordinate_unit = "condition",
   experiment_link = paste0(
       "<a href='../",
       project_id,
       "/report.php' class='btn btn-success btn-lg' role='button'>Primary Data &raquo;</a>"
   ),
   experiment_title = ""
)
## Done listing project                                                      ##
###############################################################################

##########
## Done ##
##########

###############################################################################
## Assemble documentation rmd script                                         ##

##                                                                           ##
###############################################################################


## End of template                                                           ##
###############################################################################
