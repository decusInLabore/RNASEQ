###############################################################################
# Sahai lab                                                                   #
# Experiment:                                                                 #
# Primary data anlalysis                                                      #
# Ensembl release-89                                                          #
###############################################################################

## The sample sample.ids.txt file needs to have the following columns:
# sampleID	sample.name	sample.id


rm(list = ls())

## Set database password ##
if (exists("db.pwd")){
    print("Database password is set.")
} else {
    ## Set database password ##
    library(rstudioapi)
    db.pwd <- rstudioapi::askForPassword("Please enter database password")
}

user_ids <- c(
    "project", 
    "ag.lab.all",
    "andrew.bailey",
    "demo1",
    "demo2",
    "demo3",
    "demo",
    "song.chng",
    "tiffany.heanue",
    "vp.lab.all",
    "franze.progatzky"
)

## Check if S4LvariableListing has been saved before ##

S4LvariableListing <- list(
    "folder" = "goulda/andrew.bailey/268_AGL_AB_RNASeq_ISRIB_NGF_DLKi_sensitivity_GSE95672",
    "sra.id.vector" = "SRP101348",
    "gse.id.vector" = "GSE95672",
    "experiment.type" = "rna_seq",    # any
    "species" = "mm",                 # is either mm or hs 
    "project_id" = "agl268",
    "labname" = "Gould",
    "db.user" = "boeings",
    "host" = "www.biologic-db.org",
    "timecourse.units" = "",
    "count.table.headline" = "TPM Values for all Samples",
    "count.table.sidelabel" = "TPM",
    "heamap.headline.text" = "Heatmap: Row-averaged TPM",
    "paired.end" = FALSE,
    "stranded" = TRUE,
    "release" = "release-89", #release-86, release-89
    "loadR" = "module purge;module load R/3.3.1-foss-2016b-bioc-3.3-libX11-1.6.3;R;",
    "ref.cat.db" = "reference_categories_db_new",
    "ref.cat.db.table" = "reference_categories_db_new",
    "read.length" =  "50bp",
    "webSiteDir" = "",
    "pathToSeqStorageFolder" = "/camp/stp/sequencing/inputs/instruments/fastq/180914_K00102_0256_BHWFH3BBXX/fastq/RN18087/",
    "addFullTPMtable" = FALSE
)
    
    
#library(SBwebtools)
## on camp ##
#source("/camp/stp/babs/working/boeings/Stefan/protocol_files/github/boeings/packages/20171107.package.SBwebtools.r")
source("Y:/working/boeings/Stefan/protocol_files/github/boeings/packages/20171107.package.SBwebtools.r")


S4LvariableListing <- setPathsAndVariablesEnsemblHsMmRelease89(
    S4LvariableListing
)

if (!dir.exists(outputDir)){
    dir.create(outputDir)
}

###############################################################################
## Create empty documentation vector                                         ##
## Create vector to collect all shell script lines ##
shellScriptVector <- as.vector(NULL, mode = "character")


## Done creating documentation vector                                        ##
###############################################################################

###############################################################################
## Optional module:Download FASTQ files from SRA                             ##
###############################################################################
## Download FASTQ files from SRA                                             ##
#library(SBwebtools)

## Finding SRA id ##
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE93313


## Create shell script to download FASTQfiles ##
# For this, create EMPTY FASTQ_files folder
setwd(localWorkDir)
createSRAdownloadScript(
    sra.id.vector = sra.id.vector,
    module.load.cmd = paste(
        "module use /camp/stp/babs/working/software/modules/all;",
        "module load sratoolkit/2.8.2-1"
    ),
    fastqDir = fastqDir
)

##  Done downloading FASTQ files from SRA                                    ##

###############################################################################
## Create design file                                                        ##

## Retrieve sample annotation
annotationDownloadURL <- paste0(
    "http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?",
    "save=efetch&db=sra&rettype=runinfo&term= ", 
    sra.id.vector
)

annotationFileName <- paste0(
    sra.id.vector,
    "_annotation.csv"
)

library(RCurl)

setwd(localDataDir)

download.file(
    annotationDownloadURL, 
    annotationFileName
)

dfAnno <- read.csv(
    annotationFileName
)

## Select RNA-Seq ##
dfAnno <- dfAnno[dfAnno$LibraryStrategy == "RNA-Seq",]

dfAnno <- unique(
    dfAnno[,c(
        "SRAStudy", 
        "Sample",
        "Run", 
        "avgLength", 
        "Experiment", 
        "LibraryLayout"
    )
    ]
)

## Add extensions ##
if (paired.end){
    dfAnnoR2 <- dfAnno
    dfAnno[["srrFASTQname"]] <- paste0(
        fastqDir,
        dfAnno$Run,
        "_pass_1.fastq.gz"
    )
    
    dfAnnoR2[["srrFASTQname"]] <- paste0(
        fastqDir,
        dfAnno$Run,
        "_pass_2.fastq.gz"
    )
    
    dfAnno[["srxFASTQname"]]<- paste0(
        fastqDir,
        dfAnno$Experiment,
        "_pass_1.fastq.gz"
    )
    
    dfAnnoR2[["srxFASTQname"]]<- paste0(
        fastqDir,
        dfAnno$Experiment,
        "_pass_2.fastq.gz"
    )
    
    dfAnno <- rbind(
        dfAnno, 
        dfAnnoR2
    )
    
    rm(dfAnnoR2)
    
} else {
    dfAnno[["srrFASTQname"]] <- paste0(
        fastqDir,
        dfAnno$Run,
        "_pass_1.fastq.gz"
    )
    
    dfAnno[["srxFASTQname"]]<- paste0(
        fastqDir,
        dfAnno$Experiment,
        "_pass_1.fastq.gz"
    )
}

setwd(localDataDir)

## If necessary, concatenate SRR files into SRX samples ##

## Helper function ##
createConcatenateFASTQfilesShellScript <- function(
    dfConCat = "dataframe with SRR inut and SRX output file names",
    paired.end = FALSE
){
    ## Determine need to concatenate ##
    srxVec <- as.vector(unique(dfConCat$srxFASTQname))
    srrVec <- as.vector(unique(dfConCat$srrFASTQname))
    
    if (paired.end) {
        cutOffLength <- 2 * length(srxVec)
    } else {
        cutOffLength <- 1 * length(srxVec)
    }
    
    if (length(srrVec) > cutOffLength){
        concatenate = TRUE
    } else {
        concatenate = FALSE
    }
    
    SRXvec <- as.vector(
        unique(
            dfConCat$srxFASTQname
        )
    )
    
    ## Create shell script ##
    sink("concatFASTQfiles.sh")
    cat('#!/bin/sh'); cat('\n');
    
    for (i in 1:length(SRXvec)){
        catFiles <- as.vector(
            dfConCat[dfConCat$srxFASTQname == SRXvec[i], "srrFASTQname"]
        )
        
        
        if (concatenate){
            catCMD <- paste0(
                "cat ",
                paste(catFiles, collapse = " "),
                " > ",
                SRXvec[i]
            )
        } else {
            catCMD <- paste0(
                "mv ",
                catFiles[1],
                " ",
                SRXvec[i]
            )
        }
        
        cat(catCMD); cat('\n');cat('\n');
    }
    
    sink()
}

## End helper function

SRR <- as.vector(unique(dfAnno$Run))
SRX <- as.vector(unique(dfAnno$Experiment))

dfConCat <- unique(dfAnno[,c("srrFASTQname", "srxFASTQname")])

setwd(localWorkDir)
createConcatenateFASTQfilesShellScript(
    dfConCat = dfConCat,
    paired.end = TRUE
)

print("Remove end of line characters: tr -d '\r' <concatFASTQfiles.sh> conv.concatFASTQfiles.sh")


## Get SRA file annotion ##
library(GEOquery)
gse <- getGEO(gse.id.vector, GSEMatrix = TRUE)
dfGSE <- pData(gse[[1]])

## In this case dfGSE can not be loaded automatically. Load manually the series_matrix file ##
seriesMatrix <- read.delim(
    paste0(
        localDataDir,
        "series_matrix.txt"
    ),
    header = FALSE,
    sep="\t",
    stringsAsFactors = FALSE
)

seriesMatrix <- data.frame(t(seriesMatrix))

names(seriesMatrix) <- t(seriesMatrix[1,])
dfGSE <- seriesMatrix[-1,]
names(dfGSE)[34] <- "relation.1"



selVec <- names(dfGSE)[names(dfGSE) %in% names(seriesMatrix)]



dfGSE <- dfGSE[,c("relation.1", "title")]
dfGSE[["srxID"]] <- sapply(
    dfGSE$relation.1,
    function(x) unlist(
        strsplit(as.character(x), "=")
    )[2]
)

dfGSE <- dfGSE[dfGSE$srxID %in% dfAnno$Experiment,]

dfGSE[["sample.id"]] <- dfGSE$title
dfGSE$sample.id <- gsub("SAM24327525 anti.ngf.dmso", "antiNGF_DMSO_1", dfGSE$sample.id)
dfGSE$sample.id <- gsub("SAM24327526 anti.ngf.dmso", "antiNGF_DMSO_2", dfGSE$sample.id)
dfGSE$sample.id <- gsub("SAM24327527 anti.ngf.dmso", "antiNGF_DMSO_3", dfGSE$sample.id)
dfGSE$sample.id <- gsub("SAM24327528 anti.ngf.dmso", "antiNGF_DMSO_4", dfGSE$sample.id)
dfGSE$sample.id <- gsub("SAM24327529 anti.ngf.dmso", "antiNGF_DMSO_5", dfGSE$sample.id)


dfGSE$sample.id <- gsub("SAM24327534 anti.ngf.DLKi", "antiNGF_DLKi_1", dfGSE$sample.id)
dfGSE$sample.id <- gsub("SAM24327535 anti.ngf.DLKi", "antiNGF_DLKi_2", dfGSE$sample.id)
dfGSE$sample.id <- gsub("SAM24327536 anti.ngf.DLKi", "antiNGF_DLKi_3", dfGSE$sample.id)
dfGSE$sample.id <- gsub("SAM24327537 anti.ngf.DLKi", "antiNGF_DLKi_4", dfGSE$sample.id)

dfGSE$sample.id <- gsub("SAM24327520 ngf.dmso", "NGF_DMSO_1", dfGSE$sample.id)
dfGSE$sample.id <- gsub("SAM24327521 ngf.dmso", "NGF_DMSO_2", dfGSE$sample.id)
dfGSE$sample.id <- gsub("SAM24327522 ngf.dmso", "NGF_DMSO_3", dfGSE$sample.id)
dfGSE$sample.id <- gsub("SAM24327523 ngf.dmso", "NGF_DMSO_4", dfGSE$sample.id)
dfGSE$sample.id <- gsub("SAM24327524 ngf.dmso", "NGF_DMSO_5", dfGSE$sample.id)

dfGSE$sample.id <- gsub("SAM24327530 anti.ngf.isrib", "antiNGF_ISRIB_1", dfGSE$sample.id)
dfGSE$sample.id <- gsub("SAM24327531 anti.ngf.isrib", "antiNGF_ISRIB_2", dfGSE$sample.id)
dfGSE$sample.id <- gsub("SAM24327532 anti.ngf.isrib", "antiNGF_ISRIB_3", dfGSE$sample.id)
dfGSE$sample.id <- gsub("SAM24327533 anti.ngf.isrib", "antiNGF_ISRIB_4", dfGSE$sample.id)


dfGSE <- dfGSE[,c("srxID", "sample.id")]
dfDesign <- unique(dfAnno[,c("srxFASTQname", "Experiment", "avgLength")])

dfDesign <- merge(dfDesign, dfGSE, by.x = "Experiment", by.y = "srxID")

## helper function ##
completeDesignBasedOnSampleID <- function(
    dfBasedesign,
    fastqDir = ""
){
    ## Create sample group ##
    dfBasedesign[["sample.group"]] <- as.vector(
        sapply(
            dfBasedesign$sample.id,
            function(x) paste(
                unlist(
                    strsplit(x, "_")
                )[c(1:2)],
                collapse = "_"
            )
        )
    )
    
    ## Create dataseries ##
    dfBasedesign[["dataseries"]] <- as.vector(sapply(
        dfBasedesign$sample.group,
        function(x) unlist(
            strsplit(x, "_") 
        )[1]
    ))
    
    ## Create dataseries color ##
    library(RColorBrewer)
    
    nCol <- length(unique(dfBasedesign$dataseries))
    
    selcol <- colorRampPalette(brewer.pal(9,"YlOrBr"))
    groupCols <- selcol(nCol)
    
    dfColor <- data.frame(unique(dfBasedesign$dataseries), groupCols)
    names(dfColor) <- c("dataseries", "dataseries_color")
    
    dfDesign <- merge(dfBasedesign, dfColor, by.x = "dataseries", by.y = "dataseries")
    
    ## Set original.NGS ##
    names(dfDesign) <- gsub("srxFASTQname", "original.NGS", names(dfDesign))
    
    ## Set NGS ##
    dfDesign[["NGS"]] <- paste0(
        fastqDir,
        dfDesign$sample.id,
        "_R1.fastq.gz"
    )
    
    return(dfDesign)
    
}

## end helper function

## Add dataseries and sample group columns based on sample.id ##
dfDesign <- completeDesignBasedOnSampleID(
    dfBasedesign = dfDesign,
    fastqDir = fastqDir
)

dfDesign[grep("pass_2", dfDesign$original.NGS), "NGS"] <- gsub(
    "_R1",
    "_R2",
    dfDesign[grep("pass_2", dfDesign$original.NGS), "NGS"]
)

## Add timepoint column ##
dfDesign[["timepoint"]] <- 4.5



## Set NGS ##

## Write fastq to file to add comparisons comp_x

dfDesign[["comp_1"]] <- ""
dfDesign[dfDesign$sample.group == "antiNGF_DMSO","comp_1"] <- "1_antiNGF_DMSO"
dfDesign[dfDesign$sample.group == "NGF_DMSO","comp_1"] <- "2_NGF_DMSO"

dfDesign[["comp_2"]] <- ""
dfDesign[dfDesign$sample.group == "antiNGF_ISRIB","comp_2"] <- "1_antiNGF_ISRIB"
dfDesign[dfDesign$sample.group == "antiNGF_DMSO","comp_2"] <- "2_antiNGF_DMSO"

dfDesign[["comp_3"]] <- ""
dfDesign[dfDesign$sample.group == "antiNGF_DLKi","comp_3"] <- "1_antiNGF_DLKi"
dfDesign[dfDesign$sample.group == "antiNGF_DMSO","comp_3"] <- "2_antiNGF_DMSO"



setwd(localWorkDir)
#write.table(dfDesign, design.file, row.names = FALSE, sep = "\t")


## Done creating design file                                                 ##
###############################################################################




###############################################################################
## End: Optional module:Download FASTQ files from SRA                        ##
###############################################################################

## If necessary, concatenate SRR files into SRX samples ##




## set the timepoint parameter according to condition:
## Reorder df design ##
# order by sample group - sample name
dfDesign <- dfDesign[order(
    dfDesign$dataseries, 
   dfDesign$sample.group, 
    dfDesign$sample.id),
]


setwd(localWorkDir)
#write.table(dfDesign, design.file, row.names = FALSE, sep = "\t")


## Done creating design file                                                 ##
###############################################################################


###############################################################################
## Load   design file                                                        ##

setwd(localWorkDir)

dfDesign <- read.delim(
    design.file, 
    header = TRUE, 
    sep = "\t",
    stringsAsFactors = FALSE
)

## Reorder df design ##
# order by sample group - sample name
dfDesign <- dfDesign[order(
    dfDesign$dataseries, 
    dfDesign$timepoint, 
    dfDesign$sample.group, 
    dfDesign$sample.id),
]

print("Design file ordered correctly? If in doubt, check!")
## Done creating/loading design file                                         ##
###############################################################################

## Create vector to collect all shell script lines ##
#shellScriptVector <- as.vector(NULL, mode = "character")

###############################################################################
## Create shell script to rename files                                       ##
tempShellScriptVector <- as.vector(NULL, mode = "character")

tempShellScriptVector <- c(
    tempShellScriptVector,
    "###############################################################################",
    "\n",
    "## Creating softlinks for fastq files in ASF seq storage                     ##",
    "\n"
)
    
for (i in 1:nrow(dfDesign)){
    string <- paste0(
        "ln -s ",
        dfDesign$original.NGS[i],
        " ",
        dfDesign$NGS[i]
    )
        
    tempShellScriptVector <- c(
        tempShellScriptVector,
        string,
        "\n"
    )
}
    
tempShellScriptVector <- c(
    tempShellScriptVector,
    "## Done creating softlinks to ASF seq storage                                ##",
    "\n",
    "###############################################################################",
    "\n",
    "\n",
    "\n",
    "\n"
)


## Write to shell script ##
setwd(localWorkDir)

sink("create.fastq.softlinks.sh")        

for (i in 1:length(tempShellScriptVector)){
    cat(tempShellScriptVector[i])
}

sink()

## Add to overall shell script documentation
shellScriptVector <- append(
    shellScriptVector,
    tempShellScriptVector
)

## Remove end of line characters ##
print("Remove end of line characters: tr -d '\r' <create.fastq.softlinks.sh> conv.create.fastq.softlinks.sh")
## Run ##
print("sh conv.create.fastq.softlinks.sh")

# awk '{ sub("\r$", ""); print }' windows.txt > unix.txt
#awk 'sub("$", "\r")' unixfile.txt > winfile.txt

###############################################################################
## Create strings for documentation                                          ##
## fastq folder                                                              ##
## Create documentation module here                                          ##
documentationVector <- as.vector(
    NULL, 
    mode = "character"
)

documentationVector <- c(
    documentationVector,
    paste0(
        "FASTQ file location:", 
        pathToSeqStorageFolder
    ),   
    "FASTQ file names:",
    dfDesign$original.NGS,
    "",
    "Project FASTQ file names:", 
    dfDesign$NGS,
    ""
)



if (paired.end){
    documentationVector <- c(
        documentationVector,
        "Alignment mode: paired-end"
    )
} else {
    documentationVector <- c(
        documentationVector,
        "Alignment mode: single-end"
    )
}

if (stranded){
    documentationVector <- c(
        documentationVector,
        "This is a stranded dataset"
    )
} else {
    documentationVector <- c(
        documentationVector,
        "This is not a stranded dataset"
    )
}

documentationVector <- c(
    documentationVector,
    paste0("Reference genome/transcriptome: ",
           ref.genome
    )
)


## Add here: automated creation of powerpoint slide.
## Powerpoint slides to be generated:
## Documentation slide
## MA plot for each comparison slide
## PCA plot slide
## Heatmap slide
## Sample names/specifications slide
## Cluster dendrogram slide
## RNASeqQC slide

## End create strings for documentation                                      ##
###############################################################################

###############################################################################
# Create trim galore shell script                                             #
###############################################################################
setwd(localWorkDir)

TrimGaloreReturn <- create.trim.galore.shell.script(
    df.design = dfDesign, #design file with FASTQ column
    project.code = project.code,
    fastqDir  = fastqDir,
    localWorkDir = localWorkDir,
    paired.end = paired.end, 
    module.to.load = "Trim_Galore/0.4.2-foss-2016b",
    min.length = 25,
    quality = 20
)

dfDesign <- TrimGaloreReturn$df.design

shellScriptVector <- append(
    shellScriptVector,
    TrimGaloreReturn$ShellScriptVector
)


# If this shell script is created on a windows machine, don't forget to remove the end of line '\r' characters
print(
    paste0(
        "tr -d '\r' <", 
        project.code, 
        ".trimgalore.script.sh> conv.", 
        project.code, 
        ".trimgalore.script.sh"
    )
)


# Write dfDesign to file so it can be re-read once the alignment is done
setwd(localWorkDir)
#write.table(dfDesign, design.file, row.names= FALSE, sep="\t")


###############################################################################
# Align                                                                       #
###############################################################################
setwd(localWorkDir)

if (stranded){
    # Stranded
    forward_prob = "--forward-prob 0.0"
} else {
    # Unstranded
    forward_prob = "--forward-prob 0.5"
}

tempShellScriptVector <- create.alignment.shell.script(
    localWorkDir = localWorkDir,
    df.design = dfDesign,
    df.design.fastq.column =  "FASTQ_trimgalore",
    workDir = workdir,
    fastqDir =fastqDir,
    sampleDir = paste0(
        "/camp/stp/babs/working/boeings/Projects/",
        folder,
        "/workdir/RSEM/Ensembl"
    ),
    outputEnsDir = paste0(
        workdir, 
        "RSEM/Ensembl/"
    ),
    FASTQC = paste0(
        workdir, 
        "FASTQC/"
    ),
    logDir = paste0(
        workdir, 
        "logs/"
    ),
  
    msubDir = paste0(
        workdir, 
        "msub/"
    ),
  
    localMsubDir = paste0(
        localWorkDir, 
        "msub/"
    ),
  
    msubFile = paste0(
        project.code, 
        ".fastqc.and.rsem.all.samples.sh"
    ),
  
    sample.order = "",
    forward_prob = forward_prob,
    ENSG.index.human = ENSG.index.human,
    ENSG.index.mouse = ENSG.index.mouse,
    species = species,
    paired.end = paired.end,
    project.code = project.code
)

# If this shell script is created on a windows machine, don't forget to remove the end of line '\r' characters
print(
    paste0(
        "tr -d '\r' <",
        project.code,
        ".fastqc.and.rsem.all.samples.sh> conv.fastqc.and.rsem.all.samples.sh"
    )
)

# Run
print(
    paste0(
        "sh ", 
        "conv.fastqc.and.rsem.all.samples.sh"
    )
)

## Add to shell script documentation vector ##
shellScriptVector <- c(
    shellScriptVector,
    tempShellScriptVector
)

###############################################################################
# Prepare RNAseQC script                                                      #
###############################################################################
setwd(localWorkDir)
tempShellScriptVector <- create.rnaseqc.script(
    df.design = dfDesign,
    sample.column = "sample.id",
    project.code = project.code,
    project="SB_RNAseqQC",
    basedataDir=paste0(workdir, "RSEM/Ensembl"),
    bam.suffix = "STAR.genome.bam",
    GTFfile=GTFfile,
    rRNAfile=rRNAfile,
    genome.fa=genome.fa,
    paired.end = paired.end,
    refFlatFile=refFlatFile,
    ribosomalIntervalList=ribosomalIntervalList,
    bedFile = bedFile,
    strandSpecific = stranded
)


# If this shell script is created on a windows machine, don't forget to remove the end of line '\r' characters
print(
    paste0(
        "tr -d '\r' <", 
        project.code, 
        ".rnaseqc.script.sh> conv.", 
        project.code, 
        ".rnaseqc.script.sh"
    )
)

## Add to shell script documentation vector ##
shellScriptVector <- c(
    shellScriptVector,
    tempShellScriptVector
)

###############################################################################
## Produce shell script                                                      ##

fn <- paste0(
    localWorkDir,
    project_id,
    "documentationShell.script.sh"
)

sink(fn)

for (i in 1:length(shellScriptVector)){
    cat(shellScriptVector[i])
}

sink()

## use http://hilite.me/ for code conversion and display in script.php

## Done producing documentatin shell script                                  ##
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


samples = as.vector(
    unique(
        dfDesign$sample.id
        )
)

files <- paste(
    workdir, 
    "RSEM/Ensembl/", 
    samples, 
    ".genes.results", 
    sep=""
)

files = paste(
    files, 
    collapse = " "
)

RSEM.CMD = paste(
    "rsem-generate-data-matrix ", 
    files, 
    " > ./RSEM/", 
    count.data.file, 
    sep=""
)

setwd(localWorkDir)
sink("rsem.cmd.sh")
    cat("module load RSEM/1.2.31-foss-2016b"); cat("\n")
    cat(RSEM.CMD); cat("\n")
sink()

print("Remove end of line characters: tr -d '\r' <rsem.cmd.sh> conv.rsem.cmd.sh")


###############################################################################
# Run RSEM command on command line                                            #
###############################################################################
# Load RSEM module
# >srun -c1 --mem 10000 --pty bash
#print("module load RSEM/1.2.31-foss-2016b")
# Copy and paste the RSEM.CMD to the command line and run
# This will consolidate all individual RSEM results in a single file.
#print(RSEM.CMD)

## To do: incorpoate RSEM command into the alignment script ##

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
df.fpkm <- list.tpm.fpkm$df.fpkm

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

## First, run DESeq2 in normal mode ##
dfResList <- do.differential.expression.analyis(
    raw.count.dir = paste0(localWorkDir, "RSEM"),
    count.data.file = count.data.file,             #count data filename
    DEseq2Dir = paste0(localWorkDir, "DESeq2"),    #directory for results
    df.design = dfDesign,                         #dfDesign
    gene.id = primary.alignment.gene.id,                           #primary gene id after alignment 
    batch.mode = FALSE, #if true, dfDesign needs to contain a 'replicate' column 
    string.to.be.deleted.in.raw.counts.columns = paste0("X", gsub("/", ".",paste0(workdir, "RSEM/Ensembl/")))
)

df.summary <- dfResList$df.summary

shellScriptVector <- c(
    shellScriptVector,
    "##DESeq2 Differential Gene Expression Analysis ##",
    dfResList$docuVector,
    "##"
)


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
setwd(localWorkDir)
pca.dir <- paste0(
    localWorkDir, "DESeq2", "/pca.table.txt")
df.pca <- read.delim(pca.dir, header = TRUE, sep="\t", stringsAsFactors = FALSE)
names(df.pca) <- gsub("[.]", "_", names(df.pca))

## Documentation ##
# For each sample_group_[sample_group_name] 
# 

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
    cutOff = 1,
    zeroOneCol = "logFC_cut_off",
    selCol = "logFC"
)



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
## Create Excel output table                                                 ##
outCols <- c(
    gene.id.column,
    primary.alignment.gene.id,
    "gene_description",
    "gene_type",
    names(database.table)[grep("contrast_", names(database.table))],
    names(database.table)[grep("norm_counts_", names(database.table))],
    names(database.table)[grep("raw_counts_", names(database.table))],
    "count_cut_off", 
    "CoVar"
)

dfOutput <- unique(database.table[,outCols])

## Rename columns ##
names(dfOutput) <- gsub("norm_counts_", "TPM_", names(dfOutput))
comparisons <- names(dfDesign)[grep("comp_", names(dfDesign))]

for (i in 1:length(comparisons)){
    names(dfOutput) <- gsub(
        paste0("contrast_", i, "_"), 
        "", 
        names(dfOutput)
    )
}


outPutDir <- paste0(webSiteDir,project_id)
if (!dir.exists(outPutDir)){
    dir.create(outPutDir)
}

outPutDir <- paste0(webSiteDir,project_id, "/outputs")
if (!dir.exists(outPutDir)){
    dir.create(outPutDir)
}

outPutFile <- paste0(outPutDir, "/",project.code,".result.table.txt")

write.table(
    dfOutput, 
    outPutFile, 
    row.names = FALSE, 
    sep="\t"
)

outPutFile <- paste0(outputDir,project.code,".result.table.txt")
write.table(
    dfOutput, 
    outPutFile, 
    row.names = FALSE, 
    sep="\t"
)

print(paste0("Check ", outPutFile, " and convert to .xlsx."))

## Done creating Excel output table                                          ##
###############################################################################


###############################################################################
## Set cut-off for subsequent analyses                                       ##
dim(database.table)
database.table1 <- database.table[database.table$count_cut_off >= 1,]
database.table5 <- database.table[database.table$count_cut_off >= 5,]

###############################################################################
## Create metacore table                                                     ##
df.metacore <- database.table[database.table$count_cut_off >= 1,]

## select padj and logFC columns only ##
sel.vec <- names(df.metacore)[grep("contrast_", names(df.metacore))]
sel.vec <- sel.vec[-grep("stat", sel.vec)]
sel.vec <- sel.vec[-grep("lg10p", sel.vec)]
## Remove contrast_x_ prefix > remove front 11 characters
sel.vec <- append(gene.id.column, sel.vec)
df.metacore <- unique(df.metacore[,sel.vec])

## Remove contrast_x_ tag from column labels ##
comparisons <- names(dfDesign)[grep("comp_", names(dfDesign))]
for (i in 1:length(comparisons)){
    names(df.metacore) <- gsub(
        paste0("contrast_", i, "_"), 
        "", 
        names(df.metacore)
    )
}

outPutDir <- paste0(webSiteDir,project_id)
if (!dir.exists(outPutDir)){
    dir.create(outPutDir)
}

outPutDir <- paste0(webSiteDir,project_id, "/outputs")
if (!dir.exists(outPutDir)){
    dir.create(outPutDir)
}

outPutFile <- paste0(outPutDir, "/",project.code,".metacore.input.file.txt")

write.table(
    df.metacore, 
    outPutFile, 
    row.names = FALSE, 
    sep="\t"
)

outPutFile <- paste0(outputDir, project.code,".metacore.input.file.txt")

write.table(
    df.metacore, 
    outPutFile, 
    row.names = FALSE, 
    sep="\t"
)

print(
    paste0(
        "Upload ",   
        project.code, 
        ".metacore.input.file.txt"
        , 
        " to metacore and convert to .xls"
    )
)

## Done creating Metacore input file                                         ##
###############################################################################

###############################################################################
## Metacore analysis                                                         ##

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

minimalCols <- names(database.table)

###############################################################################
## Add additional plot columns from database                                 ##


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
        "DECIMAL(6,3) NULL DEFAULT NULL" = c("^CoVar$","NES", "logFC", "lg2_avg", "intensity", "^int", "iBAQ","^localization_prob$", "stat", "lg10p"),
        "DECIMAL(6,5) NULL DEFAULT NULL" = c("padj", "pvalue","^pep$"),
        "DECIMAL(6,1) NULL DEFAULT NULL" = c("norm_counts")
    ),
    new.table = TRUE
)

killDbConnections()

  
###############################################################################
# Do GSEA                                                                     #
###############################################################################
#Create rnk files for all logFC comparisons
#Remove lowly expressed genes from table. Cut off: Less than 1/counts per sample. Should be changed to TPM cut-off
# In this case, an average of 1 TPM is required for a sample to qualify


## If necessary, get database table from database ##
database.table <- import.db.table.from.db(
    password = db.pwd,
    dbname = prim.data.db,
    dbtable = rnaseqdbTableName
)

if (exists("minimalCols")) {
    database.table2 <- database.table[, minimalCols %in% names(database.table)]
} else {
    database.table2 <- database.table
}

# rmVec <- c(
#         grep("contrast_1", names(database.table2)),
#         grep("contrast_2", names(database.table2)),
#         grep("contrast_7", names(database.table2)),
#         grep("contrast_8", names(database.table2))
#     )
# 
# database.table2 <- database.table2[,-rmVec]

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

write.table(
    dfRefGmt,
    paste0(localGmtDir, gmtFileName),
    col.names = FALSE,
    row.names = FALSE,
    sep="\t"
)

## At present it is necessary to open, save and close the gmt file once in Excel ##

## Function to run GSEA ##
# parameters: GSEA.rnk file names
# gmt file > made from list 
# chip file (optional)

## This appears to work ##


## Create 
#java -Xmx2512m -cp /camp/stp/babs/working/boeings/Projects/software/gsea-3.0.jar xtools.gsea.GseaPreranked -gmx /camp/stp/babs/working/boeings/Projects/reference_data/GSEA.gmt.files/20160508.rna.seq.txn.analysis.gmt -rnk /camp/stp/babs/working/boeings/Projects/126_SL_MP_RNA_seq_macrophage_stimulation_timecourse/workdir/GSEA/contrast_1_logFC_IFNARLA_vs_WTLA30.rnk -rpt_label contrast_1_rnaSeqTxnTest -out /camp/stp/babs/working/boeings/Projects/126_SL_MP_RNA_seq_macrophage_stimulation_timecourse/workdir/GSEA -collapse false -mode Max_probe -norm meandiv -nperm 10 -scoring_scheme classic -include_only_symbols true -make_sets true -plot_top_x 20 -rnd_seed timestamp -set_max 2500 -set_min 10 -zip_report false -gui false
## end of working

contrasts <- names(database.table2)
contrasts <- contrasts[grep("logFC", contrasts)]
contrasts <- contrasts[contrasts != "logFC_cut_off"]



setwd("GSEA")
sink("GSEAcommands.sh")
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

# Did MysigDB C2, C3, c6, ChEA, TRANSFAC/JASPAR #


###############################################################################
## Optional Module: Perform hypergeometric test                              ##

## If performed on CAMP include:
## include relevant functions ##
#funPath <- paste0(
#    "working/boeings/Stefan/protocol_files/github/boeings/packages/20180115.functions.SBwebtools.r"
#)

#source(funPath)

## Set parameters ##
# workdir <- paste0(
#     "/camp/stp/babs/working/working/boeings/Projects/",
#     "143_jsl_lw_IRDS_genes_cat_enrichments/workdir"
#)

# setwd(workdir)
#host = "www.biologic-db.org"
#dbname = "reference_categories_db_new"
#db.user = "boeings"
list.db.tables.in.db(
    password = db.pwd
)

dbtable.vec <- c("es_lab_categories", "SGP_from_GEO_up_down_combined")



## Get test gene set ##
testGeneSet <- as.vector(
    unique(
    database.table[
        database.table$contrast_2_logFC_siTFAP2C1Only_vs_Ctrl < -0.5 &
            database.table$contrast_4_logFC_nonAligned_vs_aligned < -1,
        gene.id.column
    ]
    )
)

## Get background gene set ##

backgroundGeneVec <- as.vector(
    unique(database.table[, gene.id.column])
)


backgroundGeneVec <- as.vector(
    unique(append(testGeneSet, backgroundGeneVec))
)

## Done retrieving background genes


## Creating enrichment dataframe ##
setwd(localWorkDir)

dfRes <- surveyDbCategories(
    genesOfInterestVec = testGeneSet,
    backgroundGeneVec = backgroundGeneVec,
    host = "www.biologic-db.org",
    dbname = ref.cat.db,
    db.pwd = db.pwd,
    db.user = db.user,
    dbtable.vec = dbtable.vec,
    gene.id.column =  "hgnc_symbol",
    project.code = paste0(project.code, ".TFAP2Cdown.algnmt.down"),
    NminGenesInCat = 1
)

write.table(
    dfRes,
    paste0(project.code, ".hypergeometric.test.result.txt"),
    row.names = FALSE,
    sep = "\t"
)
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


##############
# Mop up GSEA files and put into database
##############
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


## Create GSEA plot directory by executing farm command ##
## Check [project.code.enrichment.file and copy to outputs]

###############################################################################
## Create reference category selection for this project                      ##
## This procedure will create a tailored category reference table for Ley lab 
## RNA-Seq experiments VTL87


## Get base dataset VTL87 ##
## Identify all genes that are significantly regulated in at least one 
## comparison ##
## Require at least five of those genes to be present in the above for 
## inclusion, with the exception of lab categories

#list.db.tables.in.db(password = db.pwd, dbname = prim.data.db)
###############################################################################
## Variables and parameters                                                  ##
###############################################################################
#prim.data.db <- "vtl_data"



###############################################################################
## Optional Module: Correlation analysis                                     ##

## If required, run the module 
# templateRNAseqAnalysis_pluginModule_correlationAnalysis.R
# here.

## DoneOptional Module: Correlation analysis                                 ##
###############################################################################


if (exists("tables")){
    print(paste0(
        "Do you want to use GSEA tables for reference selection(",
        paste(
            tables, collapse = ", "),
        ")?"
    )
    )
    refTableSelection <- tables
} else {
    print("Define refTableSelection variable for reference data selection.")
}


list.db.tables.in.db(
    password = db.pwd,
    dbname = "reference_categories_db_new"
)

tables <- c(
    #tables,
    "mysigdb_h_hallmarks",
    "TRANSFAC_and_JASPAR_PWMs",
    #"ag_lab_categories",
    "mysigdb_c2_1329_canonical_pathways",
    "mysigdb_c5_BP",                               
    "mysigdb_c5_MF",
    "SGP_from_GEO_up_down_combined",
    "networkcategories"
)

refTableSelection <- tables



print(
    paste0(
        "Table selection to carry forward:",
        paste0(
            "Do you want to use GSEA tables for reference selection(",
            paste(
                refTableSelection, collapse = ", "),
            ")?"
        )
    )
)

###############################################################################
## Optional Module: Find best correlation of timepoints with single gene     ##
##permutations                                                               ##       
## Download single gene mutations as gmt file                                ##

# Include optional module from plug in selection
# find best correlated category       

## Done optional module                                                      ##
###############################################################################

##Define table names


tableNameVector <- refTableSelection


tableNameVector <- c(
    'Hallmarks',
    'TF binding sites',
    #'Gould lab',
    'Pathways',
    'GO-BP',
    'GO-MF',
    'Single Gene Pertubations',
    'Protein-protein interactions'
    
)

dfCheck <- data.frame(
    tableNameVector, 
    tables
)

## Check table name to tables assignment ##
dfCheck


## Define relevant genes ##

dfData <- import.db.table.from.db(password = db.pwd,
                                  dbname = prim.data.db,
                                  dbtable = rnaseqdbTableName)

dfData <- dfData[dfData$cluster_order > 0, ]

relevant.genes <- as.vector(unique(dfData[dfData$cluster_order > 0, gene.id.column]))
    
   

###############################################################################
## Lab category                                                              ##
###############################################################################

for (i in 1:length(refTableSelection)) {
    df.ref <- import.db.table.from.db(dbname = ref.cat.db,
                                      dbtable = refTableSelection[i],
                                      password = db.pwd)
    
    ## Remove temp categories ##
    temPos <- grep("temp_", df.ref$cat_type)
    
    if (length(temPos) > 0) {
        df.ref <- df.ref[-temPos,]
    }
    
    
    
    df.temp <- add2labCatSelectionDBtable(
        df.ref = df.ref,
        cat_group_name = tableNameVector[i],
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

## Created all categories ##

## Filter for temp categories ##
#temPos <- grep("temp_", df.db.table$cat_id)    
    
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
# timecourse.cat.lines <- create.timecourse.cat.lines(
#     dfDesign
# )

###############################################################################
# Transfer GSEA plot files to server                                          #
###############################################################################
# Before running the GSEA shell scripts, do the following replacements:
# "(" > "\("
# ")" > "\)"
# "$" > "\$"
# "&" > "\&"
# Run GSEA shell scripts to transfer GSEA image files
# Transfer folder to server
            
###############################################################################
# Create microwebsite                                                         #
###############################################################################

## Color samples by sample group ##
# Get sample order from database.table
sample.order <- sort(names(database.table)[grep("norm_counts_", names(database.table))])

sample.order <- c(
    'norm_counts_NGF_DMSO_1',
    'norm_counts_NGF_DMSO_2',
    'norm_counts_NGF_DMSO_3',
    'norm_counts_NGF_DMSO_4',
    'norm_counts_NGF_DMSO_5',
    'norm_counts_antiNGF_DMSO_1',
    'norm_counts_antiNGF_DMSO_2',
    'norm_counts_antiNGF_DMSO_3',
    'norm_counts_antiNGF_DMSO_4',
    'norm_counts_antiNGF_DMSO_5',
    'norm_counts_antiNGF_ISRIB_1',
    'norm_counts_antiNGF_ISRIB_2',
    'norm_counts_antiNGF_ISRIB_3',
    'norm_counts_antiNGF_ISRIB_4',
    'norm_counts_antiNGF_DLKi_1',
    'norm_counts_antiNGF_DLKi_2',
    'norm_counts_antiNGF_DLKi_3',
    'norm_counts_antiNGF_DLKi_4'
)

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
    venn.slider.selector.strings = c("_logFC", "_padj"),
    plot.selection.strings = c("_logFC", "_lg10p"), #NA, #strings to grep from col names for plot display
    plate.view.db.table = NA,
    plate.view.column.vec = NA,
    cat.seletion.table.vec = c(cat.ref.db, cat.ref.db.table)
)
            
            
###########################################
# Do multiQC and add to QC section 
###########################################
# run 
# >> module use /camp/stp/babs/working/software/modules/all
# >> module load multiqc/0.9-2016b-Python-2.7.12

## Copy mga.tsv file into the working directory ##
print(
    paste0(
        "cp ",
        dataFolder,
        "/mga.tsv ",
        workdir
    )
)

# run in top project directory:
# >> multiqc .
# multiqc /camp/stp/babs/working/boeings/Projects/103_VTL_ES_RNA_seq_BAFF_timecourse_hs/workdir/RSEM/Ensembl /camp/stp/babs/working/boeings/Projects/103_VTL_ES_RNA_seq_BAFF_timecourse_hs/workdir/RSEM/Ensembl/RNAseQC /camp/stp/babs/working/boeings/Projects/103_VTL_ES_RNA_seq_BAFF_timecourse_hs/workdir/logs /camp/stp/babs/working/boeings/Projects/103_VTL_ES_RNA_seq_BAFF_timecourse_hs/workdir/FASTQC_stranded
print("module purge; module use /camp/stp/babs/working/software/modules/all; module load multiqc/1.3-2016b-Python-2.7.12")
print("multiqc --config /camp/stp/babs/working/escudem/bin/multiqc_config.yaml .")

print(paste0("multiqc . -o ", webSiteDir))
print(paste0("multiqc . -o ", workdir, "outputs"))

## Transfer multiqc_report.html to outputs ##

## list directories to include in multiqc search
#paste0(workdir, "logs")
paste0(workdir, "RSEM/Ensembl")
paste0(workdir, "RSEM/Ensembl/RNAseQC")
paste0(workdir, "FASTQC")


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

# createNewProject(
#     dbname = "reference_categories_db_new",
#     dbtable = "project_description_table",
#     password = db.pwd,
#     project_name = "Tumor_Disc_Development",
#     project_lab = ";Gould;",
#     project_description = "Characterization of turmor disc development in glioblastoma."
# )

## Done creating new project                                                 ##
###############################################################################

###############################################################################
## List project in projects table                                            ##

addProject2ProjectTable(
   dbname = "reference_categories_db_new",
   dbtable = "project_db_table",
   password = db.pwd,
   experiment_id= project_id,
   experiment_question = "Which genes are regulated by the ISRIB inhibitor in primary murine neurons?",
   experiment_description = "<b>RNA Seq analysis</b>. Overall design RNA sequencing from primary DRG neurons grown for 4 days in vitro and deprived of NGF for 4.5 h in the presence of DMSO, DLKi 400 nM or ISRIB 400 nM. NGF with DMSO was used as control.",
   experiment_owner = ";Genentech;",
   experiment_lab = paste0(";Genentech;"),
   experiment_viewers = paste0(";", paste(user_ids, collapse = ";"),";"),
   experiment_project = ";Tumor_Disc_Development;",
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

## End of template                                                           ##
###############################################################################