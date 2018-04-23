###############################################################################
## Create TCGA download scripts for tumor/normal pairs                       ##
## Stefan Boeing                                                             ##
## CRICK Bioinformatics                                                      ##

#Solution
# Step 1: Create list of all available TCGA WXS files
# Step 2: Add TCGA-ID in order to distinguish tumor/normal samples
# Step 3: Create download manifests by TCGA cancer type

###############################################################################
## Step 1: Create list of all available TCGA WXS files                       ##
###############################################################################

## get default display fields ##
library(RCurl)
library(jsonlite)

request <- "https://api.gdc.cancer.gov/files/_mapping"

resp <- getURL(
    request
)

json_file <- fromJSON(
    resp
)

df.resp <- as.data.frame(json_file["defaults"])

defaultFields <- as.vector(df.resp[,1])
## Done with default field creation ##

#cases.samples.tumor_descriptor,cases.samples.tissue_type,cases.samples.sample_type,cases.samples.submitter_id,cases.samples.sample_id,cases.samples.portions.analytes.aliquots.aliquot_id,cases.samples.portions.analytes.aliquots.submitter_id","size"

## Add custom fileds ##
selFields <- c(
    defaultFields,
    "cases.project.project_id",
    "cases.submitter_id",
    "cases.case_id",
    "cases.samples.tumor_descriptor",
    "cases.samples.sample_type"
)

## Remove superflous columns ##
selFields <- selFields[selFields != "created_datetime"]
selFields <- selFields[selFields != "acl"]

operator <- paste0(
    '{
        "op":"and",
        "content":[',
            '{
                "op":"=",
                "content":{
                    "field":"files.data_format",
                    "value":"BAM"
                }
            },',
    
            '{
                "op":"=",
                "content":{
                    "field":"cases.project.program.name",
                    "value":"TCGA"
                }
            },',
    
            '{
                "op":"=",
                "content":{
                    "field":"experimental_strategy",
                    "value":"WXS"
                }
            }
        ]
    }'
)

# percent-encode #
#library(urltools)
#operator <- url_encode(operator)
# changed to RCurl package function
operator <-curlEscape(operator)

# Mapping is always compatible with filters

request <- paste0(
    "https://api.gdc.cancer.gov/files?filters=",
    operator,
    "&from=0&size=60000",
    "&pretty=true",
    "&fields=",
    paste(selFields, collapse=","),
    "&related_files=true"
)

library("jsonlite")
resp <- getURL(
    request
)
json_file <- fromJSON(resp)
df.resp <- as.data.frame(json_file$data$hits)

## Need to unlist the cases column ##
size <- length(unlist(df.resp$cases[1]))

for (i in 1:size){
    colName <- names(unlist(df.resp$cases[1])[i])
    df.resp[[colName]] <- sapply(df.resp$cases, function(x) unlist(x)[i])
}

df.resp$cases <- NULL

# Add tumor/normal columns #
df.resp[["tumor"]] <- 0
df.resp[["normal"]] <- 0

sampleNames <- unique(df.resp$samples.sample_type)

## Mark tumor samples ##
tumorSampleNames <- sampleNames[grep("Tumor", sampleNames)]
tumorSampleNames <- c(
    tumorSampleNames,
    "Metastatic",
    "Primary Blood Derived Cancer - Peripheral Blood"
)

## Mark normal samples ##
normalSampleNames <- sampleNames[grep("Normal", sampleNames)]

df.resp[df.resp$samples.sample_type %in% tumorSampleNames, "tumor"] <- 1
df.resp[df.resp$samples.sample_type %in% normalSampleNames, "normal"] <- 1

## Mark tumor-normal pairs ##
df.resp[["TNpaired"]] <- ""
ids <- unique(df.resp$submitter_id)

for (i in 1:length(ids)){
    t <- sum(df.resp[df.resp$submitter_id == ids[i], "tumor"])
    n <- sum(df.resp[df.resp$submitter_id == ids[i], "normal"])
    if ((t > 0) & (n > 0)){
        df.resp[df.resp$submitter_id == ids[i], "TNpaired"] <- "+"
    }
}

dfDownload <- df.resp[df.resp$TNpaired == "+",]


## Done with Step 1                                                          ##
###############################################################################


###############################################################################
## Identify downloaded bam files                                             ##
cancerTypes <- sort(unique(dfDownload$project.project_id))
sink("modify.sh")
for (i in 1:length(cancerTypes)){
    cmd <- paste0("find /srv/shared/vanloo/TCGA/",cancerTypes[i]," -type f -execdir basename {} \\; >> /home/sboeing/downloaded.files.txt")
    cat(cmd);cat("\n")
}

sink()

## Transfer files to local computer ##
print("scp vanloo-download-1:/home/sboeing/downloaded.files.txt ./")

FN <- "/Users/boeings/Desktop/downloaded.files.txt"

dfPresent <- read.delim(
    FN,
    header=FALSE,
    sep="\t",
    stringsAsFactors = FALSE
)
    

dfDownloadComplete <- dfDownload[dfDownload$file_name %in% dfPresent[,1],]

dfDownloadToDoList <- dfDownload[!(dfDownload$file_name %in% dfPresent[,1]),]

dfDownload <- dfDownloadToDoList

## Remove already downloaded files from dfDownload                           ##

## Done editing dfDownload ##


###############################################################################
## Step 2: Create download manifests by TCGA cancer type                     ##
cancerTypes <- sort(unique(dfDownload$project.project_id))

## Select cancerTypes ##
shellScriptFNvec <- as.vector(NULL, mode = "character")
TCGApath <- "/srv/shared/vanloo/TCGA/"
workDir <- "/Volumes/babs/working/boeings/Projects/184_CSL_NO_TCGA_data_transfer_to_emedlab/workdir/"

## Number of files per download script ## 
increment <- 200



## Create download script for each cancer type ##

for (i in 1:length(cancerTypes)){
    dfTypeDownload <- dfDownload[dfDownload$project.project_id == cancerTypes[i],]
    cancer.type <- cancerTypes[i]
    counter <- 0
    
    while (nrow(dfTypeDownload) > 0){
        shellScriptVector <- as.vector(NULL, mode = "character")
        counter <- counter + 1
        
        if (nrow(dfTypeDownload) >= increment){
            dfTypeDownloadIncrement <-  dfTypeDownload[1:increment,]
            dfTypeDownload <- dfTypeDownload[(increment+1):nrow(dfTypeDownload),]
        } else {
            dfTypeDownloadIncrement <-  dfTypeDownload[1:nrow(dfTypeDownload),]
            dfTypeDownload <- dfTypeDownload[0,]
        }
         
        
        ## Create file name and add to file vector ##
        DLfname <- paste0(workDir, cancer.type, ".WXS",".download.script-",counter,".sh")
        shellScriptFNvec <- c(
            shellScriptFNvec,
            DLfname
        )
        
        # Get uuids #
        uuIDs <- as.vector(unique(dfTypeDownloadIncrement$file_id))
        uuIDs[1] <- paste0('analysis_ids="', uuIDs[1])
        uuIDs <- c(
            uuIDs,
            '"'
        )
    
        ## Create shell script for the download ##
        shellScriptVector <- c(
            shellScriptVector,
            '#!/usr/bin/sh',
            '###############################################################################',
            '### Loading required modules                                                ###',
            '###############################################################################',
            'module purge',
            'module load gdc-client/1.3.0-foss-2016b-Python-2.7.12',
            '',
            '###############################################################################',
            '### VARIABLES                                                               ###',
            '###############################################################################',
            paste0('cancer_type="', cancer.type,'"'),
            paste0('cancer_dir="', TCGApath, '${cancer_type}"'),
            #paste0('md5_documentation="',TCGApath,'${cancer_type}/md5_documentation"'),
            #paste0('gto_dir="',TCGApath,'${cancer_type}/gto_dir"'),
            '### Path to token file ###',
            'token="/home/sboeing/token.txt"',
            '### Set uuids to download ###',
            uuIDs,
            '',
            '###############################################################################',
            '#### Download bam files and put into cancer.type_WXS_directory              ###',
            '###############################################################################',
            '',
            '### Create requiered directories ###',
            'if [ ! -d $cancer_dir ]; then',
            '    mkdir $cancer_dir',
            'fi',
            '',
            '### Create requiered directories ###',
            'cancer_dir="$cancer_dir/WXS"',
            'if [ ! -d $cancer_dir ]; then',
            '    mkdir $cancer_dir',
            'fi',
            '',
            #'if [ ! -d $gto_dir ]; then',
            #'    mkdir $gto_dir',
            #'fi',
            #'',
            #'if [ ! -d $md5_documentation ]; then',
            #'    mkdir $md5_documentation',
            #'fi',
            '',
            '### Download and extract all samples ###',
            #'',
            #'cd $cancer_dir',
            '',
            'for analysis_id in $analysis_ids',
            '    do',
       
            # curl --remote-name --remote-header-name --header "X-Auth-Token: $token" \\'https://api.gdc.cancer.gov/data/b119488a-0379-4e84-9d3d-f4aab253cb33\\'
            '        gdc-client download ${analysis_id} -t $token -d $cancer_dir',
            '',
            '        echo "${analysis_id}" >> $cancer_dir/download.completed.ids.txt',
            #mv *.gto ./$gto_dir/
            #'    BAMfile="https:/api.gdc.cancer.gov/data/${analysis_id}/*.bam $cancer_dir"',
            #'    BAIfile="https:/api.gdc.cancer.gov/data/${analysis_id}/*.bai"',
            '',
            #'    if [ ! -f "$BAMfile" ]; then',
            '        mv $cancer_dir/https:/api.gdc.cancer.gov/data/${analysis_id}/*.bam $cancer_dir',
            #'    fi',
            '',
            #'#    if [ ! -d ${analysis_id} ]; then',
            '        rm -r $cancer_dir/https:/api.gdc.cancer.gov/data/${analysis_id}',
            #'#    fi',
        
            '    done'
        
        )
    
        ## Write shell script to file ##
        sink(DLfname)
        for (j in 1:length(shellScriptVector)){
            cat(shellScriptVector[j]); cat("\n");
        }
        sink()
    } # end of while loop
}

## Done with step 2                                                          ##
###############################################################################

###############################################################################
## Create compounded shell scripts                                           ##

nCompoundScripts <- 40

ScriptsPerScript <- ceiling(length(shellScriptFNvec)/nCompoundScripts)

shellScriptFNvec <- gsub(
    "/Volumes/babs/working/boeings/Projects/184_CSL_NO_TCGA_data_transfer_to_emedlab/workdir/",
    "",
    shellScriptFNvec
)
shellScriptFNvec <- paste("sh", shellScriptFNvec)

script <- 0

for (i in 1:nCompoundScripts){
    assign(paste0("compoundSHscript", i), as.vector('#!/usr/bin/sh', mode="character"))
}

while (length(shellScriptFNvec) > 0){
    if (script == nCompoundScripts){
        script <- 1
    } else {
        script <- script +1
    }
    
    assign(paste0("compoundSHscript", script), c(get(paste0("compoundSHscript", script)), shellScriptFNvec[1]))
    shellScriptFNvec <- shellScriptFNvec[-1]       
    ## Add first element to the next script ##
    
}

## Write compount scripts to file ##
for (i in 1:nCompoundScripts){
    DLfname <- paste0(workDir, "compoundSHscript", i,".sh")
    ## Write shell script to file ##
    sink(DLfname)
    for (j in 1:length(get(paste0("compoundSHscript", i)))){
        cat(get(paste0("compoundSHscript", i))[j]); cat("\n");
    }
    sink()
}

## Done with step 3                                                          ##
###############################################################################

## Upload scripts to eMedLab ##
print("scp *.sh vanloo-download-1 vanloo-download-1:/home/sboeing/")

###############################################################################

