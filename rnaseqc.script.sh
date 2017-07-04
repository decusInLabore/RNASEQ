#!/bin/sh

## First part : RNAseQC
## Second part: MultiQC

###############################################################################
## First part: RNASEQC                                                       ##
###############################################################################

#Copy this shell script into the project directory and run it from there.
###############################################################################
###VARIABLES###################################################################
###############################################################################
project="SB_RNAseqQC"
projectID=""


#AlignDir gives the path to the input bam files
alignDir="/path/to/bamfiles/"
GTFfile="path/to/Mus_musculus.GRCm38.86.rnaseqc.gtf"
rRNAfile="path/to/Mus_musculus.GRCm38.86.rRNA.list"

## Input files are genome aligned BAM files.                                 ##
## Input bam files have to follow the following naming convention for this   ##
## script: [sample_name].[sample_suffix]                                     ##
## In this case, sample suffix is 'STAR.genome.bam'							 ##
## and the full file name individual_sample_name.star.genome.bam             ##

## This version of the script deals with the single end protocol.            ##
## To deal with paired end input files,                                      ##
## comment out line 331 and 340 '-singleEnd \'                               ##

## The sample variable lists all individual parts of the input bam files     ##

samplesuffix="STAR.genome.bam"

genome_fa="path/to/genome/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa"
samples="DM_0_1
DM_0_2
DM_0_3
DM_0_4
DM_0_5
DM_0_6
DM_05h_1
DM_05h_2
DM_05h_3
DM_05h_4
DM_05h_5
DM_05h_6
DM_1h_1
DM_1h_2
DM_1h_3
DM_1h_4
DM_1h_5
DM_1h_6
DM_2h_1
DM_2h_2
DM_2h_3
DM_2h_4
DM_2h_5
DM_2h_6
MT_0_1
MT_0_2
MT_0_3
MT_0_4
MT_0_5
MT_0_6
MT_05h_1
MT_05h_2
MT_05h_3
MT_05h_4
MT_05h_5
MT_05h_6
MT_1h_1
MT_1h_2
MT_1h_3
MT_1h_4
MT_1h_5
MT_1h_6
MT_2h_1
MT_2h_2
MT_2h_3
MT_2h_4
MT_2h_5
MT_2h_6
WT_0_1
WT_0_3
WT_0_4
WT_0_5
WT_0_6
WT_0_7
WT_0_8
WT_0_9
WT_05h_1
WT_05h_2
WT_05h_3
WT_05h_4
WT_05h_5
WT_05h_6
WT_05h_7
WT_05h_8
WT_05h_9
WT_1h_1
WT_1h_2
WT_1h_3
WT_1h_4
WT_1h_5
WT_1h_6
WT_1h_7
WT_1h_8
WT_1h_9
WT_2h_1
WT_2h_2
WT_2h_3
WT_2h_4
WT_2h_5
WT_2h_6
WT_2h_7
WT_2h_8
WT_2h_9"


#################################################################################
##FUNCTIONS######################################################################
#################################################################################

wait_on_lsf() { ## wait on jobs{
sleep 300
n=`squeue --name=$project  | wc -l`
while [ $n -ne 1 ]
do
n=`squeue --name=$project  | wc -l`
((z=$n-1))
#number of running
echo "$project jobs ($projectID) running: $z"
#number of pending
sleep 300
done
}



## End of function                                                             ##
#################################################################################


#################################################################################
# Prere bam files for RNASeQC                                                   #
#################################################################################
#################################################################################
# AddOrReplaceReadGroups                                                        #
#################################################################################

projectID="AddOrReplaceReadGroups"
module load R/3.3.1-foss-2016b-bioc-3.3-libX11-1.6.3

echo "module load R/3.3.1-foss-2016b-bioc-3.3-libX11-1.6.3" >> commands.txt
module load picard/2.1.1-Java-1.8.0_112

echo "module load picard/2.1.1-Java-1.8.0_112" >> commands.txt
echo "redo headers on Tophat bam output"
echo "#Submitted jobs" >> commands.txt 
#projectID="headers"
for sample in $samples
   do 
      echo "java -jar ${EBROOTPICARD}/picard.jar AddOrReplaceReadGroups \
      I=${alignDir}/${sample}.${samplesuffix} \
      O=${alignDir}/${sample}.accepted_hits.readgroups.bam \
      RGID=$sample \
      RGCN=CCCB \
      RGLB=lib1 \
      RGPL=ILLUMINA \
      RGPU=NA \
      RGSM=accepted_hits.bam" >> commands.txt
      sbatch --wrap "java -jar ${EBROOTPICARD}/picard.jar AddOrReplaceReadGroups \
      I=${alignDir}/${sample}.${samplesuffix} \
      O=${alignDir}/${sample}.accepted_hits.readgroups.bam \
      RGID=$sample \
      RGCN=TheFrancisCrickInstitute \
      RGLB=lib1 \
      RGPL=ILLUMINA \
      RGPU=NA \
      RGSM=accepted_hits.bam" --job-name=$project -c 1 --mem-per-cpu=7000 -o $sample.addorreplacereadgroups.slurm >> commands.txt 
   done
##wait on jobs
wait_on_lsf


#################################################################################
# SortSam                                                                       #
#################################################################################

projectID="SortSam"
echo "###########################################################################################" >>  commands.txt
echo "SortSam co-ordinate sort" >>  commands.txt
echo "SortSam output"
#projectID="coordinate_sort"
for sample in $samples
   do
      echo "java -Xmx10g -Djava.io.tmpdir='pwd'/tmp -jar ${EBROOTPICARD}/picard.jar SortSam \
      I=${alignDir}/${sample}.accepted_hits.readgroups.bam \
      O=${alignDir}/${sample}.accepted_hits.readgroups.sorted.bam \
      SO=coordinate \
      TMP_DIR='pwd'/tmp" >> commands.txt 
      sbatch --wrap "java -Xmx10g -Djava.io.tmpdir='pwd'/tmp -jar ${EBROOTPICARD}/picard.jar SortSam \
      I=${alignDir}/${sample}.accepted_hits.readgroups.bam \
      O=${alignDir}/${sample}.accepted_hits.readgroups.sorted.bam \
      SO=coordinate \
      TMP_DIR='pwd'/tmp" --job-name=$project -c 1 --mem-per-cpu=7000 -o $sample.sortsam.slurm >> commands.txt 
   done
##wait on jobs
wait_on_lsf

for sample in $samples
do
echo "rm ${alignDir}/${sample}.accepted_hits.readgroups.bam" >> commands.txt  
rm ${alignDir}/${sample}.accepted_hits.readgroups.bam  
done

#################################################################################
# Reorder SAM                                                                   #
#################################################################################

projectID="ReorderSam"
echo "###########################################################################################" >>  commands.txt
echo "reorder reads to match contigs in the reference" >> commands.txt
echo "Reorder reads to match contigs in the reference"
for sample in $samples
do	
echo "java -Xmx10g -Djava.io.tmpdir='pwd'/tmp -jar ${EBROOTPICARD}/picard.jar ReorderSam \
I=${alignDir}/${sample}.accepted_hits.readgroups.sorted.bam \
O=${alignDir}/${sample}.accepted_hits.readgroups.sorted.reordered.bam \
REFERENCE=${genome_fa} \
TMP_DIR='pwd'/tmp" >> commands.txt 

sbatch --wrap "java -Xmx10g -Djava.io.tmpdir='pwd'/tmp -jar ${EBROOTPICARD}/picard.jar ReorderSam \
I=${alignDir}/${sample}.accepted_hits.readgroups.sorted.bam \
O=${alignDir}/${sample}.accepted_hits.readgroups.sorted.reordered.bam \
REFERENCE=${genome_fa} \
TMP_DIR='pwd'/tmp" --job-name=$project -c 1 --mem-per-cpu=7000 -o $sample.reordersam.slurm >> commands.txt 

done
##wait on jobs
wait_on_lsf

for sample in $samples
   do
      echo "rm ${alignDir}/${sample}.accepted_hits.readgroups.sorted.bam" >> commands.txt
      rm ${alignDir}/${sample}.accepted_hits.readgroups.sorted.bam
   done

#################################################################################
# MarkDuplicates                                                                #
#################################################################################

projectID="MarkDuplicates"
echo "###########################################################################################" >>  commands.txt
projectID="markdups"
echo "mark duplicates" >> commands.txt
echo "mark duplicates"
   for sample in $samples
      do
         echo "java -Xmx10g -Djava.io.tmpdir='pwd'/tmp -jar ${EBROOTPICARD}/picard.jar MarkDuplicates \
         I=${alignDir}/${sample}.accepted_hits.readgroups.sorted.reordered.bam \
         O=${alignDir}/${sample}.accepted_hits.readgroups.sorted.reordered.dedup.bam \
         METRICS_FILE=${alignDir}/${sample}/dup_metrics.txt \
         TMP_DIR='pwd'/tmp" >> commands.txt

         sbatch --wrap "java -Xmx10g -Djava.io.tmpdir='pwd'/tmp -jar ${EBROOTPICARD}/picard.jar MarkDuplicates \
         I=${alignDir}/${sample}.accepted_hits.readgroups.sorted.reordered.bam \
         O=${alignDir}/${sample}.accepted_hits.readgroups.sorted.reordered.dedup.bam \
         METRICS_FILE=${alignDir}/${sample}/dup_metrics.txt \
         TMP_DIR='pwd'/tmp" --job-name=$project -c 2 --mem-per-cpu=7000 -o $sample.markduplicates.slurm >> commands.txt 
      done
##wait on jobs
wait_on_lsf

##############################
for sample in $samples
   do
      echo "rm ${alignDir}/${sample}.accepted_hits.readgroups.sorted.reordered.bam" >> commands.txt
      rm ${alignDir}/${sample}.accepted_hits.readgroups.sorted.reordered.bam
   done
#################################################################################
# Samtool indexing                                                              #
#################################################################################

projectID="Samtool indexing"
module load SAMtools/1.3.1-foss-2016b
echo "###########################################################################################" >>  commands.txt
echo "module load SAMtools/1.3.1-foss-2016b" >> commands.txt
echo "Samtool Indexing "
for sample in $samples
   do
      echo "samtools index ${alignDir}/${sample}.accepted_hits.readgroups.sorted.reordered.dedup.bam" >> commands.txt
      sbatch --wrap "samtools index ${alignDir}/${sample}.accepted_hits.readgroups.sorted.reordered.dedup.bam" --job-name=$project -c 1 --mem-per-cpu=7000 -o $sample.samtools.index.slurm >> commands.txt 
   done
##wait on jobs
wait_on_lsf

#################################################################################
# Run RNASeqC                                                                   #
#################################################################################

projectID="RNAseQC"
echo "###########################################################################################" >>  commands.txt
module load RNA-SeQC/1.1.8-Java-1.7.0_80

echo "module load RNA-SeQC/1.1.8-Java-1.7.0_80" >> commands.txt
echo "make RNASeqC sample list"
cd ${alignDir}

if [ ! -d ${projectID} ]; then
  mkdir ${projectID}
fi


echo "sample list" > ${alignDir}/${projectID}/sample.list
for sample in $samples
do
echo -e "$sample	${alignDir}/${sample}.accepted_hits.readgroups.sorted.reordered.dedup.bam	NA" >>${alignDir}/${projectID}/sample.list
done
#wait

echo "run RNA-seqQC"

#RNASeqQC requires Java version 1.7 and does not run on the most recent version. 

echo "java -Djava.io.tmpdir='pwd'/tmp -jar ${EBROOTRNAMINSEQC}/RNA-SeQC_v1.1.8.jar \
-singleEnd \
-o ${alignDir}/${projectID}/ \
-r ${genome_fa} \
-s ${alignDir}/${projectID}/sample.list \
-t $GTFfile \
-gatkFlags '-S SILENT -U ALLOW_SEQ_DICT_INCOMPATIBILITY' \
-rRNA $rRNAfile " >> commands.txt

sbatch --wrap "java -Djava.io.tmpdir='pwd'/tmp -jar ${EBROOTRNAMINSEQC}/RNA-SeQC_v1.1.8.jar \
-singleEnd \
-o ${alignDir}/${projectID}/ \
-r ${genome_fa} \
-s ${alignDir}/${projectID}/sample.list \
-t $GTFfile \
-gatkFlags '-S SILENT -U ALLOW_SEQ_DICT_INCOMPATIBILITY' \
-rRNA $rRNAfile " --job-name=$project -c 1 --mem-per-cpu=7000 -o ${alignDir}/rnaseqc.slurm >> ${alignDir}/commands.txt 

###############################################################################
## Second part: MultiQC                                                      ##
###############################################################################

## You may have to install the multiQC module in your lab library as we run it 
## from the BABS module collection

## Load module path                                                          ##
module use /camp/stp/babs/working/software/modules/all
## Load multiqc module                                                       ##
module load multiqc/0.9-2016b-Python-2.7.12

## If you have all FASTQC and RNASeQC and GATK output files downstream of the #
## run in that directory:
multiqc .

## In other cases you may have to run                                        ##
multiqc /path/to/first/fileset /path/to/second/fileset 

#end of file
