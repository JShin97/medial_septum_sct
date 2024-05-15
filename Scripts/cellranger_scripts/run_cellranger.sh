##### SUBMIT QIIME2 SCRIPT START #######
#!/bin/bash
#$ -cwd
# error = Merged with joblog
#$ -o log/joblog.$JOB_ID
#$ -j y
#$ -V
# RESOURCE REQUESTED: Change resource with respect to your need.
#$ -l h_rt=12:00:00,h_data=30,highp
#$ -pe shared 4
# Email address to notify
#$ -M $USER@mail
# Notify when
#$ -m bea
#$ -N cellrangerCount
#$ -t 1 # this number can change depends on how many samples you have (e.g. if you have 30 then 1-30)

source /u/local/Modules/default/init/modules.sh
module load gcc/4.9.3

module load cellranger
module load bcl2fastq

echo "Job $JOB_ID started on: " `hostname -s`
echo "Job $JOB_ID started on: " `date `
echo " "

# $SGE_TASK_ID is number you get from each run of the script given by -t parameter above. It will sequentially start from 1 to the end 
#For this particular script, you want your data to be under directory named after your sampleid and store that sample id to sampleid.txt
#Each line of sampleid.txt will have one sample id
#Command below will get you a sample id with given $SGE_TASK_ID number.
sampleid=$(sed -n ${SGE_TASK_ID}p sampleid.txt)

##change paths as needed
root_dir=/u/home/j/jshin/scratch/medial_septum_sct
fastq_dir=$root_dir/Data/Raw/Fastq
output_dir=$root_dir/Data/Raw/Counts

cd $output_dir

######################################################
#   About different parameters in cellranger count   #
######################################################
# --id : your sample id from sampleid.txt
# --transcriptome: path to your reference
# --fastqs: path to your fastq files
# --sample: sample name 
# --localcores: how many cpu to allocate
# --localmem: how many memory to allocate
# --include-introns false: only true if you are working with single nuclei data
echo "cellranger count --id $sampleid --transcriptome=/u/project/xyang123/shared/tools/10XGenomics/References/refdata-cellranger-mm10-3.0.0 --fastqs=$fastq_dir --sample=$sampleid --localcores=8 --localmem=64 --include-introns false"

cellranger count --id=$sampleid \
                --sample=$sampleid \
                --transcriptome=/u/project/xyang123/shared/tools/10XGenomics/References/refdata-cellranger-mm10-3.0.0 \
                --fastqs=$fastq_dir/XY-TriSeptum-1 \
                --project=MS_SCT \
                --localcores=8 \
                --localmem=4 \
                --include-introns=true


echo "Job $JOB_ID ended on: " `hostname -s`
echo "Job $JOB_ID ended on: " `date `
echo " "
