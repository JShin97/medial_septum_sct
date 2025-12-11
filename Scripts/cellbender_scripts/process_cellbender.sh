##### SUBMIT QIIME2 SCRIPT START #######
#!/bin/bash
#$ -cwd
# error = Merged with joblog
#$ -o log/joblog.$JOB_ID
#$ -j y
#$ -V
# RESOURCE REQUESTED: Change resource with respect to your need.
#$ -l h_rt=2:00:00,h_data=30
#$ -pe shared 4
# Email address to notify
#$ -M $USER@mail
# Notify when
#$ -m bea
#$ -N cellbender
#$ -t 3-24

source /u/local/Modules/default/init/modules.sh
module load gcc/4.9.3

module load python
module load anaconda3

conda activate cellbender

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
fastq_dir=$root_dir/Data/Raw/Counts
output_dir=$root_dir/Data/Raw/Cellbender_counts2

#mkdir $output_dir/$sampleid

cd $output_dir/$sampleid

######################################################
#   About different parameters in cellbender remove-background   #
######################################################
# --input : raw_feature_bc_matrix
# --output : output file name 
# --expected-cells : 10000 
# --total-droplets-included : 15000
# --fpr : 0.01 
# --epochs : 150
# --learning-rate

ptrepack --complevel 5 ${sampleid}_filtered.h5:/matrix ${sampleid}_filtered_seurat.h5:/matrix

echo "Job $JOB_ID ended on: " `hostname -s`
echo "Job $JOB_ID ended on: " `date `
echo " "