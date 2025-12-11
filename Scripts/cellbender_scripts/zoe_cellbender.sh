#### cellbender.sh START ####
#!/bin/bash
#$ -cwd
# error = Merged with joblog
#$ -o log/joblog.$JOB_ID.$TASK_ID
#$ -j y
## Edit the line below as needed:
#$ -l gpu,RTX2080Ti,cuda=1,h_rt=2:00:00,h_data=10G
## Modify the parallel environment
# Email address to notify
#$ -M $USER@mail
# Notify when
#$ -m bea
#$ -t 1-4
#$ -tc 3 ##error if batch size is too large

# echo job info on joblog:
echo "Job $JOB_ID started on:   " `hostname -s`
echo "Job $JOB_ID started on:   " `date `
echo " "

source /u/local/Modules/default/init/modules.sh
module load gcc/4.9.3
module load anaconda3

conda activate cellbender

## Edit the line below as needed:
#conda activate /u/project/xyang123/shared/tools/miniconda3/envs/cellbender

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

## substitute the command to run your code
## in the two lines below:
cellbender remove-background \
    --input $fastq_dir/$sampleid/outs/raw_feature_bc_matrix.h5 \
    --output $sampleid.h5 \
    --cuda

# echo job info on joblog:
echo "Job $JOB_ID ended on:   " `hostname -s`
echo "Job $JOB_ID ended on:   " `date `
echo " "

#### cellbender.sh STOP ####
