#!/bin/bash
#$ -N pathway_analysis
#$ -cwd
#$ -V
#$ -j y
#$ -o logs/$JOB_NAME.$JOB_ID.log
#$ -l h_data=5G,h_rt=1:00:00
#$ -pe shared 2
#$ -M $USER@mail
#$ -m bea

# load the job environment:
. /u/local/Modules/default/init/modules.sh
module load R/4.3.0

Rscript test.R

#Rscript ../analysis/pathway/pathway_job.R
