### submit_rstudio.sh START ###
#!/bin/bash
#$ -cwd
#$ -o  logs/joblog.$JOB_ID
#$ -j y
#  Resources requested
#  PLEASE CHANGE THE RESOURCES REQUESTED AS NEEDED:
#$ -l h_data=100G,h_rt=24:00:00
#  PLEASE CHANGE THE NUMBER OF CORES REQUESTED AS NEEDED:
#$ -pe shared 1
#  Email address to notify
#$ -M $USER@mail
#$ -m bea
#$ -t 1

#
# Output job info on joblog file:
#
echo " "
echo "Job myscript, ID no. $JOB_ID started on:   "` hostname -s `
echo "Job myscript, ID no. $JOB_ID started at:   "` date `
echo " "

#
# Set up job environment:
#
. /u/local/Modules/default/init/modules.sh
module load apptainer

# Export R Libs directory
## NOTE YOU MAY NEED TO UPDATE THE DEFINITION OF R_LIBS_USER TO
## MATCH THE VERSION OF R YOU ARE USING (IF NOT USING h2-rstudio_4.1.0.sif)
export R_LIBS_USER=/u/project/xyang123/jshin/apps/R/4.3.0

# Set trait:
trait=$(sed -n ${SGE_TASK_ID}p gwas_reprocessed_post_mdf.txt)


# Run the job:
#
## REMEMBER TO SUBSTITUTE THE NAME OF THE R SCRIPT YOU INTEND TO RUN IN THE TWO LINES BELOW &
## TO CHANGE THE RSTUDIO SERVER CONTAINER IF A DIFFERENT VERSION OF R IS NEEDED (TO SEE AVAILABLE 
## VERSIONS OF R ISSUE: module load apptainer; ls $H2_CONTAINER_LOC/h2-rstudio*sif)
echo "/usr/bin/time -v apptainer exec $H2_CONTAINER_LOC/h2-rstudio_4.3.0.sif R CMD BATCH --no-save --no-restore MSEA_correctmapping_jobarray.R $trait logs/output.$JOB_ID"
/usr/bin/time -v apptainer  exec $H2_CONTAINER_LOC/h2-rstudio_4.3.0.sif R CMD BATCH --no-save --no-restore MSEA_correctmapping_jobarray.R ${trait} logs/output.$JOB_ID

echo " "
echo "Job myscript, ID no. $JOB_ID finished at:  "` date `
echo " "
### submit_rstudio.sh STOP ###
