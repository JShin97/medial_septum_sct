##### SUBMIT GEOSUBMISSION SCRIPT START #######
#!/bin/bash
#$ -cwd
# error = Merged with joblog
#$ -o joblog.$JOB_ID
#$ -j y
#$ -V
# RESOURCE REQUESTED: Change resource with respect to your need.
#$ -l h_rt=24:00:00,h_data=100G
#$ -pe shared 2
# Email address to notify
#$ -M $USER@mail
# Notify when
#$ -m bea
#$ -N rsync

echo "Job $JOB_ID started on: " `hostname -s`
echo "Job $JOB_ID started on: " `date `
echo " "

module load rclone

origpath=/u/home/j/jshin/scratch/medial_septum_sct/GEOSubmission
#destpath=box:Hoffman2_data_backup/Current_members/Dan/data/human_retina
destpath=ftp_ncbi:uploads/julianaeshin@gmail.com_ktMpOOAU

#find $origpath -type f -size +15G > filesizeover15g.txt

#while IFS= read -r line
#do
#  echo split -b 15000m $line ${line}.part-
#  split -b 15000m $line ${line}.part-
#  rm $line
#done < filesizeover15g.txt

echo rclone sync -L $origpath $destpath
rclone sync -L $origpath $destpath
#rm filesizeover15g.txt


