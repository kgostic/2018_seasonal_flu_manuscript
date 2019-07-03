#!/bin/bash
#$ -cwd
#$ -l h_rt=04:00:00,h_data=20G
#$ -N INSIGHT-profiles
#$ -o /u/home/k/kgostic/seasonal_flu/2017_INSIGHT/.messages
#$ -e /u/home/k/kgostic/seasonal_flu/2017_INSIGHT/.errors
#$ -m abe
#$ -M kgostic
#$ -t 01-32

## Load R
source /u/local/Modules/default/init/modules.sh
module load R/3.4.0 # or whatever version you want (you can see options using 'module avail')

## Load R script and output path
script=/u/home/k/kgostic/seasonal_flu/2017_INSIGHT/02-calculate-likelihood-profiles.R
outpath=/u/home/k/kgostic/seasonal_flu/2017_INSIGHT/cluster_outputs/INSIGHT_profiles.txt

##  run your script:
Rscript $script --outpath $outpath --var $SGE_TASK_ID
