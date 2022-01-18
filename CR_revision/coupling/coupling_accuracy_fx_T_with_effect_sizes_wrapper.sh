#!/bin/bash
export QUEUE=taki
echo LSB_JOB_REPORT_MAIL=N >> ~/.bashrc
#source /project/bbl_projects/apps/default_modules.sh
#conda activate rstudio
Rscript coupling_accuracy_fx_T_with_effect_sizes.R
