#!/bin/bash

# very simple job submission script 

cd $PBS_O_WORKDIR
source /home/$USER/.bashrc
eval 'cd /home/$USER/ttZh_analysis/DeepSleep/ ; conda activate ttxenv;'
echo "$sample"
eval 'python runAna.py $sample'

