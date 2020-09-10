#!/bin/bash

# very simple job submission script 

cd $PBS_O_WORKDIR
source /home/$USER/.bashrc
eval 'cd /home/$USER/ttZh_ana/DeepSleep/ ; conda activate ttxenv;'
echo "$sample"
echo "$year"
echo "$jec"
echo "$jjec"
eval 'python runAna.py -s $sample -y $year --jjec $jjec --jec $jec'
