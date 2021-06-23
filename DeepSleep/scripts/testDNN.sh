#!/bin/bash

# very simple job submission script 
# to test dnn model archs

cd $PBS_O_WORKDIR
#export OMP_NUM_THREADS=1
source /home/$USER/.bashrc
eval 'cd /home/$USER/ttZh_ana/DeepSleep/ ; conda activate ttxenv2.0;'
args="--mode exe --inputjson $json --jnum $jobnum"

echo "$args"
eval 'python testDNN_models.py $args'
