#!/bin/bash

# very simple job submission script 

cd $PBS_O_WORKDIR
export OMP_NUM_THREADS=1
source /home/$USER/.bashrc
eval 'cd /home/$USER/ttZh_ana/DeepSleep/ ; conda activate ttxenv2.0;'
args="-s $sample -y $year"
if [ ! -z "$jec" ]; then
    args="$args -j $jec"
fi
echo "$args"
eval 'python runAna.py $args'


