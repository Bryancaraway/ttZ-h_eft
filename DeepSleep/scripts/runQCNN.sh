#!/bin/bash

# very simple job submission script 

cd $PBS_O_WORKDIR
export OMP_NUM_THREADS=1
source /home/$USER/.bashrc

eval 'cd /home/$USER/ttZh_ana/DeepSleep/ ; conda activate ttxenv2.0;'

args=""
if [ ! -z "$nn_input" ]; then
    args="$args $nn_input"
fi
echo "$args"
eval 'python qcDatacard.py $args True'
