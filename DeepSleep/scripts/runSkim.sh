#!/bin/bash

# very simple job submission script 

cd $PBS_O_WORKDIR
export OMP_NUM_THREADS=1
source /home/$USER/.bashrc

eval 'cd /home/$USER/ttZh_ana/DeepSleep/ ; conda activate ttxenv2.0;'
args="-s $sample -y $year"
if [ ! -z "$infile" ]; then
    args="$args -i $infile"
fi
if [ ! -z "$outfile" ]; then
    args="$args -o $outfile"
fi
if [ ! -z "$jec" ]; then
    args="$args -j $jec"
fi
if [ ! -z "$qsub" ]; then
    args="$args --qsub"
fi
if [ ! -z "$nopost" ]; then 
    args="$args --nopost"
fi
if [ ! -z "$noskim" ]; then
    args="$args --noskim"
fi
echo "$args"
eval 'python runSkim.py $args'
