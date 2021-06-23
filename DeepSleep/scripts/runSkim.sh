#!/bin/bash

# very simple job submission script 

cd $PBS_O_WORKDIR
export OMP_NUM_THREADS=1
source /home/$USER/.bashrc

eval 'cd /home/$USER/ttZh_ana/DeepSleep/ ; conda activate ttxenv2.0;'

if [ -z "$PBS_ARRAY_INDEX" ]; then
    export PBS_ARRAY_INDEX=0
fi
echo "$PBS_ARRAY_INDEX"
echo "$jec"
args="-s $sample -y $year"
#if [ ! -z "$infile" ]; then
#    args="$args -i $infile"
if [ ! -z "$infile" ]; then
    args="$args -i $infile$PBS_ARRAY_INDEX.list"
fi
if [ ! -z "$outfile" ]; then
    args="$args -o ${outfile}_${PBS_ARRAY_INDEX}${tag}.pkl"
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
if [ ! -z "$is4trig" ]; then
    args="$args --is4trig"
fi
echo "$args"
eval 'python runSkim.py $args'
