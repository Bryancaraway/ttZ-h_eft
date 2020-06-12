#!/bin/bash

cd $PBS_O_WORKDIR
source /home/bcaraway/.bashrc
eval 'cd /home/bcaraway/DeepSleep_v2/DeepSleep/ ; conda activate mlenv;'
echo "$sample"
eval 'python runAna.py $sample'




# THIS WORKS FOR SIMPLE USE CASE qsub -v hi=1 test.py 
# Can parralize this with a python script to do something like:
# python runmanyjobs.py --> 
# qsub -o test1.out -e test1.err -v hi=1 test.sh
# qsub -o test2.out -e test2.err -v hi=2 test.sh
