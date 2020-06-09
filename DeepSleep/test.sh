#!/bin/bash

cd $PBS_O_WORKDIR
source /home/bcaraway/.bashrc
eval 'cd /home/bcaraway/DeepSleep/ZInvisible/Tools/DeepSleep/ ; conda activate mlenv;'
echo "$hi"
eval 'python test.py $hi'




# THIS WORKS FOR SIMPLE USE CASE qsub -v hi=1 test.py 
# Can parralize this with a python script to do something like:
# python runmanyjobs.py --> 
# qsub -o test1.out -e test1.err -v hi=1 test.py
# qsub -o test2.out -e test2.err -v hi=2 test.py
