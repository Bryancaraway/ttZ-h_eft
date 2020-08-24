#!/bin/bash

# tar up data directory 
head=$(git rev-parse --show-cdup)
cdir=$PWD

cd $head/DeepSleep
cd ./data
rm -i data.tar.gz
# universal change to data dir
# and then tar contents 
tar cvfz data.tar.gz .
cd $cdir
# 
