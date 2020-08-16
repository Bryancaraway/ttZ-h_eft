#!/bin/bash
# Setup file system for analysis

mkdir -p files

YEARS=(2016 2017 2018)
for y in ${YEARS[@]}; do
    mkdir -p files/$y 
    mkdir -p files/$y/mc_files
    mkdir -p files/$y/data_files
done

mkdir -p pdf 
mkdir -p log
#mkdir NN_files
#cp /cms/data/store/user/ttxeft/NN_files/* NN_files/
#
tar -zxvf data.tar.gz