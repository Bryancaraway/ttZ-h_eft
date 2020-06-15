#!/bin/bash
# Setup file system for analysis

mkdir files

YEARS=(2016 2017 2018)
for y in ${YEARS[@]}; do
    mkdir files/$y 
    mkdir files/$y/mc_files
    mkdir files/$y/data_files
done

mkdir pdf 
mkdir log
mkdir NN_output

