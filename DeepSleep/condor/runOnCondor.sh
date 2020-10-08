#!/usr/bin/env bash
sample=$1
year=$2
estart=$3
estop=$4
tag=$5

# What sl are we running on?
lsb_release -a
#
echo "Transfering Env"
xrdcp -v root://cmseos.fnal.gov//store/user/bcaraway/AnaEnv/ttxenv.tar.gz .
mkdir -p ttxenv
tar -zxf ttxenv.tar.gz -C ttxenv
rm ttxenv.tar.gz
source ttxenv/bin/activate

tar -zxf ttzh.tar.gz
mkdir -p data
tar -zxf data.tar.gz -C data

pwd
ls -la

python runAna.py -s $sample -y $year -t $tag --estart $estart --estop $estop --condor
rm runAna.py

#do i need to move the output to ${_CONDOR_SCRATCH_DIR}?
