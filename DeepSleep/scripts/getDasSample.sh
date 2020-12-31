#!/bin/bash

das_sample=$1

source /cvmfs/cms.cern.ch/cmsset_default.sh
voms-proxy-init --rfc --voms cms
source /cvmfs/cms.cern.ch/rucio/setup.sh
export RUCIO_ACCOUNT=`whoami`

eval 'rucio add-rule cms:$das_sample 1 T3_US_Baylor'
