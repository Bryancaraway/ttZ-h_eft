#!/bin/bash

# very simple job transfer script 
cd $PBS_O_WORKDIR
export X509_USER_PROXY=/cms/data/$USER/.x509_user_proxy
#source ~cmssoft/shrc
#source /home/$USER/.bashrc
echo "$target"
echo "$destination"
echo $X509_USER_PROXY
echo "xrdcp -f root://cmseos.fnal.gov//$target $destination"
eval "xrdcp -f root://cmseos.fnal.gov//$target $destination"
