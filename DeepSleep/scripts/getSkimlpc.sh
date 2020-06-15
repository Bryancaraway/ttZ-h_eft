#!/bin/bash


#xrdcp -v root://cmseos.fnal.gov//store/user/bcaraway/skimAnaSamples/result_2016.root files/root/.
#xrdcp -v root://cmseos.fnal.gov//store/user/bcaraway/skimAnaSamples/result_2017.root files/root/.
#xrdcp -v root://cmseos.fnal.gov//store/user/bcaraway/skimAnaSamples/result_2018.root files/root/.

file=$1
file_dir=/cms/data/store/user/ttxeft/Skim_nanoAOD
lpc_dir=store/user/bcaraway/skimAnaSamples

if [ -z $file]; then
    echo "Must proved root file name to be transfered from lpc to kodiak."
    echo "Example: Data_2017.root, MC_2017.root"
    exit 1
fi

#voms-proxy-init --rfc --voms cms

echo "Moving file $file from lpc to kodiak directory: $file_dir"
echo "A backup of the previous file version (if any) can be found here: $file_dir/backup/"

mv -uv $file_dir/$file $file_dir/backup/.
rm  $file_dir/$file
echo "xrdcp -v root://cmseos.fnal.gov//$lpc_dir/$file $file_dir/$file"
xrdcp -v root://cmseos.fnal.gov//$lpc_dir/$file $file_dir/$file
