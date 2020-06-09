
mkdir files

mkdir files/val 
mkdir files/test 
mkdir files/train 
mkdir files/train_overSample 
mkdir files/skim 
mkdir files/skim_kinemFit 
mkdir files/root  

mkdir files/skim_cnn
mkdir files/val_cnn 
mkdir files/test_cnn  
mkdir files/train_cnn  

mkdir pdf 

voms-proxy-init --rfc --voms cms

xrdcp -v root://cmseos.fnal.gov//store/user/bcaraway/skimAnaSamples/result_2016.root files/root/.
xrdcp -v root://cmseos.fnal.gov//store/user/bcaraway/skimAnaSamples/result_2017.root files/root/.
xrdcp -v root://cmseos.fnal.gov//store/user/bcaraway/skimAnaSamples/result_2018.root files/root/.
