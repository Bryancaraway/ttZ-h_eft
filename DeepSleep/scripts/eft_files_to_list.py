## quick macro to make sure
## generated files from ken 
## are readable
## ---
## prevents crab jobs failing

from glob import glob
import uproot
from multiprocessing import Pool

### TTBB ###
#2018
#eft_dir  = '/cms/data/store/user/hatake/TopEFT/TTBB_13TeV_TopEFT_v12_2018/TTBB_13TeV_TopEFT_MINIAOD_v12_2018-v1/*.root'
#out_name = 'TTBB_13TeV_TopEFT_MINIAOD_v12_2018-v1.txt' 
### ---- ###

### TTBB JET ###
#2018
#eft_dir  = '/cms/data/store/user/hatake/TopEFT/TTBBJet_13TeV_TopEFT_v12_gen/TTBBJet_13TeV_TopEFT_MINIAOD_v12_gen/step_miniaodsim_*.root'
#out_name = 'TTBBJet_13TeV_TopEFT_MINIAOD_v12_2018-v1.txt' 
#
#eft_dir = '/cms/data/store/user/bcaraway/NanoAODv7/PreProcessed/2018/TTBBJet_13TeV_TopEFT_MINIAOD_v12_2018-v1/NanoHRT_v12/210401_183458/0000/prod2018MC_v7_NANO_*.root'
#out_name = 'TTbb_EFT_2018.txt'
### ---- ###

### TTZ  ###
#2016
#eft_dir  = '/cms/data/store/user/hatake/TopEFT/TTZJet_13TeV_TopEFT_v20_2016/TTZJet_13TeV_TopEFT_MINIAOD_v20_2016-v1/step_miniaodsim_*.root'
#out_name = 'TTZJet_13TeV_TopEFT_MINIAOD_2016_v20-v1.txt'
#2017
#eft_dir  = '/cms/data/store/user/hatake/TopEFT/TTZJet_13TeV_TopEFT_v20_2017/TTZJet_13TeV_TopEFT_MINIAOD_v20_2017-v1/step_miniaodsim_*.root'
#out_name = 'TTZJet_13TeV_TopEFT_MINIAOD_2017_v20-v1.txt'
#2018
#eft_dir  = '/cms/data/store/user/hatake/TopEFT/TTZJet_13TeV_TopEFT_v20/TTZJet_13TeV_TopEFT_MINIAOD_v20-v1/step_miniaodsim_*.root'
#out_name = 'TTZJet_13TeV_TopEFT_MINIAOD_2018_v20-v1.txt'
### ---- ###

### TTH  ###
#2016
#eft_dir  = '/cms/data/store/user/hatake/TopEFT/TTHJet_13TeV_TopEFT_v11_2016/TTHJet_13TeV_TopEFT_MINIAOD_v11_2016-v2/*.root'
#out_name = 'TTHJet_13TeV_TopEFT_MINIAOD_2016_v11-v2.txt'
#2017
#eft_dir  = '/cms/data/store/user/hatake/TopEFT/TTHJet_13TeV_TopEFT_v11_2017/TTHJet_13TeV_TopEFT_2017_MINIAOD_v11_2017-v1/*.root'
#out_name = 'TTHJet_13TeV_TopEFT_MINIAOD_2017_v11-v1.txt'
#2018
#eft_dir  = '/cms/data/store/user/hatake/TopEFT/TTHJet_13TeV_TopEFT_v11hs/TTHJet_13TeV_TopEFT_MINIAOD_v11hs-v3/*.root'
#out_name = 'TTHJet_13TeV_TopEFT_MINIAOD_2018_v11hs-v3.txt'
### ---- ###

#eft_dir = '/cms/data/store/user/bcaraway/NanoAODv7/PreProcessed/2018/TTZJet_13TeV_TopEFT_MINIAOD_2018_v20-v1/NanoHRT_v20/210418_201500/0000/prod2018MC_v7_NANO_*.root'
#eft_dir = '/cms/data/store/user/bcaraway/NanoAODv7/PreProcessed/2017/TTZJet_13TeV_TopEFT_MINIAOD_2017_v20-v1/NanoHRT_v20/210418_201330/0000/prod2017MC_v7_NANO_*.root'
eft_dir = '/cms/data/store/user/bcaraway/NanoAODv7/PostProcessed/2017/TTZ_EFT_2017/prod2017MC_v7_NANO_*.root'
out_name = 'ttz_eft_2017.txt'

def main():
    eft_files = glob(eft_dir)
    pool = Pool()
    out_list = pool.map(__worker, eft_files)

    with open(out_name,'w') as o_f:
        o_f.writelines(out_list)

def __worker(file_):
    try:
        with uproot.open(file_):
            return file_.replace('/cms/data/','root://kodiak-se.baylor.edu//')+'\n'
    #except IndexError: # can also get ValueError
    except:
        print('bad file',file_)
        return ''

if __name__ == '__main__':
    main()
