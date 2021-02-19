## quick macro to make sure
## generated files from ken 
## are readable
## ---
## prevents crab jobs failing

from glob import glob
import uproot
from multiprocessing import Pool

#eft_dir = '/cms/data/store/user/hatake/TopEFT/TTBB_13TeV_TopEFT_cbW_2018/TTBB_13TeV_TopEFT_MINIAOD_cbW_2018-v1/*.root'
eft_dir  = '/cms/data/store/user/hatake/TopEFT/TTBB_13TeV_TopEFT_v12_2018/TTBB_13TeV_TopEFT_MINIAOD_v12_2018-v1/*.root'
out_name = 'TTBBJet_13TeV_TopEFT_MINIAOD_2018-v12.txt'

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
    except IndexError:
        print('bad file',file_)
        return ''

if __name__ == '__main__':
    main()
