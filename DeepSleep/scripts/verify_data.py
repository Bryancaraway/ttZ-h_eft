
import sys
import os
import subprocess as sb
if __name__ == '__main__':
    sys.path.insert(1, sb.check_output('echo $(git rev-parse --show-cdup)', shell=True).decode().strip('\n')+'DeepSleep/')
import concurrent.futures
import config.ana_cff as cfg
import uproot
import multiprocessing
import re
import json
import glob
import numpy as np
executor = concurrent.futures.ThreadPoolExecutor()
#/cvmfs/cms.cern.ch/common/dasgoclient -query="dataset=/SingleElectron/Run2016B-02Apr2020_ver2-v1/NANOAOD" -json

dasgo = '/cvmfs/cms.cern.ch/common/dasgoclient'

def worker(roo):
    if 'kodiak' in roo:
        roo = roo.replace('root://kodiak-se.baylor.edu//','/cms/data/')
    t=uproot.open(roo)['Events']
    n_entries = len(t)
    return n_entries 
    
 
def process(d_sample,das,y,pool):
    # get n_entries according to Das
    das_out = sb.check_output(f"{dasgo} -query='dataset={das}' -json", shell=True).decode()
    das_entries = int(re.search(r"\"nevents\":\d*",das_out).group().split(':')[-1])
    # get pre-skimmed file list
    postpdir = f'{cfg.postproc_dir}/{y}/{d_sample}/'
    f_list = glob.glob(postpdir+'*.root')
    results = pool.map(worker, f_list)
    print(f"Data validation for  : {d_sample} --> {'GOOD' if das_entries == sum(results) else 'BAD'}")
    print(f" # entries in DAS      : {das_entries}")
    print(f" # entries in pre-Skim : {sum(results)}")
    # process results

    
        

def main(cfg_dir):
    pool = multiprocessing.Pool() # have to define this here
    for y in cfg.Years:
        jfile = f'{cfg_dir}/sampleDas_nano_{y}_v2.json'
        print(jfile)
        with open(jfile) as jf:
            s_dict = json.load(jf)
            data_sets = re.findall(r'Data\w*', ' '.join(s_dict.keys()))
            for d_sample in data_sets:
                process(d_sample,s_dict[d_sample]['DAS'][0],y,pool)


        

if __name__ == '__main__':

    cfg_dir = cfg.dataDir
    #
    main(cfg_dir)
