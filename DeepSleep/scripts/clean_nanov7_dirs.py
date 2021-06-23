'''
Script to clean postprocess directiories
'''
import os
import signal
import sys
import re
import json
from glob import glob
import uproot
import numpy as np
import subprocess as sb
sys.path.insert(1, sb.check_output('echo $(git rev-parse --show-cdup)', shell=True).decode().strip('\n')+'DeepSleep/')
from multiprocessing import Pool

job_output_dirs = '/cms/data/store/user/*/NanoAODv7/'
assembly_dir    = '/cms/data/store/user/bcaraway/NanoAODv7/PostProcessed/'

json_filelist   = {'2016': f'{sys.path[1]}/data/sampleDas_nano_2016_v2.json',
                   '2017': f'{sys.path[1]}/data/sampleDas_nano_2017_v2.json',
                   '2018': f'{sys.path[1]}/data/sampleDas_nano_2018_v2.json'        
       }

def find_and_remove_bad_files(year, json_file):
    samples = json.load(open(json_file))
    for s in samples:
        for sf in samples[s]['files']:
            sfile = sf.split('/')[-1].replace('.root','')
            files = set(glob(assembly_dir+year+f'/{s}/{sfile}*.root'))
            #print(files)
            if not files: continue
            for f_ in files:
                if not goodfile(f_):
                    print(f'Removing bad file at: {f_}')
                    os.system(f'rm {f_}')
            files = set(glob(assembly_dir+year+f'/{s}/{sfile}*.root'))
            if not files: continue
            if len(files) > 1 :
                print(files)
            
            
            
        

def goodfile(f_):
    try:
        with uproot.open(f_) as _:
            return True
    except IndexError:
        return False
        

if __name__ == '__main__':
    for y, js in json_filelist.items():
        find_and_remove_bad_files(y,js)
