'''
Script to create the postprocess architechture for ttzh analysis
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

job_output_dirs = '/cms/data/store/user/*/NanoAODv7/*'
assembly_dir    = '/cms/data/store/user/bcaraway/NanoAODv7/PostProcessed/'

json_filelist   = {'2016': f'{sys.path[1]}/data/sampleDas_nano_2016_v2.json',
                   '2017': f'{sys.path[1]}/data/sampleDas_nano_2017_v2.json',
                   '2018': f'{sys.path[1]}/data/sampleDas_nano_2018_v2.json'
               }

def get_finished_jobs(year):
    return glob(job_output_dirs+f'{year}*/*.root')

def find_and_transfer_files(year, json_file):
    samples = json.load(open(json_file))
    finished_jobs = get_finished_jobs(year)
    missing_files = {}
    for s in samples:
        missing_files[s] = {'files':[]}
        for sf in samples[s]['files']:
            sfile = sf.split('/')[-1].replace('.root','')#.replace('.root','_Skim[_]+\d+.root')
            #found_job = re.findall(rf'\w*{sfile}',' '.join(finished_jobs))
            found_job = re.findall(rf'/cms/data/store/user/\w*/NanoAODv7/\w*{year}/\w*{sfile}\w*.root',' '.join(finished_jobs))
            if not found_job:
                print(f'Missing: {sf} !!!')
            else:
                found_job = {os.stat(j).st_size:j for j in found_job}
                initial_loc = found_job[max(found_job.keys())]
                job_loc = re.search(rf'/\w*{year}/\w*{sfile}\w*.root', initial_loc).group()
                final_loc = assembly_dir+year+job_loc
                if os.path.exists(final_loc):
                    if os.stat(final_loc).st_size < os.stat(initial_loc).st_size:
                        print(f'New file better than existing one: {initial_loc}')
                        transfer_files(initial_loc, final_loc)
                else:
                    print(f'Transferring: {initial_loc}')
                    transfer_files(initial_loc, final_loc)
                

def transfer_files(initial_loc,final_loc):
    os.system(f"mkdir -p {'/'.join(final_loc.split('/')[:-1])} && cp {initial_loc} {final_loc}")

def main():
    for y, json in json_filelist.items():
        find_and_transfer_files(y,json)

if __name__ == '__main__':
    main()
