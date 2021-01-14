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
                   #'2018': f'{sys.path[1]}/data/sampleDas_nano_2018_v2.json'
               }

def get_finished_jobs(year):
    return glob(job_output_dirs+f'{year}*/*.root')

def find_and_transfer_files(year, json_file):
    samples = json.load(open(json_file))
    finished_jobs = get_finished_jobs(year)
    pp_finished_jobs = glob('/cms/data/store/user/bcaraway/NanoAODv7/PostProcessed/*/*'+f'{year}*/*.root')
    missing_files = {}
    for s in samples:
        missing_files[s] = {'files':[]}
        for sf in samples[s]['files']:
            sfile = sf.split('/')[-1].replace('.root','')#.replace('.root','_Skim[_]+\d+.root')
            #found_job = re.findall(rf'\w*{sfile}',' '.join(finished_jobs))
            #print(sf)
            #found_job = re.findall(rf'/cms/data/store/user/\w*/NanoAODv7/\w*{year}\w*/\w*{sfile}\w*.root',' '.join(finished_jobs))
            done_job = re.findall(rf'/cms/data/store/user/bcaraway/NanoAODv7/PostProcessed/{year}/\w*{year}\w*/\w*{sfile}\w*.root',' '.join(pp_finished_jobs))
            found_job = []
            if done_job:
                try:
                    with uproot.open(done_job[0]) as _:
                        continue
                except:
                    found_job = re.findall(rf'/cms/data/store/user/\w*/NanoAODv7/\w*{year}\w*/\w*{sfile}\w*.root',' '.join(finished_jobs))
                    pass
            else:
                print('Not in PostProcessed area:',sf)
                continue
                found_job = re.findall(rf'/cms/data/store/user/\w*/NanoAODv7/\w*{year}\w*/\w*{sfile}\w*.root',' '.join(finished_jobs)) # look in kens/ others area for files
            if not found_job:
                print(f'Missing: {sf} !!!')
            else:
                found_job = {os.stat(j).st_size:j for j in found_job}
                initial_loc = found_job[max(found_job.keys())]
                job_loc = re.search(rf'/\w*{year}\w*/\w*{sfile}\w*.root', initial_loc).group()
                final_loc = re.findall(rf'/cms/data/store/user/bcaraway/NanoAODv7/PostProcessed/{year}/\w*{year}\w*/\w*{sfile}\w*.root',' '.join(pp_finished_jobs)) 
                #if os.path.exists(final_loc):
                if final_loc:
                    final_loc = final_loc[0]
                    try:
                        with uproot.open(final_loc) as _ :
                            continue
                    except IndexError:
                        print(f"Bad file at {final_loc}")
                        pass
                    if os.stat(final_loc).st_size < os.stat(initial_loc).st_size:
                        print(f'New file better than existing one: {initial_loc}')
                        transfer_files(initial_loc, final_loc)
                else:
                    final_loc = assembly_dir+year+job_loc
                    print(f'Transferring: {initial_loc}')
                    transfer_files(initial_loc, final_loc)
                

def transfer_files(initial_loc,final_loc):
    try:
        with uproot.open(initial_loc) as _:
            os.system(f"mkdir -p {'/'.join(final_loc.split('/')[:-1])} && cp {initial_loc} {final_loc}")
    except IndexError:
        print(f'{initial_loc} cannot be opened!!!!')

def main():
    for y, json in json_filelist.items():
        if len(sys.argv) > 1:
            if y != str(sys.argv[1]): continue 
        find_and_transfer_files(y,json)

if __name__ == '__main__':
    main()
