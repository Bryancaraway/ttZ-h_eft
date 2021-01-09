'''
Script to add file for das entries in json output
'''
import os
import signal
import sys
import re
import json
from glob import glob
import uproot
import time
import numpy as np
import subprocess as sb
from multiprocessing import Pool

job_output_dirs = '/cms/data/store/user/*/NanoAODv7/*'


def find_missing_files(json_file):
    year = re.search(r"201\d" , json_file).group()
    samples = json.load(open(json_file))
    finished_jobs = glob(job_output_dirs+f'{year}*/*.root')
    pp_finished_jobs = glob('/cms/data/store/user/bcaraway/NanoAODv7/PostProcessed/*/*'+f'{year}*/*.root')
    missing_files = {}
    for s in samples:
        missing_files[s] = {'files':[]}
        for sf in samples[s]['files']:
            sfile = sf.split('/')[-1].replace('.root','')#.replace('.root','_Skim[_]+\d+.root')
            #found_job = re.findall(rf'\w*{sfile}',' '.join(finished_jobs))
            #found_job = re.findall(rf'/cms/data/store/user/\w*/NanoAODv7/\w*{year}\w*/\w*{sfile}\w*.root',' '.join(finished_jobs))
            found_job = re.findall(rf'/cms/data/store/user/bcaraway/NanoAODv7/PostProcessed/{year}/\w*{year}\w*/\w*{sfile}\w*.root',' '.join(pp_finished_jobs))
            if not found_job:
                found_job = re.findall(rf'/cms/data/store/user/\w*/NanoAODv7/\w*{year}\w*/\w*{sfile}\w*.root',' '.join(finished_jobs)) # look in kens/ others area for files
            #check if incomplete
            job_incomplete = False
            if found_job:
                for i_job in range(len(found_job)):
                    job = found_job[i_job]
                    try:
                        with uproot.open(job) as _:
                            break
                    except IndexError:
                        if i_job < len(found_job)-1:
                            continue
                        else:
                            print(f'All bad in: {found_job}')
                            if not attempt_transfer(found_job,year,sfile):
                                job_incomplete = True # end of for loop
                    
            if not found_job or job_incomplete:
                if not job_incomplete: print(f'Missing: {sf}!!!!!!!')
                missing_files[s]['files'].append(sf)
    out_file = '/'.join(json_file.split('/')[:-1])+f'/rerun_jobs_{year}.json'
    print(out_file)
    with open(out_file, 'w') as jf:
        json.dump(missing_files, jf, indent=4)
    

def attempt_transfer(potential_names,year,sfile):
    ### only this type ###
    #figure out target directory
    #if True: return False
    ### kodiak not working
    destination=f'/cms/data/store/user/bcaraway/NanoAODv7/PostProcessed/{year}/'
    target = ''
    #if year == '2016' or year == '2018':
    #    target = 'store/user/bcaraway/NanoAODv7/'
    #else:
    #target = 'store/user/hatake/NanoAODv7/'
    
    #with open('bryans_nanoaod_files.txt') as target_file:
    with open('kens_nanoaod_files.txt') as target_file:
        for line in target_file.readlines():
            if sfile in line:
                print(f'File found at: {line} ') 
                l = line.replace('\n','').replace('root://cmseos.fnal.gov//','')
                target = l
                destination += '/'.join(l.split('/')[-2:])
                break
            
    #for name in potential_names:
    #    if target in name:
    #        target_cand = name.replace('/cms/data/','').replace('hatake','hcal_upgrade/noreplica/hatake')
    #        try:
    #            with uproot.open(f'root://cmseos.fnal.gov//{target_cand}') as _:
    #                print(f'File found at: {target_cand} ')
    #                target = target_cand
    #        except:
    #            continue
    #        destination += '/'.join(name.split('/')[-2:])
    # see if file is remote
    if '.root' not in target:
        return False
    #
    #if 'TTbb_' not in target: return False
    # make sure to wait before trying to transfer
    num_jobs_running = lambda: int(sb.check_output('qstat -u $USER | grep transfer | wc -l', shell=True).decode())
    #num_jobs_running = lambda: int(sb.check_output('jobs | wc -l', shell=True).decode())
    print(num_jobs_running())
    while num_jobs_running() >= 30:
        time.sleep(30)

    out_name = destination.split('/')[-1]
    args = f'-v target={target},destination={destination}'
    #print(f"xrdcp -f root://cmseos.fnal.gov//{target} {destination} &>/dev/null &")
    #command = f"eval 'xrdcp -f root://cmseos.fnal.gov//{target} {destination} &>/dev/null &'"
    command = f'cp $X509_USER_PROXY /cms/data/$USER/.x509_user_proxy ;qsub -l nodes=1:ppn=1 -o {year}/{out_name}.out -e {year}/{out_name}.err {args} transfer.sh'
    #
    print(command)
    os.system(command)
    return True
    

if __name__ == '__main__':
    if not os.path.exists(sys.argv[1]) : raise NameError(sys.argv[1]) 
    i_file = sys.argv[1]
    #das_from_cfg(cfg_file)
    find_missing_files(i_file)
