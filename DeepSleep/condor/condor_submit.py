# run full analysis on HTCondor batch system
# workflow: 
# create job_dir, create custom .jdl 
# use case: python submitjobs.py
#
import sys
import os
import argparse
if __name__ == '__main__':
    import subprocess as sb
    sys.path.insert(1, sb.check_output('echo $(git rev-parse --show-cdup)', shell=True).decode().strip('\n')+'DeepSleep/')
    # this should mean that wherever the script is executed, it should work
    os.system(f'cd {sys.path[1]}/condor')
import config.ana_cff as cfg
from datetime import datetime
import numpy as np
import uproot
#
# global variables # ===========================================
ana_env_dir = 'root://cmseos.fnal.gov//store/user/bcaraway/AnaEnv/'
ana_env     = 'ttxenv.tar.gz'
data        = 'data.tar.gz'
# ==============================================================

# add parsargs at some point for minibatch
parser = argparse.ArgumentParser(description='Run analysis over many samples and years using batch system')
parser.add_argument('-s', dest='samples', type=str, choices=['all','mc','allmc','data','test'], 
                    required=True, help='analyze all, mc only, or data only')
parser.add_argument('-y', dest='year', type=str, choices=cfg.Years+['all'],
                    required=True, help='year or all years')
#parser.add_argument('-j', dest='jec',     type=str, required=False, help='Run with specified jec variation', choices=['JESUp','JESDown','JERUp','JERDown','all'], default='')
#parser.add_argument('-t', dest='tag',     type=str, required=False, help='Optional tag to add to output file', default='')
args = parser.parse_args()

class condorSubmit ():

    def __init__(self):
        sample_dict = {'all':cfg.All_MC+cfg.Data_samples,
                       'mc': cfg.MC_pow,
                       'allmc':cfg.All_MC,
                       'data':cfg.Data_samples,
                       'test':['TriBoson']}
        #
        job_dir = f"job_{datetime.now().strftime('%d%m_%H%M')}"
        os.system(f'mkdir -p {job_dir}')
        os.system(f'mkdir -p {job_dir}/logs')
        # create tarball of analysis scripts
        os.chdir('../')
        os.system(f'tar -zcf ttzh.tar.gz config lib modules runAna.py')
        os.system(f'mv ttzh.tar.gz condor/{job_dir}')
        os.chdir('condor')
        # os.system(f'cp {ana_env} {job_dir}')
        os.system(f'cp ../data/{data} {job_dir}')
        os.system(f'cp runOnCondor.sh {job_dir}')
        
        self.samples = sample_dict[args.samples]
        years = [args.year] if args.year != 'all' else cfg.Years
        # create .jdl file
        jdl = open(f'{job_dir}/{job_dir}.jdl', 'w')
        jdl.writelines(['universe = vanilla\n',
                        'Executable = runOnCondor.sh\n',
                        '\n',
                        'Should_transfer_files = YES\n',
                        'WhenToTransferOutput = ON_EXIT\n',
                        'notification = never\n',
                        '\n',
                        f'Transfer_Input_Files = ttzh.tar.gz,{data}\n',
                        '\n',
                        'notify_user = ${LOGNAME}@FNAL.GOV\n',
                        'x509userproxy = $ENV(X509_USER_PROXY)\n',
                        '\n',
                        #'request_cpus = 4\n',
                        #'request_memory = 8000\n',
                        '\n'])
               
        # create a batch map for minibatch (if needed)
        self.job_lines = []
        for year in years:
            batch_map = self.batch_map(year)
            for sample in self.samples:
                if batch_map[sample][-1] > 1:
                    self.miniBatch(sample, year, batch_map[sample])
                else:
                    self.job_lines.append(f'Arguments = {sample} {year} 0 {batch_map[sample][1]} 0')
                    self.job_lines.append(f'Output = logs/{sample}_{year}_$(Process).stdout')
                    self.job_lines.append(f'Error = logs/{sample}_{year}_$(Process).stderr')
                    self.job_lines.append(f'Log = logs/{sample}_{year}_$(Process).log')
                    self.job_lines.append('Queue')
                    #
                #
            #
        #
        jdl.writelines([line+'\n' for line in self.job_lines])
        jdl.close()
        os.chdir(f'{job_dir}')
        os.system(f'condor_submit {job_dir}.jdl')    

    def miniBatch(self, sample, year, info):
        start_stops = np.linspace(info[0],info[1],info[2]+1).astype('int') # have to add one cuase numpy
        for i,s in enumerate(start_stops[:-1]): # iterate till next to last
            start = start_stops[i]
            stop  = start_stops[i+1]
            self.job_lines.append(f'Arguments = {sample} {year} {start} {stop} {i}')
            self.job_lines.append(f'Output = logs/{sample}_{year}_$(Process).stdout')
            self.job_lines.append(f'Error = logs/{sample}_{year}_$(Process).stderr')
            self.job_lines.append(f'Log = logs/{sample}_{year}_$(Process).log')
            self.job_lines.append('Queue')
            
        
    def batch_map(self, year):
        # lets say if num entries is > 300k split into even jobs of 300k max
        max_events = 200000 # bad hardcode
        batch_map = {} # this will ultimately have the format [range(0,numentries), n_files] 
        # first get numentries from files
        # do for MC and Data file for that year
        for f_type in ['MC','Data']:
            with uproot.open(f'{cfg.file_path}/{f_type}_{year}.root') as roofile:
                trees = roofile.get(cfg.tree_dir+'_bb')
                for sample in self.samples:
                    if sample in trees:
                        n_events = trees.get(sample).numentries
                        n_files  = n_events//max_events + 1 # have to add 1 due to floor division
                        batch_map[sample] = [0,n_events,n_files]
                        #
                    #
                #
            #
        #
        return batch_map


           
if __name__ == '__main__':
    
    condorSubmit()

