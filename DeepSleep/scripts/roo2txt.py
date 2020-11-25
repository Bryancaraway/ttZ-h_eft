from glob import glob
import os
import re
import numpy as np

'''
Script to glob root file output from pre and post
processed jobs and place them in txt file
with the correct format to be read from 
lpc nodes
'''
def worker(step_str):

    nano_dir = '/cms/data/store/user/bcaraway/nanoAOD/'+f'{step_str}Processed/201*/*/'
    
    pro_files = glob(nano_dir)
    for f in pro_files:
        roofiles = np.array(glob(f+'**/*.root', recursive=True))
        # check for duplicates by checking prod201\dMC_NANO_x this part and _x.root 
        if step_str == 'post':
            roo_names = [roofile.split('/')[-1] for roofile in roofiles] # strips directory 
            prefix = re.findall(r'prod201\dMC_NANO_\d+' ,' '.join(roo_names)) 
            suffix = re.findall(r'_\d+.root' ,' '.join(roo_names))
            cat_roo_names = [i+j for i,j in zip(prefix,suffix)]
            unique_index = []
            unique_roos  = []
            for i,roo in enumerate(cat_roo_names):
                if roo not in unique_roos: 
                    unique_index.append(i)
                    unique_roos.append(roo)

            roofiles = roofiles[unique_index] # only keep the first occurance of duplicate 
            #(will have to check samples who's pre/post stats do not match)


        # need to replace '/cms/data/'  with root://kodiak-se.baylor.edu//
        roofiles = [roo.replace('/cms/data/', 'root://kodiak-se.baylor.edu//') + '\n' for roo in roofiles]
        #
     
        f2strip = f.split('_')[-1]
        year = re.search('/201\d/', f).group().strip('/')
        #
        format_map = {#'pre' : f.replace(f2strip,'')+year+'.txt',
                      'pre' : f.rstrip(f2strip)+year+'.txt',
                      'post': f.rstrip('/')+'.txt'}
        #
        txt_name = format_map[step_str]
        #
        with open(txt_name, 'w') as f_txt:
            f_txt.writelines(roofiles)
    
    #need to atr txt files for transfer
    years = ['2016','2017','2018']
    tar_dir = f'/cms/data/store/user/bcaraway/nanoAOD/{step_str}Processed/'
    for y in years:
        os.system(f'cd {tar_dir}{y}/; tar czvf sample_txt_files_{y}.tar.gz *.txt')

def main():
    # run script on both pre and post processed areas
    # its a super quick script anyways
    worker('pre')
    worker('post')

if __name__ == '__main__':
    
    main()
