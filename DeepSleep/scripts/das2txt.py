#/cvmfs/cms.cern.ch/common/dasgoclient  -query="dataset=/ST_t*/*v4*upgrade2018*/NANO* instance=/prod/global"
import os
import sys
import re
import json
import numpy as np
import subprocess as sb

def main(cfg_file):
    year = re.search(r'201\d',cfg_file).group()
    #ofile_name = '/'.join(cfg_file.split('/')[:-1])+f'/sampleDas_nano_{year}.cfg'
    ofile_name = '/'.join(cfg_file.split('/')[:-1])+f'/sampleDas_nano_{year}.json'
    out_text = {}
    with open(cfg_file,'r') as cff:
        for l in cff.readlines():
            line = l.replace('\n','')
            if '#' in line or len(line.strip()) == 0: continue
            try:
                process, mini_name, *_ = line.split()
            except:
                print(f'line.split failed for {l}: {line}')
            nano_name = getDasNano(mini_name,year)
            if len(nano_name) == 0:
                print(f"\n{process} for year:{year} has no corresponding nanoAODv7!!!\n")
                nano_name = 'N/A'
            #out_text.append(f'{process}\t{nano_name}\n')
            out_text[f'{process}_{year}'] = nano_name.split('\t')
    #
    with open(ofile_name, 'w') as outf:
        print(f"\nWriting to: {ofile_name}\n")
        #outf.writelines(out_text)
        json.dump(out_text, outf, indent=4)
        
    
        
def getDasNano(mini_name,year):
    year_to_das = {'2016': 'asymptotic',
                   '2017': 'mc2017',
                   '2018': 'upgrade2018',
               }
    dasgo = "/cvmfs/cms.cern.ch/common/dasgoclient"
    # special stupid case where someone changes 2L2Nu to 2l2nu...
    if 'TTbb_4f_TTTo2L2Nu' in mini_name and year == '2018':
        mini_name = mini_name.replace('TTTo2L2Nu','TTTo2l2nu')
    #
    try:
        request = '/'+mini_name.split('/')[1] + f"/*AODv7*{year_to_das[year]}*/NANO*"
    except:
        print(f"request construction failed: {mini_name}")
    nano_name = sb.check_output(f"{dasgo} -query='dataset={request} instance=/prod/global'",
                                shell=True).decode()
    # special case where new_pmx is in the 2017 sample
    #if year == '2017' and ('new_pmx' in nano_name or 'PSweights' in nano_name):
    if ('new_pmx' in nano_name or 'PSweights' in nano_name):
        #only keep nanoAOD with new_pmx in name
        if 'PSweights' in nano_name and 'new_pmx' in nano_name:
            nano_name= '\n'.join(re.findall(r'/.*PSweights.*/.*new_pmx.*/NANOAODSIM',nano_name))
        elif 'PSweights' in nano_name and 'new_pmx' not in nano_name:
            nano_name= '\n'.join(re.findall(r'/.*PSweights.*/.*/NANOAODSIM',nano_name))
        else:
            nano_name= '\n'.join(re.findall(r'/.*/.*new_pmx.*/NANOAODSIM',nano_name))
        #print(nano_name)            
    #
    return nano_name.replace('\n','\t').rstrip('\t')
                            

if __name__ == '__main__':
    cfg_file = sys.argv[1]
    if not os.path.exists(sys.argv[1]) : raise NameError(sys.argv[1]) 
    main(cfg_file)
