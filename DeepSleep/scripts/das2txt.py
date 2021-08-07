#/cvmfs/cms.cern.ch/common/dasgoclient  -query="dataset=/ST_t*/*v4*upgrade2018*/NANO* instance=/prod/global"
'''
Script to add file for das entries in json output
'''
import os
import signal
import sys
import re
import json
import uproot
import numpy as np
import subprocess as sb
from multiprocessing import Pool

dasgo = "/cvmfs/cms.cern.ch/common/dasgoclient"

def das_from_cfg(cfg_file):
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
    check_other_end(ofile_name)
        
    
        
def getDasNano(mini_name,year):
    year_to_das = {'2016': 'asymptotic',
                   '2017': 'mc2017',
                   '2018': 'upgrade2018',
               }
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
                            

def check_other_end(json_file):
    new_json = {}
    with open(json_file) as jf:
        dict_ = json.load(jf)
        for k,v in dict_.items():
            print(k)
            new_json[k] = {'DAS':v}
            gfile_list = []
            bfile_list = []
            file_list = []
            for das in v:
                if 'N/A' in das: continue
                nano_name = sb.check_output(f"{dasgo} -query='file dataset={das} instance=/prod/global'", shell=True).decode()
                pool = Pool()
                #results = np.array(pool.map(check_xrdcp, nano_name.splitlines()))
                results = np.array(pool.map(check_uproot, nano_name.splitlines()))
                pool.close()
                #gfile_list += list(results[:,0][results[:,3] == "False"])
                #bfile_list += list(results[:,0][results[:,3] == "True"])
                gfile_list += list(results[:,0][results[:,1] == "True"])
                bfile_list += list(results[:,0][results[:,1] == "False"])
                file_list += list(results[:,0])
                #print(results[:,1])
                #print(results[:,0][results[:,1] == True])
                #print(results[:,0][results[:,1] == False])
            #new_json[k]['gfiles'] = gfile_list
            #new_json[k]['n_gfiles'] = len(gfile_list)
            #new_json[k]['bfiles'] = bfile_list
            #new_json[k]['bfiles'] = len(bfile_list)
            new_json[k]['files'] = file_list
            new_json[k]['n_files'] = len(file_list)
            #
        #
    #
    with open(json_file.replace('.json','_v2.json'),'w') as jf:
        print(f"\nWriting to: {jf}\n")
        #outf.writelines(out_text)
        json.dump(new_json, jf, indent=4)
            

                
                    
def check_uproot(roo):
    isgood = True
    try:
        with uproot.open(f'root://cmsxrootd.fnal.gov/{roo}') as _:
            return (roo, isgood)
    except:
        isgood = False
        return (roo, isgood)

def check_xrdcp(roo):
    out = err = ''
    with sb.Popen(f"xrdcp -d 1 -f root://cmsxrootd.fnal.gov/{roo} -> /dev/null", shell=True, stdout=sb.PIPE, stderr=sb.PIPE, preexec_fn=os.setsid) as pro:
        try:
            out, err = pro.communicate(timeout=60)
        except:
            os.killpg(pro.pid, signal.SIGINT)
            out, err = pro.communicate()
    return (roo, out, err, bool(len(err) == 0 or 'ERROR' in err.decode()) )

if __name__ == '__main__':
    if not os.path.exists(sys.argv[1]) : raise NameError(sys.argv[1]) 
    i_file = sys.argv[1]
    das_from_cfg(i_file) # pass mini --> get nano
    #check_other_end(i_file) # pass nano --> get nano_v2
